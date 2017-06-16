/********************************************************
  Xpress-BCL C++ Example Problems
  ===============================

  file xbcutstk.cxx
  `````````````````
  Cutting stock problem, solved by column (= cutting 
  pattern) generation heuristic looping over the
  root node.

  (c) 2008 Fair Isaac Corporation
      author: S.Heipcke, 2001, rev. Mar. 2014
********************************************************/

#include <iostream>
#include <cmath>
#include "xprb_cpp.h"

using namespace std;
using namespace ::dashoptimization;

#define NWIDTHS 5     //number of demanded width
#define MAXWIDTH 94   //the max width of the raw material

#define EPS 1e-6
#define MAXCOL 10   //number of so-called raw material

/****DATA****/
double WIDTH[] = {17, 21, 22.5, 24, 29.5};  /* Possible widths */
int DEMAND[] = {150, 96, 48, 108, 227};     /* Demand per width */
int PATTERNS[NWIDTHS][NWIDTHS];              /* (Basic) cutting patterns, also known as the initial cutting patterns*/

XPRBvar pat[NWIDTHS + MAXCOL];               /* Rolls per pattern. Note MAXCOL here is no explicit meaning ,but to save space
 * Can prove new patterns is less than MAXCOL? */
//XPRBvar pat[100];
XPRBctr dem[NWIDTHS];                      /* Demand constraints */
XPRBctr cobj;                              /* Objective function */

XPRBprob p("CutStock");                    /* Initialize a new problem in BCL */

double knapsack(int N, double *c, double *a, double R, int *d, int *xbest);

/***********************************************************************/

void modCutStock() {
    int i, j;
    XPRBexpr le;

    for (j = 0; j < NWIDTHS; j++)
        PATTERNS[j][j] = (int) floor(MAXWIDTH / WIDTH[j]);

    /****VARIABLES****/
    for (j = 0; j < NWIDTHS; j++)
        pat[j] = p.newVar(XPRBnewname("pat_%d", j + 1), XPRB_UI, 0,
                          (int) ceil((double) DEMAND[j] / PATTERNS[j][j]));

    /****OBJECTIVE****/
    for (j = 0; j < NWIDTHS; j++)
        le += pat[j];    /* Minimize total number of rolls */
    cobj = p.newCtr("OBJ", le);
    p.setObj(cobj);

    /****CONSTRAINTS****/
    for (i = 0; i < NWIDTHS; i++)                  /* Satisfy the demand per width */
    {
        le = 0;
        for (j = 0; j < NWIDTHS; j++)
            le += PATTERNS[i][j] * pat[j];
        dem[i] = p.newCtr("Demand", le >= DEMAND[i]);
    }
}

/**************************************************************************/
/*  Column generation loop at the top node:                               */
/*    solve the LP and save the basis                                     */
/*    get the solution values                                             */
/*    generate new column(s) (=cutting pattern)                           */
/*    load the modified problem and load the saved basis                  */
/**************************************************************************/
void solveCutStock() {
    double objval;                  /* Objective value */
    int i, j;
    int starttime;
    int npatt, npass;               /* Counters for columns and passes */
    double solpat[NWIDTHS + MAXCOL];  /* Solution values for variables pat */
    double dualdem[NWIDTHS];        /* Dual values of demand constraints */
    XPRBbasis basis;
    double dw, z;
    int x[NWIDTHS];

    starttime = XPRB::getTime();
    npatt = NWIDTHS;   //initially set to the number of widths

    for (npass = 0; npass < MAXCOL; npass++) {
        p.lpOptimize("");              /* Solve the LP */
        basis = p.saveBasis();         /* Save the current basis */
        objval = p.getObjVal();        /* Get the objective value */

        /* Get the solution values: */
        for (j = 0; j < npatt; j++)
            solpat[j] = pat[j].getSol();
        for (i = 0; i < NWIDTHS; i++)
            dualdem[i] = dem[i].getDual();

        /* Solve integer knapsack problem  z = min{cx : ax<=r, x in Z^n}
           with r=MAXWIDTH, n=NWIDTHS */
        z = knapsack(NWIDTHS, dualdem, WIDTH, (double) MAXWIDTH, DEMAND, x);
        cout << "(" << (XPRB::getTime() - starttime) / 1000.0 << " sec) Pass " << npass + 1 << ": ";

        if (z < 1 + EPS) {
            cout << "no profitable column found." << endl << endl;
            basis.reset();                  /* No need to keep the basis any longer */
            break;
        } else {
            /* Print the new pattern: */
            cout << "new pattern found with marginal cost " << z - 1 << endl << "   ";
            cout << "Widths distribution: ";
            dw = 0;
            for (i = 0; i < NWIDTHS; i++) {
                cout << WIDTH[i] << ":" << x[i] << "  ";
                dw += WIDTH[i] * x[i];
            }
            cout << "Total width: " << dw << endl;

            /* Create a new variable for this pattern: */
            pat[npatt] = p.newVar(XPRBnewname("pat_%d", npatt + 1), XPRB_UI);

            cobj += pat[npatt];             /* Add new var. to the objective */
            dw = 0;
            for (i = 0; i < NWIDTHS; i++)        /* Add new var. to demand constraints*/
                if (x[i] > EPS) {
                    dem[i] += x[i] * pat[npatt];
                    if ((int) ceil((double) DEMAND[i] / x[i]) > dw)
                        dw = (int) ceil((double) DEMAND[i] / x[i]);
                }
            pat[npatt].setUB(dw);           /* Change the upper bound on the new var.*/

            npatt++;

            p.loadMat();                    /* Reload the problem */
            p.loadBasis(basis);             /* Load the saved basis */
            basis.reset();                  /* No need to keep the basis any longer */
        }
    }

    p.mipOptimize("");                /* Solve the MIP */

    cout << "(" << (XPRB::getTime() - starttime) / 1000.0 << " sec) Optimal solution: " << p.getObjVal() << " rolls, "
         << npatt << " patterns" << endl << "   ";
    cout << "Rolls per pattern: ";
    for (i = 0; i < npatt; i++)
        cout << pat[i].getSol() << ", ";
    cout << endl;
}

/**************************************************************************/
/* Integer Knapsack Algorithm for solving the integer knapsack problem    */
/*    z = max{cx : ax <= R, x <= d, x in Z^N}                             */
/*                                                                        */
/* Input data:                                                            */
/*   N:        Number of item types                                       */
/*   c[i]:     Unit profit of item type i, i=1..n                         */
/*   a[i]:     Unit resource use of item type i, i=1..n                   */
/*   R:        Total resource available                                   */
/*   d[i]:     Demand for item type i, i=1..n                             */
/* Return values:                                                         */
/*   xbest[i]: Number of items of type i used in optimal solution         */
/*   zbest:    Value of optimal solution                                  */
/**************************************************************************/
double knapsack(int N, double *c, double *a, double R, int *d, int *xbest) {
    int j;
    double zbest = 0.0;
    XPRBvar *x;
    XPRBexpr klobj, knap;
    XPRBprob pk("Knapsack");

/*
 cout << "Solving z = max{cx : ax <= b; x in Z^n}\n   c   =";
 for (j = 0; j < N; j++) cout << " " << c[j] << ",";
 cout << "\n   a   =";
 for (j = 0; j < N; j++) cout << " " << a[j] << ",";
 cout << "\n   c/a =";
 for (j = 0; j < N; j++) cout << " " << c[j]/a[j] << ",";
 cout << "\n   b   = " << R << endl;
*/

    x = new XPRBvar[N];
    if (x == NULL) cout << "Allocating memory for variables failed." << endl;
    for (j = 0; j < N; j++)
        x[j] = pk.newVar("x", XPRB_UI, 0, d[j]);

    for (j = 0; j < N; j++)
        klobj += c[j] * x[j];
    pk.setObj(pk.newCtr("OBJ", klobj));
    for (j = 0; j < N; j++)
        knap += a[j] * x[j];
    pk.newCtr("KnapXPRBctr", knap <= R);
    pk.setSense(XPRB_MAXIM);
    pk.mipOptimize("");

    zbest = pk.getObjVal();
    for (j = 0; j < N; j++)
        xbest[j] = (int) floor(x[j].getSol() + 0.5);
/*
 cout << "   z = " << zbest << endl << "   x =";
 for(j=0; j<N; j++) cout << " " << xbest[j] << ",";
 cout << endl;
*/

    delete[] x;

    return (zbest);
}

/***********************************************************************/

int main(int argc, char **argv) {
    modCutStock();                    /* Model the problem */
    solveCutStock();                  /* Solve the problem */

    return 0;
} 

