/********************************************************
  Xpress-BCL C++ Example Problems
  ===============================

  file xbels.cxx
  ``````````````
  Economic lot sizing, ELS, problem, solved by adding
  (l,S)-inequalities) in several rounds looping over 
  the root node.
  
  ELS considers production planning over a horizon
  of T periods. In period t, t=1,...,T, there is a
  given demand DEMAND[t] that must be satisfied by
  production prod[t] in period t and by inventory
  carried over from previous periods. There is a
  set-up up cost SETUPCOST[t] associated with
  production in period t. The unit production cost
  in period t is PRODCOST[t]. There is no inventory
  or stock-holding cost.

  (c) 2008 Fair Isaac Corporation
      author: S.Heipcke, 2001, rev. Mar. 2011
********************************************************/

#include <iostream>
#include "xprb_cpp.h"
#include "xprs.h"

using namespace std;
using namespace ::dashoptimization;

#define EPS    1e-6

#define T 6                             /* Number of time periods */

/****DATA****/
int DEMAND[] = {1, 3, 5, 3, 4, 2};  /* Demand per period */
int SETUPCOST[] = {17, 16, 11, 6, 9, 6};  /* Setup cost per period */
int PRODCOST[] = {5, 3, 2, 1, 3, 1};  /* Production cost per period */
//第一维是time period，第二维是从第一维开始到第二维的累加需求
int D[T][T];                            /* Total demand in periods t1 - t2 */

XPRBvar prod[T];                        /* Production in period t */
XPRBvar setup[T];                       /* Setup in period t */

XPRBprob p("Els");                      /* Initialize a new problem in BCL */

/***********************************************************************/

void modEls() {
    int s, t, k;
    XPRBexpr cobj, le;

    for (s = 0; s < T; s++)
        for (t = 0; t < T; t++)
            for (k = s; k <= t; k++)
                D[s][t] += DEMAND[k];

/****VARIABLES****/
    for (t = 0; t < T; t++) {
        prod[t] = p.newVar(XPRBnewname("prod%d", t + 1));
        setup[t] = p.newVar(XPRBnewname("setup%d", t + 1), XPRB_BV);
    }

/****OBJECTIVE****/
    for (t = 0; t < T; t++)                       /* Minimize total cost */
        cobj += SETUPCOST[t] * setup[t] + PRODCOST[t] * prod[t];
    p.setObj(cobj);

/****CONSTRAINTS****/
    /* Production in period t must not exceed the total demand for the
       remaining periods; if there is production during t then there
       is a setup in t */
    for (t = 0; t < T; t++)
        p.newCtr("Production", prod[t] <= D[t][T - 1] * setup[t]);

    /* Production in periods 0 to t must satisfy the total demand
       during this period of time */
    for (t = 0; t < T; t++) {
        le = 0;
        for (s = 0; s <= t; s++)
            le += prod[s];
        p.newCtr("Demand", le >= D[0][t]);
    }

}

/**************************************************************************/
/*  Cut generation loop at the top node:                                  */
/*    solve the LP and save the basis                                     */
/*    get the solution values                                             */
/*    identify and set up violated constraints                            */
/*    load the modified problem and load the saved basis                  */
/**************************************************************************/
void solveEls() {
    double objval;               /* Objective value */
    int t, l;
    int starttime;
    int ncut, npass, npcut;      /* Counters for cuts and passes */
    double solprod[T], solsetup[T];   /* Solution values for var.s prod & setup */
    double ds;
    XPRBbasis basis;
    XPRBexpr le;

    starttime = XPRB::getTime();
    XPRSsetintcontrol(p.getXPRSprob(), XPRS_CUTSTRATEGY, 0);
    /* Disable automatic cuts - we use our own */
    XPRSsetintcontrol(p.getXPRSprob(), XPRS_PRESOLVE, 0);
    XPRSsetintcontrol(p.getXPRSprob(), XPRS_MIPPRESOLVE, 0);
    XPRSsetintcontrol(p.getXPRSprob(), XPRS_PREPROBING, 0);
    /* Switch presolve off */
    ncut = npass = 0;

    do {
        npass++;
        npcut = 0;
        p.lpOptimize("p");          /* Solve the LP */
        basis = p.saveBasis();      /* Save the current basis */
        objval = p.getObjVal();     /* Get the objective value */

        /* Get the solution values: */
        for (t = 0; t < T; t++) {
            solprod[t] = prod[t].getSol();
            solsetup[t] = setup[t].getSol();
        }

        /* Search for violated constraints: */
        for (l = 0; l < T; l++) {
            for (ds = 0.0, t = 0; t <= l; t++) {
                if (solprod[t] < D[t][l] * solsetup[t] + EPS) ds += solprod[t];
                else ds += D[t][l] * solsetup[t];
            }

            /* Add the violated inequality: the minimum of the actual production
               prod[t] and the maximum potential production D[t][l]*setup[t]
               in periods 0 to l must at least equal the total demand in periods
               0 to l.
               sum(t=1:l) min(prod[t], D[t][l]*setup[t]) >= D[0][l]
             */
            if (ds < D[0][l] - EPS) {
                le = 0;
                for (t = 0; t <= l; t++) {
                    if (solprod[t] < D[t][l] * solsetup[t] + EPS)
                        le += prod[t];
                    else
                        le += D[t][l] * setup[t];
                }
                p.newCtr(XPRBnewname("cut%d", ncut + 1), le >= D[0][l]);
                ncut++;
                npcut++;
            }
        }

        cout << "Pass " << npass << " (" << (XPRB::getTime() - starttime) / 1000.0;
        cout << " sec), objective value " << objval << ", cuts added: " << npcut;
        cout << " (total " << ncut << ")" << endl;

        if (npcut == 0)
            cout << "Optimal integer solution found:" << endl;
        else {
            p.loadMat();                 /* Reload the problem */
            p.loadBasis(basis);          /* Load the saved basis */
            basis.reset();               /* No need to keep the basis any longer */
        }
    } while (npcut > 0);

    /* Print out the solution: */
    for (t = 0; t < T; t++) {
        cout << "Period " << t + 1 << ": prod " << prod[t].getSol() << " (demand: ";
        cout << DEMAND[t] << ", cost: " << PRODCOST[t] << "), setup ";
        cout << setup[t].getSol() << " (cost: " << SETUPCOST[t] << endl;
    }
}

/***********************************************************************/

int main(int argc, char **argv) {
    modEls();                      /* Model the problem */
    solveEls();                    /* Solve the problem */

    return 0;
} 
