#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "StateVars.h"

using namespace std;

void FinSolveUser( StateVars& Qnew, double tend, double dtv[] )
{


/* TODO - write a template here, and pull it from FinSolveRK.  See the
 * template already created in the 1D code.  (-DS)

    while (t<tend)
    {
        // initialize time step
        m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveUser.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }        

        // copy qnew into qold
        CopyQ(qnew,qold);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {
            // set current time
            told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // Set initial maximum wave speed to zero
            for (j=1-mbc; j<=(melems+mbc); j++)
            { smax.set(j, 0.0e0 ); }

            // ----------------------------------------------------------------
            //
            //    THIS IS WHERE THE USER-DEFINED TIME-STEPPING SCHEME
            //    SHOULD BE ADDED. IN THE DEFAULT FILE: DogSolveUser.cpp,
            //    THE PROGRAM WILL NOW RETURN AN ERROR MESSAGE.
            // 
            // ----------------------------------------------------------------
            cout << endl;
            cout << " No user-defined time-stepping scheme has been defined yet. " << endl;
            cout << " Copy $FINESS/lib/1d/DogSolveUser.cpp into the current " << endl;
            cout << " directory and modify as needed." << endl << endl;
            exit(1);
            // ----------------------------------------------------------------

            // do any extra work      
            AfterFullTimeStep(dt,node,prim_vol,auxstar,aux,qold,qnew);

            // compute cfl number
            cfl = GetCFL(dt,dtv[2],prim_vol,method,aux,smax);

            // output time step information
            if (method[4]>0) 
            {
                cout << setprecision(3);
                cout << "DogSolve1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
                // accept
            { m_accept = 1; }
            else 
                //reject
            {   
                t = told;
                if (method[4]>0)
                {
                    cout<<"DogSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // copy qold into qnew
                CopyQ(qold,qnew);
            }

        }

        // compute conservation and print to file
        ConSoln(method,node,aux,qnew,t,outputdir);

    }

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

*/

}
