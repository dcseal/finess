#include "tensors.h"
#include "StateVars.h"

// Zeroth order extrapolation boundary conditions
void SetBndValues(StateVars& Q)
{
    dTensorBC2&  q  = Q.ref_q  ();
    dTensorBC2& aux = Q.ref_aux();
    double t        = Q.get_t  ();

    int i,m;
    double tmp;
    int melems = q.getsize(1);
    int meqn   = q.getsize(2);
    int maux   = aux.getsize(2);
    int mbc    = q.getmbc();

    { 
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (i=0; i>=(1-mbc); i--)
        {        
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = q.get(1,m);

                q.set(i,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = aux.get(1,m);

                aux.set(i,m, tmp );
            }
        }
        // ***********************************************  


        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (i=(melems+1); i<=(melems+mbc); i++)
        {        
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = q.get(melems,m);

                q.set(i,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = aux.get(melems,m);

                aux.set(i,m, tmp );
            }
        }
        // ***********************************************

    }

}
