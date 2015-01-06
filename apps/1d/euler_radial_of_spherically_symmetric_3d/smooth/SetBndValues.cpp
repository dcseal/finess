#include "dogdefs.h"
#include "StateVars.h"

// Periodic boundary conditions
//
void SetBndValues( StateVars& Q )
{

    dTensorBC2&  q  = Q.ref_q  ();
    dTensorBC2& aux = Q.ref_aux();
    double t        = Q.get_t  ();

    double tmp;
    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();
    
    { 
        // ***********************************************
        // LEFT BOUNDARY
        // ***********************************************
        for (int i=0; i>=(1-mbc); i--)
        {        
	        q.set(i,1,  q.get(1-i, 1));
            q.set(i,2, -q.get(1-i, 2));
            q.set(i,3,  q.get(1-i, 3));
                
        }
        // ***********************************************  


        // ***********************************************
        // RIGHT BOUNDARY
        // ***********************************************
        for (int i=(mx+1); i<=(mx+mbc); i++)
        {        
	        // q values
	        for (int m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i-mx,m);    
		        q.set(i,m, tmp );
            }
                
		}
        // ***********************************************
    }
}
