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
	        // q values
	        for (int m=1; m<=meqn; m++)
	        {
	            tmp = q.get(i+mx,m);
		        q.set(i,m, tmp );
            }
                
	        // aux values
	        for (int m=1; m<=maux; m++)
            {
	            tmp = aux.get(i+mx,m);
		        aux.set(i,m, tmp );
            }
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
                
	        // aux values
	        for (int m=1; m<=maux; m++)
            {
	            tmp = aux.get(i-mx,m);
		        aux.set(i,m, tmp );
            }
		}
        // ***********************************************
    }
}
