#include "tensors.h"
#include "StateVars.h"

// Zeroth order extrapolation boundary conditions
void SetBndValues(StateVars& Q)
{
    dTensorBC2&  q  = Q.ref_q  ();
    dTensorBC2& aux = Q.ref_aux();
    double t        = Q.get_t  ();

    int i,k,m,ell;
    double tmp;
    int melems = q.getsize(1);
    int meqn   = q.getsize(2);
    int kmax   = q.getsize(3);
    int maux   = aux.getsize(2);
    int mbc    = q.getmbc();
    void L2Project(int,int,int,dTensor2,dTensorBC3,dTensorBC3,dTensorBC3&,
                   void (*Func)(dTensor1, dTensor2, dTensor2, dTensor2&));
    double hl,hr,Qflow;
    
    hl = 1.0678715170821127;   //1.0040014495278291;
    hr = 0.88362682051110220;  //1.0040014495278291; 
    Qflow = 0.25;  //0.2;
     
    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (i=0; i>=(1-mbc); i--)
    {        
        // q values - height
        q.set(i,1, hl );

	// q values - momentum
	q.set(i,2, Qflow );

               
	// aux values
	for (ell=1; ell<=1; ell++)
	{
	    for (m=1; m<=1; m++)
	    {
		tmp = aux.get(1,1);
	    
		aux.set(i,1, tmp );
	    }
	}
    }
    // ***********************************************  


    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (i=(melems+1); i<=(melems+mbc); i++)
    {        
        // q values - height
        q.set(i,1, hr );

	// q values - momentum
	q.set(i,2, Qflow );

                
	// aux values
	for (ell=1; ell<=1; ell++)
	{
	    for (m=1; m<=1; m++)
	    {
		tmp = aux.get(melems,1);
	    
		aux.set(i,1, tmp );
	    }
	}
    }
    // ***********************************************

}
