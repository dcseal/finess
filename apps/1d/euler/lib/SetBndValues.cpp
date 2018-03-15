#include "tensors.h"
#include "StateVars.h"

// Zeroth order extrapolation
void SetBndValues( StateVars& Q )
{

    dTensorBC2&  q  = Q.ref_q  ();
    dTensorBC2& aux = Q.ref_aux();
    double t        = Q.get_t  ();

    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
    const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (int i=0; i>=(1-mbc); i--)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(1,m);
            q.set(i,m, tmp );
        }

        // aux values
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(1,m);
            aux.set(i,m, tmp );
        }
    }
    // ***********************************************  


    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (int i=(melems+1); i<=(melems+mbc); i++)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(melems,m);
            q.set(i,m, tmp );
        }

        // aux values
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(melems,m);
            aux.set(i,m, tmp );
        }
    }
    // ***********************************************

}

// Zeroth order extrapolation
//
// This routine can be used to apply the same boundary conditions to anything
// passed into it
void SetBndValues( dTensorBC2& q )
{

//  dTensorBC2&  q  = Q.ref_q  ();
//  dTensorBC2& aux = Q.ref_aux();
//  double t        = Q.get_t  ();

    const int melems = q.getsize(1);
    const int meqn   = q.getsize(2);
    const int kmax   = q.getsize(3);
//  const int maux   = aux.getsize(2);
    const int mbc    = q.getmbc();

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (int i=0; i>=(1-mbc); i--)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(1,m);
            q.set(i,m, tmp );
        }

    }
    // ***********************************************  


    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (int i=(melems+1); i<=(melems+mbc); i++)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(melems,m);
            q.set(i,m, tmp );
        }

    }
    // ***********************************************

}
