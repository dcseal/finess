///@file apps/1d/euler/harten_shock_tube/SetBndValues.cpp

#include "tensors.h"

///@brief Zeroth order extrapolation boundary conditions
///
///In this problem, SetBndValues(...) simply copies q-values at each of two ends to ghost cells on that end.
void SetBndValues(
        const dTensor2& node, 
        dTensorBC2& aux, 
        dTensorBC2& q)
{

    const int mx     = q.getsize(1);
    const int meqn   = q.getsize(2);
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
    for (int i=(mx+1); i<=(mx+mbc); i++)
    {        
        // q values
        for (int m=1; m<=meqn; m++)
        {
            double tmp = q.get(mx,m);

            q.set(i,m, tmp );
        }

        // aux values
        for (int m=1; m<=maux; m++)
        {
            double tmp = aux.get(mx,m);

            aux.set(i,m, tmp );
        }
    }
    // ***********************************************

}
