#include "dogdefs.h"
#include "tensors.h"
#include "IniParams.h"

// Function that is called after a full time step (i.e., after all stages are complete)
//
// This routine exists because the domain for this problem is weird.  Here, we
// reset the values of the function to be constant inside the wedge.  Next
// time ConstructL is called, these values get overwritten.
void AfterFullTimeStep(double dt,
		       dTensorBC3& auxold, dTensorBC3& aux, 
		       dTensorBC3& qold,   dTensorBC3& q)
{
    const int   mx   = q.getsize(1);
    const int   my   = q.getsize(2);
    const int meqn   = q.getsize(3);
    const int mbc    = q.getmbc();
    const int maux   = aux.getsize(3);

    const double dx   = global_ini_params.get_dx();
    const double dy   = global_ini_params.get_dy();

    if( mx%5 != 0 || my%5 != 0 )
    {
        printf("Error: mx, my = %d, %d\n", mx, my );
        printf("In order to run this problem, please select mx and my to be a multiple of 5.\n");
        exit(1);
    }

    const int istep = (mx / 5);
    const int jstep = (my / 5)+1;

    // Reset the domain inside the wedge to be the initial conditions.
    //
    // This will be overwritten next time ConstructL is called.
    void QinitFunc(const dTensor2& xpts,
        const dTensor2& NOT_USED1, const dTensor2& NOT_USED2,
        dTensor2& qvals);
//  void SampleFunction( 
//      int istart, int iend,
//      int jstart, int jend,
//      const dTensorBC3& qin, 
//      const dTensorBC3& auxin,  
//            dTensorBC3& Fout,
//      void (*Func)(const dTensor2&, const dTensor2&, const dTensor2&, dTensor2&));
//  SampleFunction( istep+1, mx+mbc, 1-mbc, jstep-1, q, aux, q, &QinitFunc );

    for( int i = istep+1; i <= mx+mbc; i++)
    for( int j = jstep-1; j >= 1-mbc; j-- )
    {
        for (int m=1; m<=meqn; m++)
        {

            q.set(i,j,m, 0.0 );
        }
    }


}
