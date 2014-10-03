#include "tensors.h"
#include "StateVars.h"
#include "IniParams.h"

// Update the solution using the constructed Lstar
//
// The Low-Storage RK methods use a combination of two different values of q,
// together with a right hand side, L.  This routine computes the following:
//
//     qnew = alpha1*qstar + alpha2*qnew + beta*dt*Lstar.
//
// The main loop covers all elements (including boundary cells) of the arrays.
//
// See also: RKinfo.h and SetRKinfo.cpp
void UpdateSoln(double alpha1, double alpha2, double beta, double dt, 
    const StateVars& Qstar, const dTensorBC2& Lstar, StateVars& Qnew )
{

    const dTensorBC2& qstar   = Qstar.const_ref_q();
    const dTensorBC2& auxstar = Qstar.const_ref_aux();
    double tstar              = Qstar.get_t  ();

    dTensorBC2&  qnew   = Qnew.ref_q  ();
    dTensorBC2& auxnew  = Qnew.ref_aux();
    double tnew         = Qnew.get_t  ();

    const int numel = qnew.numel();

    // Update time
    Qnew.set_t( alpha1*tstar + alpha2*tnew + beta*dt );

    // Update the conserved variables
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double tmp = alpha1*qstar.vget(k)+alpha2*qnew.vget(k)+beta*dt*Lstar.vget(k);
        qnew.vset( k, tmp );
    }

//  if( aux.get_size(2) > 0 )
//  {

//      // Update the Aux arrays ( TODO - THIS NEEDS TO CHANGE! )
//  #pragma omp parallel for
//          for( int k=0; k < numel; k++ )
//          {
//              double tmp = alpha1*qstar.vget(k)+alpha2*qnew.vget(k)+beta*dt*Lstar.vget(k);
//              qnew.vset( k, tmp );
//          }
//  }

}

// Update the solution using the constructed Lstar
//
// This version of UpdateSoln is used for the fifth-order time stepping.
void UpdateSoln(
    double g1,double g2, double g3, double delta, 
    double beta,double dt,
    const dTensorBC2& aux,
    const dTensorBC2& qold, const dTensorBC2& Lstar,
    dTensorBC2& q1, dTensorBC2& q2)
{
    const int numel = q1.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double s1 = q1.vget( k );
        double s3 = qold.vget( k );

        // update q2
        double s2 = q2.vget( k ) + delta*s1;
        q2.vset(k, s2 );

        // update q
        double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.vget(k);
        q1.vset(k, tmp );

    }

}
