#include "tensors.h"
#include "StateVars.h"

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
    const StateVars& Qstar, const dTensorBC4& Lstar, StateVars& Qnew)
{

    const dTensorBC4& qstar   = Qstar.const_ref_q();
    const dTensorBC4& auxstar = Qstar.const_ref_aux();
    double tstar              = Qstar.get_t  ();

    dTensorBC4&  qnew   = Qnew.ref_q  ();
    dTensorBC4& auxnew  = Qnew.ref_aux();
    double tnew         = Qnew.get_t  ();

    // Update time
    Qnew.set_t( alpha1*tstar + alpha2*tnew + beta*dt );

    const int numel = qnew.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {
        double tmp = alpha1*qstar.vget(k)+alpha2*qnew.vget(k)+beta*dt*Lstar.vget(k);
        qnew.vset( k, tmp );
    }


}

// Update the solution using the constructed Lstar
//
// This version of UpdateSoln is used for the fifth-order time stepping.
void UpdateSoln(
    double g1, double g2, double g3, double delta, 
    double beta, double dt,
    const StateVars& Qold, const dTensorBC4& Lstar,
    StateVars& Q1, StateVars& Q2)
{

    const dTensorBC4& qold    = Qold.const_ref_q();
    const dTensorBC4& auxold  = Qold.const_ref_aux();

    dTensorBC4&  q1   = Q1.ref_q  ();
    dTensorBC4& aux1  = Q1.ref_aux();

    dTensorBC4&  q2   = Q2.ref_q  ();
    dTensorBC4& aux2  = Q2.ref_aux();

    const int     mx = q1.getsize(1);
    const int     my = q1.getsize(2);
    const int   meqn = q1.getsize(3);
    const int   maux = aux1.getsize(3);
    const int    mbc = q1.getmbc();

    const int numel = q1.numel();
#pragma omp parallel for
    for( int k=0; k < numel; k++ )
    {

        // update q2
        const double s1 = q1.vget( k );
        const double s2 = q2.vget( k ) + delta*s1;
        q2.vset(k, s2 );

        // update q
        const double s3  = qold.vget( k );
        const double tmp = g1*s1 + g2*s2 + g3*s3 + beta*dt*Lstar.vget(k);
        q1.vset(k, tmp );

    }

}
