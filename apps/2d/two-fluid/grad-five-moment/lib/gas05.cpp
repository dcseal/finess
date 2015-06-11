#include <cmath>
#include "dog_math.h"
#include "dogdefs.h"
#include "IniParams.h"
#include "gas05.h"

// five-moment gas routines

// cons and prim may reference the same array
void gas05_convert_cons_to_prim(int moffset, const dTensor1& cons, dTensor1& prim)
{

    using namespace FiveMomentComponentID;  // see gas05.h
    const int rho_ = moffset+_rho;
    const int M1_  = moffset+_M1;
    const int M2_  = moffset+_M2;
    const int M3_  = moffset+_M3;
    const int P_   = moffset+_N;

    const double rho = cons.get(rho_);
    const double rho_inv = 1./rho;
    const double M1  = cons.get(M1_);
    const double M2  = cons.get(M2_);
    const double M3  = cons.get(M3_);
    const double u1  = M1*rho_inv;
    const double u2  = M2*rho_inv;
    const double u3  = M3*rho_inv;
    const double N   = cons.get(P_);

    prim.set(rho_, rho);
    prim.set(M1_, u1);
    prim.set(M2_, u2);
    prim.set(M3_, u3);
    const double Mu = M1*u1+M2*u2+M3*u3;
    prim.set(P_, (1./3.)*(2*N - Mu));

}

// cons and prim may reference the same array
void gas05_convert_cons_to_prim(int moffset, const dTensor2& cons, dTensor2& prim)
{

    using namespace FiveMomentComponentID;
    const int rho_ = moffset+_rho;
    const int M1_  = moffset+_M1;
    const int M2_  = moffset+_M2;
    const int M3_  = moffset+_M3;
    const int P_   = moffset+_N;

    const int numpts = prim.getsize(1);
    assert_eq(cons.getsize(1),numpts);
    for(int i=1;i<=numpts;i++)
    {
        const double rho = cons.get(i, rho_);
        const double rho_inv = 1./rho;
        const double M1  = cons.get(i, M1_);
        const double M2  = cons.get(i, M2_);
        const double M3  = cons.get(i, M3_);
        const double u1  = M1*rho_inv;
        const double u2  = M2*rho_inv;
        const double u3  = M3*rho_inv;
        const double N = cons.get(i, P_);

        prim.set(i, rho_, rho);
        prim.set(i, M1_, u1);
        prim.set(i, M2_, u2);
        prim.set(i, M3_, u3);
        const double Mu = M1*u1+M2*u2+M3*u3;
        prim.set(i, P_, (1./3.)*(2*N - Mu));
    }

}

// prim and cons may reference the same array
void gas05_convert_prim_to_cons(int moffset, const dTensor1& prim, dTensor1& cons)
{

    using namespace FiveMomentComponentID;
    const int rho_ = moffset+_rho;
    const int M1_  = moffset+_M1;
    const int M2_  = moffset+_M2;
    const int M3_  = moffset+_M3;
    const int N_   = moffset+_N;

    const double rho = prim.get(rho_);
    const double u1  = prim.get(M1_);
    const double u2  = prim.get(M2_);
    const double u3  = prim.get(M3_);
    const double M1  = u1*rho;
    const double M2  = u2*rho;
    const double M3  = u3*rho;
    const double p   = prim.get(N_);

    cons.set(rho_, rho);
    cons.set(M1_, M1);
    cons.set(M2_, M2);
    cons.set(M3_, M3);

    const double Mu = M1*u1+M2*u2+M3*u3;
    cons.set(N_, 0.5*(3*p + Mu));

}

// prim and cons may reference the same array
void gas05_convert_prim_to_cons(int moffset,
        const dTensor2& prim, dTensor2& cons)
{
    using namespace FiveMomentComponentID;
    const int rho_ = moffset+_rho;
    const int M1_  = moffset+_M1;
    const int M2_  = moffset+_M2;
    const int M3_  = moffset+_M3;
    const int N_   = moffset+_N;

    const int numpts = prim.getsize(1);
    assert_eq(cons.getsize(1),numpts);
    for(int i=1;i<=numpts;i++)
    {

        const double rho = prim.get(i, rho_);
        const double u1  = prim.get(i, M1_);
        const double u2  = prim.get(i, M2_);
        const double u3  = prim.get(i, M3_);
        const double M1  = u1*rho;
        const double M2  = u2*rho;
        const double M3  = u3*rho;
        const double P = prim.get(i, N_);
        //
        cons.set(i, rho_, rho);
        cons.set(i, M1_, M1);
        cons.set(i, M2_, M2);
        cons.set(i, M3_, M3);
        const double Mu = M1*u1+M2*u2+M3*u3;
        cons.set(i, N_, 0.5*(3*P + Mu));
    }
}

void gas05_get_extreme_wave_spds( const dTensor1& nhat, const dTensor1& Ql, const dTensor1& Qr,
        int moffset, double *gas_min_speed, double *gas_max_speed)
{

    using namespace FiveMomentComponentID;
    const double& n1 = nhat.get(1);
    const double& n2 = nhat.get(2);
    const double nmag2 = n1*n1 + n2*n2;
//  assert_almost_eq(nmag2,1.);     // TODO - should probably add this check
                                    //  back! Does not compile on Mac machines -DS (6/11/2015).

    // Gas constant
    const double& gamma = global_ini_params.get_gamma();

    // Left states
    //
    const double rhol    = Ql.get(moffset+_rho);
    const double rhol_inv = 1./rhol;
    const double u1l     = Ql.get(moffset+_M1)*rhol_inv;
    const double u2l     = Ql.get(moffset+_M2)*rhol_inv;
    const double u3l     = Ql.get(moffset+_M3)*rhol_inv;
    const double energyl = Ql.get(moffset+_N);
    const double KEl = 0.5e0*rhol*(u1l*u1l+u2l*u2l+u3l*u3l);
    const double pressl  = (gamma-1.0e0)*(energyl-KEl);

    // Right states
    //
    const double rhor    = Qr.get(moffset+_rho);
    const double rhor_inv = 1./rhor;
    const double u1r     = Qr.get(moffset+_M1)*rhor_inv;
    const double u2r     = Qr.get(moffset+_M2)*rhor_inv;
    const double u3r     = Qr.get(moffset+_M3)*rhor_inv;
    const double energyr = Qr.get(moffset+_N);
    const double KEr  = 0.5e0*rhor*(u1r*u1r+u2r*u2r+u3r*u3r);
    const double pressr  = (gamma-1.0e0)*(energyr-KEr);

    // Average states
    //
    const double rho     = 0.5e0*(rhol+rhor);
    const double u1      = 0.5e0*(u1l+u1r);
    const double u2      = 0.5e0*(u2l+u2r);
    const double u3      = 0.5e0*(u3l+u3r);
    const double press   = 0.5e0*(pressl+pressr);

    // Sound speeds
    // 
    const double cl = sqrt(fabs(gamma*pressl/rhol));
    const double cr = sqrt(fabs(gamma*pressr/rhor));
    const double c  = sqrt(fabs(gamma*press/rho));

    // normal velocities
    //
    const double un  = n1*u1  + n2*u2;
    const double unl = n1*u1l + n2*u2l;
    const double unr = n1*u1r + n2*u2r;

    // Minimum speed
    (*gas_min_speed) = Min(unl-cl, un-c);

    // Maximum speed
    (*gas_max_speed) = Max(unr+cr, un+c);
}

inline double _get_entropy_per_mass05(double rho, double p)
{ return 0.5*log(pow(p,3)/pow(rho,5)); }

inline double _get_entropy_per_vol05(double rho, double p)
{ return rho*_get_entropy_per_mass05(rho,p); }

inline double _get_entropy_per_mass05( double rho,double M1,double M2,double M3,double nrg)
{
    const double u1  = M1/rho;
    const double u2  = M2/rho;
    const double u3  = M3/rho;
    const double KE  = 0.5*(u1*M1+u2*M2+u3*M3);
    const double TE = nrg - KE;
    const double p = TE*(global_ini_params.get_gamma()-1.);
    return _get_entropy_per_mass05(rho,p);
}

inline double _get_entropy_per_vol05(double rho, double M1, double M2, double M3, double nrg)
{ return rho*_get_entropy_per_mass05(rho, M1, M2, M3, nrg); }

// non-inline versions
//
double get_entropy_per_vol05(double rho, double press)
{ return _get_entropy_per_vol05(rho, press); }
double get_entropy_per_mass05(double rho,double M1,double M2,double M3,double nrg)
{ return _get_entropy_per_mass05(rho,M1,M2,M3,nrg); }
double get_entropy_per_vol05(double rho,double M1,double M2,double M3,double nrg)
{ return _get_entropy_per_vol05(rho,M1,M2,M3,nrg); }

void set_entropy_per_vol05(dTensor2& q_new, const dTensor2& qvals, int firstIdx, int _entropy_s)
{
    // Loop over each point
    for (int i=1; i<=qvals.getsize(1); i++)
    {
        int idx=1;
        const double rho = qvals.get(i, idx++);
        const double M1  = qvals.get(i, idx++);
        const double M2  = qvals.get(i, idx++);
        const double M3  = qvals.get(i, idx++);
        const double nrg = qvals.get(i, idx++);
        const double entropy = _get_entropy_per_vol05(rho,M1,M2,M3,nrg);

        q_new.set(i, _entropy_s, entropy);
    }
}

/*
 * TODO - replace this to check point values in place of DG bassis stuff
 int gas05_check_positivity_cons(int moffset, const dTensor1& Q)
 {

 const double rho_epsilon  = 0.5*LimiterParams::get_rho_epsilon();
 const double epsilon      = 0.5*LimiterParams::get_pressure_epsilon();

 using namespace FiveMomentComponentID;
 using namespace PositivityReturnType;
 const double rho = Q.get(moffset+_rho);
 if(rho<=rho_epsilon)
 {
 dprintf("negative density:"
 "\n  moffset = %d,"
 "\n  rho = %24.16e",
 moffset, rho);
 return NEG_DENSITY;
 }

// confirm positivity indicators at cell average state
double pos[4];

const double M1  = Q.get(moffset+_M1);
const double M2  = Q.get(moffset+_M2);
const double M3  = Q.get(moffset+_M3);
const double N = Q.get(moffset+_N);
const double Mu = M1*M1+M2*M2+M3*M3/rho;
const double p = (1./3.)*(2*N - Mu);

if( p<= epsilon)
{
return NEG_PRESSURE;
}
return EVALUATIONS_CONFIRMED_POSITIVITY;
}
 */


