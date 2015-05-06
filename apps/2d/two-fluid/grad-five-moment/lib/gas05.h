#ifndef _gas05_h_
#define _gas05_h_
class dTensor1;
class dTensor2;

// see also TenMomentComponentID
namespace FiveMomentComponentID
{
  enum Enum{
    _rho=1,  // rho
    _M1 =2,  // 1-momentum
    _M2 =3,  // 2-momentum
    _M3 =4,  // 3-momentum
    _N  =5,  // energy_11

    _u1  = _M1,
    _u2  = _M2,
    _u3  = _M3,
  };
}

void gas05_convert_cons_to_prim(int moffset,
  const dTensor2& cons, dTensor2& prim);
void gas05_convert_prim_to_cons(int moffset,
  const dTensor2& prim, dTensor2& cons);
void gas05_convert_cons_to_prim(int moffset,
  const dTensor1& cons, dTensor1& prim);
void gas05_convert_prim_to_cons(int moffset,
  const dTensor1& prim, dTensor1& cons);

void gas05_get_extreme_wave_spds(
  const dTensor1& nhat,
  const dTensor1& Ql,
  const dTensor1& Qr,
  int moffset,
  double *gas_min_speed,
  double *gas_max_speed);

double get_entropy_per_vol05(double rho, double press);
double get_entropy_per_vol05(double rho,double M1,double M2,double M3,double N);
void set_entropy_per_vol05( dTensor2& q_new, const dTensor2& qvals,
    int firstIdx, int _entropy_s);
double get_entropy_per_mass05(double rho,double M1,double M2,double M3,double N);

#endif
