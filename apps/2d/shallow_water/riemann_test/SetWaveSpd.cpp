#include <cmath>
#include "dog_math.h"
#include "tensors.h"
#include "stdio.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
void SetWaveSpd(const dTensor1& nvec, const dTensor1& xedge,
		const dTensor1& Ql, const dTensor1& Qr,
		const dTensor1& Auxl, const dTensor1& Auxr,
		double& s1,double& s2)
{

    double hl,u1l,u2l,unl;
    double hr,u1r,u2r,unr;
    double h,u1,u2,un;
    
    // Left states
    hl  = Ql.get(1);
    u1l = Ql.get(2)/hl;
    u2l = Ql.get(3)/hl;
    unl = u1l*nvec.get(1) + u2l*nvec.get(2);

    // Right states
    hr  = Qr.get(1);
    u1r = Qr.get(2)/hl;
    u2r = Qr.get(3)/hl;
    unr = u1r*nvec.get(1) + u2r*nvec.get(2);
        
    // Average states
    h  = 0.5*(hl+hr);
    u1 = (sqrt(hl)*u1l + sqrt(hr)*u1r)/(sqrt(hl)+sqrt(hr));
    u2 = (sqrt(hl)*u2l + sqrt(hr)*u2r)/(sqrt(hl)+sqrt(hr));
    un = u1*nvec.get(1) + u2*nvec.get(2);

    // Minimum speed
    s1 = Min(un-sqrt(h), unl-sqrt(hl));
    
    // Maximum speed
    s2 = Max(un+sqrt(h), unr+sqrt(hr));

//  printf("Printing values from SetWaveSpd\n");
//  printf("  un, u1, u2 = %f %f %f\n", un, u1,u2 );
//  printf("  n1, n2 = %f, %f\n", nvec.get(1), nvec.get(2) );
//  printf("  h, hl, hr = %f %f %f\n", h, hl, hr );
//  printf("s1 = %f, s2 = %f\n", s1, s2 );
    
}

void MaxEig( const dTensor1& nvec, const dTensor1& xedge,
    const dTensor1& Q, const dTensor1& Aux, double& eig )
{

    // Average states
    double h      = Q.get(1);
    double sqrth  = sqrt( h );
    double u1     = Q.get(2)/h;
    double u2     = Q.get(3)/h;
    double u      = u1*nvec.get(1) + u2*nvec.get(2);

    eig = Max( fabs(u-sqrth), fabs(u+sqrth) );

}
