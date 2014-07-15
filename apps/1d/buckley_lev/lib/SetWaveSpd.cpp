#include "tensors.h"
#include "dog_math.h"

// This is a user-supplied routine that sets the
// HLLE wave speeds for use in "RiemannSolve"
//
void SetWaveSpd(const dTensor1& xedge,
		const dTensor1& Ql,
		const dTensor1& Qr,
		const dTensor1& Auxl,
                const dTensor1& Auxr,
		double& s1,double& s2)
{

    // TODO - make this a user supplied parameter!
    const double M = (1./3.);

    // Left wave speed:
    const double ql  = Ql.get(1);
    const double denl = ql*ql + M*(1.0-ql)*(1.0-ql);
    const double fpl   = 2.0*M*ql*(1.0-ql)/ (denl*denl);

    // Wave speed based on Roe (not really) average:
    const double qavg     = 0.5*( Ql.get(1) + Qr.get(1) );
    const double den   = qavg*qavg + M*(1.0-qavg)*(1.0-qavg);
    const double fp    = 2.0*M*qavg*(1.0-qavg)/ (den*den);

    // Minimum speed:
    s1 = Min( fpl, fp );

    // Maximum speed:
//  s2 = Max( fpl, fp );

    // We'll use global values for the flux speeds
    //
    // The fastest wave speed in the entire system is the right-going shock,
    // which travels at speed 1.5.
    //
    // There are no wave speeds traveling to the left.
    // s1 = 0.;
    s2 = 2.205737063;

}
