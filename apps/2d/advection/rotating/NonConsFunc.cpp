#include "dogdefs.h"

// This is a user-supplied routine that sets the
// flux function for use in "RiemannSolve"
//
//     Simple advection equation
//
void NonConsFunc(const dTensor2& xpts, const dTensor2& Q, 
		 const dTensor2& Aux, dTensor2& fout)
{
    int i;
    int numpts=xpts.getsize(1);
    double x,y,qc,u,v,divu;

    for (i=1; i<=numpts; i++)
    {
        x = xpts.get(i,1);
	y = xpts.get(i,2);
	
	// Variables
	qc   = Q.get(i,1);
	u    = Aux.get(i,1);
	v    = Aux.get(i,2);
	divu = Aux.get(i,3);
	
        // (u_x + v_y)*q
        fout.set(i,1, divu * qc );
    }

}
