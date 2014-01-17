#include <cmath>
#include "DogParams.h"
#include "DogParamsCart1.h"
#include "tensors.h"
#include "dog_math.h"

// Right-hand side for hyperbolic PDE in divergence form
//
//       q_t + f(q,x,t)_x = Psi(q,x,t)
//
void ConstructL(
        const dTensor2& node,
        dTensorBC2& aux,
        dTensorBC2& q,      // setbndy conditions modifies q
        dTensorBC2& Lstar,
        dTensorBC1& smax)
{

printf("ConstructL needs to be written.  This routine is doing nothing right now.\n");

/*
    void FluxFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
    void SetWaveSpd(const dTensor1&,const dTensor1&,const dTensor1&,const dTensor1&,
            const dTensor1&,double&,double&);
    void SourceTermFunc(const dTensor1&,const dTensor2&,const dTensor2&,dTensor2&);
*/

    const int     mx = q.getsize(1);
    const int   meqn = q.getsize(2);
    const int   maux = aux.getsize(2);
    const int    mbc = q.getmbc();

    // "Plus" flux, "minus" flux and source term, respectively
    dTensorBC2  Fp(mx, meqn, mbc );
    dTensorBC2  Fm(mx, meqn, mbc );
    dTensorBC2 Psi(mx, meqn, mbc );

    // Grid spacing -- node( 1:(mx+1), 1 ) = cell edges
    const double   xlow = dogParamsCart1.get_xlow();
    const double     dx = dogParamsCart1.get_dx();

    // ---------------------------------------------------------
    // Part I: compute inter-element interaction fluxes
    // ---------------------------------------------------------

    // Boundary conditions
    //
    // TODO - this should be moved before ConstructL is called (-DS)
    void SetBndValues(const dTensor2&, dTensorBC2&, dTensorBC2&);
    SetBndValues(node, aux, q);

/*
    // Loop over interior edges and solve Riemann problems
    //#pragma omp parallel for
    //
    // This sets both smax(i) and smax(i-1), so can't be parallelized!
    //
    for (int i=(1-mbc); i<=(mx+mbc); i++)
    {
        dTensor1 Ql(meqn);
        dTensor1 Qr(meqn);
        dTensor1 Auxl(meqn);
        dTensor1 Auxr(meqn);

        // Riemann data - q
        for (int m=1; m<=meqn; m++)
        {
            Ql.set(m, 0.0 );
            Qr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Ql.set(m, Ql.get(m) + sqrt(2.0*double(k)-1.0)
                        *q.get(i-1,m,k) );
                Qr.set(m, Qr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*q.get(i,m,k) ); 
            }
        }

        // Riemann data - aux
        for (int m=1; m<=maux; m++)
        {

            Auxl.set(m, 0.0 );
            Auxr.set(m, 0.0 );

            for (int k=1; k<=method[1]; k++)
            {
                Auxl.set(m, Auxl.get(m) + sqrt(2.0*double(k)-1.0)
                        *aux.get(i-1,m,k) );
                Auxr.set(m, Auxr.get(m) + pow(-1.0,k+1)
                        *sqrt(2.0*double(k)-1.0)*aux.get(i,m,k) ); 
            }
        }

        // Solve Riemann problem
        dTensor1 xedge(1);
        xedge.set(1, xlower + (double(i)-1.0)*dx );

        dTensor1 Fl(meqn);
        dTensor1 Fr(meqn);
        double smax_edge = RiemannSolve(xedge,Ql,Qr,Auxl,Auxr,Fl,Fr,
                &FluxFunc,&SetWaveSpd);

        // This is a problem for the pragma statements! (-DS)
        smax.set(i-1, Max(smax_edge,smax.get(i-1)) );
        smax.set(i,   Max(smax_edge,smax.get(i)) );

        // Construct fluxes
        for (int m=1; m<=meqn; m++)
        {
            Fm.set(i,  m, Fr.get(m) );
            Fp.set(i-1,m, Fl.get(m) );
        }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part II: compute source term
    // --------------------------------------------------------- 
    if ( method[7]>0 )
    {        
        // Set source term on computational grid
        // Set values and apply L2-projection
        L2Project(0,1-mbc,melems+mbc,node,q,aux,Psi,&SourceTermFunc);
    }
    // ---------------------------------------------------------


    // ---------------------------------------------------------
    // Part III: construct Lstar
    // ---------------------------------------------------------
    if( method[7]==0 )  // Without Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)
                for (int k=1; k<=method[1]; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx;

                    Lstar.set(i,m,k, tmp );	      
                }
    }
    else  // With Source Term
    {
#pragma omp parallel for
        for (int i=(2-mbc); i<=(melems+mbc-1); i++)	
            for (int m=1; m<=meqn; m++)	    
                for (int k=1; k<=method[1]; k++)
                {
                    double tmp = N.get(i,m,k) - sqrt(2.0*double(k)-1.0)*
                        ( Fp.get(i,m) + pow(-1.0,k)*Fm.get(i,m) )/dx
                        + Psi.get(i,m,k);

                    Lstar.set(i,m,k, tmp );
                }
    }
    // ---------------------------------------------------------

    // ---------------------------------------------------------
    // Part IV: add extra contributions to Lstar
    // ---------------------------------------------------------
    // Call LstarExtra
    // LstarExtra(node,aux,q,Lstar);
    // ---------------------------------------------------------

*/

}
