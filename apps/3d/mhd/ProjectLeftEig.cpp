#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"
#include <iostream>
// This is a user-supplied routine that projects
// Qvals onto the left eigenvectors ofthe flux 
// Jacobian; the result is stored in Wvals
//
void ProjectLeftEig(int ixy, 
        const dTensor1& Aux_ave, 
        const dTensor1& Q_ave, 
        const dTensor2& Qvals, 
        dTensor2& Wvals)
{    

    const int meqn = Qvals.getsize(1);

    // Average states
    double const gamma  = global_ini_params.get_gamma();
    double rho    = Q_ave.get(1);
    double u1     = Q_ave.get(2)/rho;
    double u2     = Q_ave.get(3)/rho;
    double u3     = Q_ave.get(4)/rho;
    double energy = Q_ave.get(5);
    double B1     = Q_ave.get(6);
    double B2     = Q_ave.get(7);
    double B3     = Q_ave.get(8);
    double um2    = 0.5*(u1*u1 + u2*u2 + u3*u3);
    double Bm2    = 0.5*(B1*B1 + B2*B2 + B3*B3);
    double p      = (gamma-1.0)*(energy - rho*um2 - Bm2);
    
    

    dTensor1 nv(3), tv(3);
    {
        double Bmag;
        switch( ixy )
        {
            case 1:
                nv.set(1, 1.0 );
                nv.set(2, 0.0 );
                nv.set(3, 0.0 ); 

                Bmag = sqrt(B2*B2 + B3*B3);
                if (Bmag>1.0e-8)
                {
                    tv.set(1, 0.0     );
                    tv.set(2, B2/Bmag );
                    tv.set(3, B3/Bmag );
                }
                else
                {
                    tv.set(1, 0.0          );
                    tv.set(2, sin(0.25*pi) );
                    tv.set(3, cos(0.25*pi) );
                }
                break;

            case 2:
                nv.set(1, 0.0 );
                nv.set(2, 1.0 );
                nv.set(3, 0.0 );

                Bmag = sqrt(B1*B1 + B3*B3);
                if (Bmag>1.0e-8)
                {
                    tv.set(1, B1/Bmag );
                    tv.set(2, 0.0     );
                    tv.set(3, B3/Bmag );
                }
                else
                {
                    tv.set(1, sin(0.25*pi) );
                    tv.set(2, 0.0          );
                    tv.set(3, cos(0.25*pi) );
                }
                break;

            case 3:
                nv.set(1, 0.0 );
                nv.set(2, 0.0 );
                nv.set(3, 1.0 );

                Bmag = sqrt(B1*B1 + B2*B2);
                if (Bmag>1.0e-8)
                {
                    tv.set(1, B1/Bmag );
                    tv.set(2, B2/Bmag );
                    tv.set(3, 0.0     );
                }
                else
                {
                    tv.set(1, sin(0.25*pi) );
                    tv.set(2, cos(0.25*pi) );
                    tv.set(3, 0.0          );
                }
                break;
        }
    }

    double rhosq = sqrt(rho);
    double   gm1 = gamma-1.0;
    double    a2 = gamma*p/rho;
    double     a = sqrt(a2);
    double   b1s = B1/rhosq;
    double   b2s = B2/rhosq;
    double   b3s = B3/rhosq;
    double   BNs = b1s*nv.get(1)+b2s*nv.get(2)+b3s*nv.get(3);
    double  BNs2 = pow(BNs,2);
    double    BN = BNs*rhosq;
    double   BN2 = pow(BN,2);

    double   d = a2 + (b1s*b1s + b2s*b2s + b3s*b3s);
    double  d2 = pow(d,2);
    double  cf = sqrt(0.5*(d+sqrt(d2-4.0*a2*BNs2)));
    double  cs = sqrt(0.5*(d-sqrt(d2-4.0*a2*BNs2)));
    double cf2 = pow(cf,2);
    double cs2 = pow(cs,2);

    double alphaf,alphas;
    if ( fabs(cf*cf-cs*cs) <= 1.0e-12)
    {
        alphaf = 3.82683432365090e-1;
        alphas = 9.23879532511287e-1;
    }
    else
    {
        alphaf = sqrt(fabs(a2 - cs*cs))/sqrt(fabs(cf*cf-cs*cs));
        alphas = sqrt(fabs(cf*cf - a2))/sqrt(fabs(cf*cf-cs*cs));
    }
    double alphaf2 = pow(alphaf,2);
    double alphas2 = pow(alphas,2);

    double    psq = sqrt(p);
    double    sq2 = 1.4142135623730950488;
    double  beta1 = double( (BNs>=0.0) - (BNs<0.0) );
    double  sqgam = sqrt(gm1/gamma);
    double     BT = (B1*tv.get(1)+B2*tv.get(2)+B3*tv.get(3));
    double  gu1sq = 1.0/sqrt(gamma);
    double   TxN1 = nv.get(3)*tv.get(2)-nv.get(2)*tv.get(3);
    double   TxN2 = nv.get(1)*tv.get(3)-nv.get(3)*tv.get(1);
    double   TxN3 = nv.get(2)*tv.get(1)-nv.get(1)*tv.get(2);
    double     n1 = nv.get(1);
    double     n2 = nv.get(2);
    double     n3 = nv.get(3);
    double     t1 = tv.get(1);
    double     t2 = tv.get(2);
    double     t3 = tv.get(3);
    double    n1s = pow(nv.get(1),2);
    double    n2s = pow(nv.get(2),2);
    double    n3s = pow(nv.get(3),2);
    double    t1s = pow(tv.get(1),2);
    double    t2s = pow(tv.get(2),2);
    double    t3s = pow(tv.get(3),2);
    double    nen = n3s*(t1s+t2s) - 2.0*n1*n3*t1*t3 - 2.0*n2*t2*(n1*t1+n3*t3) + n2s*(t1s+t3s) + n1s*(t2s+t3s);
    double   nen2 = a*gu1sq*nen;
    double   nen3 = sq2*psq*nen;
    double  nen31 = nen3/rhosq;
    double Term51 = rho*cf*((n2*t1*t2+n3*t1*t3-n1*(t2s+t3s))*rhosq*cf2*alphaf
            -a*BN*alphas*(-n2s*t1+n1*n2*t2+n3*TxN2));
    double Term52 = rho*cf*((-t2*(n1*t1+n3*t3)+n2*(t1s+t3s))*rhosq*cf2
            *alphaf-a*BN*alphas*(-n1*n2*t1+n1s*t2+n3*TxN1));
    double Term53 = rho*cf*((n3*(t1s+t2s)-(n1*t1+n2*t2)*t3)*rhosq*cf2*
            alphaf-a*(-n1*n3*t1+n1s*t3+n2*(-n3*t2+n2*t3))*BN*alphas);
    double Term54 = alphaf/(sq2*a2*gu1sq*rhosq*(alphaf2+alphas2));
    double Term55 = (n2s*t1-n1*n2*t2+n3*(n3*t1-n1*t3))*alphas;
    double Term56 = alphas*(-n1*n2*t1+n1s*t2+n3*TxN1);
    double Term57 = (-n1*n3*t1+n1s*t3+n2*(-n3*t2+n2*t3))*alphas;
    double  nen51 = sq2*a*nen*gu1sq*(a*BN2*alphas2+rhosq*cf2*alphaf*
            (a*rhosq*alphaf+BT*alphas));
    double  nen52 = sq2*a*gu1sq*(alphaf2+alphas2)*nen;
    double Term71 = rho*cf*(a*(n2s*t1-n1*n2*t2+n3*(n3*t1-n1*t3))*rhosq*
            alphaf+((B3*n2*t1-B2*n3*t1-B3*n1*t2+B2*n1*t3)
                *(-n3*t2+n2*t3)+B1*(n2s*t1s+n3s*t1s
                    -2.0*n1*n2*t1*t2-2.0*n1*n3*t1*t3+n1s*(t2s+t3s)))*alphas);
    double Term72 = rho*cf*(((n3*t1-n1*t3)*(B3*n2*t1-B3*n1*t2+B1*n3*t2-B1*n2*t3)
                +B2*((n1s+n3s)*t2s-2.0*n2*t2*(n1*t1+n3*t3)
                    +n2s*(t1s+t3s)))*alphas+a*rhosq*alphaf*(-n1*n2*t1+n1s*t2+n3*TxN1));
    double Term73 = rho*cf*(a*(-n1*n3*t1+n1s*t3+n2*(-n3*t2+n2*t3))*rhosq
            *alphaf+alphas*(B3*(n3s*(t1s+t2s)-2.0*n3*(n1*t1+n2*t2)*t3+(n1s+n2s)*t3s)
                +B2*(n3*t1-n1*t3)*TxN3+B1*(-n3*t2+n2*t3)*TxN3));
    double Term74 = alphas/(sq2*a2*gu1sq*rhosq*(alphaf2+alphas2));
    double Term75 = -alphaf*(n2s*t1-n1*n2*t2+n3*(n3*t1-n1*t3));
    double Term76 = -alphaf*(-n1*n2*t1+n1s*t2+n3*TxN1);
    double Term77 = -alphaf*(-n1*n3*t1+n1s*t3+n2*(-n3*t2+n2*t3));
    double  nen71 = sq2*beta1*nen*gu1sq*(a*BN2*alphas2+rhosq*cf2*alphaf*(a*rhosq*alphaf+BT*alphas));
    double  nen72 = nen52;

    dTensor2 lu(8,8);

    // ---------------------------------------------------
    // PRIMITIVE E-VECTORS
    // ---------------------------------------------------
    // 1 - left eigenvector
    lu.set(1,1,  1.0/(sqgam*rhosq)     );
    lu.set(1,2,  0.0                   );
    lu.set(1,3,  0.0                   );
    lu.set(1,4,  0.0                   );
    lu.set(1,5, -1.0/(a2*sqgam*rhosq)  );
    lu.set(1,6,  0.0                   );
    lu.set(1,7,  0.0                   );
    lu.set(1,8,  0.0                   );

    // 2 - left eigenvector
    lu.set(2,1,  0.0                                    );
    lu.set(2,2,  0.0                                    );
    lu.set(2,3,  0.0                                    );
    lu.set(2,4,  0.0                                    );
    lu.set(2,5,  0.0                                    );
    lu.set(2,6,  (-n2*t1*t2-n3*t1*t3+n1*(t2s+t3s))/nen2 );
    lu.set(2,7,  (-t2*(n1*t1+n3*t3)+n2*(t1s+t3s))/nen2  );
    lu.set(2,8,  (n3*(t1s+t2s)-(n1*t1+n2*t2)*t3)/nen2   );

    // 3 - left eigenvector
    lu.set(3,1,  0.0            );
    lu.set(3,2, -rho*TxN1/nen3  );
    lu.set(3,3, -rho*TxN2/nen3  );
    lu.set(3,4, -rho*TxN3/nen3  );
    lu.set(3,5,  0.0            );
    lu.set(3,6,  TxN1/nen31     );
    lu.set(3,7,  TxN2/nen31     );
    lu.set(3,8,  TxN3/nen31     );

    // 4 - left eigenvector
    lu.set(4,1,  lu.get(3,1)  );
    lu.set(4,2, -lu.get(3,2)  );
    lu.set(4,3, -lu.get(3,3)  );
    lu.set(4,4, -lu.get(3,4)  );
    lu.set(4,5,  lu.get(3,5)  );
    lu.set(4,6,  lu.get(3,6)  );
    lu.set(4,7,  lu.get(3,7)  );
    lu.set(4,8,  lu.get(3,8)  );

    // 5 - left eigenvector
    lu.set(5,1,  0.0           );
    lu.set(5,2, -Term51/nen51  );
    lu.set(5,3,  Term52/nen51  );
    lu.set(5,4,  Term53/nen51  );
    lu.set(5,5,  Term54        );
    lu.set(5,6,  Term55/nen52  );
    lu.set(5,7,  Term56/nen52  );
    lu.set(5,8,  Term57/nen52  );

    // 6 - left eigenvector
    lu.set(6,1,  lu.get(5,1)  );
    lu.set(6,2, -lu.get(5,2)  );
    lu.set(6,3, -lu.get(5,3)  );
    lu.set(6,4, -lu.get(5,4)  );
    lu.set(6,5,  lu.get(5,5)  );
    lu.set(6,6,  lu.get(5,6)  );
    lu.set(6,7,  lu.get(5,7)  );
    lu.set(6,8,  lu.get(5,8)  );

    // 7 - left eigenvector
    lu.set(7,1,  0.0           );
    lu.set(7,2,  Term71/nen71  );
    lu.set(7,3,  Term72/nen71  );
    lu.set(7,4,  Term73/nen71  );
    lu.set(7,5,  Term74        );
    lu.set(7,6,  Term75/nen72  );
    lu.set(7,7,  Term76/nen72  );
    lu.set(7,8,  Term77/nen72  );

    // 8 - left eigenvector
    lu.set(8,1,  lu.get(7,1)  );
    lu.set(8,2, -lu.get(7,2)  );
    lu.set(8,3, -lu.get(7,3)  );
    lu.set(8,4, -lu.get(7,4)  );
    lu.set(8,5,  lu.get(7,5)  );
    lu.set(8,6,  lu.get(7,6)  );
    lu.set(8,7,  lu.get(7,7)  );
    lu.set(8,8,  lu.get(7,8)  );



    // ---------------------------------------------------
    // CONSERVATIVE E-VECTORS
    // ---------------------------------------------------

    dTensor2 Lmat(8, 8);
    for (int m=1; m<=8; m++)
    {
        Lmat.set(m,1,   lu.get(m,1)-lu.get(m,2)*u1/rho
                -lu.get(m,3)*u2/rho
                -lu.get(m,4)*u3/rho+lu.get(m,5)*um2*gm1     );
        Lmat.set(m,2,   lu.get(m,2)/rho-lu.get(m,5)*u1*gm1   );
        Lmat.set(m,3,   lu.get(m,3)/rho-lu.get(m,5)*u2*gm1   );
        Lmat.set(m,4,   lu.get(m,4)/rho-lu.get(m,5)*u3*gm1   );
        Lmat.set(m,5,   lu.get(m,5)*gm1                      );
        Lmat.set(m,6,  -lu.get(m,5)*B1*gm1+lu.get(m,6)       );
        Lmat.set(m,7,  -lu.get(m,5)*B2*gm1+lu.get(m,7)       );
        Lmat.set(m,8,  -lu.get(m,5)*B3*gm1+lu.get(m,8)       );
    }
    const int nsize = Qvals.getsize(2);
    // Project onto left eigenvectors
    for (int n=1; n<=nsize; n++)
        for (int m=1; m<=8; m++)
        {
            Wvals.set(m,n, 0.0 );

            for (int ell=1; ell<=8; ell++)
            {
                double tmp = Wvals.get(m,n);
                Wvals.set(m,n, tmp + Lmat.get(m,ell)*Qvals.get(ell,n) );
            }
        } 
}

