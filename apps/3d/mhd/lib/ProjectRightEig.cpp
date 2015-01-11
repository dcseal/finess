#include "dogdefs.h"
#include "dog_math.h"
#include "IniParams.h"

// This is a user-supplied routine that projects
// Wvals onto the right eigenvectors ofthe flux 
// Jacobian; the result is stored in Qvals
//
void ProjectRightEig(int ixy, 
        const dTensor1& Aux_ave, 
        const dTensor1& Q_ave, 
        const dTensor2& Wvals,
        dTensor2& Qvals)
{    

    const int meqn = Qvals.getsize(1);

    // Average states
    double const gamma  = global_ini_params.get_gamma();
    double rho    = fabs(Q_ave.get(1));
    double u1     = Q_ave.get(2)/rho;
    double u2     = Q_ave.get(3)/rho;
    double u3     = Q_ave.get(4)/rho;
    double energy = Q_ave.get(5);
    double B1     = Q_ave.get(6);
    double B2     = Q_ave.get(7);
    double B3     = Q_ave.get(8);
    double um2    = 0.5*(u1*u1 + u2*u2 + u3*u3);
    double Bm2    = 0.5*(B1*B1 + B2*B2 + B3*B3);
    double p      = fabs((gamma-1.0)*(energy - rho*um2 - Bm2));


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

    // Some useful quantities
    double   gm1 = gamma-1.0;
    double rhosq = sqrt(rho);
    double    a2 = (gamma*p/rho);
    double     a = sqrt(a2);
    double  sqg2 = sqrt(1.0/(2.0*gamma));
    double  sq12 = sqrt(0.5);
    double  sqpr = sqrt(p)/rho;
    double sqpor = sqrt(p/rho);
    double sq1og = sqrt(1.0/gamma);
    double   b1s = B1/rhosq;
    double   b2s = B2/rhosq;
    double   b3s = B3/rhosq;
    double    bN = b1s*nv.get(1)+b2s*nv.get(2)+b3s*nv.get(3);
    double   bN2 = pow(bN,2);
    double     d = a2 + (b1s*b1s + b2s*b2s + b3s*b3s);
    double    d2 = pow(d,2);
    double    cf = sqrt(0.5*fabs(d+sqrt(d2-4.0*a2*bN2)));
    double   cf2 = pow(cf,2);
    double    cs = sqrt(0.5*fabs(d-sqrt(d2-4.0*a2*bN2)));
    double beta1 = double( (bN>=0.0) - (bN<0.0) );
    double   bst = (b1s*tv.get(1)+b2s*tv.get(2)+b3s*tv.get(3));

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

    double TxN1 = nv.get(3)*tv.get(2)-nv.get(2)*tv.get(3);
    double TxN2 = nv.get(1)*tv.get(3)-nv.get(3)*tv.get(1);
    double TxN3 = nv.get(2)*tv.get(1)-nv.get(1)*tv.get(2);

    dTensor2 ru(8,8);

    // ---------------------------------------------------
    // PRIMITIVE E-VECTORS
    // ---------------------------------------------------
    // 1 - right eigenvector Entropy Wave
    //     lamda = V.n  
    ru.set(1,1,  sqrt(gm1/gamma)*rhosq);
    ru.set(2,1,  0.0                  );
    ru.set(3,1,  0.0                  );
    ru.set(4,1,  0.0                  );
    ru.set(5,1,  0.0                  );
    ru.set(6,1,  0.0                  );
    ru.set(7,1,  0.0                  );
    ru.set(8,1,  0.0                  );

    // 2 - right eigenvector Divergence Wave
    //     lamda = V.n
    ru.set(1,2,  0.0                  );
    ru.set(2,2,  0.0                  );
    ru.set(3,2,  0.0                  );
    ru.set(4,2,  0.0                  );
    ru.set(5,2,  0.0                  );
    ru.set(6,2,  sq1og*a*nv.get(1)    );
    ru.set(7,2,  sq1og*a*nv.get(2)    );
    ru.set(8,2,  sq1og*a*nv.get(3)    );

    // 3 - right eigenvector Alfven Wave 
    //     lamda = V.n + b.n
    ru.set(1,3,   0.0                 );
    ru.set(2,3,  -sq12*sqpr*TxN1      );
    ru.set(3,3,  -sq12*sqpr*TxN2      );
    ru.set(4,3,  -sq12*sqpr*TxN3      );
    ru.set(5,3,   0.0                 );
    ru.set(6,3,   sq12*sqpor*TxN1     );
    ru.set(7,3,   sq12*sqpor*TxN2     );
    ru.set(8,3,   sq12*sqpor*TxN3     );

    // 4 - right eigenvector Alfven Wave
    //     lamda = V.n - b.n
    ru.set(1,4,   ru.get(1,3)         );
    ru.set(2,4,  -ru.get(2,3)         );
    ru.set(3,4,  -ru.get(3,3)         );
    ru.set(4,4,  -ru.get(4,3)         );
    ru.set(5,4,   ru.get(5,3)         );
    ru.set(6,4,   ru.get(6,3)         );
    ru.set(7,4,   ru.get(7,3)         );
    ru.set(8,4,   ru.get(8,3)         );

    // 5 - right eigenvector
    //     lamda = V.n + C_f
    ru.set(1,5,   sqg2*alphaf*rhosq                                                           );
    ru.set(2,5,   sqg2*(alphaf*a2*nv.get(1)+alphas*a*(bst*nv.get(1)-bN*tv.get(1)))/(rhosq*cf) );
    ru.set(3,5,   sqg2*(alphaf*a2*nv.get(2)+alphas*a*(bst*nv.get(2)-bN*tv.get(2)))/(rhosq*cf) );
    ru.set(4,5,   sqg2*(alphaf*a2*nv.get(3)+alphas*a*(bst*nv.get(3)-bN*tv.get(3)))/(rhosq*cf) );
    ru.set(5,5,   sqg2*alphaf*rhosq*a2                                                        );
    ru.set(6,5,   sqg2*alphas*a*tv.get(1)                                                     );
    ru.set(7,5,   sqg2*alphas*a*tv.get(2)                                                     );
    ru.set(8,5,   sqg2*alphas*a*tv.get(3)                                                     );

    // 6 - right eigenvector
    //     lamda = V.n - C_f
    ru.set(1,6,   ru.get(1,5)        );
    ru.set(2,6,  -ru.get(2,5)        );
    ru.set(3,6,  -ru.get(3,5)        );
    ru.set(4,6,  -ru.get(4,5)        );
    ru.set(5,6,   ru.get(5,5)        );
    ru.set(6,6,   ru.get(6,5)        );
    ru.set(7,6,   ru.get(7,5)        );
    ru.set(8,6,   ru.get(8,5)        );

    // 7 - right eigenvector
    //     lamda = V.n + C_s
    ru.set(1,7,   sqg2*alphas*rhosq );
    ru.set(2,7,   beta1*sqg2*(alphaf*cf2*tv.get(1)+a*nv.get(1)*alphas*bN)/(rhosq*cf)   );
    ru.set(3,7,   beta1*sqg2*(alphaf*cf2*tv.get(2)+alphas*a*bN*nv.get(2))/(rhosq*cf)   );
    ru.set(4,7,   beta1*sqg2*(alphaf*cf2*tv.get(3)+alphas*a*bN*nv.get(3))/(rhosq*cf)   );
    ru.set(5,7,   a2*sqg2*alphas*rhosq                                                 );
    ru.set(6,7,  -sqg2*alphaf*a*tv.get(1)                                              );
    ru.set(7,7,  -sqg2*alphaf*a*tv.get(2)                                              );
    ru.set(8,7,  -sqg2*alphaf*a*tv.get(3)                                              );

    // 8 - right eigenvector
    //     lamda = V.n - C_s
    ru.set(1,8,   ru.get(1,7)        );
    ru.set(2,8,  -ru.get(2,7)        );
    ru.set(3,8,  -ru.get(3,7)        );
    ru.set(4,8,  -ru.get(4,7)        );
    ru.set(5,8,   ru.get(5,7)        );
    ru.set(6,8,   ru.get(6,7)        );
    ru.set(7,8,   ru.get(7,7)        );
    ru.set(8,8,   ru.get(8,7)        );


    // ---------------------------------------------------
    // CONSERVATIVE E-VECTORS
    // ---------------------------------------------------    
    dTensor2 Rmat(8, 8);
    for (int m=1; m<=8; m++)
    {
        Rmat.set(1,m,  ru.get(1,m) );
        Rmat.set(2,m,  ru.get(2,m)*rho + ru.get(1,m)*u1 );
        Rmat.set(3,m,  ru.get(3,m)*rho + ru.get(1,m)*u2 );
        Rmat.set(4,m,  ru.get(4,m)*rho + ru.get(1,m)*u3 );
        Rmat.set(5,m,  um2*ru.get(1,m) + rho*u1*ru.get(2,m)
                + rho*u2*ru.get(3,m) + rho*u3*ru.get(4,m) + ru.get(5,m)/gm1
                + B1*ru.get(6,m) + B2*ru.get(7,m) + B3*ru.get(8,m) );
        Rmat.set(6,m,  ru.get(6,m) );
        Rmat.set(7,m,  ru.get(7,m) );
        Rmat.set(8,m,  ru.get(8,m) );
    }

    const int nsize = Wvals.getsize(2);
    for (int n=1; n<=nsize; n++)
        for (int m=1; m<=8; m++)
        {
            Qvals.set(m,n, 0.0 );

            for (int ell=1; ell<=8; ell++)
            {
                double tmp = Qvals.get(m,n);
                Qvals.set(m,n, tmp + Rmat.get(m,ell)*Wvals.get(ell,n) );
            }
        }


}
