#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "ConstructHJ_L.h"

#include "StateVars.h"
#include "IniParams.h"
#include "assert.h"

#include "dog_math.h"
#include "app_defined.h"

#include <iostream>



static std::pair<double, double> max_speed_in_x_y_directions(const StateVars& Q);


// Reference: arxiv 1309.3344v2.  See equation (5.15).
// Note:  This only gives first order of accuracy in time.  To obtain
//        higher order, we must either use a high-order
//        time-integrator, such as Runge-Kutta, or add more terms to
//        this update, like what we are going to do in the other
//        Constructxxx functions in this file.
//
void ConstructHJ_L_Order1(const StateVars& Q, dTensorBC3& Lauxstar)
{
    using std::pair;
    const dTensorBC3&   q = Q.const_ref_q();
    const dTensorBC3& aux = Q.const_ref_aux();


    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Parameters for the current grid
//    const int   meqn = global_ini_params.get_meqn();
    const int   maux = global_ini_params.get_maux();
    assert(maux == 1);
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );
    
    const pair<double, double> max_speed_x_y = max_speed_in_x_y_directions(Q);
    const double max_speed_x = max_speed_x_y.first;
    const double max_speed_y = max_speed_x_y.second;
    
    assert(maux == 1);
    assert(mbc >= 2 + r);

    Lauxstar.setall(0.0);
#pragma omp parallel for
    for(int i = 1 - 2 ; i <= mx + 2 ; ++i){
         dTensor2 weno_input(maux, ws);
         dTensor2 weno_result(maux, 1);
         for(int j = 1 - 2 ; j <= my + 2; ++j){
        
            // A^{3}_{,x}^{-}
            double A3xm;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i - r + k, j, 1) - aux.get(i - r + k - 1, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xm = weno_result.get(1, 1);

            // A^{3}_{,x}^{+}
            double A3xp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i + r - (k - 1), j, 1) - aux.get(i + r - k, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xp = weno_result.get(1, 1);


            // A^{3}_{,y}^{-}
            double A3ym;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j - r + k, 1) - aux.get(i, j - r + k - 1, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3ym = weno_result.get(1, 1);

            // A^{3}_{,y}^{+}
            double A3yp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j + r - (k - 1), 1) - aux.get(i, j + r - k, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3yp = weno_result.get(1, 1);
            
            double u1 = q.get(i, j, 2) / q.get(i, j, 1);
            double u2 = q.get(i, j, 3) / q.get(i, j, 1);
            Lauxstar.set(i, j, 1,
                    -u1 * 0.5 * (A3xm + A3xp) -u2 * 0.5 * (A3ym + A3yp)
                    +max_speed_x * 0.5 * (A3xp - A3xm)
                    +max_speed_y * 0.5 * (A3yp - A3ym)   );
        }
    }
        
}


// 
//
void ConstructHJ_L_Order3(const StateVars& Q, dTensorBC3& Lauxstar, double dt)
{
    const double gamma_gas = global_ini_params.get_gamma();
    
    using std::pair;
    const dTensorBC3&   q = Q.const_ref_q();
    const dTensorBC3& aux = Q.const_ref_aux();


    // Routine for WENO reconstrution
    void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
    void (*WenoReconstruct)( const dTensor2& gin, dTensor2& diff_g ) = GetWenoReconstruct();

    // Parameters for the current grid
    const int   maux = global_ini_params.get_maux();
    assert(maux == 1);
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    const int    mbc = global_ini_params.get_mbc();
    const double dx = global_ini_params.get_dx();
    const double dy = global_ini_params.get_dy();

    // Size of the WENO stencil
    const int ws = global_ini_params.get_space_order();
    const int r = (ws + 1) / 2;
    assert_ge( mbc, r );
    
    assert(Lauxstar.getsize(1) == mx);
    assert(Lauxstar.getsize(2) == my);
    assert(Lauxstar.getsize(3) == maux);

    dTensorBC2 A3x(mx, my, mbc);
    dTensorBC2 A3y(mx, my, mbc);
    dTensorBC2 A3t(mx, my, mbc);
    dTensorBC2 A3tt(mx, my, mbc);


    const pair<double, double> max_speed_x_y = max_speed_in_x_y_directions(Q);
    const double max_speed_x = max_speed_x_y.first;
    const double max_speed_y = max_speed_x_y.second;
    
    assert(maux == 1);
    assert(mbc >= 7);
    Lauxstar.setall(0.0);
#pragma omp parallel for
    for(int i = 1 - 4 ; i <= mx + 4 ; ++i){
         dTensor2 weno_input(maux, ws);
         dTensor2 weno_result(maux, 1);
         for(int j = 1 - 4 ; j <= my + 4; ++j){
        
            // A^{3}_{,x}^{-}
            double A3xm;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i - r + k, j, 1) - aux.get(i - r + k - 1, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xm = weno_result.get(1, 1);

            // A^{3}_{,x}^{+}
            double A3xp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i + r - (k - 1), j, 1) - aux.get(i + r - k, j, 1)) / dx);
            WenoReconstruct(weno_input, weno_result);
            A3xp = weno_result.get(1, 1);


            // A^{3}_{,y}^{-}
            double A3ym;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j - r + k, 1) - aux.get(i, j - r + k - 1, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3ym = weno_result.get(1, 1);

            // A^{3}_{,y}^{+}
            double A3yp;
            for(int k = 1; k <= ws; ++k)
                weno_input.set(1, k, 
                        (aux.get(i, j + r - (k - 1), 1) - aux.get(i, j + r - k, 1)) / dy);
            WenoReconstruct(weno_input, weno_result);
            A3yp = weno_result.get(1, 1);
            
            double u1 = q.get(i, j, 2) / q.get(i, j, 1);
            double u2 = q.get(i, j, 3) / q.get(i, j, 1);
            A3t.set(i, j, 
                    -u1 * 0.5 * (A3xm + A3xp) -u2 * 0.5 * (A3ym + A3yp)
                    +max_speed_x * 0.5 * (A3xp - A3xm)
                    +max_speed_y * 0.5 * (A3yp - A3ym)   );
            A3x.set(i, j, 0.5*(A3xp + A3xm));
            A3y.set(i, j, 0.5*(A3yp + A3ym));
        }
    }

    dTensorBC2 rhot(mx, my, mbc), rhou1t(mx, my, mbc), rhou2t(mx, my, mbc), rhou3t(mx, my, mbc), 
               Et(mx, my, mbc), B1t(mx, my, mbc), B2t(mx, my, mbc), B3t(mx, my, mbc);
    dTensorBC2 xflux_rhou1tt(mx, my, mbc), yflux_rhou1tt(mx, my, mbc),
               xflux_rhou2tt(mx, my, mbc), yflux_rhou2tt(mx, my, mbc);

    dTensorBC2 u1t(mx, my, mbc), u2t(mx, my, mbc), u3t(mx, my, mbc);
#pragma omp parallel for
    for(int i = 1 - 3; i <= mx + 3; ++i){
        for(int j = 1 - 3; j <= my + 3; ++j){
            const int nonsenseInt = 0;
            dTensor2 xpts(4, nonsenseInt);
            dTensor2 q_local(4, 8);
            for(int k = 1; k <= 8; ++k){
                q_local.set(1, k, q.get(i + 1, j, k));
                q_local.set(2, k, q.get(i - 1, j, k));
                q_local.set(3, k, q.get(i, j + 1, k));
                q_local.set(4, k, q.get(i, j - 1, k));
            }
            dTensor3 flux(4, 8, 2);
            dTensor2 aux_local(nonsenseInt, nonsenseInt);
            FluxFunc(xpts, q_local, aux_local, flux);
            
            rhot.set(i, j,  -(flux.get(1, 1, 1) - flux.get(2, 1, 1)) / (2.0 * dx)
                            -(flux.get(3, 1, 2) - flux.get(4, 1, 2)) / (2.0 * dy) );
            rhou1t.set(i, j, -(flux.get(1, 2, 1) - flux.get(2, 2, 1)) / (2.0 * dx)
                             -(flux.get(3, 2, 2) - flux.get(4, 2, 2)) / (2.0 * dy) );
            rhou2t.set(i, j, -(flux.get(1, 3, 1) - flux.get(2, 3, 1)) / (2.0 * dx)
                             -(flux.get(3, 3, 2) - flux.get(4, 3, 2)) / (2.0 * dy) );
            rhou3t.set(i, j, -(flux.get(1, 4, 1) - flux.get(2, 4, 1)) / (2.0 * dx)
                             -(flux.get(3, 4, 2) - flux.get(4, 4, 2)) / (2.0 * dy) );
            Et.set(i, j, -(flux.get(1, 5, 1) - flux.get(2, 5, 1)) / (2.0 * dx)
                         -(flux.get(3, 5, 2) - flux.get(4, 5, 2)) / (2.0 * dy) );
            B1t.set(i, j, -(flux.get(1, 6, 1) - flux.get(2, 6, 1)) / (2.0 * dx)
                          -(flux.get(3, 6, 2) - flux.get(4, 6, 2)) / (2.0 * dy) );
            B2t.set(i, j, -(flux.get(1, 7, 1) - flux.get(2, 7, 1)) / (2.0 * dx)
                          -(flux.get(3, 7, 2) - flux.get(4, 7, 2)) / (2.0 * dy) );
            B3t.set(i, j, -(flux.get(1, 8, 1) - flux.get(2, 8, 1)) / (2.0 * dx)
                          -(flux.get(3, 8, 2) - flux.get(4, 8, 2)) / (2.0 * dy) );

            const double rho_ij = q.get(i, j, 1);
            const double rhou1_ij = q.get(i, j, 2);
            const double rhou2_ij = q.get(i, j, 3);
            const double rhou3_ij = q.get(i, j, 4);
            const double rhot_ij = rhot.get(i, j);
            const double rhou1t_ij = rhou1t.get(i, j);
            const double rhou2t_ij = rhou2t.get(i, j);
            const double rhou3t_ij = rhou3t.get(i, j);

            u1t.set(i, j, (rhou1t_ij * rho_ij - rhou1_ij * rhot_ij)/ (rho_ij * rho_ij));
            u2t.set(i, j, (rhou2t_ij * rho_ij - rhou2_ij * rhot_ij)/ (rho_ij * rho_ij));
            u3t.set(i, j, (rhou3t_ij * rho_ij - rhou3_ij * rhot_ij)/ (rho_ij * rho_ij));

            const double u1_ij = q.get(i, j, 2) / q.get(i, j, 1);
            const double u2_ij = q.get(i, j, 3) / q.get(i, j, 1);

            const double A3tx_ij = (A3t.get(i + 1, j) - A3t.get(i - 1, j)) / (2.0 * dx);
            const double A3ty_ij = (A3t.get(i, j + 1) - A3t.get(i, j - 1)) / (2.0 * dy);
            
            const double u1t_ij = u1t.get(i, j);
            const double u2t_ij = u2t.get(i, j);

            const double A3x_ij = A3x.get(i, j);
            const double A3y_ij = A3y.get(i, j);
            
            A3tt.set(i, j, -u1t_ij * A3x_ij - u2t_ij * A3y_ij - u1_ij * A3tx_ij - u2_ij * A3ty_ij);
            const double B1_ij = q.get(i, j, 6);
            const double B2_ij = q.get(i, j, 7);
            const double B3_ij = q.get(i, j, 8);
            const double u3_ij = q.get(i, j, 4) / q.get(i, j, 1);
            const double u3t_ij = u3t.get(i, j);
            const double Et_ij = Et.get(i, j);
            const double B1t_ij = B1t.get(i, j);
            const double B2t_ij = B2t.get(i, j);
            const double B3t_ij = B3t.get(i, j);

            xflux_rhou1tt.set(i, j, rhot_ij * u1_ij * u1_ij + rho_ij * 2.0 * u1_ij * u1t_ij
                                    + (gamma_gas - 1) * (Et_ij - 0.5 * rhot_ij * (u1_ij*u1_ij + u2_ij*u2_ij + u3_ij*u3_ij)
                                        -rho_ij * (u1_ij*u1t_ij + u2_ij*u2t_ij + u3_ij*u3t_ij)  )
                                    + (2.0 - gamma_gas) * (B1_ij*B1t_ij + B2_ij*B2t_ij + B3_ij*B3t_ij)
                                    -2.0 * B1_ij * B1t_ij  );
            yflux_rhou1tt.set(i, j, rhot_ij*u1_ij*u2_ij + rho_ij*u1t_ij*u2_ij + rho_ij*u1_ij*u2t_ij
                                    - B1t_ij*B2_ij - B1_ij*B2t_ij );

            xflux_rhou2tt.set(i, j, rhot_ij*u1_ij*u2_ij + rho_ij*u1t_ij*u2_ij + rho_ij*u1_ij*u2t_ij
                                    - B1t_ij*B2_ij - B1_ij*B2t_ij );
            yflux_rhou2tt.set(i, j, rhot_ij * u2_ij * u2_ij + rho_ij * 2.0 * u2_ij * u2t_ij
                                    + (gamma_gas - 1) * (Et_ij - 0.5 * rhot_ij * (u1_ij*u1_ij + u2_ij*u2_ij + u3_ij*u3_ij)
                                        -rho_ij * (u1_ij*u1t_ij + u2_ij*u2t_ij + u3_ij*u3t_ij)  )
                                    + (2.0 - gamma_gas) * (B1_ij*B1t_ij + B2_ij*B2t_ij + B3_ij*B3t_ij)
                                    -2.0 * B2_ij * B2t_ij  );           
        }
    }


#pragma omp parallel for
    for(int i = 1 - 2; i <= mx + 2; ++i){
        for(int j = 1 - 2; j <= my + 2; ++j){
            const double rho_ij = q.get(i, j, 1);
            const double rhou1_ij = q.get(i, j, 2);
            const double rhou2_ij = q.get(i, j, 3);
            const double u1_ij = q.get(i, j, 2) / q.get(i, j, 1);
            const double u2_ij = q.get(i, j, 3) / q.get(i, j, 1);
           
            const double A3x_ij = A3x.get(i, j);
            const double A3y_ij = A3y.get(i, j);

            const double A3tx_ij = (A3t.get(i + 1, j) - A3t.get(i - 1, j)) / (2.0 * dx);
            const double A3ty_ij = (A3t.get(i, j + 1) - A3t.get(i, j - 1)) / (2.0 * dy);
            
            const double rhot_ij = rhot.get(i, j);
            const double rhou1t_ij = rhou1t.get(i, j);
            const double rhou2t_ij = rhou2t.get(i, j);
            const double u1t_ij = u1t.get(i, j);
            const double u2t_ij = u2t.get(i, j);

            const double rhott_ij = -(rhou1t.get(i + 1, j) - rhou1t.get(i - 1, j)) / (2.0 * dx)
                                    -(rhou2t.get(i, j + 1) - rhou2t.get(i, j - 1)) / (2.0 * dy);

            const double rhou1tt_ij = -(xflux_rhou1tt.get(i + 1, j) - xflux_rhou1tt.get(i - 1, j)) / (2.0 * dx)
                                      -(yflux_rhou1tt.get(i, j + 1) - yflux_rhou1tt.get(i, j - 1)) / (2.0 * dy);
            const double rhou2tt_ij = -(xflux_rhou2tt.get(i + 1, j) - xflux_rhou2tt.get(i - 1, j)) / (2.0 * dx)
                                      -(yflux_rhou2tt.get(i, j + 1) - yflux_rhou2tt.get(i, j - 1)) / (2.0 * dy);
                                     

            // All dependencies defined
            const double u1tt_ij = (rho_ij*rho_ij * (rhou1tt_ij * rho_ij - rhou1_ij * rhott_ij) - 2*rho_ij*rhot_ij * (rhou1t_ij*rho_ij - rhou1_ij * rhot_ij)) / pow(rho_ij, 4);
            const double u2tt_ij = (rho_ij*rho_ij * (rhou2tt_ij * rho_ij - rhou2_ij * rhott_ij) - 2*rho_ij*rhot_ij * (rhou2t_ij*rho_ij - rhou2_ij * rhot_ij)) / pow(rho_ij, 4);

            const double A3ttx_ij = (A3tt.get(i + 1, j) - A3tt.get(i - 1, j)) / (2.0 * dx);
            const double A3tty_ij = (A3tt.get(i, j + 1) - A3tt.get(i, j - 1)) / (2.0 * dy);

            double A3tt_ij = A3tt.get(i, j);
            double A3ttt_ij = -u1tt_ij * A3x_ij - u2tt_ij * A3y_ij
                                    -2.0 * u1t_ij * A3tx_ij - 2.0 * u2t_ij * A3ty_ij
                                    -u1_ij * A3ttx_ij - u2_ij * A3tty_ij;

            Lauxstar.set(i, j, 1, 
                A3t.get(i, j) + 0.5*dt * A3tt_ij + 1.0/6.0*dt*dt * A3ttt_ij);
        }
    }
    
    
}

static std::pair<double, double> max_speed_in_x_y_directions(const StateVars& Q){
    using std::abs;
    using std::max;
    using std::make_pair;

    const dTensorBC3& q = Q.const_ref_q();
    const int     mx = global_ini_params.get_mx();
    const int     my = global_ini_params.get_my();
    
    double max_x_speed = 1e-15;
    double max_y_speed = 1e-15;
    for(int i = 1; i <= mx; ++i){
        for(int j = 1; j <= my; ++j){
            assert(q.get(i, j, 1) != 0);
            double x_speed = std::abs(q.get(i, j, 2) / q.get(i, j, 1));
            double y_speed = std::abs(q.get(i, j, 3) / q.get(i, j, 1));
            max_x_speed = max(max_x_speed, x_speed);
            max_y_speed = max(max_y_speed, y_speed);
        }
    }
    using std::cerr;
    using std::endl;
//    cerr << "max_x_speed: " << max_x_speed << endl;
//    cerr << "max_y_speed: " << max_y_speed << endl;
    return make_pair(max_x_speed, max_y_speed);
}




