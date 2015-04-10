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



    const pair<double, double> max_speed_x_y = max_speed_in_x_y_directions(Q);
    const double max_speed_x = max_speed_x_y.first;
    const double max_speed_y = max_speed_x_y.second;

    assert(maux == 1);
    assert(mbc >= 5);
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
            double a3t =  
                    -u1 * 0.5 * (A3xm + A3xp) -u2 * 0.5 * (A3ym + A3yp)
                    +max_speed_x * 0.5 * (A3xp - A3xm)
                    +max_speed_y * 0.5 * (A3yp - A3ym)   ;
            const dTensorBC3&   qq = q;
            using std::pow;
#define Power pow

            double a3X0Y0;
            double a3X0Ymm;
            double a3X0Ym;
            double a3X0Yp;
            double a3X0Ypp;
            double a3XmmY0;
            double a3XmY0;
            double a3XmYm;
            double a3XmYp;
            double a3XpY0;
            double a3XpYm;
            double a3XpYp;
            double a3XppY0;
            double B1X0Y0;
            double B1X0Ym;
            double B1X0Yp;
            double B1XmY0;
            double B1XmYm;
            double B1XmYp;
            double B1XpY0;
            double B1XpYm;
            double B1XpYp;
            double B2X0Y0;
            double B2X0Ym;
            double B2X0Yp;
            double B2XmY0;
            double B2XmYm;
            double B2XmYp;
            double B2XpY0;
            double B2XpYm;
            double B2XpYp;
            double B3X0Y0;
            double B3X0Ym;
            double B3X0Yp;
            double B3XmY0;
            double B3XmYm;
            double B3XmYp;
            double B3XpY0;
            double B3XpYm;
            double B3XpYp;
            double eX0Y0;
            double eX0Ym;
            double eX0Yp;
            double eXmY0;
            double eXmYm;
            double eXmYp;
            double eXpY0;
            double eXpYm;
            double eXpYp;
            double rhoX0Y0;
            double rhoX0Ym;
            double rhoX0Yp;
            double rhoXmY0;
            double rhoXmYm;
            double rhoXmYp;
            double rhoXpY0;
            double rhoXpYm;
            double rhoXpYp;
            double rhou1X0Y0;
            double rhou1X0Ym;
            double rhou1X0Yp;
            double rhou1XmY0;
            double rhou1XmYm;
            double rhou1XmYp;
            double rhou1XpY0;
            double rhou1XpYm;
            double rhou1XpYp;
            double rhou2X0Y0;
            double rhou2X0Ym;
            double rhou2X0Yp;
            double rhou2XmY0;
            double rhou2XmYm;
            double rhou2XmYp;
            double rhou2XpY0;
            double rhou2XpYm;
            double rhou2XpYp;
            double rhou3X0Y0;
            double rhou3X0Ym;
            double rhou3X0Yp;
            double rhou3XmY0;
            double rhou3XmYm;
            double rhou3XmYp;
            double rhou3XpY0;
            double rhou3XpYm;
            double rhou3XpYp;

            a3X0Y0 = aux.get(i, j, 1);
            a3X0Ymm = aux.get(i, -2 + j, 1);
            a3X0Ym = aux.get(i, -1 + j, 1);
            a3X0Yp = aux.get(i, 1 + j, 1);
            a3X0Ypp = aux.get(i, 2 + j, 1);
            a3XmmY0 = aux.get(-2 + i, j, 1);
            a3XmY0 = aux.get(-1 + i, j, 1);
            a3XmYm = aux.get(-1 + i, -1 + j, 1);
            a3XmYp = aux.get(-1 + i, 1 + j, 1);
            a3XpY0 = aux.get(1 + i, j, 1);
            a3XpYm = aux.get(1 + i, -1 + j, 1);
            a3XpYp = aux.get(1 + i, 1 + j, 1);
            a3XppY0 = aux.get(2 + i, j, 1);
            B1X0Y0 = qq.get(i, j, 6);
            B1X0Ym = qq.get(i, -1 + j, 6);
            B1X0Yp = qq.get(i, 1 + j, 6);
            B1XmY0 = qq.get(-1 + i, j, 6);
            B1XmYm = qq.get(-1 + i, -1 + j, 6);
            B1XmYp = qq.get(-1 + i, 1 + j, 6);
            B1XpY0 = qq.get(1 + i, j, 6);
            B1XpYm = qq.get(1 + i, -1 + j, 6);
            B1XpYp = qq.get(1 + i, 1 + j, 6);
            B2X0Y0 = qq.get(i, j, 7);
            B2X0Ym = qq.get(i, -1 + j, 7);
            B2X0Yp = qq.get(i, 1 + j, 7);
            B2XmY0 = qq.get(-1 + i, j, 7);
            B2XmYm = qq.get(-1 + i, -1 + j, 7);
            B2XmYp = qq.get(-1 + i, 1 + j, 7);
            B2XpY0 = qq.get(1 + i, j, 7);
            B2XpYm = qq.get(1 + i, -1 + j, 7);
            B2XpYp = qq.get(1 + i, 1 + j, 7);
            B3X0Y0 = qq.get(i, j, 8);
            B3X0Ym = qq.get(i, -1 + j, 8);
            B3X0Yp = qq.get(i, 1 + j, 8);
            B3XmY0 = qq.get(-1 + i, j, 8);
            B3XmYm = qq.get(-1 + i, -1 + j, 8);
            B3XmYp = qq.get(-1 + i, 1 + j, 8);
            B3XpY0 = qq.get(1 + i, j, 8);
            B3XpYm = qq.get(1 + i, -1 + j, 8);
            B3XpYp = qq.get(1 + i, 1 + j, 8);
            eX0Y0 = qq.get(i, j, 5);
            eX0Ym = qq.get(i, -1 + j, 5);
            eX0Yp = qq.get(i, 1 + j, 5);
            eXmY0 = qq.get(-1 + i, j, 5);
            eXmYm = qq.get(-1 + i, -1 + j, 5);
            eXmYp = qq.get(-1 + i, 1 + j, 5);
            eXpY0 = qq.get(1 + i, j, 5);
            eXpYm = qq.get(1 + i, -1 + j, 5);
            eXpYp = qq.get(1 + i, 1 + j, 5);
            rhoX0Y0 = qq.get(i, j, 1);
            rhoX0Ym = qq.get(i, -1 + j, 1);
            rhoX0Yp = qq.get(i, 1 + j, 1);
            rhoXmY0 = qq.get(-1 + i, j, 1);
            rhoXmYm = qq.get(-1 + i, -1 + j, 1);
            rhoXmYp = qq.get(-1 + i, 1 + j, 1);
            rhoXpY0 = qq.get(1 + i, j, 1);
            rhoXpYm = qq.get(1 + i, -1 + j, 1);
            rhoXpYp = qq.get(1 + i, 1 + j, 1);
            rhou1X0Y0 = qq.get(i, j, 2);
            rhou1X0Ym = qq.get(i, -1 + j, 2);
            rhou1X0Yp = qq.get(i, 1 + j, 2);
            rhou1XmY0 = qq.get(-1 + i, j, 2);
            rhou1XmYm = qq.get(-1 + i, -1 + j, 2);
            rhou1XmYp = qq.get(-1 + i, 1 + j, 2);
            rhou1XpY0 = qq.get(1 + i, j, 2);
            rhou1XpYm = qq.get(1 + i, -1 + j, 2);
            rhou1XpYp = qq.get(1 + i, 1 + j, 2);
            rhou2X0Y0 = qq.get(i, j, 3);
            rhou2X0Ym = qq.get(i, -1 + j, 3);
            rhou2X0Yp = qq.get(i, 1 + j, 3);
            rhou2XmY0 = qq.get(-1 + i, j, 3);
            rhou2XmYm = qq.get(-1 + i, -1 + j, 3);
            rhou2XmYp = qq.get(-1 + i, 1 + j, 3);
            rhou2XpY0 = qq.get(1 + i, j, 3);
            rhou2XpYm = qq.get(1 + i, -1 + j, 3);
            rhou2XpYp = qq.get(1 + i, 1 + j, 3);
            rhou3X0Y0 = qq.get(i, j, 4);
            rhou3X0Ym = qq.get(i, -1 + j, 4);
            rhou3X0Yp = qq.get(i, 1 + j, 4);
            rhou3XmY0 = qq.get(-1 + i, j, 4);
            rhou3XmYm = qq.get(-1 + i, -1 + j, 4);
            rhou3XmYp = qq.get(-1 + i, 1 + j, 4);
            rhou3XpY0 = qq.get(1 + i, j, 4);
            rhou3XpYm = qq.get(1 + i, -1 + j, 4);
            rhou3XpYp = qq.get(1 + i, 1 + j, 4);

            double a3t0x0y1z0Val =  (-a3X0Ym + a3X0Yp)/(2.*dy); //0.5 * (A3yp + A3ym);
            double B1t0x0y1z0Val = (-B1X0Ym + B1X0Yp)/(2.*dy);
            double B2t0x0y1z0Val = (-B2X0Ym + B2X0Yp)/(2.*dy);
            double B3t0x0y1z0Val = (-B3X0Ym + B3X0Yp)/(2.*dy);
            double et0x0y1z0Val = (-eX0Ym + eX0Yp)/(2.*dy);
            double rhot0x0y1z0Val = (-rhoX0Ym + rhoX0Yp)/(2.*dy);
            double rhou1t0x0y1z0Val = (-rhou1X0Ym + rhou1X0Yp)/(2.*dy);
            double rhou2t0x0y1z0Val = (-rhou2X0Ym + rhou2X0Yp)/(2.*dy);
            double rhou3t0x0y1z0Val = (-rhou3X0Ym + rhou3X0Yp)/(2.*dy);
            double a3t0x0y2z0Val = (-2*a3X0Y0 + a3X0Ym + a3X0Yp)/Power(dy,2);
            double a3t0x1y0z0Val = (-a3XmY0 + a3XpY0)/(2.*dx);  //0.5 * (A3xp + A3xm);
            double B1t0x1y0z0Val = (-B1XmY0 + B1XpY0)/(2.*dx);
            double B2t0x1y0z0Val = (-B2XmY0 + B2XpY0)/(2.*dx);
            double B3t0x1y0z0Val = (-B3XmY0 + B3XpY0)/(2.*dx);
            double et0x1y0z0Val = (-eXmY0 + eXpY0)/(2.*dx);
            double rhot0x1y0z0Val = (-rhoXmY0 + rhoXpY0)/(2.*dx);
            double rhou1t0x1y0z0Val = (-rhou1XmY0 + rhou1XpY0)/(2.*dx);
            double rhou2t0x1y0z0Val = (-rhou2XmY0 + rhou2XpY0)/(2.*dx);
            double rhou3t0x1y0z0Val = (-rhou3XmY0 + rhou3XpY0)/(2.*dx);
            double a3t0x1y1z0Val = (a3XmYm - a3XmYp - a3XpYm + a3XpYp)/(4.*dx*dy);
            double a3t0x2y0z0Val = (-2*a3X0Y0 + a3XmY0 + a3XpY0)/Power(dx,2);
            //double a3t0x0y1z0Val = (-a3X0Ym + a3X0Yp)/(2.*dy);
            //double B1t0x0y1z0Val = (-B1X0Ym + B1X0Yp)/(2.*dy);
            //double B2t0x0y1z0Val = (-B2X0Ym + B2X0Yp)/(2.*dy);
            //double B3t0x0y1z0Val = (-B3X0Ym + B3X0Yp)/(2.*dy);
            //double et0x0y1z0Val = (-eX0Ym + eX0Yp)/(2.*dy);
            //double rhot0x0y1z0Val = (-rhoX0Ym + rhoX0Yp)/(2.*dy);
            //double rhou1t0x0y1z0Val = (-rhou1X0Ym + rhou1X0Yp)/(2.*dy);
            //double rhou2t0x0y1z0Val = (-rhou2X0Ym + rhou2X0Yp)/(2.*dy);
            //double rhou3t0x0y1z0Val = (-rhou3X0Ym + rhou3X0Yp)/(2.*dy);
            //double a3t0x0y2z0Val = (-2*a3X0Y0 + a3X0Ym + a3X0Yp)/Power(dy,2);
            double B1t0x0y2z0Val = (-2*B1X0Y0 + B1X0Ym + B1X0Yp)/Power(dy,2);
            double B2t0x0y2z0Val = (-2*B2X0Y0 + B2X0Ym + B2X0Yp)/Power(dy,2);
            double B3t0x0y2z0Val = (-2*B3X0Y0 + B3X0Ym + B3X0Yp)/Power(dy,2);
            double et0x0y2z0Val = (-2*eX0Y0 + eX0Ym + eX0Yp)/Power(dy,2);
            double rhot0x0y2z0Val = (-2*rhoX0Y0 + rhoX0Ym + rhoX0Yp)/Power(dy,2);
            double rhou1t0x0y2z0Val = (-2*rhou1X0Y0 + rhou1X0Ym + rhou1X0Yp)/Power(dy,2);
            double rhou2t0x0y2z0Val = (-2*rhou2X0Y0 + rhou2X0Ym + rhou2X0Yp)/Power(dy,2);
            double rhou3t0x0y2z0Val = (-2*rhou3X0Y0 + rhou3X0Ym + rhou3X0Yp)/Power(dy,2);
            double a3t0x0y3z0Val = (2*a3X0Ym - a3X0Ymm - 2*a3X0Yp + a3X0Ypp)/(2.*Power(dy,3));
            //double a3t0x1y0z0Val = (-a3XmY0 + a3XpY0)/(2.*dx);
            //double B1t0x1y0z0Val = (-B1XmY0 + B1XpY0)/(2.*dx);
            //double B2t0x1y0z0Val = (-B2XmY0 + B2XpY0)/(2.*dx);
            //double B3t0x1y0z0Val = (-B3XmY0 + B3XpY0)/(2.*dx);
            //double et0x1y0z0Val = (-eXmY0 + eXpY0)/(2.*dx);
            //double rhot0x1y0z0Val = (-rhoXmY0 + rhoXpY0)/(2.*dx);
            //double rhou1t0x1y0z0Val = (-rhou1XmY0 + rhou1XpY0)/(2.*dx);
            //double rhou2t0x1y0z0Val = (-rhou2XmY0 + rhou2XpY0)/(2.*dx);
            //double rhou3t0x1y0z0Val = (-rhou3XmY0 + rhou3XpY0)/(2.*dx);
            //double a3t0x1y1z0Val = (a3XmYm - a3XmYp - a3XpYm + a3XpYp)/(4.*dx*dy);
            double B1t0x1y1z0Val = (B1XmYm - B1XmYp - B1XpYm + B1XpYp)/(4.*dx*dy);
            double B2t0x1y1z0Val = (B2XmYm - B2XmYp - B2XpYm + B2XpYp)/(4.*dx*dy);
            double B3t0x1y1z0Val = (B3XmYm - B3XmYp - B3XpYm + B3XpYp)/(4.*dx*dy);
            double et0x1y1z0Val = (eXmYm - eXmYp - eXpYm + eXpYp)/(4.*dx*dy);
            double rhot0x1y1z0Val = (rhoXmYm - rhoXmYp - rhoXpYm + rhoXpYp)/(4.*dx*dy);
            double rhou1t0x1y1z0Val = (rhou1XmYm - rhou1XmYp - rhou1XpYm + rhou1XpYp)/(4.*dx*dy);
            double rhou2t0x1y1z0Val = (rhou2XmYm - rhou2XmYp - rhou2XpYm + rhou2XpYp)/(4.*dx*dy);
            double rhou3t0x1y1z0Val = (rhou3XmYm - rhou3XmYp - rhou3XpYm + rhou3XpYp)/(4.*dx*dy);
            double a3t0x1y2z0Val = (2*a3XmY0 - a3XmYm - a3XmYp - 2*a3XpY0 + a3XpYm + a3XpYp)/(2.*dx*Power(dy,2));
            //double a3t0x2y0z0Val = (-2*a3X0Y0 + a3XmY0 + a3XpY0)/Power(dx,2);
            double B1t0x2y0z0Val = (-2*B1X0Y0 + B1XmY0 + B1XpY0)/Power(dx,2);
            double B2t0x2y0z0Val = (-2*B2X0Y0 + B2XmY0 + B2XpY0)/Power(dx,2);
            double B3t0x2y0z0Val = (-2*B3X0Y0 + B3XmY0 + B3XpY0)/Power(dx,2);
            double et0x2y0z0Val = (-2*eX0Y0 + eXmY0 + eXpY0)/Power(dx,2);
            double rhot0x2y0z0Val = (-2*rhoX0Y0 + rhoXmY0 + rhoXpY0)/Power(dx,2);
            double rhou1t0x2y0z0Val = (-2*rhou1X0Y0 + rhou1XmY0 + rhou1XpY0)/Power(dx,2);
            double rhou2t0x2y0z0Val = (-2*rhou2X0Y0 + rhou2XmY0 + rhou2XpY0)/Power(dx,2);
            double rhou3t0x2y0z0Val = (-2*rhou3X0Y0 + rhou3XmY0 + rhou3XpY0)/Power(dx,2);
            double a3t0x2y1z0Val = (2*a3X0Ym - 2*a3X0Yp - a3XmYm + a3XmYp - a3XpYm + a3XpYp)/(2.*Power(dx,2)*dy);
            double a3t0x3y0z0Val = (-a3XmmY0 + 2*a3XmY0 + a3XppY0 - 2*a3XpY0)/(2.*Power(dx,3));



            double rhoVal = rhoX0Y0;
            double rhou1Val = rhou1X0Y0;
            double rhou2Val = rhou2X0Y0;
            double rhou3Val = rhou3X0Y0;
            double B1Val = B1X0Y0;
            double B2Val = B2X0Y0;
            double B3Val = B3X0Y0;
            double eVal = eX0Y0;
#define gamma gamma_gas
            double a3tt = (-(a3t0x0y1z0Val*rhot0x0y1z0Val*Power(rhou1Val,2)) + a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*Power(rhou1Val,2) - 5*a3t0x1y0z0Val*rhot0x1y0z0Val*Power(rhou1Val,2) + a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*Power(rhou1Val,2) - 4*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1Val*rhou2Val - 4*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1Val*rhou2Val - 5*a3t0x0y1z0Val*rhot0x0y1z0Val*Power(rhou2Val,2) + a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*Power(rhou2Val,2) - a3t0x1y0z0Val*rhot0x1y0z0Val*Power(rhou2Val,2) + a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*Power(rhou2Val,2) - a3t0x0y1z0Val*rhot0x0y1z0Val*Power(rhou3Val,2) + a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*Power(rhou3Val,2) - a3t0x1y0z0Val*rhot0x1y0z0Val*Power(rhou3Val,2) + a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*Power(rhou3Val,2) + 2*a3t0x0y1z0Val*rhou1t0x0y1z0Val*rhou1Val*rhoVal - 2*a3t0x0y1z0Val*gamma*rhou1t0x0y1z0Val*rhou1Val*rhoVal + 6*a3t0x1y0z0Val*rhou1t0x1y0z0Val*rhou1Val*rhoVal - 2*a3t0x1y0z0Val*gamma*rhou1t0x1y0z0Val*rhou1Val*rhoVal + 2*a3t0x2y0z0Val*Power(rhou1Val,2)*rhoVal + 4*a3t0x0y1z0Val*rhou1Val*rhou2t0x1y0z0Val*rhoVal + 4*a3t0x1y0z0Val*rhou1t0x0y1z0Val*rhou2Val*rhoVal + 4*a3t0x1y1z0Val*rhou1Val*rhou2Val*rhoVal + 6*a3t0x0y1z0Val*rhou2t0x0y1z0Val*rhou2Val*rhoVal - 2*a3t0x0y1z0Val*gamma*rhou2t0x0y1z0Val*rhou2Val*rhoVal + 2*a3t0x1y0z0Val*rhou2t0x1y0z0Val*rhou2Val*rhoVal - 2*a3t0x1y0z0Val*gamma*rhou2t0x1y0z0Val*rhou2Val*rhoVal + 2*a3t0x0y2z0Val*Power(rhou2Val,2)*rhoVal + 2*a3t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val*rhoVal - 2*a3t0x0y1z0Val*gamma*rhou3t0x0y1z0Val*rhou3Val*rhoVal + 2*a3t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val*rhoVal - 2*a3t0x1y0z0Val*gamma*rhou3t0x1y0z0Val*rhou3Val*rhoVal + 4*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*Power(rhoVal,2) + 4*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*Power(rhoVal,2) - 2*a3t0x0y1z0Val*et0x0y1z0Val*Power(rhoVal,2) - 2*a3t0x1y0z0Val*et0x1y0z0Val*Power(rhoVal,2) - 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*gamma*Power(rhoVal,2) - 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*gamma*Power(rhoVal,2) + 2*a3t0x0y1z0Val*et0x0y1z0Val*gamma*Power(rhoVal,2) + 2*a3t0x1y0z0Val*et0x1y0z0Val*gamma*Power(rhoVal,2) - 2*B1Val*(a3t0x0y1z0Val*(B2t0x1y0z0Val + B1t0x0y1z0Val*(-2 + gamma)) + a3t0x1y0z0Val*(B2t0x0y1z0Val + B1t0x1y0z0Val*gamma))*Power(rhoVal,2) - 2*B2Val*(a3t0x1y0z0Val*(B1t0x0y1z0Val + B2t0x1y0z0Val*(-2 + gamma)) + a3t0x0y1z0Val*(B1t0x1y0z0Val + B2t0x0y1z0Val*gamma))*Power(rhoVal,2))/(2.*Power(rhoVal,3));

            double a3ttt = (3*((a3t0x1y0z0Val*(-1 + gamma)*Power(rhot0x0y1z0Val,2) + a3t0x0y1z0Val*(-3 + 2*gamma + Power(gamma,2))*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x1y0z0Val*(-12 + 3*gamma + Power(gamma,2))*Power(rhot0x1y0z0Val,2))*Power(rhou1Val,3) + (a3t0x1y0z0Val*(-19 + 2*gamma + Power(gamma,2))*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x0y1z0Val*((-4 + 3*gamma + Power(gamma,2))*Power(rhot0x0y1z0Val,2) + (-9 + gamma)*Power(rhot0x1y0z0Val,2)))*Power(rhou1Val,2)*rhou2Val + rhou1Val*((a3t0x1y0z0Val*(-9 + gamma)*Power(rhot0x0y1z0Val,2) + a3t0x0y1z0Val*(-19 + 2*gamma + Power(gamma,2))*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x1y0z0Val*(-4 + 3*gamma + Power(gamma,2))*Power(rhot0x1y0z0Val,2))*Power(rhou2Val,2) + (-1 + gamma)*(a3t0x1y0z0Val*Power(rhot0x0y1z0Val,2) + a3t0x0y1z0Val*(3 + gamma)*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x1y0z0Val*(4 + gamma)*Power(rhot0x1y0z0Val,2))*Power(rhou3Val,2)) + rhou2Val*((a3t0x1y0z0Val*(-3 + 2*gamma + Power(gamma,2))*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x0y1z0Val*((-12 + 3*gamma + Power(gamma,2))*Power(rhot0x0y1z0Val,2) + (-1 + gamma)*Power(rhot0x1y0z0Val,2)))*Power(rhou2Val,2) + (-1 + gamma)*(a3t0x1y0z0Val*(3 + gamma)*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x0y1z0Val*((4 + gamma)*Power(rhot0x0y1z0Val,2) + Power(rhot0x1y0z0Val,2)))*Power(rhou3Val,2))) + (-4*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou1Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou1Val + 4*a3t0x0y1z0Val*eVal*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou1Val + 2*a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou1Val - 4*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou1Val - 4*a3t0x1y0z0Val*Power(B3Val,2)*Power(rhot0x1y0z0Val,2)*rhou1Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*gamma*Power(rhot0x1y0z0Val,2)*rhou1Val + 4*a3t0x1y0z0Val*eVal*gamma*Power(rhot0x1y0z0Val,2)*rhou1Val + 2*a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*Power(rhot0x1y0z0Val,2)*rhou1Val - 4*a3t0x1y0z0Val*eVal*Power(gamma,2)*Power(rhot0x1y0z0Val,2)*rhou1Val + 9*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) - 9*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) + 11*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) - 8*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) - 3*a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) + 7*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) - 5*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) + 63*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) - 22*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) - 5*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) + 3*a3t0x1y1z0Val*rhot0x0y1z0Val*Power(rhou1Val,3) - 3*a3t0x1y1z0Val*gamma*rhot0x0y1z0Val*Power(rhou1Val,3) + 15*a3t0x2y0z0Val*rhot0x1y0z0Val*Power(rhou1Val,3) - 3*a3t0x2y0z0Val*gamma*rhot0x1y0z0Val*Power(rhou1Val,3) + 3*a3t0x0y1z0Val*rhot0x1y1z0Val*Power(rhou1Val,3) - 2*a3t0x0y1z0Val*gamma*rhot0x1y1z0Val*Power(rhou1Val,3) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y1z0Val*Power(rhou1Val,3) + 9*a3t0x1y0z0Val*rhot0x2y0z0Val*Power(rhou1Val,3) - 2*a3t0x1y0z0Val*gamma*rhot0x2y0z0Val*Power(rhou1Val,3) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x2y0z0Val*Power(rhou1Val,3) + 3*a3t0x0y1z0Val*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2t0x0y1z0Val - 3*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2t0x0y1z0Val - a3t0x1y0z0Val*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2t0x0y1z0Val + 3*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2t0x0y1z0Val + 13*a3t0x1y0z0Val*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2t0x1y0z0Val - a3t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2t0x1y0z0Val + 27*a3t0x0y1z0Val*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2t0x1y0z0Val - 3*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2t0x1y0z0Val - 4*a3t0x0y1z0Val*Power(B3Val,2)*Power(rhot0x0y1z0Val,2)*rhou2Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*gamma*Power(rhot0x0y1z0Val,2)*rhou2Val + 4*a3t0x0y1z0Val*eVal*gamma*Power(rhot0x0y1z0Val,2)*rhou2Val + 2*a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*Power(rhot0x0y1z0Val,2)*rhou2Val - 4*a3t0x0y1z0Val*eVal*Power(gamma,2)*Power(rhot0x0y1z0Val,2)*rhou2Val - 4*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou2Val + 4*a3t0x1y0z0Val*eVal*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou2Val + 2*a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou2Val - 4*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x0y1z0Val*rhot0x1y0z0Val*rhou2Val + 24*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val - 22*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val + 44*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val - 8*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val + 34*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val - 8*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val + 18*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val - 6*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val + 3*a3t0x0y2z0Val*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2Val + 12*a3t0x2y0z0Val*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2Val - 3*a3t0x0y2z0Val*gamma*rhot0x0y1z0Val*Power(rhou1Val,2)*rhou2Val + 3*a3t0x0y1z0Val*rhot0x0y2z0Val*Power(rhou1Val,2)*rhou2Val - 2*a3t0x0y1z0Val*gamma*rhot0x0y2z0Val*Power(rhou1Val,2)*rhou2Val - a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y2z0Val*Power(rhou1Val,2)*rhou2Val + 27*a3t0x1y1z0Val*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2Val - 3*a3t0x1y1z0Val*gamma*rhot0x1y0z0Val*Power(rhou1Val,2)*rhou2Val + 15*a3t0x1y0z0Val*rhot0x1y1z0Val*Power(rhou1Val,2)*rhou2Val - 2*a3t0x1y0z0Val*gamma*rhot0x1y1z0Val*Power(rhou1Val,2)*rhou2Val - a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y1z0Val*Power(rhou1Val,2)*rhou2Val + 6*a3t0x0y1z0Val*rhot0x2y0z0Val*Power(rhou1Val,2)*rhou2Val + 18*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val - 6*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val + 34*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val - 8*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val + 44*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val - 8*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val + 24*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val - 22*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val + 27*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou2Val,2) - 3*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou2Val,2) + 13*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou2Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou2Val,2) - a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou2Val,2) + 3*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou2Val,2) - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou2Val,2) + 3*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou2Val,2) - 3*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou2Val,2) + 27*a3t0x1y1z0Val*rhot0x0y1z0Val*rhou1Val*Power(rhou2Val,2) - 3*a3t0x1y1z0Val*gamma*rhot0x0y1z0Val*rhou1Val*Power(rhou2Val,2) + 6*a3t0x1y0z0Val*rhot0x0y2z0Val*rhou1Val*Power(rhou2Val,2) + 12*a3t0x0y2z0Val*rhot0x1y0z0Val*rhou1Val*Power(rhou2Val,2) + 3*a3t0x2y0z0Val*rhot0x1y0z0Val*rhou1Val*Power(rhou2Val,2) - 3*a3t0x2y0z0Val*gamma*rhot0x1y0z0Val*rhou1Val*Power(rhou2Val,2) + 15*a3t0x0y1z0Val*rhot0x1y1z0Val*rhou1Val*Power(rhou2Val,2) - 2*a3t0x0y1z0Val*gamma*rhot0x1y1z0Val*rhou1Val*Power(rhou2Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y1z0Val*rhou1Val*Power(rhou2Val,2) + 3*a3t0x1y0z0Val*rhot0x2y0z0Val*rhou1Val*Power(rhou2Val,2) - 2*a3t0x1y0z0Val*gamma*rhot0x2y0z0Val*rhou1Val*Power(rhou2Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x2y0z0Val*rhou1Val*Power(rhou2Val,2) + 63*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) - 22*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) - 5*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) + 7*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) - 5*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) + 11*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) - 8*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) - 3*a3t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) + 9*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) - 9*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) + 15*a3t0x0y2z0Val*rhot0x0y1z0Val*Power(rhou2Val,3) - 3*a3t0x0y2z0Val*gamma*rhot0x0y1z0Val*Power(rhou2Val,3) + 9*a3t0x0y1z0Val*rhot0x0y2z0Val*Power(rhou2Val,3) - 2*a3t0x0y1z0Val*gamma*rhot0x0y2z0Val*Power(rhou2Val,3) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y2z0Val*Power(rhou2Val,3) + 3*a3t0x1y1z0Val*rhot0x1y0z0Val*Power(rhou2Val,3) - 3*a3t0x1y1z0Val*gamma*rhot0x1y0z0Val*Power(rhou2Val,3) + 3*a3t0x1y0z0Val*rhot0x1y1z0Val*Power(rhou2Val,3) - 2*a3t0x1y0z0Val*gamma*rhot0x1y1z0Val*Power(rhou2Val,3) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y1z0Val*Power(rhou2Val,3) + 2*Power(B2Val,2)*(-((2*a3t0x1y0z0Val*Power(rhot0x0y1z0Val,2) - a3t0x0y1z0Val*(-1 + gamma)*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x1y0z0Val*(2 + gamma - Power(gamma,2))*Power(rhot0x1y0z0Val,2))*rhou1Val) + (-1 + gamma)*gamma*rhot0x0y1z0Val*(a3t0x0y1z0Val*rhot0x0y1z0Val + a3t0x1y0z0Val*rhot0x1y0z0Val)*rhou2Val) + 2*Power(B1Val,2)*((-1 + gamma)*gamma*rhot0x1y0z0Val*(a3t0x0y1z0Val*rhot0x0y1z0Val + a3t0x1y0z0Val*rhot0x1y0z0Val)*rhou1Val + (a3t0x1y0z0Val*(-1 + gamma)*gamma*rhot0x0y1z0Val*rhot0x1y0z0Val + a3t0x0y1z0Val*((-2 - gamma + Power(gamma,2))*Power(rhot0x0y1z0Val,2) - 2*Power(rhot0x1y0z0Val,2)))*rhou2Val) + 4*B2Val*B3Val*rhot0x0y1z0Val*(a3t0x0y1z0Val*rhot0x0y1z0Val + a3t0x1y0z0Val*rhot0x1y0z0Val)*rhou3Val + 6*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val - 6*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val + 10*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val - 8*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val + 24*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val - 22*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val + 8*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val - 8*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val + 8*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val - 8*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val + 24*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val - 22*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val + 10*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val - 8*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val - 6*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val + 3*a3t0x1y0z0Val*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou3Val,2) - 3*a3t0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1t0x0y1z0Val*Power(rhou3Val,2) + a3t0x0y1z0Val*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou3Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val*Power(rhou3Val,2) - a3t0x0y1z0Val*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou3Val,2) + 3*a3t0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou3Val,2) - 2*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val*Power(rhou3Val,2) + 3*a3t0x1y0z0Val*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou3Val,2) - 3*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val*Power(rhou3Val,2) + 3*a3t0x1y1z0Val*rhot0x0y1z0Val*rhou1Val*Power(rhou3Val,2) - 3*a3t0x1y1z0Val*gamma*rhot0x0y1z0Val*rhou1Val*Power(rhou3Val,2) + 3*a3t0x2y0z0Val*rhot0x1y0z0Val*rhou1Val*Power(rhou3Val,2) - 3*a3t0x2y0z0Val*gamma*rhot0x1y0z0Val*rhou1Val*Power(rhou3Val,2) + 3*a3t0x0y1z0Val*rhot0x1y1z0Val*rhou1Val*Power(rhou3Val,2) - 2*a3t0x0y1z0Val*gamma*rhot0x1y1z0Val*rhou1Val*Power(rhou3Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x1y1z0Val*rhou1Val*Power(rhou3Val,2) + 3*a3t0x1y0z0Val*rhot0x2y0z0Val*rhou1Val*Power(rhou3Val,2) - 2*a3t0x1y0z0Val*gamma*rhot0x2y0z0Val*rhou1Val*Power(rhou3Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x2y0z0Val*rhou1Val*Power(rhou3Val,2) + 3*a3t0x0y1z0Val*rhot0x0y1z0Val*rhou2t0x0y1z0Val*Power(rhou3Val,2) - 3*a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x0y1z0Val*Power(rhou3Val,2) - a3t0x1y0z0Val*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou3Val,2) + 3*a3t0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou3Val,2) - 2*a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val*Power(rhou3Val,2) + a3t0x1y0z0Val*rhot0x0y1z0Val*rhou2t0x1y0z0Val*Power(rhou3Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x1y0z0Val*Power(rhou3Val,2) + 3*a3t0x0y1z0Val*rhot0x1y0z0Val*rhou2t0x1y0z0Val*Power(rhou3Val,2) - 3*a3t0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou2t0x1y0z0Val*Power(rhou3Val,2) + 3*a3t0x0y2z0Val*rhot0x0y1z0Val*rhou2Val*Power(rhou3Val,2) - 3*a3t0x0y2z0Val*gamma*rhot0x0y1z0Val*rhou2Val*Power(rhou3Val,2) + 3*a3t0x0y1z0Val*rhot0x0y2z0Val*rhou2Val*Power(rhou3Val,2) - 2*a3t0x0y1z0Val*gamma*rhot0x0y2z0Val*rhou2Val*Power(rhou3Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhot0x0y2z0Val*rhou2Val*Power(rhou3Val,2) + 3*a3t0x1y1z0Val*rhot0x1y0z0Val*rhou2Val*Power(rhou3Val,2) - 3*a3t0x1y1z0Val*gamma*rhot0x1y0z0Val*rhou2Val*Power(rhou3Val,2) + 3*a3t0x1y0z0Val*rhot0x1y1z0Val*rhou2Val*Power(rhou3Val,2) - 2*a3t0x1y0z0Val*gamma*rhot0x1y1z0Val*rhou2Val*Power(rhou3Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhot0x1y1z0Val*rhou2Val*Power(rhou3Val,2) + 4*B1Val*(B2Val*(Power(rhot0x0y1z0Val,2) + Power(rhot0x1y0z0Val,2))*(a3t0x0y1z0Val*rhou1Val + a3t0x1y0z0Val*rhou2Val) + B3Val*rhot0x1y0z0Val*(a3t0x0y1z0Val*rhot0x0y1z0Val + a3t0x1y0z0Val*rhot0x1y0z0Val)*rhou3Val))*rhoVal - (-2*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x1y0z0Val*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*eVal*gamma*rhot0x1y0z0Val*rhou1t0x0y1z0Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*eVal*gamma*rhot0x0y1z0Val*rhou1t0x1y0z0Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x0y1z0Val*rhou1t0x1y0z0Val - 4*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x1y0z0Val*rhou1t0x1y0z0Val + 4*a3t0x1y0z0Val*eVal*gamma*rhot0x1y0z0Val*rhou1t0x1y0z0Val + 2*a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val - 4*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val - 12*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*rhot0x0y1z0Val*rhou1Val - 4*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*rhot0x0y1z0Val*rhou1Val + 6*a3t0x1y0z0Val*et0x0y1z0Val*rhot0x0y1z0Val*rhou1Val + 2*a3t0x0y1z0Val*et0x1y0z0Val*rhot0x0y1z0Val*rhou1Val + 6*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou1Val + 2*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou1Val - 6*a3t0x1y0z0Val*et0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou1Val - 2*a3t0x0y1z0Val*et0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou1Val - 12*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*rhot0x1y0z0Val*rhou1Val - 28*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*rhot0x1y0z0Val*rhou1Val + 4*a3t0x0y1z0Val*et0x0y1z0Val*rhot0x1y0z0Val*rhou1Val + 12*a3t0x1y0z0Val*et0x1y0z0Val*rhot0x1y0z0Val*rhou1Val + 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou1Val + 10*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou1Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou1Val - 10*a3t0x1y0z0Val*et0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou1Val + 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val + 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*Power(gamma,2)*rhot0x1y0z0Val*rhou1Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x1y1z0Val*rhou1Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x1y1z0Val*rhou1Val + 2*a3t0x0y1z0Val*eVal*gamma*rhot0x1y1z0Val*rhou1Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x1y1z0Val*rhou1Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x1y1z0Val*rhou1Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x2y0z0Val*rhou1Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x2y0z0Val*rhou1Val + 2*a3t0x1y0z0Val*eVal*gamma*rhot0x2y0z0Val*rhou1Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x2y0z0Val*rhou1Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x2y0z0Val*rhou1Val + 6*a3t0x1y0z0Val*Power(rhou1t0x0y1z0Val,2)*rhou1Val - 6*a3t0x1y0z0Val*gamma*Power(rhou1t0x0y1z0Val,2)*rhou1Val + 6*a3t0x0y1z0Val*rhou1t0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val - 4*a3t0x0y1z0Val*gamma*rhou1t0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val + 24*a3t0x1y0z0Val*Power(rhou1t0x1y0z0Val,2)*rhou1Val - 10*a3t0x1y0z0Val*gamma*Power(rhou1t0x1y0z0Val,2)*rhou1Val - 2*a3t0x1y0z0Val*Power(gamma,2)*Power(rhou1t0x1y0z0Val,2)*rhou1Val + 6*a3t0x1y1z0Val*rhou1t0x0y1z0Val*Power(rhou1Val,2) - 6*a3t0x1y1z0Val*gamma*rhou1t0x0y1z0Val*Power(rhou1Val,2) + 18*a3t0x2y0z0Val*rhou1t0x1y0z0Val*Power(rhou1Val,2) - 6*a3t0x2y0z0Val*gamma*rhou1t0x1y0z0Val*Power(rhou1Val,2) + 6*a3t0x0y1z0Val*rhou1t0x1y1z0Val*Power(rhou1Val,2) - 5*a3t0x0y1z0Val*gamma*rhou1t0x1y1z0Val*Power(rhou1Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x1y1z0Val*Power(rhou1Val,2) + 12*a3t0x1y0z0Val*rhou1t0x2y0z0Val*Power(rhou1Val,2) - 5*a3t0x1y0z0Val*gamma*rhou1t0x2y0z0Val*Power(rhou1Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x2y0z0Val*Power(rhou1Val,2) + 2*a3t0x3y0z0Val*Power(rhou1Val,3) - 4*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x0y1z0Val*rhou2t0x0y1z0Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x0y1z0Val*rhou2t0x0y1z0Val + 4*a3t0x0y1z0Val*eVal*gamma*rhot0x0y1z0Val*rhou2t0x0y1z0Val + 2*a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x0y1z0Val - 4*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x1y0z0Val*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*eVal*gamma*rhot0x1y0z0Val*rhou2t0x0y1z0Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val + 6*a3t0x0y1z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val - 4*a3t0x0y1z0Val*gamma*rhou1t0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val + 4*a3t0x1y0z0Val*gamma*rhou1t0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x1y0z0Val*rhou1Val*rhou2t0x0y1z0Val + a3t0x0y1z0Val*gamma*Power(rhou1Val,2)*rhou2t0x0y2z0Val - a3t0x0y1z0Val*Power(gamma,2)*Power(rhou1Val,2)*rhou2t0x0y2z0Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x0y1z0Val*rhou2t0x1y0z0Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x0y1z0Val*rhou2t0x1y0z0Val + 2*a3t0x1y0z0Val*eVal*gamma*rhot0x0y1z0Val*rhou2t0x1y0z0Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x1y0z0Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x0y1z0Val*rhou2t0x1y0z0Val + 14*a3t0x1y0z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2t0x1y0z0Val - 2*a3t0x1y0z0Val*gamma*rhou1t0x0y1z0Val*rhou1Val*rhou2t0x1y0z0Val + 18*a3t0x0y1z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2t0x1y0z0Val - 6*a3t0x0y1z0Val*gamma*rhou1t0x1y0z0Val*rhou1Val*rhou2t0x1y0z0Val + 12*a3t0x1y1z0Val*Power(rhou1Val,2)*rhou2t0x1y0z0Val + 18*a3t0x0y1z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2t0x1y0z0Val - 6*a3t0x0y1z0Val*gamma*rhou1Val*rhou2t0x0y1z0Val*rhou2t0x1y0z0Val + 6*a3t0x1y0z0Val*rhou1Val*Power(rhou2t0x1y0z0Val,2) - 6*a3t0x1y0z0Val*gamma*rhou1Val*Power(rhou2t0x1y0z0Val,2) + a3t0x1y0z0Val*gamma*Power(rhou1Val,2)*rhou2t0x1y1z0Val - a3t0x1y0z0Val*Power(gamma,2)*Power(rhou1Val,2)*rhou2t0x1y1z0Val + 6*a3t0x0y1z0Val*Power(rhou1Val,2)*rhou2t0x2y0z0Val - 28*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*rhot0x0y1z0Val*rhou2Val - 12*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*rhot0x0y1z0Val*rhou2Val + 12*a3t0x0y1z0Val*et0x0y1z0Val*rhot0x0y1z0Val*rhou2Val + 4*a3t0x1y0z0Val*et0x1y0z0Val*rhot0x0y1z0Val*rhou2Val + 10*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou2Val + 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou2Val - 10*a3t0x0y1z0Val*et0x0y1z0Val*gamma*rhot0x0y1z0Val*rhou2Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*gamma*rhot0x0y1z0Val*rhou2Val + 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val + 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val*rhou2Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*rhot0x0y2z0Val*rhou2Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhot0x0y2z0Val*rhou2Val + 2*a3t0x0y1z0Val*eVal*gamma*rhot0x0y2z0Val*rhou2Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x0y2z0Val*rhou2Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhot0x0y2z0Val*rhou2Val - 4*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*rhot0x1y0z0Val*rhou2Val - 12*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*rhot0x1y0z0Val*rhou2Val + 2*a3t0x1y0z0Val*et0x0y1z0Val*rhot0x1y0z0Val*rhou2Val + 6*a3t0x0y1z0Val*et0x1y0z0Val*rhot0x1y0z0Val*rhou2Val + 2*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou2Val + 6*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*et0x0y1z0Val*gamma*rhot0x1y0z0Val*rhou2Val - 6*a3t0x0y1z0Val*et0x1y0z0Val*gamma*rhot0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhot0x1y1z0Val*rhou2Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhot0x1y1z0Val*rhou2Val + 2*a3t0x1y0z0Val*eVal*gamma*rhot0x1y1z0Val*rhou2Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhot0x1y1z0Val*rhou2Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhot0x1y1z0Val*rhou2Val + 6*a3t0x0y1z0Val*Power(rhou1t0x0y1z0Val,2)*rhou2Val - 6*a3t0x0y1z0Val*gamma*Power(rhou1t0x0y1z0Val,2)*rhou2Val + 18*a3t0x1y0z0Val*rhou1t0x0y1z0Val*rhou1t0x1y0z0Val*rhou2Val - 6*a3t0x1y0z0Val*gamma*rhou1t0x0y1z0Val*rhou1t0x1y0z0Val*rhou2Val + 6*a3t0x0y2z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val + 12*a3t0x2y0z0Val*rhou1t0x0y1z0Val*rhou1Val*rhou2Val - 6*a3t0x0y2z0Val*gamma*rhou1t0x0y1z0Val*rhou1Val*rhou2Val + 6*a3t0x0y1z0Val*rhou1t0x0y2z0Val*rhou1Val*rhou2Val - 6*a3t0x0y1z0Val*gamma*rhou1t0x0y2z0Val*rhou1Val*rhou2Val + 18*a3t0x1y1z0Val*rhou1t0x1y0z0Val*rhou1Val*rhou2Val - 6*a3t0x1y1z0Val*gamma*rhou1t0x1y0z0Val*rhou1Val*rhou2Val + 18*a3t0x1y0z0Val*rhou1t0x1y1z0Val*rhou1Val*rhou2Val - 6*a3t0x1y0z0Val*gamma*rhou1t0x1y1z0Val*rhou1Val*rhou2Val + 6*a3t0x2y1z0Val*Power(rhou1Val,2)*rhou2Val + 18*a3t0x1y0z0Val*rhou1t0x0y1z0Val*rhou2t0x0y1z0Val*rhou2Val - 6*a3t0x1y0z0Val*gamma*rhou1t0x0y1z0Val*rhou2t0x0y1z0Val*rhou2Val - 2*a3t0x0y1z0Val*rhou1t0x1y0z0Val*rhou2t0x0y1z0Val*rhou2Val + 4*a3t0x0y1z0Val*gamma*rhou1t0x1y0z0Val*rhou2t0x0y1z0Val*rhou2Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x1y0z0Val*rhou2t0x0y1z0Val*rhou2Val + 18*a3t0x1y1z0Val*rhou1Val*rhou2t0x0y1z0Val*rhou2Val - 6*a3t0x1y1z0Val*gamma*rhou1Val*rhou2t0x0y1z0Val*rhou2Val + 24*a3t0x0y1z0Val*Power(rhou2t0x0y1z0Val,2)*rhou2Val - 10*a3t0x0y1z0Val*gamma*Power(rhou2t0x0y1z0Val,2)*rhou2Val - 2*a3t0x0y1z0Val*Power(gamma,2)*Power(rhou2t0x0y1z0Val,2)*rhou2Val + 14*a3t0x0y1z0Val*rhou1t0x0y1z0Val*rhou2t0x1y0z0Val*rhou2Val - 2*a3t0x0y1z0Val*gamma*rhou1t0x0y1z0Val*rhou2t0x1y0z0Val*rhou2Val + 6*a3t0x1y0z0Val*rhou1t0x1y0z0Val*rhou2t0x1y0z0Val*rhou2Val - 4*a3t0x1y0z0Val*gamma*rhou1t0x1y0z0Val*rhou2t0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x1y0z0Val*rhou2t0x1y0z0Val*rhou2Val + 12*a3t0x0y2z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val + 6*a3t0x2y0z0Val*rhou1Val*rhou2t0x1y0z0Val*rhou2Val - 6*a3t0x2y0z0Val*gamma*rhou1Val*rhou2t0x1y0z0Val*rhou2Val + 6*a3t0x1y0z0Val*rhou2t0x0y1z0Val*rhou2t0x1y0z0Val*rhou2Val - 4*a3t0x1y0z0Val*gamma*rhou2t0x0y1z0Val*rhou2t0x1y0z0Val*rhou2Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhou2t0x0y1z0Val*rhou2t0x1y0z0Val*rhou2Val + 6*a3t0x0y1z0Val*Power(rhou2t0x1y0z0Val,2)*rhou2Val - 6*a3t0x0y1z0Val*gamma*Power(rhou2t0x1y0z0Val,2)*rhou2Val + 18*a3t0x0y1z0Val*rhou1Val*rhou2t0x1y1z0Val*rhou2Val - 6*a3t0x0y1z0Val*gamma*rhou1Val*rhou2t0x1y1z0Val*rhou2Val + 6*a3t0x1y0z0Val*rhou1Val*rhou2t0x2y0z0Val*rhou2Val - 6*a3t0x1y0z0Val*gamma*rhou1Val*rhou2t0x2y0z0Val*rhou2Val + 12*a3t0x1y1z0Val*rhou1t0x0y1z0Val*Power(rhou2Val,2) + 6*a3t0x1y0z0Val*rhou1t0x0y2z0Val*Power(rhou2Val,2) + a3t0x0y1z0Val*gamma*rhou1t0x1y1z0Val*Power(rhou2Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x1y1z0Val*Power(rhou2Val,2) + a3t0x1y0z0Val*gamma*rhou1t0x2y0z0Val*Power(rhou2Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x2y0z0Val*Power(rhou2Val,2) + 6*a3t0x1y2z0Val*rhou1Val*Power(rhou2Val,2) + 18*a3t0x0y2z0Val*rhou2t0x0y1z0Val*Power(rhou2Val,2) - 6*a3t0x0y2z0Val*gamma*rhou2t0x0y1z0Val*Power(rhou2Val,2) + 12*a3t0x0y1z0Val*rhou2t0x0y2z0Val*Power(rhou2Val,2) - 5*a3t0x0y1z0Val*gamma*rhou2t0x0y2z0Val*Power(rhou2Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhou2t0x0y2z0Val*Power(rhou2Val,2) + 6*a3t0x1y1z0Val*rhou2t0x1y0z0Val*Power(rhou2Val,2) - 6*a3t0x1y1z0Val*gamma*rhou2t0x1y0z0Val*Power(rhou2Val,2) + 6*a3t0x1y0z0Val*rhou2t0x1y1z0Val*Power(rhou2Val,2) - 5*a3t0x1y0z0Val*gamma*rhou2t0x1y1z0Val*Power(rhou2Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhou2t0x1y1z0Val*Power(rhou2Val,2) + 2*a3t0x0y3z0Val*Power(rhou2Val,3) + Power(B2Val,2)*(a3t0x0y1z0Val*(-1 + gamma)*gamma*(rhot0x1y0z0Val*rhou1t0x0y1z0Val + rhot0x1y1z0Val*rhou1Val + rhot0x0y1z0Val*(rhou1t0x1y0z0Val + 2*rhou2t0x0y1z0Val) + rhot0x0y2z0Val*rhou2Val) - a3t0x1y0z0Val*(4*rhot0x1y0z0Val*rhou1t0x1y0z0Val + 2*gamma*rhot0x1y0z0Val*rhou1t0x1y0z0Val - 2*Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x1y0z0Val + (2*rhot0x0y2z0Val + (2 + gamma - Power(gamma,2))*rhot0x2y0z0Val)*rhou1Val + gamma*rhot0x1y0z0Val*rhou2t0x0y1z0Val - Power(gamma,2)*rhot0x1y0z0Val*rhou2t0x0y1z0Val + rhot0x0y1z0Val*(4*rhou1t0x0y1z0Val - (-1 + gamma)*gamma*rhou2t0x1y0z0Val) + gamma*rhot0x1y1z0Val*rhou2Val - Power(gamma,2)*rhot0x1y1z0Val*rhou2Val)) + Power(B1Val,2)*(a3t0x1y0z0Val*(-1 + gamma)*gamma*(2*rhot0x1y0z0Val*rhou1t0x1y0z0Val + rhot0x2y0z0Val*rhou1Val + rhot0x1y0z0Val*rhou2t0x0y1z0Val + rhot0x0y1z0Val*rhou2t0x1y0z0Val + rhot0x1y1z0Val*rhou2Val) + a3t0x0y1z0Val*(-(gamma*rhot0x1y0z0Val*rhou1t0x0y1z0Val) + Power(gamma,2)*rhot0x1y0z0Val*rhou1t0x0y1z0Val - gamma*rhot0x1y1z0Val*rhou1Val + Power(gamma,2)*rhot0x1y1z0Val*rhou1Val + rhot0x0y1z0Val*((-1 + gamma)*gamma*rhou1t0x1y0z0Val + 2*(-2 - gamma + Power(gamma,2))*rhou2t0x0y1z0Val) - 4*rhot0x1y0z0Val*rhou2t0x1y0z0Val + ((-2 - gamma + Power(gamma,2))*rhot0x0y2z0Val - 2*rhot0x2y0z0Val)*rhou2Val)) + 6*a3t0x0y1z0Val*rhou2Val*Power(rhou3t0x0y1z0Val,2) - 6*a3t0x0y1z0Val*gamma*rhou2Val*Power(rhou3t0x0y1z0Val,2) + 6*a3t0x0y1z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3t0x1y0z0Val - 6*a3t0x0y1z0Val*gamma*rhou1Val*rhou3t0x0y1z0Val*rhou3t0x1y0z0Val + 6*a3t0x1y0z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3t0x1y0z0Val - 6*a3t0x1y0z0Val*gamma*rhou2Val*rhou3t0x0y1z0Val*rhou3t0x1y0z0Val + 6*a3t0x1y0z0Val*rhou1Val*Power(rhou3t0x1y0z0Val,2) - 6*a3t0x1y0z0Val*gamma*rhou1Val*Power(rhou3t0x1y0z0Val,2) + 4*a3t0x0y1z0Val*B1t0x1y0z0Val*B3Val*rhot0x0y1z0Val*rhou3Val + 6*a3t0x0y1z0Val*B2t0x0y1z0Val*B3Val*rhot0x0y1z0Val*rhou3Val + 2*a3t0x1y0z0Val*B2t0x1y0z0Val*B3Val*rhot0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*B1t0x1y0z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B3Val*gamma*rhot0x0y1z0Val*rhou3Val + 2*a3t0x0y1z0Val*B1t0x0y1z0Val*B3Val*rhot0x1y0z0Val*rhou3Val + 6*a3t0x1y0z0Val*B1t0x1y0z0Val*B3Val*rhot0x1y0z0Val*rhou3Val + 4*a3t0x1y0z0Val*B2t0x0y1z0Val*B3Val*rhot0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*B2t0x0y1z0Val*B3Val*gamma*rhot0x1y0z0Val*rhou3Val + 6*a3t0x1y0z0Val*rhou1t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val - 6*a3t0x1y0z0Val*gamma*rhou1t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*rhou1t0x1y0z0Val*rhou3t0x0y1z0Val*rhou3Val + 4*a3t0x0y1z0Val*gamma*rhou1t0x1y0z0Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x1y0z0Val*rhou3t0x0y1z0Val*rhou3Val + 6*a3t0x1y1z0Val*rhou1Val*rhou3t0x0y1z0Val*rhou3Val - 6*a3t0x1y1z0Val*gamma*rhou1Val*rhou3t0x0y1z0Val*rhou3Val + 6*a3t0x0y1z0Val*rhou2t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val - 4*a3t0x0y1z0Val*gamma*rhou2t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x0y1z0Val*Power(gamma,2)*rhou2t0x0y1z0Val*rhou3t0x0y1z0Val*rhou3Val + 2*a3t0x1y0z0Val*rhou2t0x1y0z0Val*rhou3t0x0y1z0Val*rhou3Val - 2*a3t0x1y0z0Val*gamma*rhou2t0x1y0z0Val*rhou3t0x0y1z0Val*rhou3Val + 6*a3t0x0y2z0Val*rhou2Val*rhou3t0x0y1z0Val*rhou3Val - 6*a3t0x0y2z0Val*gamma*rhou2Val*rhou3t0x0y1z0Val*rhou3Val + 6*a3t0x0y1z0Val*rhou2Val*rhou3t0x0y2z0Val*rhou3Val - 6*a3t0x0y1z0Val*gamma*rhou2Val*rhou3t0x0y2z0Val*rhou3Val + 2*a3t0x0y1z0Val*rhou1t0x0y1z0Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x0y1z0Val*gamma*rhou1t0x0y1z0Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x1y0z0Val*rhou1t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val - 4*a3t0x1y0z0Val*gamma*rhou1t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x2y0z0Val*rhou1Val*rhou3t0x1y0z0Val*rhou3Val - 6*a3t0x2y0z0Val*gamma*rhou1Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*rhou2t0x0y1z0Val*rhou3t0x1y0z0Val*rhou3Val + 4*a3t0x1y0z0Val*gamma*rhou2t0x0y1z0Val*rhou3t0x1y0z0Val*rhou3Val - 2*a3t0x1y0z0Val*Power(gamma,2)*rhou2t0x0y1z0Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x0y1z0Val*rhou2t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val - 6*a3t0x0y1z0Val*gamma*rhou2t0x1y0z0Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x1y1z0Val*rhou2Val*rhou3t0x1y0z0Val*rhou3Val - 6*a3t0x1y1z0Val*gamma*rhou2Val*rhou3t0x1y0z0Val*rhou3Val + 6*a3t0x0y1z0Val*rhou1Val*rhou3t0x1y1z0Val*rhou3Val - 6*a3t0x0y1z0Val*gamma*rhou1Val*rhou3t0x1y1z0Val*rhou3Val + 6*a3t0x1y0z0Val*rhou2Val*rhou3t0x1y1z0Val*rhou3Val - 6*a3t0x1y0z0Val*gamma*rhou2Val*rhou3t0x1y1z0Val*rhou3Val + 6*a3t0x1y0z0Val*rhou1Val*rhou3t0x2y0z0Val*rhou3Val - 6*a3t0x1y0z0Val*gamma*rhou1Val*rhou3t0x2y0z0Val*rhou3Val + a3t0x0y1z0Val*gamma*rhou1t0x1y1z0Val*Power(rhou3Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhou1t0x1y1z0Val*Power(rhou3Val,2) + a3t0x1y0z0Val*gamma*rhou1t0x2y0z0Val*Power(rhou3Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhou1t0x2y0z0Val*Power(rhou3Val,2) + a3t0x0y1z0Val*gamma*rhou2t0x0y2z0Val*Power(rhou3Val,2) - a3t0x0y1z0Val*Power(gamma,2)*rhou2t0x0y2z0Val*Power(rhou3Val,2) + a3t0x1y0z0Val*gamma*rhou2t0x1y1z0Val*Power(rhou3Val,2) - a3t0x1y0z0Val*Power(gamma,2)*rhou2t0x1y1z0Val*Power(rhou3Val,2) + 2*B2Val*(2*a3t0x0y1z0Val*B1Val*rhot0x0y1z0Val*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*B1Val*rhot0x1y0z0Val*rhou1t0x1y0z0Val + (a3t0x1y0z0Val*(2*B1t0x1y0z0Val*rhot0x0y1z0Val + 3*B2t0x0y1z0Val*(-1 + gamma)*rhot0x0y1z0Val + (6*B1t0x0y1z0Val + B2t0x1y0z0Val*(-14 + 5*gamma + Power(gamma,2)))*rhot0x1y0z0Val) + a3t0x0y1z0Val*(B1t0x0y1z0Val*rhot0x0y1z0Val - 2*B2t0x1y0z0Val*rhot0x0y1z0Val + B2t0x1y0z0Val*gamma*rhot0x0y1z0Val + 3*B1t0x1y0z0Val*rhot0x1y0z0Val + B2t0x0y1z0Val*rhot0x1y0z0Val + B2t0x0y1z0Val*gamma*rhot0x1y0z0Val + B2t0x0y1z0Val*Power(gamma,2)*rhot0x1y0z0Val + B1Val*(rhot0x0y2z0Val + rhot0x2y0z0Val)))*rhou1Val + 2*a3t0x1y0z0Val*B1Val*rhot0x0y1z0Val*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B1Val*rhot0x1y0z0Val*rhou2t0x1y0z0Val + (a3t0x0y1z0Val*(-(B1t0x1y0z0Val*(-6 + gamma)*rhot0x0y1z0Val) + B2t0x0y1z0Val*gamma*(4 + gamma)*rhot0x0y1z0Val + 3*(B1t0x0y1z0Val + B2t0x1y0z0Val*(-2 + gamma))*rhot0x1y0z0Val) + a3t0x1y0z0Val*(4*B1t0x0y1z0Val*rhot0x0y1z0Val - 4*B2t0x1y0z0Val*rhot0x0y1z0Val + B2t0x1y0z0Val*gamma*rhot0x0y1z0Val + B2t0x1y0z0Val*Power(gamma,2)*rhot0x0y1z0Val + 3*B1t0x1y0z0Val*rhot0x1y0z0Val - B1t0x1y0z0Val*gamma*rhot0x1y0z0Val + B1Val*(rhot0x0y2z0Val + rhot0x2y0z0Val)))*rhou2Val + 2*a3t0x0y1z0Val*B3Val*rhot0x0y1z0Val*rhou3t0x0y1z0Val + a3t0x1y0z0Val*B3Val*rhot0x1y0z0Val*rhou3t0x0y1z0Val + a3t0x1y0z0Val*B3Val*rhot0x0y1z0Val*rhou3t0x1y0z0Val + a3t0x0y1z0Val*B3t0x0y1z0Val*rhot0x0y1z0Val*rhou3Val + a3t0x1y0z0Val*B3t0x1y0z0Val*rhot0x0y1z0Val*rhou3Val + a3t0x0y1z0Val*B3Val*rhot0x0y2z0Val*rhou3Val + a3t0x1y0z0Val*B3Val*rhot0x1y1z0Val*rhou3Val) + 2*B1Val*((a3t0x1y0z0Val*(3*B2t0x1y0z0Val*rhot0x0y1z0Val + 3*B1t0x0y1z0Val*(-2 + gamma)*rhot0x0y1z0Val + (-(B2t0x0y1z0Val*(-6 + gamma)) + B1t0x1y0z0Val*gamma*(4 + gamma))*rhot0x1y0z0Val) + a3t0x0y1z0Val*(-(B2t0x0y1z0Val*(-3 + gamma)*rhot0x0y1z0Val) + (4*B2t0x1y0z0Val + B1t0x0y1z0Val*(-4 + gamma + Power(gamma,2)))*rhot0x1y0z0Val))*rhou1Val + (a3t0x1y0z0Val*(3*B2t0x0y1z0Val*rhot0x0y1z0Val + B1t0x1y0z0Val*(1 + gamma + Power(gamma,2))*rhot0x0y1z0Val + (B2t0x1y0z0Val + B1t0x0y1z0Val*(-2 + gamma))*rhot0x1y0z0Val) + a3t0x0y1z0Val*(6*B2t0x1y0z0Val*rhot0x0y1z0Val + B1t0x0y1z0Val*(-14 + 5*gamma + Power(gamma,2))*rhot0x0y1z0Val + (2*B2t0x0y1z0Val + 3*B1t0x1y0z0Val*(-1 + gamma))*rhot0x1y0z0Val))*rhou2Val + a3t0x0y1z0Val*B3Val*rhot0x1y0z0Val*rhou3t0x0y1z0Val + a3t0x0y1z0Val*B3Val*rhot0x0y1z0Val*rhou3t0x1y0z0Val + 2*a3t0x1y0z0Val*B3Val*rhot0x1y0z0Val*rhou3t0x1y0z0Val + a3t0x0y1z0Val*B3t0x0y1z0Val*rhot0x1y0z0Val*rhou3Val + a3t0x1y0z0Val*B3t0x1y0z0Val*rhot0x1y0z0Val*rhou3Val + a3t0x0y1z0Val*B3Val*rhot0x1y1z0Val*rhou3Val + a3t0x1y0z0Val*B3Val*rhot0x2y0z0Val*rhou3Val))*Power(rhoVal,2) + (-12*a3t0x1y0z0Val*B1t0x0y1z0Val*B1Val*rhou1t0x0y1z0Val + 6*a3t0x0y1z0Val*B1Val*B2t0x0y1z0Val*rhou1t0x0y1z0Val + 6*a3t0x1y0z0Val*B1Val*B2t0x1y0z0Val*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*B1t0x0y1z0Val*B2Val*rhou1t0x0y1z0Val + 4*a3t0x1y0z0Val*B1t0x1y0z0Val*B2Val*rhou1t0x0y1z0Val - 6*a3t0x1y0z0Val*B2t0x0y1z0Val*B2Val*rhou1t0x0y1z0Val - 4*a3t0x0y1z0Val*B2t0x1y0z0Val*B2Val*rhou1t0x0y1z0Val - 12*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*rhou1t0x0y1z0Val - 4*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*rhou1t0x0y1z0Val + 6*a3t0x1y0z0Val*et0x0y1z0Val*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*et0x1y0z0Val*rhou1t0x0y1z0Val + 6*a3t0x1y0z0Val*B1t0x0y1z0Val*B1Val*gamma*rhou1t0x0y1z0Val - 2*a3t0x0y1z0Val*B1Val*B2t0x0y1z0Val*gamma*rhou1t0x0y1z0Val + 6*a3t0x1y0z0Val*B2t0x0y1z0Val*B2Val*gamma*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*B2t0x1y0z0Val*B2Val*gamma*rhou1t0x0y1z0Val + 6*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*gamma*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*gamma*rhou1t0x0y1z0Val - 6*a3t0x1y0z0Val*et0x0y1z0Val*gamma*rhou1t0x0y1z0Val - 2*a3t0x0y1z0Val*et0x1y0z0Val*gamma*rhou1t0x0y1z0Val + 2*a3t0x0y1z0Val*B1Val*B2Val*rhou1t0x0y2z0Val - 2*a3t0x1y0z0Val*Power(B2Val,2)*rhou1t0x0y2z0Val + 4*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*rhou1t0x1y0z0Val + 6*a3t0x1y0z0Val*B1Val*B2t0x0y1z0Val*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*B1Val*B2t0x1y0z0Val*rhou1t0x1y0z0Val + 6*a3t0x1y0z0Val*B1t0x0y1z0Val*B2Val*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B2Val*rhou1t0x1y0z0Val - 16*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*rhou1t0x1y0z0Val - 16*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*rhou1t0x1y0z0Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*rhou1t0x1y0z0Val + 6*a3t0x1y0z0Val*et0x1y0z0Val*rhou1t0x1y0z0Val - 4*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*gamma*rhou1t0x1y0z0Val + 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B1Val*gamma*rhou1t0x1y0z0Val - 2*a3t0x1y0z0Val*B1Val*B2t0x0y1z0Val*gamma*rhou1t0x1y0z0Val - 4*a3t0x0y1z0Val*B2t0x0y1z0Val*B2Val*gamma*rhou1t0x1y0z0Val + 4*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*gamma*rhou1t0x1y0z0Val - 4*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*gamma*rhou1t0x1y0z0Val + 4*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*gamma*rhou1t0x1y0z0Val + 4*a3t0x0y1z0Val*et0x0y1z0Val*gamma*rhou1t0x1y0z0Val - 4*a3t0x1y0z0Val*et0x1y0z0Val*gamma*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*Power(gamma,2)*rhou1t0x1y0z0Val + 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B1Val*Power(gamma,2)*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B2Val*Power(gamma,2)*rhou1t0x1y0z0Val + 2*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*Power(gamma,2)*rhou1t0x1y0z0Val + 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*Power(gamma,2)*rhou1t0x1y0z0Val + 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*Power(gamma,2)*rhou1t0x1y0z0Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*Power(gamma,2)*rhou1t0x1y0z0Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*Power(gamma,2)*rhou1t0x1y0z0Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*rhou1t0x1y1z0Val - a3t0x0y1z0Val*Power(B1Val,2)*gamma*rhou1t0x1y1z0Val - a3t0x0y1z0Val*Power(B2Val,2)*gamma*rhou1t0x1y1z0Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhou1t0x1y1z0Val + 2*a3t0x0y1z0Val*eVal*gamma*rhou1t0x1y1z0Val + a3t0x0y1z0Val*Power(B1Val,2)*Power(gamma,2)*rhou1t0x1y1z0Val + a3t0x0y1z0Val*Power(B2Val,2)*Power(gamma,2)*rhou1t0x1y1z0Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhou1t0x1y1z0Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhou1t0x1y1z0Val + 2*a3t0x0y1z0Val*B1Val*B2Val*rhou1t0x2y0z0Val - 2*a3t0x1y0z0Val*Power(B2Val,2)*rhou1t0x2y0z0Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhou1t0x2y0z0Val - a3t0x1y0z0Val*Power(B1Val,2)*gamma*rhou1t0x2y0z0Val - a3t0x1y0z0Val*Power(B2Val,2)*gamma*rhou1t0x2y0z0Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhou1t0x2y0z0Val + 2*a3t0x1y0z0Val*eVal*gamma*rhou1t0x2y0z0Val + a3t0x1y0z0Val*Power(B1Val,2)*Power(gamma,2)*rhou1t0x2y0z0Val + a3t0x1y0z0Val*Power(B2Val,2)*Power(gamma,2)*rhou1t0x2y0z0Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhou1t0x2y0z0Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhou1t0x2y0z0Val - 2*(6*a3t0x1y1z0Val*B1t0x0y1z0Val*B1Val + a3t0x1y0z0Val*Power(B2t0x0y1z0Val,2) - 3*a3t0x1y0z0Val*B1t0x0y1z0Val*B2t0x1y0z0Val - 3*a3t0x1y1z0Val*B1Val*B2t0x1y0z0Val + 6*a3t0x1y0z0Val*Power(B2t0x1y0z0Val,2) - 3*a3t0x1y0z0Val*B1Val*B2t0x1y1z0Val - 3*a3t0x2y0z0Val*B1t0x0y1z0Val*B2Val - 3*a3t0x1y1z0Val*B1t0x1y0z0Val*B2Val - 2*a3t0x1y0z0Val*B1t0x1y1z0Val*B2Val + a3t0x1y0z0Val*B2t0x0y2z0Val*B2Val + 6*a3t0x2y0z0Val*B2t0x1y0z0Val*B2Val + 6*a3t0x1y0z0Val*B2t0x2y0z0Val*B2Val + 6*a3t0x1y0z0Val*Power(B3t0x1y0z0Val,2) + 6*a3t0x1y1z0Val*B3t0x0y1z0Val*B3Val + 6*a3t0x2y0z0Val*B3t0x1y0z0Val*B3Val + 6*a3t0x1y0z0Val*B3t0x2y0z0Val*B3Val - 3*a3t0x1y1z0Val*et0x0y1z0Val - 3*a3t0x2y0z0Val*et0x1y0z0Val - 3*a3t0x1y0z0Val*et0x2y0z0Val - 2*a3t0x1y0z0Val*Power(B1t0x1y0z0Val,2)*gamma - 3*a3t0x1y1z0Val*B1t0x0y1z0Val*B1Val*gamma - 3*a3t0x2y0z0Val*B1t0x1y0z0Val*B1Val*gamma - 2*a3t0x1y0z0Val*B1t0x2y0z0Val*B1Val*gamma - 3*a3t0x1y0z0Val*Power(B2t0x1y0z0Val,2)*gamma + a3t0x1y0z0Val*B1Val*B2t0x1y1z0Val*gamma - 3*a3t0x2y0z0Val*B2t0x1y0z0Val*B2Val*gamma - 3*a3t0x1y0z0Val*B2t0x2y0z0Val*B2Val*gamma - 3*a3t0x1y0z0Val*Power(B3t0x1y0z0Val,2)*gamma - 3*a3t0x1y1z0Val*B3t0x0y1z0Val*B3Val*gamma - 3*a3t0x2y0z0Val*B3t0x1y0z0Val*B3Val*gamma - 3*a3t0x1y0z0Val*B3t0x2y0z0Val*B3Val*gamma + 3*a3t0x1y1z0Val*et0x0y1z0Val*gamma + 3*a3t0x2y0z0Val*et0x1y0z0Val*gamma + 3*a3t0x1y0z0Val*et0x2y0z0Val*gamma + a3t0x0y1z0Val*(-5*B1t0x1y0z0Val*B2t0x1y0z0Val + B2t0x0y1z0Val*B2t0x1y0z0Val - 2*B1t0x2y0z0Val*B2Val + B2t0x1y1z0Val*B2Val + 6*B3t0x0y1z0Val*B3t0x1y0z0Val + 6*B3t0x1y1z0Val*B3Val - 3*et0x1y1z0Val + B1Val*(-3*B2t0x2y0z0Val - 2*B1t0x1y1z0Val*(-2 + gamma) + B2t0x0y2z0Val*(-2 + gamma)) + B1t0x0y1z0Val*(-2*B1t0x1y0z0Val + B2t0x0y1z0Val)*(-2 + gamma) - 3*B2t0x0y1z0Val*B2t0x1y0z0Val*gamma - 3*B2t0x1y1z0Val*B2Val*gamma - 3*B3t0x0y1z0Val*B3t0x1y0z0Val*gamma - 3*B3t0x1y1z0Val*B3Val*gamma + 3*et0x1y1z0Val*gamma) + B2t0x0y1z0Val*(a3t0x1y0z0Val*B1t0x1y0z0Val*(-2 + gamma) - 3*(a3t0x2y0z0Val*B1Val + a3t0x1y1z0Val*B2Val*gamma)))*rhou1Val - 16*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B1Val*rhou2t0x0y1z0Val + 6*a3t0x0y1z0Val*B1Val*B2t0x1y0z0Val*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B1t0x0y1z0Val*B2Val*rhou2t0x0y1z0Val + 6*a3t0x0y1z0Val*B1t0x1y0z0Val*B2Val*rhou2t0x0y1z0Val + 4*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*rhou2t0x0y1z0Val - 16*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*rhou2t0x0y1z0Val + 6*a3t0x0y1z0Val*et0x0y1z0Val*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*rhou2t0x0y1z0Val + 4*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*gamma*rhou2t0x0y1z0Val - 4*a3t0x1y0z0Val*B1t0x1y0z0Val*B1Val*gamma*rhou2t0x0y1z0Val - 2*a3t0x0y1z0Val*B1t0x1y0z0Val*B2Val*gamma*rhou2t0x0y1z0Val + 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B2Val*gamma*rhou2t0x0y1z0Val - 4*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*gamma*rhou2t0x0y1z0Val + 4*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*gamma*rhou2t0x0y1z0Val - 4*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*gamma*rhou2t0x0y1z0Val - 4*a3t0x0y1z0Val*et0x0y1z0Val*gamma*rhou2t0x0y1z0Val + 4*a3t0x1y0z0Val*et0x1y0z0Val*gamma*rhou2t0x0y1z0Val + 2*a3t0x0y1z0Val*B1t0x0y1z0Val*B1Val*Power(gamma,2)*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B1Val*Power(gamma,2)*rhou2t0x0y1z0Val + 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B2Val*Power(gamma,2)*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B2t0x1y0z0Val*B2Val*Power(gamma,2)*rhou2t0x0y1z0Val + 2*a3t0x0y1z0Val*B3t0x0y1z0Val*B3Val*Power(gamma,2)*rhou2t0x0y1z0Val + 2*a3t0x1y0z0Val*B3t0x1y0z0Val*B3Val*Power(gamma,2)*rhou2t0x0y1z0Val - 2*a3t0x0y1z0Val*et0x0y1z0Val*Power(gamma,2)*rhou2t0x0y1z0Val - 2*a3t0x1y0z0Val*et0x1y0z0Val*Power(gamma,2)*rhou2t0x0y1z0Val - 2*a3t0x0y1z0Val*Power(B1Val,2)*rhou2t0x0y2z0Val + 2*a3t0x1y0z0Val*B1Val*B2Val*rhou2t0x0y2z0Val - 2*a3t0x0y1z0Val*Power(B3Val,2)*rhou2t0x0y2z0Val - a3t0x0y1z0Val*Power(B1Val,2)*gamma*rhou2t0x0y2z0Val - a3t0x0y1z0Val*Power(B2Val,2)*gamma*rhou2t0x0y2z0Val - a3t0x0y1z0Val*Power(B3Val,2)*gamma*rhou2t0x0y2z0Val + 2*a3t0x0y1z0Val*eVal*gamma*rhou2t0x0y2z0Val + a3t0x0y1z0Val*Power(B1Val,2)*Power(gamma,2)*rhou2t0x0y2z0Val + a3t0x0y1z0Val*Power(B2Val,2)*Power(gamma,2)*rhou2t0x0y2z0Val + a3t0x0y1z0Val*Power(B3Val,2)*Power(gamma,2)*rhou2t0x0y2z0Val - 2*a3t0x0y1z0Val*eVal*Power(gamma,2)*rhou2t0x0y2z0Val - 4*a3t0x1y0z0Val*B1t0x0y1z0Val*B1Val*rhou2t0x1y0z0Val - 6*a3t0x0y1z0Val*B1t0x1y0z0Val*B1Val*rhou2t0x1y0z0Val + 4*a3t0x0y1z0Val*B1Val*B2t0x0y1z0Val*rhou2t0x1y0z0Val + 2*a3t0x1y0z0Val*B1Val*B2t0x1y0z0Val*rhou2t0x1y0z0Val + 6*a3t0x0y1z0Val*B1t0x0y1z0Val*B2Val*rhou2t0x1y0z0Val + 6*a3t0x1y0z0Val*B1t0x1y0z0Val*B2Val*rhou2t0x1y0z0Val - 12*a3t0x0y1z0Val*B2t0x1y0z0Val*B2Val*rhou2t0x1y0z0Val - 4*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*rhou2t0x1y0z0Val - 12*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*rhou2t0x1y0z0Val + 2*a3t0x1y0z0Val*et0x0y1z0Val*rhou2t0x1y0z0Val + 6*a3t0x0y1z0Val*et0x1y0z0Val*rhou2t0x1y0z0Val + 2*a3t0x1y0z0Val*B1t0x0y1z0Val*B1Val*gamma*rhou2t0x1y0z0Val + 6*a3t0x0y1z0Val*B1t0x1y0z0Val*B1Val*gamma*rhou2t0x1y0z0Val - 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B2Val*gamma*rhou2t0x1y0z0Val + 6*a3t0x0y1z0Val*B2t0x1y0z0Val*B2Val*gamma*rhou2t0x1y0z0Val + 2*a3t0x1y0z0Val*B3t0x0y1z0Val*B3Val*gamma*rhou2t0x1y0z0Val + 6*a3t0x0y1z0Val*B3t0x1y0z0Val*B3Val*gamma*rhou2t0x1y0z0Val - 2*a3t0x1y0z0Val*et0x0y1z0Val*gamma*rhou2t0x1y0z0Val - 6*a3t0x0y1z0Val*et0x1y0z0Val*gamma*rhou2t0x1y0z0Val - 2*a3t0x1y0z0Val*Power(B3Val,2)*rhou2t0x1y1z0Val - a3t0x1y0z0Val*Power(B1Val,2)*gamma*rhou2t0x1y1z0Val - a3t0x1y0z0Val*Power(B2Val,2)*gamma*rhou2t0x1y1z0Val - a3t0x1y0z0Val*Power(B3Val,2)*gamma*rhou2t0x1y1z0Val + 2*a3t0x1y0z0Val*eVal*gamma*rhou2t0x1y1z0Val + a3t0x1y0z0Val*Power(B1Val,2)*Power(gamma,2)*rhou2t0x1y1z0Val + a3t0x1y0z0Val*Power(B2Val,2)*Power(gamma,2)*rhou2t0x1y1z0Val + a3t0x1y0z0Val*Power(B3Val,2)*Power(gamma,2)*rhou2t0x1y1z0Val - 2*a3t0x1y0z0Val*eVal*Power(gamma,2)*rhou2t0x1y1z0Val - 2*a3t0x0y1z0Val*Power(B1Val,2)*rhou2t0x2y0z0Val + 2*a3t0x1y0z0Val*B1Val*B2Val*rhou2t0x2y0z0Val + 2*(-(a3t0x1y0z0Val*B1t0x0y1z0Val*B1t0x1y0z0Val) + 5*a3t0x1y0z0Val*B1t0x0y1z0Val*B2t0x0y1z0Val + 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B2t0x1y0z0Val - 4*a3t0x1y0z0Val*B2t0x0y1z0Val*B2t0x1y0z0Val + 3*a3t0x1y1z0Val*B1t0x0y1z0Val*B2Val + 3*a3t0x1y0z0Val*B1t0x0y2z0Val*B2Val + 3*a3t0x0y2z0Val*B1t0x1y0z0Val*B2Val + 2*a3t0x1y0z0Val*B1t0x2y0z0Val*B2Val - 6*a3t0x1y1z0Val*B2t0x1y0z0Val*B2Val - 4*a3t0x1y0z0Val*B2t0x1y1z0Val*B2Val - 6*a3t0x1y0z0Val*B3t0x0y1z0Val*B3t0x1y0z0Val - 6*a3t0x0y2z0Val*B3t0x0y1z0Val*B3Val - 6*a3t0x1y1z0Val*B3t0x1y0z0Val*B3Val - 6*a3t0x1y0z0Val*B3t0x1y1z0Val*B3Val + 3*a3t0x0y2z0Val*et0x0y1z0Val + 3*a3t0x1y1z0Val*et0x1y0z0Val + 3*a3t0x1y0z0Val*et0x1y1z0Val + 3*a3t0x1y0z0Val*B1t0x0y1z0Val*B1t0x1y0z0Val*gamma - a3t0x1y0z0Val*B1t0x1y0z0Val*B2t0x1y0z0Val*gamma + 2*a3t0x1y0z0Val*B2t0x0y1z0Val*B2t0x1y0z0Val*gamma - a3t0x1y0z0Val*B1t0x2y0z0Val*B2Val*gamma + 3*a3t0x0y2z0Val*B2t0x0y1z0Val*B2Val*gamma + 3*a3t0x1y1z0Val*B2t0x1y0z0Val*B2Val*gamma + 2*a3t0x1y0z0Val*B2t0x1y1z0Val*B2Val*gamma + 3*a3t0x1y0z0Val*B3t0x0y1z0Val*B3t0x1y0z0Val*gamma + 3*a3t0x0y2z0Val*B3t0x0y1z0Val*B3Val*gamma + 3*a3t0x1y1z0Val*B3t0x1y0z0Val*B3Val*gamma + 3*a3t0x1y0z0Val*B3t0x1y1z0Val*B3Val*gamma - 3*a3t0x0y2z0Val*et0x0y1z0Val*gamma - 3*a3t0x1y1z0Val*et0x1y0z0Val*gamma - 3*a3t0x1y0z0Val*et0x1y1z0Val*gamma + B1Val*(-(a3t0x1y0z0Val*B1t0x1y1z0Val) + 3*a3t0x1y1z0Val*B2t0x0y1z0Val + 2*a3t0x1y0z0Val*B2t0x0y2z0Val + 3*a3t0x0y2z0Val*B2t0x1y0z0Val + 3*a3t0x0y2z0Val*B1t0x0y1z0Val*(-2 + gamma) + 3*a3t0x1y1z0Val*B1t0x1y0z0Val*gamma + 3*a3t0x1y0z0Val*B1t0x1y1z0Val*gamma) + a3t0x0y1z0Val*(-Power(B1t0x1y0z0Val,2) - 6*B1t0x0y2z0Val*B1Val - B1t0x2y0z0Val*B1Val + 3*B1t0x0y1z0Val*B2t0x1y0z0Val + 2*B1Val*B2t0x1y1z0Val + 3*B1t0x1y1z0Val*B2Val - 6*Power(B3t0x0y1z0Val,2) - 6*B3t0x0y2z0Val*B3Val + 3*et0x0y2z0Val + 3*Power(B1t0x0y1z0Val,2)*(-2 + gamma) - B1t0x1y0z0Val*B2t0x0y1z0Val*(-2 + gamma) + 3*B1t0x0y2z0Val*B1Val*gamma + 2*Power(B2t0x0y1z0Val,2)*gamma - B1t0x1y1z0Val*B2Val*gamma + 2*B2t0x0y2z0Val*B2Val*gamma + 3*Power(B3t0x0y1z0Val,2)*gamma + 3*B3t0x0y2z0Val*B3Val*gamma - 3*et0x0y2z0Val*gamma))*rhou2Val + 2*a3t0x0y1z0Val*B2Val*B3t0x0y1z0Val*rhou3t0x0y1z0Val + 2*a3t0x1y0z0Val*B2Val*B3t0x1y0z0Val*rhou3t0x0y1z0Val + 4*a3t0x0y1z0Val*B1t0x1y0z0Val*B3Val*rhou3t0x0y1z0Val + 6*a3t0x0y1z0Val*B2t0x0y1z0Val*B3Val*rhou3t0x0y1z0Val + 2*a3t0x1y0z0Val*B2t0x1y0z0Val*B3Val*rhou3t0x0y1z0Val - 2*a3t0x0y1z0Val*B1t0x1y0z0Val*B3Val*gamma*rhou3t0x0y1z0Val - 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B3Val*gamma*rhou3t0x0y1z0Val + 2*a3t0x0y1z0Val*B2Val*B3Val*rhou3t0x0y2z0Val + 2*a3t0x0y1z0Val*B1Val*B3t0x0y1z0Val*rhou3t0x1y0z0Val + 2*a3t0x1y0z0Val*B1Val*B3t0x1y0z0Val*rhou3t0x1y0z0Val + 2*a3t0x0y1z0Val*B1t0x0y1z0Val*B3Val*rhou3t0x1y0z0Val + 6*a3t0x1y0z0Val*B1t0x1y0z0Val*B3Val*rhou3t0x1y0z0Val + 4*a3t0x1y0z0Val*B2t0x0y1z0Val*B3Val*rhou3t0x1y0z0Val - 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B3Val*gamma*rhou3t0x1y0z0Val - 2*a3t0x1y0z0Val*B2t0x0y1z0Val*B3Val*gamma*rhou3t0x1y0z0Val + 2*a3t0x0y1z0Val*B1Val*B3Val*rhou3t0x1y1z0Val + 2*a3t0x1y0z0Val*B2Val*B3Val*rhou3t0x1y1z0Val + 2*a3t0x1y0z0Val*B1Val*B3Val*rhou3t0x2y0z0Val + 4*a3t0x0y1z0Val*B1t0x1y0z0Val*B3t0x0y1z0Val*rhou3Val + 4*a3t0x0y1z0Val*B2t0x0y1z0Val*B3t0x0y1z0Val*rhou3Val + 4*a3t0x1y0z0Val*B1t0x1y0z0Val*B3t0x1y0z0Val*rhou3Val + 4*a3t0x1y0z0Val*B2t0x0y1z0Val*B3t0x1y0z0Val*rhou3Val + 4*a3t0x0y1z0Val*B1t0x1y1z0Val*B3Val*rhou3Val + 4*a3t0x1y0z0Val*B1t0x2y0z0Val*B3Val*rhou3Val + 4*a3t0x0y1z0Val*B2t0x0y2z0Val*B3Val*rhou3Val + 4*a3t0x1y0z0Val*B2t0x1y1z0Val*B3Val*rhou3Val - 2*a3t0x0y1z0Val*B1t0x1y0z0Val*B3t0x0y1z0Val*gamma*rhou3Val - 2*a3t0x0y1z0Val*B2t0x0y1z0Val*B3t0x0y1z0Val*gamma*rhou3Val - 2*a3t0x1y0z0Val*B1t0x1y0z0Val*B3t0x1y0z0Val*gamma*rhou3Val - 2*a3t0x1y0z0Val*B2t0x0y1z0Val*B3t0x1y0z0Val*gamma*rhou3Val - 2*a3t0x0y1z0Val*B1t0x1y1z0Val*B3Val*gamma*rhou3Val - 2*a3t0x1y0z0Val*B1t0x2y0z0Val*B3Val*gamma*rhou3Val - 2*a3t0x0y1z0Val*B2t0x0y2z0Val*B3Val*gamma*rhou3Val - 2*a3t0x1y0z0Val*B2t0x1y1z0Val*B3Val*gamma*rhou3Val)*Power(rhoVal,3))/(2.*Power(rhoVal,5));



                    Lauxstar.set(i, j, 1, 
                    a3t + 0.5*dt * a3tt + 1.0/6.0*dt*dt * a3ttt);
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




