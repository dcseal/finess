#include "dogdefs.h"
#include "Maxwell.h"
#include "IniParams.h"

void MaxwellFluxFunc( int n_offset, const dTensor2& Q, dTensor2& flux)
{

    // Parameters
    double cs_light_squared = global_ini_params.get_cs_light_squared();

    const int n_B1  = n_offset + M_B1 ;
    const int n_B2  = n_offset + M_B2 ;
    const int n_B3  = n_offset + M_B3 ;
    const int n_E1  = n_offset + M_E1 ;
    const int n_E2  = n_offset + M_E2 ;
    const int n_E3  = n_offset + M_E3 ;
    const int n_psi = n_offset + M_psi;
    const int n_phi = n_offset + M_phi;

    for(int j=1;j<=Q.getsize(1);j++)
    {

        // Variables
        const double& B1    = Q.get(j,n_B1 );
        const double& B2    = Q.get(j,n_B2 );
        const double& B3    = Q.get(j,n_B3 );
        const double& E1    = Q.get(j,n_E1 );
        const double& E2    = Q.get(j,n_E2 );
        const double& E3    = Q.get(j,n_E3 );
    
        // 1-component of flux function
        //
        flux.set(j,n_B1 ,  0. );
        flux.set(j,n_B2 , -E3 );
        flux.set(j,n_B3 ,  E2 );

        flux.set(j,n_E1 ,  0. );
        flux.set(j,n_E2 ,  cs_light_squared*B3 );
        flux.set(j,n_E3 , -cs_light_squared*B2 );
    
        // 2-component of flux function
        //
//      flux.set(j,n_B1 ,2,  E3 );
//      flux.set(j,n_B3 ,2, -E1 );
//      flux.set(j,n_E1 ,2, -cs_light_squared*B3 );
//      flux.set(j,n_E3 ,2,  cs_light_squared*B1 );

        // Fluxes for hyperbolic divergence cleaning

        // Magnetic divergence cleaning
        const double cp_speed_squared = cs_light_squared*global_ini_params.get_cc_sqd();
        if( n_psi<=Q.getsize(2))
        {
          const double& psi   = Q.get(j,n_psi);

          // components of 1-component of flux function involving psi
          flux.set(j,n_B1 ,  psi );
          flux.set(j,n_psi,  cp_speed_squared*B1 );

          // components of 2-component of flux function involving psi
//        flux.set(j,n_B2 ,2,  psi );
//        flux.set(j,n_psi,2,  cp_speed_squared*B2 );

        }
        else
        {
          flux.set(j,n_B1 ,  0. );
//        flux.set(j,n_B2 ,2,  0. );
        }

        // Electric field cleaning
        if(n_phi<=Q.getsize(2))
        {
          if( global_ini_params.get_clean_E_field() )
          {

              const double& phi   = Q.get(j,n_phi);
              // components of 1-component of flux function involving phi
              flux.set(j,n_E1 ,  phi );
              flux.set(j,n_phi,  cp_speed_squared*E1);

              // components of 2-component of flux function involving phi
//            flux.set(j,n_E2 ,2,  phi );
//            flux.set(j,n_phi,2,  cp_speed_squared*E2);

          }
          else
          {
              flux.set(j,n_E1 , 0.);
              flux.set(j,n_phi, 0.);
//            flux.set(j,n_E2 ,2, 0.);
//            flux.set(j,n_phi,2, 0.);
          }
        }
        else
        {
              flux.set(j,n_E1 , 0.);
//            flux.set(j,n_E2 ,2, 0.);
        }
    }
}

void ProjectLeftEig_Maxwell(int n_offset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{

    const double cs_light = global_ini_params.get_cs_light();

//  const int n_E1   = n_offset + M_E1 ;
//  const int n_E2   = n_offset + M_E2 ;
//  const int n_E3   = n_offset + M_E3 ;

//  const int n_B1   = n_offset + M_B1 ;
//  const int n_B2   = n_offset + M_B2 ;
//  const int n_B3   = n_offset + M_B3 ;

//  const int n_psi  = n_offset + M_psi;
//  const int n_phi  = n_offset + M_phi;

    for (int k=1; k<=W.getsize(2); k++)
    {

        W.set(11,k, (cs_light*Q.get(12,k) + Q.get(16,k))/(2.0*cs_light) );
        W.set(12,k, (cs_light*Q.get(13,k) - Q.get(15,k))/(2.0*cs_light) );
        W.set(13,k, Q.get(11,k) );

        W.set(14,k, Q.get(14,k) );
        W.set(15,k, (cs_light*Q.get(12,k) - Q.get(16,k))/(2.0*cs_light) );
        W.set(16,k, (cs_light*Q.get(13,k) + Q.get(15,k))/(2.0*cs_light) );

    }

}

void ProjectRightEig_Maxwell(int n_offset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{

//  const int n_B3   = n_offset + M_B3 ;
//  const int n_E3   = n_offset + M_E3 ;
//  const int n_psi  = n_offset + M_psi;
//  const int n_phi  = n_offset + M_phi;

    const double cs_light = global_ini_params.get_cs_light();
    const double cp_speed = global_ini_params.get_cc_factor()*cs_light;
    for (int k=1; k<=Q.getsize(2); k++)
    {

        Q.set(11,k, W.get(13,k) );
        Q.set(12,k, W.get(11,k) + W.get(15,k) );
        Q.set(13,k, W.get(12,k) + W.get(16,k) );

        Q.set(14,k, W.get(14,k) );
        Q.set(15,k, cs_light*( W.get(16,k) - W.get(12,k) ) );
        Q.set(16,k, cs_light*( W.get(11,k) - W.get(15,k) ) );

    }
}

/*
//
// TODO - this routine is used in the p05 and p10 versions of DoGPack
//
void SymPair_MaxwellFluxFunc( int n_offset, const dTensor2& Q, dTensor3& flux)
{
    // Parameters
    const double cs_light_squared = maxwellParams.get_cs_light_squared();

    const int n_B1  = n_offset + SM_B1 ;
    const int n_B2  = n_offset + SM_B2 ;
    const int n_E3  = n_offset + SM_E3 ;
    const int n_psi = n_offset + SM_psi;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        // Variables
        const double& B1    = Q.get(j,n_B1 );
        const double& B2    = Q.get(j,n_B2 );
        const double& E3    = Q.get(j,n_E3 );
    
        // 1-component of flux function
        //
        flux.set(j,n_B2 ,1, -E3 );
        flux.set(j,n_E3 ,1, -cs_light_squared*B2 );
    
        // 2-component of flux function
        //
        flux.set(j,n_B1 ,2,  E3 );
        flux.set(j,n_E3 ,2,  cs_light_squared*B1 );

        // fluxes for hyperbolic divergence cleaning
        //
        const double cp_speed_squared = maxwellParams.get_cp_speed_squared();
        if(n_psi<=Q.getsize(2))
        {
          const double& psi   = Q.get(j,n_psi);
          // components of 1-component of flux function involving psi
          flux.set(j,n_B1 ,1,  psi );
          flux.set(j,n_psi,1,  cp_speed_squared*B1 );
          // components of 2-component of flux function involving psi
          flux.set(j,n_B2 ,2,  psi );
          flux.set(j,n_psi,2,  cp_speed_squared*B2 );
        }
        else
        {
          flux.set(j,n_B1 ,1,  0.);
          flux.set(j,n_B2 ,2,  0.);
        }
    }
}

void SymPair_MaxwellFluxFunc1( int n_offset, const dTensor2& Q, dTensor2& flux)
{
    // Parameters
    const double cs_light_squared = maxwellParams.get_cs_light_squared();

    const int n_B1  = n_offset + SM_B1 ;
    const int n_B2  = n_offset + SM_B2 ;
    const int n_E3  = n_offset + SM_E3 ;
    const int n_psi = n_offset + SM_psi;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        // Variables
        const double& B1    = Q.get(j,n_B1 );
        const double& B2    = Q.get(j,n_B2 );
        const double& E3    = Q.get(j,n_E3 );
    
        // 1-component of flux function
        //
        flux.set(j,n_B2 , -E3 );
        flux.set(j,n_E3 , -cs_light_squared*B2 );
    
        // fluxes for hyperbolic divergence cleaning
        //
        const double cp_speed_squared = maxwellParams.get_cp_speed_squared();
        if(n_psi<=Q.getsize(2))
        {
          const double& psi   = Q.get(j,n_psi);
          // components of 1-component of flux function involving psi
          flux.set(j,n_B1 ,  psi );
          flux.set(j,n_psi,  cp_speed_squared*B1 );
        }
        else
        {
          flux.set(j,n_B1 ,  0.);
        }
    }
}

void SymPair_MaxwellFluxFunc2( int n_offset, const dTensor2& Q, dTensor2& flux)
{
    // Parameters
    const double cs_light_squared = maxwellParams.get_cs_light_squared();

    const int n_B1  = n_offset + SM_B1 ;
    const int n_B2  = n_offset + SM_B2 ;
    const int n_E3  = n_offset + SM_E3 ;
    const int n_psi = n_offset + SM_psi;

    for(int j=1;j<=Q.getsize(1);j++)
    {
        // Variables
        const double& B1    = Q.get(j,n_B1 );
        const double& B2    = Q.get(j,n_B2 );
        const double& E3    = Q.get(j,n_E3 );
    
        // 2-component of flux function
        //
        flux.set(j,n_B1 ,  E3 );
        flux.set(j,n_E3 ,  cs_light_squared*B1 );

        // fluxes for hyperbolic divergence cleaning
        //
        const double cp_speed_squared = maxwellParams.get_cp_speed_squared();
        if(n_psi<=Q.getsize(2))
        {
          const double& psi   = Q.get(j,n_psi);
          // components of 2-component of flux function involving psi
          flux.set(j,n_B2 ,  psi );
          flux.set(j,n_psi,  cp_speed_squared*B2 );
        }
        else
        {
          flux.set(j,n_B2 ,  0.);
        }
    }
}

void SymPair_ProjectLeftEig_Maxwell(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
    const int n_E3   = n_offset + SM_E3 ;
    const int n_psi  = n_offset + SM_psi;
    //
    int n_B1;
    int n_B2;
    if (ixy==1)
    {
        n_B1  = n_offset + SM_B1;
        n_B2  = n_offset + SM_B2;
    }
    else
    {
        assert(ixy==2);
        n_B1  = n_offset + SM_B2 ;
        n_B2  = n_offset + SM_B1 ;
    }

    const int m_B2_E3  = n_offset + SM_m_B2_E3;
    const int p_B2_E3  = n_offset + SM_p_B2_E3;
    //
    const int m_B1_psi = n_offset + SM_m_B1_psi;
    const int p_B1_psi = n_offset + SM_p_B1_psi;

    const double cs_light = maxwellParams.get_cs_light();
    const double cp_speed = global_ini_params.get_cc_factor()*cs_light;
    const double cs2_inv = 1./(2.*cs_light);
    const double cp2_inv = 1./(2.*cp_speed);
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(m_B2_E3 ,k, (cs_light*Q.get(n_B2 ,k) - Q.get(n_E3 ,k))*cs2_inv ); // eig -c
        W.set(p_B2_E3 ,k, (cs_light*Q.get(n_B2 ,k) + Q.get(n_E3 ,k))*cs2_inv ); // eig +c
        //
        if(n_psi<=Q.getsize(1))
        {
          W.set(m_B1_psi,k, (cp_speed*Q.get(n_B1 ,k) - Q.get(n_psi,k))*cp2_inv ); // eig -cp_speed
          W.set(p_B1_psi,k, (cp_speed*Q.get(n_B1 ,k) + Q.get(n_psi,k))*cp2_inv ); // eig +cp_speed
        }
        else
        {
          // we are not using psi, so just copy the data
          W.set(m_B1_psi,k, Q.get(n_B1,k));
        }
    }
}

void SymPair_ProjectRightEig_Maxwell(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
    const int n_E3   = n_offset + SM_E3 ;
    const int n_psi  = n_offset + SM_psi;
    //
    int n_B1;
    int n_B2;
    if (ixy==1)
    {
        n_B1  = n_offset + SM_B1;
        n_B2  = n_offset + SM_B2;
    }
    else
    {
        assert(ixy==2);
        n_B1  = n_offset + SM_B2 ;
        n_B2  = n_offset + SM_B1 ;
    }

    const int m_B2_E3  = n_offset + SM_m_B2_E3 ;
    const int p_B2_E3  = n_offset + SM_p_B2_E3 ;
    //
    const int m_B1_psi = n_offset + SM_m_B1_psi;
    const int p_B1_psi = n_offset + SM_p_B1_psi;

    const double cs_light = maxwellParams.get_cs_light();
    const double cp_speed = global_ini_params.get_cc_factor()*cs_light;
    for (int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(n_B2 ,k,            W.get(m_B2_E3 ,k) + W.get(p_B2_E3, k) );
        Q.set(n_E3 ,k, cs_light*( W.get(p_B2_E3 ,k) - W.get(m_B2_E3, k)) );
        //
        if(n_psi<=Q.getsize(1))
        {
          Q.set(n_B1 ,k,            W.get(p_B1_psi,k) + W.get(m_B1_psi,k) );
          Q.set(n_psi,k, cp_speed*( W.get(p_B1_psi,k) - W.get(m_B1_psi,k)) );
        }
        else
        {
          // not using psi, so copy the data back.
          Q.set(n_B1,k, W.get(m_B1_psi,k));
        }
    }
}

*/
