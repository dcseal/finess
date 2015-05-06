#include "dogdefs.h"
#include "Maxwell.h"
#include "MaxwellParams.h"

MaxwellParams maxwellParams;

void MaxwellFluxFunc(
    int n_offset,
    const dTensor2& Q,
    dTensor3& flux)
{
    // Parameters
    double cs_light_squared = maxwellParams.get_cs_light_squared();

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
        flux.set(j,n_B2 ,1, -E3 );
        flux.set(j,n_B3 ,1,  E2 );
        flux.set(j,n_E2 ,1,  cs_light_squared*B3 );
        flux.set(j,n_E3 ,1, -cs_light_squared*B2 );
    
        // 2-component of flux function
        //
        flux.set(j,n_B1 ,2,  E3 );
        flux.set(j,n_B3 ,2, -E1 );
        flux.set(j,n_E1 ,2, -cs_light_squared*B3 );
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
          flux.set(j,n_B1 ,1,  0. );
          flux.set(j,n_B2 ,2,  0. );
        }
        //
        if(n_phi<=Q.getsize(2))
        {
          if(maxwellParams.get_clean_E_field())
          {
              const double& phi   = Q.get(j,n_phi);
              // components of 1-component of flux function involving phi
              flux.set(j,n_E1 ,1,  phi );
              flux.set(j,n_phi,1,  cp_speed_squared*E1);
              // components of 2-component of flux function involving phi
              flux.set(j,n_E2 ,2,  phi );
              flux.set(j,n_phi,2,  cp_speed_squared*E2);
          }
          else
          {
              flux.set(j,n_E1 ,1, 0.);
              flux.set(j,n_phi,1, 0.);
              flux.set(j,n_E2 ,2, 0.);
              flux.set(j,n_phi,2, 0.);
          }
        }
        else
        {
              flux.set(j,n_E1 ,1, 0.);
              flux.set(j,n_E2 ,2, 0.);
        }
    }
}

void MaxwellFluxFunc1(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    // Parameters
    double cs_light_squared = maxwellParams.get_cs_light_squared();

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
        flux.set(j,n_B2 , -E3 );
        flux.set(j,n_B3 ,  E2 );
        flux.set(j,n_E2 ,  cs_light_squared*B3 );
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
          flux.set(j,n_B1 ,  0. );
        }
        //
        if(n_phi<=Q.getsize(2))
        {
          if(maxwellParams.get_clean_E_field())
          {
              const double& phi   = Q.get(j,n_phi);
              // components of 1-component of flux function involving phi
              flux.set(j,n_E1 ,  phi );
              flux.set(j,n_phi,  cp_speed_squared*E1);
          }
          else
          {
              flux.set(j,n_E1 , 0.);
              flux.set(j,n_phi, 0.);
          }
        }
        else
        {
              flux.set(j,n_E1 , 0.);
        }
    }
}

void MaxwellFluxFunc2(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
{
    // Parameters
    double cs_light_squared = maxwellParams.get_cs_light_squared();

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
    
        // 2-component of flux function
        //
        flux.set(j,n_B1 ,  E3 );
        flux.set(j,n_B3 , -E1 );
        flux.set(j,n_E1 , -cs_light_squared*B3 );
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
          flux.set(j,n_B2 ,  0. );
        }
        //
        if(n_phi<=Q.getsize(2))
        {
          if(maxwellParams.get_clean_E_field())
          {
              const double& phi   = Q.get(j,n_phi);
              // components of 2-component of flux function involving phi
              flux.set(j,n_E2 ,  phi );
              flux.set(j,n_phi,  cp_speed_squared*E2);
          }
          else
          {
              flux.set(j,n_E2 , 0.);
              flux.set(j,n_phi, 0.);
          }
        }
        else
        {
              flux.set(j,n_E2 , 0.);
        }
    }
}

void ProjectLeftEig_Maxwell(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& Q, dTensor2& W)
{
    const int n_B3   = n_offset + M_B3 ;
    const int n_E3   = n_offset + M_E3 ;
    const int n_psi  = n_offset + M_psi;
    const int n_phi  = n_offset + M_phi;
    //
    int n_B1;
    int n_B2;
    int n_E1;
    int n_E2;
    if (ixy==1)
    {
        n_B1  = n_offset + M_B1;
        n_B2  = n_offset + M_B2;
        n_E1  = n_offset + M_E1;
        n_E2  = n_offset + M_E2;
    }
    else
    {
        assert(ixy==2);
        n_B1  = n_offset + M_B2 ;
        n_B2  = n_offset + M_B1 ;
        n_E1  = n_offset + M_E2 ;
        n_E2  = n_offset + M_E1 ;
    }

    const int m_B2_E3  = n_offset + M_m_B2_E3 ;
    const int m_B3_E2  = n_offset + M_m_B3_E2 ;
    //
    const int p_B3_E2  = n_offset + M_p_B3_E2 ;
    const int p_B2_E3  = n_offset + M_p_B2_E3 ;
    //
    const int m_B1_psi = n_offset + M_m_B1_psi;
    const int p_B1_psi = n_offset + M_p_B1_psi;
    //
    const int m_E1_phi = n_offset + M_m_E1_phi;
    const int p_E1_phi = n_offset + M_p_E1_phi;

    const double cs_light = maxwellParams.get_cs_light();
    const double cp_speed = maxwellParams.get_cp_speed();
    const double cs2_inv = 1./(2.*cs_light);
    const double cp2_inv = 1./(2.*cp_speed);
    for (int k=1; k<=W.getsize(2); k++)
    {
        W.set(m_B2_E3 ,k, (cs_light*Q.get(n_B2 ,k) - Q.get(n_E3 ,k))*cs2_inv ); // eig -c
        W.set(p_B2_E3 ,k, (cs_light*Q.get(n_B2 ,k) + Q.get(n_E3 ,k))*cs2_inv ); // eig +c
        //
        W.set(m_B3_E2 ,k, (cs_light*Q.get(n_B3 ,k) - Q.get(n_E2 ,k))*cs2_inv ); // eig -c
        W.set(p_B3_E2 ,k, (cs_light*Q.get(n_B3 ,k) + Q.get(n_E2 ,k))*cs2_inv ); // eig +c
        //
        if(n_psi<=Q.getsize(1))
        {
          W.set(m_B1_psi,k, (cp_speed*Q.get(n_B1 ,k) - Q.get(n_psi,k))*cp2_inv ); // eig -cp_speed
          W.set(p_B1_psi,k, (cp_speed*Q.get(n_B1 ,k) + Q.get(n_psi,k))*cp2_inv ); // eig +cp_speed
        }
        else
        {
          // we are not using psi so we are not using phi either, so we just copy the data
          W.set(m_B1_psi,k, Q.get(n_B1,k));
          W.set(p_B1_psi,k, Q.get(n_E1,k));
          return;
        }
        //
        if(n_phi<=Q.getsize(1))
        {
          W.set(m_E1_phi,k, (cp_speed*Q.get(n_phi,k) - Q.get(n_E1 ,k))*cp2_inv ); // eig -cp_speed
          W.set(p_E1_phi,k, (cp_speed*Q.get(n_phi,k) + Q.get(n_E1 ,k))*cp2_inv ); // eig +cp_speed
        }
        else
        {
          // not using phi so just copy the data
          W.set(m_E1_phi,k, Q.get(n_E1,k));
        }
    }
}

void ProjectRightEig_Maxwell(int ixy, int n_offset,
    const dTensor1& Q_ave, const dTensor2& W, dTensor2& Q)
{
    const int n_B3   = n_offset + M_B3 ;
    const int n_E3   = n_offset + M_E3 ;
    const int n_psi  = n_offset + M_psi;
    const int n_phi  = n_offset + M_phi;
    //
    int n_B1;
    int n_B2;
    int n_E1;
    int n_E2;
    if (ixy==1)
    {
        n_B1  = n_offset + M_B1;
        n_B2  = n_offset + M_B2;
        n_E1  = n_offset + M_E1;
        n_E2  = n_offset + M_E2;
    }
    else
    {
        assert(ixy==2);
        n_B1  = n_offset + M_B2 ;
        n_B2  = n_offset + M_B1 ;
        n_E1  = n_offset + M_E2 ;
        n_E2  = n_offset + M_E1 ;
    }

    const int m_B2_E3  = n_offset + M_m_B2_E3 ;
    const int p_B2_E3  = n_offset + M_p_B2_E3 ;
    //
    const int m_B3_E2  = n_offset + M_m_B3_E2 ;
    const int p_B3_E2  = n_offset + M_p_B3_E2 ;
    //
    const int m_B1_psi = n_offset + M_m_B1_psi;
    const int p_B1_psi = n_offset + M_p_B1_psi;
    //
    const int m_E1_phi = n_offset + M_m_E1_phi;
    const int p_E1_phi = n_offset + M_p_E1_phi;

    const double cs_light = maxwellParams.get_cs_light();
    const double cp_speed = maxwellParams.get_cp_speed();
    for (int k=1; k<=Q.getsize(2); k++)
    {
        Q.set(n_B2 ,k,            W.get(m_B2_E3 ,k) + W.get(p_B2_E3, k) );
        Q.set(n_E3 ,k, cs_light*( W.get(p_B2_E3 ,k) - W.get(m_B2_E3, k)) );
        //
        Q.set(n_B3 ,k,            W.get(p_B3_E2 ,k) + W.get(m_B3_E2, k) );
        Q.set(n_E2 ,k, cs_light*( W.get(p_B3_E2 ,k) - W.get(m_B3_E2, k)) );
        //
        if(n_psi<=Q.getsize(1))
        {
          Q.set(n_B1 ,k,            W.get(p_B1_psi,k) + W.get(m_B1_psi,k) );
          Q.set(n_psi,k, cp_speed*( W.get(p_B1_psi,k) - W.get(m_B1_psi,k)) );
        }
        else
        {
          // not using psi or phi. copy the data back.
          Q.set(n_B1,k, W.get(m_B1_psi,k));
          Q.set(n_E1,k, W.get(p_B1_psi,k));
          return;
        }
        //
        if(n_phi<=Q.getsize(1))
        {
          Q.set(n_E1 ,k, cp_speed*( W.get(p_E1_phi,k) - W.get(m_E1_phi,k)) );
          Q.set(n_phi,k,            W.get(p_E1_phi,k) + W.get(m_E1_phi,k));
        }
        else
        {
          // not using phi so copy the data back
          Q.set(n_E1,k, W.get(m_E1_phi,k));
        }
    }
}

void SymPair_MaxwellFluxFunc(
    int n_offset,
    const dTensor2& Q,
    dTensor3& flux)
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

void SymPair_MaxwellFluxFunc1(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
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

void SymPair_MaxwellFluxFunc2(
    int n_offset,
    const dTensor2& Q,
    dTensor2& flux)
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
    const double cp_speed = maxwellParams.get_cp_speed();
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
    const double cp_speed = maxwellParams.get_cp_speed();
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
