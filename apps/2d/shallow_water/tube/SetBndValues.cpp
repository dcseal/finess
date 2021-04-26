#include "dogdefs.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      REFLECTIVE BOUNDARY CONDITIONS
//
void SetBndValues( StateVars& Q )
{

  dTensorBC3&  q  = Q.ref_q  ();
  dTensorBC3& aux = Q.ref_aux();
  double t        = Q.get_t  ();

  const int mx   = q.getsize(1);
  const int my   = q.getsize(2);
  const int meqn = q.getsize(3);
  const int mbc  = q.getmbc();
  const int maux = aux.getsize(3);

  // -----------------------
  // BOUNDARY CONDITION LOOP
  // -----------------------

  // ***********************************************
  // LEFT BOUNDARY
  // ***********************************************
  for(int i=0; i>=(1-mbc); i--)
  {
    for(int j=1; j<=my; j++)
    {
      // q values
      for(int m=1; m<=meqn; m++)
      {
                double tmp = q.get(1,j,m);
//              double tmp = q.get(i+mx,j,m); // periodic
                q.set(i,j,m, tmp );
      }

      // Flip the momemntum
      q.set(i,j,2, -q.get(i,j,2) );

      // aux values
      for(int m=1; m<=maux; m++)
      {
        double tmp = aux.get(i+mx,j,m);
        aux.set(i,j,m, tmp );
      }
    }
  }
  // ***********************************************



  // ***********************************************
  // RIGHT BOUNDARY
  // ***********************************************
  for(int i=(mx+1); i<=(mx+mbc); i++)
  {
    for(int j=1; j<=my; j++)
    {
      // q values
      for (int m=1; m<=meqn; m++)
      {
          double tmp = q.get(mx,j,m);
          // double tmp = q.get(i-mx,j,m); // periodic
          q.set(i,j,m, tmp );
      }

      // Flip the momemntum
      q.set(i,j,2, -q.get(i,j,2) );

      // aux values
      for(int m=1; m<=maux; m++)
      {
        double tmp = aux.get(i-mx,j,m);
        aux.set(i,j,m, tmp );
      }
    }
  }
  // ***********************************************



  // ***********************************************
  // BOTTOM BOUNDARY
  // ***********************************************
  for(int j=0; j>=(1-mbc); j--)
  {
    for(int i=1; i<=mx; i++)
    {
      // q values
      for (int m=1; m<=meqn; m++)
        {
          double tmp = q.get(i, 1, m);
          q.set(i, j, m, tmp );
        }

      // Flip the momemntum
      q.set(i,j,3, -q.get(i,j,3) );

      // aux values
      for(int m=1; m<=maux; m++)
      {
        double tmp = aux.get(i,j+my,m);
        aux.set(i,j,m, tmp );
      }
    }
  }
  // ***********************************************



  // ***********************************************
  // TOP BOUNDARY
  // ***********************************************
  for(int j=(my+1); j<=(my+mbc); j++)
  {
    for(int i=1; i<=mx; i++)
    {

      // q values
      for(int m=1; m<=meqn; m++)
      {
        double tmp = q.get(i,my,m);
        q.set(i,j,m, tmp );
      }

      // Flip the momemntum
      q.set(i,j,3, -q.get(i,j,3) );

      // aux values
      for(int m=1; m<=maux; m++)
      {
        double tmp = aux.get(i,j-my,m);
        aux.set(i,j,m, tmp );
      }
    }
  }
  // ***********************************************

  // ***********************************************
  // BOTTOM LEFT CORNER
  // ***********************************************
  for(int i=1; i<=mbc; i++)
    for(int j=1; j<=mbc; j++)
    {
      for(int m=1; m<=meqn; m++)
      {
        q.set(1-i,1-j,m, q.get(mx+1-i,my+1-j,m) );
      }
      for(int m=1; m<=maux; m++)
      {
        aux.set(1-i,1-j,m, aux.get(mx+1-i,my+1-j,m) );
      }
    }
  // ***********************************************


  // ***********************************************
  // BOTTOM RIGHT CORNER
  // ***********************************************
  for(int i=1; i<=mbc; i++)
    for(int j=1; j<=mbc; j++)
    {
      for(int m=1; m<=meqn; m++)
      {
        q.set(mx+i,1-j,m, q.get(i,my+1-j,m) );
      }
      for(int m=1; m<=maux; m++)
      {
        aux.set(mx+i,1-j,m, aux.get(i,my+1-j,m) );
      }
    }
  // ***********************************************


  // ***********************************************
  // TOP RIGHT CORNER
  // ***********************************************
  for(int i=1; i<=mbc; i++)
    for(int j=1; j<=mbc; j++)
    {
      for(int m=1; m<=meqn; m++)
      {
        q.set(mx+i,my+j,m, q.get(i,j,m) );
      }
      for(int m=1; m<=maux; m++)
      {
        aux.set(mx+i,my+j,m, aux.get(i,j,m) );
      }
    }
  // ***********************************************


  // ***********************************************
  // TOP LEFT CORNER
  // ***********************************************
  for(int i=1; i<=mbc; i++)
    for(int j=1; j<=mbc; j++)
    {
      for(int m=1; m<=meqn; m++)
      {
        q.set(1-i,my+j,m, q.get(mx+1-i,j,m) );
      }
      for(int m=1; m<=maux; m++)
      {
        aux.set(1-i,my+j,m, aux.get(mx+1-i,j,m) );
      }
    }
  // ***********************************************

}

// Wrappers for main Euler library
void SetBndValuesX( StateVars& Q )
{ SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }
