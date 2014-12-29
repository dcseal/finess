#include "dogdefs.h"
#include "StateVars.h"

// This is a user-supplied routine that sets the the boundary conditions
//
//      ZEROTH ORDER EXTRAPOLATION BOUNDARY CONDITIONS
//
void SetBndValues( StateVars& Q )
{

    dTensorBC3&  q  = Q.ref_q  ();
    dTensorBC3& aux = Q.ref_aux();
    double t        = Q.get_t  ();

    int i,j,m;
    double tmp;
    int mx   = q.getsize(1);
    int my   = q.getsize(2);
    int meqn = q.getsize(3);
    int mbc  = q.getmbc();
    int maux = aux.getsize(3);

    // -----------------------
    // BOUNDARY CONDITION LOOP
    // -----------------------

    // ***********************************************
    // LEFT BOUNDARY
    // ***********************************************
    for (i=0; i>=(1-mbc); i--)
    {
        for (j=1; j<=my; j++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = 0.0; //q.get(1,j,m,);
                q.set(i,j,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = 0.0; //aux.get(1,j,m);
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************

    // ***********************************************
    // RIGHT BOUNDARY
    // ***********************************************
    for (i=(mx+1); i<=(mx+mbc); i++)
    {
        for (j=1; j<=my; j++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = 0.0; //q.get(mx,j,m);                    
                q.set(i,j,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = 0.0; //aux.get(mx,j,m);                    
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // BOTTOM BOUNDARY
    // ***********************************************
    for (j=0; j>=(1-mbc); j--)
    {
        for (i= 1; i<= mx; i++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp = 0.0; //q.get(i,1,m);                    
                q.set(i,j,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp =  0.0; //aux.get(i,1,m);                    
                aux.set(i,j,m, tmp );
            }
        }
    }
    // ***********************************************



    // ***********************************************
    // TOP BOUNDARY
    // ***********************************************
    for (j=(my+1); j<=(my+mbc); j++)
    {
        for (i= 1; i<= mx; i++)
        {           
            // q values
            for (m=1; m<=meqn; m++)
            {
                tmp =  0.0; //q.get(i,my,m);                    
                q.set(i,j,m, tmp );
            }

            // aux values
            for (m=1; m<=maux; m++)
            {
                tmp = 0.0; //aux.get(i,my,m);                    
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
            q.set(1-i,1-j,m, q.get(1,1,m) );
        }
        for(int m=1; m<=maux; m++)
        {     
            aux.set(1-i,1-j,m, aux.get(1,1,m) );
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
            q.set(mx+i,1-j,m, q.get(mx,1,m) );
        }
        for(int m=1; m<=maux; m++)
        {     
            aux.set(mx+i,1-j,m, aux.get(mx,1,m) );
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
            q.set(mx+i,my+j,m, q.get(mx,my,m) );
        }
        for(int m=1; m<=maux; m++)
        {     
            aux.set(mx+i,my+j,m, aux.get(mx,my,m) );
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
            q.set(1-i,my+j,m, q.get(1,my,m) );
        }
        for(int m=1; m<=maux; m++)
        {     
            aux.set(1-i,my+j,m, aux.get(1,my,m) );
        }
    }
    // ***********************************************


}


void SetBndValuesX( StateVars& Q )
{ SetBndValues( Q ); }

void SetBndValuesY( StateVars& Q )
{ SetBndValues( Q ); }
