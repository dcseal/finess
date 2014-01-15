#include "LegendrePolys1d.h"

// function that evaluates legendre polynomials and all their derivatives
// at the list of points given by spts
// 
// Pseudo duplicate
void evaluateLegendrePolys(const dTensor1& spts, dTensor2& phi )
{

    const int numpts = spts.getsize();

    for (int m=1; m<=numpts; m++)
    {
        // Legendre basis functions at each grid point
        phi.set( m,1, 1.0 );
        phi.set( m,2, sq3*spts.get(m) );
        phi.set( m,3, 0.5*sq5*( 3.0*pow(spts.get(m),2) - 1.0 ) );
        phi.set( m,4, 0.5*sq7*spts.get(m)
                *(5.0*pow(spts.get(m),2) - 3.0) );
        phi.set( m,5, (105.0/8.0)*pow(spts.get(m),4) 
                - (45.0/4.0)*pow(spts.get(m),2) + (9.0/8.0) );
        phi.set( m, 6, (63.0/8.0)*sq11 * ( pow(spts.get(m),5)  -
                    (10.0/9.0)*(pow(spts.get(m),3) ) + (5.0/21.0)*spts.get(m) ) );


    }
}

// function that evaluates legendre polynomials and all their derivatives
// at the list of points given by spts
void evaluateLegendrePolys( 
        const double dx, const dTensor1& spts, dTensor2& phi, dTensor2& phi_x)
{

    // Regular Basis functions
    evaluateLegendrePolys(spts, phi );

    // 1st Derivative of Basis functions
    const int numpts = spts.getsize();
    for (int m=1; m<=numpts; m++)
    {
        phi_x.set( m,1, 0.0 );
        phi_x.set( m,2, 2.0*sq3/dx );
        phi_x.set( m,3, 6.0*sq5*spts.get(m)/dx );
        phi_x.set( m,4, 3.0*sq7*(5.0*pow(spts.get(m),2)-1.0)/dx );
        phi_x.set( m,5, 15.0*spts.get(m)*
                (7.0*pow(spts.get(m),2)-3.0)/dx );
    }
}


void setGaussLobattoPoints1d(dTensor1& x1d, dTensor1& wgt)
{

    switch( x1d.getsize() )
    {
        case 2:
            x1d.set(1, -1.0 );
            x1d.set(2,  1.0 );

            wgt.set(1, 1.0 );
            wgt.set(2, 1.0 );
            break;

        case 3:
            x1d.set(1,-1.0);
            x1d.set(2, 0.0);
            x1d.set(3, 1.0);

            wgt.set(1, 1.0/3.0 );
            wgt.set(2, 4.0/3.0 );
            wgt.set(3, 1.0/3.0 );
            break;

        case 4:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.447213595499958 );
            x1d.set(3,  0.447213595499958 );
            x1d.set(4, 1.0 );


            wgt.set(1, 0.166666666666667 );
            wgt.set(2, 0.833333333333333 );
            wgt.set(3, wgt.get(2) );
            wgt.set(4, wgt.get(1) );
            break;

        case 5:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.654653670707977 );
            x1d.set(3,  0.0 );
            x1d.set(4, -x1d.get(2) );
            x1d.set(5, -x1d.get(1) );

            wgt.set(1, 0.10 );
            wgt.set(2, 0.544444444444444 );
            wgt.set(3, 0.711111111111111 );
            wgt.set(4, wgt.get(2) );
            wgt.set(5, wgt.get(1) );
            break;

        case 6:
            x1d.set(1, -1.0 );
            x1d.set(2, -0.765055323929465 );
            x1d.set(3, -0.285231516480645 );
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );

            wgt.set(1, 0.066666666666667 );
            wgt.set(2, 0.378474956297847 ); 
            wgt.set(3, 0.554858377035486 );
            wgt.set(4, wgt.get(3) );
            wgt.set(5, wgt.get(2) );
            wgt.set(6, wgt.get(1) );
            break;

        case 7:      
            x1d.set(1,  -1);
            x1d.set(2,  -8.302238962785670e-01);
            x1d.set(3,  -4.688487934707142e-01);
            x1d.set(4,   0);
            x1d.set(5, -x1d.get(3) );
            x1d.set(6, -x1d.get(2) );
            x1d.set(7, -x1d.get(1) );

            wgt.set(1, 0.047619047619048 );
            wgt.set(2, 0.276826047361566 );
            wgt.set(3, 0.431745381209863 );
            wgt.set(4, 0.487619047619048 );
            wgt.set(5, wgt.get(3) );
            wgt.set(6, wgt.get(2) );
            wgt.set(7, wgt.get(1) );

            break;

        default:
            exit(1);
    }
}

//////////////////////////////////////////////////////////////////////////
// Set quadrature weights and points (standard before a transformation)
//////////////////////////////////////////////////////////////////////////
void setGaussLegendrePoints1d(dTensor1& x1d, dTensor1& w1d)
{

    switch ( w1d.getsize() )
    {
        case 1:
            w1d.set(1, 2.0e0 );
            x1d.set(1, 0.0e0 );

            break;

        case 2:
            w1d.set(1,   1.0 );
            w1d.set(2,   1.0 );

            x1d.set(1, -1.0/sq3 );
            x1d.set(2,  1.0/sq3 );

            break;

        case 3:
            w1d.set(1, 5.0e0/9.0e0 );
            w1d.set(2, 8.0e0/9.0e0 );
            w1d.set(3, 5.0e0/9.0e0 );

            x1d.set(1, -sq3/sq5 );
            x1d.set(2,  0.0e0 );
            x1d.set(3,  sq3/sq5 );

            break;

        case 4:
            w1d.set(1, (18.0 - sqrt(30.0))/36.0 );
            w1d.set(2, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(3, (18.0 + sqrt(30.0))/36.0 );
            w1d.set(4, (18.0 - sqrt(30.0))/36.0 );

            x1d.set(1, -sqrt(3.0+sqrt(4.8))/sq7 );
            x1d.set(2, -sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(3,  sqrt(3.0-sqrt(4.8))/sq7 );
            x1d.set(4,  sqrt(3.0+sqrt(4.8))/sq7 );

            break;

        case 5:
            w1d.set(1, (322.0 - 13.0*sqrt(70.0))/900.0 );
            w1d.set(2, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(3,  128.0/225.0 );
            w1d.set(4, (322.0 + 13.0*sqrt(70.0))/900.0 );
            w1d.set(5, (322.0 - 13.0*sqrt(70.0))/900.0 );

            x1d.set(1, -sqrt(5.0 + sqrt(40.0/7.0))/3.0 );
            x1d.set(2, -sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(3,  0.0 );
            x1d.set(4,  sqrt(5.0 - sqrt(40.0/7.0))/3.0 );
            x1d.set(5,  sqrt(5.0 + sqrt(40.0/7.0))/3.0 );

            break;

        case 6:

            w1d.set(1,  0.171324492379170);
            w1d.set(2,  0.360761573048139);
            w1d.set(3,  0.467913934572691);
            w1d.set(4,  0.467913934572691);
            w1d.set(5,  0.360761573048139);
            w1d.set(6,  0.171324492379170);


            x1d.set(1,  0.932469514203152);
            x1d.set(2,  0.661209386466265);
            x1d.set(3,  0.238619186083197);
            x1d.set(4, -x1d.get(3) );
            x1d.set(5, -x1d.get(2) );
            x1d.set(6, -x1d.get(1) );

            break;


        case 20:

            x1d.set(1,  0.993128599185095);
            x1d.set(2,  0.963971927277914);
            x1d.set(3,  0.912234428251326);
            x1d.set(4,  0.839116971822219);
            x1d.set(5,  0.746331906460151);
            x1d.set(6,  0.636053680726515);
            x1d.set(7,  0.510867001950827);
            x1d.set(8,  0.373706088715420);
            x1d.set(9,  0.227785851141645);
            x1d.set(10, 0.076526521133497);
            for( int n=1; n <= 10; n++ ) 
            { x1d.set(10+n, -x1d.get(11-n) ); }

            w1d.set(1,  0.017614007139152 );
            w1d.set(2,  0.040601429800387);
            w1d.set(3,  0.062672048334109);
            w1d.set(4,  0.083276741576705);
            w1d.set(5,  0.101930119817240);
            w1d.set(6,  0.118194531961518);
            w1d.set(7,  0.131688638449177);
            w1d.set(8,  0.142096109318382);
            w1d.set(9,  0.149172986472604);
            w1d.set(10,  0.152753387130726);
            for( int n=1; n <= 10; n++ ) 
            { w1d.set(10+n, w1d.get(11-n) ); }

            break;

        case 50:
            x1d.set(1,  9.988664044200710e-01);
            x1d.set(2,  9.940319694320907e-01);
            x1d.set(3,  9.853540840480060e-01);
            x1d.set(4,  9.728643851066920e-01);
            x1d.set(5,  9.566109552428079e-01);
            x1d.set(6,  9.366566189448780e-01);
            x1d.set(7,  9.130785566557917e-01);
            x1d.set(8,  8.859679795236131e-01);
            x1d.set(9,  8.554297694299460e-01);
            x1d.set(10, 8.215820708593360e-01);
            x1d.set(11, 7.845558329003994e-01);
            x1d.set(12, 7.444943022260686e-01);
            x1d.set(13, 7.015524687068222e-01);
            x1d.set(14, 6.558964656854394e-01);
            x1d.set(15, 6.077029271849503e-01);
            x1d.set(16, 5.571583045146502e-01);
            x1d.set(17, 5.044581449074643e-01);
            x1d.set(18, 4.498063349740388e-01);
            x1d.set(19, 3.934143118975651e-01);
            x1d.set(20, 3.355002454194373e-01);
            x1d.set(21, 2.762881937795321e-01);
            x1d.set(22, 2.160072368760417e-01);
            x1d.set(23, 1.548905899981459e-01);
            x1d.set(24, 9.317470156008612e-02);
            x1d.set(25, 3.109833832718883e-02);
            for( int n=1; n <= 25; n++ ) 
            { x1d.set(25+n, -x1d.get(26-n) ); }

            w1d.set(1,  2.908622553155257e-03);
            w1d.set(2,  6.759799195745480e-03);
            w1d.set(3,  1.059054838365105e-02);
            w1d.set(4,  1.438082276148557e-02);
            w1d.set(5,  1.811556071348930e-02);
            w1d.set(6,  2.178024317012479e-02); 
            w1d.set(7,  2.536067357001242e-02); 
            w1d.set(8,  2.884299358053520e-02); 
            w1d.set(9,  3.221372822357796e-02); 
            w1d.set(10, 3.545983561514608e-02); 
            w1d.set(11, 3.856875661258761e-02); 
            w1d.set(12, 4.152846309014772e-02); 
            w1d.set(13, 4.432750433880325e-02); 
            w1d.set(14, 4.695505130394845e-02); 
            w1d.set(15, 4.940093844946636e-02); 
            w1d.set(16, 5.165570306958114e-02); 
            w1d.set(17, 5.371062188899628e-02); 
            w1d.set(18, 5.555774480621253e-02); 
            w1d.set(19, 5.718992564772843e-02); 
            w1d.set(20, 5.860084981322242e-02); 
            w1d.set(21, 5.978505870426547e-02); 
            w1d.set(22, 6.073797084177022e-02); 
            w1d.set(23, 6.145589959031677e-02); 
            w1d.set(24, 6.193606742068321e-02); 
            w1d.set(25, 6.217661665534725e-02); 
            for( int n=1; n <= 25; n++ ) 
            { w1d.set(25+n, w1d.get(26-n) ); }

            break;

    }

}


