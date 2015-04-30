#ifndef _CENTRAL_DIFFERENCES_H2D_H__
#define _CENTRAL_DIFFERENCES_H2D_H__

#include "diffMatrix5.h"
#include "diffMatrix7.h"
#include "diffMatrix9.h"
#include "diffMatrix11.h"
#include "Ix5.h"
#include "Ix7.h"
#include "Ix9.h"
#include "Ix11.h"
#include "Iy5.h"
#include "Iy7.h"
#include "Iy9.h"
#include "Iy11.h"

// -- Reconstructions based on linear weights -- //
//
// There's no need to call these directly.  In place of extra overhead for
// checking which method to use inside a for loop, instead use
// GetCentralDifferences one time.
void CentralDifferences2D5 ( double dx, double dy, const dTensor3& f, dTensor3& fderivs   );
void CentralDifferences2D7 ( double dx, double dy, const dTensor3& f, dTensor3& fderivs   );
void CentralDifferences2D9 ( double dx, double dy, const dTensor3& f, dTensor3& fderivs   );
void CentralDifferences2D11( double dx, double dy, const dTensor3& f, dTensor3& fderivs   );

// Wrapper function that provides access to each of the above through looking
// at the global variable wenoParams.
//void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
typedef void (*central_differences_t2d)(double, double, const dTensor3&, dTensor3&);
central_differences_t2d GetCentralDifferences2D();

// Matrices used for computing derivatives.  
//
//These matrices assume that points are organized as follows....

/* (for example for 9 points...)
(x,y)=[[-4*dx, -4*dy], [-3*dx, -4*dy], [-2*dx, -4*dy], [-dx, -4*dy], [0, -4*dy], [dx, -4*dy], [2*dx, -4*dy], [3*dx, -4*dy], [4*dx, -4*dy], [-4*dx, -3*dy], [-3*dx, -3*dy], [-2*dx, -3*dy], [-dx, -3*dy], [0, -3*dy], [dx, -3*dy], [2*dx, -3*dy], [3*dx, -3*dy], [4*dx, -3*dy], [-4*dx, -2*dy], [-3*dx, -2*dy], [-2*dx, -2*dy], [-dx, -2*dy], [0, -2*dy], [dx, -2*dy], [2*dx, -2*dy], [3*dx, -2*dy], [4*dx, -2*dy], [-4*dx, -dy], [-3*dx, -dy], [-2*dx, -dy], [-dx, -dy], [0, -dy], [dx, -dy], [2*dx, -dy], [3*dx, -dy], [4*dx, -dy], [-4*dx, 0], [-3*dx, 0], [-2*dx, 0], [-dx, 0], [0, 0], [dx, 0], [2*dx, 0], [3*dx, 0], [4*dx, 0], [-4*dx, dy], [-3*dx, dy], [-2*dx, dy], [-dx, dy], [0, dy], [dx, dy], [2*dx, dy], [3*dx, dy], [4*dx, dy], [-4*dx, 2*dy], [-3*dx, 2*dy], [-2*dx, 2*dy], [-dx, 2*dy], [0, 2*dy], [dx, 2*dy], [2*dx, 2*dy], [3*dx, 2*dy], [4*dx, 2*dy], [-4*dx, 3*dy], [-3*dx, 3*dy], [-2*dx, 3*dy], [-dx, 3*dy], [0, 3*dy], [dx, 3*dy], [2*dx, 3*dy], [3*dx, 3*dy], [4*dx, 3*dy], [-4*dx, 4*dy], [-3*dx, 4*dy], [-2*dx, 4*dy], [-dx, 4*dy], [0, 4*dy], [dx, 4*dy], [2*dx, 4*dy], [3*dx, 4*dy], [4*dx, 4*dy]]
*/

//
/*
[0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4 0 1 2 3 4]
[0 0 0 0 0 1 1 1 1 1 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4]
[f f_x f_xx f_3x f_4x f_y f_xy f_2xy f_3xy f_4xy f_2y f_x2y f_2x2y f_3x2y f_4x2y f_3y f_x3y f_2x3y f_3x3y f_4x3y f_4y f_x4y f_2x4y f_3x4y f_4x4y]
*/




// Old version of the same matrix
//  const double deriv_matrix5[4][5] = {
//     { 8.333333333333333e-02, -6.666666666666666e-01,  0, 6.666666666666666e-01,  -8.333333333333333e-02 },   // first-deriv
//     {-8.333333333333333e-02,  1.3333333333333333, -2.5, 1.3333333333333333, -8.333333333333333e-02},  // second-deriv
//     {-0.5,     1.0,    0.0,  -1.0,    0.5},     // third-deriv
//     {1.0,    -4.0,    6.0,  -4.0,    1.0}       // fourth-deriv
//  };

#endif
