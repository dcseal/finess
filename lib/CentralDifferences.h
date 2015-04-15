#ifndef _WENORECONSTRUCT_H_
#define _WENORECONSTRUCT_H_

// -- Reconstructions based on linear weights -- //
void CentralDifferences5( const dTensor2& g, dTensor2& g_reconst );
void CentralDifferences7( const dTensor2& g, dTensor2& g_reconst );
// void CentralDifferences9( const dTensor2& g, dTensor2& g_reconst );

// Wrapper function that provides access to each of the above through looking
// at the global variable wenoParams.
//void (*GetWenoReconstruct())(const dTensor2& g, dTensor2& g_reconst);
typedef void (*central_differences_t)(const dTensor2&, dTensor2&);
central_differences_t GetCentralDifferences();

#endif
