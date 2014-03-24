#ifndef _WENORECONSTRUCT_H_
#define _WENORECONSTRUCT_H_

// -- Jiang and Shu reconstructions -- //
void WenoReconstruct_JS5( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_JS7( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_JS9( const dTensor2& g, dTensor2& g_reconst );

// -- WENO-Z reconstructions -- //
void WenoReconstruct_Z5( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_Z7( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_Z9( const dTensor2& g, dTensor2& g_reconst );

// -- Reconstructions based on linear weights -- //
// These are useful for running convergence studies, and should not be used
// for problems with discontinuities.
void WenoReconstruct_FD5( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_FD7( const dTensor2& g, dTensor2& g_reconst );
void WenoReconstruct_FD9( const dTensor2& g, dTensor2& g_reconst );

// Wrapper function that provides access to each of the above through looking
// at the global variable wenoParams.
void WenoReconstruct(const dTensor2& g, dTensor2& g_reconst);

#endif
