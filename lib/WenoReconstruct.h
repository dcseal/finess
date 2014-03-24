#ifndef _WENORECONSTRUCT_H_
#define _WENORECONSTRUCT_H_

void WenoReconstruct_JS5( const dTensor2& g, dTensor2& diff_g );
void WenoReconstruct_JS7( const dTensor2& g, dTensor2& diff_g );
void WenoReconstruct_JS9( const dTensor2& g, dTensor2& diff_g );
void WenoReconstruct(const dTensor2& g, dTensor2& diff_g);

#endif