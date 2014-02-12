///@file  lib/1d/dimdefs.h
///@brief Defines #NDIMS, the dimensionality of the problem.
///
///- Use: The dimensions at which <tt>dTensorBC</tt> are padded with <tt>mbc</tt> ghost cells on two ends
///     are {1,..., #NDIMS}
///- This file precedes lib/dimdefs.h in the compile commands, thus overriding the definition #NDIMS=0 in that file.
#define NDIMS 1

