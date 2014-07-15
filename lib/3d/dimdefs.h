// This module defines NDIMS, which is used in the tensor class to define an
// appropriate number of ghost cells for, in this case, the first three indices.
// That is, valid arguments for q are q( 1-mbc:mx+mbc, ... ) in place of 
// q(1:mx, ... ).
#define NDIMS 3
