#include "stdio.h"
#include "dog_math.h"
#include "constants.h"
#include "tensors.h"

void ApplyLimiter(const dTensor2& node, const dTensorBC2& aux, dTensorBC2& q, 
        void (*ProjectRightEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&),
        void (*ProjectLeftEig)(const dTensor1&,const dTensor1&,const dTensor2&,dTensor2&))
{

    printf("WARNING: no limiter has been implemented\n");

}
