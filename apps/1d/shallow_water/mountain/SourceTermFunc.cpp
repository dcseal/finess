#include "stdio.h"      // Remove this after getting rid of print statement
#include "tensors.h"

// This is a user-supplied routine that sets the source term
void SourceTermFunc(const dTensor1& xpts, 
            const dTensor2& qvals, 
            const dTensor2& auxvals,
                    dTensor2& fvals)
{

    const int numpts=xpts.getsize();
    const int meqn=qvals.getsize(2);
    for(int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i);
        
        const double h = qvals.get(i,1);
        const double u = qvals.get(i,2)/h;
        
// TODO - this is not consistent with what is in AuxFunc!
printf("WARNING: accessing potentially awful memory here\n");
        const double bot  = auxvals.get(i,1);
        const double dbdx = auxvals.get(i,2);
        
        fvals.set(i,1,  0.0 );
        fvals.set(i,2, -h*dbdx );

    }

}
