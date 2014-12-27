#include <cmath>
#include "dogdefs.h"
#include "constants.h"

// This is a user-supplied routine that sets the
// initial conditions at all the points "xpts"
//
void QinitFunc(const dTensor2& xpts, dTensor2& qvals)
{
    using std::sin;
    const int numpts=xpts.getsize(1);
    for (int i=1; i<=numpts; i++)
    {
        const double x = xpts.get(i,1);
        const double y = xpts.get(i,2);
        const double z = xpts.get(i,3);
        qvals.set(i, 1, 0.25 + sin(M_PI * x) * sin(M_PI * y) * sin(M_PI * z));

    }
}
