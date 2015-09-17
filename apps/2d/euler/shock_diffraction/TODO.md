Notes on why this application will not compile
==============================================

The main library function, ConstructL was written with the intention that
it would not call SetBndValues, and therefore whatever values of Q get passed
in, they should be constructed as 'const' values.  However, this application
requires a total of two calls to SetBndValues:

    SetBndValuesX and SetBndValuesY.

This means that the overwritten application ConstructL will not work, because
we would like to call SetBndValuesX before defining f_x, and we would like to
call SetBndValuesY before calling g_y.

I see two obvious solutions, neither of which are satisfactory:

    1. Replace the library routine with non-const reference to StateVars Q.
    2. Copy every time-integrator in the library to this application.

Until a reasonable solution is reached, this application will not run.

Update (Qi Tang):
All the functions (except ConstructLFL_BCx2) that uses SetBndValuesX and 
SetBndValuesY have been rewritten locally. Currently this example can run. 
But it may require some optimazations if there is a similar test problem 
in the future.
