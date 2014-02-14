#ifndef _RKINFO_H_
#define _RKINFO_H_

// Runge-Kutta information

///@brief A struct which stores Runge-Kutta information
///
///How it is used in the code: (in #FinSolveRK(...) ):
/**
@code
...
RKinfo rk;
SetRKinfo(time_order, rk);
...
switch(time_order){
	case <order>:
	  rk.mstage=1;
	  BeforeStep(...);
	  ConstructL(...);
	  UpdateSoln(...);  // some of parameters are [rk.alpha1 | rk.alpha2 | rk.gamma | rk.delta] -> get(something index(-ices) involving rk.mstage
	  AfterStep(...);
	  ...
	  rk.mstage=<order>
	  BeforeStep(...);
	  ConstructL(...);
	  UpdateSoln(...);  // some of parameters are [rk.alpha1 | rk.alpha2 | rk.gamma | rk.delta] -> get(something index(-ices) involving rk.mstage
	  AfterStep(...);
	  break;
	...
}
...
DeleteRKInfo(rk);  // Releases memory for arrays in rk.
@endcode
*/
///Note the following:
///- #mstage only serves as a local iteration-like variable, thus should be removed;
///- #RKinfo is agnostic about the order of the method.
///  The first argument to SetRKinfo(int, RKinfo&) only serves 
///      as the index to an element in a set of RK methods.
///  There exists an (implicit) contract between the information (number of stages, and coefficients) set by SetRKinfo(...),
///                  and the RK method code in FinSolveRK(...).
///  This contract can be kept only by *active*, *simultaneous* maintainence 
///       of the code of both FinSolveRK(...) and SetRKinfo(...).
///
///@todo %RKinfo should contain the *algorithm* for RK methods, not only the coefficients, 
///      thus,
///      - see if this can be redesigned, while still keeping the flexibility in tuning FinSolveRK(...),
///      and an acceptable performance.
struct RKinfo
{
    ///@brief (Looks like) Index of current stage when doing Runge-Kutta
	///
	///@todo Remove it.  This is used only in FinSolveRK(...),
	///      and only serves as a local iteration variable -- why then do we store it here?
	int mstage;
	///@brief Number of stages in Runge-Kutta method
	int num_stages;
	dTensor1* alpha1;
	dTensor1* alpha2;
	dTensor1* beta;

	// These two are needed for 5th order Stepping
	dTensor2* gamma;
	dTensor1* delta;

};

#endif
