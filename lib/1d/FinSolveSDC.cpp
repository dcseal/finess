#include <iostream>
#include <iomanip>
#include "dog_math.h"
#include "stdlib.h"
#include "dogdefs.h"
#include "DogParams.h"
#include "RKinfo.h"
#include "FinSolveRK.h"

using namespace std;

void FinSolveSDC(
    dTensorBC2& aux, dTensorBC2& qold, dTensorBC2& qnew, 
    dTensorBC1& smax,
    double tstart, double tend, int nv,
    double dtv[], const double cflv[], string outputdir)
{

    // This is a bit cleaner than putting all the pre and post steps in here
    void EulerStep(const double& dt, 
            dTensorBC1& smax, dTensorBC2& Lrhs,
            dTensorBC2& aux, dTensorBC2& qin, 
            dTensorBC2& qnew);
    void TimeStepSDC(int method2, 
            double t, 
            double dt, 
            dTensor1& dtvec, 
            dTensor1& tvec);
    void ResInt(double dt, 
            const dTensorBC2& L0, 
            const dTensorBC2& L1, 
            const dTensorBC2& L2, 
            const dTensorBC2& L3, 
            const dTensorBC2& L4, 
            const dTensorBC2& L5,
            dTensorBC3& ILout);

//  void StepSDCRK2(const double& dt, const int method[], const dTensor2& node,
//          dTensorBC1& smax, dTensorBC3& Lrhs, dTensorBC3& Lstar,
//          dTensorBC3& aux, dTensorBC3& qin, dTensorBC3& qstar,
//          dTensorBC3& qnew);
//  void StepSDCdeltaRK2(const double& dt, const int method[], const dTensor2& node,
//          dTensorBC1& smax, dTensorBC3& aux, dTensorBC3& qstar, dTensorBC3& Lstar,
//          dTensorBC3& L1, dTensorBC3& L1new, dTensorBC3& L2, dTensorBC3& q1, 
//          dTensorBC3& q2, int num, dTensorBC4& IL);
    // ------------------------------------------------------------

    double t            = tstart;
    double dt           = dtv[1];   // Start with time step from last frame
    double CFL_max      = cflv[1];  // max   CFL number
    double CFL_target   = cflv[2];  // target CFL number
    double cfl          = 0.0;      // current CFL number
    double dtmin        = dt;       // Counters for max and min time step taken
    double dtmax        = dt;

    const int mx    = qold.getsize(1);
    const int meqn  = qold.getsize(2);
    const int maux  = aux.getsize(2);
    const int mbc   = qnew.getmbc();

    const int torder = dogParams.get_time_order();

    // --------------------------------------------------------------
    // Create helper arrays
    dTensor1   dtvec(torder-1);
    dTensor1    tvec(torder  );
    dTensorBC2    q1(mx,meqn,mbc);
    dTensorBC2    q2(mx,meqn,mbc);
    dTensorBC2    q3(mx,meqn,mbc);
    dTensorBC2    q4(mx,meqn,mbc);
    dTensorBC2    q5(mx,meqn,mbc);
    dTensorBC2    qs(mx,meqn,mbc);

    dTensorBC2    L0(mx,meqn,mbc);
    dTensorBC2    L1(mx,meqn,mbc);
    dTensorBC2    L2(mx,meqn,mbc);
    dTensorBC2    L3(mx,meqn,mbc);
    dTensorBC2    L4(mx,meqn,mbc);
    dTensorBC2    L5(mx,meqn,mbc);
    dTensorBC2    Ls(mx,meqn,mbc);

    dTensorBC2 L1new(mx,meqn,mbc);
    dTensorBC2 L2new(mx,meqn,mbc);
    dTensorBC2 L3new(mx,meqn,mbc);
    dTensorBC2 L4new(mx,meqn,mbc);

    dTensorBC3 IL(mx,meqn,torder-1,mbc);
    // --------------------------------------------------------------

    // -----------------------
    // Construct initial L
    // -----------------------
    if( torder > 2 )
    { 
        // create current time step vector
        TimeStepSDC(dogParams.get_time_order(), tstart, dt, dtvec, tvec);

        // Construct initial right-hand side for SDC
        SetBndValues(aux, qnew);
        BeforeStep(dtvec.get(1), aux, qnew);
        ConstructL(aux, qnew, L0, smax);
    }

    // -----------------------
    // Main time-stepping loop
    // -----------------------
    int n_step = 0;
    while (t<tend)
    {

        // initialize time step
        int m_accept = 0;      
        n_step = n_step + 1;

        // check if max number of time steps exceeded
        if (n_step>nv)
        {
            cout << " Error in DogSolveRK.cpp: "<< 
                " Exceeded allowed # of time steps " << endl;
            cout << "    n_step = " << n_step << endl;
            cout << "        nv = " << nv << endl;
            cout << endl;
            exit(1);
        }

        // Copy qnew into qold (in order to save data)
        CopyQ(qnew,qold);

        // keep trying until we get time step that doesn't violate CFL condition
        while (m_accept==0)
        {

            // set current time
            double told = t;
            if (told+dt > tend)
            { dt = tend - told; }
            t = told + dt;

            // create current time step vector
            TimeStepSDC(torder, told, dt, dtvec, tvec);

            // Set initial maximum wave speed to zero
            smax.setall(0.);

            // Choose the SDC order of accuracy
            switch( torder )
            {
//              case 2:  // 2nd order in time (this is identical to RK2)
                    // --------------------------------------------------------
//                  CopyQ( qnew, q1 );
//                  StepSDCRK2(dt, method, node, smax, L0, Ls, aux, q1, qs, qnew);
                    // --------------------------------------------------------
//                  break;

                case 3:  // 3rd order in time

                    // --------------------------------------------------------
                    // Take 1st Euler time step
                    SetBndValues(aux, qnew);
                    EulerStep(dtvec.get(1),smax,L0,aux,qnew,q1);
                    AfterStep(dtvec.get(1),aux,q1);

                    // Take 2nd Euler time step
                    SetBndValues(aux, q1);
                    BeforeStep(dtvec.get(2),aux,q1);
                    ConstructL(aux,q1,L1,smax);
                    EulerStep(dtvec.get(2),smax,L1,aux,q1,q2);
                    AfterStep(dtvec.get(2),aux,q2);

                    // Construct new right-hand side
                    SetBndValues(aux, q2);
                    BeforeStep(dtvec.get(1),aux,q2);
                    ConstructL(aux,q2,L2,smax);

                    // Iterate to construct 
                    //for(int N=1; N <= (torder-1); N++)
                    for(int N=1; N <= torder; N++)
                    {
                        // Integrate residual over time step
                        ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                        // First sub-interval    
                        //dogState->set_dt(dtvec.get(1));
                        for(int i=(1-mbc); i<=(mx+mbc); i++)
                        for(int m=1; m<=meqn; m++)
                        {    
                            // Get error in first sub-interval
                            double err1 = -(q1.get(i,m)-qnew.get(i,m)) + IL.get(i,m,1);

                            // Correct solution
                            q1.set(i,m, q1.get(i,m) + err1 );
                        }
                        AfterStep(dtvec.get(1),aux,q1);

                        // Construct new right-hand side
                        SetBndValues(aux, q1);
                        BeforeStep(dtvec.get(2),aux,q1);
                        ConstructL(aux,q1,L1new,smax);

                        // Second sub-interval            
                        for(int i=(1-mbc); i<=(mx+mbc); i++)
                        for(int m=1; m<=meqn; m++)
                        {
                            // Get error in second sub-interval
                            double err2 = dtvec.get(2)*
                                (L1new.get(i,m) - L1.get(i,m)) 
                                - (q2.get(i,m)-q1.get(i,m)) 
                                + IL.get(i,m,2);

                            // Correct solution
                            q2.set(i,m, q2.get(i,m) + err2 );
                            L1.set(i,m, L1new.get(i,m) );
                        }        
                        AfterStep(dtvec.get(2),aux,q2);

                        // Construct new right-hand side
                        SetBndValues(aux, q2);
                        BeforeStep(dtvec.get(1),aux,q2);
                        ConstructL(aux,q2,L2,smax);
                    }  

                    // Update solution
                    qnew.copyfrom(q2);
                    // --------------------------------------------------------

                    break;

/*
                case 4:  // 4th order in time

                    // RK2 Time Steps on Q
                    StepSDCRK2(dtvec.get(1), method, node, smax, L0, Ls, aux, qnew, qs, q1);
                    StepSDCRK2(dtvec.get(2), method, node, smax, L1, Ls, aux, q1, qs, q2);
                    StepSDCRK2(dtvec.get(3), method, node, smax, L2, Ls, aux, q2, qs, q3);

                    // Construct new right-hand side at final time point
                    SetBndValues(aux, q3);
                    BeforeStep(dtvec.get(1),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3,smax);

                    // Add in single correction, with RK2 time steps:
                    // Integrate residual over time step
                    ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                    //                  CopyQ(L0,L5);
                    //                  StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                    //                      qs, Ls, L0, L5, L1, qnew, q1, 1, IL);

                    StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                            qs, Ls, L0, L0, L1, qnew, q1, 1, IL);

                    // Second Interval
                    SetBndValues(aux, q1);
                    BeforeStep(dt, node, aux, q1);
                    ConstructL(method,node,aux,q1,L1new,smax);
                    StepSDCdeltaRK2(dtvec.get(2), method, node, smax, aux, 
                            qs, Ls, L1, L1new, L2, q1, q2, 2, IL );

                    // Third Interval
                    SetBndValues(aux, q2);
                    BeforeStep(dt, node, aux, q2);
                    ConstructL(method,node,aux,q2,L2new,smax);
                    StepSDCdeltaRK2(dtvec.get(3), method, node, smax, aux, 
                            qs, Ls, L2, L2new, L3, q2, q3, 3, IL );

                    // Update solution
                    CopyQ(q3,qnew);

                    break;


                     //---------------------------------------------------
                     // Start of Euler Updates ---------------------------
                     //---------------------------------------------------
                     for(N=1; N<= 3; N++ )
                     {

                    // First sub-interval    
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {    
                    // Get error in first sub-interval
                    err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                    + IL.get(i,m,k,1);

                    // Correct solution
                    q1.set(i,m,k, q1.get(i,m,k) + err1 );
                    }
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q1);  }
                    AfterStep(dtvec.get(1),node,aux,q1);

                    // Construct new right-hand side
                    SetBndValues(aux, q1);
                    BeforeStep(dtvec.get(2),node,aux,q1);
                    ConstructL(method,node,aux,q1,L1new,smax);

                    // Second sub-interval    
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {
                    // Get error in second sub-interval
                    err2 = dtvec.get(2)*
                    (L1new.get(i,m,k) - L1.get(i,m,k)) 
                    - (q2.get(i,m,k)-q1.get(i,m,k)) 
                    + IL.get(i,m,k,2);

                    // Correct solution
                    q2.set(i,m,k, q2.get(i,m,k) + err2 );
                    L1.set(i,m,k, L1new.get(i,m,k) );
                    }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q2);  }
                    AfterStep(dtvec.get(2),node,aux,q2);

                    // Construct new right-hand side
                    SetBndValues(aux, q2);
                    BeforeStep(dtvec.get(3),node,aux,q2);
                    ConstructL(method,node,aux,q2,L2new,smax);

                    // Third sub-interval    
                    //dogState->set_dt(dtvec.get(3));
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                    for (m=1; m<=meqn; m++)
                    for (k=1; k<=method[1]; k++)
                    {
                    // Get error in third sub-interval
                    err3 = dtvec.get(3)*
                    (L2new.get(i,m,k) - L2.get(i,m,k)) 
                    - (q3.get(i,m,k)-q2.get(i,m,k)) 
                    + IL.get(i,m,k,3);

                    // Correct solution
                    q3.set(i,m,k, q3.get(i,m,k) + err3 );
                    L2.set(i,m,k, L2new.get(i,m,k) );
                    }    
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q3);  }
                    AfterStep(dtvec.get(3),node,aux,q3);

                    // Construct new right-hand side
                    SetBndValues(aux, q3);
                    BeforeStep(dtvec.get(4),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3new,smax);

                    // Fourth sub-interval
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                        for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                                // Get error in fourth sub-interval
                                err4 = dtvec.get(4)*
                                    (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                    - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                    + IL.get(i,m,k,4);

                                // Correct solution
                                q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                L3.set(i,m,k, L3new.get(i,m,k) );
                            }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q4);  }
                    AfterStep(dtvec.get(4),node,aux,q4);

                    // Construct new right-hand side
                    SetBndValues(aux, q4);
                    BeforeStep(dtvec.get(5),node,aux,q4);
                    ConstructL(method,node,aux,q4,L4new,smax);

                    // Fifth sub-interval
                    //dogState->set_dt(dtvec.get(5));
                    for (i=(1-mbc); i<=(melems+mbc); i++)
                        for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                                // Get error in fourth sub-interval
                                err5 = dtvec.get(5)*
                                    (L4new.get(i,m,k) - L4.get(i,m,k)) 
                                    - (q5.get(i,m,k)-q4.get(i,m,k)) 
                                    + IL.get(i,m,k,5);

                                // Correct solution
                                q5.set(i,m,k, q5.get(i,m,k) + err5 );
                                L4.set(i,m,k, L4new.get(i,m,k) );
                            }        
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q5,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q5);  }
                    AfterStep(dtvec.get(5),node,aux,q5);

                    // Construct new right-hand side
                    SetBndValues(aux, q5);
                    BeforeStep(dtvec.get(1),node,aux,q5);
                    ConstructL(method,node,aux,q5,L5,smax);

                    // --------------------------------------------------------
                    // Take 1st Euler time step            
                    //dogState->set_dt(dtvec.get(1));
                    EulerStepSDC(dtvec.get(1),aux,qnew,L0,q1);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q1);  }
                    AfterStep(dtvec.get(1),node,aux,q1);

                    // Take 2nd Euler time step
                    //dogState->set_dt(dtvec.get(2));
                    SetBndValues(aux, q1);
                    BeforeStep(dtvec.get(2),node,aux,q1);
                    ConstructL(method,node,aux,q1,L1,smax);
                    EulerStepSDC(dtvec.get(2),aux,q1,L1,q2);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q2);  }
                    AfterStep(dtvec.get(2),node,aux,q2);

                    // Take 3rd Euler time step
                    //dogState->set_dt(dtvec.get(3));
                    SetBndValues(aux, q2);
                    BeforeStep(dtvec.get(3),node,aux,q2);
                    ConstructL(method,node,aux,q2,L2,smax);
                    EulerStepSDC(dtvec.get(3),aux,q2,L2,q3);
                    if (dogParams.using_moment_limiter())
                    { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                    else if (dogParams.using_relax_limiter())
                    {  RelaxLimiter(node,aux,q3);  }
                    AfterStep(dtvec.get(3),node,aux,q3);

                    // Construct new right-hand side
                    SetBndValues(aux, q3);
                    BeforeStep(dtvec.get(1),node,aux,q3);
                    ConstructL(method,node,aux,q3,L3,smax);

                    // Iterate to construct 
                    for (N=1; N<=(method[2]-1); N++)
                    {
                        // Integrate residual over time step
                        ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                        // First sub-interval    
                        //dogState->set_dt(dtvec.get(1));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {    
                                    // Get error in first sub-interval
                                    err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                                        + IL.get(i,m,k,1);

                                    // Correct solution
                                    q1.set(i,m,k, q1.get(i,m,k) + err1 );
                                }
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q1);  }
                        AfterStep(dtvec.get(1),node,aux,q1);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(2),node,aux,q1);
                        ConstructL(method,node,aux,q1,L1new,smax);

                        // Second sub-interval    
                        //dogState->set_dt(dtvec.get(2));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {
                                    // Get error in second sub-interval
                                    err2 = dtvec.get(2)*
                                        (L1new.get(i,m,k) - L1.get(i,m,k)) 
                                        - (q2.get(i,m,k)-q1.get(i,m,k)) 
                                        + IL.get(i,m,k,2);

                                    // Correct solution
                                    q2.set(i,m,k, q2.get(i,m,k) + err2 );
                                    L1.set(i,m,k, L1new.get(i,m,k) );
                                }        
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q2);  }
                        AfterStep(dtvec.get(2),node,aux,q2);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(3),node,aux,q2);
                        ConstructL(method,node,aux,q2,L2new,smax);

                        // Third sub-interval    
                        //dogState->set_dt(dtvec.get(3));
                        for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                                for (k=1; k<=method[1]; k++)
                                {
                                    // Get error in third sub-interval
                                    err3 = dtvec.get(3)*
                                        (L2new.get(i,m,k) - L2.get(i,m,k)) 
                                        - (q3.get(i,m,k)-q2.get(i,m,k)) 
                                        + IL.get(i,m,k,3);

                                    // Correct solution
                                    q3.set(i,m,k, q3.get(i,m,k) + err3 );
                                    L2.set(i,m,k, L2new.get(i,m,k) );
                                }        
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q3);  }
                        AfterStep(dtvec.get(3),node,aux,q3);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(1),node,aux,q3);
                        ConstructL(method,node,aux,q3,L3,smax);
                    }  

                    // Update solution
                    CopyQ(q3,qnew);
                    // --------------------------------------------------------

                    break;

                    // --------------------------------
                    // ---- End of Euler updates ------
                    // --------------------------------

                case 5:  // 5th order in time

                        // --------------------------------------------------------
                        // Take 1st Euler time step
                        //dogState->set_dt(dtvec.get(1));
                        EulerStepSDC(dtvec.get(1),aux,qnew,L0,q1);
                        if (dogParams.using_moment_limiter())
                        { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                        else if (dogParams.using_relax_limiter())
                        {  RelaxLimiter(node,aux,q1);  }
                        AfterStep(dtvec.get(1),node,aux,q1);

                        // Take 2nd Euler time step
                        StepSDC(dtvec.get(2), method, node, smax, L1, aux, q1, q2);

                        // Take 3rd Euler time step
                        StepSDC(dtvec.get(3), method, node, smax, L2, aux, q2, q3);

                        // Take 4th Euler time step
                        StepSDC(dtvec.get(4), method, node, smax, L3, aux, q3, q4);

                        // Construct new right-hand side
                        BeforeStep(dtvec.get(1),node,aux,q4);
                        ConstructL(method,node,aux,q4,L4,smax);

                        // Iterate to construct 
                        for (N=1; N<=(method[2]-1); N++)
                        {
                            // Integrate residual over time step
                            ResInt(dt,L0,L1,L2,L3,L4,L5,IL);

                            // First sub-interval    
                            //dogState->set_dt(dtvec.get(1));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {    
                                        // Get error in first sub-interval
                                        err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                                            + IL.get(i,m,k,1);

                                        // Correct solution
                                        q1.set(i,m,k, q1.get(i,m,k) + err1 );
                                    }
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q1);  }
                            AfterStep(dtvec.get(1),node,aux,q1);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(2),node,aux,q1);
                            ConstructL(method,node,aux,q1,L1new,smax);

                            // Second sub-interval    
                            //dogState->set_dt(dtvec.get(2));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in second sub-interval
                                        err2 = dtvec.get(2)*
                                            (L1new.get(i,m,k) - L1.get(i,m,k)) 
                                            - (q2.get(i,m,k)-q1.get(i,m,k)) 
                                            + IL.get(i,m,k,2);

                                        // Correct solution
                                        q2.set(i,m,k, q2.get(i,m,k) + err2 );
                                        L1.set(i,m,k, L1new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q2);  }
                            AfterStep(dtvec.get(2),node,aux,q2);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(3),node,aux,q2);
                            ConstructL(method,node,aux,q2,L2new,smax);

                            // Third sub-interval    
                            //dogState->set_dt(dtvec.get(3));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in third sub-interval
                                        err3 = dtvec.get(3)*
                                            (L2new.get(i,m,k) - L2.get(i,m,k)) 
                                            - (q3.get(i,m,k)-q2.get(i,m,k)) 
                                            + IL.get(i,m,k,3);

                                        // Correct solution
                                        q3.set(i,m,k, q3.get(i,m,k) + err3 );
                                        L2.set(i,m,k, L2new.get(i,m,k) );
                                    }    
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q3);  }
                            AfterStep(dtvec.get(3),node,aux,q3);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(4),node,aux,q3);
                            ConstructL(method,node,aux,q3,L3new,smax);

                            // Fourth sub-interval
                            //dogState->set_dt(dtvec.get(4));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err4 = dtvec.get(4)*
                                            (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                            - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                            + IL.get(i,m,k,4);

                                        // Correct solution
                                        q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                        L3.set(i,m,k, L3new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q4);  }
                            AfterStep(dtvec.get(4),node,aux,q4);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q4);
                            ConstructL(method,node,aux,q4,L4,smax);
                        }  

                        // Update solution
                        CopyQ(q4,qnew);
                        // --------------------------------------------------------

                        break;            

                case 6:  // 6th order in time

                        // RK2 Time Steps on Q
                        StepSDCRK2(dtvec.get(1), method, node, smax, L0, Ls, aux, qnew, qs, q1);
                        StepSDCRK2(dtvec.get(2), method, node, smax, L1, Ls, aux, q1, qs, q2);
                        StepSDCRK2(dtvec.get(3), method, node, smax, L2, Ls, aux, q2, qs, q3);
                        StepSDCRK2(dtvec.get(4), method, node, smax, L3, Ls, aux, q3, qs, q4);
                        StepSDCRK2(dtvec.get(5), method, node, smax, L4, Ls, aux, q4, qs, q5);

                        // Construct new right-hand side at final time point
                        BeforeStep(dtvec.get(1),node,aux,q5);
                        ConstructL(method,node,aux,q5,L5,smax);

                      //StepSDC(dtvec.get(1), method, node, smax, L0, aux, qnew, q1);

                      //// Take 2nd Euler time step
                      //StepSDC(dtvec.get(2), method, node, smax, L1, aux, q1, q2);

                      //// Take 3rd Euler time step
                      //StepSDC(dtvec.get(3), method, node, smax, L2, aux, q2, q3);

                      //// Take 4th Euler time step
                      //StepSDC(dtvec.get(4), method, node, smax, L3, aux, q3, q4);

                      //// Take 5th Euler time step
                      //StepSDC(dtvec.get(5), method, node, smax, L4, aux, q4, q5);

                      //// Construct new right-hand side at final time point
                      //BeforeStep(dtvec.get(1),node,aux,q5);
                      //ConstructL(method,node,aux,q5,L5,smax);

                        for (N=1; N<=3; N++)
                        {
                            // Integrate residual over time step
                            ResInt(dt,L0,L1,L2,L3,L4,L5,IL);            

                            StepSDCdeltaRK2(dtvec.get(1), method, node, smax, aux,
                                    qs, Ls, L0, L0, L1, qnew, q1, 1, IL);

                            // Second Interval
                            BeforeStep(dt, node, aux, q1);
                            ConstructL(method,node,aux,q1,L1new,smax);
                            StepSDCdeltaRK2(dtvec.get(2), method, node, smax, aux, 
                                    qs, Ls, L1, L1new, L2, q1, q2, 2, IL );

                            // Third Interval
                            BeforeStep(dt, node, aux, q2);
                            ConstructL(method,node,aux,q2,L2new,smax);
                            StepSDCdeltaRK2(dtvec.get(3), method, node, smax, aux, 
                                    qs, Ls, L2, L2new, L3, q2, q3, 3, IL );

                            // Fourth Interval
                            BeforeStep(dt, node, aux, q3);
                            ConstructL(method,node,aux,q3,L3new,smax);
                            StepSDCdeltaRK2(dtvec.get(4), method, node, smax, aux, 
                                    qs, Ls, L3, L3new, L4, q3, q4, 4, IL );

                            // Final Interval
                            BeforeStep(dt, node, aux, q4);
                            ConstructL(method,node,aux,q4,L4new,smax);
                            StepSDCdeltaRK2(dtvec.get(5), method, node, smax, aux, 
                                    qs, Ls, L4, L4new, L5, q4, q5, 5, IL );

                            // Save the L1s for the next integration
                            CopyQ( L1new, L1 );
                            CopyQ( L2new, L2 );
                            CopyQ( L3new, L3 );
                            CopyQ( L4new, L4 );

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q5);
                            ConstructL(method,node,aux,q5,L5,smax);

                            // ---------------------------------------------------
                            // -- Start of Euler Updates -------------------------
                            // ---------------------------------------------------

                            // First sub-interval    
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {    
                            // Get error in first sub-interval
                            err1 = -(q1.get(i,m,k)-qnew.get(i,m,k)) 
                            + IL.get(i,m,k,1);

                            // Correct solution
                            q1.set(i,m,k, q1.get(i,m,k) + err1 );
                            }
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q1,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q1);  }
                            AfterStep(dtvec.get(1),node,aux,q1);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(2),node,aux,q1);
                            ConstructL(method,node,aux,q1,L1new,smax);

                            // Second sub-interval    
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                            // Get error in second sub-interval
                            err2 = dtvec.get(2)*
                            (L1new.get(i,m,k) - L1.get(i,m,k)) 
                            - (q2.get(i,m,k)-q1.get(i,m,k)) 
                            + IL.get(i,m,k,2);

                            // Correct solution
                            q2.set(i,m,k, q2.get(i,m,k) + err2 );
                            L1.set(i,m,k, L1new.get(i,m,k) );
                            }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q2,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q2);  }
                            AfterStep(dtvec.get(2),node,aux,q2);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(3),node,aux,q2);
                            ConstructL(method,node,aux,q2,L2new,smax);

                            // Third sub-interval    
                            //dogState->set_dt(dtvec.get(3));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                            for (m=1; m<=meqn; m++)
                            for (k=1; k<=method[1]; k++)
                            {
                            // Get error in third sub-interval
                            err3 = dtvec.get(3)*
                            (L2new.get(i,m,k) - L2.get(i,m,k)) 
                            - (q3.get(i,m,k)-q2.get(i,m,k)) 
                            + IL.get(i,m,k,3);

                            // Correct solution
                            q3.set(i,m,k, q3.get(i,m,k) + err3 );
                            L2.set(i,m,k, L2new.get(i,m,k) );
                            }    
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q3,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q3);  }
                            AfterStep(dtvec.get(3),node,aux,q3);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(4),node,aux,q3);
                            ConstructL(method,node,aux,q3,L3new,smax);

                            // Fourth sub-interval
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err4 = dtvec.get(4)*
                                            (L3new.get(i,m,k) - L3.get(i,m,k)) 
                                            - (q4.get(i,m,k)-q3.get(i,m,k)) 
                                            + IL.get(i,m,k,4);

                                        // Correct solution
                                        q4.set(i,m,k, q4.get(i,m,k) + err4 );
                                        L3.set(i,m,k, L3new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q4,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q4);  }
                            AfterStep(dtvec.get(4),node,aux,q4);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(5),node,aux,q4);
                            ConstructL(method,node,aux,q4,L4new,smax);

                            // Fifth sub-interval
                            //dogState->set_dt(dtvec.get(5));
                            for (i=(1-mbc); i<=(melems+mbc); i++)
                                for (m=1; m<=meqn; m++)
                                    for (k=1; k<=method[1]; k++)
                                    {
                                        // Get error in fourth sub-interval
                                        err5 = dtvec.get(5)*
                                            (L4new.get(i,m,k) - L4.get(i,m,k)) 
                                            - (q5.get(i,m,k)-q4.get(i,m,k)) 
                                            + IL.get(i,m,k,5);

                                        // Correct solution
                                        q5.set(i,m,k, q5.get(i,m,k) + err5 );
                                        L4.set(i,m,k, L4new.get(i,m,k) );
                                    }        
                            if (dogParams.using_moment_limiter())
                            { ApplyLimiter(node,aux,q5,&ProjectRightEig,&ProjectLeftEig); }
                            else if (dogParams.using_relax_limiter())
                            {  RelaxLimiter(node,aux,q5);  }
                            AfterStep(dtvec.get(5),node,aux,q5);

                            // Construct new right-hand side
                            BeforeStep(dtvec.get(1),node,aux,q5);
                            ConstructL(method,node,aux,q5,L5,smax);

                    // --------------------------------
                    // ---- End of Euler updates ------
                    // --------------------------------

                        }

                        // Update solution
                        CopyQ(q5,qnew);

                        break;            
*/
            }

            // compute cfl number
            cfl = GetCFL(dt, dtv[2], aux, smax);

            // output time step information
            if( dogParams.get_verbosity() )
            {
                cout << setprecision(3);
                cout << "FinSolveSDC1D ... Step" << setw(5) << n_step;
                cout << "   CFL =" << setw(6) << fixed << cfl;
                cout << "   dt =" << setw(11) << scientific << dt;
                cout << "   t =" << setw(11) << scientific << t <<endl;
            }

            // choose new time step
            if (cfl>0.0)
            {   
                dt = Min(dtv[2],dt*CFL_target/cfl);
                dtmin = Min(dt,dtmin);
                dtmax = Max(dt,dtmax);
            }
            else
            {
                dt = dtv[2];
            }

            // see whether to accept or reject this step
            if (cfl<=CFL_max)
            // accept
            {  
                m_accept = 1;  

                // Copy old L into new L0
                if( torder==3 )
                {  L0.copyfrom(L2); }

                if (torder==4)
                {  L0.copyfrom(L3); }

                if (torder==5)
                {  L0.copyfrom(L4); }

            }
            else //reject
            {   
                t = told;
                if( dogParams.get_verbosity() )
                {
                    cout<<"FinSolve1D rejecting step...";
                    cout<<"CFL number too large";
                    cout<<endl;
                }

                // Copy old data back into qnew
                qnew.copyfrom( qold );

            }

        } // End of m_accept loop

        // compute conservation and print to file
        ConSoln(aux, qnew, t, outputdir);

    } // End of while loop

    // set initial time step for next call to DogSolve
    dtv[1] = dt;

}
