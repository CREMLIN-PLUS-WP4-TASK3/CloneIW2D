/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses

>>> END OF LICENSE >>>
*************************************************************************/

#ifndef _minlm_h
#define _minlm_h

#include "ap.h"
#include "amp.h"
#include "blas.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
namespace minlm
{
    template<unsigned int Precision>
    class minlmstate
    {
    public:
        bool wrongparams;
        int n;
        int m;
        amp::ampf<Precision> epsg;
        amp::ampf<Precision> epsf;
        amp::ampf<Precision> epsx;
        int maxits;
        bool xrep;
        amp::ampf<Precision> stpmax;
        int flags;
        int usermode;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> f;
        ap::template_1d_array< amp::ampf<Precision> > fi;
        ap::template_2d_array< amp::ampf<Precision> > j;
        ap::template_2d_array< amp::ampf<Precision> > h;
        ap::template_1d_array< amp::ampf<Precision> > g;
        bool needf;
        bool needfg;
        bool needfgh;
        bool needfij;
        bool xupdated;
        minlbfgs::minlbfgsstate<Precision> internalstate;
        minlbfgs::minlbfgsreport<Precision> internalrep;
        ap::template_1d_array< amp::ampf<Precision> > xprec;
        ap::template_1d_array< amp::ampf<Precision> > xbase;
        ap::template_1d_array< amp::ampf<Precision> > xdir;
        ap::template_1d_array< amp::ampf<Precision> > gbase;
        ap::template_1d_array< amp::ampf<Precision> > xprev;
        amp::ampf<Precision> fprev;
        ap::template_2d_array< amp::ampf<Precision> > rawmodel;
        ap::template_2d_array< amp::ampf<Precision> > model;
        ap::template_1d_array< amp::ampf<Precision> > work;
        amp::rcommstate<Precision> rstate;
        int repiterationscount;
        int repterminationtype;
        int repnfunc;
        int repnjac;
        int repngrad;
        int repnhess;
        int repncholesky;
        int solverinfo;
        densesolver::densesolverreport<Precision> solverrep;
        int invinfo;
        matinv::matinvreport<Precision> invrep;
    };


    template<unsigned int Precision>
    class minlmreport
    {
    public:
        int iterationscount;
        int terminationtype;
        int nfunc;
        int njac;
        int ngrad;
        int nhess;
        int ncholesky;
    };




    template<unsigned int Precision>
    void minlmcreatefgh(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state);
    template<unsigned int Precision>
    void minlmcreatefgj(const int& n,
        const int& m,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state);
    template<unsigned int Precision>
    void minlmcreatefj(const int& n,
        const int& m,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state);
    template<unsigned int Precision>
    void minlmsetcond(minlmstate<Precision>& state,
        amp::ampf<Precision> epsg,
        amp::ampf<Precision> epsf,
        amp::ampf<Precision> epsx,
        int maxits);
    template<unsigned int Precision>
    void minlmsetxrep(minlmstate<Precision>& state,
        bool needxrep);
    template<unsigned int Precision>
    void minlmsetstpmax(minlmstate<Precision>& state,
        amp::ampf<Precision> stpmax);
    template<unsigned int Precision>
    bool minlmiteration(minlmstate<Precision>& state);
    template<unsigned int Precision>
    void minlmresults(const minlmstate<Precision>& state,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmreport<Precision>& rep);
    template<unsigned int Precision>
    void lmprepare(int n,
        int m,
        bool havegrad,
        minlmstate<Precision>& state);
    template<unsigned int Precision>
    void lmclearrequestfields(minlmstate<Precision>& state);
    template<unsigned int Precision>
    bool increaselambda(amp::ampf<Precision>& lambda,
        amp::ampf<Precision>& nu,
        amp::ampf<Precision> lambdaup);
    template<unsigned int Precision>
    void decreaselambda(amp::ampf<Precision>& lambda,
        amp::ampf<Precision>& nu,
        amp::ampf<Precision> lambdadown);


    static const int lmmodefj = 0;
    static const int lmmodefgj = 1;
    static const int lmmodefgh = 2;
    static const int lmflagnoprelbfgs = 1;
    static const int lmflagnointlbfgs = 2;
    static const int lmprelbfgsm = 5;
    static const int lmintlbfgsits = 5;
    static const int lbfgsnorealloc = 1;


    /*************************************************************************
        LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

    Optimization using function gradient and Hessian.  Algorithm -  Levenberg-
    Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
    pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

    Function F has general form (not "sum-of-squares"):

        F = F(x[0], ..., x[n-1])

    EXAMPLE

    See HTML-documentation.

    INPUT PARAMETERS:
        N       -   dimension, N>1
        X       -   initial solution, array[0..N-1]

    OUTPUT PARAMETERS:
        State   -   structure which stores algorithm state between subsequent
                    calls of MinLMIteration. Used for reverse communication.
                    This structure should be passed to MinLMIteration subroutine.

    See also MinLMIteration, MinLMResults.

    NOTES:

    1. you may tune stopping conditions with MinLMSetCond() function
    2. if target function contains exp() or other fast growing functions,  and
       optimization algorithm makes too large steps which leads  to  overflow,
       use MinLMSetStpMax() function to bound algorithm's steps.

      -- ALGLIB --
         Copyright 30.03.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmcreatefgh(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state)
    {
        
        //
        // Prepare RComm
        //
        state.rstate.ia.setbounds(0, 3);
        state.rstate.ba.setbounds(0, 0);
        state.rstate.ra.setbounds(0, 7);
        state.rstate.stage = -1;
        
        //
        // prepare internal structures
        //
        lmprepare<Precision>(n, 0, true, state);
        
        //
        // initialize, check parameters
        //
        minlmsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        minlmsetxrep<Precision>(state, false);
        minlmsetstpmax<Precision>(state, amp::ampf<Precision>(0));
        state.n = n;
        state.m = 0;
        state.flags = 0;
        state.usermode = lmmodefgh;
        state.wrongparams = false;
        if( n<1 )
        {
            state.wrongparams = true;
            return;
        }
        amp::vmove(state.x.getvector(0, n-1), x.getvector(0, n-1));
    }


    /*************************************************************************
        LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

    Optimization using function gradient and Jacobian.  Algorithm -  Levenberg-
    Marquardt   modification   with   L-BFGS   pre-optimization  and  internal
    pre-conditioned L-BFGS optimization after each Levenberg-Marquardt step.

    Function F is represented as sum of squares:

        F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

    EXAMPLE

    See HTML-documentation.

    INPUT PARAMETERS:
        N       -   dimension, N>1
        M       -   number of functions f[i]
        X       -   initial solution, array[0..N-1]

    OUTPUT PARAMETERS:
        State   -   structure which stores algorithm state between subsequent
                    calls of MinLMIteration. Used for reverse communication.
                    This structure should be passed to MinLMIteration subroutine.

    See also MinLMIteration, MinLMResults.

    NOTES:

    1. you may tune stopping conditions with MinLMSetCond() function
    2. if target function contains exp() or other fast growing functions,  and
       optimization algorithm makes too large steps which leads  to  overflow,
       use MinLMSetStpMax() function to bound algorithm's steps.

      -- ALGLIB --
         Copyright 30.03.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmcreatefgj(const int& n,
        const int& m,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state)
    {
        
        //
        // Prepare RComm
        //
        state.rstate.ia.setbounds(0, 3);
        state.rstate.ba.setbounds(0, 0);
        state.rstate.ra.setbounds(0, 7);
        state.rstate.stage = -1;
        
        //
        // prepare internal structures
        //
        lmprepare<Precision>(n, m, true, state);
        
        //
        // initialize, check parameters
        //
        minlmsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        minlmsetxrep<Precision>(state, false);
        minlmsetstpmax<Precision>(state, amp::ampf<Precision>(0));
        state.n = n;
        state.m = m;
        state.flags = 0;
        state.usermode = lmmodefgj;
        state.wrongparams = false;
        if( n<1 )
        {
            state.wrongparams = true;
            return;
        }
        amp::vmove(state.x.getvector(0, n-1), x.getvector(0, n-1));
    }


    /*************************************************************************
        CLASSIC LEVENBERG-MARQUARDT METHOD FOR NON-LINEAR OPTIMIZATION

    Optimization using Jacobi matrix. Algorithm  -  classic Levenberg-Marquardt
    method.

    Function F is represented as sum of squares:

        F = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])

    EXAMPLE

    See HTML-documentation.

    INPUT PARAMETERS:
        N       -   dimension, N>1
        M       -   number of functions f[i]
        X       -   initial solution, array[0..N-1]

    OUTPUT PARAMETERS:
        State   -   structure which stores algorithm state between subsequent
                    calls of MinLMIteration. Used for reverse communication.
                    This structure should be passed to MinLMIteration subroutine.

    See also MinLMIteration, MinLMResults.

    NOTES:

    1. you may tune stopping conditions with MinLMSetCond() function
    2. if target function contains exp() or other fast growing functions,  and
       optimization algorithm makes too large steps which leads  to  overflow,
       use MinLMSetStpMax() function to bound algorithm's steps.

      -- ALGLIB --
         Copyright 30.03.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmcreatefj(const int& n,
        const int& m,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmstate<Precision>& state)
    {
        
        //
        // Prepare RComm
        //
        state.rstate.ia.setbounds(0, 3);
        state.rstate.ba.setbounds(0, 0);
        state.rstate.ra.setbounds(0, 7);
        state.rstate.stage = -1;
        
        //
        // prepare internal structures
        //
        lmprepare<Precision>(n, m, true, state);
        
        //
        // initialize, check parameters
        //
        minlmsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        minlmsetxrep<Precision>(state, false);
        minlmsetstpmax<Precision>(state, amp::ampf<Precision>(0));
        state.n = n;
        state.m = m;
        state.flags = 0;
        state.usermode = lmmodefj;
        state.wrongparams = false;
        if( n<1 )
        {
            state.wrongparams = true;
            return;
        }
        amp::vmove(state.x.getvector(0, n-1), x.getvector(0, n-1));
    }


    /*************************************************************************
    This function sets stopping conditions for Levenberg-Marquardt optimization
    algorithm.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be initialized
                    with MinLMCreate???()
        EpsG    -   >=0
                    The  subroutine  finishes  its  work   if   the  condition
                    ||G||<EpsG is satisfied, where ||.|| means Euclidian norm,
                    G - gradient.
        EpsF    -   >=0
                    The  subroutine  finishes  its work if on k+1-th iteration
                    the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                    is satisfied.
        EpsX    -   >=0
                    The subroutine finishes its work if  on  k+1-th  iteration
                    the condition |X(k+1)-X(k)| <= EpsX is fulfilled.
        MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                    iterations   is    unlimited.   Only   Levenberg-Marquardt
                    iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                    counted  because their cost is very low copared to that of
                    LM).

    Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
    automatic stopping criterion selection (small EpsX).

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmsetcond(minlmstate<Precision>& state,
        amp::ampf<Precision> epsg,
        amp::ampf<Precision> epsf,
        amp::ampf<Precision> epsx,
        int maxits)
    {
        ap::ap_error::make_assertion(epsg>=0);
        ap::ap_error::make_assertion(epsf>=0);
        ap::ap_error::make_assertion(epsx>=0);
        ap::ap_error::make_assertion(maxits>=0);
        if( epsg==0 && epsf==0 && epsx==0 && maxits==0 )
        {
            epsx = amp::ampf<Precision>("1.0E-6");
        }
        state.epsg = epsg;
        state.epsf = epsf;
        state.epsx = epsx;
        state.maxits = maxits;
    }


    /*************************************************************************
    This function turns on/off reporting.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be
                    initialized with MinLMCreate???()
        NeedXRep-   whether iteration reports are needed or not

    Usually  algorithm  returns  from  MinLMIteration()  only  when  it  needs
    function/gradient/Hessian. However, with this function we can let it  stop
    after  each  iteration  (one iteration may include  more than one function
    evaluation), which is indicated by XUpdated field.

    Both Levenberg-Marquardt and L-BFGS iterations are reported.


      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmsetxrep(minlmstate<Precision>& state,
        bool needxrep)
    {
        state.xrep = needxrep;
    }


    /*************************************************************************
    This function sets maximum step length

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be
                    initialized with MinCGCreate???()
        StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                    want to limit step length.

    Use this subroutine when you optimize target function which contains exp()
    or  other  fast  growing  functions,  and optimization algorithm makes too
    large  steps  which  leads  to overflow. This function allows us to reject
    steps  that  are  too  large  (and  therefore  expose  us  to the possible
    overflow) without actually calculating function value at the x+stp*d.

    NOTE: non-zero StpMax leads to moderate  performance  degradation  because
    intermediate  step  of  preconditioned L-BFGS optimization is incompatible
    with limits on step size.

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmsetstpmax(minlmstate<Precision>& state,
        amp::ampf<Precision> stpmax)
    {
        ap::ap_error::make_assertion(stpmax>=0);
        state.stpmax = stpmax;
    }


    /*************************************************************************
    One Levenberg-Marquardt iteration.

    Called after inialization of State structure with MinLMXXX subroutine.
    See HTML docs for examples.

    Input parameters:
        State   -   structure which stores algorithm state between subsequent
                    calls and which is used for reverse communication. Must be
                    initialized with MinLMXXX call first.

    If subroutine returned False, iterative algorithm has converged.

    If subroutine returned True, then:
    * if State.NeedF=True,      -   function value F at State.X[0..N-1]
                                    is required
    * if State.NeedFG=True      -   function value F and gradient G
                                    are required
    * if State.NeedFiJ=True     -   function vector f[i] and Jacobi matrix J
                                    are required
    * if State.NeedFGH=True     -   function value F, gradient G and Hesian H
                                    are required
    * if State.XUpdated=True    -   algorithm reports about new iteration,
                                    State.X contains current point,
                                    State.F contains function value.

    One and only one of this fields can be set at time.

    Results are stored:
    * function value            -   in MinLMState.F
    * gradient                  -   in MinLMState.G[0..N-1]
    * Jacobi matrix             -   in MinLMState.J[0..M-1,0..N-1]
    * Hessian                   -   in MinLMState.H[0..N-1,0..N-1]

      -- ALGLIB --
         Copyright 10.03.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool minlmiteration(minlmstate<Precision>& state)
    {
        bool result;
        int n;
        int m;
        int i;
        amp::ampf<Precision> stepnorm;
        bool spd;
        amp::ampf<Precision> fbase;
        amp::ampf<Precision> fnew;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> nu;
        amp::ampf<Precision> lambdaup;
        amp::ampf<Precision> lambdadown;
        int lbfgsflags;
        amp::ampf<Precision> v;


        
        //
        // Reverse communication preparations
        // I know it looks ugly, but it works the same way
        // anywhere from C++ to Python.
        //
        // This code initializes locals by:
        // * random values determined during code
        //   generation - on first subroutine call
        // * values from previous call - on subsequent calls
        //
        if( state.rstate.stage>=0 )
        {
            n = state.rstate.ia(0);
            m = state.rstate.ia(1);
            i = state.rstate.ia(2);
            lbfgsflags = state.rstate.ia(3);
            spd = state.rstate.ba(0);
            stepnorm = state.rstate.ra(0);
            fbase = state.rstate.ra(1);
            fnew = state.rstate.ra(2);
            lambda = state.rstate.ra(3);
            nu = state.rstate.ra(4);
            lambdaup = state.rstate.ra(5);
            lambdadown = state.rstate.ra(6);
            v = state.rstate.ra(7);
        }
        else
        {
            n = -983;
            m = -989;
            i = -834;
            lbfgsflags = 900;
            spd = true;
            stepnorm = 364;
            fbase = 214;
            fnew = -338;
            lambda = -686;
            nu = 912;
            lambdaup = 585;
            lambdadown = 497;
            v = -271;
        }
        if( state.rstate.stage==0 )
        {
            goto lbl_0;
        }
        if( state.rstate.stage==1 )
        {
            goto lbl_1;
        }
        if( state.rstate.stage==2 )
        {
            goto lbl_2;
        }
        if( state.rstate.stage==3 )
        {
            goto lbl_3;
        }
        if( state.rstate.stage==4 )
        {
            goto lbl_4;
        }
        if( state.rstate.stage==5 )
        {
            goto lbl_5;
        }
        if( state.rstate.stage==6 )
        {
            goto lbl_6;
        }
        if( state.rstate.stage==7 )
        {
            goto lbl_7;
        }
        if( state.rstate.stage==8 )
        {
            goto lbl_8;
        }
        if( state.rstate.stage==9 )
        {
            goto lbl_9;
        }
        if( state.rstate.stage==10 )
        {
            goto lbl_10;
        }
        if( state.rstate.stage==11 )
        {
            goto lbl_11;
        }
        if( state.rstate.stage==12 )
        {
            goto lbl_12;
        }
        if( state.rstate.stage==13 )
        {
            goto lbl_13;
        }
        if( state.rstate.stage==14 )
        {
            goto lbl_14;
        }
        if( state.rstate.stage==15 )
        {
            goto lbl_15;
        }
        
        //
        // Routine body
        //
        ap::ap_error::make_assertion(state.usermode==lmmodefj || state.usermode==lmmodefgj || state.usermode==lmmodefgh);
        if( state.wrongparams )
        {
            state.repterminationtype = -1;
            result = false;
            return result;
        }
        
        //
        // prepare params
        //
        n = state.n;
        m = state.m;
        lambdaup = 20;
        lambdadown = amp::ampf<Precision>("0.5");
        nu = 1;
        lbfgsflags = 0;
        
        //
        // if we have F and G
        //
        if( ! ((state.usermode==lmmodefgj || state.usermode==lmmodefgh) && state.flags/lmflagnoprelbfgs%2==0) )
        {
            goto lbl_16;
        }
        
        //
        // First stage of the hybrid algorithm: LBFGS
        //
        minlbfgs::minlbfgscreate<Precision>(n, ap::minint(n, lmprelbfgsm), state.x, state.internalstate);
        minlbfgs::minlbfgssetcond<Precision>(state.internalstate, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), ap::maxint(5, n));
        minlbfgs::minlbfgssetxrep<Precision>(state.internalstate, state.xrep);
        minlbfgs::minlbfgssetstpmax<Precision>(state.internalstate, state.stpmax);
    lbl_18:
        if( ! minlbfgs::minlbfgsiteration<Precision>(state.internalstate) )
        {
            goto lbl_19;
        }
        if( ! state.internalstate.needfg )
        {
            goto lbl_20;
        }
        
        //
        // RComm
        //
        amp::vmove(state.x.getvector(0, n-1), state.internalstate.x.getvector(0, n-1));
        lmclearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        state.repnfunc = state.repnfunc+1;
        state.repngrad = state.repngrad+1;
        
        //
        // Call LBFGS
        //
        state.internalstate.f = state.f;
        amp::vmove(state.internalstate.g.getvector(0, n-1), state.g.getvector(0, n-1));
    lbl_20:
        if( ! (state.internalstate.xupdated && state.xrep) )
        {
            goto lbl_22;
        }
        lmclearrequestfields<Precision>(state);
        state.f = state.internalstate.f;
        amp::vmove(state.x.getvector(0, n-1), state.internalstate.x.getvector(0, n-1));
        state.xupdated = true;
        state.rstate.stage = 1;
        goto lbl_rcomm;
    lbl_1:
    lbl_22:
        goto lbl_18;
    lbl_19:
        minlbfgs::minlbfgsresults<Precision>(state.internalstate, state.x, state.internalrep);
        goto lbl_17;
    lbl_16:
        
        //
        // No first stage.
        // However, we may need to report initial point
        //
        if( ! state.xrep )
        {
            goto lbl_24;
        }
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 2;
        goto lbl_rcomm;
    lbl_2:
        lmclearrequestfields<Precision>(state);
        state.xupdated = true;
        state.rstate.stage = 3;
        goto lbl_rcomm;
    lbl_3:
    lbl_24:
    lbl_17:
        
        //
        // Second stage of the hybrid algorithm: LM
        // Initialize quadratic model.
        //
        if( state.usermode!=lmmodefgh )
        {
            goto lbl_26;
        }
        
        //
        // RComm
        //
        lmclearrequestfields<Precision>(state);
        state.needfgh = true;
        state.rstate.stage = 4;
        goto lbl_rcomm;
    lbl_4:
        state.repnfunc = state.repnfunc+1;
        state.repngrad = state.repngrad+1;
        state.repnhess = state.repnhess+1;
        
        //
        // generate raw quadratic model
        //
        ablas::rmatrixcopy<Precision>(n, n, state.h, 0, 0, state.rawmodel, 0, 0);
        amp::vmove(state.gbase.getvector(0, n-1), state.g.getvector(0, n-1));
        fbase = state.f;
    lbl_26:
        if( ! (state.usermode==lmmodefgj || state.usermode==lmmodefj) )
        {
            goto lbl_28;
        }
        
        //
        // RComm
        //
        lmclearrequestfields<Precision>(state);
        state.needfij = true;
        state.rstate.stage = 5;
        goto lbl_rcomm;
    lbl_5:
        state.repnfunc = state.repnfunc+1;
        state.repnjac = state.repnjac+1;
        
        //
        // generate raw quadratic model
        //
        ablas::rmatrixgemm<Precision>(n, n, m, amp::ampf<Precision>("2.0"), state.j, 0, 0, 1, state.j, 0, 0, 0, amp::ampf<Precision>("0.0"), state.rawmodel, 0, 0);
        ablas::rmatrixmv<Precision>(n, m, state.j, 0, 0, 1, state.fi, 0, state.gbase, 0);
        amp::vmul(state.gbase.getvector(0, n-1), 2);
        fbase = amp::vdotproduct(state.fi.getvector(0, m-1), state.fi.getvector(0, m-1));
    lbl_28:
        lambda = amp::ampf<Precision>("0.001");
    lbl_30:
        if( false )
        {
            goto lbl_31;
        }
        
        //
        // 1. Model = RawModel+lambda*I
        // 2. Try to solve (RawModel+Lambda*I)*dx = -g.
        //    Increase lambda if left part is not positive definite.
        //
        for(i=0; i<=n-1; i++)
        {
            amp::vmove(state.model.getrow(i, 0, n-1), state.rawmodel.getrow(i, 0, n-1));
            state.model(i,i) = state.model(i,i)+lambda;
        }
        spd = trfac::spdmatrixcholesky<Precision>(state.model, n, true);
        state.repncholesky = state.repncholesky+1;
        if( spd )
        {
            goto lbl_32;
        }
        if( ! increaselambda<Precision>(lambda, nu, lambdaup) )
        {
            goto lbl_34;
        }
        goto lbl_30;
        goto lbl_35;
    lbl_34:
        state.repterminationtype = 7;
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 6;
        goto lbl_rcomm;
    lbl_6:
        goto lbl_31;
    lbl_35:
    lbl_32:
        densesolver::spdmatrixcholeskysolve<Precision>(state.model, n, true, state.gbase, state.solverinfo, state.solverrep, state.xdir);
        if( state.solverinfo>=0 )
        {
            goto lbl_36;
        }
        if( ! increaselambda<Precision>(lambda, nu, lambdaup) )
        {
            goto lbl_38;
        }
        goto lbl_30;
        goto lbl_39;
    lbl_38:
        state.repterminationtype = 7;
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 7;
        goto lbl_rcomm;
    lbl_7:
        goto lbl_31;
    lbl_39:
    lbl_36:
        amp::vmul(state.xdir.getvector(0, n-1), -1);
        
        //
        // Candidate lambda is found.
        // 1. Save old w in WBase
        // 1. Test some stopping criterions
        // 2. If error(w+wdir)>error(w), increase lambda
        //
        amp::vmove(state.xprev.getvector(0, n-1), state.x.getvector(0, n-1));
        state.fprev = state.f;
        amp::vmove(state.xbase.getvector(0, n-1), state.x.getvector(0, n-1));
        amp::vadd(state.x.getvector(0, n-1), state.xdir.getvector(0, n-1));
        stepnorm = amp::vdotproduct(state.xdir.getvector(0, n-1), state.xdir.getvector(0, n-1));
        stepnorm = amp::sqrt<Precision>(stepnorm);
        if( ! (state.stpmax>0 && stepnorm>state.stpmax) )
        {
            goto lbl_40;
        }
        
        //
        // Step is larger than the limit,
        // larger lambda is needed
        //
        amp::vmove(state.x.getvector(0, n-1), state.xbase.getvector(0, n-1));
        if( ! increaselambda<Precision>(lambda, nu, lambdaup) )
        {
            goto lbl_42;
        }
        goto lbl_30;
        goto lbl_43;
    lbl_42:
        state.repterminationtype = 7;
        amp::vmove(state.x.getvector(0, n-1), state.xprev.getvector(0, n-1));
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 8;
        goto lbl_rcomm;
    lbl_8:
        goto lbl_31;
    lbl_43:
    lbl_40:
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 9;
        goto lbl_rcomm;
    lbl_9:
        state.repnfunc = state.repnfunc+1;
        fnew = state.f;
        if( fnew<=fbase )
        {
            goto lbl_44;
        }
        
        //
        // restore state and continue search for lambda
        //
        amp::vmove(state.x.getvector(0, n-1), state.xbase.getvector(0, n-1));
        if( ! increaselambda<Precision>(lambda, nu, lambdaup) )
        {
            goto lbl_46;
        }
        goto lbl_30;
        goto lbl_47;
    lbl_46:
        state.repterminationtype = 7;
        amp::vmove(state.x.getvector(0, n-1), state.xprev.getvector(0, n-1));
        lmclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 10;
        goto lbl_rcomm;
    lbl_10:
        goto lbl_31;
    lbl_47:
    lbl_44:
        if( ! (state.stpmax==0 && (state.usermode==lmmodefgj || state.usermode==lmmodefgh) && state.flags/lmflagnointlbfgs%2==0) )
        {
            goto lbl_48;
        }
        
        //
        // Optimize using LBFGS, with inv(cholesky(H)) as preconditioner.
        //
        // It is possible only when StpMax=0, because we can't guarantee
        // that step remains bounded when preconditioner is used (we need
        // SVD decomposition to do that, which is too slow).
        //
        matinv::rmatrixtrinverse<Precision>(state.model, n, true, false, state.invinfo, state.invrep);
        if( state.invinfo<=0 )
        {
            goto lbl_50;
        }
        
        //
        // if matrix can be inverted, use it.
        // just silently move to next iteration otherwise.
        // (will be very, very rare, mostly for specially
        // designed near-degenerate tasks)
        //
        amp::vmove(state.xbase.getvector(0, n-1), state.x.getvector(0, n-1));
        for(i=0; i<=n-1; i++)
        {
            state.xprec(i) = 0;
        }
        minlbfgs::minlbfgscreatex<Precision>(n, ap::minint(n, lmintlbfgsits), state.xprec, lbfgsflags, state.internalstate);
        minlbfgs::minlbfgssetcond<Precision>(state.internalstate, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), lmintlbfgsits);
    lbl_52:
        if( ! minlbfgs::minlbfgsiteration<Precision>(state.internalstate) )
        {
            goto lbl_53;
        }
        
        //
        // convert XPrec to unpreconditioned form, then call RComm.
        //
        for(i=0; i<=n-1; i++)
        {
            v = amp::vdotproduct(state.internalstate.x.getvector(i, n-1), state.model.getrow(i, i, n-1));
            state.x(i) = state.xbase(i)+v;
        }
        lmclearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 11;
        goto lbl_rcomm;
    lbl_11:
        state.repnfunc = state.repnfunc+1;
        state.repngrad = state.repngrad+1;
        
        //
        // 1. pass State.F to State.InternalState.F
        // 2. convert gradient back to preconditioned form
        //
        state.internalstate.f = state.f;
        for(i=0; i<=n-1; i++)
        {
            state.internalstate.g(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            v = state.g(i);
            amp::vadd(state.internalstate.g.getvector(i, n-1), state.model.getrow(i, i, n-1), v);
        }
        
        //
        // next iteration
        //
        goto lbl_52;
    lbl_53:
        
        //
        // change LBFGS flags to NoRealloc.
        // L-BFGS subroutine will use memory allocated from previous run.
        // it is possible since all subsequent calls will be with same N/M.
        //
        lbfgsflags = lbfgsnorealloc;
        
        //
        // back to unpreconditioned X
        //
        minlbfgs::minlbfgsresults<Precision>(state.internalstate, state.xprec, state.internalrep);
        for(i=0; i<=n-1; i++)
        {
            v = amp::vdotproduct(state.xprec.getvector(i, n-1), state.model.getrow(i, i, n-1));
            state.x(i) = state.xbase(i)+v;
        }
    lbl_50:
    lbl_48:
        
        //
        // Composite iteration is almost over:
        // * accept new position.
        // * rebuild quadratic model
        //
        state.repiterationscount = state.repiterationscount+1;
        if( state.usermode!=lmmodefgh )
        {
            goto lbl_54;
        }
        lmclearrequestfields<Precision>(state);
        state.needfgh = true;
        state.rstate.stage = 12;
        goto lbl_rcomm;
    lbl_12:
        state.repnfunc = state.repnfunc+1;
        state.repngrad = state.repngrad+1;
        state.repnhess = state.repnhess+1;
        ablas::rmatrixcopy<Precision>(n, n, state.h, 0, 0, state.rawmodel, 0, 0);
        amp::vmove(state.gbase.getvector(0, n-1), state.g.getvector(0, n-1));
        fnew = state.f;
    lbl_54:
        if( ! (state.usermode==lmmodefgj || state.usermode==lmmodefj) )
        {
            goto lbl_56;
        }
        lmclearrequestfields<Precision>(state);
        state.needfij = true;
        state.rstate.stage = 13;
        goto lbl_rcomm;
    lbl_13:
        state.repnfunc = state.repnfunc+1;
        state.repnjac = state.repnjac+1;
        ablas::rmatrixgemm<Precision>(n, n, m, amp::ampf<Precision>("2.0"), state.j, 0, 0, 1, state.j, 0, 0, 0, amp::ampf<Precision>("0.0"), state.rawmodel, 0, 0);
        ablas::rmatrixmv<Precision>(n, m, state.j, 0, 0, 1, state.fi, 0, state.gbase, 0);
        amp::vmul(state.gbase.getvector(0, n-1), 2);
        fnew = amp::vdotproduct(state.fi.getvector(0, m-1), state.fi.getvector(0, m-1));
    lbl_56:
        
        //
        // Stopping conditions
        //
        amp::vmove(state.work.getvector(0, n-1), state.xprev.getvector(0, n-1));
        amp::vsub(state.work.getvector(0, n-1), state.x.getvector(0, n-1));
        stepnorm = amp::vdotproduct(state.work.getvector(0, n-1), state.work.getvector(0, n-1));
        stepnorm = amp::sqrt<Precision>(stepnorm);
        if( stepnorm<=state.epsx )
        {
            state.repterminationtype = 2;
            goto lbl_31;
        }
        if( state.repiterationscount>=state.maxits && state.maxits>0 )
        {
            state.repterminationtype = 5;
            goto lbl_31;
        }
        v = amp::vdotproduct(state.gbase.getvector(0, n-1), state.gbase.getvector(0, n-1));
        v = amp::sqrt<Precision>(v);
        if( v<=state.epsg )
        {
            state.repterminationtype = 4;
            goto lbl_31;
        }
        if( amp::abs<Precision>(fnew-fbase)<=state.epsf*amp::maximum<Precision>(amp::ampf<Precision>(1), amp::maximum<Precision>(amp::abs<Precision>(fnew), amp::abs<Precision>(fbase))) )
        {
            state.repterminationtype = 1;
            goto lbl_31;
        }
        
        //
        // Now, iteration is finally over:
        // * update FBase
        // * decrease lambda
        // * report new iteration
        //
        if( ! state.xrep )
        {
            goto lbl_58;
        }
        lmclearrequestfields<Precision>(state);
        state.xupdated = true;
        state.f = fnew;
        state.rstate.stage = 14;
        goto lbl_rcomm;
    lbl_14:
    lbl_58:
        fbase = fnew;
        decreaselambda<Precision>(lambda, nu, lambdadown);
        goto lbl_30;
    lbl_31:
        
        //
        // final point is reported
        //
        if( ! state.xrep )
        {
            goto lbl_60;
        }
        lmclearrequestfields<Precision>(state);
        state.xupdated = true;
        state.f = fnew;
        state.rstate.stage = 15;
        goto lbl_rcomm;
    lbl_15:
    lbl_60:
        result = false;
        return result;
        
        //
        // Saving state
        //
    lbl_rcomm:
        result = true;
        state.rstate.ia(0) = n;
        state.rstate.ia(1) = m;
        state.rstate.ia(2) = i;
        state.rstate.ia(3) = lbfgsflags;
        state.rstate.ba(0) = spd;
        state.rstate.ra(0) = stepnorm;
        state.rstate.ra(1) = fbase;
        state.rstate.ra(2) = fnew;
        state.rstate.ra(3) = lambda;
        state.rstate.ra(4) = nu;
        state.rstate.ra(5) = lambdaup;
        state.rstate.ra(6) = lambdadown;
        state.rstate.ra(7) = v;
        return result;
    }


    /*************************************************************************
    Levenberg-Marquardt algorithm results

    Called after MinLMIteration returned False.

    Input parameters:
        State   -   algorithm state (used by MinLMIteration).

    Output parameters:
        X       -   array[0..N-1], solution
        Rep     -   optimization report:
                    * Rep.TerminationType completetion code:
                        * -1    incorrect parameters were specified
                        *  1    relative function improvement is no more than
                                EpsF.
                        *  2    relative step is no more than EpsX.
                        *  4    gradient is no more than EpsG.
                        *  5    MaxIts steps was taken
                        *  7    stopping conditions are too stringent,
                                further improvement is impossible
                    * Rep.IterationsCount contains iterations count
                    * Rep.NFunc     - number of function calculations
                    * Rep.NJac      - number of Jacobi matrix calculations
                    * Rep.NGrad     - number of gradient calculations
                    * Rep.NHess     - number of Hessian calculations
                    * Rep.NCholesky - number of Cholesky decomposition calculations

      -- ALGLIB --
         Copyright 10.03.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void minlmresults(const minlmstate<Precision>& state,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        minlmreport<Precision>& rep)
    {
        x.setbounds(0, state.n-1);
        amp::vmove(x.getvector(0, state.n-1), state.x.getvector(0, state.n-1));
        rep.iterationscount = state.repiterationscount;
        rep.terminationtype = state.repterminationtype;
        rep.nfunc = state.repnfunc;
        rep.njac = state.repnjac;
        rep.ngrad = state.repngrad;
        rep.nhess = state.repnhess;
        rep.ncholesky = state.repncholesky;
    }


    /*************************************************************************
    Prepare internal structures (except for RComm).

    Note: M must be zero for FGH mode, non-zero for FJ/FGJ mode.
    *************************************************************************/
    template<unsigned int Precision>
    void lmprepare(int n,
        int m,
        bool havegrad,
        minlmstate<Precision>& state)
    {
        state.repiterationscount = 0;
        state.repterminationtype = 0;
        state.repnfunc = 0;
        state.repnjac = 0;
        state.repngrad = 0;
        state.repnhess = 0;
        state.repncholesky = 0;
        if( n<=0 || m<0 )
        {
            return;
        }
        if( havegrad )
        {
            state.g.setbounds(0, n-1);
        }
        if( m!=0 )
        {
            state.j.setbounds(0, m-1, 0, n-1);
            state.fi.setbounds(0, m-1);
            state.h.setbounds(0, 0, 0, 0);
        }
        else
        {
            state.j.setbounds(0, 0, 0, 0);
            state.fi.setbounds(0, 0);
            state.h.setbounds(0, n-1, 0, n-1);
        }
        state.x.setbounds(0, n-1);
        state.rawmodel.setbounds(0, n-1, 0, n-1);
        state.model.setbounds(0, n-1, 0, n-1);
        state.xbase.setbounds(0, n-1);
        state.xprec.setbounds(0, n-1);
        state.gbase.setbounds(0, n-1);
        state.xdir.setbounds(0, n-1);
        state.xprev.setbounds(0, n-1);
        state.work.setbounds(0, ap::maxint(n, m));
    }


    /*************************************************************************
    Clears request fileds (to be sure that we don't forgot to clear something)
    *************************************************************************/
    template<unsigned int Precision>
    void lmclearrequestfields(minlmstate<Precision>& state)
    {
        state.needf = false;
        state.needfg = false;
        state.needfgh = false;
        state.needfij = false;
        state.xupdated = false;
    }


    /*************************************************************************
    Increases lambda, returns False when there is a danger of overflow
    *************************************************************************/
    template<unsigned int Precision>
    bool increaselambda(amp::ampf<Precision>& lambda,
        amp::ampf<Precision>& nu,
        amp::ampf<Precision> lambdaup)
    {
        bool result;
        amp::ampf<Precision> lnlambda;
        amp::ampf<Precision> lnnu;
        amp::ampf<Precision> lnlambdaup;
        amp::ampf<Precision> lnmax;


        result = false;
        lnlambda = amp::log<Precision>(lambda);
        lnlambdaup = amp::log<Precision>(lambdaup);
        lnnu = amp::log<Precision>(nu);
        lnmax = amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
        if( lnlambda+lnlambdaup+lnnu>lnmax )
        {
            return result;
        }
        if( lnnu+amp::log<Precision>(amp::ampf<Precision>(2))>lnmax )
        {
            return result;
        }
        lambda = lambda*lambdaup*nu;
        nu = nu*2;
        result = true;
        return result;
    }


    /*************************************************************************
    Decreases lambda, but leaves it unchanged when there is danger of underflow.
    *************************************************************************/
    template<unsigned int Precision>
    void decreaselambda(amp::ampf<Precision>& lambda,
        amp::ampf<Precision>& nu,
        amp::ampf<Precision> lambdadown)
    {
        nu = 1;
        if( amp::log<Precision>(lambda)+amp::log<Precision>(lambdadown)<amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber()) )
        {
            lambda = amp::ampf<Precision>::getAlgoPascalMinNumber();
        }
        else
        {
            lambda = lambda*lambdadown;
        }
    }
} // namespace

#endif
