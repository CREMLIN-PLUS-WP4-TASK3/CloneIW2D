/*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _mincg_h
#define _mincg_h

#include "ap.h"
#include "amp.h"
#include "linmin.h"
namespace mincg
{
    template<unsigned int Precision>
    class mincgstate
    {
    public:
        int n;
        amp::ampf<Precision> epsg;
        amp::ampf<Precision> epsf;
        amp::ampf<Precision> epsx;
        int maxits;
        amp::ampf<Precision> stpmax;
        bool xrep;
        int cgtype;
        int nfev;
        int mcstage;
        int k;
        ap::template_1d_array< amp::ampf<Precision> > xk;
        ap::template_1d_array< amp::ampf<Precision> > dk;
        ap::template_1d_array< amp::ampf<Precision> > xn;
        ap::template_1d_array< amp::ampf<Precision> > dn;
        ap::template_1d_array< amp::ampf<Precision> > d;
        amp::ampf<Precision> fold;
        amp::ampf<Precision> stp;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_1d_array< amp::ampf<Precision> > yk;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> f;
        ap::template_1d_array< amp::ampf<Precision> > g;
        bool needfg;
        bool xupdated;
        amp::rcommstate<Precision> rstate;
        int repiterationscount;
        int repnfev;
        int repterminationtype;
        int debugrestartscount;
        linmin::linminstate<Precision> lstate;
        amp::ampf<Precision> betahs;
        amp::ampf<Precision> betady;
    };


    template<unsigned int Precision>
    class mincgreport
    {
    public:
        int iterationscount;
        int nfev;
        int terminationtype;
    };




    template<unsigned int Precision>
    void mincgcreate(int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        mincgstate<Precision>& state);
    template<unsigned int Precision>
    void mincgsetcond(mincgstate<Precision>& state,
        amp::ampf<Precision> epsg,
        amp::ampf<Precision> epsf,
        amp::ampf<Precision> epsx,
        int maxits);
    template<unsigned int Precision>
    void mincgsetxrep(mincgstate<Precision>& state,
        bool needxrep);
    template<unsigned int Precision>
    void mincgsetcgtype(mincgstate<Precision>& state,
        int cgtype);
    template<unsigned int Precision>
    void mincgsetstpmax(mincgstate<Precision>& state,
        amp::ampf<Precision> stpmax);
    template<unsigned int Precision>
    bool mincgiteration(mincgstate<Precision>& state);
    template<unsigned int Precision>
    void mincgresults(const mincgstate<Precision>& state,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        mincgreport<Precision>& rep);
    template<unsigned int Precision>
    void clearrequestfields(mincgstate<Precision>& state);


    /*************************************************************************
            NONLINEAR CONJUGATE GRADIENT METHOD

    The subroutine minimizes function F(x) of N arguments by using one of  the
    nonlinear conjugate gradient methods.

    These CG methods are globally convergent (even on non-convex functions) as
    long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
    L = { x : f(x)<=f(x0) }.

    INPUT PARAMETERS:
        N       -   problem dimension. N>0
        X       -   initial solution approximation, array[0..N-1].
        EpsG    -   positive number which  defines  a  precision  of  search.  The
                    subroutine finishes its work if the condition ||G|| < EpsG  is
                    satisfied, where ||.|| means Euclidian norm, G - gradient, X -
                    current approximation.
        EpsF    -   positive number which  defines  a  precision  of  search.  The
                    subroutine finishes its work if on iteration  number  k+1  the
                    condition |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}    is
                    satisfied.
        EpsX    -   positive number which  defines  a  precision  of  search.  The
                    subroutine finishes its work if on iteration number k+1    the
                    condition |X(k+1)-X(k)| <= EpsX is fulfilled.
        MaxIts  -   maximum number of iterations. If MaxIts=0, the number of
                    iterations is unlimited.

    OUTPUT PARAMETERS:
        State - structure used for reverse communication.

    See also MinCGIteration, MinCGResults

    NOTE:

    Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
    automatic stopping criterion selection (small EpsX).

      -- ALGLIB --
         Copyright 25.03.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgcreate(int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        mincgstate<Precision>& state)
    {
        ap::ap_error::make_assertion(n>=1);
        
        //
        // Initialize
        //
        state.n = n;
        mincgsetcond<Precision>(state, amp::ampf<Precision>(0), amp::ampf<Precision>(0), amp::ampf<Precision>(0), 0);
        mincgsetxrep<Precision>(state, false);
        mincgsetstpmax<Precision>(state, amp::ampf<Precision>(0));
        mincgsetcgtype<Precision>(state, -1);
        state.xk.setlength(n);
        state.dk.setlength(n);
        state.xn.setlength(n);
        state.dn.setlength(n);
        state.x.setlength(n);
        state.d.setlength(n);
        state.g.setlength(n);
        state.work.setlength(n);
        state.yk.setlength(n);
        
        //
        // Prepare first run
        //
        amp::vmove(state.x.getvector(0, n-1), x.getvector(0, n-1));
        state.rstate.ia.setbounds(0, 2);
        state.rstate.ra.setbounds(0, 2);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    This function sets stopping conditions for CG optimization algorithm.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be initialized
                    with MinCGCreate()
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
                    iterations is unlimited.

    Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
    automatic stopping criterion selection (small EpsX).

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgsetcond(mincgstate<Precision>& state,
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
                    initialized with MinCGCreate()
        NeedXRep-   whether iteration reports are needed or not

    Usually  algorithm  returns  from  MinCGIteration()  only  when  it  needs
    function/gradient. However, with this function we can let  it  stop  after
    each  iteration  (one  iteration  may  include   more  than  one  function
    evaluation), which is indicated by XUpdated field.

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgsetxrep(mincgstate<Precision>& state,
        bool needxrep)
    {
        state.xrep = needxrep;
    }


    /*************************************************************************
    This function sets CG algorithm.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be
                    initialized with MinCGCreate()
        CGType  -   algorithm type:
                    * -1    automatic selection of the best algorithm
                    * 0     DY (Dai and Yuan) algorithm
                    * 1     Hybrid DY-HS algorithm

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgsetcgtype(mincgstate<Precision>& state,
        int cgtype)
    {
        ap::ap_error::make_assertion(cgtype>=-1 && cgtype<=1);
        if( cgtype==-1 )
        {
            cgtype = 1;
        }
        state.cgtype = cgtype;
    }


    /*************************************************************************
    This function sets maximum step length

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be
                    initialized with MinCGCreate()
        StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                    want to limit step length.

    Use this subroutine when you optimize target function which contains exp()
    or  other  fast  growing  functions,  and optimization algorithm makes too
    large  steps  which  leads  to overflow. This function allows us to reject
    steps  that  are  too  large  (and  therefore  expose  us  to the possible
    overflow) without actually calculating function value at the x+stp*d.

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgsetstpmax(mincgstate<Precision>& state,
        amp::ampf<Precision> stpmax)
    {
        ap::ap_error::make_assertion(stpmax>=0);
        state.stpmax = stpmax;
    }


    /*************************************************************************
    One conjugate gradient iteration

    Called after initialization with MinCG.
    See HTML documentation for examples.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be initialized
                    with MinCG.

    RESULT:
    * if function returned False, iterative proces has converged.
      Use MinLBFGSResults() to obtain optimization results.
    * if subroutine returned True, then, depending on structure fields, we
      have one of the following situations


    === FUNC/GRAD REQUEST ===
    State.NeedFG is True => function value/gradient are needed.
    Caller should calculate function value State.F and gradient
    State.G[0..N-1] at State.X[0..N-1] and call MinLBFGSIteration() again.

    === NEW INTERATION IS REPORTED ===
    State.XUpdated is True => one more iteration was made.
    State.X contains current position, State.F contains function value at X.
    You can read info from these fields, but never modify  them  because  they
    contain the only copy of optimization algorithm state.

    One and only one of these fields (NeedFG, XUpdated) is true on return. New
    iterations are reported only when reports  are  explicitly  turned  on  by
    MinLBFGSSetXRep() function, so if you never called it, you can expect that
    NeedFG is always True.


      -- ALGLIB --
         Copyright 20.04.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool mincgiteration(mincgstate<Precision>& state)
    {
        bool result;
        int n;
        int i;
        amp::ampf<Precision> betak;
        amp::ampf<Precision> v;
        amp::ampf<Precision> vv;
        int mcinfo;


        
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
            i = state.rstate.ia(1);
            mcinfo = state.rstate.ia(2);
            betak = state.rstate.ra(0);
            v = state.rstate.ra(1);
            vv = state.rstate.ra(2);
        }
        else
        {
            n = -983;
            i = -989;
            mcinfo = -834;
            betak = 900;
            v = -287;
            vv = 364;
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
        
        //
        // Routine body
        //
        
        //
        // Prepare
        //
        n = state.n;
        state.repterminationtype = 0;
        state.repiterationscount = 0;
        state.repnfev = 0;
        state.debugrestartscount = 0;
        
        //
        // Calculate F/G, initialize algorithm
        //
        clearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        if( ! state.xrep )
        {
            goto lbl_4;
        }
        clearrequestfields<Precision>(state);
        state.xupdated = true;
        state.rstate.stage = 1;
        goto lbl_rcomm;
    lbl_1:
    lbl_4:
        v = amp::vdotproduct(state.g.getvector(0, n-1), state.g.getvector(0, n-1));
        v = amp::sqrt<Precision>(v);
        if( v==0 )
        {
            state.repterminationtype = 4;
            result = false;
            return result;
        }
        state.repnfev = 1;
        state.k = 0;
        state.fold = state.f;
        amp::vmove(state.xk.getvector(0, n-1), state.x.getvector(0, n-1));
        amp::vmoveneg(state.dk.getvector(0, n-1), state.g.getvector(0, n-1));
        
        //
        // Main cycle
        //
    lbl_6:
        if( false )
        {
            goto lbl_7;
        }
        
        //
        // Store G[k] for later calculation of Y[k]
        //
        amp::vmoveneg(state.yk.getvector(0, n-1), state.g.getvector(0, n-1));
        
        //
        // Calculate X(k+1): minimize F(x+alpha*d)
        //
        amp::vmove(state.d.getvector(0, n-1), state.dk.getvector(0, n-1));
        amp::vmove(state.x.getvector(0, n-1), state.xk.getvector(0, n-1));
        state.mcstage = 0;
        state.stp = amp::ampf<Precision>("1.0");
        linmin::linminnormalized<Precision>(state.d, state.stp, n);
        linmin::mcsrch<Precision>(n, state.x, state.f, state.g, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
    lbl_8:
        if( state.mcstage==0 )
        {
            goto lbl_9;
        }
        clearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 2;
        goto lbl_rcomm;
    lbl_2:
        linmin::mcsrch<Precision>(n, state.x, state.f, state.g, state.d, state.stp, state.stpmax, mcinfo, state.nfev, state.work, state.lstate, state.mcstage);
        goto lbl_8;
    lbl_9:
        if( ! state.xrep )
        {
            goto lbl_10;
        }
        clearrequestfields<Precision>(state);
        state.xupdated = true;
        state.rstate.stage = 3;
        goto lbl_rcomm;
    lbl_3:
    lbl_10:
        amp::vmove(state.xn.getvector(0, n-1), state.x.getvector(0, n-1));
        if( mcinfo==1 )
        {
            
            //
            // Standard Wolfe conditions hold
            // Calculate Y[K] and BetaK
            //
            amp::vadd(state.yk.getvector(0, n-1), state.g.getvector(0, n-1));
            vv = amp::vdotproduct(state.yk.getvector(0, n-1), state.dk.getvector(0, n-1));
            v = amp::vdotproduct(state.g.getvector(0, n-1), state.g.getvector(0, n-1));
            state.betady = v/vv;
            v = amp::vdotproduct(state.g.getvector(0, n-1), state.yk.getvector(0, n-1));
            state.betahs = v/vv;
            if( state.cgtype==0 )
            {
                betak = state.betady;
            }
            if( state.cgtype==1 )
            {
                betak = amp::maximum<Precision>(amp::ampf<Precision>(0), amp::minimum<Precision>(state.betady, state.betahs));
            }
        }
        else
        {
            
            //
            // Something is wrong (may be function is too wild or too flat).
            //
            // We'll set BetaK=0, which will restart CG algorithm.
            // We can stop later (during normal checks) if stopping conditions are met.
            //
            betak = 0;
            state.debugrestartscount = state.debugrestartscount+1;
        }
        
        //
        // Calculate D(k+1)
        //
        amp::vmoveneg(state.dn.getvector(0, n-1), state.g.getvector(0, n-1));
        amp::vadd(state.dn.getvector(0, n-1), state.dk.getvector(0, n-1), betak);
        
        //
        // Update information and Hessian.
        // Check stopping conditions.
        //
        state.repnfev = state.repnfev+state.nfev;
        state.repiterationscount = state.repiterationscount+1;
        if( state.repiterationscount>=state.maxits && state.maxits>0 )
        {
            
            //
            // Too many iterations
            //
            state.repterminationtype = 5;
            result = false;
            return result;
        }
        v = amp::vdotproduct(state.g.getvector(0, n-1), state.g.getvector(0, n-1));
        if( amp::sqrt<Precision>(v)<=state.epsg )
        {
            
            //
            // Gradient is small enough
            //
            state.repterminationtype = 4;
            result = false;
            return result;
        }
        if( state.fold-state.f<=state.epsf*amp::maximum<Precision>(amp::abs<Precision>(state.fold), amp::maximum<Precision>(amp::abs<Precision>(state.f), amp::ampf<Precision>("1.0"))) )
        {
            
            //
            // F(k+1)-F(k) is small enough
            //
            state.repterminationtype = 1;
            result = false;
            return result;
        }
        v = amp::vdotproduct(state.d.getvector(0, n-1), state.d.getvector(0, n-1));
        if( amp::sqrt<Precision>(v)*state.stp<=state.epsx )
        {
            
            //
            // X(k+1)-X(k) is small enough
            //
            state.repterminationtype = 2;
            result = false;
            return result;
        }
        
        //
        // Shift Xk/Dk, update other information
        //
        amp::vmove(state.xk.getvector(0, n-1), state.xn.getvector(0, n-1));
        amp::vmove(state.dk.getvector(0, n-1), state.dn.getvector(0, n-1));
        state.fold = state.f;
        state.k = state.k+1;
        goto lbl_6;
    lbl_7:
        result = false;
        return result;
        
        //
        // Saving state
        //
    lbl_rcomm:
        result = true;
        state.rstate.ia(0) = n;
        state.rstate.ia(1) = i;
        state.rstate.ia(2) = mcinfo;
        state.rstate.ra(0) = betak;
        state.rstate.ra(1) = v;
        state.rstate.ra(2) = vv;
        return result;
    }


    /*************************************************************************
    Conjugate gradient results

    Called after MinCG returned False.

    INPUT PARAMETERS:
        State   -   algorithm state (used by MinCGIteration).

    OUTPUT PARAMETERS:
        X       -   array[0..N-1], solution
        Rep     -   optimization report:
                    * Rep.TerminationType completetion code:
                        * -2    rounding errors prevent further improvement.
                                X contains best point found.
                        * -1    incorrect parameters were specified
                        *  1    relative function improvement is no more than
                                EpsF.
                        *  2    relative step is no more than EpsX.
                        *  4    gradient norm is no more than EpsG
                        *  5    MaxIts steps was taken
                        *  7    stopping conditions are too stringent,
                                further improvement is impossible
                    * Rep.IterationsCount contains iterations count
                    * NFEV countains number of function calculations

      -- ALGLIB --
         Copyright 20.04.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void mincgresults(const mincgstate<Precision>& state,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        mincgreport<Precision>& rep)
    {
        x.setbounds(0, state.n-1);
        amp::vmove(x.getvector(0, state.n-1), state.xn.getvector(0, state.n-1));
        rep.iterationscount = state.repiterationscount;
        rep.nfev = state.repnfev;
        rep.terminationtype = state.repterminationtype;
    }


    /*************************************************************************
    Clears request fileds (to be sure that we don't forgot to clear something)
    *************************************************************************/
    template<unsigned int Precision>
    void clearrequestfields(mincgstate<Precision>& state)
    {
        state.needfg = false;
        state.xupdated = false;
    }
} // namespace

#endif
