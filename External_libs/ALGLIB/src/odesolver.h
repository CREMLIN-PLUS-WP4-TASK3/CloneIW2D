/*************************************************************************
Copyright 2009 by Sergey Bochkanov (ALGLIB project).

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

#ifndef _odesolver_h
#define _odesolver_h

#include "ap.h"
#include "amp.h"
namespace odesolver
{
    template<unsigned int Precision>
    class odesolverstate
    {
    public:
        int n;
        int m;
        amp::ampf<Precision> xscale;
        amp::ampf<Precision> h;
        amp::ampf<Precision> eps;
        bool fraceps;
        ap::template_1d_array< amp::ampf<Precision> > yc;
        ap::template_1d_array< amp::ampf<Precision> > escale;
        ap::template_1d_array< amp::ampf<Precision> > xg;
        int solvertype;
        amp::ampf<Precision> x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > dy;
        ap::template_2d_array< amp::ampf<Precision> > ytbl;
        int repterminationtype;
        int repnfev;
        ap::template_1d_array< amp::ampf<Precision> > yn;
        ap::template_1d_array< amp::ampf<Precision> > yns;
        ap::template_1d_array< amp::ampf<Precision> > rka;
        ap::template_1d_array< amp::ampf<Precision> > rkc;
        ap::template_1d_array< amp::ampf<Precision> > rkcs;
        ap::template_2d_array< amp::ampf<Precision> > rkb;
        ap::template_2d_array< amp::ampf<Precision> > rkk;
        amp::rcommstate<Precision> rstate;
    };


    template<unsigned int Precision>
    class odesolverreport
    {
    public:
        int nfev;
        int terminationtype;
    };




    template<unsigned int Precision>
    void odesolverrkck(const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int m,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> h,
        odesolverstate<Precision>& state);
    template<unsigned int Precision>
    bool odesolveriteration(odesolverstate<Precision>& state);
    template<unsigned int Precision>
    void odesolverresults(const odesolverstate<Precision>& state,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& xtbl,
        ap::template_2d_array< amp::ampf<Precision> >& ytbl,
        odesolverreport<Precision>& rep);
    template<unsigned int Precision>
    void odesolverinit(int solvertype,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int m,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> h,
        odesolverstate<Precision>& state);


    template<unsigned int Precision>
    const amp::ampf<Precision>& odesolvermaxgrow()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("3.0");
        return v;
    }
    template<unsigned int Precision>
    const amp::ampf<Precision>& odesolvermaxshrink()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("10.0");
        return v;
    }


    /*************************************************************************
    Cash-Karp adaptive ODE solver.

    This subroutine solves ODE  Y'=f(Y,x)  with  initial  conditions  Y(xs)=Ys
    (here Y may be single variable or vector of N variables).

    INPUT PARAMETERS:
        Y       -   initial conditions, array[0..N-1].
                    contains values of Y[] at X[0]
        N       -   system size
        X       -   points at which Y should be tabulated, array[0..M-1]
                    integrations starts at X[0], ends at X[M-1],  intermediate
                    values at X[i] are returned too.
                    SHOULD BE ORDERED BY ASCENDING OR BY DESCENDING!!!!
        M       -   number of intermediate points + first point + last point:
                    * M>2 means that you need both Y(X[M-1]) and M-2 values at
                      intermediate points
                    * M=2 means that you want just to integrate from  X[0]  to
                      X[1] and don't interested in intermediate values.
                    * M=1 means that you don't want to integrate :)
                      it is degenerate case, but it will be handled correctly.
                    * M<1 means error
        Eps     -   tolerance (absolute/relative error on each  step  will  be
                    less than Eps). When passing:
                    * Eps>0, it means desired ABSOLUTE error
                    * Eps<0, it means desired RELATIVE error.  Relative errors
                      are calculated with respect to maximum values of  Y seen
                      so far. Be careful to use this criterion  when  starting
                      from Y[] that are close to zero.
        H       -   initial  step  lenth,  it  will  be adjusted automatically
                    after the first  step.  If  H=0,  step  will  be  selected
                    automatically  (usualy  it  will  be  equal  to  0.001  of
                    min(x[i]-x[j])).

    OUTPUT PARAMETERS
        State   -   structure which stores algorithm state between  subsequent
                    calls of OdeSolverIteration. Used for reverse communication.
                    This structure should be passed  to the OdeSolverIteration
                    subroutine.

    SEE ALSO
        AutoGKSmoothW, AutoGKSingular, AutoGKIteration, AutoGKResults.


      -- ALGLIB --
         Copyright 01.09.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void odesolverrkck(const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int m,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> h,
        odesolverstate<Precision>& state)
    {
        odesolverinit<Precision>(0, y, n, x, m, eps, h, state);
    }


    /*************************************************************************
    One iteration of ODE solver.

    Called after inialization of State structure with OdeSolverXXX subroutine.
    See HTML docs for examples.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between subsequent
                    calls and which is used for reverse communication. Must be
                    initialized with OdeSolverXXX() call first.

    If subroutine returned False, algorithm have finished its work.
    If subroutine returned True, then user should:
    * calculate F(State.X, State.Y)
    * store it in State.DY
    Here State.X is real, State.Y and State.DY are arrays[0..N-1] of reals.

      -- ALGLIB --
         Copyright 01.09.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool odesolveriteration(odesolverstate<Precision>& state)
    {
        bool result;
        int n;
        int m;
        int i;
        int j;
        int k;
        amp::ampf<Precision> xc;
        amp::ampf<Precision> v;
        amp::ampf<Precision> h;
        amp::ampf<Precision> h2;
        bool gridpoint;
        amp::ampf<Precision> err;
        amp::ampf<Precision> maxgrowpow;
        int klimit;


        
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
            j = state.rstate.ia(3);
            k = state.rstate.ia(4);
            klimit = state.rstate.ia(5);
            gridpoint = state.rstate.ba(0);
            xc = state.rstate.ra(0);
            v = state.rstate.ra(1);
            h = state.rstate.ra(2);
            h2 = state.rstate.ra(3);
            err = state.rstate.ra(4);
            maxgrowpow = state.rstate.ra(5);
        }
        else
        {
            n = -983;
            m = -989;
            i = -834;
            j = 900;
            k = -287;
            klimit = 364;
            gridpoint = false;
            xc = -338;
            v = -686;
            h = 912;
            h2 = 585;
            err = 497;
            maxgrowpow = -271;
        }
        if( state.rstate.stage==0 )
        {
            goto lbl_0;
        }
        
        //
        // Routine body
        //
        
        //
        // prepare
        //
        if( state.repterminationtype!=0 )
        {
            result = false;
            return result;
        }
        n = state.n;
        m = state.m;
        h = state.h;
        state.y.setlength(n);
        state.dy.setlength(n);
        maxgrowpow = amp::pow<Precision>(odesolvermaxgrow<Precision>(), amp::ampf<Precision>(5));
        state.repnfev = 0;
        
        //
        // some preliminary checks for internal errors
        // after this we assume that H>0 and M>1
        //
        ap::ap_error::make_assertion(state.h>0);
        ap::ap_error::make_assertion(m>1);
        
        //
        // choose solver
        //
        if( state.solvertype!=0 )
        {
            goto lbl_1;
        }
        
        //
        // Cask-Karp solver
        // Prepare coefficients table.
        // Check it for errors
        //
        state.rka.setlength(6);
        state.rka(0) = 0;
        state.rka(1) = amp::ampf<Precision>(1)/amp::ampf<Precision>(5);
        state.rka(2) = amp::ampf<Precision>(3)/amp::ampf<Precision>(10);
        state.rka(3) = amp::ampf<Precision>(3)/amp::ampf<Precision>(5);
        state.rka(4) = 1;
        state.rka(5) = amp::ampf<Precision>(7)/amp::ampf<Precision>(8);
        state.rkb.setlength(6, 5);
        state.rkb(1,0) = amp::ampf<Precision>(1)/amp::ampf<Precision>(5);
        state.rkb(2,0) = amp::ampf<Precision>(3)/amp::ampf<Precision>(40);
        state.rkb(2,1) = amp::ampf<Precision>(9)/amp::ampf<Precision>(40);
        state.rkb(3,0) = amp::ampf<Precision>(3)/amp::ampf<Precision>(10);
        state.rkb(3,1) = -amp::ampf<Precision>(9)/amp::ampf<Precision>(10);
        state.rkb(3,2) = amp::ampf<Precision>(6)/amp::ampf<Precision>(5);
        state.rkb(4,0) = -amp::ampf<Precision>(11)/amp::ampf<Precision>(54);
        state.rkb(4,1) = amp::ampf<Precision>(5)/amp::ampf<Precision>(2);
        state.rkb(4,2) = -amp::ampf<Precision>(70)/amp::ampf<Precision>(27);
        state.rkb(4,3) = amp::ampf<Precision>(35)/amp::ampf<Precision>(27);
        state.rkb(5,0) = amp::ampf<Precision>(1631)/amp::ampf<Precision>(55296);
        state.rkb(5,1) = amp::ampf<Precision>(175)/amp::ampf<Precision>(512);
        state.rkb(5,2) = amp::ampf<Precision>(575)/amp::ampf<Precision>(13824);
        state.rkb(5,3) = amp::ampf<Precision>(44275)/amp::ampf<Precision>(110592);
        state.rkb(5,4) = amp::ampf<Precision>(253)/amp::ampf<Precision>(4096);
        state.rkc.setlength(6);
        state.rkc(0) = amp::ampf<Precision>(37)/amp::ampf<Precision>(378);
        state.rkc(1) = 0;
        state.rkc(2) = amp::ampf<Precision>(250)/amp::ampf<Precision>(621);
        state.rkc(3) = amp::ampf<Precision>(125)/amp::ampf<Precision>(594);
        state.rkc(4) = 0;
        state.rkc(5) = amp::ampf<Precision>(512)/amp::ampf<Precision>(1771);
        state.rkcs.setlength(6);
        state.rkcs(0) = amp::ampf<Precision>(2825)/amp::ampf<Precision>(27648);
        state.rkcs(1) = 0;
        state.rkcs(2) = amp::ampf<Precision>(18575)/amp::ampf<Precision>(48384);
        state.rkcs(3) = amp::ampf<Precision>(13525)/amp::ampf<Precision>(55296);
        state.rkcs(4) = amp::ampf<Precision>(277)/amp::ampf<Precision>(14336);
        state.rkcs(5) = amp::ampf<Precision>(1)/amp::ampf<Precision>(4);
        state.rkk.setlength(6, n);
        
        //
        // Main cycle consists of two iterations:
        // * outer where we travel from X[i-1] to X[i]
        // * inner where we travel inside [X[i-1],X[i]]
        //
        state.ytbl.setlength(m, n);
        state.escale.setlength(n);
        state.yn.setlength(n);
        state.yns.setlength(n);
        xc = state.xg(0);
        amp::vmove(state.ytbl.getrow(0, 0, n-1), state.yc.getvector(0, n-1));
        for(j=0; j<=n-1; j++)
        {
            state.escale(j) = 0;
        }
        i = 1;
    lbl_3:
        if( i>m-1 )
        {
            goto lbl_5;
        }
        
        //
        // begin inner iteration
        //
    lbl_6:
        if( false )
        {
            goto lbl_7;
        }
        
        //
        // truncate step if needed (beyond right boundary).
        // determine should we store X or not
        //
        if( xc+h>=state.xg(i) )
        {
            h = state.xg(i)-xc;
            gridpoint = true;
        }
        else
        {
            gridpoint = false;
        }
        
        //
        // Update error scale maximums
        //
        // These maximums are initialized by zeros,
        // then updated every iterations.
        //
        for(j=0; j<=n-1; j++)
        {
            state.escale(j) = amp::maximum<Precision>(state.escale(j), amp::abs<Precision>(state.yc(j)));
        }
        
        //
        // make one step:
        // 1. calculate all info needed to do step
        // 2. update errors scale maximums using values/derivatives
        //    obtained during (1)
        //
        // Take into account that we use scaling of X to reduce task
        // to the form where x[0] < x[1] < ... < x[n-1]. So X is
        // replaced by x=xscale*t, and dy/dx=f(y,x) is replaced
        // by dy/dt=xscale*f(y,xscale*t).
        //
        amp::vmove(state.yn.getvector(0, n-1), state.yc.getvector(0, n-1));
        amp::vmove(state.yns.getvector(0, n-1), state.yc.getvector(0, n-1));
        k = 0;
    lbl_8:
        if( k>5 )
        {
            goto lbl_10;
        }
        
        //
        // prepare data for the next update of YN/YNS
        //
        state.x = state.xscale*(xc+state.rka(k)*h);
        amp::vmove(state.y.getvector(0, n-1), state.yc.getvector(0, n-1));
        for(j=0; j<=k-1; j++)
        {
            v = state.rkb(k,j);
            amp::vadd(state.y.getvector(0, n-1), state.rkk.getrow(j, 0, n-1), v);
        }
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        state.repnfev = state.repnfev+1;
        v = h*state.xscale;
        amp::vmove(state.rkk.getrow(k, 0, n-1), state.dy.getvector(0, n-1), v);
        
        //
        // update YN/YNS
        //
        v = state.rkc(k);
        amp::vadd(state.yn.getvector(0, n-1), state.rkk.getrow(k, 0, n-1), v);
        v = state.rkcs(k);
        amp::vadd(state.yns.getvector(0, n-1), state.rkk.getrow(k, 0, n-1), v);
        k = k+1;
        goto lbl_8;
    lbl_10:
        
        //
        // estimate error
        //
        err = 0;
        for(j=0; j<=n-1; j++)
        {
            if( !state.fraceps )
            {
                
                //
                // absolute error is estimated
                //
                err = amp::maximum<Precision>(err, amp::abs<Precision>(state.yn(j)-state.yns(j)));
            }
            else
            {
                
                //
                // Relative error is estimated
                //
                v = state.escale(j);
                if( v==0 )
                {
                    v = 1;
                }
                err = amp::maximum<Precision>(err, amp::abs<Precision>(state.yn(j)-state.yns(j))/v);
            }
        }
        
        //
        // calculate new step, restart if necessary
        //
        if( maxgrowpow*err<=state.eps )
        {
            h2 = odesolvermaxgrow<Precision>()*h;
        }
        else
        {
            h2 = h*amp::pow<Precision>(state.eps/err, amp::ampf<Precision>("0.2"));
        }
        if( h2<h/odesolvermaxshrink<Precision>() )
        {
            h2 = h/odesolvermaxshrink<Precision>();
        }
        if( err>state.eps )
        {
            h = h2;
            goto lbl_6;
        }
        
        //
        // advance position
        //
        xc = xc+h;
        amp::vmove(state.yc.getvector(0, n-1), state.yn.getvector(0, n-1));
        
        //
        // update H
        //
        h = h2;
        
        //
        // break on grid point
        //
        if( gridpoint )
        {
            goto lbl_7;
        }
        goto lbl_6;
    lbl_7:
        
        //
        // save result
        //
        amp::vmove(state.ytbl.getrow(i, 0, n-1), state.yc.getvector(0, n-1));
        i = i+1;
        goto lbl_3;
    lbl_5:
        state.repterminationtype = 1;
        result = false;
        return result;
    lbl_1:
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
        state.rstate.ia(3) = j;
        state.rstate.ia(4) = k;
        state.rstate.ia(5) = klimit;
        state.rstate.ba(0) = gridpoint;
        state.rstate.ra(0) = xc;
        state.rstate.ra(1) = v;
        state.rstate.ra(2) = h;
        state.rstate.ra(3) = h2;
        state.rstate.ra(4) = err;
        state.rstate.ra(5) = maxgrowpow;
        return result;
    }


    /*************************************************************************
    ODE solver results

    Called after OdeSolverIteration returned False.

    INPUT PARAMETERS:
        State   -   algorithm state (used by OdeSolverIteration).

    OUTPUT PARAMETERS:
        M       -   number of tabulated values, M>=1
        XTbl    -   array[0..M-1], values of X
        YTbl    -   array[0..M-1,0..N-1], values of Y in X[i]
        Rep     -   solver report:
                    * Rep.TerminationType completetion code:
                        * -2    X is not ordered  by  ascending/descending  or
                                there are non-distinct X[],  i.e.  X[i]=X[i+1]
                        * -1    incorrect parameters were specified
                        *  1    task has been solved
                    * Rep.NFEV contains number of function calculations

      -- ALGLIB --
         Copyright 01.09.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void odesolverresults(const odesolverstate<Precision>& state,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& xtbl,
        ap::template_2d_array< amp::ampf<Precision> >& ytbl,
        odesolverreport<Precision>& rep)
    {
        amp::ampf<Precision> v;
        int i;


        rep.terminationtype = state.repterminationtype;
        if( rep.terminationtype>0 )
        {
            m = state.m;
            rep.nfev = state.repnfev;
            xtbl.setlength(state.m);
            v = state.xscale;
            amp::vmove(xtbl.getvector(0, state.m-1), state.xg.getvector(0, state.m-1), v);
            ytbl.setlength(state.m, state.n);
            for(i=0; i<=state.m-1; i++)
            {
                amp::vmove(ytbl.getrow(i, 0, state.n-1), state.ytbl.getrow(i, 0, state.n-1));
            }
        }
        else
        {
            rep.nfev = 0;
        }
    }


    /*************************************************************************
    Internal initialization subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void odesolverinit(int solvertype,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int m,
        amp::ampf<Precision> eps,
        amp::ampf<Precision> h,
        odesolverstate<Precision>& state)
    {
        int i;
        amp::ampf<Precision> v;


        
        //
        // Prepare RComm
        //
        state.rstate.ia.setbounds(0, 5);
        state.rstate.ba.setbounds(0, 0);
        state.rstate.ra.setbounds(0, 5);
        state.rstate.stage = -1;
        
        //
        // check parameters.
        //
        if( n<=0 || m<1 || eps==0 )
        {
            state.repterminationtype = -1;
            return;
        }
        if( h<0 )
        {
            h = -h;
        }
        
        //
        // quick exit if necessary.
        // after this block we assume that M>1
        //
        if( m==1 )
        {
            state.repnfev = 0;
            state.repterminationtype = 1;
            state.ytbl.setlength(1, n);
            amp::vmove(state.ytbl.getrow(0, 0, n-1), y.getvector(0, n-1));
            state.xg.setlength(m);
            amp::vmove(state.xg.getvector(0, m-1), x.getvector(0, m-1));
            return;
        }
        
        //
        // check again: correct order of X[]
        //
        if( x(1)==x(0) )
        {
            state.repterminationtype = -2;
            return;
        }
        for(i=1; i<=m-1; i++)
        {
            if( x(1)>x(0) && x(i)<=x(i-1) || x(1)<x(0) && x(i)>=x(i-1) )
            {
                state.repterminationtype = -2;
                return;
            }
        }
        
        //
        // auto-select H if necessary
        //
        if( h==0 )
        {
            v = amp::abs<Precision>(x(1)-x(0));
            for(i=2; i<=m-1; i++)
            {
                v = amp::minimum<Precision>(v, amp::abs<Precision>(x(i)-x(i-1)));
            }
            h = amp::ampf<Precision>("0.001")*v;
        }
        
        //
        // store parameters
        //
        state.n = n;
        state.m = m;
        state.h = h;
        state.eps = amp::abs<Precision>(eps);
        state.fraceps = eps<0;
        state.xg.setlength(m);
        amp::vmove(state.xg.getvector(0, m-1), x.getvector(0, m-1));
        if( x(1)>x(0) )
        {
            state.xscale = 1;
        }
        else
        {
            state.xscale = -1;
            amp::vmul(state.xg.getvector(0, m-1), -1);
        }
        state.yc.setlength(n);
        amp::vmove(state.yc.getvector(0, n-1), y.getvector(0, n-1));
        state.solvertype = solvertype;
        state.repterminationtype = 0;
    }
} // namespace

#endif
