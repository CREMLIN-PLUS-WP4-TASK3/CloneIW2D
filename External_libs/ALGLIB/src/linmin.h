/*************************************************************************
ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE

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

#ifndef _linmin_h
#define _linmin_h

#include "ap.h"
#include "amp.h"
namespace linmin
{
    template<unsigned int Precision>
    class linminstate
    {
    public:
        bool brackt;
        bool stage1;
        int infoc;
        amp::ampf<Precision> dg;
        amp::ampf<Precision> dgm;
        amp::ampf<Precision> dginit;
        amp::ampf<Precision> dgtest;
        amp::ampf<Precision> dgx;
        amp::ampf<Precision> dgxm;
        amp::ampf<Precision> dgy;
        amp::ampf<Precision> dgym;
        amp::ampf<Precision> finit;
        amp::ampf<Precision> ftest1;
        amp::ampf<Precision> fm;
        amp::ampf<Precision> fx;
        amp::ampf<Precision> fxm;
        amp::ampf<Precision> fy;
        amp::ampf<Precision> fym;
        amp::ampf<Precision> stx;
        amp::ampf<Precision> sty;
        amp::ampf<Precision> stmin;
        amp::ampf<Precision> stmax;
        amp::ampf<Precision> width;
        amp::ampf<Precision> width1;
        amp::ampf<Precision> xtrapf;
    };




    template<unsigned int Precision>
    void linminnormalized(ap::template_1d_array< amp::ampf<Precision> >& d,
        amp::ampf<Precision>& stp,
        int n);
    template<unsigned int Precision>
    void mcsrch(const int& n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& f,
        ap::template_1d_array< amp::ampf<Precision> >& g,
        const ap::template_1d_array< amp::ampf<Precision> >& s,
        amp::ampf<Precision>& stp,
        amp::ampf<Precision> stpmax,
        int& info,
        int& nfev,
        ap::template_1d_array< amp::ampf<Precision> >& wa,
        linminstate<Precision>& state,
        int& stage);
    template<unsigned int Precision>
    void mcstep(amp::ampf<Precision>& stx,
        amp::ampf<Precision>& fx,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& sty,
        amp::ampf<Precision>& fy,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& stp,
        const amp::ampf<Precision>& fp,
        const amp::ampf<Precision>& dp,
        bool& brackt,
        const amp::ampf<Precision>& stmin,
        const amp::ampf<Precision>& stmax,
        int& info);


    template<unsigned int Precision>
    const amp::ampf<Precision>& ftol()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("0.001");
        return v;
    }
    template<unsigned int Precision>
    const amp::ampf<Precision>& xtol()
    {
        static amp::ampf<Precision> v = 100*amp::ampf<Precision>::getAlgoPascalEpsilon();
        return v;
    }
    template<unsigned int Precision>
    const amp::ampf<Precision>& gtol()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("0.3");
        return v;
    }
    static const int maxfev = 20;
    template<unsigned int Precision>
    const amp::ampf<Precision>& stpmin()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("1.0E-50");
        return v;
    }
    template<unsigned int Precision>
    const amp::ampf<Precision>& defstpmax()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("1.0E+50");
        return v;
    }


    /*************************************************************************
    Normalizes direction/step pair: makes |D|=1, scales Stp.
    If |D|=0, it returns, leavind D/Stp unchanged.

      -- ALGLIB --
         Copyright 01.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void linminnormalized(ap::template_1d_array< amp::ampf<Precision> >& d,
        amp::ampf<Precision>& stp,
        int n)
    {
        amp::ampf<Precision> mx;
        amp::ampf<Precision> s;
        int i;


        
        //
        // first, scale D to avoid underflow/overflow durng squaring
        //
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(d(i)));
        }
        if( mx==0 )
        {
            return;
        }
        s = 1/mx;
        amp::vmul(d.getvector(0, n-1), s);
        stp = stp/s;
        
        //
        // normalize D
        //
        s = amp::vdotproduct(d.getvector(0, n-1), d.getvector(0, n-1));
        s = 1/amp::sqrt<Precision>(s);
        amp::vmul(d.getvector(0, n-1), s);
        stp = stp/s;
    }


    /*************************************************************************
    THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
    DECREASE CONDITION AND A CURVATURE CONDITION.

    AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
    ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
    SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

        F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

    IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
    FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
    UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

    THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
    DECREASE CONDITION

        F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

    AND THE CURVATURE CONDITION

        ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

    IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
    BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
    IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
    ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
    IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.

    PARAMETERS DESCRIPRION

    N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

    X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
    THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

    F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
    IT CONTAINS THE VALUE OF F AT X + STP*S.

    G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
    ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

    S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

    STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
    OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

    FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
    SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
    SATISFIED.

    XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
    WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

    STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
    UPPER BOUNDS FOR THE STEP.

    MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
    NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

    INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
        INFO = 0  IMPROPER INPUT PARAMETERS.

        INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
                  DIRECTIONAL DERIVATIVE CONDITION HOLD.

        INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
                  IS AT MOST XTOL.

        INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

        INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

        INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

        INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
                  THERE MAY NOT BE A STEP WHICH SATISFIES THE
                  SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
                  TOLERANCES MAY BE TOO SMALL.

    NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

    WA IS A WORK ARRAY OF LENGTH N.

    ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
    JORGE J. MORE', DAVID J. THUENTE
    *************************************************************************/
    template<unsigned int Precision>
    void mcsrch(const int& n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& f,
        ap::template_1d_array< amp::ampf<Precision> >& g,
        const ap::template_1d_array< amp::ampf<Precision> >& s,
        amp::ampf<Precision>& stp,
        amp::ampf<Precision> stpmax,
        int& info,
        int& nfev,
        ap::template_1d_array< amp::ampf<Precision> >& wa,
        linminstate<Precision>& state,
        int& stage)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> p5;
        amp::ampf<Precision> p66;
        amp::ampf<Precision> zero;


        
        //
        // init
        //
        p5 = amp::ampf<Precision>("0.5");
        p66 = amp::ampf<Precision>("0.66");
        state.xtrapf = amp::ampf<Precision>("4.0");
        zero = 0;
        if( stpmax==0 )
        {
            stpmax = defstpmax<Precision>();
        }
        if( stp<stpmin<Precision>() )
        {
            stp = stpmin<Precision>();
        }
        if( stp>stpmax )
        {
            stp = stpmax;
        }
        
        //
        // Main cycle
        //
        while( true )
        {
            if( stage==0 )
            {
                
                //
                // NEXT
                //
                stage = 2;
                continue;
            }
            if( stage==2 )
            {
                state.infoc = 1;
                info = 0;
                
                //
                //     CHECK THE INPUT PARAMETERS FOR ERRORS.
                //
                if( n<=0 || stp<=0 || ftol<Precision>()<0 || gtol<Precision>()<zero || xtol<Precision>()<zero || stpmin<Precision>()<zero || stpmax<stpmin<Precision>() || maxfev<=0 )
                {
                    stage = 0;
                    return;
                }
                
                //
                //     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
                //     AND CHECK THAT S IS A DESCENT DIRECTION.
                //
                v = amp::vdotproduct(g.getvector(0, n-1), s.getvector(0, n-1));
                state.dginit = v;
                if( state.dginit>=0 )
                {
                    stage = 0;
                    return;
                }
                
                //
                //     INITIALIZE LOCAL VARIABLES.
                //
                state.brackt = false;
                state.stage1 = true;
                nfev = 0;
                state.finit = f;
                state.dgtest = ftol<Precision>()*state.dginit;
                state.width = stpmax-stpmin<Precision>();
                state.width1 = state.width/p5;
                amp::vmove(wa.getvector(0, n-1), x.getvector(0, n-1));
                
                //
                //     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
                //     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
                //     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
                //     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
                //     THE INTERVAL OF UNCERTAINTY.
                //     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
                //     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
                //
                state.stx = 0;
                state.fx = state.finit;
                state.dgx = state.dginit;
                state.sty = 0;
                state.fy = state.finit;
                state.dgy = state.dginit;
                
                //
                // NEXT
                //
                stage = 3;
                continue;
            }
            if( stage==3 )
            {
                
                //
                //     START OF ITERATION.
                //
                //     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
                //     TO THE PRESENT INTERVAL OF UNCERTAINTY.
                //
                if( state.brackt )
                {
                    if( state.stx<state.sty )
                    {
                        state.stmin = state.stx;
                        state.stmax = state.sty;
                    }
                    else
                    {
                        state.stmin = state.sty;
                        state.stmax = state.stx;
                    }
                }
                else
                {
                    state.stmin = state.stx;
                    state.stmax = stp+state.xtrapf*(stp-state.stx);
                }
                
                //
                //        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
                //
                if( stp>stpmax )
                {
                    stp = stpmax;
                }
                if( stp<stpmin<Precision>() )
                {
                    stp = stpmin<Precision>();
                }
                
                //
                //        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
                //        STP BE THE LOWEST POINT OBTAINED SO FAR.
                //
                if( state.brackt && (stp<=state.stmin || stp>=state.stmax) || nfev>=maxfev-1 || state.infoc==0 || state.brackt && state.stmax-state.stmin<=xtol<Precision>()*state.stmax )
                {
                    stp = state.stx;
                }
                
                //
                //        EVALUATE THE FUNCTION AND GRADIENT AT STP
                //        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
                //
                amp::vmove(x.getvector(0, n-1), wa.getvector(0, n-1));
                amp::vadd(x.getvector(0, n-1), s.getvector(0, n-1), stp);
                
                //
                // NEXT
                //
                stage = 4;
                return;
            }
            if( stage==4 )
            {
                info = 0;
                nfev = nfev+1;
                v = amp::vdotproduct(g.getvector(0, n-1), s.getvector(0, n-1));
                state.dg = v;
                state.ftest1 = state.finit+stp*state.dgtest;
                
                //
                //        TEST FOR CONVERGENCE.
                //
                if( state.brackt && (stp<=state.stmin || stp>=state.stmax) || state.infoc==0 )
                {
                    info = 6;
                }
                if( stp==stpmax && f<=state.ftest1 && state.dg<=state.dgtest )
                {
                    info = 5;
                }
                if( stp==stpmin<Precision>() && (f>state.ftest1 || state.dg>=state.dgtest) )
                {
                    info = 4;
                }
                if( nfev>=maxfev )
                {
                    info = 3;
                }
                if( state.brackt && state.stmax-state.stmin<=xtol<Precision>()*state.stmax )
                {
                    info = 2;
                }
                if( f<=state.ftest1 && amp::abs<Precision>(state.dg)<=-gtol<Precision>()*state.dginit )
                {
                    info = 1;
                }
                
                //
                //        CHECK FOR TERMINATION.
                //
                if( info!=0 )
                {
                    stage = 0;
                    return;
                }
                
                //
                //        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
                //        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
                //
                if( state.stage1 && f<=state.ftest1 && state.dg>=amp::minimum<Precision>(ftol<Precision>(), gtol<Precision>())*state.dginit )
                {
                    state.stage1 = false;
                }
                
                //
                //        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
                //        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
                //        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
                //        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
                //        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
                //
                if( state.stage1 && f<=state.fx && f>state.ftest1 )
                {
                    
                    //
                    //           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                    //
                    state.fm = f-stp*state.dgtest;
                    state.fxm = state.fx-state.stx*state.dgtest;
                    state.fym = state.fy-state.sty*state.dgtest;
                    state.dgm = state.dg-state.dgtest;
                    state.dgxm = state.dgx-state.dgtest;
                    state.dgym = state.dgy-state.dgtest;
                    
                    //
                    //           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                    //           AND TO COMPUTE THE NEW STEP.
                    //
                    mcstep<Precision>(state.stx, state.fxm, state.dgxm, state.sty, state.fym, state.dgym, stp, state.fm, state.dgm, state.brackt, state.stmin, state.stmax, state.infoc);
                    
                    //
                    //           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                    //
                    state.fx = state.fxm+state.stx*state.dgtest;
                    state.fy = state.fym+state.sty*state.dgtest;
                    state.dgx = state.dgxm+state.dgtest;
                    state.dgy = state.dgym+state.dgtest;
                }
                else
                {
                    
                    //
                    //           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                    //           AND TO COMPUTE THE NEW STEP.
                    //
                    mcstep<Precision>(state.stx, state.fx, state.dgx, state.sty, state.fy, state.dgy, stp, f, state.dg, state.brackt, state.stmin, state.stmax, state.infoc);
                }
                
                //
                //        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
                //        INTERVAL OF UNCERTAINTY.
                //
                if( state.brackt )
                {
                    if( amp::abs<Precision>(state.sty-state.stx)>=p66*state.width1 )
                    {
                        stp = state.stx+p5*(state.sty-state.stx);
                    }
                    state.width1 = state.width;
                    state.width = amp::abs<Precision>(state.sty-state.stx);
                }
                
                //
                //  NEXT.
                //
                stage = 3;
                continue;
            }
        }
    }


    template<unsigned int Precision>
    void mcstep(amp::ampf<Precision>& stx,
        amp::ampf<Precision>& fx,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& sty,
        amp::ampf<Precision>& fy,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& stp,
        const amp::ampf<Precision>& fp,
        const amp::ampf<Precision>& dp,
        bool& brackt,
        const amp::ampf<Precision>& stmin,
        const amp::ampf<Precision>& stmax,
        int& info)
    {
        bool bound;
        amp::ampf<Precision> gamma;
        amp::ampf<Precision> p;
        amp::ampf<Precision> q;
        amp::ampf<Precision> r;
        amp::ampf<Precision> s;
        amp::ampf<Precision> sgnd;
        amp::ampf<Precision> stpc;
        amp::ampf<Precision> stpf;
        amp::ampf<Precision> stpq;
        amp::ampf<Precision> theta;


        info = 0;
        
        //
        //     CHECK THE INPUT PARAMETERS FOR ERRORS.
        //
        if( brackt && (stp<=amp::minimum<Precision>(stx, sty) || stp>=amp::maximum<Precision>(stx, sty)) || dx*(stp-stx)>=0 || stmax<stmin )
        {
            return;
        }
        
        //
        //     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
        //
        sgnd = dp*(dx/amp::abs<Precision>(dx));
        
        //
        //     FIRST CASE. A HIGHER FUNCTION VALUE.
        //     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
        //     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
        //     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
        //
        if( fp>fx )
        {
            info = 1;
            bound = true;
            theta = 3*(fx-fp)/(stp-stx)+dx+dp;
            s = amp::maximum<Precision>(amp::abs<Precision>(theta), amp::maximum<Precision>(amp::abs<Precision>(dx), amp::abs<Precision>(dp)));
            gamma = s*amp::sqrt<Precision>(amp::sqr<Precision>(theta/s)-dx/s*(dp/s));
            if( stp<stx )
            {
                gamma = -gamma;
            }
            p = gamma-dx+theta;
            q = gamma-dx+gamma+dp;
            r = p/q;
            stpc = stx+r*(stp-stx);
            stpq = stx+dx/((fx-fp)/(stp-stx)+dx)/2*(stp-stx);
            if( amp::abs<Precision>(stpc-stx)<amp::abs<Precision>(stpq-stx) )
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpc+(stpq-stpc)/2;
            }
            brackt = true;
        }
        else
        {
            if( sgnd<0 )
            {
                
                //
                //     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
                //     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
                //     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
                //     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
                //
                info = 2;
                bound = false;
                theta = 3*(fx-fp)/(stp-stx)+dx+dp;
                s = amp::maximum<Precision>(amp::abs<Precision>(theta), amp::maximum<Precision>(amp::abs<Precision>(dx), amp::abs<Precision>(dp)));
                gamma = s*amp::sqrt<Precision>(amp::sqr<Precision>(theta/s)-dx/s*(dp/s));
                if( stp>stx )
                {
                    gamma = -gamma;
                }
                p = gamma-dp+theta;
                q = gamma-dp+gamma+dx;
                r = p/q;
                stpc = stp+r*(stx-stp);
                stpq = stp+dp/(dp-dx)*(stx-stp);
                if( amp::abs<Precision>(stpc-stp)>amp::abs<Precision>(stpq-stp) )
                {
                    stpf = stpc;
                }
                else
                {
                    stpf = stpq;
                }
                brackt = true;
            }
            else
            {
                if( amp::abs<Precision>(dp)<amp::abs<Precision>(dx) )
                {
                    
                    //
                    //     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                    //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                    //     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                    //     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                    //     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                    //     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                    //     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                    //     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                    //
                    info = 3;
                    bound = true;
                    theta = 3*(fx-fp)/(stp-stx)+dx+dp;
                    s = amp::maximum<Precision>(amp::abs<Precision>(theta), amp::maximum<Precision>(amp::abs<Precision>(dx), amp::abs<Precision>(dp)));
                    
                    //
                    //        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                    //        TO INFINITY IN THE DIRECTION OF THE STEP.
                    //
                    gamma = s*amp::sqrt<Precision>(amp::maximum<Precision>(amp::ampf<Precision>(0), amp::sqr<Precision>(theta/s)-dx/s*(dp/s)));
                    if( stp>stx )
                    {
                        gamma = -gamma;
                    }
                    p = gamma-dp+theta;
                    q = gamma+(dx-dp)+gamma;
                    r = p/q;
                    if( r<0 && gamma!=0 )
                    {
                        stpc = stp+r*(stx-stp);
                    }
                    else
                    {
                        if( stp>stx )
                        {
                            stpc = stmax;
                        }
                        else
                        {
                            stpc = stmin;
                        }
                    }
                    stpq = stp+dp/(dp-dx)*(stx-stp);
                    if( brackt )
                    {
                        if( amp::abs<Precision>(stp-stpc)<amp::abs<Precision>(stp-stpq) )
                        {
                            stpf = stpc;
                        }
                        else
                        {
                            stpf = stpq;
                        }
                    }
                    else
                    {
                        if( amp::abs<Precision>(stp-stpc)>amp::abs<Precision>(stp-stpq) )
                        {
                            stpf = stpc;
                        }
                        else
                        {
                            stpf = stpq;
                        }
                    }
                }
                else
                {
                    
                    //
                    //     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                    //     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                    //     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                    //     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                    //
                    info = 4;
                    bound = false;
                    if( brackt )
                    {
                        theta = 3*(fp-fy)/(sty-stp)+dy+dp;
                        s = amp::maximum<Precision>(amp::abs<Precision>(theta), amp::maximum<Precision>(amp::abs<Precision>(dy), amp::abs<Precision>(dp)));
                        gamma = s*amp::sqrt<Precision>(amp::sqr<Precision>(theta/s)-dy/s*(dp/s));
                        if( stp>sty )
                        {
                            gamma = -gamma;
                        }
                        p = gamma-dp+theta;
                        q = gamma-dp+gamma+dy;
                        r = p/q;
                        stpc = stp+r*(sty-stp);
                        stpf = stpc;
                    }
                    else
                    {
                        if( stp>stx )
                        {
                            stpf = stmax;
                        }
                        else
                        {
                            stpf = stmin;
                        }
                    }
                }
            }
        }
        
        //
        //     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
        //     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
        //
        if( fp>fx )
        {
            sty = stp;
            fy = fp;
            dy = dp;
        }
        else
        {
            if( sgnd<amp::ampf<Precision>("0.0") )
            {
                sty = stx;
                fy = fx;
                dy = dx;
            }
            stx = stp;
            fx = fp;
            dx = dp;
        }
        
        //
        //     COMPUTE THE NEW STEP AND SAFEGUARD IT.
        //
        stpf = amp::minimum<Precision>(stmax, stpf);
        stpf = amp::maximum<Precision>(stmin, stpf);
        stp = stpf;
        if( brackt && bound )
        {
            if( sty>stx )
            {
                stp = amp::minimum<Precision>(stx+amp::ampf<Precision>("0.66")*(sty-stx), stp);
            }
            else
            {
                stp = amp::maximum<Precision>(stx+amp::ampf<Precision>("0.66")*(sty-stx), stp);
            }
        }
    }
} // namespace

#endif
