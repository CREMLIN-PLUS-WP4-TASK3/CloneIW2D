/*************************************************************************
Copyright (c) 2007-2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _ratint_h
#define _ratint_h

#include "ap.h"
#include "amp.h"
#include "tsort.h"
#include "ratinterpolation.h"
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
#include "minlm.h"
#include "lsfit.h"
namespace ratint
{
    /*************************************************************************
    Barycentric interpolant.
    *************************************************************************/
    template<unsigned int Precision>
    class barycentricinterpolant
    {
    public:
        int n;
        amp::ampf<Precision> sy;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > w;
    };


    /*************************************************************************
    Barycentric fitting report:
        TaskRCond       reciprocal of task's condition number
        RMSError        RMS error
        AvgError        average error
        AvgRelError     average relative error (for non-zero Y[I])
        MaxError        maximum error
    *************************************************************************/
    template<unsigned int Precision>
    class barycentricfitreport
    {
    public:
        amp::ampf<Precision> taskrcond;
        int dbest;
        amp::ampf<Precision> rmserror;
        amp::ampf<Precision> avgerror;
        amp::ampf<Precision> avgrelerror;
        amp::ampf<Precision> maxerror;
    };




    template<unsigned int Precision>
    amp::ampf<Precision> barycentriccalc(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t);
    template<unsigned int Precision>
    void barycentricdiff1(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& df);
    template<unsigned int Precision>
    void barycentricdiff2(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& df,
        amp::ampf<Precision>& d2f);
    template<unsigned int Precision>
    void barycentriclintransx(barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> ca,
        amp::ampf<Precision> cb);
    template<unsigned int Precision>
    void barycentriclintransy(barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> ca,
        amp::ampf<Precision> cb);
    template<unsigned int Precision>
    void barycentricunpack(const barycentricinterpolant<Precision>& b,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& w);
    template<unsigned int Precision>
    void barycentricserialize(const barycentricinterpolant<Precision>& b,
        ap::template_1d_array< amp::ampf<Precision> >& ra,
        int& ralen);
    template<unsigned int Precision>
    void barycentricunserialize(const ap::template_1d_array< amp::ampf<Precision> >& ra,
        barycentricinterpolant<Precision>& b);
    template<unsigned int Precision>
    void barycentriccopy(const barycentricinterpolant<Precision>& b,
        barycentricinterpolant<Precision>& b2);
    template<unsigned int Precision>
    void barycentricbuildxyw(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        barycentricinterpolant<Precision>& b);
    template<unsigned int Precision>
    void barycentricbuildfloaterhormann(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int d,
        barycentricinterpolant<Precision>& b);
    template<unsigned int Precision>
    void barycentricfitfloaterhormannwc(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& xc,
        const ap::template_1d_array< amp::ampf<Precision> >& yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep);
    template<unsigned int Precision>
    void barycentricfitfloaterhormann(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int m,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep);
    template<unsigned int Precision>
    void barycentricnormalize(barycentricinterpolant<Precision>& b);
    template<unsigned int Precision>
    void barycentriccalcbasis(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    void barycentricfitwcfixedd(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        ap::template_1d_array< amp::ampf<Precision> > xc,
        ap::template_1d_array< amp::ampf<Precision> > yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int d,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep);


    static const int brcvnum = 10;


    /*************************************************************************
    Rational interpolation using barycentric formula

    F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

    Input parameters:
        B   -   barycentric interpolant built with one of model building
                subroutines.
        T   -   interpolation point

    Result:
        barycentric interpolant F(t)

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> barycentriccalc(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;
        int i;


        
        //
        // special case: N=1
        //
        if( b.n==1 )
        {
            result = b.sy*b.y(0);
            return result;
        }
        
        //
        // Here we assume that task is normalized, i.e.:
        // 1. abs(Y[i])<=1
        // 2. abs(W[i])<=1
        // 3. X[] is ordered
        //
        s = amp::abs<Precision>(t-b.x(0));
        for(i=0; i<=b.n-1; i++)
        {
            v = b.x(i);
            if( v==t )
            {
                result = b.sy*b.y(i);
                return result;
            }
            v = amp::abs<Precision>(t-v);
            if( v<s )
            {
                s = v;
            }
        }
        s1 = 0;
        s2 = 0;
        for(i=0; i<=b.n-1; i++)
        {
            v = s/(t-b.x(i));
            v = v*b.w(i);
            s1 = s1+v*b.y(i);
            s2 = s2+v;
        }
        result = b.sy*s1/s2;
        return result;
    }


    /*************************************************************************
    Differentiation of barycentric interpolant: first derivative.

    Algorithm used in this subroutine is very robust and should not fail until
    provided with values too close to MaxRealNumber  (usually  MaxRealNumber/N
    or greater will overflow).

    INPUT PARAMETERS:
        B   -   barycentric interpolant built with one of model building
                subroutines.
        T   -   interpolation point

    OUTPUT PARAMETERS:
        F   -   barycentric interpolant at T
        DF  -   first derivative
        
    NOTE


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricdiff1(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& df)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> vv;
        int i;
        int k;
        amp::ampf<Precision> n0;
        amp::ampf<Precision> n1;
        amp::ampf<Precision> d0;
        amp::ampf<Precision> d1;
        amp::ampf<Precision> s0;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> xk;
        amp::ampf<Precision> xi;
        amp::ampf<Precision> xmin;
        amp::ampf<Precision> xmax;
        amp::ampf<Precision> xscale1;
        amp::ampf<Precision> xoffs1;
        amp::ampf<Precision> xscale2;
        amp::ampf<Precision> xoffs2;
        amp::ampf<Precision> xprev;


        
        //
        // special case: N=1
        //
        if( b.n==1 )
        {
            f = b.sy*b.y(0);
            df = 0;
            return;
        }
        if( b.sy==0 )
        {
            f = 0;
            df = 0;
            return;
        }
        ap::ap_error::make_assertion(b.sy>0);
        
        //
        // We assume than N>1 and B.SY>0. Find:
        // 1. pivot point (X[i] closest to T)
        // 2. width of interval containing X[i]
        //
        v = amp::abs<Precision>(b.x(0)-t);
        k = 0;
        xmin = b.x(0);
        xmax = b.x(0);
        for(i=1; i<=b.n-1; i++)
        {
            vv = b.x(i);
            if( amp::abs<Precision>(vv-t)<v )
            {
                v = amp::abs<Precision>(vv-t);
                k = i;
            }
            xmin = amp::minimum<Precision>(xmin, vv);
            xmax = amp::maximum<Precision>(xmax, vv);
        }
        
        //
        // pivot point found, calculate dNumerator and dDenominator
        //
        xscale1 = 1/(xmax-xmin);
        xoffs1 = -xmin/(xmax-xmin)+1;
        xscale2 = 2;
        xoffs2 = -3;
        t = t*xscale1+xoffs1;
        t = t*xscale2+xoffs2;
        xk = b.x(k);
        xk = xk*xscale1+xoffs1;
        xk = xk*xscale2+xoffs2;
        v = t-xk;
        n0 = 0;
        n1 = 0;
        d0 = 0;
        d1 = 0;
        xprev = -2;
        for(i=0; i<=b.n-1; i++)
        {
            xi = b.x(i);
            xi = xi*xscale1+xoffs1;
            xi = xi*xscale2+xoffs2;
            ap::ap_error::make_assertion(xi>xprev);
            xprev = xi;
            if( i!=k )
            {
                vv = amp::sqr<Precision>(t-xi);
                s0 = (t-xk)/(t-xi);
                s1 = (xk-xi)/vv;
            }
            else
            {
                s0 = 1;
                s1 = 0;
            }
            vv = b.w(i)*b.y(i);
            n0 = n0+s0*vv;
            n1 = n1+s1*vv;
            vv = b.w(i);
            d0 = d0+s0*vv;
            d1 = d1+s1*vv;
        }
        f = b.sy*n0/d0;
        df = (n1*d0-n0*d1)/amp::sqr<Precision>(d0);
        if( df!=0 )
        {
            df = amp::sign<Precision>(df)*amp::exp<Precision>(amp::log<Precision>(amp::abs<Precision>(df))+amp::log<Precision>(b.sy)+amp::log<Precision>(xscale1)+amp::log<Precision>(xscale2));
        }
    }


    /*************************************************************************
    Differentiation of barycentric interpolant: first/second derivatives.

    INPUT PARAMETERS:
        B   -   barycentric interpolant built with one of model building
                subroutines.
        T   -   interpolation point

    OUTPUT PARAMETERS:
        F   -   barycentric interpolant at T
        DF  -   first derivative
        D2F -   second derivative

    NOTE: this algorithm may fail due to overflow/underflor if  used  on  data
    whose values are close to MaxRealNumber or MinRealNumber.  Use more robust
    BarycentricDiff1() subroutine in such cases.


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricdiff2(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& f,
        amp::ampf<Precision>& df,
        amp::ampf<Precision>& d2f)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> vv;
        int i;
        int k;
        amp::ampf<Precision> n0;
        amp::ampf<Precision> n1;
        amp::ampf<Precision> n2;
        amp::ampf<Precision> d0;
        amp::ampf<Precision> d1;
        amp::ampf<Precision> d2;
        amp::ampf<Precision> s0;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> xk;
        amp::ampf<Precision> xi;


        f = 0;
        df = 0;
        d2f = 0;
        
        //
        // special case: N=1
        //
        if( b.n==1 )
        {
            f = b.sy*b.y(0);
            df = 0;
            d2f = 0;
            return;
        }
        if( b.sy==0 )
        {
            f = 0;
            df = 0;
            d2f = 0;
            return;
        }
        ap::ap_error::make_assertion(b.sy>0);
        
        //
        // We assume than N>1 and B.SY>0. Find:
        // 1. pivot point (X[i] closest to T)
        // 2. width of interval containing X[i]
        //
        v = amp::abs<Precision>(b.x(0)-t);
        k = 0;
        for(i=1; i<=b.n-1; i++)
        {
            vv = b.x(i);
            if( amp::abs<Precision>(vv-t)<v )
            {
                v = amp::abs<Precision>(vv-t);
                k = i;
            }
        }
        
        //
        // pivot point found, calculate dNumerator and dDenominator
        //
        xk = b.x(k);
        v = t-xk;
        n0 = 0;
        n1 = 0;
        n2 = 0;
        d0 = 0;
        d1 = 0;
        d2 = 0;
        for(i=0; i<=b.n-1; i++)
        {
            if( i!=k )
            {
                xi = b.x(i);
                vv = amp::sqr<Precision>(t-xi);
                s0 = (t-xk)/(t-xi);
                s1 = (xk-xi)/vv;
                s2 = -2*(xk-xi)/(vv*(t-xi));
            }
            else
            {
                s0 = 1;
                s1 = 0;
                s2 = 0;
            }
            vv = b.w(i)*b.y(i);
            n0 = n0+s0*vv;
            n1 = n1+s1*vv;
            n2 = n2+s2*vv;
            vv = b.w(i);
            d0 = d0+s0*vv;
            d1 = d1+s1*vv;
            d2 = d2+s2*vv;
        }
        f = b.sy*n0/d0;
        df = b.sy*(n1*d0-n0*d1)/amp::sqr<Precision>(d0);
        d2f = b.sy*((n2*d0-n0*d2)*amp::sqr<Precision>(d0)-(n1*d0-n0*d1)*2*d0*d1)/amp::sqr<Precision>(amp::sqr<Precision>(d0));
    }


    /*************************************************************************
    This subroutine performs linear transformation of the argument.

    INPUT PARAMETERS:
        B       -   rational interpolant in barycentric form
        CA, CB  -   transformation coefficients: x = CA*t + CB

    OUTPUT PARAMETERS:
        B       -   transformed interpolant with X replaced by T

      -- ALGLIB PROJECT --
         Copyright 19.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentriclintransx(barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> ca,
        amp::ampf<Precision> cb)
    {
        int i;
        int j;
        amp::ampf<Precision> v;


        
        //
        // special case, replace by constant F(CB)
        //
        if( ca==0 )
        {
            b.sy = barycentriccalc<Precision>(b, cb);
            v = 1;
            for(i=0; i<=b.n-1; i++)
            {
                b.y(i) = 1;
                b.w(i) = v;
                v = -v;
            }
            return;
        }
        
        //
        // general case: CA<>0
        //
        for(i=0; i<=b.n-1; i++)
        {
            b.x(i) = (b.x(i)-cb)/ca;
        }
        if( ca<0 )
        {
            for(i=0; i<=b.n-1; i++)
            {
                if( i<b.n-1-i )
                {
                    j = b.n-1-i;
                    v = b.x(i);
                    b.x(i) = b.x(j);
                    b.x(j) = v;
                    v = b.y(i);
                    b.y(i) = b.y(j);
                    b.y(j) = v;
                    v = b.w(i);
                    b.w(i) = b.w(j);
                    b.w(j) = v;
                }
                else
                {
                    break;
                }
            }
        }
    }


    /*************************************************************************
    This  subroutine   performs   linear  transformation  of  the  barycentric
    interpolant.

    INPUT PARAMETERS:
        B       -   rational interpolant in barycentric form
        CA, CB  -   transformation coefficients: B2(x) = CA*B(x) + CB

    OUTPUT PARAMETERS:
        B       -   transformed interpolant

      -- ALGLIB PROJECT --
         Copyright 19.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentriclintransy(barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> ca,
        amp::ampf<Precision> cb)
    {
        int i;
        amp::ampf<Precision> v;


        for(i=0; i<=b.n-1; i++)
        {
            b.y(i) = ca*b.sy*b.y(i)+cb;
        }
        b.sy = 0;
        for(i=0; i<=b.n-1; i++)
        {
            b.sy = amp::maximum<Precision>(b.sy, amp::abs<Precision>(b.y(i)));
        }
        if( b.sy>0 )
        {
            v = 1/b.sy;
            amp::vmul(b.y.getvector(0, b.n-1), v);
        }
    }


    /*************************************************************************
    Extracts X/Y/W arrays from rational interpolant

    INPUT PARAMETERS:
        B   -   barycentric interpolant

    OUTPUT PARAMETERS:
        N   -   nodes count, N>0
        X   -   interpolation nodes, array[0..N-1]
        F   -   function values, array[0..N-1]
        W   -   barycentric weights, array[0..N-1]

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricunpack(const barycentricinterpolant<Precision>& b,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        amp::ampf<Precision> v;


        n = b.n;
        x.setlength(n);
        y.setlength(n);
        w.setlength(n);
        v = b.sy;
        amp::vmove(x.getvector(0, n-1), b.x.getvector(0, n-1));
        amp::vmove(y.getvector(0, n-1), b.y.getvector(0, n-1), v);
        amp::vmove(w.getvector(0, n-1), b.w.getvector(0, n-1));
    }


    /*************************************************************************
    Serialization of the barycentric interpolant

    INPUT PARAMETERS:
        B   -   barycentric interpolant

    OUTPUT PARAMETERS:
        RA      -   array of real numbers which contains interpolant,
                    array[0..RLen-1]
        RLen    -   RA lenght

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricserialize(const barycentricinterpolant<Precision>& b,
        ap::template_1d_array< amp::ampf<Precision> >& ra,
        int& ralen)
    {
        ralen = 2+2+3*b.n;
        ra.setlength(ralen);
        ra(0) = ralen;
        ra(1) = brcvnum;
        ra(2) = b.n;
        ra(3) = b.sy;
        amp::vmove(ra.getvector(4, 4+b.n-1), b.x.getvector(0, b.n-1));
        amp::vmove(ra.getvector(4+b.n, 4+2*b.n-1), b.y.getvector(0, b.n-1));
        amp::vmove(ra.getvector(4+2*b.n, 4+3*b.n-1), b.w.getvector(0, b.n-1));
    }


    /*************************************************************************
    Unserialization of the barycentric interpolant

    INPUT PARAMETERS:
        RA  -   array of real numbers which contains interpolant,

    OUTPUT PARAMETERS:
        B   -   barycentric interpolant

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricunserialize(const ap::template_1d_array< amp::ampf<Precision> >& ra,
        barycentricinterpolant<Precision>& b)
    {
        ap::ap_error::make_assertion(amp::round<Precision>(ra(1))==brcvnum);
        b.n = amp::round<Precision>(ra(2));
        b.sy = ra(3);
        b.x.setlength(b.n);
        b.y.setlength(b.n);
        b.w.setlength(b.n);
        amp::vmove(b.x.getvector(0, b.n-1), ra.getvector(4, 4+b.n-1));
        amp::vmove(b.y.getvector(0, b.n-1), ra.getvector(4+b.n, 4+2*b.n-1));
        amp::vmove(b.w.getvector(0, b.n-1), ra.getvector(4+2*b.n, 4+3*b.n-1));
    }


    /*************************************************************************
    Copying of the barycentric interpolant

    INPUT PARAMETERS:
        B   -   barycentric interpolant

    OUTPUT PARAMETERS:
        B2  -   copy(B1)

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentriccopy(const barycentricinterpolant<Precision>& b,
        barycentricinterpolant<Precision>& b2)
    {
        b2.n = b.n;
        b2.sy = b.sy;
        b2.x.setlength(b2.n);
        b2.y.setlength(b2.n);
        b2.w.setlength(b2.n);
        amp::vmove(b2.x.getvector(0, b2.n-1), b.x.getvector(0, b2.n-1));
        amp::vmove(b2.y.getvector(0, b2.n-1), b.y.getvector(0, b2.n-1));
        amp::vmove(b2.w.getvector(0, b2.n-1), b.w.getvector(0, b2.n-1));
    }


    /*************************************************************************
    Rational interpolant from X/Y/W arrays

    F(t) = SUM(i=0,n-1,w[i]*f[i]/(t-x[i])) / SUM(i=0,n-1,w[i]/(t-x[i]))

    INPUT PARAMETERS:
        X   -   interpolation nodes, array[0..N-1]
        F   -   function values, array[0..N-1]
        W   -   barycentric weights, array[0..N-1]
        N   -   nodes count, N>0

    OUTPUT PARAMETERS:
        B   -   barycentric interpolant built from (X, Y, W)

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricbuildxyw(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        barycentricinterpolant<Precision>& b)
    {
        ap::ap_error::make_assertion(n>0);
        
        //
        // fill X/Y/W
        //
        b.x.setlength(n);
        b.y.setlength(n);
        b.w.setlength(n);
        amp::vmove(b.x.getvector(0, n-1), x.getvector(0, n-1));
        amp::vmove(b.y.getvector(0, n-1), y.getvector(0, n-1));
        amp::vmove(b.w.getvector(0, n-1), w.getvector(0, n-1));
        b.n = n;
        
        //
        // Normalize
        //
        barycentricnormalize<Precision>(b);
    }


    /*************************************************************************
    Rational interpolant without poles

    The subroutine constructs the rational interpolating function without real
    poles  (see  'Barycentric rational interpolation with no  poles  and  high
    rates of approximation', Michael S. Floater. and  Kai  Hormann,  for  more
    information on this subject).

    Input parameters:
        X   -   interpolation nodes, array[0..N-1].
        Y   -   function values, array[0..N-1].
        N   -   number of nodes, N>0.
        D   -   order of the interpolation scheme, 0 <= D <= N-1.
                D<0 will cause an error.
                D>=N it will be replaced with D=N-1.
                if you don't know what D to choose, use small value about 3-5.

    Output parameters:
        B   -   barycentric interpolant.

    Note:
        this algorithm always succeeds and calculates the weights  with  close
        to machine precision.

      -- ALGLIB PROJECT --
         Copyright 17.06.2007 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricbuildfloaterhormann(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int d,
        barycentricinterpolant<Precision>& b)
    {
        amp::ampf<Precision> s0;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;
        int i;
        int j;
        int k;
        ap::template_1d_array< int > perm;
        ap::template_1d_array< amp::ampf<Precision> > wtemp;


        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(d>=0);
        
        //
        // Prepare
        //
        if( d>n-1 )
        {
            d = n-1;
        }
        b.n = n;
        
        //
        // special case: N=1
        //
        if( n==1 )
        {
            b.x.setlength(n);
            b.y.setlength(n);
            b.w.setlength(n);
            b.x(0) = x(0);
            b.y(0) = y(0);
            b.w(0) = 1;
            barycentricnormalize<Precision>(b);
            return;
        }
        
        //
        // Fill X/Y
        //
        b.x.setlength(n);
        b.y.setlength(n);
        amp::vmove(b.x.getvector(0, n-1), x.getvector(0, n-1));
        amp::vmove(b.y.getvector(0, n-1), y.getvector(0, n-1));
        tsort::tagsortfastr<Precision>(b.x, b.y, n);
        
        //
        // Calculate Wk
        //
        b.w.setlength(n);
        s0 = 1;
        for(k=1; k<=d; k++)
        {
            s0 = -s0;
        }
        for(k=0; k<=n-1; k++)
        {
            
            //
            // Wk
            //
            s = 0;
            for(i=ap::maxint(k-d, 0); i<=ap::minint(k, n-1-d); i++)
            {
                v = 1;
                for(j=i; j<=i+d; j++)
                {
                    if( j!=k )
                    {
                        v = v/amp::abs<Precision>(b.x(k)-b.x(j));
                    }
                }
                s = s+v;
            }
            b.w(k) = s0*s;
            
            //
            // Next S0
            //
            s0 = -s0;
        }
        
        //
        // Normalize
        //
        barycentricnormalize<Precision>(b);
    }


    /*************************************************************************
    Weghted rational least  squares  fitting  using  Floater-Hormann  rational
    functions  with  optimal  D  chosen  from  [0,9],  with  constraints   and
    individual weights.

    Equidistant  grid  with M node on [min(x),max(x)]  is  used to build basis
    functions. Different values of D are tried, optimal D (least WEIGHTED root
    mean square error) is chosen.  Task  is  linear,  so  linear least squares
    solver  is  used.  Complexity  of  this  computational  scheme is O(N*M^2)
    (mostly dominated by the least squares solver).

    SEE ALSO
    * BarycentricFitFloaterHormann(), "lightweight" fitting without invididual
      weights and constraints.

    INPUT PARAMETERS:
        X   -   points, array[0..N-1].
        Y   -   function values, array[0..N-1].
        W   -   weights, array[0..N-1]
                Each summand in square  sum  of  approximation deviations from
                given  values  is  multiplied  by  the square of corresponding
                weight. Fill it by 1's if you don't  want  to  solve  weighted
                task.
        N   -   number of points, N>0.
        XC  -   points where function values/derivatives are constrained,
                array[0..K-1].
        YC  -   values of constraints, array[0..K-1]
        DC  -   array[0..K-1], types of constraints:
                * DC[i]=0   means that S(XC[i])=YC[i]
                * DC[i]=1   means that S'(XC[i])=YC[i]
                SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
        K   -   number of constraints, 0<=K<M.
                K=0 means no constraints (XC/YC/DC are not used in such cases)
        M   -   number of basis functions ( = number_of_nodes), M>=2.

    OUTPUT PARAMETERS:
        Info-   same format as in LSFitLinearWC() subroutine.
                * Info>0    task is solved
                * Info<=0   an error occured:
                            -4 means inconvergence of internal SVD
                            -3 means inconsistent constraints
                            -1 means another errors in parameters passed
                               (N<=0, for example)
        B   -   barycentric interpolant.
        Rep -   report, same format as in LSFitLinearWC() subroutine.
                Following fields are set:
                * DBest         best value of the D parameter
                * RMSError      rms error on the (X,Y).
                * AvgError      average error on the (X,Y).
                * AvgRelError   average relative error on the non-zero Y
                * MaxError      maximum error
                                NON-WEIGHTED ERRORS ARE CALCULATED

    IMPORTANT:
        this subroitine doesn't calculate task's condition number for K<>0.

    SETTING CONSTRAINTS - DANGERS AND OPPORTUNITIES:

    Setting constraints can lead  to undesired  results,  like ill-conditioned
    behavior, or inconsistency being detected. From the other side,  it allows
    us to improve quality of the fit. Here we summarize  our  experience  with
    constrained barycentric interpolants:
    * excessive  constraints  can  be  inconsistent.   Floater-Hormann   basis
      functions aren't as flexible as splines (although they are very smooth).
    * the more evenly constraints are spread across [min(x),max(x)],  the more
      chances that they will be consistent
    * the  greater  is  M (given  fixed  constraints),  the  more chances that
      constraints will be consistent
    * in the general case, consistency of constraints IS NOT GUARANTEED.
    * in the several special cases, however, we CAN guarantee consistency.
    * one of this cases is constraints on the function  VALUES at the interval
      boundaries. Note that consustency of the  constraints  on  the  function
      DERIVATIVES is NOT guaranteed (you can use in such cases  cubic  splines
      which are more flexible).
    * another  special  case  is ONE constraint on the function value (OR, but
      not AND, derivative) anywhere in the interval

    Our final recommendation is to use constraints  WHEN  AND  ONLY  WHEN  you
    can't solve your task without them. Anything beyond  special  cases  given
    above is not guaranteed and may result in inconsistency.

      -- ALGLIB PROJECT --
         Copyright 18.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricfitfloaterhormannwc(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& xc,
        const ap::template_1d_array< amp::ampf<Precision> >& yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep)
    {
        int d;
        int i;
        amp::ampf<Precision> wrmscur;
        amp::ampf<Precision> wrmsbest;
        barycentricinterpolant<Precision> locb;
        barycentricfitreport<Precision> locrep;
        int locinfo;


        if( n<1 || m<2 || k<0 || k>=m )
        {
            info = -1;
            return;
        }
        
        //
        // Find optimal D
        //
        // Info is -3 by default (degenerate constraints).
        // If LocInfo will always be equal to -3, Info will remain equal to -3.
        // If at least once LocInfo will be -4, Info will be -4.
        //
        wrmsbest = amp::ampf<Precision>::getAlgoPascalMaxNumber();
        rep.dbest = -1;
        info = -3;
        for(d=0; d<=ap::minint(9, n-1); d++)
        {
            barycentricfitwcfixedd<Precision>(x, y, w, n, xc, yc, dc, k, m, d, locinfo, locb, locrep);
            ap::ap_error::make_assertion(locinfo==-4 || locinfo==-3 || locinfo>0);
            if( locinfo>0 )
            {
                
                //
                // Calculate weghted RMS
                //
                wrmscur = 0;
                for(i=0; i<=n-1; i++)
                {
                    wrmscur = wrmscur+amp::sqr<Precision>(w(i)*(y(i)-barycentriccalc<Precision>(locb, x(i))));
                }
                wrmscur = amp::sqrt<Precision>(wrmscur/n);
                if( wrmscur<wrmsbest || rep.dbest<0 )
                {
                    barycentriccopy<Precision>(locb, b);
                    rep.dbest = d;
                    info = 1;
                    rep.rmserror = locrep.rmserror;
                    rep.avgerror = locrep.avgerror;
                    rep.avgrelerror = locrep.avgrelerror;
                    rep.maxerror = locrep.maxerror;
                    rep.taskrcond = locrep.taskrcond;
                    wrmsbest = wrmscur;
                }
            }
            else
            {
                if( locinfo!=-3 && info<0 )
                {
                    info = locinfo;
                }
            }
        }
    }


    /*************************************************************************
    Rational least squares fitting, without weights and constraints.

    See BarycentricFitFloaterHormannWC() for more information.

      -- ALGLIB PROJECT --
         Copyright 18.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricfitfloaterhormann(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int m,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep)
    {
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > yc;
        ap::template_1d_array< int > dc;
        int i;


        if( n<1 )
        {
            info = -1;
            return;
        }
        w.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = 1;
        }
        barycentricfitfloaterhormannwc<Precision>(x, y, w, n, xc, yc, dc, 0, m, info, b, rep);
    }


    /*************************************************************************
    Normalization of barycentric interpolant:
    * B.N, B.X, B.Y and B.W are initialized
    * B.SY is NOT initialized
    * Y[] is normalized, scaling coefficient is stored in B.SY
    * W[] is normalized, no scaling coefficient is stored
    * X[] is sorted

    Internal subroutine.
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricnormalize(barycentricinterpolant<Precision>& b)
    {
        ap::template_1d_array< int > p1;
        ap::template_1d_array< int > p2;
        int i;
        int j;
        int j2;
        amp::ampf<Precision> v;


        
        //
        // Normalize task: |Y|<=1, |W|<=1, sort X[]
        //
        b.sy = 0;
        for(i=0; i<=b.n-1; i++)
        {
            b.sy = amp::maximum<Precision>(b.sy, amp::abs<Precision>(b.y(i)));
        }
        if( b.sy>0 && amp::abs<Precision>(b.sy-1)>10*amp::ampf<Precision>::getAlgoPascalEpsilon() )
        {
            v = 1/b.sy;
            amp::vmul(b.y.getvector(0, b.n-1), v);
        }
        v = 0;
        for(i=0; i<=b.n-1; i++)
        {
            v = amp::maximum<Precision>(v, amp::abs<Precision>(b.w(i)));
        }
        if( v>0 && amp::abs<Precision>(v-1)>10*amp::ampf<Precision>::getAlgoPascalEpsilon() )
        {
            v = 1/v;
            amp::vmul(b.w.getvector(0, b.n-1), v);
        }
        for(i=0; i<=b.n-2; i++)
        {
            if( b.x(i+1)<b.x(i) )
            {
                tsort::tagsort<Precision>(b.x, b.n, p1, p2);
                for(j=0; j<=b.n-1; j++)
                {
                    j2 = p2(j);
                    v = b.y(j);
                    b.y(j) = b.y(j2);
                    b.y(j2) = v;
                    v = b.w(j);
                    b.w(j) = b.w(j2);
                    b.w(j2) = v;
                }
                break;
            }
        }
    }


    /*************************************************************************
    Internal subroutine, calculates barycentric basis functions.
    Used for efficient simultaneous calculation of N basis functions.

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void barycentriccalcbasis(const barycentricinterpolant<Precision>& b,
        amp::ampf<Precision> t,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        amp::ampf<Precision> s2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;
        int i;
        int j;


        
        //
        // special case: N=1
        //
        if( b.n==1 )
        {
            y(0) = 1;
            return;
        }
        
        //
        // Here we assume that task is normalized, i.e.:
        // 1. abs(Y[i])<=1
        // 2. abs(W[i])<=1
        // 3. X[] is ordered
        //
        // First, we decide: should we use "safe" formula (guarded
        // against overflow) or fast one?
        //
        s = amp::abs<Precision>(t-b.x(0));
        for(i=0; i<=b.n-1; i++)
        {
            v = b.x(i);
            if( v==t )
            {
                for(j=0; j<=b.n-1; j++)
                {
                    y(j) = 0;
                }
                y(i) = 1;
                return;
            }
            v = amp::abs<Precision>(t-v);
            if( v<s )
            {
                s = v;
            }
        }
        s2 = 0;
        for(i=0; i<=b.n-1; i++)
        {
            v = s/(t-b.x(i));
            v = v*b.w(i);
            y(i) = v;
            s2 = s2+v;
        }
        v = 1/s2;
        amp::vmul(y.getvector(0, b.n-1), v);
    }


    /*************************************************************************
    Internal Floater-Hormann fitting subroutine for fixed D
    *************************************************************************/
    template<unsigned int Precision>
    void barycentricfitwcfixedd(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        ap::template_1d_array< amp::ampf<Precision> > xc,
        ap::template_1d_array< amp::ampf<Precision> > yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int d,
        int& info,
        barycentricinterpolant<Precision>& b,
        barycentricfitreport<Precision>& rep)
    {
        ap::template_2d_array< amp::ampf<Precision> > fmatrix;
        ap::template_2d_array< amp::ampf<Precision> > cmatrix;
        ap::template_1d_array< amp::ampf<Precision> > y2;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > sx;
        ap::template_1d_array< amp::ampf<Precision> > sy;
        ap::template_1d_array< amp::ampf<Precision> > sbf;
        ap::template_1d_array< amp::ampf<Precision> > xoriginal;
        ap::template_1d_array< amp::ampf<Precision> > yoriginal;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        lsfit::lsfitreport<Precision> lrep;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> mx;
        barycentricinterpolant<Precision> b2;
        int i;
        int j;
        int relcnt;
        amp::ampf<Precision> xa;
        amp::ampf<Precision> xb;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> sb;
        amp::ampf<Precision> decay;


        if( n<1 || m<2 || k<0 || k>=m )
        {
            info = -1;
            return;
        }
        for(i=0; i<=k-1; i++)
        {
            info = 0;
            if( dc(i)<0 )
            {
                info = -1;
            }
            if( dc(i)>1 )
            {
                info = -1;
            }
            if( info<0 )
            {
                return;
            }
        }
        
        //
        // weight decay for correct handling of task which becomes
        // degenerate after constraints are applied
        //
        decay = 10000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Scale X, Y, XC, YC
        //
        lsfit::lsfitscalexy<Precision>(x, y, n, xc, yc, dc, k, xa, xb, sa, sb, xoriginal, yoriginal);
        
        //
        // allocate space, initialize:
        // * FMatrix-   values of basis functions at X[]
        // * CMatrix-   values (derivatives) of basis functions at XC[]
        //
        y2.setlength(n+m);
        w2.setlength(n+m);
        fmatrix.setlength(n+m, m);
        if( k>0 )
        {
            cmatrix.setlength(k, m+1);
        }
        y2.setlength(n+m);
        w2.setlength(n+m);
        
        //
        // Prepare design and constraints matrices:
        // * fill constraints matrix
        // * fill first N rows of design matrix with values
        // * fill next M rows of design matrix with regularizing term
        // * append M zeros to Y
        // * append M elements, mean(abs(W)) each, to W
        //
        sx.setlength(m);
        sy.setlength(m);
        sbf.setlength(m);
        for(j=0; j<=m-1; j++)
        {
            sx(j) = amp::ampf<Precision>(2*j)/(amp::ampf<Precision>(m-1))-1;
        }
        for(i=0; i<=m-1; i++)
        {
            sy(i) = 1;
        }
        barycentricbuildfloaterhormann<Precision>(sx, sy, m, d, b2);
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            barycentriccalcbasis<Precision>(b2, x(i), sbf);
            amp::vmove(fmatrix.getrow(i, 0, m-1), sbf.getvector(0, m-1));
            y2(i) = y(i);
            w2(i) = w(i);
            mx = mx+amp::abs<Precision>(w(i))/n;
        }
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=m-1; j++)
            {
                if( i==j )
                {
                    fmatrix(n+i,j) = decay;
                }
                else
                {
                    fmatrix(n+i,j) = 0;
                }
            }
            y2(n+i) = 0;
            w2(n+i) = mx;
        }
        if( k>0 )
        {
            for(j=0; j<=m-1; j++)
            {
                for(i=0; i<=m-1; i++)
                {
                    sy(i) = 0;
                }
                sy(j) = 1;
                barycentricbuildfloaterhormann<Precision>(sx, sy, m, d, b2);
                for(i=0; i<=k-1; i++)
                {
                    ap::ap_error::make_assertion(dc(i)>=0 && dc(i)<=1);
                    barycentricdiff1<Precision>(b2, xc(i), v0, v1);
                    if( dc(i)==0 )
                    {
                        cmatrix(i,j) = v0;
                    }
                    if( dc(i)==1 )
                    {
                        cmatrix(i,j) = v1;
                    }
                }
            }
            for(i=0; i<=k-1; i++)
            {
                cmatrix(i,m) = yc(i);
            }
        }
        
        //
        // Solve constrained task
        //
        if( k>0 )
        {
            
            //
            // solve using regularization
            //
            lsfit::lsfitlinearwc<Precision>(y2, w2, fmatrix, cmatrix, n+m, m, k, info, tmp, lrep);
        }
        else
        {
            
            //
            // no constraints, no regularization needed
            //
            lsfit::lsfitlinearwc<Precision>(y, w, fmatrix, cmatrix, n, m, k, info, tmp, lrep);
        }
        if( info<0 )
        {
            return;
        }
        
        //
        // Generate interpolant and scale it
        //
        amp::vmove(sy.getvector(0, m-1), tmp.getvector(0, m-1));
        barycentricbuildfloaterhormann<Precision>(sx, sy, m, d, b);
        barycentriclintransx<Precision>(b, 2/(xb-xa), -(xa+xb)/(xb-xa));
        barycentriclintransy<Precision>(b, sb-sa, sa);
        
        //
        // Scale absolute errors obtained from LSFitLinearW.
        // Relative error should be calculated separately
        // (because of shifting/scaling of the task)
        //
        rep.taskrcond = lrep.taskrcond;
        rep.rmserror = lrep.rmserror*(sb-sa);
        rep.avgerror = lrep.avgerror*(sb-sa);
        rep.maxerror = lrep.maxerror*(sb-sa);
        rep.avgrelerror = 0;
        relcnt = 0;
        for(i=0; i<=n-1; i++)
        {
            if( yoriginal(i)!=0 )
            {
                rep.avgrelerror = rep.avgrelerror+amp::abs<Precision>(barycentriccalc<Precision>(b, xoriginal(i))-yoriginal(i))/amp::abs<Precision>(yoriginal(i));
                relcnt = relcnt+1;
            }
        }
        if( relcnt!=0 )
        {
            rep.avgrelerror = rep.avgrelerror/relcnt;
        }
    }
} // namespace

#endif
