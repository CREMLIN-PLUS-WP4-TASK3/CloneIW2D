/*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _polint_h
#define _polint_h

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
#include "ratint.h"
#include "apserv.h"
namespace polint
{
    /*************************************************************************
    Polynomial fitting report:
        TaskRCond       reciprocal of task's condition number
        RMSError        RMS error
        AvgError        average error
        AvgRelError     average relative error (for non-zero Y[I])
        MaxError        maximum error
    *************************************************************************/
    template<unsigned int Precision>
    class polynomialfitreport
    {
    public:
        amp::ampf<Precision> taskrcond;
        amp::ampf<Precision> rmserror;
        amp::ampf<Precision> avgerror;
        amp::ampf<Precision> avgrelerror;
        amp::ampf<Precision> maxerror;
    };




    template<unsigned int Precision>
    void polynomialbuild(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p);
    template<unsigned int Precision>
    void polynomialbuildeqdist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p);
    template<unsigned int Precision>
    void polynomialbuildcheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p);
    template<unsigned int Precision>
    void polynomialbuildcheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p);
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalceqdist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t);
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalccheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t);
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalccheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t);
    template<unsigned int Precision>
    void polynomialfit(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int m,
        int& info,
        ratint::barycentricinterpolant<Precision>& p,
        polynomialfitreport<Precision>& rep);
    template<unsigned int Precision>
    void polynomialfitwc(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        ap::template_1d_array< amp::ampf<Precision> > xc,
        ap::template_1d_array< amp::ampf<Precision> > yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int& info,
        ratint::barycentricinterpolant<Precision>& p,
        polynomialfitreport<Precision>& rep);


    /*************************************************************************
    Lagrange intepolant: generation of the model on the general grid.
    This function has O(N^2) complexity.

    INPUT PARAMETERS:
        X   -   abscissas, array[0..N-1]
        Y   -   function values, array[0..N-1]
        N   -   number of points, N>=1

    OIYTPUT PARAMETERS
        P   -   barycentric model which represents Lagrange interpolant
                (see ratint unit info and BarycentricCalc() description for
                more information).

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialbuild(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p)
    {
        int j;
        int k;
        ap::template_1d_array< amp::ampf<Precision> > w;
        amp::ampf<Precision> b;
        amp::ampf<Precision> a;
        amp::ampf<Precision> v;
        amp::ampf<Precision> mx;


        ap::ap_error::make_assertion(n>0);
        
        //
        // calculate W[j]
        // multi-pass algorithm is used to avoid overflow
        //
        w.setlength(n);
        a = x(0);
        b = x(0);
        for(j=0; j<=n-1; j++)
        {
            w(j) = 1;
            a = amp::minimum<Precision>(a, x(j));
            b = amp::maximum<Precision>(b, x(j));
        }
        for(k=0; k<=n-1; k++)
        {
            
            //
            // W[K] is used instead of 0.0 because
            // cycle on J does not touch K-th element
            // and we MUST get maximum from ALL elements
            //
            mx = amp::abs<Precision>(w(k));
            for(j=0; j<=n-1; j++)
            {
                if( j!=k )
                {
                    v = (b-a)/(x(j)-x(k));
                    w(j) = w(j)*v;
                    mx = amp::maximum<Precision>(mx, amp::abs<Precision>(w(j)));
                }
            }
            if( k%5==0 )
            {
                
                //
                // every 5-th run we renormalize W[]
                //
                v = 1/mx;
                amp::vmul(w.getvector(0, n-1), v);
            }
        }
        ratint::barycentricbuildxyw<Precision>(x, y, w, n, p);
    }


    /*************************************************************************
    Lagrange intepolant: generation of the model on equidistant grid.
    This function has O(N) complexity.

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        Y   -   function values at the nodes, array[0..N-1]
        N   -   number of points, N>=1
                for N=1 a constant model is constructed.

    OIYTPUT PARAMETERS
        P   -   barycentric model which represents Lagrange interpolant
                (see ratint unit info and BarycentricCalc() description for
                more information).

      -- ALGLIB --
         Copyright 03.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialbuildeqdist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>0);
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            x.setlength(1);
            w.setlength(1);
            x(0) = amp::ampf<Precision>("0.5")*(b+a);
            w(0) = 1;
            ratint::barycentricbuildxyw<Precision>(x, y, w, 1, p);
            return;
        }
        
        //
        // general case
        //
        x.setlength(n);
        w.setlength(n);
        v = 1;
        for(i=0; i<=n-1; i++)
        {
            w(i) = v;
            x(i) = a+(b-a)*i/(n-1);
            v = -v*(n-1-i);
            v = v/(i+1);
        }
        ratint::barycentricbuildxyw<Precision>(x, y, w, n, p);
    }


    /*************************************************************************
    Lagrange intepolant on Chebyshev grid (first kind).
    This function has O(N) complexity.

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        Y   -   function values at the nodes, array[0..N-1],
                Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n)))
        N   -   number of points, N>=1
                for N=1 a constant model is constructed.

    OIYTPUT PARAMETERS
        P   -   barycentric model which represents Lagrange interpolant
                (see ratint unit info and BarycentricCalc() description for
                more information).

      -- ALGLIB --
         Copyright 03.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialbuildcheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> v;
        amp::ampf<Precision> t;


        ap::ap_error::make_assertion(n>0);
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            x.setlength(1);
            w.setlength(1);
            x(0) = amp::ampf<Precision>("0.5")*(b+a);
            w(0) = 1;
            ratint::barycentricbuildxyw<Precision>(x, y, w, 1, p);
            return;
        }
        
        //
        // general case
        //
        x.setlength(n);
        w.setlength(n);
        v = 1;
        for(i=0; i<=n-1; i++)
        {
            t = amp::tan<Precision>(amp::ampf<Precision>("0.5")*amp::pi<Precision>()*(2*i+1)/(2*n));
            w(i) = 2*v*t/(1+amp::sqr<Precision>(t));
            x(i) = amp::ampf<Precision>("0.5")*(b+a)+amp::ampf<Precision>("0.5")*(b-a)*(1-amp::sqr<Precision>(t))/(1+amp::sqr<Precision>(t));
            v = -v;
        }
        ratint::barycentricbuildxyw<Precision>(x, y, w, n, p);
    }


    /*************************************************************************
    Lagrange intepolant on Chebyshev grid (second kind).
    This function has O(N) complexity.

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        Y   -   function values at the nodes, array[0..N-1],
                Y[I] = Y(0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1)))
        N   -   number of points, N>=1
                for N=1 a constant model is constructed.

    OIYTPUT PARAMETERS
        P   -   barycentric model which represents Lagrange interpolant
                (see ratint unit info and BarycentricCalc() description for
                more information).

      -- ALGLIB --
         Copyright 03.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialbuildcheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ratint::barycentricinterpolant<Precision>& p)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > x;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>0);
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            x.setlength(1);
            w.setlength(1);
            x(0) = amp::ampf<Precision>("0.5")*(b+a);
            w(0) = 1;
            ratint::barycentricbuildxyw<Precision>(x, y, w, 1, p);
            return;
        }
        
        //
        // general case
        //
        x.setlength(n);
        w.setlength(n);
        v = 1;
        for(i=0; i<=n-1; i++)
        {
            if( i==0 || i==n-1 )
            {
                w(i) = v*amp::ampf<Precision>("0.5");
            }
            else
            {
                w(i) = v;
            }
            x(i) = amp::ampf<Precision>("0.5")*(b+a)+amp::ampf<Precision>("0.5")*(b-a)*amp::cos<Precision>(amp::pi<Precision>()*i/(n-1));
            v = -v;
        }
        ratint::barycentricbuildxyw<Precision>(x, y, w, n, p);
    }


    /*************************************************************************
    Fast equidistant polynomial interpolation function with O(N) complexity

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        F   -   function values, array[0..N-1]
        N   -   number of points on equidistant grid, N>=1
                for N=1 a constant model is constructed.
        T   -   position where P(x) is calculated

    RESULT
        value of the Lagrange interpolant at T
        
    IMPORTANT
        this function provides fast interface which is not overflow-safe
        nor it is very precise.
        the best option is to use  PolynomialBuildEqDist()/BarycentricCalc()
        subroutines unless you are pretty sure that your data will not result
        in overflow.

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalceqdist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> s;
        amp::ampf<Precision> h;
        int i;
        int j;
        amp::ampf<Precision> w;
        amp::ampf<Precision> x;


        ap::ap_error::make_assertion(n>0);
        threshold = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber());
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            result = f(0);
            return result;
        }
        
        //
        // First, decide: should we use "safe" formula (guarded
        // against overflow) or fast one?
        //
        j = 0;
        s = t-a;
        for(i=1; i<=n-1; i++)
        {
            x = a+amp::ampf<Precision>(i)/(amp::ampf<Precision>(n-1))*(b-a);
            if( amp::abs<Precision>(t-x)<amp::abs<Precision>(s) )
            {
                s = t-x;
                j = i;
            }
        }
        if( s==0 )
        {
            result = f(j);
            return result;
        }
        if( amp::abs<Precision>(s)>threshold )
        {
            
            //
            // use fast formula
            //
            j = -1;
            s = amp::ampf<Precision>("1.0");
        }
        
        //
        // Calculate using safe or fast barycentric formula
        //
        s1 = 0;
        s2 = 0;
        w = amp::ampf<Precision>("1.0");
        h = (b-a)/(n-1);
        for(i=0; i<=n-1; i++)
        {
            if( i!=j )
            {
                v = s*w/(t-(a+i*h));
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            else
            {
                v = w;
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            w = -w*(n-1-i);
            w = w/(i+1);
        }
        result = s1/s2;
        return result;
    }


    /*************************************************************************
    Fast polynomial interpolation function on Chebyshev points (first kind)
    with O(N) complexity.

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        F   -   function values, array[0..N-1]
        N   -   number of points on Chebyshev grid (first kind),
                X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*(2*i+1)/(2*n))
                for N=1 a constant model is constructed.
        T   -   position where P(x) is calculated

    RESULT
        value of the Lagrange interpolant at T

    IMPORTANT
        this function provides fast interface which is not overflow-safe
        nor it is very precise.
        the best option is to use  PolIntBuildCheb1()/BarycentricCalc()
        subroutines unless you are pretty sure that your data will not result
        in overflow.

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalccheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> s;
        int i;
        int j;
        amp::ampf<Precision> a0;
        amp::ampf<Precision> delta;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> ca;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> tempc;
        amp::ampf<Precision> temps;
        amp::ampf<Precision> x;
        amp::ampf<Precision> w;
        amp::ampf<Precision> p1;


        ap::ap_error::make_assertion(n>0);
        threshold = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber());
        t = (t-amp::ampf<Precision>("0.5")*(a+b))/(amp::ampf<Precision>("0.5")*(b-a));
        
        //
        // Fast exit
        //
        if( n==1 )
        {
            result = f(0);
            return result;
        }
        
        //
        // Prepare information for the recurrence formula
        // used to calculate sin(pi*(2j+1)/(2n+2)) and
        // cos(pi*(2j+1)/(2n+2)):
        //
        // A0    = pi/(2n+2)
        // Delta = pi/(n+1)
        // Alpha = 2 sin^2 (Delta/2)
        // Beta  = sin(Delta)
        //
        // so that sin(..) = sin(A0+j*delta) and cos(..) = cos(A0+j*delta).
        // Then we use
        //
        // sin(x+delta) = sin(x) - (alpha*sin(x) - beta*cos(x))
        // cos(x+delta) = cos(x) - (alpha*cos(x) - beta*sin(x))
        //
        // to repeatedly calculate sin(..) and cos(..).
        //
        a0 = amp::pi<Precision>()/(2*(n-1)+2);
        delta = 2*amp::pi<Precision>()/(2*(n-1)+2);
        alpha = 2*amp::sqr<Precision>(amp::sin<Precision>(delta/2));
        beta = amp::sin<Precision>(delta);
        
        //
        // First, decide: should we use "safe" formula (guarded
        // against overflow) or fast one?
        //
        ca = amp::cos<Precision>(a0);
        sa = amp::sin<Precision>(a0);
        j = 0;
        x = ca;
        s = t-x;
        for(i=1; i<=n-1; i++)
        {
            
            //
            // Next X[i]
            //
            temps = sa-(alpha*sa-beta*ca);
            tempc = ca-(alpha*ca+beta*sa);
            sa = temps;
            ca = tempc;
            x = ca;
            
            //
            // Use X[i]
            //
            if( amp::abs<Precision>(t-x)<amp::abs<Precision>(s) )
            {
                s = t-x;
                j = i;
            }
        }
        if( s==0 )
        {
            result = f(j);
            return result;
        }
        if( amp::abs<Precision>(s)>threshold )
        {
            
            //
            // use fast formula
            //
            j = -1;
            s = amp::ampf<Precision>("1.0");
        }
        
        //
        // Calculate using safe or fast barycentric formula
        //
        s1 = 0;
        s2 = 0;
        ca = amp::cos<Precision>(a0);
        sa = amp::sin<Precision>(a0);
        p1 = amp::ampf<Precision>("1.0");
        for(i=0; i<=n-1; i++)
        {
            
            //
            // Calculate X[i], W[i]
            //
            x = ca;
            w = p1*sa;
            
            //
            // Proceed
            //
            if( i!=j )
            {
                v = s*w/(t-x);
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            else
            {
                v = w;
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            
            //
            // Next CA, SA, P1
            //
            temps = sa-(alpha*sa-beta*ca);
            tempc = ca-(alpha*ca+beta*sa);
            sa = temps;
            ca = tempc;
            p1 = -p1;
        }
        result = s1/s2;
        return result;
    }


    /*************************************************************************
    Fast polynomial interpolation function on Chebyshev points (second kind)
    with O(N) complexity.

    INPUT PARAMETERS:
        A   -   left boundary of [A,B]
        B   -   right boundary of [A,B]
        F   -   function values, array[0..N-1]
        N   -   number of points on Chebyshev grid (second kind),
                X[i] = 0.5*(B+A) + 0.5*(B-A)*Cos(PI*i/(n-1))
                for N=1 a constant model is constructed.
        T   -   position where P(x) is calculated

    RESULT
        value of the Lagrange interpolant at T

    IMPORTANT
        this function provides fast interface which is not overflow-safe
        nor it is very precise.
        the best option is to use PolIntBuildCheb2()/BarycentricCalc()
        subroutines unless you are pretty sure that your data will not result
        in overflow.

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> polynomialcalccheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        const ap::template_1d_array< amp::ampf<Precision> >& f,
        int n,
        amp::ampf<Precision> t)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> s1;
        amp::ampf<Precision> s2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> threshold;
        amp::ampf<Precision> s;
        int i;
        int j;
        amp::ampf<Precision> a0;
        amp::ampf<Precision> delta;
        amp::ampf<Precision> alpha;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> ca;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> tempc;
        amp::ampf<Precision> temps;
        amp::ampf<Precision> x;
        amp::ampf<Precision> w;
        amp::ampf<Precision> p1;


        ap::ap_error::make_assertion(n>0);
        threshold = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber());
        t = (t-amp::ampf<Precision>("0.5")*(a+b))/(amp::ampf<Precision>("0.5")*(b-a));
        
        //
        // Fast exit
        //
        if( n==1 )
        {
            result = f(0);
            return result;
        }
        
        //
        // Prepare information for the recurrence formula
        // used to calculate sin(pi*i/n) and
        // cos(pi*i/n):
        //
        // A0    = 0
        // Delta = pi/n
        // Alpha = 2 sin^2 (Delta/2)
        // Beta  = sin(Delta)
        //
        // so that sin(..) = sin(A0+j*delta) and cos(..) = cos(A0+j*delta).
        // Then we use
        //
        // sin(x+delta) = sin(x) - (alpha*sin(x) - beta*cos(x))
        // cos(x+delta) = cos(x) - (alpha*cos(x) - beta*sin(x))
        //
        // to repeatedly calculate sin(..) and cos(..).
        //
        a0 = amp::ampf<Precision>("0.0");
        delta = amp::pi<Precision>()/(n-1);
        alpha = 2*amp::sqr<Precision>(amp::sin<Precision>(delta/2));
        beta = amp::sin<Precision>(delta);
        
        //
        // First, decide: should we use "safe" formula (guarded
        // against overflow) or fast one?
        //
        ca = amp::cos<Precision>(a0);
        sa = amp::sin<Precision>(a0);
        j = 0;
        x = ca;
        s = t-x;
        for(i=1; i<=n-1; i++)
        {
            
            //
            // Next X[i]
            //
            temps = sa-(alpha*sa-beta*ca);
            tempc = ca-(alpha*ca+beta*sa);
            sa = temps;
            ca = tempc;
            x = ca;
            
            //
            // Use X[i]
            //
            if( amp::abs<Precision>(t-x)<amp::abs<Precision>(s) )
            {
                s = t-x;
                j = i;
            }
        }
        if( s==0 )
        {
            result = f(j);
            return result;
        }
        if( amp::abs<Precision>(s)>threshold )
        {
            
            //
            // use fast formula
            //
            j = -1;
            s = amp::ampf<Precision>("1.0");
        }
        
        //
        // Calculate using safe or fast barycentric formula
        //
        s1 = 0;
        s2 = 0;
        ca = amp::cos<Precision>(a0);
        sa = amp::sin<Precision>(a0);
        p1 = amp::ampf<Precision>("1.0");
        for(i=0; i<=n-1; i++)
        {
            
            //
            // Calculate X[i], W[i]
            //
            x = ca;
            if( i==0 || i==n-1 )
            {
                w = amp::ampf<Precision>("0.5")*p1;
            }
            else
            {
                w = amp::ampf<Precision>("1.0")*p1;
            }
            
            //
            // Proceed
            //
            if( i!=j )
            {
                v = s*w/(t-x);
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            else
            {
                v = w;
                s1 = s1+v*f(i);
                s2 = s2+v;
            }
            
            //
            // Next CA, SA, P1
            //
            temps = sa-(alpha*sa-beta*ca);
            tempc = ca-(alpha*ca+beta*sa);
            sa = temps;
            ca = tempc;
            p1 = -p1;
        }
        result = s1/s2;
        return result;
    }


    /*************************************************************************
    Least squares fitting by polynomial.

    This subroutine is "lightweight" alternative for more complex and feature-
    rich PolynomialFitWC().  See  PolynomialFitWC() for more information about
    subroutine parameters (we don't duplicate it here because of length)

      -- ALGLIB PROJECT --
         Copyright 12.10.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialfit(const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        int m,
        int& info,
        ratint::barycentricinterpolant<Precision>& p,
        polynomialfitreport<Precision>& rep)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > yc;
        ap::template_1d_array< int > dc;


        if( n>0 )
        {
            w.setlength(n);
            for(i=0; i<=n-1; i++)
            {
                w(i) = 1;
            }
        }
        polynomialfitwc<Precision>(x, y, w, n, xc, yc, dc, 0, m, info, p, rep);
    }


    /*************************************************************************
    Weighted  fitting  by  Chebyshev  polynomial  in  barycentric  form,  with
    constraints on function values or first derivatives.

    Small regularizing term is used when solving constrained tasks (to improve
    stability).

    Task is linear, so linear least squares solver is used. Complexity of this
    computational scheme is O(N*M^2), mostly dominated by least squares solver

    SEE ALSO:
        PolynomialFit()

    INPUT PARAMETERS:
        X   -   points, array[0..N-1].
        Y   -   function values, array[0..N-1].
        W   -   weights, array[0..N-1]
                Each summand in square  sum  of  approximation deviations from
                given  values  is  multiplied  by  the square of corresponding
                weight. Fill it by 1's if you don't  want  to  solve  weighted
                task.
        N   -   number of points, N>0.
        XC  -   points where polynomial values/derivatives are constrained,
                array[0..K-1].
        YC  -   values of constraints, array[0..K-1]
        DC  -   array[0..K-1], types of constraints:
                * DC[i]=0   means that P(XC[i])=YC[i]
                * DC[i]=1   means that P'(XC[i])=YC[i]
                SEE BELOW FOR IMPORTANT INFORMATION ON CONSTRAINTS
        K   -   number of constraints, 0<=K<M.
                K=0 means no constraints (XC/YC/DC are not used in such cases)
        M   -   number of basis functions (= polynomial_degree + 1), M>=1

    OUTPUT PARAMETERS:
        Info-   same format as in LSFitLinearW() subroutine:
                * Info>0    task is solved
                * Info<=0   an error occured:
                            -4 means inconvergence of internal SVD
                            -3 means inconsistent constraints
                            -1 means another errors in parameters passed
                               (N<=0, for example)
        P   -   interpolant in barycentric form.
        Rep -   report, same format as in LSFitLinearW() subroutine.
                Following fields are set:
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
    constrained regression splines:
    * even simple constraints can be inconsistent, see  Wikipedia  article  on
      this subject: http://en.wikipedia.org/wiki/Birkhoff_interpolation
    * the  greater  is  M (given  fixed  constraints),  the  more chances that
      constraints will be consistent
    * in the general case, consistency of constraints is NOT GUARANTEED.
    * in the one special cases, however, we can  guarantee  consistency.  This
      case  is:  M>1  and constraints on the function values (NOT DERIVATIVES)

    Our final recommendation is to use constraints  WHEN  AND  ONLY  when  you
    can't solve your task without them. Anything beyond  special  cases  given
    above is not guaranteed and may result in inconsistency.

      -- ALGLIB PROJECT --
         Copyright 10.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void polynomialfitwc(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        int n,
        ap::template_1d_array< amp::ampf<Precision> > xc,
        ap::template_1d_array< amp::ampf<Precision> > yc,
        const ap::template_1d_array< int >& dc,
        int k,
        int m,
        int& info,
        ratint::barycentricinterpolant<Precision>& p,
        polynomialfitreport<Precision>& rep)
    {
        amp::ampf<Precision> xa;
        amp::ampf<Precision> xb;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> sb;
        ap::template_1d_array< amp::ampf<Precision> > xoriginal;
        ap::template_1d_array< amp::ampf<Precision> > yoriginal;
        ap::template_1d_array< amp::ampf<Precision> > y2;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        ap::template_1d_array< amp::ampf<Precision> > tmp2;
        ap::template_1d_array< amp::ampf<Precision> > tmpdiff;
        ap::template_1d_array< amp::ampf<Precision> > bx;
        ap::template_1d_array< amp::ampf<Precision> > by;
        ap::template_1d_array< amp::ampf<Precision> > bw;
        ap::template_2d_array< amp::ampf<Precision> > fmatrix;
        ap::template_2d_array< amp::ampf<Precision> > cmatrix;
        int i;
        int j;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> decay;
        amp::ampf<Precision> u;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;
        int relcnt;
        lsfit::lsfitreport<Precision> lrep;


        if( m<1 || n<1 || k<0 || k>=m )
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
        // allocate space, initialize/fill:
        // * FMatrix-   values of basis functions at X[]
        // * CMatrix-   values (derivatives) of basis functions at XC[]
        // * fill constraints matrix
        // * fill first N rows of design matrix with values
        // * fill next M rows of design matrix with regularizing term
        // * append M zeros to Y
        // * append M elements, mean(abs(W)) each, to W
        //
        y2.setlength(n+m);
        w2.setlength(n+m);
        tmp.setlength(m);
        tmpdiff.setlength(m);
        fmatrix.setlength(n+m, m);
        if( k>0 )
        {
            cmatrix.setlength(k, m+1);
        }
        
        //
        // Fill design matrix, Y2, W2:
        // * first N rows with basis functions for original points
        // * next M rows with decay terms
        //
        for(i=0; i<=n-1; i++)
        {
            
            //
            // prepare Ith row
            // use Tmp for calculations to avoid multidimensional arrays overhead
            //
            for(j=0; j<=m-1; j++)
            {
                if( j==0 )
                {
                    tmp(j) = 1;
                }
                else
                {
                    if( j==1 )
                    {
                        tmp(j) = x(i);
                    }
                    else
                    {
                        tmp(j) = 2*x(i)*tmp(j-1)-tmp(j-2);
                    }
                }
            }
            amp::vmove(fmatrix.getrow(i, 0, m-1), tmp.getvector(0, m-1));
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
        }
        amp::vmove(y2.getvector(0, n-1), y.getvector(0, n-1));
        amp::vmove(w2.getvector(0, n-1), w.getvector(0, n-1));
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            mx = mx+amp::abs<Precision>(w(i));
        }
        mx = mx/n;
        for(i=0; i<=m-1; i++)
        {
            y2(n+i) = 0;
            w2(n+i) = mx;
        }
        
        //
        // fill constraints matrix
        //
        for(i=0; i<=k-1; i++)
        {
            
            //
            // prepare Ith row
            // use Tmp for basis function values,
            // TmpDiff for basos function derivatives
            //
            for(j=0; j<=m-1; j++)
            {
                if( j==0 )
                {
                    tmp(j) = 1;
                    tmpdiff(j) = 0;
                }
                else
                {
                    if( j==1 )
                    {
                        tmp(j) = xc(i);
                        tmpdiff(j) = 1;
                    }
                    else
                    {
                        tmp(j) = 2*xc(i)*tmp(j-1)-tmp(j-2);
                        tmpdiff(j) = 2*(tmp(j-1)+xc(i)*tmpdiff(j-1))-tmpdiff(j-2);
                    }
                }
            }
            if( dc(i)==0 )
            {
                amp::vmove(cmatrix.getrow(i, 0, m-1), tmp.getvector(0, m-1));
            }
            if( dc(i)==1 )
            {
                amp::vmove(cmatrix.getrow(i, 0, m-1), tmpdiff.getvector(0, m-1));
            }
            cmatrix(i,m) = yc(i);
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
            lsfit::lsfitlinearwc<Precision>(y, w, fmatrix, cmatrix, n, m, 0, info, tmp, lrep);
        }
        if( info<0 )
        {
            return;
        }
        
        //
        // Generate barycentric model and scale it
        // * BX, BY store barycentric model nodes
        // * FMatrix is reused (remember - it is at least MxM, what we need)
        //
        // Model intialization is done in O(M^2). In principle, it can be
        // done in O(M*log(M)), but before it we solved task with O(N*M^2)
        // complexity, so it is only a small amount of total time spent.
        //
        bx.setlength(m);
        by.setlength(m);
        bw.setlength(m);
        tmp2.setlength(m);
        s = 1;
        for(i=0; i<=m-1; i++)
        {
            if( m!=1 )
            {
                u = amp::cos<Precision>(amp::pi<Precision>()*i/(m-1));
            }
            else
            {
                u = 0;
            }
            v = 0;
            for(j=0; j<=m-1; j++)
            {
                if( j==0 )
                {
                    tmp2(j) = 1;
                }
                else
                {
                    if( j==1 )
                    {
                        tmp2(j) = u;
                    }
                    else
                    {
                        tmp2(j) = 2*u*tmp2(j-1)-tmp2(j-2);
                    }
                }
                v = v+tmp(j)*tmp2(j);
            }
            bx(i) = u;
            by(i) = v;
            bw(i) = s;
            if( i==0 || i==m-1 )
            {
                bw(i) = amp::ampf<Precision>("0.5")*bw(i);
            }
            s = -s;
        }
        ratint::barycentricbuildxyw<Precision>(bx, by, bw, m, p);
        ratint::barycentriclintransx<Precision>(p, 2/(xb-xa), -(xa+xb)/(xb-xa));
        ratint::barycentriclintransy<Precision>(p, sb-sa, sa);
        
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
                rep.avgrelerror = rep.avgrelerror+amp::abs<Precision>(ratint::barycentriccalc<Precision>(p, xoriginal(i))-yoriginal(i))/amp::abs<Precision>(yoriginal(i));
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
