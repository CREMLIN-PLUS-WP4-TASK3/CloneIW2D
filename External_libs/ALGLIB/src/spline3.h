/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _spline3_h
#define _spline3_h

#include "ap.h"
#include "amp.h"
namespace spline3
{
    template<unsigned int Precision>
    void buildlinearspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c);
    template<unsigned int Precision>
    void buildcubicspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        int boundltype,
        amp::ampf<Precision> boundl,
        int boundrtype,
        amp::ampf<Precision> boundr,
        ap::template_1d_array< amp::ampf<Precision> >& c);
    template<unsigned int Precision>
    void buildhermitespline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        ap::template_1d_array< amp::ampf<Precision> > d,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c);
    template<unsigned int Precision>
    void buildakimaspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c);
    template<unsigned int Precision>
    amp::ampf<Precision> splineinterpolation(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x);
    template<unsigned int Precision>
    void splinedifferentiation(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x,
        amp::ampf<Precision>& s,
        amp::ampf<Precision>& ds,
        amp::ampf<Precision>& d2s);
    template<unsigned int Precision>
    void splinecopy(const ap::template_1d_array< amp::ampf<Precision> >& c,
        ap::template_1d_array< amp::ampf<Precision> >& cc);
    template<unsigned int Precision>
    void splineunpack(const ap::template_1d_array< amp::ampf<Precision> >& c,
        int& n,
        ap::template_2d_array< amp::ampf<Precision> >& tbl);
    template<unsigned int Precision>
    void splinelintransx(ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    void splinelintransy(ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    amp::ampf<Precision> splineintegration(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x);
    template<unsigned int Precision>
    void spline3buildtable(int n,
        const int& diffn,
        ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const amp::ampf<Precision>& boundl,
        const amp::ampf<Precision>& boundr,
        ap::template_2d_array< amp::ampf<Precision> >& ctbl);
    template<unsigned int Precision>
    amp::ampf<Precision> spline3interpolate(int n,
        const ap::template_2d_array< amp::ampf<Precision> >& c,
        const amp::ampf<Precision>& x);
    template<unsigned int Precision>
    void heapsortpoints(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int n);
    template<unsigned int Precision>
    void heapsortdpoints(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        int n);
    template<unsigned int Precision>
    void solvetridiagonal(ap::template_1d_array< amp::ampf<Precision> > a,
        ap::template_1d_array< amp::ampf<Precision> > b,
        ap::template_1d_array< amp::ampf<Precision> > c,
        ap::template_1d_array< amp::ampf<Precision> > d,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    amp::ampf<Precision> diffthreepoint(amp::ampf<Precision> t,
        amp::ampf<Precision> x0,
        amp::ampf<Precision> f0,
        amp::ampf<Precision> x1,
        amp::ampf<Precision> f1,
        amp::ampf<Precision> x2,
        amp::ampf<Precision> f2);


    template<unsigned int Precision>
    void buildlinearspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        int i;
        int tblsize;


        ap::ap_error::make_assertion(n>=2);
        
        //
        // Sort points
        //
        heapsortpoints<Precision>(x, y, n);
        
        //
        // Fill C:
        //  C[0]            -   length(C)
        //  C[1]            -   type(C):
        //                      3 - general cubic spline
        //  C[2]            -   N
        //  C[3]...C[3+N-1] -   x[i], i = 0...N-1
        //  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
        //
        tblsize = 3+n+(n-1)*4;
        c.setbounds(0, tblsize-1);
        c(0) = tblsize;
        c(1) = 3;
        c(2) = n;
        for(i=0; i<=n-1; i++)
        {
            c(3+i) = x(i);
        }
        for(i=0; i<=n-2; i++)
        {
            c(3+n+4*i+0) = y(i);
            c(3+n+4*i+1) = (y(i+1)-y(i))/(x(i+1)-x(i));
            c(3+n+4*i+2) = 0;
            c(3+n+4*i+3) = 0;
        }
    }


    template<unsigned int Precision>
    void buildcubicspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        int boundltype,
        amp::ampf<Precision> boundl,
        int boundrtype,
        amp::ampf<Precision> boundr,
        ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        ap::template_1d_array< amp::ampf<Precision> > a1;
        ap::template_1d_array< amp::ampf<Precision> > a2;
        ap::template_1d_array< amp::ampf<Precision> > a3;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > d;
        int i;
        int tblsize;
        amp::ampf<Precision> delta;
        amp::ampf<Precision> delta2;
        amp::ampf<Precision> delta3;


        ap::ap_error::make_assertion(n>=2);
        ap::ap_error::make_assertion(boundltype==0 || boundltype==1 || boundltype==2);
        ap::ap_error::make_assertion(boundrtype==0 || boundrtype==1 || boundrtype==2);
        a1.setbounds(0, n-1);
        a2.setbounds(0, n-1);
        a3.setbounds(0, n-1);
        b.setbounds(0, n-1);
        
        //
        // Special case:
        // * N=2
        // * parabolic terminated boundary condition on both ends
        //
        if( n==2 && boundltype==0 && boundrtype==0 )
        {
            
            //
            // Change task type
            //
            boundltype = 2;
            boundl = 0;
            boundrtype = 2;
            boundr = 0;
        }
        
        //
        //
        // Sort points
        //
        heapsortpoints<Precision>(x, y, n);
        
        //
        // Left boundary conditions
        //
        if( boundltype==0 )
        {
            a1(0) = 0;
            a2(0) = 1;
            a3(0) = 1;
            b(0) = 2*(y(1)-y(0))/(x(1)-x(0));
        }
        if( boundltype==1 )
        {
            a1(0) = 0;
            a2(0) = 1;
            a3(0) = 0;
            b(0) = boundl;
        }
        if( boundltype==2 )
        {
            a1(0) = 0;
            a2(0) = 2;
            a3(0) = 1;
            b(0) = 3*(y(1)-y(0))/(x(1)-x(0))-amp::ampf<Precision>("0.5")*boundl*(x(1)-x(0));
        }
        
        //
        // Central conditions
        //
        for(i=1; i<=n-2; i++)
        {
            a1(i) = x(i+1)-x(i);
            a2(i) = 2*(x(i+1)-x(i-1));
            a3(i) = x(i)-x(i-1);
            b(i) = 3*(y(i)-y(i-1))/(x(i)-x(i-1))*(x(i+1)-x(i))+3*(y(i+1)-y(i))/(x(i+1)-x(i))*(x(i)-x(i-1));
        }
        
        //
        // Right boundary conditions
        //
        if( boundrtype==0 )
        {
            a1(n-1) = 1;
            a2(n-1) = 1;
            a3(n-1) = 0;
            b(n-1) = 2*(y(n-1)-y(n-2))/(x(n-1)-x(n-2));
        }
        if( boundrtype==1 )
        {
            a1(n-1) = 0;
            a2(n-1) = 1;
            a3(n-1) = 0;
            b(n-1) = boundr;
        }
        if( boundrtype==2 )
        {
            a1(n-1) = 1;
            a2(n-1) = 2;
            a3(n-1) = 0;
            b(n-1) = 3*(y(n-1)-y(n-2))/(x(n-1)-x(n-2))+amp::ampf<Precision>("0.5")*boundr*(x(n-1)-x(n-2));
        }
        
        //
        // Solve
        //
        solvetridiagonal<Precision>(a1, a2, a3, b, n, d);
        
        //
        // Now problem is reduced to the cubic Hermite spline
        //
        buildhermitespline<Precision>(x, y, d, n, c);
    }


    template<unsigned int Precision>
    void buildhermitespline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        ap::template_1d_array< amp::ampf<Precision> > d,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        int i;
        int tblsize;
        amp::ampf<Precision> delta;
        amp::ampf<Precision> delta2;
        amp::ampf<Precision> delta3;


        ap::ap_error::make_assertion(n>=2);
        
        //
        // Sort points
        //
        heapsortdpoints<Precision>(x, y, d, n);
        
        //
        // Fill C:
        //  C[0]            -   length(C)
        //  C[1]            -   type(C):
        //                      3 - general cubic spline
        //  C[2]            -   N
        //  C[3]...C[3+N-1] -   x[i], i = 0...N-1
        //  C[3+N]...C[3+N+(N-1)*4-1] - coefficients table
        //
        tblsize = 3+n+(n-1)*4;
        c.setbounds(0, tblsize-1);
        c(0) = tblsize;
        c(1) = 3;
        c(2) = n;
        for(i=0; i<=n-1; i++)
        {
            c(3+i) = x(i);
        }
        for(i=0; i<=n-2; i++)
        {
            delta = x(i+1)-x(i);
            delta2 = amp::sqr<Precision>(delta);
            delta3 = delta*delta2;
            c(3+n+4*i+0) = y(i);
            c(3+n+4*i+1) = d(i);
            c(3+n+4*i+2) = (3*(y(i+1)-y(i))-2*d(i)*delta-d(i+1)*delta)/delta2;
            c(3+n+4*i+3) = (2*(y(i)-y(i+1))+d(i)*delta+d(i+1)*delta)/delta3;
        }
    }


    template<unsigned int Precision>
    void buildakimaspline(ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        int i;
        ap::template_1d_array< amp::ampf<Precision> > d;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > diff;


        ap::ap_error::make_assertion(n>=5);
        
        //
        // Sort points
        //
        heapsortpoints<Precision>(x, y, n);
        
        //
        // Prepare W (weights), Diff (divided differences)
        //
        w.setbounds(1, n-2);
        diff.setbounds(0, n-2);
        for(i=0; i<=n-2; i++)
        {
            diff(i) = (y(i+1)-y(i))/(x(i+1)-x(i));
        }
        for(i=1; i<=n-2; i++)
        {
            w(i) = amp::abs<Precision>(diff(i)-diff(i-1));
        }
        
        //
        // Prepare Hermite interpolation scheme
        //
        d.setbounds(0, n-1);
        for(i=2; i<=n-3; i++)
        {
            if( amp::abs<Precision>(w(i-1))+amp::abs<Precision>(w(i+1))!=0 )
            {
                d(i) = (w(i+1)*diff(i-1)+w(i-1)*diff(i))/(w(i+1)+w(i-1));
            }
            else
            {
                d(i) = ((x(i+1)-x(i))*diff(i-1)+(x(i)-x(i-1))*diff(i))/(x(i+1)-x(i-1));
            }
        }
        d(0) = diffthreepoint<Precision>(x(0), x(0), y(0), x(1), y(1), x(2), y(2));
        d(1) = diffthreepoint<Precision>(x(1), x(0), y(0), x(1), y(1), x(2), y(2));
        d(n-2) = diffthreepoint<Precision>(x(n-2), x(n-3), y(n-3), x(n-2), y(n-2), x(n-1), y(n-1));
        d(n-1) = diffthreepoint<Precision>(x(n-1), x(n-3), y(n-3), x(n-2), y(n-2), x(n-1), y(n-1));
        
        //
        // Build Akima spline using Hermite interpolation scheme
        //
        buildhermitespline<Precision>(x, y, d, n, c);
    }


    template<unsigned int Precision>
    amp::ampf<Precision> splineinterpolation(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x)
    {
        amp::ampf<Precision> result;
        int n;
        int l;
        int r;
        int m;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        
        //
        // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
        //
        l = 3;
        r = 3+n-2+1;
        while( l!=r-1 )
        {
            m = (l+r)/2;
            if( c(m)>=x )
            {
                r = m;
            }
            else
            {
                l = m;
            }
        }
        
        //
        // Interpolation
        //
        x = x-c(l);
        m = 3+n+4*(l-3);
        result = c(m)+x*(c(m+1)+x*(c(m+2)+x*c(m+3)));
        return result;
    }


    template<unsigned int Precision>
    void splinedifferentiation(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x,
        amp::ampf<Precision>& s,
        amp::ampf<Precision>& ds,
        amp::ampf<Precision>& d2s)
    {
        int n;
        int l;
        int r;
        int m;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        
        //
        // Binary search
        //
        l = 3;
        r = 3+n-2+1;
        while( l!=r-1 )
        {
            m = (l+r)/2;
            if( c(m)>=x )
            {
                r = m;
            }
            else
            {
                l = m;
            }
        }
        
        //
        // Differentiation
        //
        x = x-c(l);
        m = 3+n+4*(l-3);
        s = c(m)+x*(c(m+1)+x*(c(m+2)+x*c(m+3)));
        ds = c(m+1)+2*x*c(m+2)+3*amp::sqr<Precision>(x)*c(m+3);
        d2s = 2*c(m+2)+6*x*c(m+3);
    }


    template<unsigned int Precision>
    void splinecopy(const ap::template_1d_array< amp::ampf<Precision> >& c,
        ap::template_1d_array< amp::ampf<Precision> >& cc)
    {
        int s;


        s = amp::round<Precision>(c(0));
        cc.setbounds(0, s-1);
        amp::vmove(cc.getvector(0, s-1), c.getvector(0, s-1));
    }


    template<unsigned int Precision>
    void splineunpack(const ap::template_1d_array< amp::ampf<Precision> >& c,
        int& n,
        ap::template_2d_array< amp::ampf<Precision> >& tbl)
    {
        int i;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        tbl.setbounds(0, n-2, 0, 5);
        
        //
        // Fill
        //
        for(i=0; i<=n-2; i++)
        {
            tbl(i,0) = c(3+i);
            tbl(i,1) = c(3+i+1);
            tbl(i,2) = c(3+n+4*i);
            tbl(i,3) = c(3+n+4*i+1);
            tbl(i,4) = c(3+n+4*i+2);
            tbl(i,5) = c(3+n+4*i+3);
        }
    }


    template<unsigned int Precision>
    void splinelintransx(ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        int i;
        int n;
        amp::ampf<Precision> v;
        amp::ampf<Precision> dv;
        amp::ampf<Precision> d2v;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > d;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        
        //
        // Special case: A=0
        //
        if( a==0 )
        {
            v = splineinterpolation<Precision>(c, b);
            for(i=0; i<=n-2; i++)
            {
                c(3+n+4*i) = v;
                c(3+n+4*i+1) = 0;
                c(3+n+4*i+2) = 0;
                c(3+n+4*i+3) = 0;
            }
            return;
        }
        
        //
        // General case: A<>0.
        // Unpack, X, Y, dY/dX.
        // Scale and pack again.
        //
        x.setbounds(0, n-1);
        y.setbounds(0, n-1);
        d.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            x(i) = c(3+i);
            splinedifferentiation<Precision>(c, x(i), v, dv, d2v);
            x(i) = (x(i)-b)/a;
            y(i) = v;
            d(i) = a*dv;
        }
        buildhermitespline<Precision>(x, y, d, n, c);
    }


    template<unsigned int Precision>
    void splinelintransy(ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        int i;
        int n;
        amp::ampf<Precision> v;
        amp::ampf<Precision> dv;
        amp::ampf<Precision> d2v;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > d;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        
        //
        // Special case: A=0
        //
        for(i=0; i<=n-2; i++)
        {
            c(3+n+4*i) = a*c(3+n+4*i)+b;
            c(3+n+4*i+1) = a*c(3+n+4*i+1);
            c(3+n+4*i+2) = a*c(3+n+4*i+2);
            c(3+n+4*i+3) = a*c(3+n+4*i+3);
        }
    }


    template<unsigned int Precision>
    amp::ampf<Precision> splineintegration(const ap::template_1d_array< amp::ampf<Precision> >& c,
        amp::ampf<Precision> x)
    {
        amp::ampf<Precision> result;
        int n;
        int i;
        int l;
        int r;
        int m;
        amp::ampf<Precision> w;


        ap::ap_error::make_assertion(amp::round<Precision>(c(1))==3);
        n = amp::round<Precision>(c(2));
        
        //
        // Binary search in the [ x[0], ..., x[n-2] ] (x[n-1] is not included)
        //
        l = 3;
        r = 3+n-2+1;
        while( l!=r-1 )
        {
            m = (l+r)/2;
            if( c(m)>=x )
            {
                r = m;
            }
            else
            {
                l = m;
            }
        }
        
        //
        // Integration
        //
        result = 0;
        for(i=3; i<=l-1; i++)
        {
            w = c(i+1)-c(i);
            m = 3+n+4*(i-3);
            result = result+c(m)*w;
            result = result+c(m+1)*amp::sqr<Precision>(w)/2;
            result = result+c(m+2)*amp::sqr<Precision>(w)*w/3;
            result = result+c(m+3)*amp::sqr<Precision>(amp::sqr<Precision>(w))/4;
        }
        w = x-c(l);
        m = 3+n+4*(l-3);
        result = result+c(m)*w;
        result = result+c(m+1)*amp::sqr<Precision>(w)/2;
        result = result+c(m+2)*amp::sqr<Precision>(w)*w/3;
        result = result+c(m+3)*amp::sqr<Precision>(amp::sqr<Precision>(w))/4;
        return result;
    }


    template<unsigned int Precision>
    void spline3buildtable(int n,
        const int& diffn,
        ap::template_1d_array< amp::ampf<Precision> > x,
        ap::template_1d_array< amp::ampf<Precision> > y,
        const amp::ampf<Precision>& boundl,
        const amp::ampf<Precision>& boundr,
        ap::template_2d_array< amp::ampf<Precision> >& ctbl)
    {
        bool c;
        int e;
        int g;
        amp::ampf<Precision> tmp;
        int nxm1;
        int i;
        int j;
        amp::ampf<Precision> dx;
        amp::ampf<Precision> dxj;
        amp::ampf<Precision> dyj;
        amp::ampf<Precision> dxjp1;
        amp::ampf<Precision> dyjp1;
        amp::ampf<Precision> dxp;
        amp::ampf<Precision> dyp;
        amp::ampf<Precision> yppa;
        amp::ampf<Precision> yppb;
        amp::ampf<Precision> pj;
        amp::ampf<Precision> b1;
        amp::ampf<Precision> b2;
        amp::ampf<Precision> b3;
        amp::ampf<Precision> b4;


        n = n-1;
        g = (n+1)/2;
        do
        {
            i = g;
            do
            {
                j = i-g;
                c = true;
                do
                {
                    if( x(j)<=x(j+g) )
                    {
                        c = false;
                    }
                    else
                    {
                        tmp = x(j);
                        x(j) = x(j+g);
                        x(j+g) = tmp;
                        tmp = y(j);
                        y(j) = y(j+g);
                        y(j+g) = tmp;
                    }
                    j = j-1;
                }
                while( j>=0 && c );
                i = i+1;
            }
            while( i<=n );
            g = g/2;
        }
        while( g>0 );
        ctbl.setbounds(0, 4, 0, n);
        n = n+1;
        if( diffn==1 )
        {
            b1 = 1;
            b2 = 6/(x(1)-x(0))*((y(1)-y(0))/(x(1)-x(0))-boundl);
            b3 = 1;
            b4 = 6/(x(n-1)-x(n-2))*(boundr-(y(n-1)-y(n-2))/(x(n-1)-x(n-2)));
        }
        else
        {
            b1 = 0;
            b2 = 2*boundl;
            b3 = 0;
            b4 = 2*boundr;
        }
        nxm1 = n-1;
        if( n>=2 )
        {
            if( n>2 )
            {
                dxj = x(1)-x(0);
                dyj = y(1)-y(0);
                j = 2;
                while( j<=nxm1 )
                {
                    dxjp1 = x(j)-x(j-1);
                    dyjp1 = y(j)-y(j-1);
                    dxp = dxj+dxjp1;
                    ctbl(1,j-1) = dxjp1/dxp;
                    ctbl(2,j-1) = 1-ctbl(1,j-1);
                    ctbl(3,j-1) = 6*(dyjp1/dxjp1-dyj/dxj)/dxp;
                    dxj = dxjp1;
                    dyj = dyjp1;
                    j = j+1;
                }
            }
            ctbl(1,0) = -b1/2;
            ctbl(2,0) = b2/2;
            if( n!=2 )
            {
                j = 2;
                while( j<=nxm1 )
                {
                    pj = ctbl(2,j-1)*ctbl(1,j-2)+2;
                    ctbl(1,j-1) = -ctbl(1,j-1)/pj;
                    ctbl(2,j-1) = (ctbl(3,j-1)-ctbl(2,j-1)*ctbl(2,j-2))/pj;
                    j = j+1;
                }
            }
            yppb = (b4-b3*ctbl(2,nxm1-1))/(b3*ctbl(1,nxm1-1)+2);
            i = 1;
            while( i<=nxm1 )
            {
                j = n-i;
                yppa = ctbl(1,j-1)*yppb+ctbl(2,j-1);
                dx = x(j)-x(j-1);
                ctbl(3,j-1) = (yppb-yppa)/dx/6;
                ctbl(2,j-1) = yppa/2;
                ctbl(1,j-1) = (y(j)-y(j-1))/dx-(ctbl(2,j-1)+ctbl(3,j-1)*dx)*dx;
                yppb = yppa;
                i = i+1;
            }
            for(i=1; i<=n; i++)
            {
                ctbl(0,i-1) = y(i-1);
                ctbl(4,i-1) = x(i-1);
            }
        }
    }


    template<unsigned int Precision>
    amp::ampf<Precision> spline3interpolate(int n,
        const ap::template_2d_array< amp::ampf<Precision> >& c,
        const amp::ampf<Precision>& x)
    {
        amp::ampf<Precision> result;
        int i;
        int l;
        int half;
        int first;
        int middle;


        n = n-1;
        l = n;
        first = 0;
        while( l>0 )
        {
            half = l/2;
            middle = first+half;
            if( c(4,middle)<x )
            {
                first = middle+1;
                l = l-half-1;
            }
            else
            {
                l = half;
            }
        }
        i = first-1;
        if( i<0 )
        {
            i = 0;
        }
        result = c(0,i)+(x-c(4,i))*(c(1,i)+(x-c(4,i))*(c(2,i)+c(3,i)*(x-c(4,i))));
        return result;
    }


    template<unsigned int Precision>
    void heapsortpoints(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int n)
    {
        int i;
        int j;
        int k;
        int t;
        amp::ampf<Precision> tmp;
        bool isascending;
        bool isdescending;


        
        //
        // Test for already sorted set
        //
        isascending = true;
        isdescending = true;
        for(i=1; i<=n-1; i++)
        {
            isascending = isascending && x(i)>x(i-1);
            isdescending = isdescending && x(i)<x(i-1);
        }
        if( isascending )
        {
            return;
        }
        if( isdescending )
        {
            for(i=0; i<=n-1; i++)
            {
                j = n-1-i;
                if( j<=i )
                {
                    break;
                }
                tmp = x(i);
                x(i) = x(j);
                x(j) = tmp;
                tmp = y(i);
                y(i) = y(j);
                y(j) = tmp;
            }
            return;
        }
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            return;
        }
        
        //
        // General case
        //
        i = 2;
        do
        {
            t = i;
            while( t!=1 )
            {
                k = t/2;
                if( x(k-1)>=x(t-1) )
                {
                    t = 1;
                }
                else
                {
                    tmp = x(k-1);
                    x(k-1) = x(t-1);
                    x(t-1) = tmp;
                    tmp = y(k-1);
                    y(k-1) = y(t-1);
                    y(t-1) = tmp;
                    t = k;
                }
            }
            i = i+1;
        }
        while( i<=n );
        i = n-1;
        do
        {
            tmp = x(i);
            x(i) = x(0);
            x(0) = tmp;
            tmp = y(i);
            y(i) = y(0);
            y(0) = tmp;
            t = 1;
            while( t!=0 )
            {
                k = 2*t;
                if( k>i )
                {
                    t = 0;
                }
                else
                {
                    if( k<i )
                    {
                        if( x(k)>x(k-1) )
                        {
                            k = k+1;
                        }
                    }
                    if( x(t-1)>=x(k-1) )
                    {
                        t = 0;
                    }
                    else
                    {
                        tmp = x(k-1);
                        x(k-1) = x(t-1);
                        x(t-1) = tmp;
                        tmp = y(k-1);
                        y(k-1) = y(t-1);
                        y(t-1) = tmp;
                        t = k;
                    }
                }
            }
            i = i-1;
        }
        while( i>=1 );
    }


    template<unsigned int Precision>
    void heapsortdpoints(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        int n)
    {
        int i;
        int j;
        int k;
        int t;
        amp::ampf<Precision> tmp;
        bool isascending;
        bool isdescending;


        
        //
        // Test for already sorted set
        //
        isascending = true;
        isdescending = true;
        for(i=1; i<=n-1; i++)
        {
            isascending = isascending && x(i)>x(i-1);
            isdescending = isdescending && x(i)<x(i-1);
        }
        if( isascending )
        {
            return;
        }
        if( isdescending )
        {
            for(i=0; i<=n-1; i++)
            {
                j = n-1-i;
                if( j<=i )
                {
                    break;
                }
                tmp = x(i);
                x(i) = x(j);
                x(j) = tmp;
                tmp = y(i);
                y(i) = y(j);
                y(j) = tmp;
                tmp = d(i);
                d(i) = d(j);
                d(j) = tmp;
            }
            return;
        }
        
        //
        // Special case: N=1
        //
        if( n==1 )
        {
            return;
        }
        
        //
        // General case
        //
        i = 2;
        do
        {
            t = i;
            while( t!=1 )
            {
                k = t/2;
                if( x(k-1)>=x(t-1) )
                {
                    t = 1;
                }
                else
                {
                    tmp = x(k-1);
                    x(k-1) = x(t-1);
                    x(t-1) = tmp;
                    tmp = y(k-1);
                    y(k-1) = y(t-1);
                    y(t-1) = tmp;
                    tmp = d(k-1);
                    d(k-1) = d(t-1);
                    d(t-1) = tmp;
                    t = k;
                }
            }
            i = i+1;
        }
        while( i<=n );
        i = n-1;
        do
        {
            tmp = x(i);
            x(i) = x(0);
            x(0) = tmp;
            tmp = y(i);
            y(i) = y(0);
            y(0) = tmp;
            tmp = d(i);
            d(i) = d(0);
            d(0) = tmp;
            t = 1;
            while( t!=0 )
            {
                k = 2*t;
                if( k>i )
                {
                    t = 0;
                }
                else
                {
                    if( k<i )
                    {
                        if( x(k)>x(k-1) )
                        {
                            k = k+1;
                        }
                    }
                    if( x(t-1)>=x(k-1) )
                    {
                        t = 0;
                    }
                    else
                    {
                        tmp = x(k-1);
                        x(k-1) = x(t-1);
                        x(t-1) = tmp;
                        tmp = y(k-1);
                        y(k-1) = y(t-1);
                        y(t-1) = tmp;
                        tmp = d(k-1);
                        d(k-1) = d(t-1);
                        d(t-1) = tmp;
                        t = k;
                    }
                }
            }
            i = i-1;
        }
        while( i>=1 );
    }


    template<unsigned int Precision>
    void solvetridiagonal(ap::template_1d_array< amp::ampf<Precision> > a,
        ap::template_1d_array< amp::ampf<Precision> > b,
        ap::template_1d_array< amp::ampf<Precision> > c,
        ap::template_1d_array< amp::ampf<Precision> > d,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        int k;
        amp::ampf<Precision> t;


        x.setbounds(0, n-1);
        a(0) = 0;
        c(n-1) = 0;
        for(k=1; k<=n-1; k++)
        {
            t = a(k)/b(k-1);
            b(k) = b(k)-t*c(k-1);
            d(k) = d(k)-t*d(k-1);
        }
        x(n-1) = d(n-1)/b(n-1);
        for(k=n-2; k>=0; k--)
        {
            x(k) = (d(k)-c(k)*x(k+1))/b(k);
        }
    }


    template<unsigned int Precision>
    amp::ampf<Precision> diffthreepoint(amp::ampf<Precision> t,
        amp::ampf<Precision> x0,
        amp::ampf<Precision> f0,
        amp::ampf<Precision> x1,
        amp::ampf<Precision> f1,
        amp::ampf<Precision> x2,
        amp::ampf<Precision> f2)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;


        t = t-x0;
        x1 = x1-x0;
        x2 = x2-x0;
        a = (f2-f0-x2/x1*(f1-f0))/(amp::sqr<Precision>(x2)-x1*x2);
        b = (f1-f0-a*amp::sqr<Precision>(x1))/x1;
        result = 2*a*t+b;
        return result;
    }
} // namespace

#endif
