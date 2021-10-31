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

#ifndef _apserv_h
#define _apserv_h

#include "ap.h"
#include "amp.h"
namespace apserv
{
    template<unsigned int Precision>
    void taskgenint1d(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    void taskgenint1dequidist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    void taskgenint1dcheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    void taskgenint1dcheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    bool apservaredistinct(ap::template_1d_array< amp::ampf<Precision> > x,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> safepythag2(amp::ampf<Precision> x,
        amp::ampf<Precision> y);
    template<unsigned int Precision>
    amp::ampf<Precision> safepythag3(amp::ampf<Precision> x,
        amp::ampf<Precision> y,
        amp::ampf<Precision> z);
    template<unsigned int Precision>
    void apperiodicmap(amp::ampf<Precision>& x,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision>& k);


    /*************************************************************************
    This  function  generates  1-dimensional  general  interpolation task with
    moderate Lipshitz constant (close to 1.0)

    If N=1 then suborutine generates only one point at the middle of [A,B]

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void taskgenint1d(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        int i;
        amp::ampf<Precision> h;


        ap::ap_error::make_assertion(n>=1);
        x.setlength(n);
        y.setlength(n);
        if( n>1 )
        {
            x(0) = a;
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
            h = (b-a)/(n-1);
            for(i=1; i<=n-1; i++)
            {
                if( i!=n-1 )
                {
                    x(i) = a+(i+amp::ampf<Precision>("0.2")*(2*amp::ampf<Precision>::getRandom()-1))*h;
                }
                else
                {
                    x(i) = b;
                }
                y(i) = y(i-1)+(2*amp::ampf<Precision>::getRandom()-1)*(x(i)-x(i-1));
            }
        }
        else
        {
            x(0) = amp::ampf<Precision>("0.5")*(a+b);
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
        }
    }


    /*************************************************************************
    This function generates  1-dimensional equidistant interpolation task with
    moderate Lipshitz constant (close to 1.0)

    If N=1 then suborutine generates only one point at the middle of [A,B]

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void taskgenint1dequidist(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        int i;
        amp::ampf<Precision> h;


        ap::ap_error::make_assertion(n>=1);
        x.setlength(n);
        y.setlength(n);
        if( n>1 )
        {
            x(0) = a;
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
            h = (b-a)/(n-1);
            for(i=1; i<=n-1; i++)
            {
                x(i) = a+i*h;
                y(i) = y(i-1)+(2*amp::ampf<Precision>::getRandom()-1)*h;
            }
        }
        else
        {
            x(0) = amp::ampf<Precision>("0.5")*(a+b);
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
        }
    }


    /*************************************************************************
    This function generates  1-dimensional Chebyshev-1 interpolation task with
    moderate Lipshitz constant (close to 1.0)

    If N=1 then suborutine generates only one point at the middle of [A,B]

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void taskgenint1dcheb1(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        int i;


        ap::ap_error::make_assertion(n>=1);
        x.setlength(n);
        y.setlength(n);
        if( n>1 )
        {
            for(i=0; i<=n-1; i++)
            {
                x(i) = amp::ampf<Precision>("0.5")*(b+a)+amp::ampf<Precision>("0.5")*(b-a)*amp::cos<Precision>(amp::pi<Precision>()*(2*i+1)/(2*n));
                if( i==0 )
                {
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                else
                {
                    y(i) = y(i-1)+(2*amp::ampf<Precision>::getRandom()-1)*(x(i)-x(i-1));
                }
            }
        }
        else
        {
            x(0) = amp::ampf<Precision>("0.5")*(a+b);
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
        }
    }


    /*************************************************************************
    This function generates  1-dimensional Chebyshev-2 interpolation task with
    moderate Lipshitz constant (close to 1.0)

    If N=1 then suborutine generates only one point at the middle of [A,B]

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void taskgenint1dcheb2(amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        int i;


        ap::ap_error::make_assertion(n>=1);
        x.setlength(n);
        y.setlength(n);
        if( n>1 )
        {
            for(i=0; i<=n-1; i++)
            {
                x(i) = amp::ampf<Precision>("0.5")*(b+a)+amp::ampf<Precision>("0.5")*(b-a)*amp::cos<Precision>(amp::pi<Precision>()*i/(n-1));
                if( i==0 )
                {
                    y(i) = 2*amp::ampf<Precision>::getRandom()-1;
                }
                else
                {
                    y(i) = y(i-1)+(2*amp::ampf<Precision>::getRandom()-1)*(x(i)-x(i-1));
                }
            }
        }
        else
        {
            x(0) = amp::ampf<Precision>("0.5")*(a+b);
            y(0) = 2*amp::ampf<Precision>::getRandom()-1;
        }
    }


    /*************************************************************************
    This function checks that all values from X[] are distinct. It does more
    than just usual floating point comparison:
    * first, it calculates max(X) and min(X)
    * second, it maps X[] from [min,max] to [1,2]
    * only at this stage actual comparison is done

    The meaning of such check is to ensure that all values are "distinct enough"
    and will not cause interpolation subroutine to fail.

    NOTE:
        X[] must be sorted by ascending (subroutine ASSERT's it)

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool apservaredistinct(ap::template_1d_array< amp::ampf<Precision> > x,
        int n)
    {
        bool result;
        bool issorted;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        int i;


        ap::ap_error::make_assertion(n>=1);
        if( n==1 )
        {
            
            //
            // everything is alright, it is up to caller to decide whether it
            // can interpolate something with just one point
            //
            result = true;
            return result;
        }
        a = x(0);
        b = x(0);
        for(i=1; i<=n-1; i++)
        {
            a = amp::minimum<Precision>(a, x(i));
            b = amp::maximum<Precision>(b, x(i));
            ap::ap_error::make_assertion(x(i)>=x(i-1));
        }
        for(i=0; i<=n-1; i++)
        {
            x(i) = (x(i)-a)/(b-a)+1;
        }
        for(i=1; i<=n-1; i++)
        {
            if( x(i)==x(i-1) )
            {
                result = false;
                return result;
            }
        }
        result = true;
        return result;
    }


    /*************************************************************************
    Safe sqrt(x^2+y^2)

      -- ALGLIB --
         Copyright by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> safepythag2(amp::ampf<Precision> x,
        amp::ampf<Precision> y)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> w;
        amp::ampf<Precision> xabs;
        amp::ampf<Precision> yabs;
        amp::ampf<Precision> z;


        xabs = amp::abs<Precision>(x);
        yabs = amp::abs<Precision>(y);
        w = amp::maximum<Precision>(xabs, yabs);
        z = amp::minimum<Precision>(xabs, yabs);
        if( z==0 )
        {
            result = w;
        }
        else
        {
            result = w*amp::sqrt<Precision>(1+amp::sqr<Precision>(z/w));
        }
        return result;
    }


    /*************************************************************************
    Safe sqrt(x^2+y^2)

      -- ALGLIB --
         Copyright by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> safepythag3(amp::ampf<Precision> x,
        amp::ampf<Precision> y,
        amp::ampf<Precision> z)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> w;


        w = amp::maximum<Precision>(amp::abs<Precision>(x), amp::maximum<Precision>(amp::abs<Precision>(y), amp::abs<Precision>(z)));
        if( w==0 )
        {
            result = 0;
            return result;
        }
        x = x/w;
        y = y/w;
        z = z/w;
        result = w*amp::sqrt<Precision>(amp::sqr<Precision>(x)+amp::sqr<Precision>(y)+amp::sqr<Precision>(z));
        return result;
    }


    /*************************************************************************
    This function makes periodic mapping of X to [A,B].

    It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
    such that X = T + K*(B-A).

    NOTES:
    * K is represented as real value, although actually it is integer
    * T is guaranteed to be in [A,B]
    * T replaces X

      -- ALGLIB --
         Copyright by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void apperiodicmap(amp::ampf<Precision>& x,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        amp::ampf<Precision>& k)
    {
        ap::ap_error::make_assertion(a<b);
        k = amp::floor<Precision>((x-a)/(b-a));
        x = x-k*(b-a);
        while( x<a )
        {
            x = x+(b-a);
            k = k-1;
        }
        while( x>b )
        {
            x = x-(b-a);
            k = k+1;
        }
        x = amp::maximum<Precision>(x, a);
        x = amp::minimum<Precision>(x, b);
    }
} // namespace

#endif
