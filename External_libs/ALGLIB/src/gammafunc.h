/*************************************************************************
Cephes Math Library Release 2.8:  June, 2000
Copyright by Stephen L. Moshier

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from C to
      pseudocode.

See subroutines comments for additional copyrights.

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

#ifndef _gammafunc_h
#define _gammafunc_h

#include "ap.h"
#include "amp.h"
namespace gammafunc
{
    template<unsigned int Precision>
    amp::ampf<Precision> gamma(amp::ampf<Precision> x);
    template<unsigned int Precision>
    amp::ampf<Precision> lngamma(amp::ampf<Precision> x,
        amp::ampf<Precision>& sgngam);
    template<unsigned int Precision>
    amp::ampf<Precision> gammastirf(amp::ampf<Precision> x);


    /*************************************************************************
    Gamma function

    Input parameters:
        X   -   argument

    Domain:
        0 < X < 171.6
        -170 < X < 0, X is not an integer.

    Relative error:
     arithmetic   domain     # trials      peak         rms
        IEEE    -170,-33      20000       2.3e-15     3.3e-16
        IEEE     -33,  33     20000       9.4e-16     2.2e-16
        IEEE      33, 171.6   20000       2.3e-15     3.2e-16

    Cephes Math Library Release 2.8:  June, 2000
    Original copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
    Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> gamma(amp::ampf<Precision> x)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_SPECFUNCS
        amp::ampf<Precision> result;
        amp::ampf<Precision> p;
        amp::ampf<Precision> pp;
        amp::ampf<Precision> q;
        amp::ampf<Precision> qq;
        amp::ampf<Precision> z;
        int i;
        amp::ampf<Precision> sgngam;


        sgngam = 1;
        q = amp::abs<Precision>(x);
        if( q>amp::ampf<Precision>("33.0") )
        {
            if( x<amp::ampf<Precision>("0.0") )
            {
                p = amp::floor<Precision>(q);
                i = amp::round<Precision>(p);
                if( i%2==0 )
                {
                    sgngam = -1;
                }
                z = q-p;
                if( z>amp::ampf<Precision>("0.5") )
                {
                    p = p+1;
                    z = q-p;
                }
                z = q*amp::sin<Precision>(amp::pi<Precision>()*z);
                z = amp::abs<Precision>(z);
                z = amp::pi<Precision>()/(z*gammastirf<Precision>(q));
            }
            else
            {
                z = gammastirf<Precision>(x);
            }
            result = sgngam*z;
            return result;
        }
        z = 1;
        while( x>=3 )
        {
            x = x-1;
            z = z*x;
        }
        while( x<0 )
        {
            if( x>-amp::ampf<Precision>("0.000000001") )
            {
                result = z/((1+amp::ampf<Precision>("0.5772156649015329")*x)*x);
                return result;
            }
            z = z/x;
            x = x+1;
        }
        while( x<2 )
        {
            if( x<amp::ampf<Precision>("0.000000001") )
            {
                result = z/((1+amp::ampf<Precision>("0.5772156649015329")*x)*x);
                return result;
            }
            z = z/x;
            x = x+amp::ampf<Precision>("1.0");
        }
        if( x==2 )
        {
            result = z;
            return result;
        }
        x = x-amp::ampf<Precision>("2.0");
        pp = amp::ampf<Precision>("1.60119522476751861407E-4");
        pp = amp::ampf<Precision>("1.19135147006586384913E-3")+x*pp;
        pp = amp::ampf<Precision>("1.04213797561761569935E-2")+x*pp;
        pp = amp::ampf<Precision>("4.76367800457137231464E-2")+x*pp;
        pp = amp::ampf<Precision>("2.07448227648435975150E-1")+x*pp;
        pp = amp::ampf<Precision>("4.94214826801497100753E-1")+x*pp;
        pp = amp::ampf<Precision>("9.99999999999999996796E-1")+x*pp;
        qq = -amp::ampf<Precision>("2.31581873324120129819E-5");
        qq = amp::ampf<Precision>("5.39605580493303397842E-4")+x*qq;
        qq = -amp::ampf<Precision>("4.45641913851797240494E-3")+x*qq;
        qq = amp::ampf<Precision>("1.18139785222060435552E-2")+x*qq;
        qq = amp::ampf<Precision>("3.58236398605498653373E-2")+x*qq;
        qq = -amp::ampf<Precision>("2.34591795718243348568E-1")+x*qq;
        qq = amp::ampf<Precision>("7.14304917030273074085E-2")+x*qq;
        qq = amp::ampf<Precision>("1.00000000000000000320")+x*qq;
        result = z*pp/qq;
        return result;
        return result;
    #else
        return amp::_i_gamma<Precision>(x);
    #endif
    }


    /*************************************************************************
    Natural logarithm of gamma function

    Input parameters:
        X       -   argument

    Result:
        logarithm of the absolute value of the Gamma(X).

    Output parameters:
        SgnGam  -   sign(Gamma(X))

    Domain:
        0 < X < 2.55e305
        -2.55e305 < X < 0, X is not an integer.

    ACCURACY:
    arithmetic      domain        # trials     peak         rms
       IEEE    0, 3                 28000     5.4e-16     1.1e-16
       IEEE    2.718, 2.556e305     40000     3.5e-16     8.3e-17
    The error criterion was relative when the function magnitude
    was greater than one but absolute when it was less than one.

    The following test used the relative error criterion, though
    at certain points the relative error could be much higher than
    indicated.
       IEEE    -200, -4             10000     4.8e-16     1.3e-16

    Cephes Math Library Release 2.8:  June, 2000
    Copyright 1984, 1987, 1989, 1992, 2000 by Stephen L. Moshier
    Translated to AlgoPascal by Bochkanov Sergey (2005, 2006, 2007).
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> lngamma(amp::ampf<Precision> x,
        amp::ampf<Precision>& sgngam)
    {
    #ifndef ALGLIB_MP_INTERCEPTS_SPECFUNCS
        amp::ampf<Precision> result;
        amp::ampf<Precision> a;
        amp::ampf<Precision> b;
        amp::ampf<Precision> c;
        amp::ampf<Precision> p;
        amp::ampf<Precision> q;
        amp::ampf<Precision> u;
        amp::ampf<Precision> w;
        amp::ampf<Precision> z;
        int i;
        amp::ampf<Precision> logpi;
        amp::ampf<Precision> ls2pi;
        amp::ampf<Precision> tmp;


        sgngam = 1;
        logpi = amp::ampf<Precision>("1.14472988584940017414");
        ls2pi = amp::ampf<Precision>("0.91893853320467274178");
        if( x<-amp::ampf<Precision>("34.0") )
        {
            q = -x;
            w = lngamma<Precision>(q, tmp);
            p = amp::floor<Precision>(q);
            i = amp::round<Precision>(p);
            if( i%2==0 )
            {
                sgngam = -1;
            }
            else
            {
                sgngam = 1;
            }
            z = q-p;
            if( z>amp::ampf<Precision>("0.5") )
            {
                p = p+1;
                z = p-q;
            }
            z = q*amp::sin<Precision>(amp::pi<Precision>()*z);
            result = logpi-amp::log<Precision>(z)-w;
            return result;
        }
        if( x<13 )
        {
            z = 1;
            p = 0;
            u = x;
            while( u>=3 )
            {
                p = p-1;
                u = x+p;
                z = z*u;
            }
            while( u<2 )
            {
                z = z/u;
                p = p+1;
                u = x+p;
            }
            if( z<0 )
            {
                sgngam = -1;
                z = -z;
            }
            else
            {
                sgngam = 1;
            }
            if( u==2 )
            {
                result = amp::log<Precision>(z);
                return result;
            }
            p = p-2;
            x = x+p;
            b = -amp::ampf<Precision>("1378.25152569120859100");
            b = -amp::ampf<Precision>("38801.6315134637840924")+x*b;
            b = -amp::ampf<Precision>("331612.992738871184744")+x*b;
            b = -amp::ampf<Precision>("1162370.97492762307383")+x*b;
            b = -amp::ampf<Precision>("1721737.00820839662146")+x*b;
            b = -amp::ampf<Precision>("853555.664245765465627")+x*b;
            c = 1;
            c = -amp::ampf<Precision>("351.815701436523470549")+x*c;
            c = -amp::ampf<Precision>("17064.2106651881159223")+x*c;
            c = -amp::ampf<Precision>("220528.590553854454839")+x*c;
            c = -amp::ampf<Precision>("1139334.44367982507207")+x*c;
            c = -amp::ampf<Precision>("2532523.07177582951285")+x*c;
            c = -amp::ampf<Precision>("2018891.41433532773231")+x*c;
            p = x*b/c;
            result = amp::log<Precision>(z)+p;
            return result;
        }
        q = (x-amp::ampf<Precision>("0.5"))*amp::log<Precision>(x)-x+ls2pi;
        if( x>100000000 )
        {
            result = q;
            return result;
        }
        p = 1/(x*x);
        if( x>=amp::ampf<Precision>("1000.0") )
        {
            q = q+((amp::ampf<Precision>("7.9365079365079365079365")*amp::ampf<Precision>("0.0001")*p-amp::ampf<Precision>("2.7777777777777777777778")*amp::ampf<Precision>("0.001"))*p+amp::ampf<Precision>("0.0833333333333333333333"))/x;
        }
        else
        {
            a = amp::ampf<Precision>("8.11614167470508450300")*amp::ampf<Precision>("0.0001");
            a = -amp::ampf<Precision>("5.95061904284301438324")*amp::ampf<Precision>("0.0001")+p*a;
            a = amp::ampf<Precision>("7.93650340457716943945")*amp::ampf<Precision>("0.0001")+p*a;
            a = -amp::ampf<Precision>("2.77777777730099687205")*amp::ampf<Precision>("0.001")+p*a;
            a = amp::ampf<Precision>("8.33333333333331927722")*amp::ampf<Precision>("0.01")+p*a;
            q = q+a/x;
        }
        result = q;
        return result;
    #else
        return amp::_i_lngamma<Precision>(x, sgngam);
    #endif
    }


    template<unsigned int Precision>
    amp::ampf<Precision> gammastirf(amp::ampf<Precision> x)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> y;
        amp::ampf<Precision> w;
        amp::ampf<Precision> v;
        amp::ampf<Precision> stir;


        w = 1/x;
        stir = amp::ampf<Precision>("7.87311395793093628397E-4");
        stir = -amp::ampf<Precision>("2.29549961613378126380E-4")+w*stir;
        stir = -amp::ampf<Precision>("2.68132617805781232825E-3")+w*stir;
        stir = amp::ampf<Precision>("3.47222221605458667310E-3")+w*stir;
        stir = amp::ampf<Precision>("8.33333333333482257126E-2")+w*stir;
        w = 1+w*stir;
        y = amp::exp<Precision>(x);
        if( x>amp::ampf<Precision>("143.01608") )
        {
            v = amp::pow<Precision>(x, amp::ampf<Precision>("0.5")*x-amp::ampf<Precision>("0.25"));
            y = v*(v/y);
        }
        else
        {
            y = amp::pow<Precision>(x, x-amp::ampf<Precision>("0.5"))/y;
        }
        result = amp::ampf<Precision>("2.50662827463100050242")*y*w;
        return result;
    }
} // namespace

#endif
