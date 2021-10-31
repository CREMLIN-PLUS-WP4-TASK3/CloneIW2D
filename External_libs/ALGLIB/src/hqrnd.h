/*************************************************************************
Copyright (c)
    2007, Sergey Bochkanov (ALGLIB project).
    1988, Pierre L'Ecuyer

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

#ifndef _hqrnd_h
#define _hqrnd_h

#include "ap.h"
#include "amp.h"
namespace hqrnd
{
    /*************************************************************************
    Portable high quality random number generator state.
    Initialized with HQRNDRandomize() or HQRNDSeed().

    Fields:
        S1, S2      -   seed values
        V           -   precomputed value
        MagicV      -   'magic' value used to determine whether State structure
                        was correctly initialized.
    *************************************************************************/
    template<unsigned int Precision>
    class hqrndstate
    {
    public:
        int s1;
        int s2;
        amp::ampf<Precision> v;
        int magicv;
    };




    template<unsigned int Precision>
    void hqrndrandomize(hqrndstate<Precision>& state);
    template<unsigned int Precision>
    void hqrndseed(int s1,
        int s2,
        hqrndstate<Precision>& state);
    template<unsigned int Precision>
    amp::ampf<Precision> hqrnduniformr(hqrndstate<Precision>& state);
    template<unsigned int Precision>
    int hqrnduniformi(int n,
        hqrndstate<Precision>& state);
    template<unsigned int Precision>
    amp::ampf<Precision> hqrndnormal(hqrndstate<Precision>& state);
    template<unsigned int Precision>
    void hqrndunit2(hqrndstate<Precision>& state,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y);
    template<unsigned int Precision>
    void hqrndnormal2(hqrndstate<Precision>& state,
        amp::ampf<Precision>& x1,
        amp::ampf<Precision>& x2);
    template<unsigned int Precision>
    amp::ampf<Precision> hqrndexponential(amp::ampf<Precision> lambda,
        hqrndstate<Precision>& state);
    template<unsigned int Precision>
    int hqrndintegerbase(hqrndstate<Precision>& state);


    static const int hqrndmax = 2147483563;
    static const int hqrndm1 = 2147483563;
    static const int hqrndm2 = 2147483399;
    static const int hqrndmagic = 1634357784;


    /*************************************************************************
    HQRNDState  initialization  with  random  values  which come from standard
    RNG.

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hqrndrandomize(hqrndstate<Precision>& state)
    {
        hqrndseed<Precision>(ap::randominteger(hqrndm1), ap::randominteger(hqrndm2), state);
    }


    /*************************************************************************
    HQRNDState initialization with seed values

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hqrndseed(int s1,
        int s2,
        hqrndstate<Precision>& state)
    {
        state.s1 = s1%(hqrndm1-1)+1;
        state.s2 = s2%(hqrndm2-1)+1;
        state.v = amp::ampf<Precision>(1)/amp::ampf<Precision>(hqrndmax);
        state.magicv = hqrndmagic;
    }


    /*************************************************************************
    This function generates random real number in (0,1),
    not including interval boundaries

    State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> hqrnduniformr(hqrndstate<Precision>& state)
    {
        amp::ampf<Precision> result;


        result = state.v*hqrndintegerbase<Precision>(state);
        return result;
    }


    /*************************************************************************
    This function generates random integer number in [0, N)

    1. N must be less than HQRNDMax-1.
    2. State structure must be initialized with HQRNDRandomize() or HQRNDSeed()

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int hqrnduniformi(int n,
        hqrndstate<Precision>& state)
    {
        int result;
        int mx;


        
        //
        // Correct handling of N's close to RNDBaseMax
        // (avoiding skewed distributions for RNDBaseMax<>K*N)
        //
        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(n<hqrndmax-1);
        mx = hqrndmax-1-(hqrndmax-1)%n;
        do
        {
            result = hqrndintegerbase<Precision>(state)-1;
        }
        while( result>=mx );
        result = result%n;
        return result;
    }


    /*************************************************************************
    Random number generator: normal numbers

    This function generates one random number from normal distribution.
    Its performance is equal to that of HQRNDNormal2()

    State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> hqrndnormal(hqrndstate<Precision>& state)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;


        hqrndnormal2<Precision>(state, v1, v2);
        result = v1;
        return result;
    }


    /*************************************************************************
    Random number generator: random X and Y such that X^2+Y^2=1

    State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hqrndunit2(hqrndstate<Precision>& state,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> mn;


        do
        {
            hqrndnormal2<Precision>(state, x, y);
        }
        while( ! (x!=0 || y!=0) );
        mx = amp::maximum<Precision>(amp::abs<Precision>(x), amp::abs<Precision>(y));
        mn = amp::minimum<Precision>(amp::abs<Precision>(x), amp::abs<Precision>(y));
        v = mx*amp::sqrt<Precision>(1+amp::sqr<Precision>(mn/mx));
        x = x/v;
        y = y/v;
    }


    /*************************************************************************
    Random number generator: normal numbers

    This function generates two independent random numbers from normal
    distribution. Its performance is equal to that of HQRNDNormal()

    State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

      -- ALGLIB --
         Copyright 02.12.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hqrndnormal2(hqrndstate<Precision>& state,
        amp::ampf<Precision>& x1,
        amp::ampf<Precision>& x2)
    {
        amp::ampf<Precision> u;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;


        while( true )
        {
            u = 2*hqrnduniformr<Precision>(state)-1;
            v = 2*hqrnduniformr<Precision>(state)-1;
            s = amp::sqr<Precision>(u)+amp::sqr<Precision>(v);
            if( s>0 && s<1 )
            {
                
                //
                // two Sqrt's instead of one to
                // avoid overflow when S is too small
                //
                s = amp::sqrt<Precision>(-2*amp::log<Precision>(s))/amp::sqrt<Precision>(s);
                x1 = u*s;
                x2 = v*s;
                return;
            }
        }
    }


    /*************************************************************************
    Random number generator: exponential distribution

    State structure must be initialized with HQRNDRandomize() or HQRNDSeed().

      -- ALGLIB --
         Copyright 11.08.2007 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> hqrndexponential(amp::ampf<Precision> lambda,
        hqrndstate<Precision>& state)
    {
        amp::ampf<Precision> result;


        ap::ap_error::make_assertion(lambda>0);
        result = -amp::log<Precision>(hqrnduniformr<Precision>(state))/lambda;
        return result;
    }


    /*************************************************************************

    L'Ecuyer, Efficient and portable combined random number generators
    *************************************************************************/
    template<unsigned int Precision>
    int hqrndintegerbase(hqrndstate<Precision>& state)
    {
        int result;
        int k;


        ap::ap_error::make_assertion(state.magicv==hqrndmagic);
        k = state.s1/53668;
        state.s1 = 40014*(state.s1-k*53668)-k*12211;
        if( state.s1<0 )
        {
            state.s1 = state.s1+2147483563;
        }
        k = state.s2/52774;
        state.s2 = 40692*(state.s2-k*52774)-k*3791;
        if( state.s2<0 )
        {
            state.s2 = state.s2+2147483399;
        }
        
        //
        // Result
        //
        result = state.s1-state.s2;
        if( result<1 )
        {
            result = result+2147483562;
        }
        return result;
    }
} // namespace

#endif
