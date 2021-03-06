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

#ifndef _corr_h
#define _corr_h

#include "ap.h"
#include "amp.h"
#include "ftbase.h"
#include "fft.h"
#include "conv.h"
namespace corr
{
    template<unsigned int Precision>
    void corrc1d(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r);
    template<unsigned int Precision>
    void corrc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& c);
    template<unsigned int Precision>
    void corrr1d(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r);
    template<unsigned int Precision>
    void corrr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c);


    /*************************************************************************
    1-dimensional complex cross-correlation.

    For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

    Correlation is calculated using reduction to  convolution.  Algorithm with
    max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
    about performance).

    IMPORTANT:
        for  historical reasons subroutine accepts its parameters in  reversed
        order: CorrC1D(Signal, Pattern) = Pattern x Signal (using  traditional
        definition of cross-correlation, denoting cross-correlation as "x").

    INPUT PARAMETERS
        Signal  -   array[0..N-1] - complex function to be transformed,
                    signal containing pattern
        N       -   problem size
        Pattern -   array[0..M-1] - complex function to be transformed,
                    pattern to search withing signal
        M       -   problem size

    OUTPUT PARAMETERS
        R       -   cross-correlation, array[0..N+M-2]:
                    * positive lags are stored in R[0..N-1],
                      R[i] = sum(conj(pattern[j])*signal[i+j]
                    * negative lags are stored in R[N..N+M-2],
                      R[N+M-1-i] = sum(conj(pattern[j])*signal[-i+j]

    NOTE:
        It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
    on [-K..M-1],  you can still use this subroutine, just shift result by K.

      -- ALGLIB --
         Copyright 21.07.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void corrc1d(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::campf<Precision> >& r)
    {
        ap::template_1d_array< amp::campf<Precision> > p;
        ap::template_1d_array< amp::campf<Precision> > b;
        int i;
        int i_;
        int i1_;


        ap::ap_error::make_assertion(n>0 && m>0);
        p.setlength(m);
        for(i=0; i<=m-1; i++)
        {
            p(m-1-i) = amp::conj<Precision>(pattern(i));
        }
        conv::convc1d<Precision>(p, m, signal, n, b);
        r.setlength(m+n-1);
        i1_ = (m-1) - (0);
        for(i_=0; i_<=n-1;i_++)
        {
            r(i_) = b(i_+i1_);
        }
        if( m+n-2>=n )
        {
            i1_ = (0) - (n);
            for(i_=n; i_<=m+n-2;i_++)
            {
                r(i_) = b(i_+i1_);
            }
        }
    }


    /*************************************************************************
    1-dimensional circular complex cross-correlation.

    For given Pattern/Signal returns corr(Pattern,Signal) (circular).
    Algorithm has linearithmic complexity for any M/N.

    IMPORTANT:
        for  historical reasons subroutine accepts its parameters in  reversed
        order:   CorrC1DCircular(Signal, Pattern) = Pattern x Signal    (using
        traditional definition of cross-correlation, denoting cross-correlation
        as "x").

    INPUT PARAMETERS
        Signal  -   array[0..N-1] - complex function to be transformed,
                    periodic signal containing pattern
        N       -   problem size
        Pattern -   array[0..M-1] - complex function to be transformed,
                    non-periodic pattern to search withing signal
        M       -   problem size

    OUTPUT PARAMETERS
        R   -   convolution: A*B. array[0..M-1].


      -- ALGLIB --
         Copyright 21.07.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void corrc1dcircular(const ap::template_1d_array< amp::campf<Precision> >& signal,
        int m,
        const ap::template_1d_array< amp::campf<Precision> >& pattern,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& c)
    {
        ap::template_1d_array< amp::campf<Precision> > p;
        ap::template_1d_array< amp::campf<Precision> > b;
        int i1;
        int i2;
        int i;
        int j2;
        int i_;
        int i1_;


        ap::ap_error::make_assertion(n>0 && m>0);
        
        //
        // normalize task: make M>=N,
        // so A will be longer (at least - not shorter) that B.
        //
        if( m<n )
        {
            b.setlength(m);
            for(i1=0; i1<=m-1; i1++)
            {
                b(i1) = 0;
            }
            i1 = 0;
            while( i1<n )
            {
                i2 = ap::minint(i1+m-1, n-1);
                j2 = i2-i1;
                i1_ = (i1) - (0);
                for(i_=0; i_<=j2;i_++)
                {
                    b(i_) = b(i_) + pattern(i_+i1_);
                }
                i1 = i1+m;
            }
            corrc1dcircular<Precision>(signal, m, b, m, c);
            return;
        }
        
        //
        // Task is normalized
        //
        p.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            p(n-1-i) = amp::conj<Precision>(pattern(i));
        }
        conv::convc1dcircular<Precision>(signal, m, p, n, b);
        c.setlength(m);
        i1_ = (n-1) - (0);
        for(i_=0; i_<=m-n;i_++)
        {
            c(i_) = b(i_+i1_);
        }
        if( m-n+1<=m-1 )
        {
            i1_ = (0) - (m-n+1);
            for(i_=m-n+1; i_<=m-1;i_++)
            {
                c(i_) = b(i_+i1_);
            }
        }
    }


    /*************************************************************************
    1-dimensional real cross-correlation.

    For given Pattern/Signal returns corr(Pattern,Signal) (non-circular).

    Correlation is calculated using reduction to  convolution.  Algorithm with
    max(N,N)*log(max(N,N)) complexity is used (see  ConvC1D()  for  more  info
    about performance).

    IMPORTANT:
        for  historical reasons subroutine accepts its parameters in  reversed
        order: CorrR1D(Signal, Pattern) = Pattern x Signal (using  traditional
        definition of cross-correlation, denoting cross-correlation as "x").

    INPUT PARAMETERS
        Signal  -   array[0..N-1] - real function to be transformed,
                    signal containing pattern
        N       -   problem size
        Pattern -   array[0..M-1] - real function to be transformed,
                    pattern to search withing signal
        M       -   problem size

    OUTPUT PARAMETERS
        R       -   cross-correlation, array[0..N+M-2]:
                    * positive lags are stored in R[0..N-1],
                      R[i] = sum(pattern[j]*signal[i+j]
                    * negative lags are stored in R[N..N+M-2],
                      R[N+M-1-i] = sum(pattern[j]*signal[-i+j]

    NOTE:
        It is assumed that pattern domain is [0..M-1].  If Pattern is non-zero
    on [-K..M-1],  you can still use this subroutine, just shift result by K.

      -- ALGLIB --
         Copyright 21.07.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void corrr1d(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int m,
        ap::template_1d_array< amp::ampf<Precision> >& r)
    {
        ap::template_1d_array< amp::ampf<Precision> > p;
        ap::template_1d_array< amp::ampf<Precision> > b;
        int i;


        ap::ap_error::make_assertion(n>0 && m>0);
        p.setlength(m);
        for(i=0; i<=m-1; i++)
        {
            p(m-1-i) = pattern(i);
        }
        conv::convr1d<Precision>(p, m, signal, n, b);
        r.setlength(m+n-1);
        amp::vmove(r.getvector(0, n-1), b.getvector(m-1, m+n-2));
        if( m+n-2>=n )
        {
            amp::vmove(r.getvector(n, m+n-2), b.getvector(0, m-2));
        }
    }


    /*************************************************************************
    1-dimensional circular real cross-correlation.

    For given Pattern/Signal returns corr(Pattern,Signal) (circular).
    Algorithm has linearithmic complexity for any M/N.

    IMPORTANT:
        for  historical reasons subroutine accepts its parameters in  reversed
        order:   CorrR1DCircular(Signal, Pattern) = Pattern x Signal    (using
        traditional definition of cross-correlation, denoting cross-correlation
        as "x").

    INPUT PARAMETERS
        Signal  -   array[0..N-1] - real function to be transformed,
                    periodic signal containing pattern
        N       -   problem size
        Pattern -   array[0..M-1] - real function to be transformed,
                    non-periodic pattern to search withing signal
        M       -   problem size

    OUTPUT PARAMETERS
        R   -   convolution: A*B. array[0..M-1].


      -- ALGLIB --
         Copyright 21.07.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void corrr1dcircular(const ap::template_1d_array< amp::ampf<Precision> >& signal,
        int m,
        const ap::template_1d_array< amp::ampf<Precision> >& pattern,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& c)
    {
        ap::template_1d_array< amp::ampf<Precision> > p;
        ap::template_1d_array< amp::ampf<Precision> > b;
        int i1;
        int i2;
        int i;
        int j2;


        ap::ap_error::make_assertion(n>0 && m>0);
        
        //
        // normalize task: make M>=N,
        // so A will be longer (at least - not shorter) that B.
        //
        if( m<n )
        {
            b.setlength(m);
            for(i1=0; i1<=m-1; i1++)
            {
                b(i1) = 0;
            }
            i1 = 0;
            while( i1<n )
            {
                i2 = ap::minint(i1+m-1, n-1);
                j2 = i2-i1;
                amp::vadd(b.getvector(0, j2), pattern.getvector(i1, i2));
                i1 = i1+m;
            }
            corrr1dcircular<Precision>(signal, m, b, m, c);
            return;
        }
        
        //
        // Task is normalized
        //
        p.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            p(n-1-i) = pattern(i);
        }
        conv::convr1dcircular<Precision>(signal, m, p, n, b);
        c.setlength(m);
        amp::vmove(c.getvector(0, m-n), b.getvector(n-1, m-1));
        if( m-n+1<=m-1 )
        {
            amp::vmove(c.getvector(m-n+1, m-1), b.getvector(0, n-2));
        }
    }
} // namespace

#endif
