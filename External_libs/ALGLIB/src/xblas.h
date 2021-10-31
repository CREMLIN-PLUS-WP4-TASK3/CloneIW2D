
#ifndef _xblas_h
#define _xblas_h

#include "ap.h"
#include "amp.h"
namespace xblas
{
    template<unsigned int Precision>
    void xdot(const ap::template_1d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& temp,
        amp::ampf<Precision>& r,
        amp::ampf<Precision>& rerr);
    template<unsigned int Precision>
    void xcdot(const ap::template_1d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& temp,
        amp::campf<Precision>& r,
        amp::ampf<Precision>& rerr);
    template<unsigned int Precision>
    void xsum(ap::template_1d_array< amp::ampf<Precision> >& w,
        amp::ampf<Precision> mx,
        int n,
        amp::ampf<Precision>& r,
        amp::ampf<Precision>& rerr);
    template<unsigned int Precision>
    amp::ampf<Precision> xfastpow(amp::ampf<Precision> r,
        int n);


    /*************************************************************************
    More precise dot-product. Absolute error of  subroutine  result  is  about
    1 ulp of max(MX,V), where:
        MX = max( |a[i]*b[i]| )
        V  = |(a,b)|

    INPUT PARAMETERS
        A       -   array[0..N-1], vector 1
        B       -   array[0..N-1], vector 2
        N       -   vectors length, N<2^29.
        Temp    -   array[0..N-1], pre-allocated temporary storage

    OUTPUT PARAMETERS
        R       -   (A,B)
        RErr    -   estimate of error. This estimate accounts for both  errors
                    during  calculation  of  (A,B)  and  errors  introduced by
                    rounding of A and B to fit in double (about 1 ulp).

      -- ALGLIB --
         Copyright 24.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void xdot(const ap::template_1d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& temp,
        amp::ampf<Precision>& r,
        amp::ampf<Precision>& rerr)
    {
        int i;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;


        
        //
        // special cases:
        // * N=0
        //
        if( n==0 )
        {
            r = 0;
            rerr = 0;
            return;
        }
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            v = a(i)*b(i);
            temp(i) = v;
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(v));
        }
        if( mx==0 )
        {
            r = 0;
            rerr = 0;
            return;
        }
        xsum<Precision>(temp, mx, n, r, rerr);
    }


    /*************************************************************************
    More precise complex dot-product. Absolute error of  subroutine  result is
    about 1 ulp of max(MX,V), where:
        MX = max( |a[i]*b[i]| )
        V  = |(a,b)|

    INPUT PARAMETERS
        A       -   array[0..N-1], vector 1
        B       -   array[0..N-1], vector 2
        N       -   vectors length, N<2^29.
        Temp    -   array[0..2*N-1], pre-allocated temporary storage

    OUTPUT PARAMETERS
        R       -   (A,B)
        RErr    -   estimate of error. This estimate accounts for both  errors
                    during  calculation  of  (A,B)  and  errors  introduced by
                    rounding of A and B to fit in double (about 1 ulp).

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void xcdot(const ap::template_1d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& temp,
        amp::campf<Precision>& r,
        amp::ampf<Precision>& rerr)
    {
        int i;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> rerrx;
        amp::ampf<Precision> rerry;


        
        //
        // special cases:
        // * N=0
        //
        if( n==0 )
        {
            r = 0;
            rerr = 0;
            return;
        }
        
        //
        // calculate real part
        //
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            v = a(i).x*b(i).x;
            temp(2*i+0) = v;
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(v));
            v = -a(i).y*b(i).y;
            temp(2*i+1) = v;
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(v));
        }
        if( mx==0 )
        {
            r.x = 0;
            rerrx = 0;
        }
        else
        {
            xsum<Precision>(temp, mx, 2*n, r.x, rerrx);
        }
        
        //
        // calculate imaginary part
        //
        mx = 0;
        for(i=0; i<=n-1; i++)
        {
            v = a(i).x*b(i).y;
            temp(2*i+0) = v;
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(v));
            v = a(i).y*b(i).x;
            temp(2*i+1) = v;
            mx = amp::maximum<Precision>(mx, amp::abs<Precision>(v));
        }
        if( mx==0 )
        {
            r.y = 0;
            rerry = 0;
        }
        else
        {
            xsum<Precision>(temp, mx, 2*n, r.y, rerry);
        }
        
        //
        // total error
        //
        if( rerrx==0 && rerry==0 )
        {
            rerr = 0;
        }
        else
        {
            rerr = amp::maximum<Precision>(rerrx, rerry)*amp::sqrt<Precision>(1+amp::sqr<Precision>(amp::minimum<Precision>(rerrx, rerry)/amp::maximum<Precision>(rerrx, rerry)));
        }
    }


    /*************************************************************************
    Internal subroutine for extra-precise calculation of SUM(w[i]).

    INPUT PARAMETERS:
        W   -   array[0..N-1], values to be added
                W is modified during calculations.
        MX  -   max(W[i])
        N   -   array size
        
    OUTPUT PARAMETERS:
        R   -   SUM(w[i])
        RErr-   error estimate for R

      -- ALGLIB --
         Copyright 24.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void xsum(ap::template_1d_array< amp::ampf<Precision> >& w,
        amp::ampf<Precision> mx,
        int n,
        amp::ampf<Precision>& r,
        amp::ampf<Precision>& rerr)
    {
        int i;
        int k;
        int ks;
        amp::ampf<Precision> v;
        amp::ampf<Precision> s;
        amp::ampf<Precision> ln2;
        amp::ampf<Precision> chunk;
        amp::ampf<Precision> invchunk;
        bool allzeros;


        
        //
        // special cases:
        // * N=0
        // * N is too large to use integer arithmetics
        //
        if( n==0 )
        {
            r = 0;
            rerr = 0;
            return;
        }
        if( mx==0 )
        {
            r = 0;
            rerr = 0;
            return;
        }
        ap::ap_error::make_assertion(n<536870912);
        
        //
        // Prepare
        //
        ln2 = amp::log<Precision>(amp::ampf<Precision>(2));
        rerr = mx*amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // 1. find S such that 0.5<=S*MX<1
        // 2. multiply W by S, so task is normalized in some sense
        // 3. S:=1/S so we can obtain original vector multiplying by S
        //
        k = amp::round<Precision>(amp::log<Precision>(mx)/ln2);
        s = xfastpow<Precision>(amp::ampf<Precision>(2), -k);
        while( s*mx>=1 )
        {
            s = amp::ampf<Precision>("0.5")*s;
        }
        while( s*mx<amp::ampf<Precision>("0.5") )
        {
            s = 2*s;
        }
        amp::vmul(w.getvector(0, n-1), s);
        s = 1/s;
        
        //
        // find Chunk=2^M such that N*Chunk<2^29
        //
        // we have chosen upper limit (2^29) with enough space left
        // to tolerate possible problems with rounding and N's close
        // to the limit, so we don't want to be very strict here.
        //
        k = amp::trunc<Precision>(amp::log<Precision>(amp::ampf<Precision>(536870912)/amp::ampf<Precision>(n))/ln2);
        chunk = xfastpow<Precision>(amp::ampf<Precision>(2), k);
        if( chunk<2 )
        {
            chunk = 2;
        }
        invchunk = 1/chunk;
        
        //
        // calculate result
        //
        r = 0;
        amp::vmul(w.getvector(0, n-1), chunk);
        while( true )
        {
            s = s*invchunk;
            allzeros = true;
            ks = 0;
            for(i=0; i<=n-1; i++)
            {
                v = w(i);
                k = amp::trunc<Precision>(v);
                if( v!=k )
                {
                    allzeros = false;
                }
                w(i) = chunk*(v-k);
                ks = ks+k;
            }
            r = r+s*ks;
            v = amp::abs<Precision>(r);
            if( allzeros || s*n+mx==mx )
            {
                break;
            }
        }
        
        //
        // correct error
        //
        rerr = amp::maximum<Precision>(rerr, amp::abs<Precision>(r)*amp::ampf<Precision>::getAlgoPascalEpsilon());
    }


    /*************************************************************************
    Fast Pow

      -- ALGLIB --
         Copyright 24.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> xfastpow(amp::ampf<Precision> r,
        int n)
    {
        amp::ampf<Precision> result;


        if( n>0 )
        {
            if( n%2==0 )
            {
                result = amp::sqr<Precision>(xfastpow<Precision>(r, n/2));
            }
            else
            {
                result = r*xfastpow<Precision>(r, n-1);
            }
            return result;
        }
        if( n==0 )
        {
            result = 1;
        }
        if( n<0 )
        {
            result = xfastpow<Precision>(1/r, -n);
        }
        return result;
    }
} // namespace

#endif
