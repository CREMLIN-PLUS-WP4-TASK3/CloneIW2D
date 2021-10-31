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

#ifndef _matgen_h
#define _matgen_h

#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
namespace matgen
{
    template<unsigned int Precision>
    void rmatrixrndorthogonal(int n,
        ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void rmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void cmatrixrndorthogonal(int n,
        ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void cmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void smatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void spdmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a);
    template<unsigned int Precision>
    void hmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void hpdmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a);
    template<unsigned int Precision>
    void rmatrixrndorthogonalfromtheright(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n);
    template<unsigned int Precision>
    void rmatrixrndorthogonalfromtheleft(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n);
    template<unsigned int Precision>
    void cmatrixrndorthogonalfromtheright(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n);
    template<unsigned int Precision>
    void cmatrixrndorthogonalfromtheleft(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n);
    template<unsigned int Precision>
    void smatrixrndmultiply(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    void hmatrixrndmultiply(ap::template_2d_array< amp::campf<Precision> >& a,
        int n);


    /*************************************************************************
    Generation of a random uniformly distributed (Haar) orthogonal matrix

    INPUT PARAMETERS:
        N   -   matrix size, N>=1
        
    OUTPUT PARAMETERS:
        A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrndorthogonal(int n,
        ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        int i;
        int j;


        ap::ap_error::make_assertion(n>=1);
        a.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( i==j )
                {
                    a(i,j) = 1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
        rmatrixrndorthogonalfromtheright<Precision>(a, n, n);
    }


    /*************************************************************************
    Generation of random NxN matrix with given condition number and norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;


        ap::ap_error::make_assertion(n>=1 && c>=1);
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            
            //
            // special case
            //
            a(0,0) = 2*ap::randominteger(2)-1;
            return;
        }
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        rmatrixrndorthogonalfromtheleft<Precision>(a, n, n);
        rmatrixrndorthogonalfromtheright<Precision>(a, n, n);
    }


    /*************************************************************************
    Generation of a random Haar distributed orthogonal complex matrix

    INPUT PARAMETERS:
        N   -   matrix size, N>=1

    OUTPUT PARAMETERS:
        A   -   orthogonal NxN matrix, array[0..N-1,0..N-1]

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrndorthogonal(int n,
        ap::template_2d_array< amp::campf<Precision> >& a)
    {
        int i;
        int j;


        ap::ap_error::make_assertion(n>=1);
        a.setbounds(0, n-1, 0, n-1);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                if( i==j )
                {
                    a(i,j) = 1;
                }
                else
                {
                    a(i,j) = 0;
                }
            }
        }
        cmatrixrndorthogonalfromtheright<Precision>(a, n, n);
    }


    /*************************************************************************
    Generation of random NxN complex matrix with given condition number C and
    norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;
        hqrnd::hqrndstate<Precision> state;
        amp::campf<Precision> v;


        ap::ap_error::make_assertion(n>=1 && c>=1);
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            
            //
            // special case
            //
            hqrnd::hqrndrandomize<Precision>(state);
            hqrnd::hqrndunit2<Precision>(state, v.x, v.y);
            a(0,0) = v;
            return;
        }
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        cmatrixrndorthogonalfromtheleft<Precision>(a, n, n);
        cmatrixrndorthogonalfromtheright<Precision>(a, n, n);
    }


    /*************************************************************************
    Generation of random NxN symmetric matrix with given condition number  and
    norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void smatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;


        ap::ap_error::make_assertion(n>=1 && c>=1);
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            
            //
            // special case
            //
            a(0,0) = 2*ap::randominteger(2)-1;
            return;
        }
        
        //
        // Prepare matrix
        //
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = (2*ap::randominteger(2)-1)*amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        
        //
        // Multiply
        //
        smatrixrndmultiply<Precision>(a, n);
    }


    /*************************************************************************
    Generation of random NxN symmetric positive definite matrix with given
    condition number and norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random SPD matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::ampf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;


        
        //
        // Special cases
        //
        if( n<=0 || c<1 )
        {
            return;
        }
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            a(0,0) = 1;
            return;
        }
        
        //
        // Prepare matrix
        //
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        
        //
        // Multiply
        //
        smatrixrndmultiply<Precision>(a, n);
    }


    /*************************************************************************
    Generation of random NxN Hermitian matrix with given condition number  and
    norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;


        ap::ap_error::make_assertion(n>=1 && c>=1);
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            
            //
            // special case
            //
            a(0,0) = 2*ap::randominteger(2)-1;
            return;
        }
        
        //
        // Prepare matrix
        //
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = (2*ap::randominteger(2)-1)*amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        
        //
        // Multiply
        //
        hmatrixrndmultiply<Precision>(a, n);
        
        //
        // post-process to ensure that matrix diagonal is real
        //
        for(i=0; i<=n-1; i++)
        {
            a(i,i).y = 0;
        }
    }


    /*************************************************************************
    Generation of random NxN Hermitian positive definite matrix with given
    condition number and norm2(A)=1

    INPUT PARAMETERS:
        N   -   matrix size
        C   -   condition number (in 2-norm)

    OUTPUT PARAMETERS:
        A   -   random HPD matrix with norm2(A)=1 and cond(A)=C

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixrndcond(int n,
        amp::ampf<Precision> c,
        ap::template_2d_array< amp::campf<Precision> >& a)
    {
        int i;
        int j;
        amp::ampf<Precision> l1;
        amp::ampf<Precision> l2;


        
        //
        // Special cases
        //
        if( n<=0 || c<1 )
        {
            return;
        }
        a.setbounds(0, n-1, 0, n-1);
        if( n==1 )
        {
            a(0,0) = 1;
            return;
        }
        
        //
        // Prepare matrix
        //
        l1 = 0;
        l2 = amp::log<Precision>(1/c);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                a(i,j) = 0;
            }
        }
        a(0,0) = amp::exp<Precision>(l1);
        for(i=1; i<=n-2; i++)
        {
            a(i,i) = amp::exp<Precision>(amp::ampf<Precision>::getRandom()*(l2-l1)+l1);
        }
        a(n-1,n-1) = amp::exp<Precision>(l2);
        
        //
        // Multiply
        //
        hmatrixrndmultiply<Precision>(a, n);
        
        //
        // post-process to ensure that matrix diagonal is real
        //
        for(i=0; i<=n-1; i++)
        {
            a(i,i).y = 0;
        }
    }


    /*************************************************************************
    Multiplication of MxN matrix by NxN random Haar distributed orthogonal matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..M-1, 0..N-1]
        M, N-   matrix size

    OUTPUT PARAMETERS:
        A   -   A*Q, where Q is random NxN orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrndorthogonalfromtheright(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n)
    {
        amp::ampf<Precision> tau;
        amp::ampf<Precision> lambda;
        int s;
        int i;
        amp::ampf<Precision> u1;
        amp::ampf<Precision> u2;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;


        ap::ap_error::make_assertion(n>=1 && m>=1);
        if( n==1 )
        {
            
            //
            // Special case
            //
            tau = 2*ap::randominteger(2)-1;
            for(i=0; i<=m-1; i++)
            {
                a(i,0) = a(i,0)*tau;
            }
            return;
        }
        
        //
        // General case.
        // First pass.
        //
        w.setbounds(0, m-1);
        v.setbounds(1, n);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=n; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                i = 1;
                while( i<=s )
                {
                    hqrnd::hqrndnormal2<Precision>(state, u1, u2);
                    v(i) = u1;
                    if( i+1<=s )
                    {
                        v(i+1) = u2;
                    }
                    i = i+2;
                }
                lambda = amp::vdotproduct(v.getvector(1, s), v.getvector(1, s));
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            reflections::generatereflection<Precision>(v, s, tau);
            v(1) = 1;
            reflections::applyreflectionfromtheright<Precision>(a, tau, v, 0, m-1, n-s, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=n-1; i++)
        {
            tau = 2*ap::randominteger(2)-1;
            amp::vmul(a.getcolumn(i, 0, m-1), tau);
        }
    }


    /*************************************************************************
    Multiplication of MxN matrix by MxM random Haar distributed orthogonal matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..M-1, 0..N-1]
        M, N-   matrix size

    OUTPUT PARAMETERS:
        A   -   Q*A, where Q is random MxM orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrndorthogonalfromtheleft(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n)
    {
        amp::ampf<Precision> tau;
        amp::ampf<Precision> lambda;
        int s;
        int i;
        int j;
        amp::ampf<Precision> u1;
        amp::ampf<Precision> u2;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;


        ap::ap_error::make_assertion(n>=1 && m>=1);
        if( m==1 )
        {
            
            //
            // special case
            //
            tau = 2*ap::randominteger(2)-1;
            for(j=0; j<=n-1; j++)
            {
                a(0,j) = a(0,j)*tau;
            }
            return;
        }
        
        //
        // General case.
        // First pass.
        //
        w.setbounds(0, n-1);
        v.setbounds(1, m);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=m; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                i = 1;
                while( i<=s )
                {
                    hqrnd::hqrndnormal2<Precision>(state, u1, u2);
                    v(i) = u1;
                    if( i+1<=s )
                    {
                        v(i+1) = u2;
                    }
                    i = i+2;
                }
                lambda = amp::vdotproduct(v.getvector(1, s), v.getvector(1, s));
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            reflections::generatereflection<Precision>(v, s, tau);
            v(1) = 1;
            reflections::applyreflectionfromtheleft<Precision>(a, tau, v, m-s, m-1, 0, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=m-1; i++)
        {
            tau = 2*ap::randominteger(2)-1;
            amp::vmul(a.getrow(i, 0, n-1), tau);
        }
    }


    /*************************************************************************
    Multiplication of MxN complex matrix by NxN random Haar distributed
    complex orthogonal matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..M-1, 0..N-1]
        M, N-   matrix size

    OUTPUT PARAMETERS:
        A   -   A*Q, where Q is random NxN orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrndorthogonalfromtheright(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n)
    {
        amp::campf<Precision> lambda;
        amp::campf<Precision> tau;
        int s;
        int i;
        ap::template_1d_array< amp::campf<Precision> > w;
        ap::template_1d_array< amp::campf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;
        int i_;


        ap::ap_error::make_assertion(n>=1 && m>=1);
        if( n==1 )
        {
            
            //
            // Special case
            //
            hqrnd::hqrndrandomize<Precision>(state);
            hqrnd::hqrndunit2<Precision>(state, tau.x, tau.y);
            for(i=0; i<=m-1; i++)
            {
                a(i,0) = a(i,0)*tau;
            }
            return;
        }
        
        //
        // General case.
        // First pass.
        //
        w.setbounds(0, m-1);
        v.setbounds(1, n);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=n; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                for(i=1; i<=s; i++)
                {
                    hqrnd::hqrndnormal2<Precision>(state, tau.x, tau.y);
                    v(i) = tau;
                }
                lambda = 0.0;
                for(i_=1; i_<=s;i_++)
                {
                    lambda += v(i_)*amp::conj(v(i_));
                }
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            creflections::complexgeneratereflection<Precision>(v, s, tau);
            v(1) = 1;
            creflections::complexapplyreflectionfromtheright<Precision>(a, tau, v, 0, m-1, n-s, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=n-1; i++)
        {
            hqrnd::hqrndunit2<Precision>(state, tau.x, tau.y);
            for(i_=0; i_<=m-1;i_++)
            {
                a(i_,i) = tau*a(i_,i);
            }
        }
    }


    /*************************************************************************
    Multiplication of MxN complex matrix by MxM random Haar distributed
    complex orthogonal matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..M-1, 0..N-1]
        M, N-   matrix size

    OUTPUT PARAMETERS:
        A   -   Q*A, where Q is random MxM orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrndorthogonalfromtheleft(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n)
    {
        amp::campf<Precision> tau;
        amp::campf<Precision> lambda;
        int s;
        int i;
        int j;
        ap::template_1d_array< amp::campf<Precision> > w;
        ap::template_1d_array< amp::campf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;
        int i_;


        ap::ap_error::make_assertion(n>=1 && m>=1);
        if( m==1 )
        {
            
            //
            // special case
            //
            hqrnd::hqrndrandomize<Precision>(state);
            hqrnd::hqrndunit2<Precision>(state, tau.x, tau.y);
            for(j=0; j<=n-1; j++)
            {
                a(0,j) = a(0,j)*tau;
            }
            return;
        }
        
        //
        // General case.
        // First pass.
        //
        w.setbounds(0, n-1);
        v.setbounds(1, m);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=m; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                for(i=1; i<=s; i++)
                {
                    hqrnd::hqrndnormal2<Precision>(state, tau.x, tau.y);
                    v(i) = tau;
                }
                lambda = 0.0;
                for(i_=1; i_<=s;i_++)
                {
                    lambda += v(i_)*amp::conj(v(i_));
                }
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            creflections::complexgeneratereflection<Precision>(v, s, tau);
            v(1) = 1;
            creflections::complexapplyreflectionfromtheleft<Precision>(a, tau, v, m-s, m-1, 0, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=m-1; i++)
        {
            hqrnd::hqrndunit2<Precision>(state, tau.x, tau.y);
            for(i_=0; i_<=n-1;i_++)
            {
                a(i,i_) = tau*a(i,i_);
            }
        }
    }


    /*************************************************************************
    Symmetric multiplication of NxN matrix by random Haar distributed
    orthogonal  matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..N-1, 0..N-1]
        N   -   matrix size

    OUTPUT PARAMETERS:
        A   -   Q'*A*Q, where Q is random NxN orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void smatrixrndmultiply(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        amp::ampf<Precision> tau;
        amp::ampf<Precision> lambda;
        int s;
        int i;
        amp::ampf<Precision> u1;
        amp::ampf<Precision> u2;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;


        
        //
        // General case.
        //
        w.setbounds(0, n-1);
        v.setbounds(1, n);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=n; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                i = 1;
                while( i<=s )
                {
                    hqrnd::hqrndnormal2<Precision>(state, u1, u2);
                    v(i) = u1;
                    if( i+1<=s )
                    {
                        v(i+1) = u2;
                    }
                    i = i+2;
                }
                lambda = amp::vdotproduct(v.getvector(1, s), v.getvector(1, s));
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            reflections::generatereflection<Precision>(v, s, tau);
            v(1) = 1;
            reflections::applyreflectionfromtheright<Precision>(a, tau, v, 0, n-1, n-s, n-1, w);
            reflections::applyreflectionfromtheleft<Precision>(a, tau, v, n-s, n-1, 0, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=n-1; i++)
        {
            tau = 2*ap::randominteger(2)-1;
            amp::vmul(a.getcolumn(i, 0, n-1), tau);
            amp::vmul(a.getrow(i, 0, n-1), tau);
        }
    }


    /*************************************************************************
    Hermitian multiplication of NxN matrix by random Haar distributed
    complex orthogonal matrix

    INPUT PARAMETERS:
        A   -   matrix, array[0..N-1, 0..N-1]
        N   -   matrix size

    OUTPUT PARAMETERS:
        A   -   Q^H*A*Q, where Q is random NxN orthogonal matrix

      -- ALGLIB routine --
         04.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hmatrixrndmultiply(ap::template_2d_array< amp::campf<Precision> >& a,
        int n)
    {
        amp::campf<Precision> tau;
        amp::campf<Precision> lambda;
        int s;
        int i;
        ap::template_1d_array< amp::campf<Precision> > w;
        ap::template_1d_array< amp::campf<Precision> > v;
        hqrnd::hqrndstate<Precision> state;
        int i_;


        
        //
        // General case.
        //
        w.setbounds(0, n-1);
        v.setbounds(1, n);
        hqrnd::hqrndrandomize<Precision>(state);
        for(s=2; s<=n; s++)
        {
            
            //
            // Prepare random normal v
            //
            do
            {
                for(i=1; i<=s; i++)
                {
                    hqrnd::hqrndnormal2<Precision>(state, tau.x, tau.y);
                    v(i) = tau;
                }
                lambda = 0.0;
                for(i_=1; i_<=s;i_++)
                {
                    lambda += v(i_)*amp::conj(v(i_));
                }
            }
            while( lambda==0 );
            
            //
            // Prepare and apply reflection
            //
            creflections::complexgeneratereflection<Precision>(v, s, tau);
            v(1) = 1;
            creflections::complexapplyreflectionfromtheright<Precision>(a, tau, v, 0, n-1, n-s, n-1, w);
            creflections::complexapplyreflectionfromtheleft<Precision>(a, amp::conj<Precision>(tau), v, n-s, n-1, 0, n-1, w);
        }
        
        //
        // Second pass.
        //
        for(i=0; i<=n-1; i++)
        {
            hqrnd::hqrndunit2<Precision>(state, tau.x, tau.y);
            for(i_=0; i_<=n-1;i_++)
            {
                a(i_,i) = tau*a(i_,i);
            }
            tau = amp::conj<Precision>(tau);
            for(i_=0; i_<=n-1;i_++)
            {
                a(i,i_) = tau*a(i,i_);
            }
        }
    }
} // namespace

#endif
