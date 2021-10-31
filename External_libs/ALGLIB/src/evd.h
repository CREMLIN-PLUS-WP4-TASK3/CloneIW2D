/*************************************************************************
Copyright (c) 2005-2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _evd_h
#define _evd_h

#include "ap.h"
#include "amp.h"
#include "hblas.h"
#include "reflections.h"
#include "creflections.h"
#include "sblas.h"
#include "ablasf.h"
#include "ablas.h"
#include "ortfac.h"
#include "blas.h"
#include "rotations.h"
#include "hsschur.h"
namespace evd
{
    template<unsigned int Precision>
    bool smatrixevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixevdr(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixevdi(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        int i1,
        int i2,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool hmatrixevd(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::campf<Precision> >& z);
    template<unsigned int Precision>
    bool hmatrixevdr(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::campf<Precision> >& z);
    template<unsigned int Precision>
    bool hmatrixevdi(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        int i1,
        int i2,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::campf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixtdevd(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int zneeded,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixtdevdr(ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        int zneeded,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int& m,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixtdevdi(ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        int zneeded,
        int i1,
        int i2,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool rmatrixevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int vneeded,
        ap::template_1d_array< amp::ampf<Precision> >& wr,
        ap::template_1d_array< amp::ampf<Precision> >& wi,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr);
    template<unsigned int Precision>
    bool internalbisectioneigenvalues(ap::template_1d_array< amp::ampf<Precision> > d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int irange,
        int iorder,
        amp::ampf<Precision> vl,
        amp::ampf<Precision> vu,
        int il,
        int iu,
        amp::ampf<Precision> abstol,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        int& m,
        int& nsplit,
        ap::template_1d_array< int >& iblock,
        ap::template_1d_array< int >& isplit,
        int& errorcode);
    template<unsigned int Precision>
    void internaldstein(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        const int& m,
        ap::template_1d_array< amp::ampf<Precision> > w,
        const ap::template_1d_array< int >& iblock,
        const ap::template_1d_array< int >& isplit,
        ap::template_2d_array< amp::ampf<Precision> >& z,
        ap::template_1d_array< int >& ifail,
        int& info);
    template<unsigned int Precision>
    bool tridiagonalevd(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int zneeded,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    void tdevde2(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        amp::ampf<Precision>& rt1,
        amp::ampf<Precision>& rt2);
    template<unsigned int Precision>
    void tdevdev2(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        amp::ampf<Precision>& rt1,
        amp::ampf<Precision>& rt2,
        amp::ampf<Precision>& cs1,
        amp::ampf<Precision>& sn1);
    template<unsigned int Precision>
    amp::ampf<Precision> tdevdpythag(amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    amp::ampf<Precision> tdevdextsign(amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    void tdininternaldlagtf(const int& n,
        ap::template_1d_array< amp::ampf<Precision> >& a,
        const amp::ampf<Precision>& lambda,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        const amp::ampf<Precision>& tol,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< int >& iin,
        int& info);
    template<unsigned int Precision>
    void tdininternaldlagts(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< int >& iin,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        amp::ampf<Precision>& tol,
        int& info);
    template<unsigned int Precision>
    void internaldlaebz(const int& ijob,
        const int& nitmax,
        const int& n,
        const int& mmax,
        const int& minp,
        const amp::ampf<Precision>& abstol,
        const amp::ampf<Precision>& reltol,
        const amp::ampf<Precision>& pivmin,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        const ap::template_1d_array< amp::ampf<Precision> >& e2,
        ap::template_1d_array< int >& nval,
        ap::template_2d_array< amp::ampf<Precision> >& ab,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        int& mout,
        ap::template_2d_array< int >& nab,
        ap::template_1d_array< amp::ampf<Precision> >& work,
        ap::template_1d_array< int >& iwork,
        int& info);
    template<unsigned int Precision>
    void internaltrevc(const ap::template_2d_array< amp::ampf<Precision> >& t,
        int n,
        int side,
        int howmny,
        ap::template_1d_array< bool > vselect,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr,
        int& m,
        int& info);
    template<unsigned int Precision>
    void internalhsevdlaln2(const bool& ltrans,
        const int& na,
        const int& nw,
        const amp::ampf<Precision>& smin,
        const amp::ampf<Precision>& ca,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        const amp::ampf<Precision>& d1,
        const amp::ampf<Precision>& d2,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        const amp::ampf<Precision>& wr,
        const amp::ampf<Precision>& wi,
        ap::template_1d_array< bool >& rswap4,
        ap::template_1d_array< bool >& zswap4,
        ap::template_2d_array< int >& ipivot44,
        ap::template_1d_array< amp::ampf<Precision> >& civ4,
        ap::template_1d_array< amp::ampf<Precision> >& crv4,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& scl,
        amp::ampf<Precision>& xnorm,
        int& info);
    template<unsigned int Precision>
    void internalhsevdladiv(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        const amp::ampf<Precision>& d,
        amp::ampf<Precision>& p,
        amp::ampf<Precision>& q);
    template<unsigned int Precision>
    bool nonsymmetricevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int vneeded,
        ap::template_1d_array< amp::ampf<Precision> >& wr,
        ap::template_1d_array< amp::ampf<Precision> >& wi,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr);
    template<unsigned int Precision>
    void toupperhessenberg(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& tau);
    template<unsigned int Precision>
    void unpackqfromupperhessenberg(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& tau,
        ap::template_2d_array< amp::ampf<Precision> >& q);
    template<unsigned int Precision>
    void unpackhfromupperhessenberg(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& tau,
        ap::template_2d_array< amp::ampf<Precision> >& h);


    /*************************************************************************
    Finding the eigenvalues and eigenvectors of a symmetric matrix

    The algorithm finds eigen pairs of a symmetric matrix by reducing it to
    tridiagonal form and using the QL/QR algorithm.

    Input parameters:
        A       -   symmetric matrix which is given by its upper or lower
                    triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.

    Output parameters:
        D       -   eigenvalues in ascending order.
                    Array whose index ranges within [0..N-1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains the eigenvectors.
                    Array whose indexes range within [0..N-1, 0..N-1].
                    The eigenvectors are stored in the matrix columns.

    Result:
        True, if the algorithm has converged.
        False, if the algorithm hasn't converged (rare case).

      -- ALGLIB --
         Copyright 2005-2008 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        ortfac::smatrixtd<Precision>(a, n, isupper, tau, d, e);
        if( zneeded==1 )
        {
            ortfac::smatrixtdunpackq<Precision>(a, n, isupper, tau, z);
        }
        result = smatrixtdevd<Precision>(d, e, n, zneeded, z);
        return result;
    }


    /*************************************************************************
    Subroutine for finding the eigenvalues (and eigenvectors) of  a  symmetric
    matrix  in  a  given half open interval (A, B] by using  a  bisection  and
    inverse iteration

    Input parameters:
        A       -   symmetric matrix which is given by its upper or lower
                    triangular part. Array [0..N-1, 0..N-1].
        N       -   size of matrix A.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.
        IsUpperA -  storage format of matrix A.
        B1, B2 -    half open interval (B1, B2] to search eigenvalues in.

    Output parameters:
        M       -   number of eigenvalues found in a given half-interval (M>=0).
        W       -   array of the eigenvalues found.
                    Array whose index ranges within [0..M-1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains eigenvectors.
                    Array whose indexes range within [0..N-1, 0..M-1].
                    The eigenvectors are stored in the matrix columns.

    Result:
        True, if successful. M contains the number of eigenvalues in the given
        half-interval (could be equal to 0), W contains the eigenvalues,
        Z contains the eigenvectors (if needed).

        False, if the bisection method subroutine wasn't able to find the
        eigenvalues in the given interval or if the inverse iteration subroutine
        wasn't able to find all the corresponding eigenvectors.
        In that case, the eigenvalues and eigenvectors are not returned,
        M is equal to 0.

      -- ALGLIB --
         Copyright 07.01.2006 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixevdr(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        ortfac::smatrixtd<Precision>(a, n, isupper, tau, w, e);
        if( zneeded==1 )
        {
            ortfac::smatrixtdunpackq<Precision>(a, n, isupper, tau, z);
        }
        result = smatrixtdevdr<Precision>(w, e, n, zneeded, b1, b2, m, z);
        return result;
    }


    /*************************************************************************
    Subroutine for finding the eigenvalues and  eigenvectors  of  a  symmetric
    matrix with given indexes by using bisection and inverse iteration methods.

    Input parameters:
        A       -   symmetric matrix which is given by its upper or lower
                    triangular part. Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.
        IsUpperA -  storage format of matrix A.
        I1, I2 -    index interval for searching (from I1 to I2).
                    0 <= I1 <= I2 <= N-1.

    Output parameters:
        W       -   array of the eigenvalues found.
                    Array whose index ranges within [0..I2-I1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains eigenvectors.
                    Array whose indexes range within [0..N-1, 0..I2-I1].
                    In that case, the eigenvectors are stored in the matrix columns.

    Result:
        True, if successful. W contains the eigenvalues, Z contains the
        eigenvectors (if needed).

        False, if the bisection method subroutine wasn't able to find the
        eigenvalues in the given interval or if the inverse iteration subroutine
        wasn't able to find all the corresponding eigenvectors.
        In that case, the eigenvalues and eigenvectors are not returned.

      -- ALGLIB --
         Copyright 07.01.2006 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixevdi(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        int i1,
        int i2,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        ortfac::smatrixtd<Precision>(a, n, isupper, tau, w, e);
        if( zneeded==1 )
        {
            ortfac::smatrixtdunpackq<Precision>(a, n, isupper, tau, z);
        }
        result = smatrixtdevdi<Precision>(w, e, n, zneeded, i1, i2, z);
        return result;
    }


    /*************************************************************************
    Finding the eigenvalues and eigenvectors of a Hermitian matrix

    The algorithm finds eigen pairs of a Hermitian matrix by  reducing  it  to
    real tridiagonal form and using the QL/QR algorithm.

    Input parameters:
        A       -   Hermitian matrix which is given  by  its  upper  or  lower
                    triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.
        ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                    not. If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.

    Output parameters:
        D       -   eigenvalues in ascending order.
                    Array whose index ranges within [0..N-1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains the eigenvectors.
                    Array whose indexes range within [0..N-1, 0..N-1].
                    The eigenvectors are stored in the matrix columns.

    Result:
        True, if the algorithm has converged.
        False, if the algorithm hasn't converged (rare case).

    Note:
        eigenvectors of Hermitian matrix are defined up to  multiplication  by
        a complex number L, such that |L|=1.

      -- ALGLIB --
         Copyright 2005, 23 March 2007 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool hmatrixevd(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::campf<Precision> >& z)
    {
        bool result;
        ap::template_1d_array< amp::campf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_2d_array< amp::campf<Precision> > q;
        int i;
        int k;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        
        //
        // Reduce to tridiagonal form
        //
        ortfac::hmatrixtd<Precision>(a, n, isupper, tau, d, e);
        if( zneeded==1 )
        {
            ortfac::hmatrixtdunpackq<Precision>(a, n, isupper, tau, q);
            zneeded = 2;
        }
        
        //
        // TDEVD
        //
        result = smatrixtdevd<Precision>(d, e, n, zneeded, t);
        
        //
        // Eigenvectors are needed
        // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
        //
        if( result && zneeded!=0 )
        {
            work.setbounds(0, n-1);
            z.setbounds(0, n-1, 0, n-1);
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Calculate real part
                //
                for(k=0; k<=n-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).x;
                    amp::vadd(work.getvector(0, n-1), t.getrow(k, 0, n-1), v);
                }
                for(k=0; k<=n-1; k++)
                {
                    z(i,k).x = work(k);
                }
                
                //
                // Calculate imaginary part
                //
                for(k=0; k<=n-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).y;
                    amp::vadd(work.getvector(0, n-1), t.getrow(k, 0, n-1), v);
                }
                for(k=0; k<=n-1; k++)
                {
                    z(i,k).y = work(k);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Subroutine for finding the eigenvalues (and eigenvectors) of  a  Hermitian
    matrix  in  a  given half-interval (A, B] by using a bisection and inverse
    iteration

    Input parameters:
        A       -   Hermitian matrix which is given  by  its  upper  or  lower
                    triangular  part.  Array  whose   indexes   range   within
                    [0..N-1, 0..N-1].
        N       -   size of matrix A.
        ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                    not. If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.
        IsUpperA -  storage format of matrix A.
        B1, B2 -    half-interval (B1, B2] to search eigenvalues in.

    Output parameters:
        M       -   number of eigenvalues found in a given half-interval, M>=0
        W       -   array of the eigenvalues found.
                    Array whose index ranges within [0..M-1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains eigenvectors.
                    Array whose indexes range within [0..N-1, 0..M-1].
                    The eigenvectors are stored in the matrix columns.

    Result:
        True, if successful. M contains the number of eigenvalues in the given
        half-interval (could be equal to 0), W contains the eigenvalues,
        Z contains the eigenvectors (if needed).

        False, if the bisection method subroutine  wasn't  able  to  find  the
        eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
        subroutine  wasn't  able  to  find all the corresponding eigenvectors.
        In that case, the eigenvalues and eigenvectors are not returned, M  is
        equal to 0.

    Note:
        eigen vectors of Hermitian matrix are defined up to multiplication  by
        a complex number L, such as |L|=1.

      -- ALGLIB --
         Copyright 07.01.2006, 24.03.2007 by Bochkanov Sergey.
    *************************************************************************/
    template<unsigned int Precision>
    bool hmatrixevdr(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        amp::ampf<Precision> b1,
        amp::ampf<Precision> b2,
        int& m,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::campf<Precision> >& z)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::campf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_1d_array< amp::ampf<Precision> > work;
        int i;
        int k;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        
        //
        // Reduce to tridiagonal form
        //
        ortfac::hmatrixtd<Precision>(a, n, isupper, tau, w, e);
        if( zneeded==1 )
        {
            ortfac::hmatrixtdunpackq<Precision>(a, n, isupper, tau, q);
            zneeded = 2;
        }
        
        //
        // Bisection and inverse iteration
        //
        result = smatrixtdevdr<Precision>(w, e, n, zneeded, b1, b2, m, t);
        
        //
        // Eigenvectors are needed
        // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
        //
        if( result && zneeded!=0 && m!=0 )
        {
            work.setbounds(0, m-1);
            z.setbounds(0, n-1, 0, m-1);
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Calculate real part
                //
                for(k=0; k<=m-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).x;
                    amp::vadd(work.getvector(0, m-1), t.getrow(k, 0, m-1), v);
                }
                for(k=0; k<=m-1; k++)
                {
                    z(i,k).x = work(k);
                }
                
                //
                // Calculate imaginary part
                //
                for(k=0; k<=m-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).y;
                    amp::vadd(work.getvector(0, m-1), t.getrow(k, 0, m-1), v);
                }
                for(k=0; k<=m-1; k++)
                {
                    z(i,k).y = work(k);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Subroutine for finding the eigenvalues and  eigenvectors  of  a  Hermitian
    matrix with given indexes by using bisection and inverse iteration methods

    Input parameters:
        A       -   Hermitian matrix which is given  by  its  upper  or  lower
                    triangular part.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        ZNeeded -   flag controlling whether the eigenvectors  are  needed  or
                    not. If ZNeeded is equal to:
                     * 0, the eigenvectors are not returned;
                     * 1, the eigenvectors are returned.
        IsUpperA -  storage format of matrix A.
        I1, I2 -    index interval for searching (from I1 to I2).
                    0 <= I1 <= I2 <= N-1.

    Output parameters:
        W       -   array of the eigenvalues found.
                    Array whose index ranges within [0..I2-I1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains eigenvectors.
                    Array whose indexes range within [0..N-1, 0..I2-I1].
                    In  that  case,  the eigenvectors are stored in the matrix
                    columns.

    Result:
        True, if successful. W contains the eigenvalues, Z contains the
        eigenvectors (if needed).

        False, if the bisection method subroutine  wasn't  able  to  find  the
        eigenvalues  in  the  given  interval  or  if  the  inverse  iteration
        subroutine wasn't able to find  all  the  corresponding  eigenvectors.
        In that case, the eigenvalues and eigenvectors are not returned.

    Note:
        eigen vectors of Hermitian matrix are defined up to multiplication  by
        a complex number L, such as |L|=1.

      -- ALGLIB --
         Copyright 07.01.2006, 24.03.2007 by Bochkanov Sergey.
    *************************************************************************/
    template<unsigned int Precision>
    bool hmatrixevdi(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        int zneeded,
        bool isupper,
        int i1,
        int i2,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::campf<Precision> >& z)
    {
        bool result;
        ap::template_2d_array< amp::campf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::campf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_1d_array< amp::ampf<Precision> > work;
        int i;
        int k;
        amp::ampf<Precision> v;
        int m;


        ap::ap_error::make_assertion(zneeded==0 || zneeded==1);
        
        //
        // Reduce to tridiagonal form
        //
        ortfac::hmatrixtd<Precision>(a, n, isupper, tau, w, e);
        if( zneeded==1 )
        {
            ortfac::hmatrixtdunpackq<Precision>(a, n, isupper, tau, q);
            zneeded = 2;
        }
        
        //
        // Bisection and inverse iteration
        //
        result = smatrixtdevdi<Precision>(w, e, n, zneeded, i1, i2, t);
        
        //
        // Eigenvectors are needed
        // Calculate Z = Q*T = Re(Q)*T + i*Im(Q)*T
        //
        m = i2-i1+1;
        if( result && zneeded!=0 )
        {
            work.setbounds(0, m-1);
            z.setbounds(0, n-1, 0, m-1);
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Calculate real part
                //
                for(k=0; k<=m-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).x;
                    amp::vadd(work.getvector(0, m-1), t.getrow(k, 0, m-1), v);
                }
                for(k=0; k<=m-1; k++)
                {
                    z(i,k).x = work(k);
                }
                
                //
                // Calculate imaginary part
                //
                for(k=0; k<=m-1; k++)
                {
                    work(k) = 0;
                }
                for(k=0; k<=n-1; k++)
                {
                    v = q(i,k).y;
                    amp::vadd(work.getvector(0, m-1), t.getrow(k, 0, m-1), v);
                }
                for(k=0; k<=m-1; k++)
                {
                    z(i,k).y = work(k);
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Finding the eigenvalues and eigenvectors of a tridiagonal symmetric matrix

    The algorithm finds the eigen pairs of a tridiagonal symmetric matrix by
    using an QL/QR algorithm with implicit shifts.

    Input parameters:
        D       -   the main diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-1].
        E       -   the secondary diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-2].
        N       -   size of matrix A.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not needed;
                     * 1, the eigenvectors of a tridiagonal matrix
                       are multiplied by the square matrix Z. It is used if the
                       tridiagonal matrix is obtained by the similarity
                       transformation of a symmetric matrix;
                     * 2, the eigenvectors of a tridiagonal matrix replace the
                       square matrix Z;
                     * 3, matrix Z contains the first row of the eigenvectors
                       matrix.
        Z       -   if ZNeeded=1, Z contains the square matrix by which the
                    eigenvectors are multiplied.
                    Array whose indexes range within [0..N-1, 0..N-1].

    Output parameters:
        D       -   eigenvalues in ascending order.
                    Array whose index ranges within [0..N-1].
        Z       -   if ZNeeded is equal to:
                     * 0, Z hasn’t changed;
                     * 1, Z contains the product of a given matrix (from the left)
                       and the eigenvectors matrix (from the right);
                     * 2, Z contains the eigenvectors.
                     * 3, Z contains the first row of the eigenvectors matrix.
                    If ZNeeded<3, Z is the array whose indexes range within [0..N-1, 0..N-1].
                    In that case, the eigenvectors are stored in the matrix columns.
                    If ZNeeded=3, Z is the array whose indexes range within [0..0, 0..N-1].

    Result:
        True, if the algorithm has converged.
        False, if the algorithm hasn't converged.

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         September 30, 1994
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixtdevd(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int zneeded,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > d1;
        ap::template_1d_array< amp::ampf<Precision> > e1;
        ap::template_2d_array< amp::ampf<Precision> > z1;
        int i;


        
        //
        // Prepare 1-based task
        //
        d1.setbounds(1, n);
        e1.setbounds(1, n);
        amp::vmove(d1.getvector(1, n), d.getvector(0, n-1));
        if( n>1 )
        {
            amp::vmove(e1.getvector(1, n-1), e.getvector(0, n-2));
        }
        if( zneeded==1 )
        {
            z1.setbounds(1, n, 1, n);
            for(i=1; i<=n; i++)
            {
                amp::vmove(z1.getrow(i, 1, n), z.getrow(i-1, 0, n-1));
            }
        }
        
        //
        // Solve 1-based task
        //
        result = tridiagonalevd<Precision>(d1, e1, n, zneeded, z1);
        if( !result )
        {
            return result;
        }
        
        //
        // Convert back to 0-based result
        //
        amp::vmove(d.getvector(0, n-1), d1.getvector(1, n));
        if( zneeded!=0 )
        {
            if( zneeded==1 )
            {
                for(i=1; i<=n; i++)
                {
                    amp::vmove(z.getrow(i-1, 0, n-1), z1.getrow(i, 1, n));
                }
                return result;
            }
            if( zneeded==2 )
            {
                z.setbounds(0, n-1, 0, n-1);
                for(i=1; i<=n; i++)
                {
                    amp::vmove(z.getrow(i-1, 0, n-1), z1.getrow(i, 1, n));
                }
                return result;
            }
            if( zneeded==3 )
            {
                z.setbounds(0, 0, 0, n-1);
                amp::vmove(z.getrow(0, 0, n-1), z1.getrow(1, 1, n));
                return result;
            }
            ap::ap_error::make_assertion(false);
        }
        return result;
    }


    /*************************************************************************
    Subroutine for finding the tridiagonal matrix eigenvalues/vectors in a
    given half-interval (A, B] by using bisection and inverse iteration.

    Input parameters:
        D       -   the main diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-1].
        E       -   the secondary diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-2].
        N       -   size of matrix, N>=0.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not needed;
                     * 1, the eigenvectors of a tridiagonal matrix are multiplied
                       by the square matrix Z. It is used if the tridiagonal
                       matrix is obtained by the similarity transformation
                       of a symmetric matrix.
                     * 2, the eigenvectors of a tridiagonal matrix replace matrix Z.
        A, B    -   half-interval (A, B] to search eigenvalues in.
        Z       -   if ZNeeded is equal to:
                     * 0, Z isn't used and remains unchanged;
                     * 1, Z contains the square matrix (array whose indexes range
                       within [0..N-1, 0..N-1]) which reduces the given symmetric
                       matrix to tridiagonal form;
                     * 2, Z isn't used (but changed on the exit).

    Output parameters:
        D       -   array of the eigenvalues found.
                    Array whose index ranges within [0..M-1].
        M       -   number of eigenvalues found in the given half-interval (M>=0).
        Z       -   if ZNeeded is equal to:
                     * 0, doesn't contain any information;
                     * 1, contains the product of a given NxN matrix Z (from the
                       left) and NxM matrix of the eigenvectors found (from the
                       right). Array whose indexes range within [0..N-1, 0..M-1].
                     * 2, contains the matrix of the eigenvectors found.
                       Array whose indexes range within [0..N-1, 0..M-1].

    Result:

        True, if successful. In that case, M contains the number of eigenvalues
        in the given half-interval (could be equal to 0), D contains the eigenvalues,
        Z contains the eigenvectors (if needed).
        It should be noted that the subroutine changes the size of arrays D and Z.

        False, if the bisection method subroutine wasn't able to find the
        eigenvalues in the given interval or if the inverse iteration subroutine
        wasn't able to find all the corresponding eigenvectors. In that case,
        the eigenvalues and eigenvectors are not returned, M is equal to 0.

      -- ALGLIB --
         Copyright 31.03.2008 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixtdevdr(ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        int zneeded,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b,
        int& m,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        int errorcode;
        int nsplit;
        int i;
        int j;
        int k;
        int cr;
        ap::template_1d_array< int > iblock;
        ap::template_1d_array< int > isplit;
        ap::template_1d_array< int > ifail;
        ap::template_1d_array< amp::ampf<Precision> > d1;
        ap::template_1d_array< amp::ampf<Precision> > e1;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_2d_array< amp::ampf<Precision> > z2;
        ap::template_2d_array< amp::ampf<Precision> > z3;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(zneeded>=0 && zneeded<=2);
        
        //
        // Special cases
        //
        if( b<=a )
        {
            m = 0;
            result = true;
            return result;
        }
        if( n<=0 )
        {
            m = 0;
            result = true;
            return result;
        }
        
        //
        // Copy D,E to D1, E1
        //
        d1.setbounds(1, n);
        amp::vmove(d1.getvector(1, n), d.getvector(0, n-1));
        if( n>1 )
        {
            e1.setbounds(1, n-1);
            amp::vmove(e1.getvector(1, n-1), e.getvector(0, n-2));
        }
        
        //
        // No eigen vectors
        //
        if( zneeded==0 )
        {
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 2, 1, a, b, 0, 0, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result || m==0 )
            {
                m = 0;
                return result;
            }
            d.setbounds(0, m-1);
            amp::vmove(d.getvector(0, m-1), w.getvector(1, m));
            return result;
        }
        
        //
        // Eigen vectors are multiplied by Z
        //
        if( zneeded==1 )
        {
            
            //
            // Find eigen pairs
            //
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 2, 2, a, b, 0, 0, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result || m==0 )
            {
                m = 0;
                return result;
            }
            internaldstein<Precision>(n, d1, e1, m, w, iblock, isplit, z2, ifail, cr);
            if( cr!=0 )
            {
                m = 0;
                result = false;
                return result;
            }
            
            //
            // Sort eigen values and vectors
            //
            for(i=1; i<=m; i++)
            {
                k = i;
                for(j=i; j<=m; j++)
                {
                    if( w(j)<w(k) )
                    {
                        k = j;
                    }
                }
                v = w(i);
                w(i) = w(k);
                w(k) = v;
                for(j=1; j<=n; j++)
                {
                    v = z2(j,i);
                    z2(j,i) = z2(j,k);
                    z2(j,k) = v;
                }
            }
            
            //
            // Transform Z2 and overwrite Z
            //
            z3.setbounds(1, m, 1, n);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z3.getrow(i, 1, n), z2.getcolumn(i, 1, n));
            }
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=m; j++)
                {
                    v = amp::vdotproduct(z.getrow(i-1, 0, n-1), z3.getrow(j, 1, n));
                    z2(i,j) = v;
                }
            }
            z.setbounds(0, n-1, 0, m-1);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z.getcolumn(i-1, 0, n-1), z2.getcolumn(i, 1, n));
            }
            
            //
            // Store W
            //
            d.setbounds(0, m-1);
            for(i=1; i<=m; i++)
            {
                d(i-1) = w(i);
            }
            return result;
        }
        
        //
        // Eigen vectors are stored in Z
        //
        if( zneeded==2 )
        {
            
            //
            // Find eigen pairs
            //
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 2, 2, a, b, 0, 0, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result || m==0 )
            {
                m = 0;
                return result;
            }
            internaldstein<Precision>(n, d1, e1, m, w, iblock, isplit, z2, ifail, cr);
            if( cr!=0 )
            {
                m = 0;
                result = false;
                return result;
            }
            
            //
            // Sort eigen values and vectors
            //
            for(i=1; i<=m; i++)
            {
                k = i;
                for(j=i; j<=m; j++)
                {
                    if( w(j)<w(k) )
                    {
                        k = j;
                    }
                }
                v = w(i);
                w(i) = w(k);
                w(k) = v;
                for(j=1; j<=n; j++)
                {
                    v = z2(j,i);
                    z2(j,i) = z2(j,k);
                    z2(j,k) = v;
                }
            }
            
            //
            // Store W
            //
            d.setbounds(0, m-1);
            for(i=1; i<=m; i++)
            {
                d(i-1) = w(i);
            }
            z.setbounds(0, n-1, 0, m-1);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z.getcolumn(i-1, 0, n-1), z2.getcolumn(i, 1, n));
            }
            return result;
        }
        result = false;
        return result;
    }


    /*************************************************************************
    Subroutine for finding tridiagonal matrix eigenvalues/vectors with given
    indexes (in ascending order) by using the bisection and inverse iteraion.

    Input parameters:
        D       -   the main diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-1].
        E       -   the secondary diagonal of a tridiagonal matrix.
                    Array whose index ranges within [0..N-2].
        N       -   size of matrix. N>=0.
        ZNeeded -   flag controlling whether the eigenvectors are needed or not.
                    If ZNeeded is equal to:
                     * 0, the eigenvectors are not needed;
                     * 1, the eigenvectors of a tridiagonal matrix are multiplied
                       by the square matrix Z. It is used if the
                       tridiagonal matrix is obtained by the similarity transformation
                       of a symmetric matrix.
                     * 2, the eigenvectors of a tridiagonal matrix replace
                       matrix Z.
        I1, I2  -   index interval for searching (from I1 to I2).
                    0 <= I1 <= I2 <= N-1.
        Z       -   if ZNeeded is equal to:
                     * 0, Z isn't used and remains unchanged;
                     * 1, Z contains the square matrix (array whose indexes range within [0..N-1, 0..N-1])
                       which reduces the given symmetric matrix to  tridiagonal form;
                     * 2, Z isn't used (but changed on the exit).

    Output parameters:
        D       -   array of the eigenvalues found.
                    Array whose index ranges within [0..I2-I1].
        Z       -   if ZNeeded is equal to:
                     * 0, doesn't contain any information;
                     * 1, contains the product of a given NxN matrix Z (from the left) and
                       Nx(I2-I1) matrix of the eigenvectors found (from the right).
                       Array whose indexes range within [0..N-1, 0..I2-I1].
                     * 2, contains the matrix of the eigenvalues found.
                       Array whose indexes range within [0..N-1, 0..I2-I1].


    Result:

        True, if successful. In that case, D contains the eigenvalues,
        Z contains the eigenvectors (if needed).
        It should be noted that the subroutine changes the size of arrays D and Z.

        False, if the bisection method subroutine wasn't able to find the eigenvalues
        in the given interval or if the inverse iteration subroutine wasn't able
        to find all the corresponding eigenvectors. In that case, the eigenvalues
        and eigenvectors are not returned.

      -- ALGLIB --
         Copyright 25.12.2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixtdevdi(ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        int n,
        int zneeded,
        int i1,
        int i2,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        int errorcode;
        int nsplit;
        int i;
        int j;
        int k;
        int m;
        int cr;
        ap::template_1d_array< int > iblock;
        ap::template_1d_array< int > isplit;
        ap::template_1d_array< int > ifail;
        ap::template_1d_array< amp::ampf<Precision> > w;
        ap::template_1d_array< amp::ampf<Precision> > d1;
        ap::template_1d_array< amp::ampf<Precision> > e1;
        ap::template_2d_array< amp::ampf<Precision> > z2;
        ap::template_2d_array< amp::ampf<Precision> > z3;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(0<=i1 && i1<=i2 && i2<n);
        
        //
        // Copy D,E to D1, E1
        //
        d1.setbounds(1, n);
        amp::vmove(d1.getvector(1, n), d.getvector(0, n-1));
        if( n>1 )
        {
            e1.setbounds(1, n-1);
            amp::vmove(e1.getvector(1, n-1), e.getvector(0, n-2));
        }
        
        //
        // No eigen vectors
        //
        if( zneeded==0 )
        {
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 3, 1, amp::ampf<Precision>(0), amp::ampf<Precision>(0), i1+1, i2+1, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result )
            {
                return result;
            }
            if( m!=i2-i1+1 )
            {
                result = false;
                return result;
            }
            d.setbounds(0, m-1);
            for(i=1; i<=m; i++)
            {
                d(i-1) = w(i);
            }
            return result;
        }
        
        //
        // Eigen vectors are multiplied by Z
        //
        if( zneeded==1 )
        {
            
            //
            // Find eigen pairs
            //
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 3, 2, amp::ampf<Precision>(0), amp::ampf<Precision>(0), i1+1, i2+1, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result )
            {
                return result;
            }
            if( m!=i2-i1+1 )
            {
                result = false;
                return result;
            }
            internaldstein<Precision>(n, d1, e1, m, w, iblock, isplit, z2, ifail, cr);
            if( cr!=0 )
            {
                result = false;
                return result;
            }
            
            //
            // Sort eigen values and vectors
            //
            for(i=1; i<=m; i++)
            {
                k = i;
                for(j=i; j<=m; j++)
                {
                    if( w(j)<w(k) )
                    {
                        k = j;
                    }
                }
                v = w(i);
                w(i) = w(k);
                w(k) = v;
                for(j=1; j<=n; j++)
                {
                    v = z2(j,i);
                    z2(j,i) = z2(j,k);
                    z2(j,k) = v;
                }
            }
            
            //
            // Transform Z2 and overwrite Z
            //
            z3.setbounds(1, m, 1, n);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z3.getrow(i, 1, n), z2.getcolumn(i, 1, n));
            }
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=m; j++)
                {
                    v = amp::vdotproduct(z.getrow(i-1, 0, n-1), z3.getrow(j, 1, n));
                    z2(i,j) = v;
                }
            }
            z.setbounds(0, n-1, 0, m-1);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z.getcolumn(i-1, 0, n-1), z2.getcolumn(i, 1, n));
            }
            
            //
            // Store W
            //
            d.setbounds(0, m-1);
            for(i=1; i<=m; i++)
            {
                d(i-1) = w(i);
            }
            return result;
        }
        
        //
        // Eigen vectors are stored in Z
        //
        if( zneeded==2 )
        {
            
            //
            // Find eigen pairs
            //
            result = internalbisectioneigenvalues<Precision>(d1, e1, n, 3, 2, amp::ampf<Precision>(0), amp::ampf<Precision>(0), i1+1, i2+1, amp::ampf<Precision>(-1), w, m, nsplit, iblock, isplit, errorcode);
            if( !result )
            {
                return result;
            }
            if( m!=i2-i1+1 )
            {
                result = false;
                return result;
            }
            internaldstein<Precision>(n, d1, e1, m, w, iblock, isplit, z2, ifail, cr);
            if( cr!=0 )
            {
                result = false;
                return result;
            }
            
            //
            // Sort eigen values and vectors
            //
            for(i=1; i<=m; i++)
            {
                k = i;
                for(j=i; j<=m; j++)
                {
                    if( w(j)<w(k) )
                    {
                        k = j;
                    }
                }
                v = w(i);
                w(i) = w(k);
                w(k) = v;
                for(j=1; j<=n; j++)
                {
                    v = z2(j,i);
                    z2(j,i) = z2(j,k);
                    z2(j,k) = v;
                }
            }
            
            //
            // Store Z
            //
            z.setbounds(0, n-1, 0, m-1);
            for(i=1; i<=m; i++)
            {
                amp::vmove(z.getcolumn(i-1, 0, n-1), z2.getcolumn(i, 1, n));
            }
            
            //
            // Store W
            //
            d.setbounds(0, m-1);
            for(i=1; i<=m; i++)
            {
                d(i-1) = w(i);
            }
            return result;
        }
        result = false;
        return result;
    }


    /*************************************************************************
    Finding eigenvalues and eigenvectors of a general matrix

    The algorithm finds eigenvalues and eigenvectors of a general matrix by
    using the QR algorithm with multiple shifts. The algorithm can find
    eigenvalues and both left and right eigenvectors.

    The right eigenvector is a vector x such that A*x = w*x, and the left
    eigenvector is a vector y such that y'*A = w*y' (here y' implies a complex
    conjugate transposition of vector y).

    Input parameters:
        A       -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        VNeeded -   flag controlling whether eigenvectors are needed or not.
                    If VNeeded is equal to:
                     * 0, eigenvectors are not returned;
                     * 1, right eigenvectors are returned;
                     * 2, left eigenvectors are returned;
                     * 3, both left and right eigenvectors are returned.

    Output parameters:
        WR      -   real parts of eigenvalues.
                    Array whose index ranges within [0..N-1].
        WR      -   imaginary parts of eigenvalues.
                    Array whose index ranges within [0..N-1].
        VL, VR  -   arrays of left and right eigenvectors (if they are needed).
                    If WI[i]=0, the respective eigenvalue is a real number,
                    and it corresponds to the column number I of matrices VL/VR.
                    If WI[i]>0, we have a pair of complex conjugate numbers with
                    positive and negative imaginary parts:
                        the first eigenvalue WR[i] + sqrt(-1)*WI[i];
                        the second eigenvalue WR[i+1] + sqrt(-1)*WI[i+1];
                        WI[i]>0
                        WI[i+1] = -WI[i] < 0
                    In that case, the eigenvector  corresponding to the first
                    eigenvalue is located in i and i+1 columns of matrices
                    VL/VR (the column number i contains the real part, and the
                    column number i+1 contains the imaginary part), and the vector
                    corresponding to the second eigenvalue is a complex conjugate to
                    the first vector.
                    Arrays whose indexes range within [0..N-1, 0..N-1].

    Result:
        True, if the algorithm has converged.
        False, if the algorithm has not converged.

    Note 1:
        Some users may ask the following question: what if WI[N-1]>0?
        WI[N] must contain an eigenvalue which is complex conjugate to the
        N-th eigenvalue, but the array has only size N?
        The answer is as follows: such a situation cannot occur because the
        algorithm finds a pairs of eigenvalues, therefore, if WI[i]>0, I is
        strictly less than N-1.

    Note 2:
        The algorithm performance depends on the value of the internal parameter
        NS of the InternalSchurDecomposition subroutine which defines the number
        of shifts in the QR algorithm (similarly to the block width in block-matrix
        algorithms of linear algebra). If you require maximum performance
        on your machine, it is recommended to adjust this parameter manually.


    See also the InternalTREVC subroutine.

    The algorithm is based on the LAPACK 3.0 library.
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int vneeded,
        ap::template_1d_array< amp::ampf<Precision> >& wr,
        ap::template_1d_array< amp::ampf<Precision> >& wi,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_2d_array< amp::ampf<Precision> > vl1;
        ap::template_2d_array< amp::ampf<Precision> > vr1;
        ap::template_1d_array< amp::ampf<Precision> > wr1;
        ap::template_1d_array< amp::ampf<Precision> > wi1;
        int i;


        ap::ap_error::make_assertion(vneeded>=0 && vneeded<=3);
        a1.setbounds(1, n, 1, n);
        for(i=1; i<=n; i++)
        {
            amp::vmove(a1.getrow(i, 1, n), a.getrow(i-1, 0, n-1));
        }
        result = nonsymmetricevd<Precision>(a1, n, vneeded, wr1, wi1, vl1, vr1);
        if( result )
        {
            wr.setbounds(0, n-1);
            wi.setbounds(0, n-1);
            amp::vmove(wr.getvector(0, n-1), wr1.getvector(1, n));
            amp::vmove(wi.getvector(0, n-1), wi1.getvector(1, n));
            if( vneeded==2 || vneeded==3 )
            {
                vl.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(vl.getrow(i, 0, n-1), vl1.getrow(i+1, 1, n));
                }
            }
            if( vneeded==1 || vneeded==3 )
            {
                vr.setbounds(0, n-1, 0, n-1);
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(vr.getrow(i, 0, n-1), vr1.getrow(i+1, 1, n));
                }
            }
        }
        return result;
    }


    template<unsigned int Precision>
    bool internalbisectioneigenvalues(ap::template_1d_array< amp::ampf<Precision> > d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int irange,
        int iorder,
        amp::ampf<Precision> vl,
        amp::ampf<Precision> vu,
        int il,
        int iu,
        amp::ampf<Precision> abstol,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        int& m,
        int& nsplit,
        ap::template_1d_array< int >& iblock,
        ap::template_1d_array< int >& isplit,
        int& errorcode)
    {
        bool result;
        amp::ampf<Precision> fudge;
        amp::ampf<Precision> relfac;
        bool ncnvrg;
        bool toofew;
        int ib;
        int ibegin;
        int idiscl;
        int idiscu;
        int ie;
        int iend;
        int iinfo;
        int im;
        int iin;
        int ioff;
        int iout;
        int itmax;
        int iw;
        int iwoff;
        int j;
        int itmp1;
        int jb;
        int jdisc;
        int je;
        int nwl;
        int nwu;
        amp::ampf<Precision> atoli;
        amp::ampf<Precision> bnorm;
        amp::ampf<Precision> gl;
        amp::ampf<Precision> gu;
        amp::ampf<Precision> pivmin;
        amp::ampf<Precision> rtoli;
        amp::ampf<Precision> safemn;
        amp::ampf<Precision> tmp1;
        amp::ampf<Precision> tmp2;
        amp::ampf<Precision> tnorm;
        amp::ampf<Precision> ulp;
        amp::ampf<Precision> wkill;
        amp::ampf<Precision> wl;
        amp::ampf<Precision> wlu;
        amp::ampf<Precision> wu;
        amp::ampf<Precision> wul;
        amp::ampf<Precision> scalefactor;
        amp::ampf<Precision> t;
        ap::template_1d_array< int > idumma;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_1d_array< int > iwork;
        ap::template_1d_array< int > ia1s2;
        ap::template_1d_array< amp::ampf<Precision> > ra1s2;
        ap::template_2d_array< amp::ampf<Precision> > ra1s2x2;
        ap::template_2d_array< int > ia1s2x2;
        ap::template_1d_array< amp::ampf<Precision> > ra1siin;
        ap::template_1d_array< amp::ampf<Precision> > ra2siin;
        ap::template_1d_array< amp::ampf<Precision> > ra3siin;
        ap::template_1d_array< amp::ampf<Precision> > ra4siin;
        ap::template_2d_array< amp::ampf<Precision> > ra1siinx2;
        ap::template_2d_array< int > ia1siinx2;
        ap::template_1d_array< int > iworkspace;
        ap::template_1d_array< amp::ampf<Precision> > rworkspace;
        int tmpi;


        
        //
        // Quick return if possible
        //
        m = 0;
        if( n==0 )
        {
            result = true;
            return result;
        }
        
        //
        // Get machine constants
        // NB is the minimum vector length for vector bisection, or 0
        // if only scalar is to be done.
        //
        fudge = 2;
        relfac = 2;
        safemn = amp::ampf<Precision>::getAlgoPascalMinNumber();
        ulp = 2*amp::ampf<Precision>::getAlgoPascalEpsilon();
        rtoli = ulp*relfac;
        idumma.setbounds(1, 1);
        work.setbounds(1, 4*n);
        iwork.setbounds(1, 3*n);
        w.setbounds(1, n);
        iblock.setbounds(1, n);
        isplit.setbounds(1, n);
        ia1s2.setbounds(1, 2);
        ra1s2.setbounds(1, 2);
        ra1s2x2.setbounds(1, 2, 1, 2);
        ia1s2x2.setbounds(1, 2, 1, 2);
        ra1siin.setbounds(1, n);
        ra2siin.setbounds(1, n);
        ra3siin.setbounds(1, n);
        ra4siin.setbounds(1, n);
        ra1siinx2.setbounds(1, n, 1, 2);
        ia1siinx2.setbounds(1, n, 1, 2);
        iworkspace.setbounds(1, n);
        rworkspace.setbounds(1, n);
        
        //
        // Check for Errors
        //
        result = false;
        errorcode = 0;
        if( irange<=0 || irange>=4 )
        {
            errorcode = -4;
        }
        if( iorder<=0 || iorder>=3 )
        {
            errorcode = -5;
        }
        if( n<0 )
        {
            errorcode = -3;
        }
        if( irange==2 && vl>=vu )
        {
            errorcode = -6;
        }
        if( irange==3 && (il<1 || il>ap::maxint(1, n)) )
        {
            errorcode = -8;
        }
        if( irange==3 && (iu<ap::minint(n, il) || iu>n) )
        {
            errorcode = -9;
        }
        if( errorcode!=0 )
        {
            return result;
        }
        
        //
        // Initialize error flags
        //
        ncnvrg = false;
        toofew = false;
        
        //
        // Simplifications:
        //
        if( irange==3 && il==1 && iu==n )
        {
            irange = 1;
        }
        
        //
        // Special Case when N=1
        //
        if( n==1 )
        {
            nsplit = 1;
            isplit(1) = 1;
            if( irange==2 && (vl>=d(1) || vu<d(1)) )
            {
                m = 0;
            }
            else
            {
                w(1) = d(1);
                iblock(1) = 1;
                m = 1;
            }
            result = true;
            return result;
        }
        
        //
        // Scaling
        //
        t = amp::abs<Precision>(d(n));
        for(j=1; j<=n-1; j++)
        {
            t = amp::maximum<Precision>(t, amp::abs<Precision>(d(j)));
            t = amp::maximum<Precision>(t, amp::abs<Precision>(e(j)));
        }
        scalefactor = 1;
        if( t!=0 )
        {
            if( t>amp::sqrt<Precision>(amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber()))*amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()) )
            {
                scalefactor = t;
            }
            if( t<amp::sqrt<Precision>(amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber()))*amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber()) )
            {
                scalefactor = t;
            }
            for(j=1; j<=n-1; j++)
            {
                d(j) = d(j)/scalefactor;
                e(j) = e(j)/scalefactor;
            }
            d(n) = d(n)/scalefactor;
        }
        
        //
        // Compute Splitting Points
        //
        nsplit = 1;
        work(n) = 0;
        pivmin = 1;
        for(j=2; j<=n; j++)
        {
            tmp1 = amp::sqr<Precision>(e(j-1));
            if( amp::abs<Precision>(d(j)*d(j-1))*amp::sqr<Precision>(ulp)+safemn>tmp1 )
            {
                isplit(nsplit) = j-1;
                nsplit = nsplit+1;
                work(j-1) = 0;
            }
            else
            {
                work(j-1) = tmp1;
                pivmin = amp::maximum<Precision>(pivmin, tmp1);
            }
        }
        isplit(nsplit) = n;
        pivmin = pivmin*safemn;
        
        //
        // Compute Interval and ATOLI
        //
        if( irange==3 )
        {
            
            //
            // RANGE='I': Compute the interval containing eigenvalues
            //     IL through IU.
            //
            // Compute Gershgorin interval for entire (split) matrix
            // and use it as the initial interval
            //
            gu = d(1);
            gl = d(1);
            tmp1 = 0;
            for(j=1; j<=n-1; j++)
            {
                tmp2 = amp::sqrt<Precision>(work(j));
                gu = amp::maximum<Precision>(gu, d(j)+tmp1+tmp2);
                gl = amp::minimum<Precision>(gl, d(j)-tmp1-tmp2);
                tmp1 = tmp2;
            }
            gu = amp::maximum<Precision>(gu, d(n)+tmp1);
            gl = amp::minimum<Precision>(gl, d(n)-tmp1);
            tnorm = amp::maximum<Precision>(amp::abs<Precision>(gl), amp::abs<Precision>(gu));
            gl = gl-fudge*tnorm*ulp*n-fudge*2*pivmin;
            gu = gu+fudge*tnorm*ulp*n+fudge*pivmin;
            
            //
            // Compute Iteration parameters
            //
            itmax = amp::ceil<Precision>((amp::log<Precision>(tnorm+pivmin)-amp::log<Precision>(pivmin))/amp::log<Precision>(amp::ampf<Precision>(2)))+2;
            if( abstol<=0 )
            {
                atoli = ulp*tnorm;
            }
            else
            {
                atoli = abstol;
            }
            work(n+1) = gl;
            work(n+2) = gl;
            work(n+3) = gu;
            work(n+4) = gu;
            work(n+5) = gl;
            work(n+6) = gu;
            iwork(1) = -1;
            iwork(2) = -1;
            iwork(3) = n+1;
            iwork(4) = n+1;
            iwork(5) = il-1;
            iwork(6) = iu;
            
            //
            // Calling DLAEBZ
            //
            // DLAEBZ( 3, ITMAX, N, 2, 2, NB, ATOLI, RTOLI, PIVMIN, D, E,
            //    WORK, IWORK( 5 ), WORK( N+1 ), WORK( N+5 ), IOUT,
            //    IWORK, W, IBLOCK, IINFO )
            //
            ia1s2(1) = iwork(5);
            ia1s2(2) = iwork(6);
            ra1s2(1) = work(n+5);
            ra1s2(2) = work(n+6);
            ra1s2x2(1,1) = work(n+1);
            ra1s2x2(2,1) = work(n+2);
            ra1s2x2(1,2) = work(n+3);
            ra1s2x2(2,2) = work(n+4);
            ia1s2x2(1,1) = iwork(1);
            ia1s2x2(2,1) = iwork(2);
            ia1s2x2(1,2) = iwork(3);
            ia1s2x2(2,2) = iwork(4);
            internaldlaebz<Precision>(3, itmax, n, 2, 2, atoli, rtoli, pivmin, d, e, work, ia1s2, ra1s2x2, ra1s2, iout, ia1s2x2, w, iblock, iinfo);
            iwork(5) = ia1s2(1);
            iwork(6) = ia1s2(2);
            work(n+5) = ra1s2(1);
            work(n+6) = ra1s2(2);
            work(n+1) = ra1s2x2(1,1);
            work(n+2) = ra1s2x2(2,1);
            work(n+3) = ra1s2x2(1,2);
            work(n+4) = ra1s2x2(2,2);
            iwork(1) = ia1s2x2(1,1);
            iwork(2) = ia1s2x2(2,1);
            iwork(3) = ia1s2x2(1,2);
            iwork(4) = ia1s2x2(2,2);
            if( iwork(6)==iu )
            {
                wl = work(n+1);
                wlu = work(n+3);
                nwl = iwork(1);
                wu = work(n+4);
                wul = work(n+2);
                nwu = iwork(4);
            }
            else
            {
                wl = work(n+2);
                wlu = work(n+4);
                nwl = iwork(2);
                wu = work(n+3);
                wul = work(n+1);
                nwu = iwork(3);
            }
            if( nwl<0 || nwl>=n || nwu<1 || nwu>n )
            {
                errorcode = 4;
                result = false;
                return result;
            }
        }
        else
        {
            
            //
            // RANGE='A' or 'V' -- Set ATOLI
            //
            tnorm = amp::maximum<Precision>(amp::abs<Precision>(d(1))+amp::abs<Precision>(e(1)), amp::abs<Precision>(d(n))+amp::abs<Precision>(e(n-1)));
            for(j=2; j<=n-1; j++)
            {
                tnorm = amp::maximum<Precision>(tnorm, amp::abs<Precision>(d(j))+amp::abs<Precision>(e(j-1))+amp::abs<Precision>(e(j)));
            }
            if( abstol<=0 )
            {
                atoli = ulp*tnorm;
            }
            else
            {
                atoli = abstol;
            }
            if( irange==2 )
            {
                wl = vl;
                wu = vu;
            }
            else
            {
                wl = 0;
                wu = 0;
            }
        }
        
        //
        // Find Eigenvalues -- Loop Over Blocks and recompute NWL and NWU.
        // NWL accumulates the number of eigenvalues .le. WL,
        // NWU accumulates the number of eigenvalues .le. WU
        //
        m = 0;
        iend = 0;
        errorcode = 0;
        nwl = 0;
        nwu = 0;
        for(jb=1; jb<=nsplit; jb++)
        {
            ioff = iend;
            ibegin = ioff+1;
            iend = isplit(jb);
            iin = iend-ioff;
            if( iin==1 )
            {
                
                //
                // Special Case -- IIN=1
                //
                if( irange==1 || wl>=d(ibegin)-pivmin )
                {
                    nwl = nwl+1;
                }
                if( irange==1 || wu>=d(ibegin)-pivmin )
                {
                    nwu = nwu+1;
                }
                if( irange==1 || wl<d(ibegin)-pivmin && wu>=d(ibegin)-pivmin )
                {
                    m = m+1;
                    w(m) = d(ibegin);
                    iblock(m) = jb;
                }
            }
            else
            {
                
                //
                // General Case -- IIN > 1
                //
                // Compute Gershgorin Interval
                // and use it as the initial interval
                //
                gu = d(ibegin);
                gl = d(ibegin);
                tmp1 = 0;
                for(j=ibegin; j<=iend-1; j++)
                {
                    tmp2 = amp::abs<Precision>(e(j));
                    gu = amp::maximum<Precision>(gu, d(j)+tmp1+tmp2);
                    gl = amp::minimum<Precision>(gl, d(j)-tmp1-tmp2);
                    tmp1 = tmp2;
                }
                gu = amp::maximum<Precision>(gu, d(iend)+tmp1);
                gl = amp::minimum<Precision>(gl, d(iend)-tmp1);
                bnorm = amp::maximum<Precision>(amp::abs<Precision>(gl), amp::abs<Precision>(gu));
                gl = gl-fudge*bnorm*ulp*iin-fudge*pivmin;
                gu = gu+fudge*bnorm*ulp*iin+fudge*pivmin;
                
                //
                // Compute ATOLI for the current submatrix
                //
                if( abstol<=0 )
                {
                    atoli = ulp*amp::maximum<Precision>(amp::abs<Precision>(gl), amp::abs<Precision>(gu));
                }
                else
                {
                    atoli = abstol;
                }
                if( irange>1 )
                {
                    if( gu<wl )
                    {
                        nwl = nwl+iin;
                        nwu = nwu+iin;
                        continue;
                    }
                    gl = amp::maximum<Precision>(gl, wl);
                    gu = amp::minimum<Precision>(gu, wu);
                    if( gl>=gu )
                    {
                        continue;
                    }
                }
                
                //
                // Set Up Initial Interval
                //
                work(n+1) = gl;
                work(n+iin+1) = gu;
                
                //
                // Calling DLAEBZ
                //
                // CALL DLAEBZ( 1, 0, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
                //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
                //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IM,
                //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
                //
                for(tmpi=1; tmpi<=iin; tmpi++)
                {
                    ra1siin(tmpi) = d(ibegin-1+tmpi);
                    if( ibegin-1+tmpi<n )
                    {
                        ra2siin(tmpi) = e(ibegin-1+tmpi);
                    }
                    ra3siin(tmpi) = work(ibegin-1+tmpi);
                    ra1siinx2(tmpi,1) = work(n+tmpi);
                    ra1siinx2(tmpi,2) = work(n+tmpi+iin);
                    ra4siin(tmpi) = work(n+2*iin+tmpi);
                    rworkspace(tmpi) = w(m+tmpi);
                    iworkspace(tmpi) = iblock(m+tmpi);
                    ia1siinx2(tmpi,1) = iwork(tmpi);
                    ia1siinx2(tmpi,2) = iwork(tmpi+iin);
                }
                internaldlaebz<Precision>(1, 0, iin, iin, 1, atoli, rtoli, pivmin, ra1siin, ra2siin, ra3siin, idumma, ra1siinx2, ra4siin, im, ia1siinx2, rworkspace, iworkspace, iinfo);
                for(tmpi=1; tmpi<=iin; tmpi++)
                {
                    work(n+tmpi) = ra1siinx2(tmpi,1);
                    work(n+tmpi+iin) = ra1siinx2(tmpi,2);
                    work(n+2*iin+tmpi) = ra4siin(tmpi);
                    w(m+tmpi) = rworkspace(tmpi);
                    iblock(m+tmpi) = iworkspace(tmpi);
                    iwork(tmpi) = ia1siinx2(tmpi,1);
                    iwork(tmpi+iin) = ia1siinx2(tmpi,2);
                }
                nwl = nwl+iwork(1);
                nwu = nwu+iwork(iin+1);
                iwoff = m-iwork(1);
                
                //
                // Compute Eigenvalues
                //
                itmax = amp::ceil<Precision>((amp::log<Precision>(gu-gl+pivmin)-amp::log<Precision>(pivmin))/amp::log<Precision>(amp::ampf<Precision>(2)))+2;
                
                //
                // Calling DLAEBZ
                //
                //CALL DLAEBZ( 2, ITMAX, IN, IN, 1, NB, ATOLI, RTOLI, PIVMIN,
                //    D( IBEGIN ), E( IBEGIN ), WORK( IBEGIN ),
                //    IDUMMA, WORK( N+1 ), WORK( N+2*IN+1 ), IOUT,
                //    IWORK, W( M+1 ), IBLOCK( M+1 ), IINFO )
                //
                for(tmpi=1; tmpi<=iin; tmpi++)
                {
                    ra1siin(tmpi) = d(ibegin-1+tmpi);
                    if( ibegin-1+tmpi<n )
                    {
                        ra2siin(tmpi) = e(ibegin-1+tmpi);
                    }
                    ra3siin(tmpi) = work(ibegin-1+tmpi);
                    ra1siinx2(tmpi,1) = work(n+tmpi);
                    ra1siinx2(tmpi,2) = work(n+tmpi+iin);
                    ra4siin(tmpi) = work(n+2*iin+tmpi);
                    rworkspace(tmpi) = w(m+tmpi);
                    iworkspace(tmpi) = iblock(m+tmpi);
                    ia1siinx2(tmpi,1) = iwork(tmpi);
                    ia1siinx2(tmpi,2) = iwork(tmpi+iin);
                }
                internaldlaebz<Precision>(2, itmax, iin, iin, 1, atoli, rtoli, pivmin, ra1siin, ra2siin, ra3siin, idumma, ra1siinx2, ra4siin, iout, ia1siinx2, rworkspace, iworkspace, iinfo);
                for(tmpi=1; tmpi<=iin; tmpi++)
                {
                    work(n+tmpi) = ra1siinx2(tmpi,1);
                    work(n+tmpi+iin) = ra1siinx2(tmpi,2);
                    work(n+2*iin+tmpi) = ra4siin(tmpi);
                    w(m+tmpi) = rworkspace(tmpi);
                    iblock(m+tmpi) = iworkspace(tmpi);
                    iwork(tmpi) = ia1siinx2(tmpi,1);
                    iwork(tmpi+iin) = ia1siinx2(tmpi,2);
                }
                
                //
                // Copy Eigenvalues Into W and IBLOCK
                // Use -JB for block number for unconverged eigenvalues.
                //
                for(j=1; j<=iout; j++)
                {
                    tmp1 = amp::ampf<Precision>("0.5")*(work(j+n)+work(j+iin+n));
                    
                    //
                    // Flag non-convergence.
                    //
                    if( j>iout-iinfo )
                    {
                        ncnvrg = true;
                        ib = -jb;
                    }
                    else
                    {
                        ib = jb;
                    }
                    for(je=iwork(j)+1+iwoff; je<=iwork(j+iin)+iwoff; je++)
                    {
                        w(je) = tmp1;
                        iblock(je) = ib;
                    }
                }
                m = m+im;
            }
        }
        
        //
        // If RANGE='I', then (WL,WU) contains eigenvalues NWL+1,...,NWU
        // If NWL+1 < IL or NWU > IU, discard extra eigenvalues.
        //
        if( irange==3 )
        {
            im = 0;
            idiscl = il-1-nwl;
            idiscu = nwu-iu;
            if( idiscl>0 || idiscu>0 )
            {
                for(je=1; je<=m; je++)
                {
                    if( w(je)<=wlu && idiscl>0 )
                    {
                        idiscl = idiscl-1;
                    }
                    else
                    {
                        if( w(je)>=wul && idiscu>0 )
                        {
                            idiscu = idiscu-1;
                        }
                        else
                        {
                            im = im+1;
                            w(im) = w(je);
                            iblock(im) = iblock(je);
                        }
                    }
                }
                m = im;
            }
            if( idiscl>0 || idiscu>0 )
            {
                
                //
                // Code to deal with effects of bad arithmetic:
                // Some low eigenvalues to be discarded are not in (WL,WLU],
                // or high eigenvalues to be discarded are not in (WUL,WU]
                // so just kill off the smallest IDISCL/largest IDISCU
                // eigenvalues, by simply finding the smallest/largest
                // eigenvalue(s).
                //
                // (If N(w) is monotone non-decreasing, this should never
                //  happen.)
                //
                if( idiscl>0 )
                {
                    wkill = wu;
                    for(jdisc=1; jdisc<=idiscl; jdisc++)
                    {
                        iw = 0;
                        for(je=1; je<=m; je++)
                        {
                            if( iblock(je)!=0 && (w(je)<wkill || iw==0) )
                            {
                                iw = je;
                                wkill = w(je);
                            }
                        }
                        iblock(iw) = 0;
                    }
                }
                if( idiscu>0 )
                {
                    wkill = wl;
                    for(jdisc=1; jdisc<=idiscu; jdisc++)
                    {
                        iw = 0;
                        for(je=1; je<=m; je++)
                        {
                            if( iblock(je)!=0 && (w(je)>wkill || iw==0) )
                            {
                                iw = je;
                                wkill = w(je);
                            }
                        }
                        iblock(iw) = 0;
                    }
                }
                im = 0;
                for(je=1; je<=m; je++)
                {
                    if( iblock(je)!=0 )
                    {
                        im = im+1;
                        w(im) = w(je);
                        iblock(im) = iblock(je);
                    }
                }
                m = im;
            }
            if( idiscl<0 || idiscu<0 )
            {
                toofew = true;
            }
        }
        
        //
        // If ORDER='B', do nothing -- the eigenvalues are already sorted
        //    by block.
        // If ORDER='E', sort the eigenvalues from smallest to largest
        //
        if( iorder==1 && nsplit>1 )
        {
            for(je=1; je<=m-1; je++)
            {
                ie = 0;
                tmp1 = w(je);
                for(j=je+1; j<=m; j++)
                {
                    if( w(j)<tmp1 )
                    {
                        ie = j;
                        tmp1 = w(j);
                    }
                }
                if( ie!=0 )
                {
                    itmp1 = iblock(ie);
                    w(ie) = w(je);
                    iblock(ie) = iblock(je);
                    w(je) = tmp1;
                    iblock(je) = itmp1;
                }
            }
        }
        for(j=1; j<=m; j++)
        {
            w(j) = w(j)*scalefactor;
        }
        errorcode = 0;
        if( ncnvrg )
        {
            errorcode = errorcode+1;
        }
        if( toofew )
        {
            errorcode = errorcode+2;
        }
        result = errorcode==0;
        return result;
    }


    template<unsigned int Precision>
    void internaldstein(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        const int& m,
        ap::template_1d_array< amp::ampf<Precision> > w,
        const ap::template_1d_array< int >& iblock,
        const ap::template_1d_array< int >& isplit,
        ap::template_2d_array< amp::ampf<Precision> >& z,
        ap::template_1d_array< int >& ifail,
        int& info)
    {
        int maxits;
        int extra;
        int b1;
        int blksiz;
        int bn;
        int gpind;
        int i;
        int iinfo;
        int its;
        int j;
        int j1;
        int jblk;
        int jmax;
        int nblk;
        int nrmchk;
        amp::ampf<Precision> dtpcrt;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> eps1;
        amp::ampf<Precision> nrm;
        amp::ampf<Precision> onenrm;
        amp::ampf<Precision> ortol;
        amp::ampf<Precision> pertol;
        amp::ampf<Precision> scl;
        amp::ampf<Precision> sep;
        amp::ampf<Precision> tol;
        amp::ampf<Precision> xj;
        amp::ampf<Precision> xjm;
        amp::ampf<Precision> ztr;
        ap::template_1d_array< amp::ampf<Precision> > work1;
        ap::template_1d_array< amp::ampf<Precision> > work2;
        ap::template_1d_array< amp::ampf<Precision> > work3;
        ap::template_1d_array< amp::ampf<Precision> > work4;
        ap::template_1d_array< amp::ampf<Precision> > work5;
        ap::template_1d_array< int > iwork;
        bool tmpcriterion;
        int ti;
        int i1;
        int i2;
        amp::ampf<Precision> v;


        maxits = 5;
        extra = 2;
        work1.setbounds(1, ap::maxint(n, 1));
        work2.setbounds(1, ap::maxint(n-1, 1));
        work3.setbounds(1, ap::maxint(n, 1));
        work4.setbounds(1, ap::maxint(n, 1));
        work5.setbounds(1, ap::maxint(n, 1));
        iwork.setbounds(1, ap::maxint(n, 1));
        ifail.setbounds(1, ap::maxint(m, 1));
        z.setbounds(1, ap::maxint(n, 1), 1, ap::maxint(m, 1));
        
        //
        // Test the input parameters.
        //
        info = 0;
        for(i=1; i<=m; i++)
        {
            ifail(i) = 0;
        }
        if( n<0 )
        {
            info = -1;
            return;
        }
        if( m<0 || m>n )
        {
            info = -4;
            return;
        }
        for(j=2; j<=m; j++)
        {
            if( iblock(j)<iblock(j-1) )
            {
                info = -6;
                break;
            }
            if( iblock(j)==iblock(j-1) && w(j)<w(j-1) )
            {
                info = -5;
                break;
            }
        }
        if( info!=0 )
        {
            return;
        }
        
        //
        // Quick return if possible
        //
        if( n==0 || m==0 )
        {
            return;
        }
        if( n==1 )
        {
            z(1,1) = 1;
            return;
        }
        
        //
        // Some preparations
        //
        ti = n-1;
        amp::vmove(work1.getvector(1, ti), e.getvector(1, ti));
        e.setbounds(1, n);
        amp::vmove(e.getvector(1, ti), work1.getvector(1, ti));
        amp::vmove(work1.getvector(1, m), w.getvector(1, m));
        w.setbounds(1, n);
        amp::vmove(w.getvector(1, m), work1.getvector(1, m));
        
        //
        // Get machine constants.
        //
        eps = amp::ampf<Precision>::getAlgoPascalEpsilon();
        
        //
        // Compute eigenvectors of matrix blocks.
        //
        j1 = 1;
        for(nblk=1; nblk<=iblock(m); nblk++)
        {
            
            //
            // Find starting and ending indices of block nblk.
            //
            if( nblk==1 )
            {
                b1 = 1;
            }
            else
            {
                b1 = isplit(nblk-1)+1;
            }
            bn = isplit(nblk);
            blksiz = bn-b1+1;
            if( blksiz!=1 )
            {
                
                //
                // Compute reorthogonalization criterion and stopping criterion.
                //
                gpind = b1;
                onenrm = amp::abs<Precision>(d(b1))+amp::abs<Precision>(e(b1));
                onenrm = amp::maximum<Precision>(onenrm, amp::abs<Precision>(d(bn))+amp::abs<Precision>(e(bn-1)));
                for(i=b1+1; i<=bn-1; i++)
                {
                    onenrm = amp::maximum<Precision>(onenrm, amp::abs<Precision>(d(i))+amp::abs<Precision>(e(i-1))+amp::abs<Precision>(e(i)));
                }
                ortol = amp::ampf<Precision>("0.001")*onenrm;
                dtpcrt = amp::sqrt<Precision>(amp::ampf<Precision>("0.1")/blksiz);
            }
            
            //
            // Loop through eigenvalues of block nblk.
            //
            jblk = 0;
            for(j=j1; j<=m; j++)
            {
                if( iblock(j)!=nblk )
                {
                    j1 = j;
                    break;
                }
                jblk = jblk+1;
                xj = w(j);
                if( blksiz==1 )
                {
                    
                    //
                    // Skip all the work if the block size is one.
                    //
                    work1(1) = 1;
                }
                else
                {
                    
                    //
                    // If eigenvalues j and j-1 are too close, add a relatively
                    // small perturbation.
                    //
                    if( jblk>1 )
                    {
                        eps1 = amp::abs<Precision>(eps*xj);
                        pertol = 10*eps1;
                        sep = xj-xjm;
                        if( sep<pertol )
                        {
                            xj = xjm+pertol;
                        }
                    }
                    its = 0;
                    nrmchk = 0;
                    
                    //
                    // Get random starting vector.
                    //
                    for(ti=1; ti<=blksiz; ti++)
                    {
                        work1(ti) = 2*amp::ampf<Precision>::getRandom()-1;
                    }
                    
                    //
                    // Copy the matrix T so it won't be destroyed in factorization.
                    //
                    for(ti=1; ti<=blksiz-1; ti++)
                    {
                        work2(ti) = e(b1+ti-1);
                        work3(ti) = e(b1+ti-1);
                        work4(ti) = d(b1+ti-1);
                    }
                    work4(blksiz) = d(b1+blksiz-1);
                    
                    //
                    // Compute LU factors with partial pivoting  ( PT = LU )
                    //
                    tol = 0;
                    tdininternaldlagtf<Precision>(blksiz, work4, xj, work2, work3, tol, work5, iwork, iinfo);
                    
                    //
                    // Update iteration count.
                    //
                    do
                    {
                        its = its+1;
                        if( its>maxits )
                        {
                            
                            //
                            // If stopping criterion was not satisfied, update info and
                            // store eigenvector number in array ifail.
                            //
                            info = info+1;
                            ifail(info) = j;
                            break;
                        }
                        
                        //
                        // Normalize and scale the righthand side vector Pb.
                        //
                        v = 0;
                        for(ti=1; ti<=blksiz; ti++)
                        {
                            v = v+amp::abs<Precision>(work1(ti));
                        }
                        scl = blksiz*onenrm*amp::maximum<Precision>(eps, amp::abs<Precision>(work4(blksiz)))/v;
                        amp::vmul(work1.getvector(1, blksiz), scl);
                        
                        //
                        // Solve the system LU = Pb.
                        //
                        tdininternaldlagts<Precision>(blksiz, work4, work2, work3, work5, iwork, work1, tol, iinfo);
                        
                        //
                        // Reorthogonalize by modified Gram-Schmidt if eigenvalues are
                        // close enough.
                        //
                        if( jblk!=1 )
                        {
                            if( amp::abs<Precision>(xj-xjm)>ortol )
                            {
                                gpind = j;
                            }
                            if( gpind!=j )
                            {
                                for(i=gpind; i<=j-1; i++)
                                {
                                    i1 = b1;
                                    i2 = b1+blksiz-1;
                                    ztr = amp::vdotproduct(work1.getvector(1, blksiz), z.getcolumn(i, i1, i2));
                                    amp::vsub(work1.getvector(1, blksiz), z.getcolumn(i, i1, i2), ztr);
                                }
                            }
                        }
                        
                        //
                        // Check the infinity norm of the iterate.
                        //
                        jmax = blas::vectoridxabsmax<Precision>(work1, 1, blksiz);
                        nrm = amp::abs<Precision>(work1(jmax));
                        
                        //
                        // Continue for additional iterations after norm reaches
                        // stopping criterion.
                        //
                        tmpcriterion = false;
                        if( nrm<dtpcrt )
                        {
                            tmpcriterion = true;
                        }
                        else
                        {
                            nrmchk = nrmchk+1;
                            if( nrmchk<extra+1 )
                            {
                                tmpcriterion = true;
                            }
                        }
                    }
                    while( tmpcriterion );
                    
                    //
                    // Accept iterate as jth eigenvector.
                    //
                    scl = 1/blas::vectornorm2<Precision>(work1, 1, blksiz);
                    jmax = blas::vectoridxabsmax<Precision>(work1, 1, blksiz);
                    if( work1(jmax)<0 )
                    {
                        scl = -scl;
                    }
                    amp::vmul(work1.getvector(1, blksiz), scl);
                }
                for(i=1; i<=n; i++)
                {
                    z(i,j) = 0;
                }
                for(i=1; i<=blksiz; i++)
                {
                    z(b1+i-1,j) = work1(i);
                }
                
                //
                // Save the shift to check eigenvalue spacing at next
                // iteration.
                //
                xjm = xj;
            }
        }
    }


    template<unsigned int Precision>
    bool tridiagonalevd(ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< amp::ampf<Precision> > e,
        int n,
        int zneeded,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        int maxit;
        int i;
        int ii;
        int iscale;
        int j;
        int jtot;
        int k;
        int t;
        int l;
        int l1;
        int lend;
        int lendm1;
        int lendp1;
        int lendsv;
        int lm1;
        int lsv;
        int m;
        int mm;
        int mm1;
        int nm1;
        int nmaxit;
        int tmpint;
        amp::ampf<Precision> anorm;
        amp::ampf<Precision> b;
        amp::ampf<Precision> c;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> eps2;
        amp::ampf<Precision> f;
        amp::ampf<Precision> g;
        amp::ampf<Precision> p;
        amp::ampf<Precision> r;
        amp::ampf<Precision> rt1;
        amp::ampf<Precision> rt2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> safmax;
        amp::ampf<Precision> safmin;
        amp::ampf<Precision> ssfmax;
        amp::ampf<Precision> ssfmin;
        amp::ampf<Precision> tst;
        amp::ampf<Precision> tmp;
        ap::template_1d_array< amp::ampf<Precision> > work1;
        ap::template_1d_array< amp::ampf<Precision> > work2;
        ap::template_1d_array< amp::ampf<Precision> > workc;
        ap::template_1d_array< amp::ampf<Precision> > works;
        ap::template_1d_array< amp::ampf<Precision> > wtemp;
        bool gotoflag;
        int zrows;
        bool wastranspose;


        ap::ap_error::make_assertion(zneeded>=0 && zneeded<=3);
        
        //
        // Quick return if possible
        //
        if( zneeded<0 || zneeded>3 )
        {
            result = false;
            return result;
        }
        result = true;
        if( n==0 )
        {
            return result;
        }
        if( n==1 )
        {
            if( zneeded==2 || zneeded==3 )
            {
                z.setbounds(1, 1, 1, 1);
                z(1,1) = 1;
            }
            return result;
        }
        maxit = 30;
        
        //
        // Initialize arrays
        //
        wtemp.setbounds(1, n);
        work1.setbounds(1, n-1);
        work2.setbounds(1, n-1);
        workc.setbounds(1, n);
        works.setbounds(1, n);
        
        //
        // Determine the unit roundoff and over/underflow thresholds.
        //
        eps = amp::ampf<Precision>::getAlgoPascalEpsilon();
        eps2 = amp::sqr<Precision>(eps);
        safmin = amp::ampf<Precision>::getAlgoPascalMinNumber();
        safmax = amp::ampf<Precision>::getAlgoPascalMaxNumber();
        ssfmax = amp::sqrt<Precision>(safmax)/3;
        ssfmin = amp::sqrt<Precision>(safmin)/eps2;
        
        //
        // Prepare Z
        //
        // Here we are using transposition to get rid of column operations
        //
        //
        wastranspose = false;
        if( zneeded==0 )
        {
            zrows = 0;
        }
        if( zneeded==1 )
        {
            zrows = n;
        }
        if( zneeded==2 )
        {
            zrows = n;
        }
        if( zneeded==3 )
        {
            zrows = 1;
        }
        if( zneeded==1 )
        {
            wastranspose = true;
            blas::inplacetranspose<Precision>(z, 1, n, 1, n, wtemp);
        }
        if( zneeded==2 )
        {
            wastranspose = true;
            z.setbounds(1, n, 1, n);
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++)
                {
                    if( i==j )
                    {
                        z(i,j) = 1;
                    }
                    else
                    {
                        z(i,j) = 0;
                    }
                }
            }
        }
        if( zneeded==3 )
        {
            wastranspose = false;
            z.setbounds(1, 1, 1, n);
            for(j=1; j<=n; j++)
            {
                if( j==1 )
                {
                    z(1,j) = 1;
                }
                else
                {
                    z(1,j) = 0;
                }
            }
        }
        nmaxit = n*maxit;
        jtot = 0;
        
        //
        // Determine where the matrix splits and choose QL or QR iteration
        // for each block, according to whether top or bottom diagonal
        // element is smaller.
        //
        l1 = 1;
        nm1 = n-1;
        while( true )
        {
            if( l1>n )
            {
                break;
            }
            if( l1>1 )
            {
                e(l1-1) = 0;
            }
            gotoflag = false;
            if( l1<=nm1 )
            {
                for(m=l1; m<=nm1; m++)
                {
                    tst = amp::abs<Precision>(e(m));
                    if( tst==0 )
                    {
                        gotoflag = true;
                        break;
                    }
                    if( tst<=amp::sqrt<Precision>(amp::abs<Precision>(d(m)))*amp::sqrt<Precision>(amp::abs<Precision>(d(m+1)))*eps )
                    {
                        e(m) = 0;
                        gotoflag = true;
                        break;
                    }
                }
            }
            if( !gotoflag )
            {
                m = n;
            }
            
            //
            // label 30:
            //
            l = l1;
            lsv = l;
            lend = m;
            lendsv = lend;
            l1 = m+1;
            if( lend==l )
            {
                continue;
            }
            
            //
            // Scale submatrix in rows and columns L to LEND
            //
            if( l==lend )
            {
                anorm = amp::abs<Precision>(d(l));
            }
            else
            {
                anorm = amp::maximum<Precision>(amp::abs<Precision>(d(l))+amp::abs<Precision>(e(l)), amp::abs<Precision>(e(lend-1))+amp::abs<Precision>(d(lend)));
                for(i=l+1; i<=lend-1; i++)
                {
                    anorm = amp::maximum<Precision>(anorm, amp::abs<Precision>(d(i))+amp::abs<Precision>(e(i))+amp::abs<Precision>(e(i-1)));
                }
            }
            iscale = 0;
            if( anorm==0 )
            {
                continue;
            }
            if( anorm>ssfmax )
            {
                iscale = 1;
                tmp = ssfmax/anorm;
                tmpint = lend-1;
                amp::vmul(d.getvector(l, lend), tmp);
                amp::vmul(e.getvector(l, tmpint), tmp);
            }
            if( anorm<ssfmin )
            {
                iscale = 2;
                tmp = ssfmin/anorm;
                tmpint = lend-1;
                amp::vmul(d.getvector(l, lend), tmp);
                amp::vmul(e.getvector(l, tmpint), tmp);
            }
            
            //
            // Choose between QL and QR iteration
            //
            if( amp::abs<Precision>(d(lend))<amp::abs<Precision>(d(l)) )
            {
                lend = lsv;
                l = lendsv;
            }
            if( lend>l )
            {
                
                //
                // QL Iteration
                //
                // Look for small subdiagonal element.
                //
                while( true )
                {
                    gotoflag = false;
                    if( l!=lend )
                    {
                        lendm1 = lend-1;
                        for(m=l; m<=lendm1; m++)
                        {
                            tst = amp::sqr<Precision>(amp::abs<Precision>(e(m)));
                            if( tst<=eps2*amp::abs<Precision>(d(m))*amp::abs<Precision>(d(m+1))+safmin )
                            {
                                gotoflag = true;
                                break;
                            }
                        }
                    }
                    if( !gotoflag )
                    {
                        m = lend;
                    }
                    if( m<lend )
                    {
                        e(m) = 0;
                    }
                    p = d(l);
                    if( m!=l )
                    {
                        
                        //
                        // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                        // to compute its eigensystem.
                        //
                        if( m==l+1 )
                        {
                            if( zneeded>0 )
                            {
                                tdevdev2<Precision>(d(l), e(l), d(l+1), rt1, rt2, c, s);
                                work1(l) = c;
                                work2(l) = s;
                                workc(1) = work1(l);
                                works(1) = work2(l);
                                if( !wastranspose )
                                {
                                    rotations::applyrotationsfromtheright<Precision>(false, 1, zrows, l, l+1, workc, works, z, wtemp);
                                }
                                else
                                {
                                    rotations::applyrotationsfromtheleft<Precision>(false, l, l+1, 1, zrows, workc, works, z, wtemp);
                                }
                            }
                            else
                            {
                                tdevde2<Precision>(d(l), e(l), d(l+1), rt1, rt2);
                            }
                            d(l) = rt1;
                            d(l+1) = rt2;
                            e(l) = 0;
                            l = l+2;
                            if( l<=lend )
                            {
                                continue;
                            }
                            
                            //
                            // GOTO 140
                            //
                            break;
                        }
                        if( jtot==nmaxit )
                        {
                            
                            //
                            // GOTO 140
                            //
                            break;
                        }
                        jtot = jtot+1;
                        
                        //
                        // Form shift.
                        //
                        g = (d(l+1)-p)/(2*e(l));
                        r = tdevdpythag<Precision>(g, amp::ampf<Precision>(1));
                        g = d(m)-p+e(l)/(g+tdevdextsign<Precision>(r, g));
                        s = 1;
                        c = 1;
                        p = 0;
                        
                        //
                        // Inner loop
                        //
                        mm1 = m-1;
                        for(i=mm1; i>=l; i--)
                        {
                            f = s*e(i);
                            b = c*e(i);
                            rotations::generaterotation<Precision>(g, f, c, s, r);
                            if( i!=m-1 )
                            {
                                e(i+1) = r;
                            }
                            g = d(i+1)-p;
                            r = (d(i)-g)*s+2*c*b;
                            p = s*r;
                            d(i+1) = g+p;
                            g = c*r-b;
                            
                            //
                            // If eigenvectors are desired, then save rotations.
                            //
                            if( zneeded>0 )
                            {
                                work1(i) = c;
                                work2(i) = -s;
                            }
                        }
                        
                        //
                        // If eigenvectors are desired, then apply saved rotations.
                        //
                        if( zneeded>0 )
                        {
                            for(i=l; i<=m-1; i++)
                            {
                                workc(i-l+1) = work1(i);
                                works(i-l+1) = work2(i);
                            }
                            if( !wastranspose )
                            {
                                rotations::applyrotationsfromtheright<Precision>(false, 1, zrows, l, m, workc, works, z, wtemp);
                            }
                            else
                            {
                                rotations::applyrotationsfromtheleft<Precision>(false, l, m, 1, zrows, workc, works, z, wtemp);
                            }
                        }
                        d(l) = d(l)-p;
                        e(l) = g;
                        continue;
                    }
                    
                    //
                    // Eigenvalue found.
                    //
                    d(l) = p;
                    l = l+1;
                    if( l<=lend )
                    {
                        continue;
                    }
                    break;
                }
            }
            else
            {
                
                //
                // QR Iteration
                //
                // Look for small superdiagonal element.
                //
                while( true )
                {
                    gotoflag = false;
                    if( l!=lend )
                    {
                        lendp1 = lend+1;
                        for(m=l; m>=lendp1; m--)
                        {
                            tst = amp::sqr<Precision>(amp::abs<Precision>(e(m-1)));
                            if( tst<=eps2*amp::abs<Precision>(d(m))*amp::abs<Precision>(d(m-1))+safmin )
                            {
                                gotoflag = true;
                                break;
                            }
                        }
                    }
                    if( !gotoflag )
                    {
                        m = lend;
                    }
                    if( m>lend )
                    {
                        e(m-1) = 0;
                    }
                    p = d(l);
                    if( m!=l )
                    {
                        
                        //
                        // If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
                        // to compute its eigensystem.
                        //
                        if( m==l-1 )
                        {
                            if( zneeded>0 )
                            {
                                tdevdev2<Precision>(d(l-1), e(l-1), d(l), rt1, rt2, c, s);
                                work1(m) = c;
                                work2(m) = s;
                                workc(1) = c;
                                works(1) = s;
                                if( !wastranspose )
                                {
                                    rotations::applyrotationsfromtheright<Precision>(true, 1, zrows, l-1, l, workc, works, z, wtemp);
                                }
                                else
                                {
                                    rotations::applyrotationsfromtheleft<Precision>(true, l-1, l, 1, zrows, workc, works, z, wtemp);
                                }
                            }
                            else
                            {
                                tdevde2<Precision>(d(l-1), e(l-1), d(l), rt1, rt2);
                            }
                            d(l-1) = rt1;
                            d(l) = rt2;
                            e(l-1) = 0;
                            l = l-2;
                            if( l>=lend )
                            {
                                continue;
                            }
                            break;
                        }
                        if( jtot==nmaxit )
                        {
                            break;
                        }
                        jtot = jtot+1;
                        
                        //
                        // Form shift.
                        //
                        g = (d(l-1)-p)/(2*e(l-1));
                        r = tdevdpythag<Precision>(g, amp::ampf<Precision>(1));
                        g = d(m)-p+e(l-1)/(g+tdevdextsign<Precision>(r, g));
                        s = 1;
                        c = 1;
                        p = 0;
                        
                        //
                        // Inner loop
                        //
                        lm1 = l-1;
                        for(i=m; i<=lm1; i++)
                        {
                            f = s*e(i);
                            b = c*e(i);
                            rotations::generaterotation<Precision>(g, f, c, s, r);
                            if( i!=m )
                            {
                                e(i-1) = r;
                            }
                            g = d(i)-p;
                            r = (d(i+1)-g)*s+2*c*b;
                            p = s*r;
                            d(i) = g+p;
                            g = c*r-b;
                            
                            //
                            // If eigenvectors are desired, then save rotations.
                            //
                            if( zneeded>0 )
                            {
                                work1(i) = c;
                                work2(i) = s;
                            }
                        }
                        
                        //
                        // If eigenvectors are desired, then apply saved rotations.
                        //
                        if( zneeded>0 )
                        {
                            mm = l-m+1;
                            for(i=m; i<=l-1; i++)
                            {
                                workc(i-m+1) = work1(i);
                                works(i-m+1) = work2(i);
                            }
                            if( !wastranspose )
                            {
                                rotations::applyrotationsfromtheright<Precision>(true, 1, zrows, m, l, workc, works, z, wtemp);
                            }
                            else
                            {
                                rotations::applyrotationsfromtheleft<Precision>(true, m, l, 1, zrows, workc, works, z, wtemp);
                            }
                        }
                        d(l) = d(l)-p;
                        e(lm1) = g;
                        continue;
                    }
                    
                    //
                    // Eigenvalue found.
                    //
                    d(l) = p;
                    l = l-1;
                    if( l>=lend )
                    {
                        continue;
                    }
                    break;
                }
            }
            
            //
            // Undo scaling if necessary
            //
            if( iscale==1 )
            {
                tmp = anorm/ssfmax;
                tmpint = lendsv-1;
                amp::vmul(d.getvector(lsv, lendsv), tmp);
                amp::vmul(e.getvector(lsv, tmpint), tmp);
            }
            if( iscale==2 )
            {
                tmp = anorm/ssfmin;
                tmpint = lendsv-1;
                amp::vmul(d.getvector(lsv, lendsv), tmp);
                amp::vmul(e.getvector(lsv, tmpint), tmp);
            }
            
            //
            // Check for no convergence to an eigenvalue after a total
            // of N*MAXIT iterations.
            //
            if( jtot>=nmaxit )
            {
                result = false;
                if( wastranspose )
                {
                    blas::inplacetranspose<Precision>(z, 1, n, 1, n, wtemp);
                }
                return result;
            }
        }
        
        //
        // Order eigenvalues and eigenvectors.
        //
        if( zneeded==0 )
        {
            
            //
            // Sort
            //
            if( n==1 )
            {
                return result;
            }
            if( n==2 )
            {
                if( d(1)>d(2) )
                {
                    tmp = d(1);
                    d(1) = d(2);
                    d(2) = tmp;
                }
                return result;
            }
            i = 2;
            do
            {
                t = i;
                while( t!=1 )
                {
                    k = t/2;
                    if( d(k)>=d(t) )
                    {
                        t = 1;
                    }
                    else
                    {
                        tmp = d(k);
                        d(k) = d(t);
                        d(t) = tmp;
                        t = k;
                    }
                }
                i = i+1;
            }
            while( i<=n );
            i = n-1;
            do
            {
                tmp = d(i+1);
                d(i+1) = d(1);
                d(+1) = tmp;
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
                            if( d(k+1)>d(k) )
                            {
                                k = k+1;
                            }
                        }
                        if( d(t)>=d(k) )
                        {
                            t = 0;
                        }
                        else
                        {
                            tmp = d(k);
                            d(k) = d(t);
                            d(t) = tmp;
                            t = k;
                        }
                    }
                }
                i = i-1;
            }
            while( i>=1 );
        }
        else
        {
            
            //
            // Use Selection Sort to minimize swaps of eigenvectors
            //
            for(ii=2; ii<=n; ii++)
            {
                i = ii-1;
                k = i;
                p = d(i);
                for(j=ii; j<=n; j++)
                {
                    if( d(j)<p )
                    {
                        k = j;
                        p = d(j);
                    }
                }
                if( k!=i )
                {
                    d(k) = d(i);
                    d(i) = p;
                    if( wastranspose )
                    {
                        amp::vmove(wtemp.getvector(1, n), z.getrow(i, 1, n));
                        amp::vmove(z.getrow(i, 1, n), z.getrow(k, 1, n));
                        amp::vmove(z.getrow(k, 1, n), wtemp.getvector(1, n));
                    }
                    else
                    {
                        amp::vmove(wtemp.getvector(1, zrows), z.getcolumn(i, 1, zrows));
                        amp::vmove(z.getcolumn(i, 1, zrows), z.getcolumn(k, 1, zrows));
                        amp::vmove(z.getcolumn(k, 1, zrows), wtemp.getvector(1, zrows));
                    }
                }
            }
            if( wastranspose )
            {
                blas::inplacetranspose<Precision>(z, 1, n, 1, n, wtemp);
            }
        }
        return result;
    }


    /*************************************************************************
    DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
       [  A   B  ]
       [  B   C  ].
    On return, RT1 is the eigenvalue of larger absolute value, and RT2
    is the eigenvalue of smaller absolute value.

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         October 31, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void tdevde2(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        amp::ampf<Precision>& rt1,
        amp::ampf<Precision>& rt2)
    {
        amp::ampf<Precision> ab;
        amp::ampf<Precision> acmn;
        amp::ampf<Precision> acmx;
        amp::ampf<Precision> adf;
        amp::ampf<Precision> df;
        amp::ampf<Precision> rt;
        amp::ampf<Precision> sm;
        amp::ampf<Precision> tb;


        sm = a+c;
        df = a-c;
        adf = amp::abs<Precision>(df);
        tb = b+b;
        ab = amp::abs<Precision>(tb);
        if( amp::abs<Precision>(a)>amp::abs<Precision>(c) )
        {
            acmx = a;
            acmn = c;
        }
        else
        {
            acmx = c;
            acmn = a;
        }
        if( adf>ab )
        {
            rt = adf*amp::sqrt<Precision>(1+amp::sqr<Precision>(ab/adf));
        }
        else
        {
            if( adf<ab )
            {
                rt = ab*amp::sqrt<Precision>(1+amp::sqr<Precision>(adf/ab));
            }
            else
            {
                
                //
                // Includes case AB=ADF=0
                //
                rt = ab*amp::sqrt<Precision>(amp::ampf<Precision>(2));
            }
        }
        if( sm<0 )
        {
            rt1 = amp::ampf<Precision>("0.5")*(sm-rt);
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            rt2 = acmx/rt1*acmn-b/rt1*b;
        }
        else
        {
            if( sm>0 )
            {
                rt1 = amp::ampf<Precision>("0.5")*(sm+rt);
                
                //
                // Order of execution important.
                // To get fully accurate smaller eigenvalue,
                // next line needs to be executed in higher precision.
                //
                rt2 = acmx/rt1*acmn-b/rt1*b;
            }
            else
            {
                
                //
                // Includes case RT1 = RT2 = 0
                //
                rt1 = amp::ampf<Precision>("0.5")*rt;
                rt2 = -amp::ampf<Precision>("0.5")*rt;
            }
        }
    }


    /*************************************************************************
    DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix

       [  A   B  ]
       [  B   C  ].

    On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
    eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
    eigenvector for RT1, giving the decomposition

       [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
       [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].


      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         October 31, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void tdevdev2(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        amp::ampf<Precision>& rt1,
        amp::ampf<Precision>& rt2,
        amp::ampf<Precision>& cs1,
        amp::ampf<Precision>& sn1)
    {
        int sgn1;
        int sgn2;
        amp::ampf<Precision> ab;
        amp::ampf<Precision> acmn;
        amp::ampf<Precision> acmx;
        amp::ampf<Precision> acs;
        amp::ampf<Precision> adf;
        amp::ampf<Precision> cs;
        amp::ampf<Precision> ct;
        amp::ampf<Precision> df;
        amp::ampf<Precision> rt;
        amp::ampf<Precision> sm;
        amp::ampf<Precision> tb;
        amp::ampf<Precision> tn;


        
        //
        // Compute the eigenvalues
        //
        sm = a+c;
        df = a-c;
        adf = amp::abs<Precision>(df);
        tb = b+b;
        ab = amp::abs<Precision>(tb);
        if( amp::abs<Precision>(a)>amp::abs<Precision>(c) )
        {
            acmx = a;
            acmn = c;
        }
        else
        {
            acmx = c;
            acmn = a;
        }
        if( adf>ab )
        {
            rt = adf*amp::sqrt<Precision>(1+amp::sqr<Precision>(ab/adf));
        }
        else
        {
            if( adf<ab )
            {
                rt = ab*amp::sqrt<Precision>(1+amp::sqr<Precision>(adf/ab));
            }
            else
            {
                
                //
                // Includes case AB=ADF=0
                //
                rt = ab*amp::sqrt<Precision>(amp::ampf<Precision>(2));
            }
        }
        if( sm<0 )
        {
            rt1 = amp::ampf<Precision>("0.5")*(sm-rt);
            sgn1 = -1;
            
            //
            // Order of execution important.
            // To get fully accurate smaller eigenvalue,
            // next line needs to be executed in higher precision.
            //
            rt2 = acmx/rt1*acmn-b/rt1*b;
        }
        else
        {
            if( sm>0 )
            {
                rt1 = amp::ampf<Precision>("0.5")*(sm+rt);
                sgn1 = 1;
                
                //
                // Order of execution important.
                // To get fully accurate smaller eigenvalue,
                // next line needs to be executed in higher precision.
                //
                rt2 = acmx/rt1*acmn-b/rt1*b;
            }
            else
            {
                
                //
                // Includes case RT1 = RT2 = 0
                //
                rt1 = amp::ampf<Precision>("0.5")*rt;
                rt2 = -amp::ampf<Precision>("0.5")*rt;
                sgn1 = 1;
            }
        }
        
        //
        // Compute the eigenvector
        //
        if( df>=0 )
        {
            cs = df+rt;
            sgn2 = 1;
        }
        else
        {
            cs = df-rt;
            sgn2 = -1;
        }
        acs = amp::abs<Precision>(cs);
        if( acs>ab )
        {
            ct = -tb/cs;
            sn1 = 1/amp::sqrt<Precision>(1+ct*ct);
            cs1 = ct*sn1;
        }
        else
        {
            if( ab==0 )
            {
                cs1 = 1;
                sn1 = 0;
            }
            else
            {
                tn = -cs/tb;
                cs1 = 1/amp::sqrt<Precision>(1+tn*tn);
                sn1 = tn*cs1;
            }
        }
        if( sgn1==sgn2 )
        {
            tn = cs1;
            cs1 = -sn1;
            sn1 = tn;
        }
    }


    /*************************************************************************
    Internal routine
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> tdevdpythag(amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( amp::abs<Precision>(a)<amp::abs<Precision>(b) )
        {
            result = amp::abs<Precision>(b)*amp::sqrt<Precision>(1+amp::sqr<Precision>(a/b));
        }
        else
        {
            result = amp::abs<Precision>(a)*amp::sqrt<Precision>(1+amp::sqr<Precision>(b/a));
        }
        return result;
    }


    /*************************************************************************
    Internal routine
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> tdevdextsign(amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( b>=0 )
        {
            result = amp::abs<Precision>(a);
        }
        else
        {
            result = -amp::abs<Precision>(a);
        }
        return result;
    }


    template<unsigned int Precision>
    void tdininternaldlagtf(const int& n,
        ap::template_1d_array< amp::ampf<Precision> >& a,
        const amp::ampf<Precision>& lambda,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        const amp::ampf<Precision>& tol,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_1d_array< int >& iin,
        int& info)
    {
        int k;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> mult;
        amp::ampf<Precision> piv1;
        amp::ampf<Precision> piv2;
        amp::ampf<Precision> scale1;
        amp::ampf<Precision> scale2;
        amp::ampf<Precision> temp;
        amp::ampf<Precision> tl;


        info = 0;
        if( n<0 )
        {
            info = -1;
            return;
        }
        if( n==0 )
        {
            return;
        }
        a(1) = a(1)-lambda;
        iin(n) = 0;
        if( n==1 )
        {
            if( a(1)==0 )
            {
                iin(1) = 1;
            }
            return;
        }
        eps = amp::ampf<Precision>::getAlgoPascalEpsilon();
        tl = amp::maximum<Precision>(tol, eps);
        scale1 = amp::abs<Precision>(a(1))+amp::abs<Precision>(b(1));
        for(k=1; k<=n-1; k++)
        {
            a(k+1) = a(k+1)-lambda;
            scale2 = amp::abs<Precision>(c(k))+amp::abs<Precision>(a(k+1));
            if( k<n-1 )
            {
                scale2 = scale2+amp::abs<Precision>(b(k+1));
            }
            if( a(k)==0 )
            {
                piv1 = 0;
            }
            else
            {
                piv1 = amp::abs<Precision>(a(k))/scale1;
            }
            if( c(k)==0 )
            {
                iin(k) = 0;
                piv2 = 0;
                scale1 = scale2;
                if( k<n-1 )
                {
                    d(k) = 0;
                }
            }
            else
            {
                piv2 = amp::abs<Precision>(c(k))/scale2;
                if( piv2<=piv1 )
                {
                    iin(k) = 0;
                    scale1 = scale2;
                    c(k) = c(k)/a(k);
                    a(k+1) = a(k+1)-c(k)*b(k);
                    if( k<n-1 )
                    {
                        d(k) = 0;
                    }
                }
                else
                {
                    iin(k) = 1;
                    mult = a(k)/c(k);
                    a(k) = c(k);
                    temp = a(k+1);
                    a(k+1) = b(k)-mult*temp;
                    if( k<n-1 )
                    {
                        d(k) = b(k+1);
                        b(k+1) = -mult*d(k);
                    }
                    b(k) = temp;
                    c(k) = mult;
                }
            }
            if( amp::maximum<Precision>(piv1, piv2)<=tl && iin(n)==0 )
            {
                iin(n) = k;
            }
        }
        if( amp::abs<Precision>(a(n))<=scale1*tl && iin(n)==0 )
        {
            iin(n) = n;
        }
    }


    template<unsigned int Precision>
    void tdininternaldlagts(const int& n,
        const ap::template_1d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< int >& iin,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        amp::ampf<Precision>& tol,
        int& info)
    {
        int k;
        amp::ampf<Precision> absak;
        amp::ampf<Precision> ak;
        amp::ampf<Precision> bignum;
        amp::ampf<Precision> eps;
        amp::ampf<Precision> pert;
        amp::ampf<Precision> sfmin;
        amp::ampf<Precision> temp;


        info = 0;
        if( n<0 )
        {
            info = -1;
            return;
        }
        if( n==0 )
        {
            return;
        }
        eps = amp::ampf<Precision>::getAlgoPascalEpsilon();
        sfmin = amp::ampf<Precision>::getAlgoPascalMinNumber();
        bignum = 1/sfmin;
        if( tol<=0 )
        {
            tol = amp::abs<Precision>(a(1));
            if( n>1 )
            {
                tol = amp::maximum<Precision>(tol, amp::maximum<Precision>(amp::abs<Precision>(a(2)), amp::abs<Precision>(b(1))));
            }
            for(k=3; k<=n; k++)
            {
                tol = amp::maximum<Precision>(tol, amp::maximum<Precision>(amp::abs<Precision>(a(k)), amp::maximum<Precision>(amp::abs<Precision>(b(k-1)), amp::abs<Precision>(d(k-2)))));
            }
            tol = tol*eps;
            if( tol==0 )
            {
                tol = eps;
            }
        }
        for(k=2; k<=n; k++)
        {
            if( iin(k-1)==0 )
            {
                y(k) = y(k)-c(k-1)*y(k-1);
            }
            else
            {
                temp = y(k-1);
                y(k-1) = y(k);
                y(k) = temp-c(k-1)*y(k);
            }
        }
        for(k=n; k>=1; k--)
        {
            if( k<=n-2 )
            {
                temp = y(k)-b(k)*y(k+1)-d(k)*y(k+2);
            }
            else
            {
                if( k==n-1 )
                {
                    temp = y(k)-b(k)*y(k+1);
                }
                else
                {
                    temp = y(k);
                }
            }
            ak = a(k);
            pert = amp::abs<Precision>(tol);
            if( ak<0 )
            {
                pert = -pert;
            }
            while( true )
            {
                absak = amp::abs<Precision>(ak);
                if( absak<1 )
                {
                    if( absak<sfmin )
                    {
                        if( absak==0 || amp::abs<Precision>(temp)*sfmin>absak )
                        {
                            ak = ak+pert;
                            pert = 2*pert;
                            continue;
                        }
                        else
                        {
                            temp = temp*bignum;
                            ak = ak*bignum;
                        }
                    }
                    else
                    {
                        if( amp::abs<Precision>(temp)>absak*bignum )
                        {
                            ak = ak+pert;
                            pert = 2*pert;
                            continue;
                        }
                    }
                }
                break;
            }
            y(k) = temp/ak;
        }
    }


    template<unsigned int Precision>
    void internaldlaebz(const int& ijob,
        const int& nitmax,
        const int& n,
        const int& mmax,
        const int& minp,
        const amp::ampf<Precision>& abstol,
        const amp::ampf<Precision>& reltol,
        const amp::ampf<Precision>& pivmin,
        const ap::template_1d_array< amp::ampf<Precision> >& d,
        const ap::template_1d_array< amp::ampf<Precision> >& e,
        const ap::template_1d_array< amp::ampf<Precision> >& e2,
        ap::template_1d_array< int >& nval,
        ap::template_2d_array< amp::ampf<Precision> >& ab,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        int& mout,
        ap::template_2d_array< int >& nab,
        ap::template_1d_array< amp::ampf<Precision> >& work,
        ap::template_1d_array< int >& iwork,
        int& info)
    {
        int itmp1;
        int itmp2;
        int j;
        int ji;
        int jit;
        int jp;
        int kf;
        int kfnew;
        int kl;
        int klnew;
        amp::ampf<Precision> tmp1;
        amp::ampf<Precision> tmp2;


        info = 0;
        if( ijob<1 || ijob>3 )
        {
            info = -1;
            return;
        }
        
        //
        // Initialize NAB
        //
        if( ijob==1 )
        {
            
            //
            // Compute the number of eigenvalues in the initial intervals.
            //
            mout = 0;
            
            //
            //DIR$ NOVECTOR
            //
            for(ji=1; ji<=minp; ji++)
            {
                for(jp=1; jp<=2; jp++)
                {
                    tmp1 = d(1)-ab(ji,jp);
                    if( amp::abs<Precision>(tmp1)<pivmin )
                    {
                        tmp1 = -pivmin;
                    }
                    nab(ji,jp) = 0;
                    if( tmp1<=0 )
                    {
                        nab(ji,jp) = 1;
                    }
                    for(j=2; j<=n; j++)
                    {
                        tmp1 = d(j)-e2(j-1)/tmp1-ab(ji,jp);
                        if( amp::abs<Precision>(tmp1)<pivmin )
                        {
                            tmp1 = -pivmin;
                        }
                        if( tmp1<=0 )
                        {
                            nab(ji,jp) = nab(ji,jp)+1;
                        }
                    }
                }
                mout = mout+nab(ji,2)-nab(ji,1);
            }
            return;
        }
        
        //
        // Initialize for loop
        //
        // KF and KL have the following meaning:
        //   Intervals 1,...,KF-1 have converged.
        //   Intervals KF,...,KL  still need to be refined.
        //
        kf = 1;
        kl = minp;
        
        //
        // If IJOB=2, initialize C.
        // If IJOB=3, use the user-supplied starting point.
        //
        if( ijob==2 )
        {
            for(ji=1; ji<=minp; ji++)
            {
                c(ji) = amp::ampf<Precision>("0.5")*(ab(ji,1)+ab(ji,2));
            }
        }
        
        //
        // Iteration loop
        //
        for(jit=1; jit<=nitmax; jit++)
        {
            
            //
            // Loop over intervals
            //
            //
            // Serial Version of the loop
            //
            klnew = kl;
            for(ji=kf; ji<=kl; ji++)
            {
                
                //
                // Compute N(w), the number of eigenvalues less than w
                //
                tmp1 = c(ji);
                tmp2 = d(1)-tmp1;
                itmp1 = 0;
                if( tmp2<=pivmin )
                {
                    itmp1 = 1;
                    tmp2 = amp::minimum<Precision>(tmp2, -pivmin);
                }
                
                //
                // A series of compiler directives to defeat vectorization
                // for the next loop
                //
                //*$PL$ CMCHAR=' '
                //CDIR$          NEXTSCALAR
                //C$DIR          SCALAR
                //CDIR$          NEXT SCALAR
                //CVD$L          NOVECTOR
                //CDEC$          NOVECTOR
                //CVD$           NOVECTOR
                //*VDIR          NOVECTOR
                //*VOCL          LOOP,SCALAR
                //CIBM           PREFER SCALAR
                //*$PL$ CMCHAR='*'
                //
                for(j=2; j<=n; j++)
                {
                    tmp2 = d(j)-e2(j-1)/tmp2-tmp1;
                    if( tmp2<=pivmin )
                    {
                        itmp1 = itmp1+1;
                        tmp2 = amp::minimum<Precision>(tmp2, -pivmin);
                    }
                }
                if( ijob<=2 )
                {
                    
                    //
                    // IJOB=2: Choose all intervals containing eigenvalues.
                    //
                    // Insure that N(w) is monotone
                    //
                    itmp1 = ap::minint(nab(ji,2), ap::maxint(nab(ji,1), itmp1));
                    
                    //
                    // Update the Queue -- add intervals if both halves
                    // contain eigenvalues.
                    //
                    if( itmp1==nab(ji,2) )
                    {
                        
                        //
                        // No eigenvalue in the upper interval:
                        // just use the lower interval.
                        //
                        ab(ji,2) = tmp1;
                    }
                    else
                    {
                        if( itmp1==nab(ji,1) )
                        {
                            
                            //
                            // No eigenvalue in the lower interval:
                            // just use the upper interval.
                            //
                            ab(ji,1) = tmp1;
                        }
                        else
                        {
                            if( klnew<mmax )
                            {
                                
                                //
                                // Eigenvalue in both intervals -- add upper to queue.
                                //
                                klnew = klnew+1;
                                ab(klnew,2) = ab(ji,2);
                                nab(klnew,2) = nab(ji,2);
                                ab(klnew,1) = tmp1;
                                nab(klnew,1) = itmp1;
                                ab(ji,2) = tmp1;
                                nab(ji,2) = itmp1;
                            }
                            else
                            {
                                info = mmax+1;
                                return;
                            }
                        }
                    }
                }
                else
                {
                    
                    //
                    // IJOB=3: Binary search.  Keep only the interval
                    // containing  w  s.t. N(w) = NVAL
                    //
                    if( itmp1<=nval(ji) )
                    {
                        ab(ji,1) = tmp1;
                        nab(ji,1) = itmp1;
                    }
                    if( itmp1>=nval(ji) )
                    {
                        ab(ji,2) = tmp1;
                        nab(ji,2) = itmp1;
                    }
                }
            }
            kl = klnew;
            
            //
            // Check for convergence
            //
            kfnew = kf;
            for(ji=kf; ji<=kl; ji++)
            {
                tmp1 = amp::abs<Precision>(ab(ji,2)-ab(ji,1));
                tmp2 = amp::maximum<Precision>(amp::abs<Precision>(ab(ji,2)), amp::abs<Precision>(ab(ji,1)));
                if( tmp1<amp::maximum<Precision>(abstol, amp::maximum<Precision>(pivmin, reltol*tmp2)) || nab(ji,1)>=nab(ji,2) )
                {
                    
                    //
                    // Converged -- Swap with position KFNEW,
                    // then increment KFNEW
                    //
                    if( ji>kfnew )
                    {
                        tmp1 = ab(ji,1);
                        tmp2 = ab(ji,2);
                        itmp1 = nab(ji,1);
                        itmp2 = nab(ji,2);
                        ab(ji,1) = ab(kfnew,1);
                        ab(ji,2) = ab(kfnew,2);
                        nab(ji,1) = nab(kfnew,1);
                        nab(ji,2) = nab(kfnew,2);
                        ab(kfnew,1) = tmp1;
                        ab(kfnew,2) = tmp2;
                        nab(kfnew,1) = itmp1;
                        nab(kfnew,2) = itmp2;
                        if( ijob==3 )
                        {
                            itmp1 = nval(ji);
                            nval(ji) = nval(kfnew);
                            nval(kfnew) = itmp1;
                        }
                    }
                    kfnew = kfnew+1;
                }
            }
            kf = kfnew;
            
            //
            // Choose Midpoints
            //
            for(ji=kf; ji<=kl; ji++)
            {
                c(ji) = amp::ampf<Precision>("0.5")*(ab(ji,1)+ab(ji,2));
            }
            
            //
            // If no more intervals to refine, quit.
            //
            if( kf>kl )
            {
                break;
            }
        }
        
        //
        // Converged
        //
        info = ap::maxint(kl+1-kf, 0);
        mout = kl;
    }


    /*************************************************************************
    Internal subroutine

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         June 30, 1999
    *************************************************************************/
    template<unsigned int Precision>
    void internaltrevc(const ap::template_2d_array< amp::ampf<Precision> >& t,
        int n,
        int side,
        int howmny,
        ap::template_1d_array< bool > vselect,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr,
        int& m,
        int& info)
    {
        bool allv;
        bool bothv;
        bool leftv;
        bool over;
        bool pair;
        bool rightv;
        bool somev;
        int i;
        int ierr;
        int ii;
        int ip;
        int iis;
        int j;
        int j1;
        int j2;
        int jnxt;
        int k;
        int ki;
        int n2;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> bignum;
        amp::ampf<Precision> emax;
        amp::ampf<Precision> ovfl;
        amp::ampf<Precision> rec;
        amp::ampf<Precision> remax;
        amp::ampf<Precision> scl;
        amp::ampf<Precision> smin;
        amp::ampf<Precision> smlnum;
        amp::ampf<Precision> ulp;
        amp::ampf<Precision> unfl;
        amp::ampf<Precision> vcrit;
        amp::ampf<Precision> vmax;
        amp::ampf<Precision> wi;
        amp::ampf<Precision> wr;
        amp::ampf<Precision> xnorm;
        ap::template_2d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_1d_array< amp::ampf<Precision> > temp;
        ap::template_2d_array< amp::ampf<Precision> > temp11;
        ap::template_2d_array< amp::ampf<Precision> > temp22;
        ap::template_2d_array< amp::ampf<Precision> > temp11b;
        ap::template_2d_array< amp::ampf<Precision> > temp21b;
        ap::template_2d_array< amp::ampf<Precision> > temp12b;
        ap::template_2d_array< amp::ampf<Precision> > temp22b;
        bool skipflag;
        int k1;
        int k2;
        int k3;
        int k4;
        amp::ampf<Precision> vt;
        ap::template_1d_array< bool > rswap4;
        ap::template_1d_array< bool > zswap4;
        ap::template_2d_array< int > ipivot44;
        ap::template_1d_array< amp::ampf<Precision> > civ4;
        ap::template_1d_array< amp::ampf<Precision> > crv4;


        x.setbounds(1, 2, 1, 2);
        temp11.setbounds(1, 1, 1, 1);
        temp11b.setbounds(1, 1, 1, 1);
        temp21b.setbounds(1, 2, 1, 1);
        temp12b.setbounds(1, 1, 1, 2);
        temp22b.setbounds(1, 2, 1, 2);
        temp22.setbounds(1, 2, 1, 2);
        work.setbounds(1, 3*n);
        temp.setbounds(1, n);
        rswap4.setbounds(1, 4);
        zswap4.setbounds(1, 4);
        ipivot44.setbounds(1, 4, 1, 4);
        civ4.setbounds(1, 4);
        crv4.setbounds(1, 4);
        if( howmny!=1 )
        {
            if( side==1 || side==3 )
            {
                vr.setbounds(1, n, 1, n);
            }
            if( side==2 || side==3 )
            {
                vl.setbounds(1, n, 1, n);
            }
        }
        
        //
        // Decode and test the input parameters
        //
        bothv = side==3;
        rightv = side==1 || bothv;
        leftv = side==2 || bothv;
        allv = howmny==2;
        over = howmny==1;
        somev = howmny==3;
        info = 0;
        if( n<0 )
        {
            info = -2;
            return;
        }
        if( !rightv && !leftv )
        {
            info = -3;
            return;
        }
        if( !allv && !over && !somev )
        {
            info = -4;
            return;
        }
        
        //
        // Set M to the number of columns required to store the selected
        // eigenvectors, standardize the array SELECT if necessary, and
        // test MM.
        //
        if( somev )
        {
            m = 0;
            pair = false;
            for(j=1; j<=n; j++)
            {
                if( pair )
                {
                    pair = false;
                    vselect(j) = false;
                }
                else
                {
                    if( j<n )
                    {
                        if( t(j+1,j)==0 )
                        {
                            if( vselect(j) )
                            {
                                m = m+1;
                            }
                        }
                        else
                        {
                            pair = true;
                            if( vselect(j) || vselect(j+1) )
                            {
                                vselect(j) = true;
                                m = m+2;
                            }
                        }
                    }
                    else
                    {
                        if( vselect(n) )
                        {
                            m = m+1;
                        }
                    }
                }
            }
        }
        else
        {
            m = n;
        }
        
        //
        // Quick return if possible.
        //
        if( n==0 )
        {
            return;
        }
        
        //
        // Set the constants to control overflow.
        //
        unfl = amp::ampf<Precision>::getAlgoPascalMinNumber();
        ovfl = 1/unfl;
        ulp = amp::ampf<Precision>::getAlgoPascalEpsilon();
        smlnum = unfl*(n/ulp);
        bignum = (1-ulp)/smlnum;
        
        //
        // Compute 1-norm of each column of strictly upper triangular
        // part of T to control overflow in triangular solver.
        //
        work(1) = 0;
        for(j=2; j<=n; j++)
        {
            work(j) = 0;
            for(i=1; i<=j-1; i++)
            {
                work(j) = work(j)+amp::abs<Precision>(t(i,j));
            }
        }
        
        //
        // Index IP is used to specify the real or complex eigenvalue:
        // IP = 0, real eigenvalue,
        //      1, first of conjugate complex pair: (wr,wi)
        //     -1, second of conjugate complex pair: (wr,wi)
        //
        n2 = 2*n;
        if( rightv )
        {
            
            //
            // Compute right eigenvectors.
            //
            ip = 0;
            iis = m;
            for(ki=n; ki>=1; ki--)
            {
                skipflag = false;
                if( ip==1 )
                {
                    skipflag = true;
                }
                else
                {
                    if( ki!=1 )
                    {
                        if( t(ki,ki-1)!=0 )
                        {
                            ip = -1;
                        }
                    }
                    if( somev )
                    {
                        if( ip==0 )
                        {
                            if( !vselect(ki) )
                            {
                                skipflag = true;
                            }
                        }
                        else
                        {
                            if( !vselect(ki-1) )
                            {
                                skipflag = true;
                            }
                        }
                    }
                }
                if( !skipflag )
                {
                    
                    //
                    // Compute the KI-th eigenvalue (WR,WI).
                    //
                    wr = t(ki,ki);
                    wi = 0;
                    if( ip!=0 )
                    {
                        wi = amp::sqrt<Precision>(amp::abs<Precision>(t(ki,ki-1)))*amp::sqrt<Precision>(amp::abs<Precision>(t(ki-1,ki)));
                    }
                    smin = amp::maximum<Precision>(ulp*(amp::abs<Precision>(wr)+amp::abs<Precision>(wi)), smlnum);
                    if( ip==0 )
                    {
                        
                        //
                        // Real right eigenvector
                        //
                        work(ki+n) = 1;
                        
                        //
                        // Form right-hand side
                        //
                        for(k=1; k<=ki-1; k++)
                        {
                            work(k+n) = -t(k,ki);
                        }
                        
                        //
                        // Solve the upper quasi-triangular system:
                        //   (T(1:KI-1,1:KI-1) - WR)*X = SCALE*WORK.
                        //
                        jnxt = ki-1;
                        for(j=ki-1; j>=1; j--)
                        {
                            if( j>jnxt )
                            {
                                continue;
                            }
                            j1 = j;
                            j2 = j;
                            jnxt = j-1;
                            if( j>1 )
                            {
                                if( t(j,j-1)!=0 )
                                {
                                    j1 = j-1;
                                    jnxt = j-2;
                                }
                            }
                            if( j1==j2 )
                            {
                                
                                //
                                // 1-by-1 diagonal block
                                //
                                temp11(1,1) = t(j,j);
                                temp11b(1,1) = work(j+n);
                                internalhsevdlaln2<Precision>(false, 1, 1, smin, amp::ampf<Precision>(1), temp11, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp11b, wr, amp::ampf<Precision>("0.0"), rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale X(1,1) to avoid overflow when updating
                                // the right-hand side.
                                //
                                if( xnorm>1 )
                                {
                                    if( work(j)>bignum/xnorm )
                                    {
                                        x(1,1) = x(1,1)/xnorm;
                                        scl = scl/xnorm;
                                    }
                                }
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    k1 = n+1;
                                    k2 = n+ki;
                                    amp::vmul(work.getvector(k1, k2), scl);
                                }
                                work(j+n) = x(1,1);
                                
                                //
                                // Update right-hand side
                                //
                                k1 = 1+n;
                                k2 = j-1+n;
                                k3 = j-1;
                                vt = -x(1,1);
                                amp::vadd(work.getvector(k1, k2), t.getcolumn(j, 1, k3), vt);
                            }
                            else
                            {
                                
                                //
                                // 2-by-2 diagonal block
                                //
                                temp22(1,1) = t(j-1,j-1);
                                temp22(1,2) = t(j-1,j);
                                temp22(2,1) = t(j,j-1);
                                temp22(2,2) = t(j,j);
                                temp21b(1,1) = work(j-1+n);
                                temp21b(2,1) = work(j+n);
                                internalhsevdlaln2<Precision>(false, 2, 1, smin, amp::ampf<Precision>("1.0"), temp22, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp21b, wr, amp::ampf<Precision>(0), rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale X(1,1) and X(2,1) to avoid overflow when
                                // updating the right-hand side.
                                //
                                if( xnorm>1 )
                                {
                                    beta = amp::maximum<Precision>(work(j-1), work(j));
                                    if( beta>bignum/xnorm )
                                    {
                                        x(1,1) = x(1,1)/xnorm;
                                        x(2,1) = x(2,1)/xnorm;
                                        scl = scl/xnorm;
                                    }
                                }
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    k1 = 1+n;
                                    k2 = ki+n;
                                    amp::vmul(work.getvector(k1, k2), scl);
                                }
                                work(j-1+n) = x(1,1);
                                work(j+n) = x(2,1);
                                
                                //
                                // Update right-hand side
                                //
                                k1 = 1+n;
                                k2 = j-2+n;
                                k3 = j-2;
                                k4 = j-1;
                                vt = -x(1,1);
                                amp::vadd(work.getvector(k1, k2), t.getcolumn(k4, 1, k3), vt);
                                vt = -x(2,1);
                                amp::vadd(work.getvector(k1, k2), t.getcolumn(j, 1, k3), vt);
                            }
                        }
                        
                        //
                        // Copy the vector x or Q*x to VR and normalize.
                        //
                        if( !over )
                        {
                            k1 = 1+n;
                            k2 = ki+n;
                            amp::vmove(vr.getcolumn(iis, 1, ki), work.getvector(k1, k2));
                            ii = blas::columnidxabsmax<Precision>(vr, 1, ki, iis);
                            remax = 1/amp::abs<Precision>(vr(ii,iis));
                            amp::vmul(vr.getcolumn(iis, 1, ki), remax);
                            for(k=ki+1; k<=n; k++)
                            {
                                vr(k,iis) = 0;
                            }
                        }
                        else
                        {
                            if( ki>1 )
                            {
                                amp::vmove(temp.getvector(1, n), vr.getcolumn(ki, 1, n));
                                blas::matrixvectormultiply<Precision>(vr, 1, n, 1, ki-1, false, work, 1+n, ki-1+n, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki+n));
                                amp::vmove(vr.getcolumn(ki, 1, n), temp.getvector(1, n));
                            }
                            ii = blas::columnidxabsmax<Precision>(vr, 1, n, ki);
                            remax = 1/amp::abs<Precision>(vr(ii,ki));
                            amp::vmul(vr.getcolumn(ki, 1, n), remax);
                        }
                    }
                    else
                    {
                        
                        //
                        // Complex right eigenvector.
                        //
                        // Initial solve
                        //     [ (T(KI-1,KI-1) T(KI-1,KI) ) - (WR + I* WI)]*X = 0.
                        //     [ (T(KI,KI-1)   T(KI,KI)   )               ]
                        //
                        if( amp::abs<Precision>(t(ki-1,ki))>=amp::abs<Precision>(t(ki,ki-1)) )
                        {
                            work(ki-1+n) = 1;
                            work(ki+n2) = wi/t(ki-1,ki);
                        }
                        else
                        {
                            work(ki-1+n) = -wi/t(ki,ki-1);
                            work(ki+n2) = 1;
                        }
                        work(ki+n) = 0;
                        work(ki-1+n2) = 0;
                        
                        //
                        // Form right-hand side
                        //
                        for(k=1; k<=ki-2; k++)
                        {
                            work(k+n) = -work(ki-1+n)*t(k,ki-1);
                            work(k+n2) = -work(ki+n2)*t(k,ki);
                        }
                        
                        //
                        // Solve upper quasi-triangular system:
                        // (T(1:KI-2,1:KI-2) - (WR+i*WI))*X = SCALE*(WORK+i*WORK2)
                        //
                        jnxt = ki-2;
                        for(j=ki-2; j>=1; j--)
                        {
                            if( j>jnxt )
                            {
                                continue;
                            }
                            j1 = j;
                            j2 = j;
                            jnxt = j-1;
                            if( j>1 )
                            {
                                if( t(j,j-1)!=0 )
                                {
                                    j1 = j-1;
                                    jnxt = j-2;
                                }
                            }
                            if( j1==j2 )
                            {
                                
                                //
                                // 1-by-1 diagonal block
                                //
                                temp11(1,1) = t(j,j);
                                temp12b(1,1) = work(j+n);
                                temp12b(1,2) = work(j+n+n);
                                internalhsevdlaln2<Precision>(false, 1, 2, smin, amp::ampf<Precision>("1.0"), temp11, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp12b, wr, wi, rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale X(1,1) and X(1,2) to avoid overflow when
                                // updating the right-hand side.
                                //
                                if( xnorm>1 )
                                {
                                    if( work(j)>bignum/xnorm )
                                    {
                                        x(1,1) = x(1,1)/xnorm;
                                        x(1,2) = x(1,2)/xnorm;
                                        scl = scl/xnorm;
                                    }
                                }
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    k1 = 1+n;
                                    k2 = ki+n;
                                    amp::vmul(work.getvector(k1, k2), scl);
                                    k1 = 1+n2;
                                    k2 = ki+n2;
                                    amp::vmul(work.getvector(k1, k2), scl);
                                }
                                work(j+n) = x(1,1);
                                work(j+n2) = x(1,2);
                                
                                //
                                // Update the right-hand side
                                //
                                k1 = 1+n;
                                k2 = j-1+n;
                                k3 = 1;
                                k4 = j-1;
                                vt = -x(1,1);
                                amp::vadd(work.getvector(k1, k2), t.getcolumn(j, k3, k4), vt);
                                k1 = 1+n2;
                                k2 = j-1+n2;
                                k3 = 1;
                                k4 = j-1;
                                vt = -x(1,2);
                                amp::vadd(work.getvector(k1, k2), t.getcolumn(j, k3, k4), vt);
                            }
                            else
                            {
                                
                                //
                                // 2-by-2 diagonal block
                                //
                                temp22(1,1) = t(j-1,j-1);
                                temp22(1,2) = t(j-1,j);
                                temp22(2,1) = t(j,j-1);
                                temp22(2,2) = t(j,j);
                                temp22b(1,1) = work(j-1+n);
                                temp22b(1,2) = work(j-1+n+n);
                                temp22b(2,1) = work(j+n);
                                temp22b(2,2) = work(j+n+n);
                                internalhsevdlaln2<Precision>(false, 2, 2, smin, amp::ampf<Precision>("1.0"), temp22, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp22b, wr, wi, rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale X to avoid overflow when updating
                                // the right-hand side.
                                //
                                if( xnorm>1 )
                                {
                                    beta = amp::maximum<Precision>(work(j-1), work(j));
                                    if( beta>bignum/xnorm )
                                    {
                                        rec = 1/xnorm;
                                        x(1,1) = x(1,1)*rec;
                                        x(1,2) = x(1,2)*rec;
                                        x(2,1) = x(2,1)*rec;
                                        x(2,2) = x(2,2)*rec;
                                        scl = scl*rec;
                                    }
                                }
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    amp::vmul(work.getvector(1+n, ki+n), scl);
                                    amp::vmul(work.getvector(1+n2, ki+n2), scl);
                                }
                                work(j-1+n) = x(1,1);
                                work(j+n) = x(2,1);
                                work(j-1+n2) = x(1,2);
                                work(j+n2) = x(2,2);
                                
                                //
                                // Update the right-hand side
                                //
                                vt = -x(1,1);
                                amp::vadd(work.getvector(n+1, n+j-2), t.getcolumn(j-1, 1, j-2), vt);
                                vt = -x(2,1);
                                amp::vadd(work.getvector(n+1, n+j-2), t.getcolumn(j, 1, j-2), vt);
                                vt = -x(1,2);
                                amp::vadd(work.getvector(n2+1, n2+j-2), t.getcolumn(j-1, 1, j-2), vt);
                                vt = -x(2,2);
                                amp::vadd(work.getvector(n2+1, n2+j-2), t.getcolumn(j, 1, j-2), vt);
                            }
                        }
                        
                        //
                        // Copy the vector x or Q*x to VR and normalize.
                        //
                        if( !over )
                        {
                            amp::vmove(vr.getcolumn(iis-1, 1, ki), work.getvector(n+1, n+ki));
                            amp::vmove(vr.getcolumn(iis, 1, ki), work.getvector(n2+1, n2+ki));
                            emax = 0;
                            for(k=1; k<=ki; k++)
                            {
                                emax = amp::maximum<Precision>(emax, amp::abs<Precision>(vr(k,iis-1))+amp::abs<Precision>(vr(k,iis)));
                            }
                            remax = 1/emax;
                            amp::vmul(vr.getcolumn(iis-1, 1, ki), remax);
                            amp::vmul(vr.getcolumn(iis, 1, ki), remax);
                            for(k=ki+1; k<=n; k++)
                            {
                                vr(k,iis-1) = 0;
                                vr(k,iis) = 0;
                            }
                        }
                        else
                        {
                            if( ki>2 )
                            {
                                amp::vmove(temp.getvector(1, n), vr.getcolumn(ki-1, 1, n));
                                blas::matrixvectormultiply<Precision>(vr, 1, n, 1, ki-2, false, work, 1+n, ki-2+n, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki-1+n));
                                amp::vmove(vr.getcolumn(ki-1, 1, n), temp.getvector(1, n));
                                amp::vmove(temp.getvector(1, n), vr.getcolumn(ki, 1, n));
                                blas::matrixvectormultiply<Precision>(vr, 1, n, 1, ki-2, false, work, 1+n2, ki-2+n2, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki+n2));
                                amp::vmove(vr.getcolumn(ki, 1, n), temp.getvector(1, n));
                            }
                            else
                            {
                                vt = work(ki-1+n);
                                amp::vmul(vr.getcolumn(ki-1, 1, n), vt);
                                vt = work(ki+n2);
                                amp::vmul(vr.getcolumn(ki, 1, n), vt);
                            }
                            emax = 0;
                            for(k=1; k<=n; k++)
                            {
                                emax = amp::maximum<Precision>(emax, amp::abs<Precision>(vr(k,ki-1))+amp::abs<Precision>(vr(k,ki)));
                            }
                            remax = 1/emax;
                            amp::vmul(vr.getcolumn(ki-1, 1, n), remax);
                            amp::vmul(vr.getcolumn(ki, 1, n), remax);
                        }
                    }
                    iis = iis-1;
                    if( ip!=0 )
                    {
                        iis = iis-1;
                    }
                }
                if( ip==1 )
                {
                    ip = 0;
                }
                if( ip==-1 )
                {
                    ip = 1;
                }
            }
        }
        if( leftv )
        {
            
            //
            // Compute left eigenvectors.
            //
            ip = 0;
            iis = 1;
            for(ki=1; ki<=n; ki++)
            {
                skipflag = false;
                if( ip==-1 )
                {
                    skipflag = true;
                }
                else
                {
                    if( ki!=n )
                    {
                        if( t(ki+1,ki)!=0 )
                        {
                            ip = 1;
                        }
                    }
                    if( somev )
                    {
                        if( !vselect(ki) )
                        {
                            skipflag = true;
                        }
                    }
                }
                if( !skipflag )
                {
                    
                    //
                    // Compute the KI-th eigenvalue (WR,WI).
                    //
                    wr = t(ki,ki);
                    wi = 0;
                    if( ip!=0 )
                    {
                        wi = amp::sqrt<Precision>(amp::abs<Precision>(t(ki,ki+1)))*amp::sqrt<Precision>(amp::abs<Precision>(t(ki+1,ki)));
                    }
                    smin = amp::maximum<Precision>(ulp*(amp::abs<Precision>(wr)+amp::abs<Precision>(wi)), smlnum);
                    if( ip==0 )
                    {
                        
                        //
                        // Real left eigenvector.
                        //
                        work(ki+n) = 1;
                        
                        //
                        // Form right-hand side
                        //
                        for(k=ki+1; k<=n; k++)
                        {
                            work(k+n) = -t(ki,k);
                        }
                        
                        //
                        // Solve the quasi-triangular system:
                        // (T(KI+1:N,KI+1:N) - WR)'*X = SCALE*WORK
                        //
                        vmax = 1;
                        vcrit = bignum;
                        jnxt = ki+1;
                        for(j=ki+1; j<=n; j++)
                        {
                            if( j<jnxt )
                            {
                                continue;
                            }
                            j1 = j;
                            j2 = j;
                            jnxt = j+1;
                            if( j<n )
                            {
                                if( t(j+1,j)!=0 )
                                {
                                    j2 = j+1;
                                    jnxt = j+2;
                                }
                            }
                            if( j1==j2 )
                            {
                                
                                //
                                // 1-by-1 diagonal block
                                //
                                // Scale if necessary to avoid overflow when forming
                                // the right-hand side.
                                //
                                if( work(j)>vcrit )
                                {
                                    rec = 1/vmax;
                                    amp::vmul(work.getvector(ki+n, n+n), rec);
                                    vmax = 1;
                                    vcrit = bignum;
                                }
                                vt = amp::vdotproduct(t.getcolumn(j, ki+1, j-1), work.getvector(ki+1+n, j-1+n));
                                work(j+n) = work(j+n)-vt;
                                
                                //
                                // Solve (T(J,J)-WR)'*X = WORK
                                //
                                temp11(1,1) = t(j,j);
                                temp11b(1,1) = work(j+n);
                                internalhsevdlaln2<Precision>(false, 1, 1, smin, amp::ampf<Precision>("1.0"), temp11, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp11b, wr, amp::ampf<Precision>(0), rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    amp::vmul(work.getvector(ki+n, n+n), scl);
                                }
                                work(j+n) = x(1,1);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(work(j+n)), vmax);
                                vcrit = bignum/vmax;
                            }
                            else
                            {
                                
                                //
                                // 2-by-2 diagonal block
                                //
                                // Scale if necessary to avoid overflow when forming
                                // the right-hand side.
                                //
                                beta = amp::maximum<Precision>(work(j), work(j+1));
                                if( beta>vcrit )
                                {
                                    rec = 1/vmax;
                                    amp::vmul(work.getvector(ki+n, n+n), rec);
                                    vmax = 1;
                                    vcrit = bignum;
                                }
                                vt = amp::vdotproduct(t.getcolumn(j, ki+1, j-1), work.getvector(ki+1+n, j-1+n));
                                work(j+n) = work(j+n)-vt;
                                vt = amp::vdotproduct(t.getcolumn(j+1, ki+1, j-1), work.getvector(ki+1+n, j-1+n));
                                work(j+1+n) = work(j+1+n)-vt;
                                
                                //
                                // Solve
                                //    [T(J,J)-WR   T(J,J+1)     ]'* X = SCALE*( WORK1 )
                                //    [T(J+1,J)    T(J+1,J+1)-WR]             ( WORK2 )
                                //
                                temp22(1,1) = t(j,j);
                                temp22(1,2) = t(j,j+1);
                                temp22(2,1) = t(j+1,j);
                                temp22(2,2) = t(j+1,j+1);
                                temp21b(1,1) = work(j+n);
                                temp21b(2,1) = work(j+1+n);
                                internalhsevdlaln2<Precision>(true, 2, 1, smin, amp::ampf<Precision>("1.0"), temp22, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp21b, wr, amp::ampf<Precision>(0), rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    amp::vmul(work.getvector(ki+n, n+n), scl);
                                }
                                work(j+n) = x(1,1);
                                work(j+1+n) = x(2,1);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(work(j+n)), amp::maximum<Precision>(amp::abs<Precision>(work(j+1+n)), vmax));
                                vcrit = bignum/vmax;
                            }
                        }
                        
                        //
                        // Copy the vector x or Q*x to VL and normalize.
                        //
                        if( !over )
                        {
                            amp::vmove(vl.getcolumn(iis, ki, n), work.getvector(ki+n, n+n));
                            ii = blas::columnidxabsmax<Precision>(vl, ki, n, iis);
                            remax = 1/amp::abs<Precision>(vl(ii,iis));
                            amp::vmul(vl.getcolumn(iis, ki, n), remax);
                            for(k=1; k<=ki-1; k++)
                            {
                                vl(k,iis) = 0;
                            }
                        }
                        else
                        {
                            if( ki<n )
                            {
                                amp::vmove(temp.getvector(1, n), vl.getcolumn(ki, 1, n));
                                blas::matrixvectormultiply<Precision>(vl, 1, n, ki+1, n, false, work, ki+1+n, n+n, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki+n));
                                amp::vmove(vl.getcolumn(ki, 1, n), temp.getvector(1, n));
                            }
                            ii = blas::columnidxabsmax<Precision>(vl, 1, n, ki);
                            remax = 1/amp::abs<Precision>(vl(ii,ki));
                            amp::vmul(vl.getcolumn(ki, 1, n), remax);
                        }
                    }
                    else
                    {
                        
                        //
                        // Complex left eigenvector.
                        //
                        // Initial solve:
                        //   ((T(KI,KI)    T(KI,KI+1) )' - (WR - I* WI))*X = 0.
                        //   ((T(KI+1,KI) T(KI+1,KI+1))                )
                        //
                        if( amp::abs<Precision>(t(ki,ki+1))>=amp::abs<Precision>(t(ki+1,ki)) )
                        {
                            work(ki+n) = wi/t(ki,ki+1);
                            work(ki+1+n2) = 1;
                        }
                        else
                        {
                            work(ki+n) = 1;
                            work(ki+1+n2) = -wi/t(ki+1,ki);
                        }
                        work(ki+1+n) = 0;
                        work(ki+n2) = 0;
                        
                        //
                        // Form right-hand side
                        //
                        for(k=ki+2; k<=n; k++)
                        {
                            work(k+n) = -work(ki+n)*t(ki,k);
                            work(k+n2) = -work(ki+1+n2)*t(ki+1,k);
                        }
                        
                        //
                        // Solve complex quasi-triangular system:
                        // ( T(KI+2,N:KI+2,N) - (WR-i*WI) )*X = WORK1+i*WORK2
                        //
                        vmax = 1;
                        vcrit = bignum;
                        jnxt = ki+2;
                        for(j=ki+2; j<=n; j++)
                        {
                            if( j<jnxt )
                            {
                                continue;
                            }
                            j1 = j;
                            j2 = j;
                            jnxt = j+1;
                            if( j<n )
                            {
                                if( t(j+1,j)!=0 )
                                {
                                    j2 = j+1;
                                    jnxt = j+2;
                                }
                            }
                            if( j1==j2 )
                            {
                                
                                //
                                // 1-by-1 diagonal block
                                //
                                // Scale if necessary to avoid overflow when
                                // forming the right-hand side elements.
                                //
                                if( work(j)>vcrit )
                                {
                                    rec = 1/vmax;
                                    amp::vmul(work.getvector(ki+n, n+n), rec);
                                    amp::vmul(work.getvector(ki+n2, n+n2), rec);
                                    vmax = 1;
                                    vcrit = bignum;
                                }
                                vt = amp::vdotproduct(t.getcolumn(j, ki+2, j-1), work.getvector(ki+2+n, j-1+n));
                                work(j+n) = work(j+n)-vt;
                                vt = amp::vdotproduct(t.getcolumn(j, ki+2, j-1), work.getvector(ki+2+n2, j-1+n2));
                                work(j+n2) = work(j+n2)-vt;
                                
                                //
                                // Solve (T(J,J)-(WR-i*WI))*(X11+i*X12)= WK+I*WK2
                                //
                                temp11(1,1) = t(j,j);
                                temp12b(1,1) = work(j+n);
                                temp12b(1,2) = work(j+n+n);
                                internalhsevdlaln2<Precision>(false, 1, 2, smin, amp::ampf<Precision>("1.0"), temp11, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp12b, wr, -wi, rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    amp::vmul(work.getvector(ki+n, n+n), scl);
                                    amp::vmul(work.getvector(ki+n2, n+n2), scl);
                                }
                                work(j+n) = x(1,1);
                                work(j+n2) = x(1,2);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(work(j+n)), amp::maximum<Precision>(amp::abs<Precision>(work(j+n2)), vmax));
                                vcrit = bignum/vmax;
                            }
                            else
                            {
                                
                                //
                                // 2-by-2 diagonal block
                                //
                                // Scale if necessary to avoid overflow when forming
                                // the right-hand side elements.
                                //
                                beta = amp::maximum<Precision>(work(j), work(j+1));
                                if( beta>vcrit )
                                {
                                    rec = 1/vmax;
                                    amp::vmul(work.getvector(ki+n, n+n), rec);
                                    amp::vmul(work.getvector(ki+n2, n+n2), rec);
                                    vmax = 1;
                                    vcrit = bignum;
                                }
                                vt = amp::vdotproduct(t.getcolumn(j, ki+2, j-1), work.getvector(ki+2+n, j-1+n));
                                work(j+n) = work(j+n)-vt;
                                vt = amp::vdotproduct(t.getcolumn(j, ki+2, j-1), work.getvector(ki+2+n2, j-1+n2));
                                work(j+n2) = work(j+n2)-vt;
                                vt = amp::vdotproduct(t.getcolumn(j+1, ki+2, j-1), work.getvector(ki+2+n, j-1+n));
                                work(j+1+n) = work(j+1+n)-vt;
                                vt = amp::vdotproduct(t.getcolumn(j+1, ki+2, j-1), work.getvector(ki+2+n2, j-1+n2));
                                work(j+1+n2) = work(j+1+n2)-vt;
                                
                                //
                                // Solve 2-by-2 complex linear equation
                                //   ([T(j,j)   T(j,j+1)  ]'-(wr-i*wi)*I)*X = SCALE*B
                                //   ([T(j+1,j) T(j+1,j+1)]             )
                                //
                                temp22(1,1) = t(j,j);
                                temp22(1,2) = t(j,j+1);
                                temp22(2,1) = t(j+1,j);
                                temp22(2,2) = t(j+1,j+1);
                                temp22b(1,1) = work(j+n);
                                temp22b(1,2) = work(j+n+n);
                                temp22b(2,1) = work(j+1+n);
                                temp22b(2,2) = work(j+1+n+n);
                                internalhsevdlaln2<Precision>(true, 2, 2, smin, amp::ampf<Precision>("1.0"), temp22, amp::ampf<Precision>("1.0"), amp::ampf<Precision>("1.0"), temp22b, wr, -wi, rswap4, zswap4, ipivot44, civ4, crv4, x, scl, xnorm, ierr);
                                
                                //
                                // Scale if necessary
                                //
                                if( scl!=1 )
                                {
                                    amp::vmul(work.getvector(ki+n, n+n), scl);
                                    amp::vmul(work.getvector(ki+n2, n+n2), scl);
                                }
                                work(j+n) = x(1,1);
                                work(j+n2) = x(1,2);
                                work(j+1+n) = x(2,1);
                                work(j+1+n2) = x(2,2);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(x(1,1)), vmax);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(x(1,2)), vmax);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(x(2,1)), vmax);
                                vmax = amp::maximum<Precision>(amp::abs<Precision>(x(2,2)), vmax);
                                vcrit = bignum/vmax;
                            }
                        }
                        
                        //
                        // Copy the vector x or Q*x to VL and normalize.
                        //
                        if( !over )
                        {
                            amp::vmove(vl.getcolumn(iis, ki, n), work.getvector(ki+n, n+n));
                            amp::vmove(vl.getcolumn(iis+1, ki, n), work.getvector(ki+n2, n+n2));
                            emax = 0;
                            for(k=ki; k<=n; k++)
                            {
                                emax = amp::maximum<Precision>(emax, amp::abs<Precision>(vl(k,iis))+amp::abs<Precision>(vl(k,iis+1)));
                            }
                            remax = 1/emax;
                            amp::vmul(vl.getcolumn(iis, ki, n), remax);
                            amp::vmul(vl.getcolumn(iis+1, ki, n), remax);
                            for(k=1; k<=ki-1; k++)
                            {
                                vl(k,iis) = 0;
                                vl(k,iis+1) = 0;
                            }
                        }
                        else
                        {
                            if( ki<n-1 )
                            {
                                amp::vmove(temp.getvector(1, n), vl.getcolumn(ki, 1, n));
                                blas::matrixvectormultiply<Precision>(vl, 1, n, ki+2, n, false, work, ki+2+n, n+n, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki+n));
                                amp::vmove(vl.getcolumn(ki, 1, n), temp.getvector(1, n));
                                amp::vmove(temp.getvector(1, n), vl.getcolumn(ki+1, 1, n));
                                blas::matrixvectormultiply<Precision>(vl, 1, n, ki+2, n, false, work, ki+2+n2, n+n2, amp::ampf<Precision>("1.0"), temp, 1, n, work(ki+1+n2));
                                amp::vmove(vl.getcolumn(ki+1, 1, n), temp.getvector(1, n));
                            }
                            else
                            {
                                vt = work(ki+n);
                                amp::vmul(vl.getcolumn(ki, 1, n), vt);
                                vt = work(ki+1+n2);
                                amp::vmul(vl.getcolumn(ki+1, 1, n), vt);
                            }
                            emax = 0;
                            for(k=1; k<=n; k++)
                            {
                                emax = amp::maximum<Precision>(emax, amp::abs<Precision>(vl(k,ki))+amp::abs<Precision>(vl(k,ki+1)));
                            }
                            remax = 1/emax;
                            amp::vmul(vl.getcolumn(ki, 1, n), remax);
                            amp::vmul(vl.getcolumn(ki+1, 1, n), remax);
                        }
                    }
                    iis = iis+1;
                    if( ip!=0 )
                    {
                        iis = iis+1;
                    }
                }
                if( ip==-1 )
                {
                    ip = 0;
                }
                if( ip==1 )
                {
                    ip = -1;
                }
            }
        }
    }


    /*************************************************************************
    DLALN2 solves a system of the form  (ca A - w D ) X = s B
    or (ca A' - w D) X = s B   with possible scaling ("s") and
    perturbation of A.  (A' means A-transpose.)

    A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
    real diagonal matrix, w is a real or complex value, and X and B are
    NA x 1 matrices -- real if w is real, complex if w is complex.  NA
    may be 1 or 2.

    If w is complex, X and B are represented as NA x 2 matrices,
    the first column of each being the real part and the second
    being the imaginary part.

    "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
    so chosen that X can be computed without overflow.  X is further
    scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
    than overflow.

    If both singular values of (ca A - w D) are less than SMIN,
    SMIN*identity will be used instead of (ca A - w D).  If only one
    singular value is less than SMIN, one element of (ca A - w D) will be
    perturbed enough to make the smallest singular value roughly SMIN.
    If both singular values are at least SMIN, (ca A - w D) will not be
    perturbed.  In any case, the perturbation will be at most some small
    multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
    are computed by infinity-norm approximations, and thus will only be
    correct to a factor of 2 or so.

    Note: all input quantities are assumed to be smaller than overflow
    by a reasonable factor.  (See BIGNUM.)

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         October 31, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void internalhsevdlaln2(const bool& ltrans,
        const int& na,
        const int& nw,
        const amp::ampf<Precision>& smin,
        const amp::ampf<Precision>& ca,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        const amp::ampf<Precision>& d1,
        const amp::ampf<Precision>& d2,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        const amp::ampf<Precision>& wr,
        const amp::ampf<Precision>& wi,
        ap::template_1d_array< bool >& rswap4,
        ap::template_1d_array< bool >& zswap4,
        ap::template_2d_array< int >& ipivot44,
        ap::template_1d_array< amp::ampf<Precision> >& civ4,
        ap::template_1d_array< amp::ampf<Precision> >& crv4,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& scl,
        amp::ampf<Precision>& xnorm,
        int& info)
    {
        int icmax;
        int j;
        amp::ampf<Precision> bbnd;
        amp::ampf<Precision> bi1;
        amp::ampf<Precision> bi2;
        amp::ampf<Precision> bignum;
        amp::ampf<Precision> bnorm;
        amp::ampf<Precision> br1;
        amp::ampf<Precision> br2;
        amp::ampf<Precision> ci21;
        amp::ampf<Precision> ci22;
        amp::ampf<Precision> cmax;
        amp::ampf<Precision> cnorm;
        amp::ampf<Precision> cr21;
        amp::ampf<Precision> cr22;
        amp::ampf<Precision> csi;
        amp::ampf<Precision> csr;
        amp::ampf<Precision> li21;
        amp::ampf<Precision> lr21;
        amp::ampf<Precision> smini;
        amp::ampf<Precision> smlnum;
        amp::ampf<Precision> temp;
        amp::ampf<Precision> u22abs;
        amp::ampf<Precision> ui11;
        amp::ampf<Precision> ui11r;
        amp::ampf<Precision> ui12;
        amp::ampf<Precision> ui12s;
        amp::ampf<Precision> ui22;
        amp::ampf<Precision> ur11;
        amp::ampf<Precision> ur11r;
        amp::ampf<Precision> ur12;
        amp::ampf<Precision> ur12s;
        amp::ampf<Precision> ur22;
        amp::ampf<Precision> xi1;
        amp::ampf<Precision> xi2;
        amp::ampf<Precision> xr1;
        amp::ampf<Precision> xr2;
        amp::ampf<Precision> tmp1;
        amp::ampf<Precision> tmp2;


        zswap4(1) = false;
        zswap4(2) = false;
        zswap4(3) = true;
        zswap4(4) = true;
        rswap4(1) = false;
        rswap4(2) = true;
        rswap4(3) = false;
        rswap4(4) = true;
        ipivot44(1,1) = 1;
        ipivot44(2,1) = 2;
        ipivot44(3,1) = 3;
        ipivot44(4,1) = 4;
        ipivot44(1,2) = 2;
        ipivot44(2,2) = 1;
        ipivot44(3,2) = 4;
        ipivot44(4,2) = 3;
        ipivot44(1,3) = 3;
        ipivot44(2,3) = 4;
        ipivot44(3,3) = 1;
        ipivot44(4,3) = 2;
        ipivot44(1,4) = 4;
        ipivot44(2,4) = 3;
        ipivot44(3,4) = 2;
        ipivot44(4,4) = 1;
        smlnum = 2*amp::ampf<Precision>::getAlgoPascalMinNumber();
        bignum = 1/smlnum;
        smini = amp::maximum<Precision>(smin, smlnum);
        
        //
        // Don't check for input errors
        //
        info = 0;
        
        //
        // Standard Initializations
        //
        scl = 1;
        if( na==1 )
        {
            
            //
            // 1 x 1  (i.e., scalar) system   C X = B
            //
            if( nw==1 )
            {
                
                //
                // Real 1x1 system.
                //
                // C = ca A - w D
                //
                csr = ca*a(1,1)-wr*d1;
                cnorm = amp::abs<Precision>(csr);
                
                //
                // If | C | < SMINI, use C = SMINI
                //
                if( cnorm<smini )
                {
                    csr = smini;
                    cnorm = smini;
                    info = 1;
                }
                
                //
                // Check scaling for  X = B / C
                //
                bnorm = amp::abs<Precision>(b(1,1));
                if( cnorm<1 && bnorm>1 )
                {
                    if( bnorm>bignum*cnorm )
                    {
                        scl = 1/bnorm;
                    }
                }
                
                //
                // Compute X
                //
                x(1,1) = b(1,1)*scl/csr;
                xnorm = amp::abs<Precision>(x(1,1));
            }
            else
            {
                
                //
                // Complex 1x1 system (w is complex)
                //
                // C = ca A - w D
                //
                csr = ca*a(1,1)-wr*d1;
                csi = -wi*d1;
                cnorm = amp::abs<Precision>(csr)+amp::abs<Precision>(csi);
                
                //
                // If | C | < SMINI, use C = SMINI
                //
                if( cnorm<smini )
                {
                    csr = smini;
                    csi = 0;
                    cnorm = smini;
                    info = 1;
                }
                
                //
                // Check scaling for  X = B / C
                //
                bnorm = amp::abs<Precision>(b(1,1))+amp::abs<Precision>(b(1,2));
                if( cnorm<1 && bnorm>1 )
                {
                    if( bnorm>bignum*cnorm )
                    {
                        scl = 1/bnorm;
                    }
                }
                
                //
                // Compute X
                //
                internalhsevdladiv<Precision>(scl*b(1,1), scl*b(1,2), csr, csi, tmp1, tmp2);
                x(1,1) = tmp1;
                x(1,2) = tmp2;
                xnorm = amp::abs<Precision>(x(1,1))+amp::abs<Precision>(x(1,2));
            }
        }
        else
        {
            
            //
            // 2x2 System
            //
            // Compute the real part of  C = ca A - w D  (or  ca A' - w D )
            //
            crv4(1+0) = ca*a(1,1)-wr*d1;
            crv4(2+2) = ca*a(2,2)-wr*d2;
            if( ltrans )
            {
                crv4(1+2) = ca*a(2,1);
                crv4(2+0) = ca*a(1,2);
            }
            else
            {
                crv4(2+0) = ca*a(2,1);
                crv4(1+2) = ca*a(1,2);
            }
            if( nw==1 )
            {
                
                //
                // Real 2x2 system  (w is real)
                //
                // Find the largest element in C
                //
                cmax = 0;
                icmax = 0;
                for(j=1; j<=4; j++)
                {
                    if( amp::abs<Precision>(crv4(j))>cmax )
                    {
                        cmax = amp::abs<Precision>(crv4(j));
                        icmax = j;
                    }
                }
                
                //
                // If norm(C) < SMINI, use SMINI*identity.
                //
                if( cmax<smini )
                {
                    bnorm = amp::maximum<Precision>(amp::abs<Precision>(b(1,1)), amp::abs<Precision>(b(2,1)));
                    if( smini<1 && bnorm>1 )
                    {
                        if( bnorm>bignum*smini )
                        {
                            scl = 1/bnorm;
                        }
                    }
                    temp = scl/smini;
                    x(1,1) = temp*b(1,1);
                    x(2,1) = temp*b(2,1);
                    xnorm = temp*bnorm;
                    info = 1;
                    return;
                }
                
                //
                // Gaussian elimination with complete pivoting.
                //
                ur11 = crv4(icmax);
                cr21 = crv4(ipivot44(2,icmax));
                ur12 = crv4(ipivot44(3,icmax));
                cr22 = crv4(ipivot44(4,icmax));
                ur11r = 1/ur11;
                lr21 = ur11r*cr21;
                ur22 = cr22-ur12*lr21;
                
                //
                // If smaller pivot < SMINI, use SMINI
                //
                if( amp::abs<Precision>(ur22)<smini )
                {
                    ur22 = smini;
                    info = 1;
                }
                if( rswap4(icmax) )
                {
                    br1 = b(2,1);
                    br2 = b(1,1);
                }
                else
                {
                    br1 = b(1,1);
                    br2 = b(2,1);
                }
                br2 = br2-lr21*br1;
                bbnd = amp::maximum<Precision>(amp::abs<Precision>(br1*(ur22*ur11r)), amp::abs<Precision>(br2));
                if( bbnd>1 && amp::abs<Precision>(ur22)<1 )
                {
                    if( bbnd>=bignum*amp::abs<Precision>(ur22) )
                    {
                        scl = 1/bbnd;
                    }
                }
                xr2 = br2*scl/ur22;
                xr1 = scl*br1*ur11r-xr2*(ur11r*ur12);
                if( zswap4(icmax) )
                {
                    x(1,1) = xr2;
                    x(2,1) = xr1;
                }
                else
                {
                    x(1,1) = xr1;
                    x(2,1) = xr2;
                }
                xnorm = amp::maximum<Precision>(amp::abs<Precision>(xr1), amp::abs<Precision>(xr2));
                
                //
                // Further scaling if  norm(A) norm(X) > overflow
                //
                if( xnorm>1 && cmax>1 )
                {
                    if( xnorm>bignum/cmax )
                    {
                        temp = cmax/bignum;
                        x(1,1) = temp*x(1,1);
                        x(2,1) = temp*x(2,1);
                        xnorm = temp*xnorm;
                        scl = temp*scl;
                    }
                }
            }
            else
            {
                
                //
                // Complex 2x2 system  (w is complex)
                //
                // Find the largest element in C
                //
                civ4(1+0) = -wi*d1;
                civ4(2+0) = 0;
                civ4(1+2) = 0;
                civ4(2+2) = -wi*d2;
                cmax = 0;
                icmax = 0;
                for(j=1; j<=4; j++)
                {
                    if( amp::abs<Precision>(crv4(j))+amp::abs<Precision>(civ4(j))>cmax )
                    {
                        cmax = amp::abs<Precision>(crv4(j))+amp::abs<Precision>(civ4(j));
                        icmax = j;
                    }
                }
                
                //
                // If norm(C) < SMINI, use SMINI*identity.
                //
                if( cmax<smini )
                {
                    bnorm = amp::maximum<Precision>(amp::abs<Precision>(b(1,1))+amp::abs<Precision>(b(1,2)), amp::abs<Precision>(b(2,1))+amp::abs<Precision>(b(2,2)));
                    if( smini<1 && bnorm>1 )
                    {
                        if( bnorm>bignum*smini )
                        {
                            scl = 1/bnorm;
                        }
                    }
                    temp = scl/smini;
                    x(1,1) = temp*b(1,1);
                    x(2,1) = temp*b(2,1);
                    x(1,2) = temp*b(1,2);
                    x(2,2) = temp*b(2,2);
                    xnorm = temp*bnorm;
                    info = 1;
                    return;
                }
                
                //
                // Gaussian elimination with complete pivoting.
                //
                ur11 = crv4(icmax);
                ui11 = civ4(icmax);
                cr21 = crv4(ipivot44(2,icmax));
                ci21 = civ4(ipivot44(2,icmax));
                ur12 = crv4(ipivot44(3,icmax));
                ui12 = civ4(ipivot44(3,icmax));
                cr22 = crv4(ipivot44(4,icmax));
                ci22 = civ4(ipivot44(4,icmax));
                if( icmax==1 || icmax==4 )
                {
                    
                    //
                    // Code when off-diagonals of pivoted C are real
                    //
                    if( amp::abs<Precision>(ur11)>amp::abs<Precision>(ui11) )
                    {
                        temp = ui11/ur11;
                        ur11r = 1/(ur11*(1+amp::sqr<Precision>(temp)));
                        ui11r = -temp*ur11r;
                    }
                    else
                    {
                        temp = ur11/ui11;
                        ui11r = -1/(ui11*(1+amp::sqr<Precision>(temp)));
                        ur11r = -temp*ui11r;
                    }
                    lr21 = cr21*ur11r;
                    li21 = cr21*ui11r;
                    ur12s = ur12*ur11r;
                    ui12s = ur12*ui11r;
                    ur22 = cr22-ur12*lr21;
                    ui22 = ci22-ur12*li21;
                }
                else
                {
                    
                    //
                    // Code when diagonals of pivoted C are real
                    //
                    ur11r = 1/ur11;
                    ui11r = 0;
                    lr21 = cr21*ur11r;
                    li21 = ci21*ur11r;
                    ur12s = ur12*ur11r;
                    ui12s = ui12*ur11r;
                    ur22 = cr22-ur12*lr21+ui12*li21;
                    ui22 = -ur12*li21-ui12*lr21;
                }
                u22abs = amp::abs<Precision>(ur22)+amp::abs<Precision>(ui22);
                
                //
                // If smaller pivot < SMINI, use SMINI
                //
                if( u22abs<smini )
                {
                    ur22 = smini;
                    ui22 = 0;
                    info = 1;
                }
                if( rswap4(icmax) )
                {
                    br2 = b(1,1);
                    br1 = b(2,1);
                    bi2 = b(1,2);
                    bi1 = b(2,2);
                }
                else
                {
                    br1 = b(1,1);
                    br2 = b(2,1);
                    bi1 = b(1,2);
                    bi2 = b(2,2);
                }
                br2 = br2-lr21*br1+li21*bi1;
                bi2 = bi2-li21*br1-lr21*bi1;
                bbnd = amp::maximum<Precision>((amp::abs<Precision>(br1)+amp::abs<Precision>(bi1))*(u22abs*(amp::abs<Precision>(ur11r)+amp::abs<Precision>(ui11r))), amp::abs<Precision>(br2)+amp::abs<Precision>(bi2));
                if( bbnd>1 && u22abs<1 )
                {
                    if( bbnd>=bignum*u22abs )
                    {
                        scl = 1/bbnd;
                        br1 = scl*br1;
                        bi1 = scl*bi1;
                        br2 = scl*br2;
                        bi2 = scl*bi2;
                    }
                }
                internalhsevdladiv<Precision>(br2, bi2, ur22, ui22, xr2, xi2);
                xr1 = ur11r*br1-ui11r*bi1-ur12s*xr2+ui12s*xi2;
                xi1 = ui11r*br1+ur11r*bi1-ui12s*xr2-ur12s*xi2;
                if( zswap4(icmax) )
                {
                    x(1,1) = xr2;
                    x(2,1) = xr1;
                    x(1,2) = xi2;
                    x(2,2) = xi1;
                }
                else
                {
                    x(1,1) = xr1;
                    x(2,1) = xr2;
                    x(1,2) = xi1;
                    x(2,2) = xi2;
                }
                xnorm = amp::maximum<Precision>(amp::abs<Precision>(xr1)+amp::abs<Precision>(xi1), amp::abs<Precision>(xr2)+amp::abs<Precision>(xi2));
                
                //
                // Further scaling if  norm(A) norm(X) > overflow
                //
                if( xnorm>1 && cmax>1 )
                {
                    if( xnorm>bignum/cmax )
                    {
                        temp = cmax/bignum;
                        x(1,1) = temp*x(1,1);
                        x(2,1) = temp*x(2,1);
                        x(1,2) = temp*x(1,2);
                        x(2,2) = temp*x(2,2);
                        xnorm = temp*xnorm;
                        scl = temp*scl;
                    }
                }
            }
        }
    }


    /*************************************************************************
    performs complex division in  real arithmetic

                            a + i*b
                 p + i*q = ---------
                            c + i*d

    The algorithm is due to Robert L. Smith and can be found
    in D. Knuth, The art of Computer Programming, Vol.2, p.195

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         October 31, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void internalhsevdladiv(const amp::ampf<Precision>& a,
        const amp::ampf<Precision>& b,
        const amp::ampf<Precision>& c,
        const amp::ampf<Precision>& d,
        amp::ampf<Precision>& p,
        amp::ampf<Precision>& q)
    {
        amp::ampf<Precision> e;
        amp::ampf<Precision> f;


        if( amp::abs<Precision>(d)<amp::abs<Precision>(c) )
        {
            e = d/c;
            f = c+d*e;
            p = (a+b*e)/f;
            q = (b-a*e)/f;
        }
        else
        {
            e = c/d;
            f = d+c*e;
            p = (b+a*e)/f;
            q = (-a+b*e)/f;
        }
    }


    template<unsigned int Precision>
    bool nonsymmetricevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        int vneeded,
        ap::template_1d_array< amp::ampf<Precision> >& wr,
        ap::template_1d_array< amp::ampf<Precision> >& wi,
        ap::template_2d_array< amp::ampf<Precision> >& vl,
        ap::template_2d_array< amp::ampf<Precision> >& vr)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > s;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< bool > sel;
        int i;
        int info;
        int m;


        ap::ap_error::make_assertion(vneeded>=0 && vneeded<=3);
        if( vneeded==0 )
        {
            
            //
            // Eigen values only
            //
            toupperhessenberg<Precision>(a, n, tau);
            hsschur::internalschurdecomposition<Precision>(a, n, 0, 0, wr, wi, s, info);
            result = info==0;
            return result;
        }
        
        //
        // Eigen values and vectors
        //
        toupperhessenberg<Precision>(a, n, tau);
        unpackqfromupperhessenberg<Precision>(a, n, tau, s);
        hsschur::internalschurdecomposition<Precision>(a, n, 1, 1, wr, wi, s, info);
        result = info==0;
        if( !result )
        {
            return result;
        }
        if( vneeded==1 || vneeded==3 )
        {
            vr.setbounds(1, n, 1, n);
            for(i=1; i<=n; i++)
            {
                amp::vmove(vr.getrow(i, 1, n), s.getrow(i, 1, n));
            }
        }
        if( vneeded==2 || vneeded==3 )
        {
            vl.setbounds(1, n, 1, n);
            for(i=1; i<=n; i++)
            {
                amp::vmove(vl.getrow(i, 1, n), s.getrow(i, 1, n));
            }
        }
        internaltrevc<Precision>(a, n, vneeded, 1, sel, vl, vr, m, info);
        result = info==0;
        return result;
    }


    template<unsigned int Precision>
    void toupperhessenberg(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& tau)
    {
        int i;
        int ip1;
        int nmi;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::ampf<Precision> > work;


        ap::ap_error::make_assertion(n>=0);
        
        //
        // Quick return if possible
        //
        if( n<=1 )
        {
            return;
        }
        tau.setbounds(1, n-1);
        t.setbounds(1, n);
        work.setbounds(1, n);
        for(i=1; i<=n-1; i++)
        {
            
            //
            // Compute elementary reflector H(i) to annihilate A(i+2:ihi,i)
            //
            ip1 = i+1;
            nmi = n-i;
            amp::vmove(t.getvector(1, nmi), a.getcolumn(i, ip1, n));
            reflections::generatereflection<Precision>(t, nmi, v);
            amp::vmove(a.getcolumn(i, ip1, n), t.getvector(1, nmi));
            tau(i) = v;
            t(1) = 1;
            
            //
            // Apply H(i) to A(1:ihi,i+1:ihi) from the right
            //
            reflections::applyreflectionfromtheright<Precision>(a, v, t, 1, n, i+1, n, work);
            
            //
            // Apply H(i) to A(i+1:ihi,i+1:n) from the left
            //
            reflections::applyreflectionfromtheleft<Precision>(a, v, t, i+1, n, i+1, n, work);
        }
    }


    template<unsigned int Precision>
    void unpackqfromupperhessenberg(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& tau,
        ap::template_2d_array< amp::ampf<Precision> >& q)
    {
        int i;
        int j;
        ap::template_1d_array< amp::ampf<Precision> > v;
        ap::template_1d_array< amp::ampf<Precision> > work;
        int ip1;
        int nmi;


        if( n==0 )
        {
            return;
        }
        
        //
        // init
        //
        q.setbounds(1, n, 1, n);
        v.setbounds(1, n);
        work.setbounds(1, n);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                if( i==j )
                {
                    q(i,j) = 1;
                }
                else
                {
                    q(i,j) = 0;
                }
            }
        }
        
        //
        // unpack Q
        //
        for(i=1; i<=n-1; i++)
        {
            
            //
            // Apply H(i)
            //
            ip1 = i+1;
            nmi = n-i;
            amp::vmove(v.getvector(1, nmi), a.getcolumn(i, ip1, n));
            v(1) = 1;
            reflections::applyreflectionfromtheright<Precision>(q, tau(i), v, 1, n, i+1, n, work);
        }
    }


    template<unsigned int Precision>
    void unpackhfromupperhessenberg(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& tau,
        ap::template_2d_array< amp::ampf<Precision> >& h)
    {
        int i;
        int j;
        ap::template_1d_array< amp::ampf<Precision> > v;
        ap::template_1d_array< amp::ampf<Precision> > work;


        if( n==0 )
        {
            return;
        }
        h.setbounds(1, n, 1, n);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=i-2; j++)
            {
                h(i,j) = 0;
            }
            j = ap::maxint(1, i-1);
            amp::vmove(h.getrow(i, j, n), a.getrow(i, j, n));
        }
    }
} // namespace

#endif
