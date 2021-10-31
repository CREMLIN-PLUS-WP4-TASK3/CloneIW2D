/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee.  All rights reserved.

Contributors:
    * Sergey Bochkanov (ALGLIB project). Translation from FORTRAN to
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

#ifndef _rcond_h
#define _rcond_h

#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
namespace rcond
{
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixrcond1(ap::template_2d_array< amp::ampf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixrcondinf(ap::template_2d_array< amp::ampf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixrcond(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixtrrcond1(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixtrrcondinf(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit);
    template<unsigned int Precision>
    amp::ampf<Precision> hpdmatrixrcond(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixrcond1(ap::template_2d_array< amp::campf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixrcondinf(ap::template_2d_array< amp::campf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixlurcond1(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixlurcondinf(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixcholeskyrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    amp::ampf<Precision> hpdmatrixcholeskyrcond(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixlurcond1(const ap::template_2d_array< amp::campf<Precision> >& lua,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixlurcondinf(const ap::template_2d_array< amp::campf<Precision> >& lua,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixtrrcond1(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit);
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixtrrcondinf(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit);
    template<unsigned int Precision>
    amp::ampf<Precision> rcondthreshold();
    template<unsigned int Precision>
    void rmatrixrcondtrinternal(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        bool onenorm,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void cmatrixrcondtrinternal(const ap::template_2d_array< amp::campf<Precision> >& a,
        const int& n,
        bool isupper,
        bool isunit,
        bool onenorm,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void spdmatrixrcondcholeskyinternal(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        bool isnormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void hpdmatrixrcondcholeskyinternal(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        bool isnormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void rmatrixrcondluinternal(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n,
        bool onenorm,
        bool isanormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void cmatrixrcondluinternal(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const int& n,
        bool onenorm,
        bool isanormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc);
    template<unsigned int Precision>
    void rmatrixestimatenorm(int n,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< int >& isgn,
        amp::ampf<Precision>& est,
        int& kase);
    template<unsigned int Precision>
    void cmatrixestimatenorm(const int& n,
        ap::template_1d_array< amp::campf<Precision> >& v,
        ap::template_1d_array< amp::campf<Precision> >& x,
        amp::ampf<Precision>& est,
        int& kase,
        ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave);
    template<unsigned int Precision>
    amp::ampf<Precision> internalcomplexrcondscsum1(const ap::template_1d_array< amp::campf<Precision> >& x,
        int n);
    template<unsigned int Precision>
    int internalcomplexrcondicmax1(const ap::template_1d_array< amp::campf<Precision> >& x,
        int n);
    template<unsigned int Precision>
    void internalcomplexrcondsaveall(ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave,
        int& i,
        int& iter,
        int& j,
        int& jlast,
        int& jump,
        amp::ampf<Precision>& absxi,
        amp::ampf<Precision>& altsgn,
        amp::ampf<Precision>& estold,
        amp::ampf<Precision>& temp);
    template<unsigned int Precision>
    void internalcomplexrcondloadall(ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave,
        int& i,
        int& iter,
        int& j,
        int& jlast,
        int& jump,
        amp::ampf<Precision>& absxi,
        amp::ampf<Precision>& altsgn,
        amp::ampf<Precision>& estold,
        amp::ampf<Precision>& temp);


    /*************************************************************************
    Estimate of a matrix condition number (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixrcond1(ap::template_2d_array< amp::ampf<Precision> > a,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        ap::template_1d_array< amp::ampf<Precision> > t;


        ap::ap_error::make_assertion(n>=1);
        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(j) = t(j)+amp::abs<Precision>(a(i,j));
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        trfac::rmatrixlu<Precision>(a, n, n, pivots);
        rmatrixrcondluinternal<Precision>(a, n, true, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of a matrix condition number (infinity-norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixrcondinf(ap::template_2d_array< amp::ampf<Precision> > a,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;


        ap::ap_error::make_assertion(n>=1);
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=n-1; j++)
            {
                v = v+amp::abs<Precision>(a(i,j));
            }
            nrm = amp::maximum<Precision>(nrm, v);
        }
        trfac::rmatrixlu<Precision>(a, n, n, pivots);
        rmatrixrcondluinternal<Precision>(a, n, false, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Condition number estimate of a symmetric positive definite matrix.

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    It should be noted that 1-norm and inf-norm of condition numbers of symmetric
    matrices are equal, so the algorithm doesn't take into account the
    differences between these types of norms.

    Input parameters:
        A       -   symmetric positive definite matrix which is given by its
                    upper or lower triangle depending on the value of
                    IsUpper. Array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.

    Result:
        1/LowerBound(cond(A)), if matrix A is positive definite,
       -1, if matrix A is not positive definite, and its condition number
        could not be found by this algorithm.

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixrcond(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        int j1;
        int j2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< amp::ampf<Precision> > t;


        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            for(j=j1; j<=j2; j++)
            {
                if( i==j )
                {
                    t(i) = t(i)+amp::abs<Precision>(a(i,i));
                }
                else
                {
                    t(i) = t(i)+amp::abs<Precision>(a(i,j));
                    t(j) = t(j)+amp::abs<Precision>(a(i,j));
                }
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        if( trfac::spdmatrixcholesky<Precision>(a, n, isupper) )
        {
            spdmatrixrcondcholeskyinternal<Precision>(a, n, isupper, true, nrm, v);
            result = v;
        }
        else
        {
            result = -1;
        }
        return result;
    }


    /*************************************************************************
    Triangular matrix: estimate of a condition number (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A       -   matrix. Array[0..N-1, 0..N-1].
        N       -   size of A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixtrrcond1(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        ap::template_1d_array< amp::ampf<Precision> > t;
        int j1;
        int j2;


        ap::ap_error::make_assertion(n>=1);
        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            for(j=j1; j<=j2; j++)
            {
                t(j) = t(j)+amp::abs<Precision>(a(i,j));
            }
            if( isunit )
            {
                t(i) = t(i)+1;
            }
            else
            {
                t(i) = t(i)+amp::abs<Precision>(a(i,i));
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        rmatrixrcondtrinternal<Precision>(a, n, isupper, isunit, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Triangular matrix: estimate of a matrix condition number (infinity-norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixtrrcondinf(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        int j1;
        int j2;


        ap::ap_error::make_assertion(n>=1);
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            v = 0;
            for(j=j1; j<=j2; j++)
            {
                v = v+amp::abs<Precision>(a(i,j));
            }
            if( isunit )
            {
                v = v+1;
            }
            else
            {
                v = v+amp::abs<Precision>(a(i,i));
            }
            nrm = amp::maximum<Precision>(nrm, v);
        }
        rmatrixrcondtrinternal<Precision>(a, n, isupper, isunit, false, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Condition number estimate of a Hermitian positive definite matrix.

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    It should be noted that 1-norm and inf-norm of condition numbers of symmetric
    matrices are equal, so the algorithm doesn't take into account the
    differences between these types of norms.

    Input parameters:
        A       -   Hermitian positive definite matrix which is given by its
                    upper or lower triangle depending on the value of
                    IsUpper. Array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.

    Result:
        1/LowerBound(cond(A)), if matrix A is positive definite,
       -1, if matrix A is not positive definite, and its condition number
        could not be found by this algorithm.

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> hpdmatrixrcond(ap::template_2d_array< amp::campf<Precision> > a,
        int n,
        bool isupper)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        int j1;
        int j2;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< amp::ampf<Precision> > t;


        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i;
            }
            for(j=j1; j<=j2; j++)
            {
                if( i==j )
                {
                    t(i) = t(i)+amp::abscomplex<Precision>(a(i,i));
                }
                else
                {
                    t(i) = t(i)+amp::abscomplex<Precision>(a(i,j));
                    t(j) = t(j)+amp::abscomplex<Precision>(a(i,j));
                }
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        if( trfac::hpdmatrixcholesky<Precision>(a, n, isupper) )
        {
            hpdmatrixrcondcholeskyinternal<Precision>(a, n, isupper, true, nrm, v);
            result = v;
        }
        else
        {
            result = -1;
        }
        return result;
    }


    /*************************************************************************
    Estimate of a matrix condition number (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixrcond1(ap::template_2d_array< amp::campf<Precision> > a,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        ap::template_1d_array< amp::ampf<Precision> > t;


        ap::ap_error::make_assertion(n>=1);
        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                t(j) = t(j)+amp::abscomplex<Precision>(a(i,j));
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        trfac::cmatrixlu<Precision>(a, n, n, pivots);
        cmatrixrcondluinternal<Precision>(a, n, true, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of a matrix condition number (infinity-norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixrcondinf(ap::template_2d_array< amp::campf<Precision> > a,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;


        ap::ap_error::make_assertion(n>=1);
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=n-1; j++)
            {
                v = v+amp::abscomplex<Precision>(a(i,j));
            }
            nrm = amp::maximum<Precision>(nrm, v);
        }
        trfac::cmatrixlu<Precision>(a, n, n, pivots);
        cmatrixrcondluinternal<Precision>(a, n, false, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        LUA         -   LU decomposition of a matrix in compact form. Output of
                        the RMatrixLU subroutine.
        N           -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixlurcond1(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        rmatrixrcondluinternal<Precision>(lua, n, true, false, amp::ampf<Precision>(0), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of the condition number of a matrix given by its LU decomposition
    (infinity norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        LUA     -   LU decomposition of a matrix in compact form. Output of
                    the RMatrixLU subroutine.
        N       -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixlurcondinf(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        rmatrixrcondluinternal<Precision>(lua, n, false, false, amp::ampf<Precision>(0), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Condition number estimate of a symmetric positive definite matrix given by
    Cholesky decomposition.

    The algorithm calculates a lower bound of the condition number. In this
    case, the algorithm does not return a lower bound of the condition number,
    but an inverse number (to avoid an overflow in case of a singular matrix).

    It should be noted that 1-norm and inf-norm condition numbers of symmetric
    matrices are equal, so the algorithm doesn't take into account the
    differences between these types of norms.

    Input parameters:
        CD  - Cholesky decomposition of matrix A,
              output of SMatrixCholesky subroutine.
        N   - size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixcholeskyrcond(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        spdmatrixrcondcholeskyinternal<Precision>(a, n, isupper, false, amp::ampf<Precision>(0), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Condition number estimate of a Hermitian positive definite matrix given by
    Cholesky decomposition.

    The algorithm calculates a lower bound of the condition number. In this
    case, the algorithm does not return a lower bound of the condition number,
    but an inverse number (to avoid an overflow in case of a singular matrix).

    It should be noted that 1-norm and inf-norm condition numbers of symmetric
    matrices are equal, so the algorithm doesn't take into account the
    differences between these types of norms.

    Input parameters:
        CD  - Cholesky decomposition of matrix A,
              output of SMatrixCholesky subroutine.
        N   - size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> hpdmatrixcholeskyrcond(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        hpdmatrixrcondcholeskyinternal<Precision>(a, n, isupper, false, amp::ampf<Precision>(0), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of the condition number of a matrix given by its LU decomposition (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        LUA         -   LU decomposition of a matrix in compact form. Output of
                        the CMatrixLU subroutine.
        N           -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixlurcond1(const ap::template_2d_array< amp::campf<Precision> >& lua,
        int n)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>=1);
        cmatrixrcondluinternal<Precision>(lua, n, true, false, amp::ampf<Precision>("0.0"), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Estimate of the condition number of a matrix given by its LU decomposition
    (infinity norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        LUA     -   LU decomposition of a matrix in compact form. Output of
                    the CMatrixLU subroutine.
        N       -   size of matrix A.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixlurcondinf(const ap::template_2d_array< amp::campf<Precision> >& lua,
        int n)
    {
        amp::ampf<Precision> result;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>=1);
        cmatrixrcondluinternal<Precision>(lua, n, false, false, amp::ampf<Precision>("0.0"), v);
        result = v;
        return result;
    }


    /*************************************************************************
    Triangular matrix: estimate of a condition number (1-norm)

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A       -   matrix. Array[0..N-1, 0..N-1].
        N       -   size of A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixtrrcond1(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        ap::template_1d_array< amp::ampf<Precision> > t;
        int j1;
        int j2;


        ap::ap_error::make_assertion(n>=1);
        t.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            t(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            for(j=j1; j<=j2; j++)
            {
                t(j) = t(j)+amp::abscomplex<Precision>(a(i,j));
            }
            if( isunit )
            {
                t(i) = t(i)+1;
            }
            else
            {
                t(i) = t(i)+amp::abscomplex<Precision>(a(i,i));
            }
        }
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            nrm = amp::maximum<Precision>(nrm, t(i));
        }
        cmatrixrcondtrinternal<Precision>(a, n, isupper, isunit, true, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Triangular matrix: estimate of a matrix condition number (infinity-norm).

    The algorithm calculates a lower bound of the condition number. In this case,
    the algorithm does not return a lower bound of the condition number, but an
    inverse number (to avoid an overflow in case of a singular matrix).

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Result: 1/LowerBound(cond(A))

    NOTE:
        if k(A) is very large, then matrix is  assumed  degenerate,  k(A)=INF,
        0.0 is returned in such cases.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> cmatrixtrrcondinf(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit)
    {
        amp::ampf<Precision> result;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> nrm;
        ap::template_1d_array< int > pivots;
        int j1;
        int j2;


        ap::ap_error::make_assertion(n>=1);
        nrm = 0;
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            v = 0;
            for(j=j1; j<=j2; j++)
            {
                v = v+amp::abscomplex<Precision>(a(i,j));
            }
            if( isunit )
            {
                v = v+1;
            }
            else
            {
                v = v+amp::abscomplex<Precision>(a(i,i));
            }
            nrm = amp::maximum<Precision>(nrm, v);
        }
        cmatrixrcondtrinternal<Precision>(a, n, isupper, isunit, false, nrm, v);
        result = v;
        return result;
    }


    /*************************************************************************
    Threshold for rcond: matrices with condition number beyond this  threshold
    are considered singular.

    Threshold must be far enough from underflow, at least Sqr(Threshold)  must
    be greater than underflow.
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rcondthreshold()
    {
        amp::ampf<Precision> result;


        result = amp::sqrt<Precision>(amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber()));
        return result;
    }


    /*************************************************************************
    Internal subroutine for condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrcondtrinternal(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        bool onenorm,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        ap::template_1d_array< amp::ampf<Precision> > ex;
        ap::template_1d_array< amp::ampf<Precision> > ev;
        ap::template_1d_array< int > iwork;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;
        int j;
        int kase;
        int kase1;
        int j1;
        int j2;
        amp::ampf<Precision> ainvnm;
        amp::ampf<Precision> maxgrowth;
        amp::ampf<Precision> s;
        bool mupper;
        bool mtrans;
        bool munit;


        
        //
        // RC=0 if something happens
        //
        rc = 0;
        
        //
        // init
        //
        if( onenorm )
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        mupper = true;
        mtrans = true;
        munit = true;
        iwork.setlength(n+1);
        tmp.setlength(n);
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        s = 0;
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            for(j=j1; j<=j2; j++)
            {
                s = amp::maximum<Precision>(s, amp::abs<Precision>(a(i,j)));
            }
            if( isunit )
            {
                s = amp::maximum<Precision>(s, amp::ampf<Precision>(1));
            }
            else
            {
                s = amp::maximum<Precision>(s, amp::abs<Precision>(a(i,i)));
            }
        }
        if( s==0 )
        {
            s = 1;
        }
        s = 1/s;
        
        //
        // Scale according to S
        //
        anorm = anorm*s;
        
        //
        // Quick return if possible
        // We assume that ANORM<>0 after this block
        //
        if( anorm==0 )
        {
            return;
        }
        if( n==1 )
        {
            rc = 1;
            return;
        }
        
        //
        // Estimate the norm of inv(A).
        //
        ainvnm = 0;
        kase = 0;
        while( true )
        {
            rmatrixestimatenorm<Precision>(n, ev, ex, iwork, ainvnm, kase);
            if( kase==0 )
            {
                break;
            }
            
            //
            // from 1-based array to 0-based
            //
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            
            //
            // multiply by inv(A) or inv(A')
            //
            if( kase==kase1 )
            {
                
                //
                // multiply by inv(A)
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(a, s, n, ex, isupper, 0, isunit, maxgrowth) )
                {
                    return;
                }
            }
            else
            {
                
                //
                // multiply by inv(A')
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(a, s, n, ex, isupper, 1, isunit, maxgrowth) )
                {
                    return;
                }
            }
            
            //
            // from 0-based array to 1-based
            //
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            rc = 1/ainvnm;
            rc = rc/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         March 31, 1993
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrcondtrinternal(const ap::template_2d_array< amp::campf<Precision> >& a,
        const int& n,
        bool isupper,
        bool isunit,
        bool onenorm,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        ap::template_1d_array< amp::campf<Precision> > ex;
        ap::template_1d_array< amp::campf<Precision> > cwork2;
        ap::template_1d_array< amp::campf<Precision> > cwork3;
        ap::template_1d_array< amp::campf<Precision> > cwork4;
        ap::template_1d_array< int > isave;
        ap::template_1d_array< amp::ampf<Precision> > rsave;
        int kase;
        int kase1;
        amp::ampf<Precision> ainvnm;
        amp::campf<Precision> v;
        int i;
        int j;
        int j1;
        int j2;
        amp::ampf<Precision> s;
        amp::ampf<Precision> maxgrowth;


        
        //
        // RC=0 if something happens
        //
        rc = 0;
        
        //
        // init
        //
        if( n<=0 )
        {
            return;
        }
        if( n==0 )
        {
            rc = 1;
            return;
        }
        cwork2.setlength(n+1);
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        s = 0;
        for(i=0; i<=n-1; i++)
        {
            if( isupper )
            {
                j1 = i+1;
                j2 = n-1;
            }
            else
            {
                j1 = 0;
                j2 = i-1;
            }
            for(j=j1; j<=j2; j++)
            {
                s = amp::maximum<Precision>(s, amp::abscomplex<Precision>(a(i,j)));
            }
            if( isunit )
            {
                s = amp::maximum<Precision>(s, amp::ampf<Precision>(1));
            }
            else
            {
                s = amp::maximum<Precision>(s, amp::abscomplex<Precision>(a(i,i)));
            }
        }
        if( s==0 )
        {
            s = 1;
        }
        s = 1/s;
        
        //
        // Scale according to S
        //
        anorm = anorm*s;
        
        //
        // Quick return if possible
        //
        if( anorm==0 )
        {
            return;
        }
        
        //
        // Estimate the norm of inv(A).
        //
        ainvnm = 0;
        if( onenorm )
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        kase = 0;
        while( true )
        {
            cmatrixestimatenorm<Precision>(n, cwork4, ex, ainvnm, kase, isave, rsave);
            if( kase==0 )
            {
                break;
            }
            
            //
            // From 1-based to 0-based
            //
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            
            //
            // multiply by inv(A) or inv(A')
            //
            if( kase==kase1 )
            {
                
                //
                // multiply by inv(A)
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(a, s, n, ex, isupper, 0, isunit, maxgrowth) )
                {
                    return;
                }
            }
            else
            {
                
                //
                // multiply by inv(A')
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(a, s, n, ex, isupper, 2, isunit, maxgrowth) )
                {
                    return;
                }
            }
            
            //
            // from 0-based to 1-based
            //
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            rc = 1/ainvnm;
            rc = rc/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Internal subroutine for condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixrcondcholeskyinternal(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        bool isnormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        int i;
        int j;
        int kase;
        amp::ampf<Precision> ainvnm;
        ap::template_1d_array< amp::ampf<Precision> > ex;
        ap::template_1d_array< amp::ampf<Precision> > ev;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        ap::template_1d_array< int > iwork;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> v;
        amp::ampf<Precision> maxgrowth;


        ap::ap_error::make_assertion(n>=1);
        tmp.setlength(n);
        
        //
        // RC=0 if something happens
        //
        rc = 0;
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        sa = 0;
        if( isupper )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    sa = amp::maximum<Precision>(sa, amp::abscomplex<Precision>(cha(i,j)));
                }
            }
        }
        else
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    sa = amp::maximum<Precision>(sa, amp::abscomplex<Precision>(cha(i,j)));
                }
            }
        }
        if( sa==0 )
        {
            sa = 1;
        }
        sa = 1/sa;
        
        //
        // Estimate the norm of A.
        //
        if( !isnormprovided )
        {
            kase = 0;
            anorm = 0;
            while( true )
            {
                rmatrixestimatenorm<Precision>(n, ev, ex, iwork, anorm, kase);
                if( kase==0 )
                {
                    break;
                }
                if( isupper )
                {
                    
                    //
                    // Multiply by U
                    //
                    for(i=1; i<=n; i++)
                    {
                        v = amp::vdotproduct(cha.getrow(i-1, i-1, n-1), ex.getvector(i, n));
                        ex(i) = v;
                    }
                    amp::vmul(ex.getvector(1, n), sa);
                    
                    //
                    // Multiply by U'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        amp::vadd(tmp.getvector(i, n-1), cha.getrow(i, i, n-1), v);
                    }
                    amp::vmove(ex.getvector(1, n), tmp.getvector(0, n-1));
                    amp::vmul(ex.getvector(1, n), sa);
                }
                else
                {
                    
                    //
                    // Multiply by L'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        amp::vadd(tmp.getvector(0, i), cha.getrow(i, 0, i), v);
                    }
                    amp::vmove(ex.getvector(1, n), tmp.getvector(0, n-1));
                    amp::vmul(ex.getvector(1, n), sa);
                    
                    //
                    // Multiply by L
                    //
                    for(i=n; i>=1; i--)
                    {
                        v = amp::vdotproduct(cha.getrow(i-1, 0, i-1), ex.getvector(1, i));
                        ex(i) = v;
                    }
                    amp::vmul(ex.getvector(1, n), sa);
                }
            }
        }
        
        //
        // Quick return if possible
        //
        if( anorm==0 )
        {
            return;
        }
        if( n==1 )
        {
            rc = 1;
            return;
        }
        
        //
        // Estimate the 1-norm of inv(A).
        //
        kase = 0;
        while( true )
        {
            rmatrixestimatenorm<Precision>(n, ev, ex, iwork, ainvnm, kase);
            if( kase==0 )
            {
                break;
            }
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            if( isupper )
            {
                
                //
                // Multiply by inv(U').
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 1, false, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(U).
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
                {
                    return;
                }
            }
            else
            {
                
                //
                // Multiply by inv(L).
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(L').
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 1, false, maxgrowth) )
                {
                    return;
                }
            }
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            v = 1/ainvnm;
            rc = v/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Internal subroutine for condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixrcondcholeskyinternal(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        bool isnormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        ap::template_1d_array< int > isave;
        ap::template_1d_array< amp::ampf<Precision> > rsave;
        ap::template_1d_array< amp::campf<Precision> > ex;
        ap::template_1d_array< amp::campf<Precision> > ev;
        ap::template_1d_array< amp::campf<Precision> > tmp;
        int kase;
        amp::ampf<Precision> ainvnm;
        amp::campf<Precision> v;
        int i;
        int j;
        amp::ampf<Precision> sa;
        amp::ampf<Precision> maxgrowth;
        int i_;
        int i1_;


        ap::ap_error::make_assertion(n>=1);
        tmp.setlength(n);
        
        //
        // RC=0 if something happens
        //
        rc = 0;
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        sa = 0;
        if( isupper )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=i; j<=n-1; j++)
                {
                    sa = amp::maximum<Precision>(sa, amp::abscomplex<Precision>(cha(i,j)));
                }
            }
        }
        else
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i; j++)
                {
                    sa = amp::maximum<Precision>(sa, amp::abscomplex<Precision>(cha(i,j)));
                }
            }
        }
        if( sa==0 )
        {
            sa = 1;
        }
        sa = 1/sa;
        
        //
        // Estimate the norm of A
        //
        if( !isnormprovided )
        {
            anorm = 0;
            kase = 0;
            while( true )
            {
                cmatrixestimatenorm<Precision>(n, ev, ex, anorm, kase, isave, rsave);
                if( kase==0 )
                {
                    break;
                }
                if( isupper )
                {
                    
                    //
                    // Multiply by U
                    //
                    for(i=1; i<=n; i++)
                    {
                        i1_ = (i)-(i-1);
                        v = 0.0;
                        for(i_=i-1; i_<=n-1;i_++)
                        {
                            v += cha(i-1,i_)*ex(i_+i1_);
                        }
                        ex(i) = v;
                    }
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = sa*ex(i_);
                    }
                    
                    //
                    // Multiply by U'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        for(i_=i; i_<=n-1;i_++)
                        {
                            tmp(i_) = tmp(i_) + v*amp::conj(cha(i,i_));
                        }
                    }
                    i1_ = (0) - (1);
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = tmp(i_+i1_);
                    }
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = sa*ex(i_);
                    }
                }
                else
                {
                    
                    //
                    // Multiply by L'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        for(i_=0; i_<=i;i_++)
                        {
                            tmp(i_) = tmp(i_) + v*amp::conj(cha(i,i_));
                        }
                    }
                    i1_ = (0) - (1);
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = tmp(i_+i1_);
                    }
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = sa*ex(i_);
                    }
                    
                    //
                    // Multiply by L
                    //
                    for(i=n; i>=1; i--)
                    {
                        i1_ = (1)-(0);
                        v = 0.0;
                        for(i_=0; i_<=i-1;i_++)
                        {
                            v += cha(i-1,i_)*ex(i_+i1_);
                        }
                        ex(i) = v;
                    }
                    for(i_=1; i_<=n;i_++)
                    {
                        ex(i_) = sa*ex(i_);
                    }
                }
            }
        }
        
        //
        // Quick return if possible
        // After this block we assume that ANORM<>0
        //
        if( anorm==0 )
        {
            return;
        }
        if( n==1 )
        {
            rc = 1;
            return;
        }
        
        //
        // Estimate the norm of inv(A).
        //
        ainvnm = 0;
        kase = 0;
        while( true )
        {
            cmatrixestimatenorm<Precision>(n, ev, ex, ainvnm, kase, isave, rsave);
            if( kase==0 )
            {
                break;
            }
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            if( isupper )
            {
                
                //
                // Multiply by inv(U').
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 2, false, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(U).
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
                {
                    return;
                }
            }
            else
            {
                
                //
                // Multiply by inv(L).
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 0, false, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(L').
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(cha, sa, n, ex, isupper, 2, false, maxgrowth) )
                {
                    return;
                }
            }
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            rc = 1/ainvnm;
            rc = rc/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Internal subroutine for condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixrcondluinternal(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        int n,
        bool onenorm,
        bool isanormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        ap::template_1d_array< amp::ampf<Precision> > ex;
        ap::template_1d_array< amp::ampf<Precision> > ev;
        ap::template_1d_array< int > iwork;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;
        int j;
        int kase;
        int kase1;
        amp::ampf<Precision> ainvnm;
        amp::ampf<Precision> maxgrowth;
        amp::ampf<Precision> su;
        amp::ampf<Precision> sl;
        bool mupper;
        bool mtrans;
        bool munit;


        
        //
        // RC=0 if something happens
        //
        rc = 0;
        
        //
        // init
        //
        if( onenorm )
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        mupper = true;
        mtrans = true;
        munit = true;
        iwork.setlength(n+1);
        tmp.setlength(n);
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        su = 0;
        sl = 1;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                sl = amp::maximum<Precision>(sl, amp::abs<Precision>(lua(i,j)));
            }
            for(j=i; j<=n-1; j++)
            {
                su = amp::maximum<Precision>(su, amp::abs<Precision>(lua(i,j)));
            }
        }
        if( su==0 )
        {
            su = 1;
        }
        su = 1/su;
        sl = 1/sl;
        
        //
        // Estimate the norm of A.
        //
        if( !isanormprovided )
        {
            kase = 0;
            anorm = 0;
            while( true )
            {
                rmatrixestimatenorm<Precision>(n, ev, ex, iwork, anorm, kase);
                if( kase==0 )
                {
                    break;
                }
                if( kase==kase1 )
                {
                    
                    //
                    // Multiply by U
                    //
                    for(i=1; i<=n; i++)
                    {
                        v = amp::vdotproduct(lua.getrow(i-1, i-1, n-1), ex.getvector(i, n));
                        ex(i) = v;
                    }
                    
                    //
                    // Multiply by L
                    //
                    for(i=n; i>=1; i--)
                    {
                        if( i>1 )
                        {
                            v = amp::vdotproduct(lua.getrow(i-1, 0, i-2), ex.getvector(1, i-1));
                        }
                        else
                        {
                            v = 0;
                        }
                        ex(i) = ex(i)+v;
                    }
                }
                else
                {
                    
                    //
                    // Multiply by L'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        if( i>=1 )
                        {
                            amp::vadd(tmp.getvector(0, i-1), lua.getrow(i, 0, i-1), v);
                        }
                        tmp(i) = tmp(i)+v;
                    }
                    amp::vmove(ex.getvector(1, n), tmp.getvector(0, n-1));
                    
                    //
                    // Multiply by U'
                    //
                    for(i=0; i<=n-1; i++)
                    {
                        tmp(i) = 0;
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        v = ex(i+1);
                        amp::vadd(tmp.getvector(i, n-1), lua.getrow(i, i, n-1), v);
                    }
                    amp::vmove(ex.getvector(1, n), tmp.getvector(0, n-1));
                }
            }
        }
        
        //
        // Scale according to SU/SL
        //
        anorm = anorm*su*sl;
        
        //
        // Quick return if possible
        // We assume that ANORM<>0 after this block
        //
        if( anorm==0 )
        {
            return;
        }
        if( n==1 )
        {
            rc = 1;
            return;
        }
        
        //
        // Estimate the norm of inv(A).
        //
        ainvnm = 0;
        kase = 0;
        while( true )
        {
            rmatrixestimatenorm<Precision>(n, ev, ex, iwork, ainvnm, kase);
            if( kase==0 )
            {
                break;
            }
            
            //
            // from 1-based array to 0-based
            //
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            
            //
            // multiply by inv(A) or inv(A')
            //
            if( kase==kase1 )
            {
                
                //
                // Multiply by inv(L).
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(lua, sl, n, ex, !mupper, 0, munit, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(U).
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(lua, su, n, ex, mupper, 0, !munit, maxgrowth) )
                {
                    return;
                }
            }
            else
            {
                
                //
                // Multiply by inv(U').
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(lua, su, n, ex, mupper, 1, !munit, maxgrowth) )
                {
                    return;
                }
                
                //
                // Multiply by inv(L').
                //
                if( !safesolve::rmatrixscaledtrsafesolve<Precision>(lua, sl, n, ex, !mupper, 1, munit, maxgrowth) )
                {
                    return;
                }
            }
            
            //
            // from 0-based array to 1-based
            //
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            rc = 1/ainvnm;
            rc = rc/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Condition number estimation

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         March 31, 1993
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixrcondluinternal(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const int& n,
        bool onenorm,
        bool isanormprovided,
        amp::ampf<Precision> anorm,
        amp::ampf<Precision>& rc)
    {
        ap::template_1d_array< amp::campf<Precision> > ex;
        ap::template_1d_array< amp::campf<Precision> > cwork2;
        ap::template_1d_array< amp::campf<Precision> > cwork3;
        ap::template_1d_array< amp::campf<Precision> > cwork4;
        ap::template_1d_array< int > isave;
        ap::template_1d_array< amp::ampf<Precision> > rsave;
        int kase;
        int kase1;
        amp::ampf<Precision> ainvnm;
        amp::campf<Precision> v;
        int i;
        int j;
        amp::ampf<Precision> su;
        amp::ampf<Precision> sl;
        amp::ampf<Precision> maxgrowth;
        int i_;
        int i1_;


        if( n<=0 )
        {
            return;
        }
        cwork2.setlength(n+1);
        rc = 0;
        if( n==0 )
        {
            rc = 1;
            return;
        }
        
        //
        // prepare parameters for triangular solver
        //
        maxgrowth = 1/rcondthreshold<Precision>();
        su = 0;
        sl = 1;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=i-1; j++)
            {
                sl = amp::maximum<Precision>(sl, amp::abscomplex<Precision>(lua(i,j)));
            }
            for(j=i; j<=n-1; j++)
            {
                su = amp::maximum<Precision>(su, amp::abscomplex<Precision>(lua(i,j)));
            }
        }
        if( su==0 )
        {
            su = 1;
        }
        su = 1/su;
        sl = 1/sl;
        
        //
        // Estimate the norm of SU*SL*A.
        //
        if( !isanormprovided )
        {
            anorm = 0;
            if( onenorm )
            {
                kase1 = 1;
            }
            else
            {
                kase1 = 2;
            }
            kase = 0;
            do
            {
                cmatrixestimatenorm<Precision>(n, cwork4, ex, anorm, kase, isave, rsave);
                if( kase!=0 )
                {
                    if( kase==kase1 )
                    {
                        
                        //
                        // Multiply by U
                        //
                        for(i=1; i<=n; i++)
                        {
                            i1_ = (i)-(i-1);
                            v = 0.0;
                            for(i_=i-1; i_<=n-1;i_++)
                            {
                                v += lua(i-1,i_)*ex(i_+i1_);
                            }
                            ex(i) = v;
                        }
                        
                        //
                        // Multiply by L
                        //
                        for(i=n; i>=1; i--)
                        {
                            v = 0;
                            if( i>1 )
                            {
                                i1_ = (1)-(0);
                                v = 0.0;
                                for(i_=0; i_<=i-2;i_++)
                                {
                                    v += lua(i-1,i_)*ex(i_+i1_);
                                }
                            }
                            ex(i) = v+ex(i);
                        }
                    }
                    else
                    {
                        
                        //
                        // Multiply by L'
                        //
                        for(i=1; i<=n; i++)
                        {
                            cwork2(i) = 0;
                        }
                        for(i=1; i<=n; i++)
                        {
                            v = ex(i);
                            if( i>1 )
                            {
                                i1_ = (0) - (1);
                                for(i_=1; i_<=i-1;i_++)
                                {
                                    cwork2(i_) = cwork2(i_) + v*amp::conj(lua(i-1,i_+i1_));
                                }
                            }
                            cwork2(i) = cwork2(i)+v;
                        }
                        
                        //
                        // Multiply by U'
                        //
                        for(i=1; i<=n; i++)
                        {
                            ex(i) = 0;
                        }
                        for(i=1; i<=n; i++)
                        {
                            v = cwork2(i);
                            i1_ = (i-1) - (i);
                            for(i_=i; i_<=n;i_++)
                            {
                                ex(i_) = ex(i_) + v*amp::conj(lua(i-1,i_+i1_));
                            }
                        }
                    }
                }
            }
            while( kase!=0 );
        }
        
        //
        // Scale according to SU/SL
        //
        anorm = anorm*su*sl;
        
        //
        // Quick return if possible
        //
        if( anorm==0 )
        {
            return;
        }
        
        //
        // Estimate the norm of inv(A).
        //
        ainvnm = 0;
        if( onenorm )
        {
            kase1 = 1;
        }
        else
        {
            kase1 = 2;
        }
        kase = 0;
        while( true )
        {
            cmatrixestimatenorm<Precision>(n, cwork4, ex, ainvnm, kase, isave, rsave);
            if( kase==0 )
            {
                break;
            }
            
            //
            // From 1-based to 0-based
            //
            for(i=0; i<=n-1; i++)
            {
                ex(i) = ex(i+1);
            }
            
            //
            // multiply by inv(A) or inv(A')
            //
            if( kase==kase1 )
            {
                
                //
                // Multiply by inv(L).
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(lua, sl, n, ex, false, 0, true, maxgrowth) )
                {
                    rc = 0;
                    return;
                }
                
                //
                // Multiply by inv(U).
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(lua, su, n, ex, true, 0, false, maxgrowth) )
                {
                    rc = 0;
                    return;
                }
            }
            else
            {
                
                //
                // Multiply by inv(U').
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(lua, su, n, ex, true, 2, false, maxgrowth) )
                {
                    rc = 0;
                    return;
                }
                
                //
                // Multiply by inv(L').
                //
                if( !safesolve::cmatrixscaledtrsafesolve<Precision>(lua, sl, n, ex, false, 2, true, maxgrowth) )
                {
                    rc = 0;
                    return;
                }
            }
            
            //
            // from 0-based to 1-based
            //
            for(i=n-1; i>=0; i--)
            {
                ex(i+1) = ex(i);
            }
        }
        
        //
        // Compute the estimate of the reciprocal condition number.
        //
        if( ainvnm!=0 )
        {
            rc = 1/ainvnm;
            rc = rc/anorm;
            if( rc<rcondthreshold<Precision>() )
            {
                rc = 0;
            }
        }
    }


    /*************************************************************************
    Internal subroutine for matrix norm estimation

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixestimatenorm(int n,
        ap::template_1d_array< amp::ampf<Precision> >& v,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< int >& isgn,
        amp::ampf<Precision>& est,
        int& kase)
    {
        int itmax;
        int i;
        amp::ampf<Precision> t;
        bool flg;
        int positer;
        int posj;
        int posjlast;
        int posjump;
        int posaltsgn;
        int posestold;
        int postemp;


        itmax = 5;
        posaltsgn = n+1;
        posestold = n+2;
        postemp = n+3;
        positer = n+1;
        posj = n+2;
        posjlast = n+3;
        posjump = n+4;
        if( kase==0 )
        {
            v.setlength(n+4);
            x.setlength(n+1);
            isgn.setlength(n+5);
            t = amp::ampf<Precision>(1)/amp::ampf<Precision>(n);
            for(i=1; i<=n; i++)
            {
                x(i) = t;
            }
            kase = 1;
            isgn(posjump) = 1;
            return;
        }
        
        //
        //     ................ ENTRY   (JUMP = 1)
        //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( isgn(posjump)==1 )
        {
            if( n==1 )
            {
                v(1) = x(1);
                est = amp::abs<Precision>(v(1));
                kase = 0;
                return;
            }
            est = 0;
            for(i=1; i<=n; i++)
            {
                est = est+amp::abs<Precision>(x(i));
            }
            for(i=1; i<=n; i++)
            {
                if( x(i)>=0 )
                {
                    x(i) = 1;
                }
                else
                {
                    x(i) = -1;
                }
                isgn(i) = amp::sign<Precision>(x(i));
            }
            kase = 2;
            isgn(posjump) = 2;
            return;
        }
        
        //
        //     ................ ENTRY   (JUMP = 2)
        //     FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
        //
        if( isgn(posjump)==2 )
        {
            isgn(posj) = 1;
            for(i=2; i<=n; i++)
            {
                if( amp::abs<Precision>(x(i))>amp::abs<Precision>(x(isgn(posj))) )
                {
                    isgn(posj) = i;
                }
            }
            isgn(positer) = 2;
            
            //
            // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
            //
            for(i=1; i<=n; i++)
            {
                x(i) = 0;
            }
            x(isgn(posj)) = 1;
            kase = 1;
            isgn(posjump) = 3;
            return;
        }
        
        //
        //     ................ ENTRY   (JUMP = 3)
        //     X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( isgn(posjump)==3 )
        {
            amp::vmove(v.getvector(1, n), x.getvector(1, n));
            v(posestold) = est;
            est = 0;
            for(i=1; i<=n; i++)
            {
                est = est+amp::abs<Precision>(v(i));
            }
            flg = false;
            for(i=1; i<=n; i++)
            {
                if( x(i)>=0 && isgn(i)<0 || x(i)<0 && isgn(i)>=0 )
                {
                    flg = true;
                }
            }
            
            //
            // REPEATED SIGN VECTOR DETECTED, HENCE ALGORITHM HAS CONVERGED.
            // OR MAY BE CYCLING.
            //
            if( !flg || est<=v(posestold) )
            {
                v(posaltsgn) = 1;
                for(i=1; i<=n; i++)
                {
                    x(i) = v(posaltsgn)*(1+(amp::ampf<Precision>(i-1))/(amp::ampf<Precision>(n-1)));
                    v(posaltsgn) = -v(posaltsgn);
                }
                kase = 1;
                isgn(posjump) = 5;
                return;
            }
            for(i=1; i<=n; i++)
            {
                if( x(i)>=0 )
                {
                    x(i) = 1;
                    isgn(i) = 1;
                }
                else
                {
                    x(i) = -1;
                    isgn(i) = -1;
                }
            }
            kase = 2;
            isgn(posjump) = 4;
            return;
        }
        
        //
        //     ................ ENTRY   (JUMP = 4)
        //     X HAS BEEN OVERWRITTEN BY TRANDPOSE(A)*X.
        //
        if( isgn(posjump)==4 )
        {
            isgn(posjlast) = isgn(posj);
            isgn(posj) = 1;
            for(i=2; i<=n; i++)
            {
                if( amp::abs<Precision>(x(i))>amp::abs<Precision>(x(isgn(posj))) )
                {
                    isgn(posj) = i;
                }
            }
            if( x(isgn(posjlast))!=amp::abs<Precision>(x(isgn(posj))) && isgn(positer)<itmax )
            {
                isgn(positer) = isgn(positer)+1;
                for(i=1; i<=n; i++)
                {
                    x(i) = 0;
                }
                x(isgn(posj)) = 1;
                kase = 1;
                isgn(posjump) = 3;
                return;
            }
            
            //
            // ITERATION COMPLETE.  FINAL STAGE.
            //
            v(posaltsgn) = 1;
            for(i=1; i<=n; i++)
            {
                x(i) = v(posaltsgn)*(1+(amp::ampf<Precision>(i-1))/(amp::ampf<Precision>(n-1)));
                v(posaltsgn) = -v(posaltsgn);
            }
            kase = 1;
            isgn(posjump) = 5;
            return;
        }
        
        //
        //     ................ ENTRY   (JUMP = 5)
        //     X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( isgn(posjump)==5 )
        {
            v(postemp) = 0;
            for(i=1; i<=n; i++)
            {
                v(postemp) = v(postemp)+amp::abs<Precision>(x(i));
            }
            v(postemp) = 2*v(postemp)/(3*n);
            if( v(postemp)>est )
            {
                amp::vmove(v.getvector(1, n), x.getvector(1, n));
                est = v(postemp);
            }
            kase = 0;
            return;
        }
    }


    template<unsigned int Precision>
    void cmatrixestimatenorm(const int& n,
        ap::template_1d_array< amp::campf<Precision> >& v,
        ap::template_1d_array< amp::campf<Precision> >& x,
        amp::ampf<Precision>& est,
        int& kase,
        ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave)
    {
        int itmax;
        int i;
        int iter;
        int j;
        int jlast;
        int jump;
        amp::ampf<Precision> absxi;
        amp::ampf<Precision> altsgn;
        amp::ampf<Precision> estold;
        amp::ampf<Precision> safmin;
        amp::ampf<Precision> temp;
        int i_;


        
        //
        //Executable Statements ..
        //
        itmax = 5;
        safmin = amp::ampf<Precision>::getAlgoPascalMinNumber();
        if( kase==0 )
        {
            v.setlength(n+1);
            x.setlength(n+1);
            isave.setlength(5);
            rsave.setlength(4);
            for(i=1; i<=n; i++)
            {
                x(i) = amp::ampf<Precision>(1)/amp::ampf<Precision>(n);
            }
            kase = 1;
            jump = 1;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        internalcomplexrcondloadall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
        
        //
        // ENTRY   (JUMP = 1)
        // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( jump==1 )
        {
            if( n==1 )
            {
                v(1) = x(1);
                est = amp::abscomplex<Precision>(v(1));
                kase = 0;
                internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
                return;
            }
            est = internalcomplexrcondscsum1<Precision>(x, n);
            for(i=1; i<=n; i++)
            {
                absxi = amp::abscomplex<Precision>(x(i));
                if( absxi>safmin )
                {
                    x(i) = x(i)/absxi;
                }
                else
                {
                    x(i) = 1;
                }
            }
            kase = 2;
            jump = 2;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        
        //
        // ENTRY   (JUMP = 2)
        // FIRST ITERATION.  X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
        //
        if( jump==2 )
        {
            j = internalcomplexrcondicmax1<Precision>(x, n);
            iter = 2;
            
            //
            // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
            //
            for(i=1; i<=n; i++)
            {
                x(i) = 0;
            }
            x(j) = 1;
            kase = 1;
            jump = 3;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        
        //
        // ENTRY   (JUMP = 3)
        // X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( jump==3 )
        {
            for(i_=1; i_<=n;i_++)
            {
                v(i_) = x(i_);
            }
            estold = est;
            est = internalcomplexrcondscsum1<Precision>(v, n);
            
            //
            // TEST FOR CYCLING.
            //
            if( est<=estold )
            {
                
                //
                // ITERATION COMPLETE.  FINAL STAGE.
                //
                altsgn = 1;
                for(i=1; i<=n; i++)
                {
                    x(i) = altsgn*(1+(amp::ampf<Precision>(i-1))/(amp::ampf<Precision>(n-1)));
                    altsgn = -altsgn;
                }
                kase = 1;
                jump = 5;
                internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
                return;
            }
            for(i=1; i<=n; i++)
            {
                absxi = amp::abscomplex<Precision>(x(i));
                if( absxi>safmin )
                {
                    x(i) = x(i)/absxi;
                }
                else
                {
                    x(i) = 1;
                }
            }
            kase = 2;
            jump = 4;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        
        //
        // ENTRY   (JUMP = 4)
        // X HAS BEEN OVERWRITTEN BY CTRANS(A)*X.
        //
        if( jump==4 )
        {
            jlast = j;
            j = internalcomplexrcondicmax1<Precision>(x, n);
            if( amp::abscomplex<Precision>(x(jlast))!=amp::abscomplex<Precision>(x(j)) && iter<itmax )
            {
                iter = iter+1;
                
                //
                // MAIN LOOP - ITERATIONS 2,3,...,ITMAX.
                //
                for(i=1; i<=n; i++)
                {
                    x(i) = 0;
                }
                x(j) = 1;
                kase = 1;
                jump = 3;
                internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
                return;
            }
            
            //
            // ITERATION COMPLETE.  FINAL STAGE.
            //
            altsgn = 1;
            for(i=1; i<=n; i++)
            {
                x(i) = altsgn*(1+(amp::ampf<Precision>(i-1))/(amp::ampf<Precision>(n-1)));
                altsgn = -altsgn;
            }
            kase = 1;
            jump = 5;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
        
        //
        // ENTRY   (JUMP = 5)
        // X HAS BEEN OVERWRITTEN BY A*X.
        //
        if( jump==5 )
        {
            temp = 2*(internalcomplexrcondscsum1<Precision>(x, n)/(3*n));
            if( temp>est )
            {
                for(i_=1; i_<=n;i_++)
                {
                    v(i_) = x(i_);
                }
                est = temp;
            }
            kase = 0;
            internalcomplexrcondsaveall<Precision>(isave, rsave, i, iter, j, jlast, jump, absxi, altsgn, estold, temp);
            return;
        }
    }


    template<unsigned int Precision>
    amp::ampf<Precision> internalcomplexrcondscsum1(const ap::template_1d_array< amp::campf<Precision> >& x,
        int n)
    {
        amp::ampf<Precision> result;
        int i;


        result = 0;
        for(i=1; i<=n; i++)
        {
            result = result+amp::abscomplex<Precision>(x(i));
        }
        return result;
    }


    template<unsigned int Precision>
    int internalcomplexrcondicmax1(const ap::template_1d_array< amp::campf<Precision> >& x,
        int n)
    {
        int result;
        int i;
        amp::ampf<Precision> m;


        result = 1;
        m = amp::abscomplex<Precision>(x(1));
        for(i=2; i<=n; i++)
        {
            if( amp::abscomplex<Precision>(x(i))>m )
            {
                result = i;
                m = amp::abscomplex<Precision>(x(i));
            }
        }
        return result;
    }


    template<unsigned int Precision>
    void internalcomplexrcondsaveall(ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave,
        int& i,
        int& iter,
        int& j,
        int& jlast,
        int& jump,
        amp::ampf<Precision>& absxi,
        amp::ampf<Precision>& altsgn,
        amp::ampf<Precision>& estold,
        amp::ampf<Precision>& temp)
    {
        isave(0) = i;
        isave(1) = iter;
        isave(2) = j;
        isave(3) = jlast;
        isave(4) = jump;
        rsave(0) = absxi;
        rsave(1) = altsgn;
        rsave(2) = estold;
        rsave(3) = temp;
    }


    template<unsigned int Precision>
    void internalcomplexrcondloadall(ap::template_1d_array< int >& isave,
        ap::template_1d_array< amp::ampf<Precision> >& rsave,
        int& i,
        int& iter,
        int& j,
        int& jlast,
        int& jump,
        amp::ampf<Precision>& absxi,
        amp::ampf<Precision>& altsgn,
        amp::ampf<Precision>& estold,
        amp::ampf<Precision>& temp)
    {
        i = isave(0);
        iter = isave(1);
        j = isave(2);
        jlast = isave(3);
        jump = isave(4);
        absxi = rsave(0);
        altsgn = rsave(1);
        estold = rsave(2);
        temp = rsave(3);
    }
} // namespace

#endif
