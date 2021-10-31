/*************************************************************************
Copyright (c) 1992-2007 The University of Tennessee. All rights reserved.

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

#ifndef _trfac_h
#define _trfac_h

#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
namespace trfac
{
    template<unsigned int Precision>
    void rmatrixlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void cmatrixlu(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    bool hpdmatrixcholesky(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    bool spdmatrixcholesky(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    void rmatrixlup(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void cmatrixlup(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void rmatrixplu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void cmatrixplu(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots);
    template<unsigned int Precision>
    void cmatrixluprec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    void rmatrixluprec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void cmatrixplurec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    void rmatrixplurec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void cmatrixlup2(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    void rmatrixlup2(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void cmatrixplu2(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    void rmatrixplu2(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    bool hpdmatrixcholeskyrec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    bool spdmatrixcholeskyrec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    bool hpdmatrixcholesky2(ap::template_2d_array< amp::campf<Precision> >& aaa,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    bool spdmatrixcholesky2(ap::template_2d_array< amp::ampf<Precision> >& aaa,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);


    /*************************************************************************
    LU decomposition of a general real matrix with row pivoting

    A is represented as A = P*L*U, where:
    * L is lower unitriangular matrix
    * U is upper triangular matrix
    * P = P0*P1*...*PK, K=min(M,N)-1,
      Pi - permutation matrix for I and Pivots[I]

    This is cache-oblivous implementation of LU decomposition.
    It is optimized for square matrices. As for rectangular matrices:
    * best case - M>>N
    * worst case - N>>M, small M, large N, matrix does not fit in CPU cache

    INPUT PARAMETERS:
        A       -   array[0..M-1, 0..N-1].
        M       -   number of rows in matrix A.
        N       -   number of columns in matrix A.


    OUTPUT PARAMETERS:
        A       -   matrices L and U in compact form:
                    * L is stored under main diagonal
                    * U is stored on and above main diagonal
        Pivots  -   permutation matrix in compact form.
                    array[0..Min(M-1,N-1)].

      -- ALGLIB routine --
         10.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        rmatrixplu<Precision>(a, m, n, pivots);
    }


    /*************************************************************************
    LU decomposition of a general complex matrix with row pivoting

    A is represented as A = P*L*U, where:
    * L is lower unitriangular matrix
    * U is upper triangular matrix
    * P = P0*P1*...*PK, K=min(M,N)-1,
      Pi - permutation matrix for I and Pivots[I]

    This is cache-oblivous implementation of LU decomposition. It is optimized
    for square matrices. As for rectangular matrices:
    * best case - M>>N
    * worst case - N>>M, small M, large N, matrix does not fit in CPU cache

    INPUT PARAMETERS:
        A       -   array[0..M-1, 0..N-1].
        M       -   number of rows in matrix A.
        N       -   number of columns in matrix A.


    OUTPUT PARAMETERS:
        A       -   matrices L and U in compact form:
                    * L is stored under main diagonal
                    * U is stored on and above main diagonal
        Pivots  -   permutation matrix in compact form.
                    array[0..Min(M-1,N-1)].

      -- ALGLIB routine --
         10.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlu(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        cmatrixplu<Precision>(a, m, n, pivots);
    }


    /*************************************************************************
    Cache-oblivious Cholesky decomposition

    The algorithm computes Cholesky decomposition  of  a  Hermitian  positive-
    definite matrix. The result of an algorithm is a representation  of  A  as
    A=U'*U  or A=L*L' (here X' detones conj(X^T)).

    INPUT PARAMETERS:
        A       -   upper or lower triangle of a factorized matrix.
                    array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   if IsUpper=True, then A contains an upper triangle of
                    a symmetric matrix, otherwise A contains a lower one.

    OUTPUT PARAMETERS:
        A       -   the result of factorization. If IsUpper=True, then
                    the upper triangle contains matrix U, so that A = U'*U,
                    and the elements below the main diagonal are not modified.
                    Similarly, if IsUpper = False.

    RESULT:
        If  the  matrix  is  positive-definite,  the  function  returns  True.
        Otherwise, the function returns False. Contents of A is not determined
        in such case.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool hpdmatrixcholesky(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< amp::campf<Precision> > tmp;


        if( n<1 )
        {
            result = false;
            return result;
        }
        tmp.setlength(2*n);
        result = hpdmatrixcholeskyrec<Precision>(a, 0, n, isupper, tmp);
        return result;
    }


    /*************************************************************************
    Cache-oblivious Cholesky decomposition

    The algorithm computes Cholesky decomposition  of  a  symmetric  positive-
    definite matrix. The result of an algorithm is a representation  of  A  as
    A=U^T*U  or A=L*L^T

    INPUT PARAMETERS:
        A       -   upper or lower triangle of a factorized matrix.
                    array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   if IsUpper=True, then A contains an upper triangle of
                    a symmetric matrix, otherwise A contains a lower one.

    OUTPUT PARAMETERS:
        A       -   the result of factorization. If IsUpper=True, then
                    the upper triangle contains matrix U, so that A = U^T*U,
                    and the elements below the main diagonal are not modified.
                    Similarly, if IsUpper = False.

    RESULT:
        If  the  matrix  is  positive-definite,  the  function  returns  True.
        Otherwise, the function returns False. Contents of A is not determined
        in such case.

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool spdmatrixcholesky(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > tmp;


        if( n<1 )
        {
            result = false;
            return result;
        }
        tmp.setlength(2*n);
        result = spdmatrixcholeskyrec<Precision>(a, 0, n, isupper, tmp);
        return result;
    }


    template<unsigned int Precision>
    void rmatrixlup(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        int i;
        int j;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;


        
        //
        // Internal LU decomposition subroutine.
        // Never call it directly.
        //
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        
        //
        // Scale matrix to avoid overflows,
        // decompose it, then scale back.
        //
        mx = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abs<Precision>(a(i,j)));
            }
        }
        if( mx!=0 )
        {
            v = 1/mx;
            for(i=0; i<=m-1; i++)
            {
                amp::vmul(a.getrow(i, 0, n-1), v);
            }
        }
        pivots.setlength(ap::minint(m, n));
        tmp.setlength(2*ap::maxint(m, n));
        rmatrixluprec<Precision>(a, 0, m, n, pivots, tmp);
        if( mx!=0 )
        {
            v = mx;
            for(i=0; i<=m-1; i++)
            {
                amp::vmul(a.getrow(i, 0, ap::minint(i, n-1)), v);
            }
        }
    }


    template<unsigned int Precision>
    void cmatrixlup(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::template_1d_array< amp::campf<Precision> > tmp;
        int i;
        int j;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;
        int i_;


        
        //
        // Internal LU decomposition subroutine.
        // Never call it directly.
        //
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        
        //
        // Scale matrix to avoid overflows,
        // decompose it, then scale back.
        //
        mx = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( mx!=0 )
        {
            v = 1/mx;
            for(i=0; i<=m-1; i++)
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    a(i,i_) = v*a(i,i_);
                }
            }
        }
        pivots.setlength(ap::minint(m, n));
        tmp.setlength(2*ap::maxint(m, n));
        cmatrixluprec<Precision>(a, 0, m, n, pivots, tmp);
        if( mx!=0 )
        {
            v = mx;
            for(i=0; i<=m-1; i++)
            {
                for(i_=0; i_<=ap::minint(i, n-1);i_++)
                {
                    a(i,i_) = v*a(i,i_);
                }
            }
        }
    }


    template<unsigned int Precision>
    void rmatrixplu(ap::template_2d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        int i;
        int j;
        amp::ampf<Precision> mx;
        amp::ampf<Precision> v;


        
        //
        // Internal LU decomposition subroutine.
        // Never call it directly.
        //
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        tmp.setlength(2*ap::maxint(m, n));
        pivots.setlength(ap::minint(m, n));
        
        //
        // Scale matrix to avoid overflows,
        // decompose it, then scale back.
        //
        mx = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abs<Precision>(a(i,j)));
            }
        }
        if( mx!=0 )
        {
            v = 1/mx;
            for(i=0; i<=m-1; i++)
            {
                amp::vmul(a.getrow(i, 0, n-1), v);
            }
        }
        rmatrixplurec<Precision>(a, 0, m, n, pivots, tmp);
        if( mx!=0 )
        {
            v = mx;
            for(i=0; i<=ap::minint(m, n)-1; i++)
            {
                amp::vmul(a.getrow(i, i, n-1), v);
            }
        }
    }


    template<unsigned int Precision>
    void cmatrixplu(ap::template_2d_array< amp::campf<Precision> >& a,
        int m,
        int n,
        ap::template_1d_array< int >& pivots)
    {
        ap::template_1d_array< amp::campf<Precision> > tmp;
        int i;
        int j;
        amp::ampf<Precision> mx;
        amp::campf<Precision> v;
        int i_;


        
        //
        // Internal LU decomposition subroutine.
        // Never call it directly.
        //
        ap::ap_error::make_assertion(m>0);
        ap::ap_error::make_assertion(n>0);
        tmp.setlength(2*ap::maxint(m, n));
        pivots.setlength(ap::minint(m, n));
        
        //
        // Scale matrix to avoid overflows,
        // decompose it, then scale back.
        //
        mx = 0;
        for(i=0; i<=m-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                mx = amp::maximum<Precision>(mx, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( mx!=0 )
        {
            v = 1/mx;
            for(i=0; i<=m-1; i++)
            {
                for(i_=0; i_<=n-1;i_++)
                {
                    a(i,i_) = v*a(i,i_);
                }
            }
        }
        cmatrixplurec<Precision>(a, 0, m, n, pivots, tmp);
        if( mx!=0 )
        {
            v = mx;
            for(i=0; i<=ap::minint(m, n)-1; i++)
            {
                for(i_=i; i_<=n-1;i_++)
                {
                    a(i,i_) = v*a(i,i_);
                }
            }
        }
    }


    /*************************************************************************
    Recurrent complex LU subroutine.
    Never call it directly.

      -- ALGLIB routine --
         04.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixluprec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        int m1;
        int m2;
        int i_;
        int i1_;


        
        //
        // Kernel case
        //
        if( ap::minint(m, n)<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            cmatrixlup2<Precision>(a, offs, m, n, pivots, tmp);
            return;
        }
        
        //
        // Preliminary step, make N>=M
        //
        //     ( A1 )
        // A = (    ), where A1 is square
        //     ( A2 )
        //
        // Factorize A1, update A2
        //
        if( m>n )
        {
            cmatrixluprec<Precision>(a, offs, n, n, pivots, tmp);
            for(i=0; i<=n-1; i++)
            {
                i1_ = (offs+n) - (0);
                for(i_=0; i_<=m-n-1;i_++)
                {
                    tmp(i_) = a(i_+i1_,offs+i);
                }
                for(i_=offs+n; i_<=offs+m-1;i_++)
                {
                    a(i_,offs+i) = a(i_,pivots(offs+i));
                }
                i1_ = (0) - (offs+n);
                for(i_=offs+n; i_<=offs+m-1;i_++)
                {
                    a(i_,pivots(offs+i)) = tmp(i_+i1_);
                }
            }
            ablas::cmatrixrighttrsm<Precision>(m-n, n, a, offs, offs, true, true, 0, a, offs+n, offs);
            return;
        }
        
        //
        // Non-kernel case
        //
        ablas::ablascomplexsplitlength<Precision>(a, m, m1, m2);
        cmatrixluprec<Precision>(a, offs, m1, n, pivots, tmp);
        if( m2>0 )
        {
            for(i=0; i<=m1-1; i++)
            {
                if( offs+i!=pivots(offs+i) )
                {
                    i1_ = (offs+m1) - (0);
                    for(i_=0; i_<=m2-1;i_++)
                    {
                        tmp(i_) = a(i_+i1_,offs+i);
                    }
                    for(i_=offs+m1; i_<=offs+m-1;i_++)
                    {
                        a(i_,offs+i) = a(i_,pivots(offs+i));
                    }
                    i1_ = (0) - (offs+m1);
                    for(i_=offs+m1; i_<=offs+m-1;i_++)
                    {
                        a(i_,pivots(offs+i)) = tmp(i_+i1_);
                    }
                }
            }
            ablas::cmatrixrighttrsm<Precision>(m2, m1, a, offs, offs, true, true, 0, a, offs+m1, offs);
            ablas::cmatrixgemm<Precision>(m-m1, n-m1, m1, -amp::ampf<Precision>("1.0"), a, offs+m1, offs, 0, a, offs, offs+m1, 0, +amp::ampf<Precision>("1.0"), a, offs+m1, offs+m1);
            cmatrixluprec<Precision>(a, offs+m1, m-m1, n-m1, pivots, tmp);
            for(i=0; i<=m2-1; i++)
            {
                if( offs+m1+i!=pivots(offs+m1+i) )
                {
                    i1_ = (offs) - (0);
                    for(i_=0; i_<=m1-1;i_++)
                    {
                        tmp(i_) = a(i_+i1_,offs+m1+i);
                    }
                    for(i_=offs; i_<=offs+m1-1;i_++)
                    {
                        a(i_,offs+m1+i) = a(i_,pivots(offs+m1+i));
                    }
                    i1_ = (0) - (offs);
                    for(i_=offs; i_<=offs+m1-1;i_++)
                    {
                        a(i_,pivots(offs+m1+i)) = tmp(i_+i1_);
                    }
                }
            }
        }
    }


    /*************************************************************************
    Recurrent real LU subroutine.
    Never call it directly.

      -- ALGLIB routine --
         04.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixluprec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        int m1;
        int m2;


        
        //
        // Kernel case
        //
        if( ap::minint(m, n)<=ablas::ablasblocksize<Precision>(a) )
        {
            rmatrixlup2<Precision>(a, offs, m, n, pivots, tmp);
            return;
        }
        
        //
        // Preliminary step, make N>=M
        //
        //     ( A1 )
        // A = (    ), where A1 is square
        //     ( A2 )
        //
        // Factorize A1, update A2
        //
        if( m>n )
        {
            rmatrixluprec<Precision>(a, offs, n, n, pivots, tmp);
            for(i=0; i<=n-1; i++)
            {
                if( offs+i!=pivots(offs+i) )
                {
                    amp::vmove(tmp.getvector(0, m-n-1), a.getcolumn(offs+i, offs+n, offs+m-1));
                    amp::vmove(a.getcolumn(offs+i, offs+n, offs+m-1), a.getcolumn(pivots(offs+i), offs+n, offs+m-1));
                    amp::vmove(a.getcolumn(pivots(offs+i), offs+n, offs+m-1), tmp.getvector(0, m-n-1));
                }
            }
            ablas::rmatrixrighttrsm<Precision>(m-n, n, a, offs, offs, true, true, 0, a, offs+n, offs);
            return;
        }
        
        //
        // Non-kernel case
        //
        ablas::ablassplitlength<Precision>(a, m, m1, m2);
        rmatrixluprec<Precision>(a, offs, m1, n, pivots, tmp);
        if( m2>0 )
        {
            for(i=0; i<=m1-1; i++)
            {
                if( offs+i!=pivots(offs+i) )
                {
                    amp::vmove(tmp.getvector(0, m2-1), a.getcolumn(offs+i, offs+m1, offs+m-1));
                    amp::vmove(a.getcolumn(offs+i, offs+m1, offs+m-1), a.getcolumn(pivots(offs+i), offs+m1, offs+m-1));
                    amp::vmove(a.getcolumn(pivots(offs+i), offs+m1, offs+m-1), tmp.getvector(0, m2-1));
                }
            }
            ablas::rmatrixrighttrsm<Precision>(m2, m1, a, offs, offs, true, true, 0, a, offs+m1, offs);
            ablas::rmatrixgemm<Precision>(m-m1, n-m1, m1, -amp::ampf<Precision>("1.0"), a, offs+m1, offs, 0, a, offs, offs+m1, 0, +amp::ampf<Precision>("1.0"), a, offs+m1, offs+m1);
            rmatrixluprec<Precision>(a, offs+m1, m-m1, n-m1, pivots, tmp);
            for(i=0; i<=m2-1; i++)
            {
                if( offs+m1+i!=pivots(offs+m1+i) )
                {
                    amp::vmove(tmp.getvector(0, m1-1), a.getcolumn(offs+m1+i, offs, offs+m1-1));
                    amp::vmove(a.getcolumn(offs+m1+i, offs, offs+m1-1), a.getcolumn(pivots(offs+m1+i), offs, offs+m1-1));
                    amp::vmove(a.getcolumn(pivots(offs+m1+i), offs, offs+m1-1), tmp.getvector(0, m1-1));
                }
            }
        }
    }


    /*************************************************************************
    Recurrent complex LU subroutine.
    Never call it directly.

      -- ALGLIB routine --
         04.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixplurec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        int n1;
        int n2;
        int i_;
        int i1_;


        
        //
        // Kernel case
        //
        if( ap::minint(m, n)<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            cmatrixplu2<Precision>(a, offs, m, n, pivots, tmp);
            return;
        }
        
        //
        // Preliminary step, make M>=N.
        //
        // A = (A1 A2), where A1 is square
        // Factorize A1, update A2
        //
        if( n>m )
        {
            cmatrixplurec<Precision>(a, offs, m, m, pivots, tmp);
            for(i=0; i<=m-1; i++)
            {
                i1_ = (offs+m) - (0);
                for(i_=0; i_<=n-m-1;i_++)
                {
                    tmp(i_) = a(offs+i,i_+i1_);
                }
                for(i_=offs+m; i_<=offs+n-1;i_++)
                {
                    a(offs+i,i_) = a(pivots(offs+i),i_);
                }
                i1_ = (0) - (offs+m);
                for(i_=offs+m; i_<=offs+n-1;i_++)
                {
                    a(pivots(offs+i),i_) = tmp(i_+i1_);
                }
            }
            ablas::cmatrixlefttrsm<Precision>(m, n-m, a, offs, offs, false, true, 0, a, offs, offs+m);
            return;
        }
        
        //
        // Non-kernel case
        //
        ablas::ablascomplexsplitlength<Precision>(a, n, n1, n2);
        cmatrixplurec<Precision>(a, offs, m, n1, pivots, tmp);
        if( n2>0 )
        {
            for(i=0; i<=n1-1; i++)
            {
                if( offs+i!=pivots(offs+i) )
                {
                    i1_ = (offs+n1) - (0);
                    for(i_=0; i_<=n2-1;i_++)
                    {
                        tmp(i_) = a(offs+i,i_+i1_);
                    }
                    for(i_=offs+n1; i_<=offs+n-1;i_++)
                    {
                        a(offs+i,i_) = a(pivots(offs+i),i_);
                    }
                    i1_ = (0) - (offs+n1);
                    for(i_=offs+n1; i_<=offs+n-1;i_++)
                    {
                        a(pivots(offs+i),i_) = tmp(i_+i1_);
                    }
                }
            }
            ablas::cmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, false, true, 0, a, offs, offs+n1);
            ablas::cmatrixgemm<Precision>(m-n1, n-n1, n1, -amp::ampf<Precision>("1.0"), a, offs+n1, offs, 0, a, offs, offs+n1, 0, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1);
            cmatrixplurec<Precision>(a, offs+n1, m-n1, n-n1, pivots, tmp);
            for(i=0; i<=n2-1; i++)
            {
                if( offs+n1+i!=pivots(offs+n1+i) )
                {
                    i1_ = (offs) - (0);
                    for(i_=0; i_<=n1-1;i_++)
                    {
                        tmp(i_) = a(offs+n1+i,i_+i1_);
                    }
                    for(i_=offs; i_<=offs+n1-1;i_++)
                    {
                        a(offs+n1+i,i_) = a(pivots(offs+n1+i),i_);
                    }
                    i1_ = (0) - (offs);
                    for(i_=offs; i_<=offs+n1-1;i_++)
                    {
                        a(pivots(offs+n1+i),i_) = tmp(i_+i1_);
                    }
                }
            }
        }
    }


    /*************************************************************************
    Recurrent real LU subroutine.
    Never call it directly.

      -- ALGLIB routine --
         04.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixplurec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        int n1;
        int n2;


        
        //
        // Kernel case
        //
        if( ap::minint(m, n)<=ablas::ablasblocksize<Precision>(a) )
        {
            rmatrixplu2<Precision>(a, offs, m, n, pivots, tmp);
            return;
        }
        
        //
        // Preliminary step, make M>=N.
        //
        // A = (A1 A2), where A1 is square
        // Factorize A1, update A2
        //
        if( n>m )
        {
            rmatrixplurec<Precision>(a, offs, m, m, pivots, tmp);
            for(i=0; i<=m-1; i++)
            {
                amp::vmove(tmp.getvector(0, n-m-1), a.getrow(offs+i, offs+m, offs+n-1));
                amp::vmove(a.getrow(offs+i, offs+m, offs+n-1), a.getrow(pivots(offs+i), offs+m, offs+n-1));
                amp::vmove(a.getrow(pivots(offs+i), offs+m, offs+n-1), tmp.getvector(0, n-m-1));
            }
            ablas::rmatrixlefttrsm<Precision>(m, n-m, a, offs, offs, false, true, 0, a, offs, offs+m);
            return;
        }
        
        //
        // Non-kernel case
        //
        ablas::ablassplitlength<Precision>(a, n, n1, n2);
        rmatrixplurec<Precision>(a, offs, m, n1, pivots, tmp);
        if( n2>0 )
        {
            for(i=0; i<=n1-1; i++)
            {
                if( offs+i!=pivots(offs+i) )
                {
                    amp::vmove(tmp.getvector(0, n2-1), a.getrow(offs+i, offs+n1, offs+n-1));
                    amp::vmove(a.getrow(offs+i, offs+n1, offs+n-1), a.getrow(pivots(offs+i), offs+n1, offs+n-1));
                    amp::vmove(a.getrow(pivots(offs+i), offs+n1, offs+n-1), tmp.getvector(0, n2-1));
                }
            }
            ablas::rmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, false, true, 0, a, offs, offs+n1);
            ablas::rmatrixgemm<Precision>(m-n1, n-n1, n1, -amp::ampf<Precision>("1.0"), a, offs+n1, offs, 0, a, offs, offs+n1, 0, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1);
            rmatrixplurec<Precision>(a, offs+n1, m-n1, n-n1, pivots, tmp);
            for(i=0; i<=n2-1; i++)
            {
                if( offs+n1+i!=pivots(offs+n1+i) )
                {
                    amp::vmove(tmp.getvector(0, n1-1), a.getrow(offs+n1+i, offs, offs+n1-1));
                    amp::vmove(a.getrow(offs+n1+i, offs, offs+n1-1), a.getrow(pivots(offs+n1+i), offs, offs+n1-1));
                    amp::vmove(a.getrow(pivots(offs+n1+i), offs, offs+n1-1), tmp.getvector(0, n1-1));
                }
            }
        }
    }


    /*************************************************************************
    Complex LUP kernel

      -- ALGLIB routine --
         10.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlup2(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        int j;
        int jp;
        amp::campf<Precision> s;
        int i_;
        int i1_;


        
        //
        // Quick return if possible
        //
        if( m==0 || n==0 )
        {
            return;
        }
        
        //
        // main cycle
        //
        for(j=0; j<=ap::minint(m-1, n-1); j++)
        {
            
            //
            // Find pivot, swap columns
            //
            jp = j;
            for(i=j+1; i<=n-1; i++)
            {
                if( amp::abscomplex<Precision>(a(offs+j,offs+i))>amp::abscomplex<Precision>(a(offs+j,offs+jp)) )
                {
                    jp = i;
                }
            }
            pivots(offs+j) = offs+jp;
            if( jp!=j )
            {
                i1_ = (offs) - (0);
                for(i_=0; i_<=m-1;i_++)
                {
                    tmp(i_) = a(i_+i1_,offs+j);
                }
                for(i_=offs; i_<=offs+m-1;i_++)
                {
                    a(i_,offs+j) = a(i_,offs+jp);
                }
                i1_ = (0) - (offs);
                for(i_=offs; i_<=offs+m-1;i_++)
                {
                    a(i_,offs+jp) = tmp(i_+i1_);
                }
            }
            
            //
            // LU decomposition of 1x(N-J) matrix
            //
            if( a(offs+j,offs+j)!=0 && j+1<=n-1 )
            {
                s = 1/a(offs+j,offs+j);
                for(i_=offs+j+1; i_<=offs+n-1;i_++)
                {
                    a(offs+j,i_) = s*a(offs+j,i_);
                }
            }
            
            //
            // Update trailing (M-J-1)x(N-J-1) matrix
            //
            if( j<ap::minint(m-1, n-1) )
            {
                i1_ = (offs+j+1) - (0);
                for(i_=0; i_<=m-j-2;i_++)
                {
                    tmp(i_) = a(i_+i1_,offs+j);
                }
                i1_ = (offs+j+1) - (m);
                for(i_=m; i_<=m+n-j-2;i_++)
                {
                    tmp(i_) = -a(offs+j,i_+i1_);
                }
                ablas::cmatrixrank1<Precision>(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
            }
        }
    }


    /*************************************************************************
    Real LUP kernel

      -- ALGLIB routine --
         10.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlup2(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        int j;
        int jp;
        amp::ampf<Precision> s;


        
        //
        // Quick return if possible
        //
        if( m==0 || n==0 )
        {
            return;
        }
        
        //
        // main cycle
        //
        for(j=0; j<=ap::minint(m-1, n-1); j++)
        {
            
            //
            // Find pivot, swap columns
            //
            jp = j;
            for(i=j+1; i<=n-1; i++)
            {
                if( amp::abs<Precision>(a(offs+j,offs+i))>amp::abs<Precision>(a(offs+j,offs+jp)) )
                {
                    jp = i;
                }
            }
            pivots(offs+j) = offs+jp;
            if( jp!=j )
            {
                amp::vmove(tmp.getvector(0, m-1), a.getcolumn(offs+j, offs, offs+m-1));
                amp::vmove(a.getcolumn(offs+j, offs, offs+m-1), a.getcolumn(offs+jp, offs, offs+m-1));
                amp::vmove(a.getcolumn(offs+jp, offs, offs+m-1), tmp.getvector(0, m-1));
            }
            
            //
            // LU decomposition of 1x(N-J) matrix
            //
            if( a(offs+j,offs+j)!=0 && j+1<=n-1 )
            {
                s = 1/a(offs+j,offs+j);
                amp::vmul(a.getrow(offs+j, offs+j+1, offs+n-1), s);
            }
            
            //
            // Update trailing (M-J-1)x(N-J-1) matrix
            //
            if( j<ap::minint(m-1, n-1) )
            {
                amp::vmove(tmp.getvector(0, m-j-2), a.getcolumn(offs+j, offs+j+1, offs+m-1));
                amp::vmoveneg(tmp.getvector(m, m+n-j-2), a.getrow(offs+j, offs+j+1, offs+n-1));
                ablas::rmatrixrank1<Precision>(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
            }
        }
    }


    /*************************************************************************
    Complex PLU kernel

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         June 30, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixplu2(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        int j;
        int jp;
        amp::campf<Precision> s;
        int i_;
        int i1_;


        
        //
        // Quick return if possible
        //
        if( m==0 || n==0 )
        {
            return;
        }
        for(j=0; j<=ap::minint(m-1, n-1); j++)
        {
            
            //
            // Find pivot and test for singularity.
            //
            jp = j;
            for(i=j+1; i<=m-1; i++)
            {
                if( amp::abscomplex<Precision>(a(offs+i,offs+j))>amp::abscomplex<Precision>(a(offs+jp,offs+j)) )
                {
                    jp = i;
                }
            }
            pivots(offs+j) = offs+jp;
            if( a(offs+jp,offs+j)!=0 )
            {
                
                //
                //Apply the interchange to rows
                //
                if( jp!=j )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        s = a(offs+j,offs+i);
                        a(offs+j,offs+i) = a(offs+jp,offs+i);
                        a(offs+jp,offs+i) = s;
                    }
                }
                
                //
                //Compute elements J+1:M of J-th column.
                //
                if( j+1<=m-1 )
                {
                    s = 1/a(offs+j,offs+j);
                    for(i_=offs+j+1; i_<=offs+m-1;i_++)
                    {
                        a(i_,offs+j) = s*a(i_,offs+j);
                    }
                }
            }
            if( j<ap::minint(m, n)-1 )
            {
                
                //
                //Update trailing submatrix.
                //
                i1_ = (offs+j+1) - (0);
                for(i_=0; i_<=m-j-2;i_++)
                {
                    tmp(i_) = a(i_+i1_,offs+j);
                }
                i1_ = (offs+j+1) - (m);
                for(i_=m; i_<=m+n-j-2;i_++)
                {
                    tmp(i_) = -a(offs+j,i_+i1_);
                }
                ablas::cmatrixrank1<Precision>(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
            }
        }
    }


    /*************************************************************************
    Real PLU kernel

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         June 30, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixplu2(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int m,
        int n,
        ap::template_1d_array< int >& pivots,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        int j;
        int jp;
        amp::ampf<Precision> s;


        
        //
        // Quick return if possible
        //
        if( m==0 || n==0 )
        {
            return;
        }
        for(j=0; j<=ap::minint(m-1, n-1); j++)
        {
            
            //
            // Find pivot and test for singularity.
            //
            jp = j;
            for(i=j+1; i<=m-1; i++)
            {
                if( amp::abs<Precision>(a(offs+i,offs+j))>amp::abs<Precision>(a(offs+jp,offs+j)) )
                {
                    jp = i;
                }
            }
            pivots(offs+j) = offs+jp;
            if( a(offs+jp,offs+j)!=0 )
            {
                
                //
                //Apply the interchange to rows
                //
                if( jp!=j )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        s = a(offs+j,offs+i);
                        a(offs+j,offs+i) = a(offs+jp,offs+i);
                        a(offs+jp,offs+i) = s;
                    }
                }
                
                //
                //Compute elements J+1:M of J-th column.
                //
                if( j+1<=m-1 )
                {
                    s = 1/a(offs+j,offs+j);
                    amp::vmul(a.getcolumn(offs+j, offs+j+1, offs+m-1), s);
                }
            }
            if( j<ap::minint(m, n)-1 )
            {
                
                //
                //Update trailing submatrix.
                //
                amp::vmove(tmp.getvector(0, m-j-2), a.getcolumn(offs+j, offs+j+1, offs+m-1));
                amp::vmoveneg(tmp.getvector(m, m+n-j-2), a.getrow(offs+j, offs+j+1, offs+n-1));
                ablas::rmatrixrank1<Precision>(m-j-1, n-j-1, a, offs+j+1, offs+j+1, tmp, 0, tmp, m);
            }
        }
    }


    /*************************************************************************
    Recursive computational subroutine for HPDMatrixCholesky

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool hpdmatrixcholeskyrec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        bool result;
        int n1;
        int n2;


        
        //
        // check N
        //
        if( n<1 )
        {
            result = false;
            return result;
        }
        
        //
        // special cases
        //
        if( n==1 )
        {
            if( a(offs,offs).x>0 )
            {
                a(offs,offs) = amp::sqrt<Precision>(a(offs,offs).x);
                result = true;
            }
            else
            {
                result = false;
            }
            return result;
        }
        if( n<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            result = hpdmatrixcholesky2<Precision>(a, offs, n, isupper, tmp);
            return result;
        }
        
        //
        // general case: split task in cache-oblivious manner
        //
        result = true;
        ablas::ablascomplexsplitlength<Precision>(a, n, n1, n2);
        result = hpdmatrixcholeskyrec<Precision>(a, offs, n1, isupper, tmp);
        if( !result )
        {
            return result;
        }
        if( n2>0 )
        {
            if( isupper )
            {
                ablas::cmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, false, 2, a, offs, offs+n1);
                ablas::cmatrixsyrk<Precision>(n2, n1, -amp::ampf<Precision>("1.0"), a, offs, offs+n1, 2, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1, isupper);
            }
            else
            {
                ablas::cmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, false, 2, a, offs+n1, offs);
                ablas::cmatrixsyrk<Precision>(n2, n1, -amp::ampf<Precision>("1.0"), a, offs+n1, offs, 0, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1, isupper);
            }
            result = hpdmatrixcholeskyrec<Precision>(a, offs+n1, n2, isupper, tmp);
            if( !result )
            {
                return result;
            }
        }
        return result;
    }


    /*************************************************************************
    Recursive computational subroutine for SPDMatrixCholesky

      -- ALGLIB routine --
         15.12.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool spdmatrixcholeskyrec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        bool result;
        int n1;
        int n2;


        
        //
        // check N
        //
        if( n<1 )
        {
            result = false;
            return result;
        }
        
        //
        // special cases
        //
        if( n==1 )
        {
            if( a(offs,offs)>0 )
            {
                a(offs,offs) = amp::sqrt<Precision>(a(offs,offs));
                result = true;
            }
            else
            {
                result = false;
            }
            return result;
        }
        if( n<=ablas::ablasblocksize<Precision>(a) )
        {
            result = spdmatrixcholesky2<Precision>(a, offs, n, isupper, tmp);
            return result;
        }
        
        //
        // general case: split task in cache-oblivious manner
        //
        result = true;
        ablas::ablassplitlength<Precision>(a, n, n1, n2);
        result = spdmatrixcholeskyrec<Precision>(a, offs, n1, isupper, tmp);
        if( !result )
        {
            return result;
        }
        if( n2>0 )
        {
            if( isupper )
            {
                ablas::rmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, false, 1, a, offs, offs+n1);
                ablas::rmatrixsyrk<Precision>(n2, n1, -amp::ampf<Precision>("1.0"), a, offs, offs+n1, 1, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1, isupper);
            }
            else
            {
                ablas::rmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, false, 1, a, offs+n1, offs);
                ablas::rmatrixsyrk<Precision>(n2, n1, -amp::ampf<Precision>("1.0"), a, offs+n1, offs, 0, +amp::ampf<Precision>("1.0"), a, offs+n1, offs+n1, isupper);
            }
            result = spdmatrixcholeskyrec<Precision>(a, offs+n1, n2, isupper, tmp);
            if( !result )
            {
                return result;
            }
        }
        return result;
    }


    /*************************************************************************
    Level-2 Hermitian Cholesky subroutine.

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    bool hpdmatrixcholesky2(ap::template_2d_array< amp::campf<Precision> >& aaa,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        bool result;
        int i;
        int j;
        int k;
        int j1;
        int j2;
        amp::ampf<Precision> ajj;
        amp::campf<Precision> v;
        amp::ampf<Precision> r;
        int i_;
        int i1_;


        result = true;
        if( n<0 )
        {
            result = false;
            return result;
        }
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        if( isupper )
        {
            
            //
            // Compute the Cholesky factorization A = U'*U.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute U(J,J) and test for non-positive-definiteness.
                //
                v = 0.0;
                for(i_=offs; i_<=offs+j-1;i_++)
                {
                    v += amp::conj(aaa(i_,offs+j))*aaa(i_,offs+j);
                }
                ajj = (aaa(offs+j,offs+j)-v).x;
                if( ajj<=0 )
                {
                    aaa(offs+j,offs+j) = ajj;
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                aaa(offs+j,offs+j) = ajj;
                
                //
                // Compute elements J+1:N-1 of row J.
                //
                if( j<n-1 )
                {
                    if( j>0 )
                    {
                        i1_ = (offs) - (0);
                        for(i_=0; i_<=j-1;i_++)
                        {
                            tmp(i_) = -amp::conj(aaa(i_+i1_,offs+j));
                        }
                        ablas::cmatrixmv<Precision>(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n);
                        i1_ = (n) - (offs+j+1);
                        for(i_=offs+j+1; i_<=offs+n-1;i_++)
                        {
                            aaa(offs+j,i_) = aaa(offs+j,i_) + tmp(i_+i1_);
                        }
                    }
                    r = 1/ajj;
                    for(i_=offs+j+1; i_<=offs+n-1;i_++)
                    {
                        aaa(offs+j,i_) = r*aaa(offs+j,i_);
                    }
                }
            }
        }
        else
        {
            
            //
            // Compute the Cholesky factorization A = L*L'.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute L(J+1,J+1) and test for non-positive-definiteness.
                //
                v = 0.0;
                for(i_=offs; i_<=offs+j-1;i_++)
                {
                    v += amp::conj(aaa(offs+j,i_))*aaa(offs+j,i_);
                }
                ajj = (aaa(offs+j,offs+j)-v).x;
                if( ajj<=0 )
                {
                    aaa(offs+j,offs+j) = ajj;
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                aaa(offs+j,offs+j) = ajj;
                
                //
                // Compute elements J+1:N of column J.
                //
                if( j<n-1 )
                {
                    if( j>0 )
                    {
                        i1_ = (offs) - (0);
                        for(i_=0; i_<=j-1;i_++)
                        {
                            tmp(i_) = amp::conj(aaa(offs+j,i_+i1_));
                        }
                        ablas::cmatrixmv<Precision>(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n);
                        for(i=0; i<=n-j-2; i++)
                        {
                            aaa(offs+j+1+i,offs+j) = (aaa(offs+j+1+i,offs+j)-tmp(n+i))/ajj;
                        }
                    }
                    else
                    {
                        for(i=0; i<=n-j-2; i++)
                        {
                            aaa(offs+j+1+i,offs+j) = aaa(offs+j+1+i,offs+j)/ajj;
                        }
                    }
                }
            }
        }
        return result;
    }


    /*************************************************************************
    Level-2 Cholesky subroutine

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992
    *************************************************************************/
    template<unsigned int Precision>
    bool spdmatrixcholesky2(ap::template_2d_array< amp::ampf<Precision> >& aaa,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        bool result;
        int i;
        int j;
        int k;
        int j1;
        int j2;
        amp::ampf<Precision> ajj;
        amp::ampf<Precision> v;
        amp::ampf<Precision> r;


        result = true;
        if( n<0 )
        {
            result = false;
            return result;
        }
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        if( isupper )
        {
            
            //
            // Compute the Cholesky factorization A = U'*U.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute U(J,J) and test for non-positive-definiteness.
                //
                v = amp::vdotproduct(aaa.getcolumn(offs+j, offs, offs+j-1), aaa.getcolumn(offs+j, offs, offs+j-1));
                ajj = aaa(offs+j,offs+j)-v;
                if( ajj<=0 )
                {
                    aaa(offs+j,offs+j) = ajj;
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                aaa(offs+j,offs+j) = ajj;
                
                //
                // Compute elements J+1:N-1 of row J.
                //
                if( j<n-1 )
                {
                    if( j>0 )
                    {
                        amp::vmoveneg(tmp.getvector(0, j-1), aaa.getcolumn(offs+j, offs, offs+j-1));
                        ablas::rmatrixmv<Precision>(n-j-1, j, aaa, offs, offs+j+1, 1, tmp, 0, tmp, n);
                        amp::vadd(aaa.getrow(offs+j, offs+j+1, offs+n-1), tmp.getvector(n, 2*n-j-2));
                    }
                    r = 1/ajj;
                    amp::vmul(aaa.getrow(offs+j, offs+j+1, offs+n-1), r);
                }
            }
        }
        else
        {
            
            //
            // Compute the Cholesky factorization A = L*L'.
            //
            for(j=0; j<=n-1; j++)
            {
                
                //
                // Compute L(J+1,J+1) and test for non-positive-definiteness.
                //
                v = amp::vdotproduct(aaa.getrow(offs+j, offs, offs+j-1), aaa.getrow(offs+j, offs, offs+j-1));
                ajj = aaa(offs+j,offs+j)-v;
                if( ajj<=0 )
                {
                    aaa(offs+j,offs+j) = ajj;
                    result = false;
                    return result;
                }
                ajj = amp::sqrt<Precision>(ajj);
                aaa(offs+j,offs+j) = ajj;
                
                //
                // Compute elements J+1:N of column J.
                //
                if( j<n-1 )
                {
                    if( j>0 )
                    {
                        amp::vmove(tmp.getvector(0, j-1), aaa.getrow(offs+j, offs, offs+j-1));
                        ablas::rmatrixmv<Precision>(n-j-1, j, aaa, offs+j+1, offs, 0, tmp, 0, tmp, n);
                        for(i=0; i<=n-j-2; i++)
                        {
                            aaa(offs+j+1+i,offs+j) = (aaa(offs+j+1+i,offs+j)-tmp(n+i))/ajj;
                        }
                    }
                    else
                    {
                        for(i=0; i<=n-j-2; i++)
                        {
                            aaa(offs+j+1+i,offs+j) = aaa(offs+j+1+i,offs+j)/ajj;
                        }
                    }
                }
            }
        }
        return result;
    }
} // namespace

#endif
