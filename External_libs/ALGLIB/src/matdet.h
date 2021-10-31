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

#ifndef _matdet_h
#define _matdet_h

#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
namespace matdet
{
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixludet(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixdet(ap::template_2d_array< amp::ampf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::campf<Precision> cmatrixludet(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n);
    template<unsigned int Precision>
    amp::campf<Precision> cmatrixdet(ap::template_2d_array< amp::campf<Precision> > a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixcholeskydet(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n);
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixdet(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper);


    /*************************************************************************
    Determinant calculation of the matrix given by its LU decomposition.

    Input parameters:
        A       -   LU decomposition of the matrix (output of
                    RMatrixLU subroutine).
        Pivots  -   table of permutations which were made during
                    the LU decomposition.
                    Output of RMatrixLU subroutine.
        N       -   size of matrix A.

    Result: matrix determinant.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixludet(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n)
    {
        amp::ampf<Precision> result;
        int i;
        int s;


        result = 1;
        s = 1;
        for(i=0; i<=n-1; i++)
        {
            result = result*a(i,i);
            if( pivots(i)!=i )
            {
                s = -s;
            }
        }
        result = result*s;
        return result;
    }


    /*************************************************************************
    Calculation of the determinant of a general matrix

    Input parameters:
        A       -   matrix, array[0..N-1, 0..N-1]
        N       -   size of matrix A.

    Result: determinant of matrix A.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> rmatrixdet(ap::template_2d_array< amp::ampf<Precision> > a,
        int n)
    {
        amp::ampf<Precision> result;
        ap::template_1d_array< int > pivots;


        trfac::rmatrixlu<Precision>(a, n, n, pivots);
        result = rmatrixludet<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    Determinant calculation of the matrix given by its LU decomposition.

    Input parameters:
        A       -   LU decomposition of the matrix (output of
                    RMatrixLU subroutine).
        Pivots  -   table of permutations which were made during
                    the LU decomposition.
                    Output of RMatrixLU subroutine.
        N       -   size of matrix A.

    Result: matrix determinant.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::campf<Precision> cmatrixludet(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n)
    {
        amp::campf<Precision> result;
        int i;
        int s;


        result = 1;
        s = 1;
        for(i=0; i<=n-1; i++)
        {
            result = result*a(i,i);
            if( pivots(i)!=i )
            {
                s = -s;
            }
        }
        result = result*s;
        return result;
    }


    /*************************************************************************
    Calculation of the determinant of a general matrix

    Input parameters:
        A       -   matrix, array[0..N-1, 0..N-1]
        N       -   size of matrix A.

    Result: determinant of matrix A.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::campf<Precision> cmatrixdet(ap::template_2d_array< amp::campf<Precision> > a,
        int n)
    {
        amp::campf<Precision> result;
        ap::template_1d_array< int > pivots;


        trfac::cmatrixlu<Precision>(a, n, n, pivots);
        result = cmatrixludet<Precision>(a, pivots, n);
        return result;
    }


    /*************************************************************************
    Determinant calculation of the matrix given by the Cholesky decomposition.

    Input parameters:
        A   -   Cholesky decomposition,
                output of SMatrixCholesky subroutine.
        N   -   size of matrix A.

    As the determinant is equal to the product of squares of diagonal elements,
    it’s not necessary to specify which triangle - lower or upper - the matrix
    is stored in.

    Result:
        matrix determinant.

      -- ALGLIB --
         Copyright 2005-2008 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixcholeskydet(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n)
    {
        amp::ampf<Precision> result;
        int i;


        result = 1;
        for(i=0; i<=n-1; i++)
        {
            result = result*amp::sqr<Precision>(a(i,i));
        }
        return result;
    }


    /*************************************************************************
    Determinant calculation of the symmetric positive definite matrix.

    Input parameters:
        A       -   matrix. Array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   if IsUpper = True, then the symmetric matrix A is given by
                    its upper triangle, and the lower triangle isn’t used by
                    subroutine. Similarly, if IsUpper = False, then A is given
                    by its lower triangle.

    Result:
        determinant of matrix A.
        If matrix A is not positive definite, then subroutine returns -1.

      -- ALGLIB --
         Copyright 2005-2008 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> spdmatrixdet(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isupper)
    {
        amp::ampf<Precision> result;


        if( !trfac::spdmatrixcholesky<Precision>(a, n, isupper) )
        {
            result = -1;
        }
        else
        {
            result = spdmatrixcholeskydet<Precision>(a, n);
        }
        return result;
    }
} // namespace

#endif
