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

#ifndef _matinv_h
#define _matinv_h

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
#include "rcond.h"
namespace matinv
{
    template<unsigned int Precision>
    class matinvreport
    {
    public:
        amp::ampf<Precision> r1;
        amp::ampf<Precision> rinf;
    };




    template<unsigned int Precision>
    void rmatrixluinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void rmatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void cmatrixluinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void cmatrixinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void spdmatrixcholeskyinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void spdmatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void hpdmatrixcholeskyinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void hpdmatrixinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void rmatrixtrinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void cmatrixtrinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void rmatrixtrinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        bool isunit,
        ap::template_1d_array< amp::ampf<Precision> >& tmp,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void cmatrixtrinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        bool isunit,
        ap::template_1d_array< amp::campf<Precision> >& tmp,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void rmatrixluinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& work,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void cmatrixluinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& work,
        int& info,
        matinvreport<Precision>& rep);
    template<unsigned int Precision>
    void spdmatrixcholeskyinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void hpdmatrixcholeskyinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp);


    /*************************************************************************
    Inversion of a matrix given by its LU decomposition.

    INPUT PARAMETERS:
        A       -   LU decomposition of the matrix (output of RMatrixLU subroutine).
        Pivots  -   table of permutations which were made during the LU decomposition
                    (the output of RMatrixLU subroutine).
        N       -   size of matrix A.

    OUTPUT PARAMETERS:
        Info    -   return code:
                    * -3    A is singular, or VERY close to singular.
                            it is filled by zeros in such cases.
                    * -1    N<=0 was passed, or incorrect Pivots was passed
                    *  1    task is solved (but matrix A may be ill-conditioned,
                            check R1/RInf parameters for condition numbers).
        Rep     -   solver report, see below for more info
        A       -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].

    SOLVER REPORT

    Subroutine sets following fields of the Rep structure:
    * R1        reciprocal of condition number: 1/cond(A), 1-norm.
    * RInf      reciprocal of condition number: 1/cond(A), inf-norm.

      -- ALGLIB routine --
         05.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixluinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        int& info,
        matinvreport<Precision>& rep)
    {
        ap::template_1d_array< amp::ampf<Precision> > work;
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;


        info = 1;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            info = -1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            if( pivots(i)>n-1 || pivots(i)<i )
            {
                info = -1;
                return;
            }
        }
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::rmatrixlurcond1<Precision>(a, n);
        rep.rinf = rcond::rmatrixlurcondinf<Precision>(a, n);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Call cache-oblivious code
        //
        work.setlength(n);
        rmatrixluinverserec<Precision>(a, 0, n, work, info, rep);
        
        //
        // apply permutations
        //
        for(i=0; i<=n-1; i++)
        {
            for(j=n-2; j>=0; j--)
            {
                k = pivots(j);
                v = a(i,j);
                a(i,j) = a(i,k);
                a(i,k) = v;
            }
        }
    }


    /*************************************************************************
    Inversion of a general matrix.

    Input parameters:
        A   -   matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N   -   size of matrix A.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

    Result:
        True, if the matrix is not singular.
        False, if the matrix is singular.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        int& info,
        matinvreport<Precision>& rep)
    {
        ap::template_1d_array< int > pivots;


        trfac::rmatrixlu<Precision>(a, n, n, pivots);
        rmatrixluinverse<Precision>(a, pivots, n, info, rep);
    }


    /*************************************************************************
    Inversion of a matrix given by its LU decomposition.

    INPUT PARAMETERS:
        A       -   LU decomposition of the matrix (output of CMatrixLU subroutine).
        Pivots  -   table of permutations which were made during the LU decomposition
                    (the output of CMatrixLU subroutine).
        N       -   size of matrix A.

    OUTPUT PARAMETERS:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB routine --
         05.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixluinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        int& info,
        matinvreport<Precision>& rep)
    {
        ap::template_1d_array< amp::campf<Precision> > work;
        int i;
        int j;
        int k;
        amp::campf<Precision> v;


        info = 1;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            info = -1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            if( pivots(i)>n-1 || pivots(i)<i )
            {
                info = -1;
                return;
            }
        }
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::cmatrixlurcond1<Precision>(a, n);
        rep.rinf = rcond::cmatrixlurcondinf<Precision>(a, n);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Call cache-oblivious code
        //
        work.setlength(n);
        cmatrixluinverserec<Precision>(a, 0, n, work, info, rep);
        
        //
        // apply permutations
        //
        for(i=0; i<=n-1; i++)
        {
            for(j=n-2; j>=0; j--)
            {
                k = pivots(j);
                v = a(i,j);
                a(i,j) = a(i,k);
                a(i,k) = v;
            }
        }
    }


    /*************************************************************************
    Inversion of a general matrix.

    Input parameters:
        A   -   matrix, array[0..N-1,0..N-1].
        N   -   size of A.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        int& info,
        matinvreport<Precision>& rep)
    {
        ap::template_1d_array< int > pivots;


        trfac::cmatrixlu<Precision>(a, n, n, pivots);
        cmatrixluinverse<Precision>(a, pivots, n, info, rep);
    }


    /*************************************************************************
    Inversion of a symmetric positive definite matrix which is given
    by Cholesky decomposition.

    Input parameters:
        A       -   Cholesky decomposition of the matrix to be inverted:
                    A=U’*U or A = L*L'.
                    Output of  SPDMatrixCholesky subroutine.
        N       -   size of matrix A.
        IsUpper –   storage format.
                    If IsUpper = True, then matrix A is given as A = U'*U
                    (matrix contains upper triangle).
                    Similarly, if IsUpper = False, then A = L*L'.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixcholeskyinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int j;
        int k;
        amp::ampf<Precision> v;
        amp::ampf<Precision> ajj;
        amp::ampf<Precision> aii;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        int info2;
        matinvreport<Precision> rep2;


        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::spdmatrixcholeskyrcond<Precision>(a, n, isupper);
        rep.rinf = rep.r1;
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            if( isupper )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=i; j<=n-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=i; j++)
                    {
                        a(i,j) = 0;
                    }
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Inverse
        //
        tmp.setlength(n);
        spdmatrixcholeskyinverserec<Precision>(a, 0, n, isupper, tmp);
    }


    /*************************************************************************
    Inversion of a symmetric positive definite matrix.

    Given an upper or lower triangle of a symmetric positive definite matrix,
    the algorithm generates matrix A^-1 and saves the upper or lower triangle
    depending on the input.

    Input parameters:
        A       -   matrix to be inverted (upper or lower triangle).
                    Array with elements [0..N-1,0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.
                    If IsUpper = True, then the upper triangle of matrix A is
                    given, otherwise the lower triangle is given.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep)
    {
        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        if( trfac::spdmatrixcholesky<Precision>(a, n, isupper) )
        {
            spdmatrixcholeskyinverse<Precision>(a, n, isupper, info, rep);
        }
        else
        {
            info = -3;
        }
    }


    /*************************************************************************
    Inversion of a Hermitian positive definite matrix which is given
    by Cholesky decomposition.

    Input parameters:
        A       -   Cholesky decomposition of the matrix to be inverted:
                    A=U’*U or A = L*L'.
                    Output of  HPDMatrixCholesky subroutine.
        N       -   size of matrix A.
        IsUpper –   storage format.
                    If IsUpper = True, then matrix A is given as A = U'*U
                    (matrix contains upper triangle).
                    Similarly, if IsUpper = False, then A = L*L'.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixcholeskyinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int j;
        int info2;
        matinvreport<Precision> rep2;
        ap::template_1d_array< amp::campf<Precision> > tmp;
        amp::campf<Precision> v;


        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::hpdmatrixcholeskyrcond<Precision>(a, n, isupper);
        rep.rinf = rep.r1;
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            if( isupper )
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=i; j<=n-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=i; j++)
                    {
                        a(i,j) = 0;
                    }
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Inverse
        //
        tmp.setlength(n);
        hpdmatrixcholeskyinverserec<Precision>(a, 0, n, isupper, tmp);
    }


    /*************************************************************************
    Inversion of a Hermitian positive definite matrix.

    Given an upper or lower triangle of a Hermitian positive definite matrix,
    the algorithm generates matrix A^-1 and saves the upper or lower triangle
    depending on the input.

    Input parameters:
        A       -   matrix to be inverted (upper or lower triangle).
                    Array with elements [0..N-1,0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format.
                    If IsUpper = True, then the upper triangle of matrix A is
                    given, otherwise the lower triangle is given.

    Output parameters:
        Info    -   return code, same as in RMatrixLUInverse
        Rep     -   solver report, same as in RMatrixLUInverse
        A       -   inverse of matrix A, same as in RMatrixLUInverse

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        int& info,
        matinvreport<Precision>& rep)
    {
        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        if( trfac::hpdmatrixcholesky<Precision>(a, n, isupper) )
        {
            hpdmatrixcholeskyinverse<Precision>(a, n, isupper, info, rep);
        }
        else
        {
            info = -3;
        }
    }


    /*************************************************************************
    Triangular matrix inverse (real)

    The subroutine inverts the following types of matrices:
        * upper triangular
        * upper triangular with unit diagonal
        * lower triangular
        * lower triangular with unit diagonal

    In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
    also be upper (lower) triangular, and after the end of the algorithm,  the
    inverse matrix replaces the source matrix. The elements  below (above) the
    main diagonal are not changed by the algorithm.

    If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
    diagonal, and the diagonal elements are not passed to the algorithm.

    Input parameters:
        A       -   matrix, array[0..N-1, 0..N-1].
        N       -   size of A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Output parameters:
        Info    -   same as for RMatrixLUInverse
        Rep     -   same as for RMatrixLUInverse
        A       -   same as for RMatrixLUInverse.

      -- ALGLIB --
         Copyright 05.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixtrinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int j;
        ap::template_1d_array< amp::ampf<Precision> > tmp;


        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::rmatrixtrrcond1<Precision>(a, n, isupper, isunit);
        rep.rinf = rcond::rmatrixtrrcondinf<Precision>(a, n, isupper, isunit);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Invert
        //
        tmp.setlength(n);
        rmatrixtrinverserec<Precision>(a, 0, n, isupper, isunit, tmp, info, rep);
    }


    /*************************************************************************
    Triangular matrix inverse (complex)

    The subroutine inverts the following types of matrices:
        * upper triangular
        * upper triangular with unit diagonal
        * lower triangular
        * lower triangular with unit diagonal

    In case of an upper (lower) triangular matrix,  the  inverse  matrix  will
    also be upper (lower) triangular, and after the end of the algorithm,  the
    inverse matrix replaces the source matrix. The elements  below (above) the
    main diagonal are not changed by the algorithm.

    If  the matrix  has a unit diagonal, the inverse matrix also  has  a  unit
    diagonal, and the diagonal elements are not passed to the algorithm.

    Input parameters:
        A       -   matrix, array[0..N-1, 0..N-1].
        N       -   size of A.
        IsUpper -   True, if the matrix is upper triangular.
        IsUnit  -   True, if the matrix has a unit diagonal.

    Output parameters:
        Info    -   same as for RMatrixLUInverse
        Rep     -   same as for RMatrixLUInverse
        A       -   same as for RMatrixLUInverse.

      -- ALGLIB --
         Copyright 05.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixtrinverse(ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        bool isunit,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int j;
        ap::template_1d_array< amp::campf<Precision> > tmp;


        if( n<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        
        //
        // calculate condition numbers
        //
        rep.r1 = rcond::cmatrixtrrcond1<Precision>(a, n, isupper, isunit);
        rep.rinf = rcond::cmatrixtrrcondinf<Precision>(a, n, isupper, isunit);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    a(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        
        //
        // Invert
        //
        tmp.setlength(n);
        cmatrixtrinverserec<Precision>(a, 0, n, isupper, isunit, tmp, info, rep);
    }


    /*************************************************************************
    Triangular matrix inversion, recursive subroutine

      -- ALGLIB --
         05.02.2010, Bochkanov Sergey.
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992.
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixtrinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        bool isunit,
        ap::template_1d_array< amp::ampf<Precision> >& tmp,
        int& info,
        matinvreport<Precision>& rep)
    {
        int n1;
        int n2;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> ajj;


        if( n<1 )
        {
            info = -1;
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablasblocksize<Precision>(a) )
        {
            if( isupper )
            {
                
                //
                // Compute inverse of upper triangular matrix.
                //
                for(j=0; j<=n-1; j++)
                {
                    if( !isunit )
                    {
                        if( a(offs+j,offs+j)==0 )
                        {
                            info = -3;
                            return;
                        }
                        a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                        ajj = -a(offs+j,offs+j);
                    }
                    else
                    {
                        ajj = -1;
                    }
                    
                    //
                    // Compute elements 1:j-1 of j-th column.
                    //
                    if( j>0 )
                    {
                        amp::vmove(tmp.getvector(0, j-1), a.getcolumn(offs+j, offs+0, offs+j-1));
                        for(i=0; i<=j-1; i++)
                        {
                            if( i<j-1 )
                            {
                                v = amp::vdotproduct(a.getrow(offs+i, offs+i+1, offs+j-1), tmp.getvector(i+1, j-1));
                            }
                            else
                            {
                                v = 0;
                            }
                            if( !isunit )
                            {
                                a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                            }
                            else
                            {
                                a(offs+i,offs+j) = v+tmp(i);
                            }
                        }
                        amp::vmul(a.getcolumn(offs+j, offs+0, offs+j-1), ajj);
                    }
                }
            }
            else
            {
                
                //
                // Compute inverse of lower triangular matrix.
                //
                for(j=n-1; j>=0; j--)
                {
                    if( !isunit )
                    {
                        if( a(offs+j,offs+j)==0 )
                        {
                            info = -3;
                            return;
                        }
                        a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                        ajj = -a(offs+j,offs+j);
                    }
                    else
                    {
                        ajj = -1;
                    }
                    if( j<n-1 )
                    {
                        
                        //
                        // Compute elements j+1:n of j-th column.
                        //
                        amp::vmove(tmp.getvector(j+1, n-1), a.getcolumn(offs+j, offs+j+1, offs+n-1));
                        for(i=j+1; i<=n-1; i++)
                        {
                            if( i>j+1 )
                            {
                                v = amp::vdotproduct(a.getrow(offs+i, offs+j+1, offs+i-1), tmp.getvector(j+1, i-1));
                            }
                            else
                            {
                                v = 0;
                            }
                            if( !isunit )
                            {
                                a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                            }
                            else
                            {
                                a(offs+i,offs+j) = v+tmp(i);
                            }
                        }
                        amp::vmul(a.getcolumn(offs+j, offs+j+1, offs+n-1), ajj);
                    }
                }
            }
            return;
        }
        
        //
        // Recursive case
        //
        ablas::ablassplitlength<Precision>(a, n, n1, n2);
        if( n2>0 )
        {
            if( isupper )
            {
                for(i=0; i<=n1-1; i++)
                {
                    amp::vmul(a.getrow(offs+i, offs+n1, offs+n-1), -1);
                }
                ablas::rmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs+n1);
                ablas::rmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs, offs+n1);
            }
            else
            {
                for(i=0; i<=n2-1; i++)
                {
                    amp::vmul(a.getrow(offs+n1+i, offs, offs+n1-1), -1);
                }
                ablas::rmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs+n1, offs);
                ablas::rmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs+n1, offs);
            }
            rmatrixtrinverserec<Precision>(a, offs+n1, n2, isupper, isunit, tmp, info, rep);
        }
        rmatrixtrinverserec<Precision>(a, offs, n1, isupper, isunit, tmp, info, rep);
    }


    /*************************************************************************
    Triangular matrix inversion, recursive subroutine

      -- ALGLIB --
         05.02.2010, Bochkanov Sergey.
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         February 29, 1992.
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixtrinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        bool isunit,
        ap::template_1d_array< amp::campf<Precision> >& tmp,
        int& info,
        matinvreport<Precision>& rep)
    {
        int n1;
        int n2;
        int i;
        int j;
        amp::campf<Precision> v;
        amp::campf<Precision> ajj;
        int i_;
        int i1_;


        if( n<1 )
        {
            info = -1;
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            if( isupper )
            {
                
                //
                // Compute inverse of upper triangular matrix.
                //
                for(j=0; j<=n-1; j++)
                {
                    if( !isunit )
                    {
                        if( a(offs+j,offs+j)==0 )
                        {
                            info = -3;
                            return;
                        }
                        a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                        ajj = -a(offs+j,offs+j);
                    }
                    else
                    {
                        ajj = -1;
                    }
                    
                    //
                    // Compute elements 1:j-1 of j-th column.
                    //
                    if( j>0 )
                    {
                        i1_ = (offs+0) - (0);
                        for(i_=0; i_<=j-1;i_++)
                        {
                            tmp(i_) = a(i_+i1_,offs+j);
                        }
                        for(i=0; i<=j-1; i++)
                        {
                            if( i<j-1 )
                            {
                                i1_ = (i+1)-(offs+i+1);
                                v = 0.0;
                                for(i_=offs+i+1; i_<=offs+j-1;i_++)
                                {
                                    v += a(offs+i,i_)*tmp(i_+i1_);
                                }
                            }
                            else
                            {
                                v = 0;
                            }
                            if( !isunit )
                            {
                                a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                            }
                            else
                            {
                                a(offs+i,offs+j) = v+tmp(i);
                            }
                        }
                        for(i_=offs+0; i_<=offs+j-1;i_++)
                        {
                            a(i_,offs+j) = ajj*a(i_,offs+j);
                        }
                    }
                }
            }
            else
            {
                
                //
                // Compute inverse of lower triangular matrix.
                //
                for(j=n-1; j>=0; j--)
                {
                    if( !isunit )
                    {
                        if( a(offs+j,offs+j)==0 )
                        {
                            info = -3;
                            return;
                        }
                        a(offs+j,offs+j) = 1/a(offs+j,offs+j);
                        ajj = -a(offs+j,offs+j);
                    }
                    else
                    {
                        ajj = -1;
                    }
                    if( j<n-1 )
                    {
                        
                        //
                        // Compute elements j+1:n of j-th column.
                        //
                        i1_ = (offs+j+1) - (j+1);
                        for(i_=j+1; i_<=n-1;i_++)
                        {
                            tmp(i_) = a(i_+i1_,offs+j);
                        }
                        for(i=j+1; i<=n-1; i++)
                        {
                            if( i>j+1 )
                            {
                                i1_ = (j+1)-(offs+j+1);
                                v = 0.0;
                                for(i_=offs+j+1; i_<=offs+i-1;i_++)
                                {
                                    v += a(offs+i,i_)*tmp(i_+i1_);
                                }
                            }
                            else
                            {
                                v = 0;
                            }
                            if( !isunit )
                            {
                                a(offs+i,offs+j) = v+a(offs+i,offs+i)*tmp(i);
                            }
                            else
                            {
                                a(offs+i,offs+j) = v+tmp(i);
                            }
                        }
                        for(i_=offs+j+1; i_<=offs+n-1;i_++)
                        {
                            a(i_,offs+j) = ajj*a(i_,offs+j);
                        }
                    }
                }
            }
            return;
        }
        
        //
        // Recursive case
        //
        ablas::ablascomplexsplitlength<Precision>(a, n, n1, n2);
        if( n2>0 )
        {
            if( isupper )
            {
                for(i=0; i<=n1-1; i++)
                {
                    for(i_=offs+n1; i_<=offs+n-1;i_++)
                    {
                        a(offs+i,i_) = -1*a(offs+i,i_);
                    }
                }
                ablas::cmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, isunit, 0, a, offs, offs+n1);
                ablas::cmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs, offs+n1);
            }
            else
            {
                for(i=0; i<=n2-1; i++)
                {
                    for(i_=offs; i_<=offs+n1-1;i_++)
                    {
                        a(offs+n1+i,i_) = -1*a(offs+n1+i,i_);
                    }
                }
                ablas::cmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, isunit, 0, a, offs+n1, offs);
                ablas::cmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, isunit, 0, a, offs+n1, offs);
            }
            cmatrixtrinverserec<Precision>(a, offs+n1, n2, isupper, isunit, tmp, info, rep);
        }
        cmatrixtrinverserec<Precision>(a, offs, n1, isupper, isunit, tmp, info, rep);
    }


    template<unsigned int Precision>
    void rmatrixluinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& work,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int iws;
        int j;
        int jb;
        int jj;
        int jp;
        int k;
        amp::ampf<Precision> v;
        int n1;
        int n2;


        if( n<1 )
        {
            info = -1;
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablasblocksize<Precision>(a) )
        {
            
            //
            // Form inv(U)
            //
            rmatrixtrinverserec<Precision>(a, offs, n, true, false, work, info, rep);
            if( info<=0 )
            {
                return;
            }
            
            //
            // Solve the equation inv(A)*L = inv(U) for inv(A).
            //
            for(j=n-1; j>=0; j--)
            {
                
                //
                // Copy current column of L to WORK and replace with zeros.
                //
                for(i=j+1; i<=n-1; i++)
                {
                    work(i) = a(offs+i,offs+j);
                    a(offs+i,offs+j) = 0;
                }
                
                //
                // Compute current column of inv(A).
                //
                if( j<n-1 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        v = amp::vdotproduct(a.getrow(offs+i, offs+j+1, offs+n-1), work.getvector(j+1, n-1));
                        a(offs+i,offs+j) = a(offs+i,offs+j)-v;
                    }
                }
            }
            return;
        }
        
        //
        // Recursive code:
        //
        //         ( L1      )   ( U1  U12 )
        // A    =  (         ) * (         )
        //         ( L12  L2 )   (     U2  )
        //
        //         ( W   X )
        // A^-1 =  (       )
        //         ( Y   Z )
        //
        ablas::ablassplitlength<Precision>(a, n, n1, n2);
        ap::ap_error::make_assertion(n2>0);
        
        //
        // X := inv(U1)*U12*inv(U2)
        //
        ablas::rmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, true, false, 0, a, offs, offs+n1);
        ablas::rmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, true, false, 0, a, offs, offs+n1);
        
        //
        // Y := inv(L2)*L12*inv(L1)
        //
        ablas::rmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, false, true, 0, a, offs+n1, offs);
        ablas::rmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, false, true, 0, a, offs+n1, offs);
        
        //
        // W := inv(L1*U1)+X*Y
        //
        rmatrixluinverserec<Precision>(a, offs, n1, work, info, rep);
        if( info<=0 )
        {
            return;
        }
        ablas::rmatrixgemm<Precision>(n1, n1, n2, amp::ampf<Precision>("1.0"), a, offs, offs+n1, 0, a, offs+n1, offs, 0, amp::ampf<Precision>("1.0"), a, offs, offs);
        
        //
        // X := -X*inv(L2)
        // Y := -inv(U2)*Y
        //
        ablas::rmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, false, true, 0, a, offs, offs+n1);
        for(i=0; i<=n1-1; i++)
        {
            amp::vmul(a.getrow(offs+i, offs+n1, offs+n-1), -1);
        }
        ablas::rmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, true, false, 0, a, offs+n1, offs);
        for(i=0; i<=n2-1; i++)
        {
            amp::vmul(a.getrow(offs+n1+i, offs, offs+n1-1), -1);
        }
        
        //
        // Z := inv(L2*U2)
        //
        rmatrixluinverserec<Precision>(a, offs+n1, n2, work, info, rep);
    }


    template<unsigned int Precision>
    void cmatrixluinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& work,
        int& info,
        matinvreport<Precision>& rep)
    {
        int i;
        int iws;
        int j;
        int jb;
        int jj;
        int jp;
        int k;
        amp::campf<Precision> v;
        int n1;
        int n2;
        int i_;
        int i1_;


        if( n<1 )
        {
            info = -1;
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            
            //
            // Form inv(U)
            //
            cmatrixtrinverserec<Precision>(a, offs, n, true, false, work, info, rep);
            if( info<=0 )
            {
                return;
            }
            
            //
            // Solve the equation inv(A)*L = inv(U) for inv(A).
            //
            for(j=n-1; j>=0; j--)
            {
                
                //
                // Copy current column of L to WORK and replace with zeros.
                //
                for(i=j+1; i<=n-1; i++)
                {
                    work(i) = a(offs+i,offs+j);
                    a(offs+i,offs+j) = 0;
                }
                
                //
                // Compute current column of inv(A).
                //
                if( j<n-1 )
                {
                    for(i=0; i<=n-1; i++)
                    {
                        i1_ = (j+1)-(offs+j+1);
                        v = 0.0;
                        for(i_=offs+j+1; i_<=offs+n-1;i_++)
                        {
                            v += a(offs+i,i_)*work(i_+i1_);
                        }
                        a(offs+i,offs+j) = a(offs+i,offs+j)-v;
                    }
                }
            }
            return;
        }
        
        //
        // Recursive code:
        //
        //         ( L1      )   ( U1  U12 )
        // A    =  (         ) * (         )
        //         ( L12  L2 )   (     U2  )
        //
        //         ( W   X )
        // A^-1 =  (       )
        //         ( Y   Z )
        //
        ablas::ablascomplexsplitlength<Precision>(a, n, n1, n2);
        ap::ap_error::make_assertion(n2>0);
        
        //
        // X := inv(U1)*U12*inv(U2)
        //
        ablas::cmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, true, false, 0, a, offs, offs+n1);
        ablas::cmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, true, false, 0, a, offs, offs+n1);
        
        //
        // Y := inv(L2)*L12*inv(L1)
        //
        ablas::cmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, false, true, 0, a, offs+n1, offs);
        ablas::cmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, false, true, 0, a, offs+n1, offs);
        
        //
        // W := inv(L1*U1)+X*Y
        //
        cmatrixluinverserec<Precision>(a, offs, n1, work, info, rep);
        if( info<=0 )
        {
            return;
        }
        ablas::cmatrixgemm<Precision>(n1, n1, n2, amp::ampf<Precision>("1.0"), a, offs, offs+n1, 0, a, offs+n1, offs, 0, amp::ampf<Precision>("1.0"), a, offs, offs);
        
        //
        // X := -X*inv(L2)
        // Y := -inv(U2)*Y
        //
        ablas::cmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, false, true, 0, a, offs, offs+n1);
        for(i=0; i<=n1-1; i++)
        {
            for(i_=offs+n1; i_<=offs+n-1;i_++)
            {
                a(offs+i,i_) = -1*a(offs+i,i_);
            }
        }
        ablas::cmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, true, false, 0, a, offs+n1, offs);
        for(i=0; i<=n2-1; i++)
        {
            for(i_=offs; i_<=offs+n1-1;i_++)
            {
                a(offs+n1+i,i_) = -1*a(offs+n1+i,i_);
            }
        }
        
        //
        // Z := inv(L2*U2)
        //
        cmatrixluinverserec<Precision>(a, offs+n1, n2, work, info, rep);
    }


    /*************************************************************************
    Recursive subroutine for SPD inversion.

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixcholeskyinverserec(ap::template_2d_array< amp::ampf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        int j;
        amp::ampf<Precision> v;
        int n1;
        int n2;
        int info2;
        matinvreport<Precision> rep2;


        if( n<1 )
        {
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablasblocksize<Precision>(a) )
        {
            rmatrixtrinverserec<Precision>(a, offs, n, isupper, false, tmp, info2, rep2);
            if( isupper )
            {
                
                //
                // Compute the product U * U'.
                // NOTE: we never assume that diagonal of U is real
                //
                for(i=0; i<=n-1; i++)
                {
                    if( i==0 )
                    {
                        
                        //
                        // 1x1 matrix
                        //
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i));
                    }
                    else
                    {
                        
                        //
                        // (I+1)x(I+1) matrix,
                        //
                        // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
                        // (          ) * (              ) = (                                )
                        // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
                        //
                        // A11 is IxI, A22 is 1x1.
                        //
                        amp::vmove(tmp.getvector(0, i-1), a.getcolumn(offs+i, offs, offs+i-1));
                        for(j=0; j<=i-1; j++)
                        {
                            v = a(offs+j,offs+i);
                            amp::vadd(a.getrow(offs+j, offs+j, offs+i-1), tmp.getvector(j, i-1), v);
                        }
                        v = a(offs+i,offs+i);
                        amp::vmul(a.getcolumn(offs+i, offs, offs+i-1), v);
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i));
                    }
                }
            }
            else
            {
                
                //
                // Compute the product L' * L
                // NOTE: we never assume that diagonal of L is real
                //
                for(i=0; i<=n-1; i++)
                {
                    if( i==0 )
                    {
                        
                        //
                        // 1x1 matrix
                        //
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i));
                    }
                    else
                    {
                        
                        //
                        // (I+1)x(I+1) matrix,
                        //
                        // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
                        // (              ) * (          ) = (                                )
                        // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
                        //
                        // A11 is IxI, A22 is 1x1.
                        //
                        amp::vmove(tmp.getvector(0, i-1), a.getrow(offs+i, offs, offs+i-1));
                        for(j=0; j<=i-1; j++)
                        {
                            v = a(offs+i,offs+j);
                            amp::vadd(a.getrow(offs+j, offs, offs+j), tmp.getvector(0, j), v);
                        }
                        v = a(offs+i,offs+i);
                        amp::vmul(a.getrow(offs+i, offs, offs+i-1), v);
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i));
                    }
                }
            }
            return;
        }
        
        //
        // Recursive code: triangular factor inversion merged with
        // UU' or L'L multiplication
        //
        ablas::ablassplitlength<Precision>(a, n, n1, n2);
        
        //
        // form off-diagonal block of trangular inverse
        //
        if( isupper )
        {
            for(i=0; i<=n1-1; i++)
            {
                amp::vmul(a.getrow(offs+i, offs+n1, offs+n-1), -1);
            }
            ablas::rmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs+n1);
            ablas::rmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, false, 0, a, offs, offs+n1);
        }
        else
        {
            for(i=0; i<=n2-1; i++)
            {
                amp::vmul(a.getrow(offs+n1+i, offs, offs+n1-1), -1);
            }
            ablas::rmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, false, 0, a, offs+n1, offs);
            ablas::rmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, false, 0, a, offs+n1, offs);
        }
        
        //
        // invert first diagonal block
        //
        spdmatrixcholeskyinverserec<Precision>(a, offs, n1, isupper, tmp);
        
        //
        // update first diagonal block with off-diagonal block,
        // update off-diagonal block
        //
        if( isupper )
        {
            ablas::rmatrixsyrk<Precision>(n1, n2, amp::ampf<Precision>("1.0"), a, offs, offs+n1, 0, amp::ampf<Precision>("1.0"), a, offs, offs, isupper);
            ablas::rmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, false, 1, a, offs, offs+n1);
        }
        else
        {
            ablas::rmatrixsyrk<Precision>(n1, n2, amp::ampf<Precision>("1.0"), a, offs+n1, offs, 1, amp::ampf<Precision>("1.0"), a, offs, offs, isupper);
            ablas::rmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, false, 1, a, offs+n1, offs);
        }
        
        //
        // invert second diagonal block
        //
        spdmatrixcholeskyinverserec<Precision>(a, offs+n1, n2, isupper, tmp);
    }


    /*************************************************************************
    Recursive subroutine for HPD inversion.

      -- ALGLIB routine --
         10.02.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixcholeskyinverserec(ap::template_2d_array< amp::campf<Precision> >& a,
        int offs,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        int j;
        amp::campf<Precision> v;
        int n1;
        int n2;
        int info2;
        matinvreport<Precision> rep2;
        int i_;
        int i1_;


        if( n<1 )
        {
            return;
        }
        
        //
        // Base case
        //
        if( n<=ablas::ablascomplexblocksize<Precision>(a) )
        {
            cmatrixtrinverserec<Precision>(a, offs, n, isupper, false, tmp, info2, rep2);
            if( isupper )
            {
                
                //
                // Compute the product U * U'.
                // NOTE: we never assume that diagonal of U is real
                //
                for(i=0; i<=n-1; i++)
                {
                    if( i==0 )
                    {
                        
                        //
                        // 1x1 matrix
                        //
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i).x)+amp::sqr<Precision>(a(offs+i,offs+i).y);
                    }
                    else
                    {
                        
                        //
                        // (I+1)x(I+1) matrix,
                        //
                        // ( A11  A12 )   ( A11^H        )   ( A11*A11^H+A12*A12^H  A12*A22^H )
                        // (          ) * (              ) = (                                )
                        // (      A22 )   ( A12^H  A22^H )   ( A22*A12^H            A22*A22^H )
                        //
                        // A11 is IxI, A22 is 1x1.
                        //
                        i1_ = (offs) - (0);
                        for(i_=0; i_<=i-1;i_++)
                        {
                            tmp(i_) = amp::conj(a(i_+i1_,offs+i));
                        }
                        for(j=0; j<=i-1; j++)
                        {
                            v = a(offs+j,offs+i);
                            i1_ = (j) - (offs+j);
                            for(i_=offs+j; i_<=offs+i-1;i_++)
                            {
                                a(offs+j,i_) = a(offs+j,i_) + v*tmp(i_+i1_);
                            }
                        }
                        v = amp::conj<Precision>(a(offs+i,offs+i));
                        for(i_=offs; i_<=offs+i-1;i_++)
                        {
                            a(i_,offs+i) = v*a(i_,offs+i);
                        }
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i).x)+amp::sqr<Precision>(a(offs+i,offs+i).y);
                    }
                }
            }
            else
            {
                
                //
                // Compute the product L' * L
                // NOTE: we never assume that diagonal of L is real
                //
                for(i=0; i<=n-1; i++)
                {
                    if( i==0 )
                    {
                        
                        //
                        // 1x1 matrix
                        //
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i).x)+amp::sqr<Precision>(a(offs+i,offs+i).y);
                    }
                    else
                    {
                        
                        //
                        // (I+1)x(I+1) matrix,
                        //
                        // ( A11^H  A21^H )   ( A11      )   ( A11^H*A11+A21^H*A21  A21^H*A22 )
                        // (              ) * (          ) = (                                )
                        // (        A22^H )   ( A21  A22 )   ( A22^H*A21            A22^H*A22 )
                        //
                        // A11 is IxI, A22 is 1x1.
                        //
                        i1_ = (offs) - (0);
                        for(i_=0; i_<=i-1;i_++)
                        {
                            tmp(i_) = a(offs+i,i_+i1_);
                        }
                        for(j=0; j<=i-1; j++)
                        {
                            v = amp::conj<Precision>(a(offs+i,offs+j));
                            i1_ = (0) - (offs);
                            for(i_=offs; i_<=offs+j;i_++)
                            {
                                a(offs+j,i_) = a(offs+j,i_) + v*tmp(i_+i1_);
                            }
                        }
                        v = amp::conj<Precision>(a(offs+i,offs+i));
                        for(i_=offs; i_<=offs+i-1;i_++)
                        {
                            a(offs+i,i_) = v*a(offs+i,i_);
                        }
                        a(offs+i,offs+i) = amp::sqr<Precision>(a(offs+i,offs+i).x)+amp::sqr<Precision>(a(offs+i,offs+i).y);
                    }
                }
            }
            return;
        }
        
        //
        // Recursive code: triangular factor inversion merged with
        // UU' or L'L multiplication
        //
        ablas::ablascomplexsplitlength<Precision>(a, n, n1, n2);
        
        //
        // form off-diagonal block of trangular inverse
        //
        if( isupper )
        {
            for(i=0; i<=n1-1; i++)
            {
                for(i_=offs+n1; i_<=offs+n-1;i_++)
                {
                    a(offs+i,i_) = -1*a(offs+i,i_);
                }
            }
            ablas::cmatrixlefttrsm<Precision>(n1, n2, a, offs, offs, isupper, false, 0, a, offs, offs+n1);
            ablas::cmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, false, 0, a, offs, offs+n1);
        }
        else
        {
            for(i=0; i<=n2-1; i++)
            {
                for(i_=offs; i_<=offs+n1-1;i_++)
                {
                    a(offs+n1+i,i_) = -1*a(offs+n1+i,i_);
                }
            }
            ablas::cmatrixrighttrsm<Precision>(n2, n1, a, offs, offs, isupper, false, 0, a, offs+n1, offs);
            ablas::cmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, false, 0, a, offs+n1, offs);
        }
        
        //
        // invert first diagonal block
        //
        hpdmatrixcholeskyinverserec<Precision>(a, offs, n1, isupper, tmp);
        
        //
        // update first diagonal block with off-diagonal block,
        // update off-diagonal block
        //
        if( isupper )
        {
            ablas::cmatrixsyrk<Precision>(n1, n2, amp::ampf<Precision>("1.0"), a, offs, offs+n1, 0, amp::ampf<Precision>("1.0"), a, offs, offs, isupper);
            ablas::cmatrixrighttrsm<Precision>(n1, n2, a, offs+n1, offs+n1, isupper, false, 2, a, offs, offs+n1);
        }
        else
        {
            ablas::cmatrixsyrk<Precision>(n1, n2, amp::ampf<Precision>("1.0"), a, offs+n1, offs, 2, amp::ampf<Precision>("1.0"), a, offs, offs, isupper);
            ablas::cmatrixlefttrsm<Precision>(n2, n1, a, offs+n1, offs+n1, isupper, false, 2, a, offs+n1, offs);
        }
        
        //
        // invert second diagonal block
        //
        hpdmatrixcholeskyinverserec<Precision>(a, offs+n1, n2, isupper, tmp);
    }
} // namespace

#endif
