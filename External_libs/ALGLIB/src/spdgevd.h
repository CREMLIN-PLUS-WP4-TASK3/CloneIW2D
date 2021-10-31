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

#ifndef _spdgevd_h
#define _spdgevd_h

#include "ap.h"
#include "amp.h"
#include "reflections.h"
#include "creflections.h"
#include "hqrnd.h"
#include "matgen.h"
#include "ablasf.h"
#include "ablas.h"
#include "trfac.h"
#include "sblas.h"
#include "blas.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "matinv.h"
#include "hblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "hsschur.h"
#include "evd.h"
namespace spdgevd
{
    template<unsigned int Precision>
    bool smatrixgevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isuppera,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        bool isupperb,
        int zneeded,
        int problemtype,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::ampf<Precision> >& z);
    template<unsigned int Precision>
    bool smatrixgevdreduce(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isuppera,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        bool isupperb,
        int problemtype,
        ap::template_2d_array< amp::ampf<Precision> >& r,
        bool& isupperr);


    /*************************************************************************
    Algorithm for solving the following generalized symmetric positive-definite
    eigenproblem:
        A*x = lambda*B*x (1) or
        A*B*x = lambda*x (2) or
        B*A*x = lambda*x (3).
    where A is a symmetric matrix, B - symmetric positive-definite matrix.
    The problem is solved by reducing it to an ordinary  symmetric  eigenvalue
    problem.

    Input parameters:
        A           -   symmetric matrix which is given by its upper or lower
                        triangular part.
                        Array whose indexes range within [0..N-1, 0..N-1].
        N           -   size of matrices A and B.
        IsUpperA    -   storage format of matrix A.
        B           -   symmetric positive-definite matrix which is given by
                        its upper or lower triangular part.
                        Array whose indexes range within [0..N-1, 0..N-1].
        IsUpperB    -   storage format of matrix B.
        ZNeeded     -   if ZNeeded is equal to:
                         * 0, the eigenvectors are not returned;
                         * 1, the eigenvectors are returned.
        ProblemType -   if ProblemType is equal to:
                         * 1, the following problem is solved: A*x = lambda*B*x;
                         * 2, the following problem is solved: A*B*x = lambda*x;
                         * 3, the following problem is solved: B*A*x = lambda*x.

    Output parameters:
        D           -   eigenvalues in ascending order.
                        Array whose index ranges within [0..N-1].
        Z           -   if ZNeeded is equal to:
                         * 0, Z hasn’t changed;
                         * 1, Z contains eigenvectors.
                        Array whose indexes range within [0..N-1, 0..N-1].
                        The eigenvectors are stored in matrix columns. It should
                        be noted that the eigenvectors in such problems do not
                        form an orthogonal system.

    Result:
        True, if the problem was solved successfully.
        False, if the error occurred during the Cholesky decomposition of matrix
        B (the matrix isn’t positive-definite) or during the work of the iterative
        algorithm for solving the symmetric eigenproblem.

    See also the GeneralizedSymmetricDefiniteEVDReduce subroutine.

      -- ALGLIB --
         Copyright 1.28.2006 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixgevd(ap::template_2d_array< amp::ampf<Precision> > a,
        int n,
        bool isuppera,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        bool isupperb,
        int zneeded,
        int problemtype,
        ap::template_1d_array< amp::ampf<Precision> >& d,
        ap::template_2d_array< amp::ampf<Precision> >& z)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > r;
        ap::template_2d_array< amp::ampf<Precision> > t;
        bool isupperr;
        int j1;
        int j2;
        int j1inc;
        int j2inc;
        int i;
        int j;
        amp::ampf<Precision> v;


        
        //
        // Reduce and solve
        //
        result = smatrixgevdreduce<Precision>(a, n, isuppera, b, isupperb, problemtype, r, isupperr);
        if( !result )
        {
            return result;
        }
        result = evd::smatrixevd<Precision>(a, n, zneeded, isuppera, d, t);
        if( !result )
        {
            return result;
        }
        
        //
        // Transform eigenvectors if needed
        //
        if( zneeded!=0 )
        {
            
            //
            // fill Z with zeros
            //
            z.setbounds(0, n-1, 0, n-1);
            for(j=0; j<=n-1; j++)
            {
                z(0,j) = amp::ampf<Precision>("0.0");
            }
            for(i=1; i<=n-1; i++)
            {
                amp::vmove(z.getrow(i, 0, n-1), z.getrow(0, 0, n-1));
            }
            
            //
            // Setup R properties
            //
            if( isupperr )
            {
                j1 = 0;
                j2 = n-1;
                j1inc = +1;
                j2inc = 0;
            }
            else
            {
                j1 = 0;
                j2 = 0;
                j1inc = 0;
                j2inc = +1;
            }
            
            //
            // Calculate R*Z
            //
            for(i=0; i<=n-1; i++)
            {
                for(j=j1; j<=j2; j++)
                {
                    v = r(i,j);
                    amp::vadd(z.getrow(i, 0, n-1), t.getrow(j, 0, n-1), v);
                }
                j1 = j1+j1inc;
                j2 = j2+j2inc;
            }
        }
        return result;
    }


    /*************************************************************************
    Algorithm for reduction of the following generalized symmetric positive-
    definite eigenvalue problem:
        A*x = lambda*B*x (1) or
        A*B*x = lambda*x (2) or
        B*A*x = lambda*x (3)
    to the symmetric eigenvalues problem C*y = lambda*y (eigenvalues of this and
    the given problems are the same, and the eigenvectors of the given problem
    could be obtained by multiplying the obtained eigenvectors by the
    transformation matrix x = R*y).

    Here A is a symmetric matrix, B - symmetric positive-definite matrix.

    Input parameters:
        A           -   symmetric matrix which is given by its upper or lower
                        triangular part.
                        Array whose indexes range within [0..N-1, 0..N-1].
        N           -   size of matrices A and B.
        IsUpperA    -   storage format of matrix A.
        B           -   symmetric positive-definite matrix which is given by
                        its upper or lower triangular part.
                        Array whose indexes range within [0..N-1, 0..N-1].
        IsUpperB    -   storage format of matrix B.
        ProblemType -   if ProblemType is equal to:
                         * 1, the following problem is solved: A*x = lambda*B*x;
                         * 2, the following problem is solved: A*B*x = lambda*x;
                         * 3, the following problem is solved: B*A*x = lambda*x.

    Output parameters:
        A           -   symmetric matrix which is given by its upper or lower
                        triangle depending on IsUpperA. Contains matrix C.
                        Array whose indexes range within [0..N-1, 0..N-1].
        R           -   upper triangular or low triangular transformation matrix
                        which is used to obtain the eigenvectors of a given problem
                        as the product of eigenvectors of C (from the right) and
                        matrix R (from the left). If the matrix is upper
                        triangular, the elements below the main diagonal
                        are equal to 0 (and vice versa). Thus, we can perform
                        the multiplication without taking into account the
                        internal structure (which is an easier though less
                        effective way).
                        Array whose indexes range within [0..N-1, 0..N-1].
        IsUpperR    -   type of matrix R (upper or lower triangular).

    Result:
        True, if the problem was reduced successfully.
        False, if the error occurred during the Cholesky decomposition of
            matrix B (the matrix is not positive-definite).

      -- ALGLIB --
         Copyright 1.28.2006 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixgevdreduce(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isuppera,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        bool isupperb,
        int problemtype,
        ap::template_2d_array< amp::ampf<Precision> >& r,
        bool& isupperr)
    {
        bool result;
        ap::template_2d_array< amp::ampf<Precision> > t;
        ap::template_1d_array< amp::ampf<Precision> > w1;
        ap::template_1d_array< amp::ampf<Precision> > w2;
        ap::template_1d_array< amp::ampf<Precision> > w3;
        int i;
        int j;
        amp::ampf<Precision> v;
        matinv::matinvreport<Precision> rep;
        int info;


        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(problemtype==1 || problemtype==2 || problemtype==3);
        result = true;
        
        //
        // Problem 1:  A*x = lambda*B*x
        //
        // Reducing to:
        //     C*y = lambda*y
        //     C = L^(-1) * A * L^(-T)
        //     x = L^(-T) * y
        //
        if( problemtype==1 )
        {
            
            //
            // Factorize B in T: B = LL'
            //
            t.setbounds(0, n-1, 0, n-1);
            if( isupperb )
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(t.getcolumn(i, i, n-1), b.getrow(i, i, n-1));
                }
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(t.getrow(i, 0, i), b.getrow(i, 0, i));
                }
            }
            if( !trfac::spdmatrixcholesky<Precision>(t, n, false) )
            {
                result = false;
                return result;
            }
            
            //
            // Invert L in T
            //
            matinv::rmatrixtrinverse<Precision>(t, n, false, false, info, rep);
            if( info<=0 )
            {
                result = false;
                return result;
            }
            
            //
            // Build L^(-1) * A * L^(-T) in R
            //
            w1.setbounds(1, n);
            w2.setbounds(1, n);
            r.setbounds(0, n-1, 0, n-1);
            for(j=1; j<=n; j++)
            {
                
                //
                // Form w2 = A * l'(j) (here l'(j) is j-th column of L^(-T))
                //
                amp::vmove(w1.getvector(1, j), t.getrow(j-1, 0, j-1));
                sblas::symmetricmatrixvectormultiply<Precision>(a, isuppera, 0, j-1, w1, amp::ampf<Precision>("1.0"), w2);
                if( isuppera )
                {
                    blas::matrixvectormultiply<Precision>(a, 0, j-1, j, n-1, true, w1, 1, j, amp::ampf<Precision>("1.0"), w2, j+1, n, amp::ampf<Precision>("0.0"));
                }
                else
                {
                    blas::matrixvectormultiply<Precision>(a, j, n-1, 0, j-1, false, w1, 1, j, amp::ampf<Precision>("1.0"), w2, j+1, n, amp::ampf<Precision>("0.0"));
                }
                
                //
                // Form l(i)*w2 (here l(i) is i-th row of L^(-1))
                //
                for(i=1; i<=n; i++)
                {
                    v = amp::vdotproduct(t.getrow(i-1, 0, i-1), w2.getvector(1, i));
                    r(i-1,j-1) = v;
                }
            }
            
            //
            // Copy R to A
            //
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(a.getrow(i, 0, n-1), r.getrow(i, 0, n-1));
            }
            
            //
            // Copy L^(-1) from T to R and transpose
            //
            isupperr = true;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=i-1; j++)
                {
                    r(i,j) = 0;
                }
            }
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(r.getrow(i, i, n-1), t.getcolumn(i, i, n-1));
            }
            return result;
        }
        
        //
        // Problem 2:  A*B*x = lambda*x
        // or
        // problem 3:  B*A*x = lambda*x
        //
        // Reducing to:
        //     C*y = lambda*y
        //     C = U * A * U'
        //     B = U'* U
        //
        if( problemtype==2 || problemtype==3 )
        {
            
            //
            // Factorize B in T: B = U'*U
            //
            t.setbounds(0, n-1, 0, n-1);
            if( isupperb )
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(t.getrow(i, i, n-1), b.getrow(i, i, n-1));
                }
            }
            else
            {
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(t.getrow(i, i, n-1), b.getcolumn(i, i, n-1));
                }
            }
            if( !trfac::spdmatrixcholesky<Precision>(t, n, true) )
            {
                result = false;
                return result;
            }
            
            //
            // Build U * A * U' in R
            //
            w1.setbounds(1, n);
            w2.setbounds(1, n);
            w3.setbounds(1, n);
            r.setbounds(0, n-1, 0, n-1);
            for(j=1; j<=n; j++)
            {
                
                //
                // Form w2 = A * u'(j) (here u'(j) is j-th column of U')
                //
                amp::vmove(w1.getvector(1, n-j+1), t.getrow(j-1, j-1, n-1));
                sblas::symmetricmatrixvectormultiply<Precision>(a, isuppera, j-1, n-1, w1, amp::ampf<Precision>("1.0"), w3);
                amp::vmove(w2.getvector(j, n), w3.getvector(1, n-j+1));
                amp::vmove(w1.getvector(j, n), t.getrow(j-1, j-1, n-1));
                if( isuppera )
                {
                    blas::matrixvectormultiply<Precision>(a, 0, j-2, j-1, n-1, false, w1, j, n, amp::ampf<Precision>("1.0"), w2, 1, j-1, amp::ampf<Precision>("0.0"));
                }
                else
                {
                    blas::matrixvectormultiply<Precision>(a, j-1, n-1, 0, j-2, true, w1, j, n, amp::ampf<Precision>("1.0"), w2, 1, j-1, amp::ampf<Precision>("0.0"));
                }
                
                //
                // Form u(i)*w2 (here u(i) is i-th row of U)
                //
                for(i=1; i<=n; i++)
                {
                    v = amp::vdotproduct(t.getrow(i-1, i-1, n-1), w2.getvector(i, n));
                    r(i-1,j-1) = v;
                }
            }
            
            //
            // Copy R to A
            //
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(a.getrow(i, 0, n-1), r.getrow(i, 0, n-1));
            }
            if( problemtype==2 )
            {
                
                //
                // Invert U in T
                //
                matinv::rmatrixtrinverse<Precision>(t, n, true, false, info, rep);
                if( info<=0 )
                {
                    result = false;
                    return result;
                }
                
                //
                // Copy U^-1 from T to R
                //
                isupperr = true;
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=i-1; j++)
                    {
                        r(i,j) = 0;
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(r.getrow(i, i, n-1), t.getrow(i, i, n-1));
                }
            }
            else
            {
                
                //
                // Copy U from T to R and transpose
                //
                isupperr = false;
                for(i=0; i<=n-1; i++)
                {
                    for(j=i+1; j<=n-1; j++)
                    {
                        r(i,j) = 0;
                    }
                }
                for(i=0; i<=n-1; i++)
                {
                    amp::vmove(r.getcolumn(i, i, n-1), t.getrow(i, i, n-1));
                }
            }
        }
        return result;
    }
} // namespace

#endif
