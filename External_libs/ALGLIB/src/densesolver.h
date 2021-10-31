/*************************************************************************
Copyright (c) 2007-2008, Sergey Bochkanov (ALGLIB project).

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

#ifndef _densesolver_h
#define _densesolver_h

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
#include "bdsvd.h"
#include "svd.h"
#include "hqrnd.h"
#include "matgen.h"
#include "trfac.h"
#include "trlinsolve.h"
#include "safesolve.h"
#include "rcond.h"
#include "xblas.h"
namespace densesolver
{
    template<unsigned int Precision>
    class densesolverreport
    {
    public:
        amp::ampf<Precision> r1;
        amp::ampf<Precision> rinf;
    };


    template<unsigned int Precision>
    class densesolverlsreport
    {
    public:
        amp::ampf<Precision> r2;
        ap::template_2d_array< amp::ampf<Precision> > cx;
        int n;
        int k;
    };




    template<unsigned int Precision>
    void rmatrixsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        bool rfs,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixlusolve(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixlusolvem(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixmixedsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixmixedsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        bool rfs,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixlusolvem(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixlusolve(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixmixedsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixmixedsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void spdmatrixsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void spdmatrixsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void spdmatrixcholeskysolvem(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void spdmatrixcholeskysolve(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void hpdmatrixsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void hpdmatrixsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void hpdmatrixcholeskysolvem(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void hpdmatrixcholeskysolve(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixsolvels(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int nrows,
        int ncols,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        amp::ampf<Precision> threshold,
        int& info,
        densesolverlsreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void rmatrixlusolveinternal(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        const amp::ampf<Precision>& scalea,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void spdmatrixcholeskysolveinternal(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        const amp::ampf<Precision>& sqrtscalea,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    void cmatrixlusolveinternal(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        const amp::ampf<Precision>& scalea,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    void hpdmatrixcholeskysolveinternal(const ap::template_2d_array< amp::campf<Precision> >& cha,
        const amp::ampf<Precision>& sqrtscalea,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x);
    template<unsigned int Precision>
    int densesolverrfsmax(int n,
        amp::ampf<Precision> r1,
        amp::ampf<Precision> rinf);
    template<unsigned int Precision>
    int densesolverrfsmaxv2(int n,
        amp::ampf<Precision> r2);
    template<unsigned int Precision>
    void rbasiclusolve(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        amp::ampf<Precision> scalea,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& xb,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void spdbasiccholeskysolve(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        amp::ampf<Precision> sqrtscalea,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& xb,
        ap::template_1d_array< amp::ampf<Precision> >& tmp);
    template<unsigned int Precision>
    void cbasiclusolve(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        amp::ampf<Precision> scalea,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& xb,
        ap::template_1d_array< amp::campf<Precision> >& tmp);
    template<unsigned int Precision>
    void hpdbasiccholeskysolve(const ap::template_2d_array< amp::campf<Precision> >& cha,
        amp::ampf<Precision> sqrtscalea,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& xb,
        ap::template_1d_array< amp::campf<Precision> >& tmp);


    /*************************************************************************
    Dense solver.

    This  subroutine  solves  a  system  A*x=b,  where A is NxN non-denegerate
    real matrix, x and b are vectors.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(N^3) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   return code:
                    * -3    A is singular, or VERY close to singular.
                            X is filled by zeros in such cases.
                    * -1    N<=0 was passed
                    *  1    task is solved (but matrix A may be ill-conditioned,
                            check R1/RInf parameters for condition numbers).
        Rep     -   solver report, see below for more info
        X       -   array[0..N-1], it contains:
                    * solution of A*x=b if A is non-singular (well-conditioned
                      or ill-conditioned, but not very close to singular)
                    * zeros,  if  A  is  singular  or  VERY  close to singular
                      (in this case Info=-3).

    SOLVER REPORT

    Subroutine sets following fields of the Rep structure:
    * R1        reciprocal of condition number: 1/cond(A), 1-norm.
    * RInf      reciprocal of condition number: 1/cond(A), inf-norm.

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > bm;
        ap::template_2d_array< amp::ampf<Precision> > xm;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        amp::vmove(bm.getcolumn(0, 0, n-1), b.getvector(0, n-1));
        rmatrixsolvem<Precision>(a, n, bm, 1, true, info, rep, xm);
        x.setlength(n);
        amp::vmove(x.getvector(0, n-1), xm.getcolumn(0, 0, n-1));
    }


    /*************************************************************************
    Dense solver.

    Similar to RMatrixSolve() but solves task with multiple right parts (where
    b and x are NxM matrices).

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * optional iterative refinement
    * O(N^3+M*N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size
        RFS     -   iterative refinement switch:
                    * True - refinement is used.
                      Less performance, more precision.
                    * False - refinement is not used.
                      More performance, less precision.

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        bool rfs,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > da;
        ap::template_2d_array< amp::ampf<Precision> > emptya;
        ap::template_1d_array< int > p;
        amp::ampf<Precision> scalea;
        int i;
        int j;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        da.setlength(n, n);
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abs<Precision>(a(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        for(i=0; i<=n-1; i++)
        {
            amp::vmove(da.getrow(i, 0, n-1), a.getrow(i, 0, n-1));
        }
        trfac::rmatrixlu<Precision>(da, n, n, p);
        if( rfs )
        {
            rmatrixlusolveinternal<Precision>(da, p, scalea, n, a, true, b, m, info, rep, x);
        }
        else
        {
            rmatrixlusolveinternal<Precision>(da, p, scalea, n, emptya, false, b, m, info, rep, x);
        }
    }


    /*************************************************************************
    Dense solver.

    This  subroutine  solves  a  system  A*X=B,  where A is NxN non-denegerate
    real matrix given by its LU decomposition, X and B are NxM real matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(N^2) complexity
    * condition number estimation

    No iterative refinement  is provided because exact form of original matrix
    is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
        P       -   array[0..N-1], pivots array, RMatrixLU result
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve
        
      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlusolve(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > bm;
        ap::template_2d_array< amp::ampf<Precision> > xm;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        amp::vmove(bm.getcolumn(0, 0, n-1), b.getvector(0, n-1));
        rmatrixlusolvem<Precision>(lua, p, n, bm, 1, info, rep, xm);
        x.setlength(n);
        amp::vmove(x.getvector(0, n-1), xm.getcolumn(0, 0, n-1));
    }


    /*************************************************************************
    Dense solver.

    Similar to RMatrixLUSolve() but solves task with multiple right parts
    (where b and x are NxM matrices).

    Algorithm features:
    * automatic detection of degenerate cases
    * O(M*N^2) complexity
    * condition number estimation

    No iterative refinement  is provided because exact form of original matrix
    is not known to subroutine. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
        P       -   array[0..N-1], pivots array, RMatrixLU result
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlusolvem(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > emptya;
        int i;
        int j;
        amp::ampf<Precision> scalea;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|U[i,j]|)
        //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
        // 2. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=i; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abs<Precision>(lua(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        rmatrixlusolveinternal<Precision>(lua, p, scalea, n, emptya, false, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver.

    This  subroutine  solves  a  system  A*x=b,  where BOTH ORIGINAL A AND ITS
    LU DECOMPOSITION ARE KNOWN. You can use it if for some  reasons  you  have
    both A and its LU decomposition.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
        P       -   array[0..N-1], pivots array, RMatrixLU result
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolveM
        Rep     -   same as in RMatrixSolveM
        X       -   same as in RMatrixSolveM

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixmixedsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > bm;
        ap::template_2d_array< amp::ampf<Precision> > xm;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        amp::vmove(bm.getcolumn(0, 0, n-1), b.getvector(0, n-1));
        rmatrixmixedsolvem<Precision>(a, lua, p, n, bm, 1, info, rep, xm);
        x.setlength(n);
        amp::vmove(x.getvector(0, n-1), xm.getcolumn(0, 0, n-1));
    }


    /*************************************************************************
    Dense solver.

    Similar to RMatrixMixedSolve() but  solves task with multiple right  parts
    (where b and x are NxM matrices).

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(M*N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
        P       -   array[0..N-1], pivots array, RMatrixLU result
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolveM
        Rep     -   same as in RMatrixSolveM
        X       -   same as in RMatrixSolveM

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixmixedsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        amp::ampf<Precision> scalea;
        int i;
        int j;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abs<Precision>(a(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        rmatrixlusolveinternal<Precision>(lua, p, scalea, n, a, true, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolveM(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(N^3+M*N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size
        RFS     -   iterative refinement switch:
                    * True - refinement is used.
                      Less performance, more precision.
                    * False - refinement is not used.
                      More performance, less precision.

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        bool rfs,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > da;
        ap::template_2d_array< amp::campf<Precision> > emptya;
        ap::template_1d_array< int > p;
        amp::ampf<Precision> scalea;
        int i;
        int j;
        int i_;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        da.setlength(n, n);
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        for(i=0; i<=n-1; i++)
        {
            for(i_=0; i_<=n-1;i_++)
            {
                da(i,i_) = a(i,i_);
            }
        }
        trfac::cmatrixlu<Precision>(da, n, n, p);
        if( rfs )
        {
            cmatrixlusolveinternal<Precision>(da, p, scalea, n, a, true, b, m, info, rep, x);
        }
        else
        {
            cmatrixlusolveinternal<Precision>(da, p, scalea, n, emptya, false, b, m, info, rep, x);
        }
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolve(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(N^3) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > bm;
        ap::template_2d_array< amp::campf<Precision> > xm;
        int i_;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            bm(i_,0) = b(i_);
        }
        cmatrixsolvem<Precision>(a, n, bm, 1, true, info, rep, xm);
        x.setlength(n);
        for(i_=0; i_<=n-1;i_++)
        {
            x(i_) = xm(i_,0);
        }
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolveM(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(M*N^2) complexity
    * condition number estimation

    No iterative refinement  is provided because exact form of original matrix
    is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        LUA     -   array[0..N-1,0..N-1], LU decomposition, RMatrixLU result
        P       -   array[0..N-1], pivots array, RMatrixLU result
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlusolvem(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > emptya;
        int i;
        int j;
        amp::ampf<Precision> scalea;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|U[i,j]|)
        //    we assume that LU is in its normal form, i.e. |L[i,j]|<=1
        // 2. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=i; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abscomplex<Precision>(lua(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        cmatrixlusolveinternal<Precision>(lua, p, scalea, n, emptya, false, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolve(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(N^2) complexity
    * condition number estimation

    No iterative refinement is provided because exact form of original matrix
    is not known to subroutine. Use CMatrixSolve or CMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
        P       -   array[0..N-1], pivots array, CMatrixLU result
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlusolve(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > bm;
        ap::template_2d_array< amp::campf<Precision> > xm;
        int i_;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            bm(i_,0) = b(i_);
        }
        cmatrixlusolvem<Precision>(lua, p, n, bm, 1, info, rep, xm);
        x.setlength(n);
        for(i_=0; i_<=n-1;i_++)
        {
            x(i_) = xm(i_,0);
        }
    }


    /*************************************************************************
    Dense solver. Same as RMatrixMixedSolveM(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(M*N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
        P       -   array[0..N-1], pivots array, CMatrixLU result
        N       -   size of A
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolveM
        Rep     -   same as in RMatrixSolveM
        X       -   same as in RMatrixSolveM

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixmixedsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        amp::ampf<Precision> scalea;
        int i;
        int j;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        scalea = 0;
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                scalea = amp::maximum<Precision>(scalea, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( scalea==0 )
        {
            scalea = 1;
        }
        scalea = 1/scalea;
        cmatrixlusolveinternal<Precision>(lua, p, scalea, n, a, true, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixMixedSolve(), but for complex matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * iterative refinement
    * O(N^2) complexity

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        LUA     -   array[0..N-1,0..N-1], LU decomposition, CMatrixLU result
        P       -   array[0..N-1], pivots array, CMatrixLU result
        N       -   size of A
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolveM
        Rep     -   same as in RMatrixSolveM
        X       -   same as in RMatrixSolveM

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixmixedsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        int n,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > bm;
        ap::template_2d_array< amp::campf<Precision> > xm;
        int i_;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            bm(i_,0) = b(i_);
        }
        cmatrixmixedsolvem<Precision>(a, lua, p, n, bm, 1, info, rep, xm);
        x.setlength(n);
        for(i_=0; i_<=n-1;i_++)
        {
            x(i_) = xm(i_,0);
        }
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolveM(), but for symmetric positive definite
    matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * O(N^3+M*N^2) complexity
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        IsUpper -   what half of A is provided
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve.
                    Returns -3 for non-SPD matrices.
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixsolvem(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > da;
        amp::ampf<Precision> sqrtscalea;
        int i;
        int j;
        int j1;
        int j2;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        da.setlength(n, n);
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        sqrtscalea = 0;
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
                sqrtscalea = amp::maximum<Precision>(sqrtscalea, amp::abs<Precision>(a(i,j)));
            }
        }
        if( sqrtscalea==0 )
        {
            sqrtscalea = 1;
        }
        sqrtscalea = 1/sqrtscalea;
        sqrtscalea = amp::sqrt<Precision>(sqrtscalea);
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
            amp::vmove(da.getrow(i, j1, j2), a.getrow(i, j1, j2));
        }
        if( !trfac::spdmatrixcholesky<Precision>(da, n, isupper) )
        {
            x.setlength(n, m);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        spdmatrixcholeskysolveinternal<Precision>(da, sqrtscalea, n, isupper, a, true, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolve(), but for SPD matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * O(N^3) complexity
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        IsUpper -   what half of A is provided
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
                    Returns -3 for non-SPD matrices.
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixsolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > bm;
        ap::template_2d_array< amp::ampf<Precision> > xm;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        amp::vmove(bm.getcolumn(0, 0, n-1), b.getvector(0, n-1));
        spdmatrixsolvem<Precision>(a, n, isupper, bm, 1, info, rep, xm);
        x.setlength(n);
        amp::vmove(x.getvector(0, n-1), xm.getcolumn(0, 0, n-1));
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolveM(), but for SPD matrices  represented
    by their Cholesky decomposition.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(M*N^2) complexity
    * condition number estimation
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                    SPDMatrixCholesky result
        N       -   size of CHA
        IsUpper -   what half of CHA is provided
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixcholeskysolvem(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > emptya;
        amp::ampf<Precision> sqrtscalea;
        int i;
        int j;
        int j1;
        int j2;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|U[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        sqrtscalea = 0;
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
                sqrtscalea = amp::maximum<Precision>(sqrtscalea, amp::abs<Precision>(cha(i,j)));
            }
        }
        if( sqrtscalea==0 )
        {
            sqrtscalea = 1;
        }
        sqrtscalea = 1/sqrtscalea;
        spdmatrixcholeskysolveinternal<Precision>(cha, sqrtscalea, n, isupper, emptya, false, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolve(), but for  SPD matrices  represented
    by their Cholesky decomposition.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(N^2) complexity
    * condition number estimation
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                    SPDMatrixCholesky result
        N       -   size of A
        IsUpper -   what half of CHA is provided
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixcholeskysolve(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_2d_array< amp::ampf<Precision> > bm;
        ap::template_2d_array< amp::ampf<Precision> > xm;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        amp::vmove(bm.getcolumn(0, 0, n-1), b.getvector(0, n-1));
        spdmatrixcholeskysolvem<Precision>(cha, n, isupper, bm, 1, info, rep, xm);
        x.setlength(n);
        amp::vmove(x.getvector(0, n-1), xm.getcolumn(0, 0, n-1));
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolveM(), but for Hermitian positive definite
    matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * O(N^3+M*N^2) complexity
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        IsUpper -   what half of A is provided
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve.
                    Returns -3 for non-HPD matrices.
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixsolvem(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > da;
        amp::ampf<Precision> sqrtscalea;
        int i;
        int j;
        int j1;
        int j2;
        int i_;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        da.setlength(n, n);
        
        //
        // 1. scale matrix, max(|A[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        sqrtscalea = 0;
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
                sqrtscalea = amp::maximum<Precision>(sqrtscalea, amp::abscomplex<Precision>(a(i,j)));
            }
        }
        if( sqrtscalea==0 )
        {
            sqrtscalea = 1;
        }
        sqrtscalea = 1/sqrtscalea;
        sqrtscalea = amp::sqrt<Precision>(sqrtscalea);
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
            for(i_=j1; i_<=j2;i_++)
            {
                da(i,i_) = a(i,i_);
            }
        }
        if( !trfac::hpdmatrixcholesky<Precision>(da, n, isupper) )
        {
            x.setlength(n, m);
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        hpdmatrixcholeskysolveinternal<Precision>(da, sqrtscalea, n, isupper, a, true, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixSolve(),  but for Hermitian positive definite
    matrices.

    Algorithm features:
    * automatic detection of degenerate cases
    * condition number estimation
    * O(N^3) complexity
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        A       -   array[0..N-1,0..N-1], system matrix
        N       -   size of A
        IsUpper -   what half of A is provided
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
                    Returns -3 for non-HPD matrices.
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixsolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > bm;
        ap::template_2d_array< amp::campf<Precision> > xm;
        int i_;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            bm(i_,0) = b(i_);
        }
        hpdmatrixsolvem<Precision>(a, n, isupper, bm, 1, info, rep, xm);
        x.setlength(n);
        for(i_=0; i_<=n-1;i_++)
        {
            x(i_) = xm(i_,0);
        }
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolveM(), but for HPD matrices  represented
    by their Cholesky decomposition.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(M*N^2) complexity
    * condition number estimation
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                    HPDMatrixCholesky result
        N       -   size of CHA
        IsUpper -   what half of CHA is provided
        B       -   array[0..N-1,0..M-1], right part
        M       -   right part size

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixcholeskysolvem(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > emptya;
        amp::ampf<Precision> sqrtscalea;
        int i;
        int j;
        int j1;
        int j2;


        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        
        //
        // 1. scale matrix, max(|U[i,j]|)
        // 2. factorize scaled matrix
        // 3. solve
        //
        sqrtscalea = 0;
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
                sqrtscalea = amp::maximum<Precision>(sqrtscalea, amp::abscomplex<Precision>(cha(i,j)));
            }
        }
        if( sqrtscalea==0 )
        {
            sqrtscalea = 1;
        }
        sqrtscalea = 1/sqrtscalea;
        hpdmatrixcholeskysolveinternal<Precision>(cha, sqrtscalea, n, isupper, emptya, false, b, m, info, rep, x);
    }


    /*************************************************************************
    Dense solver. Same as RMatrixLUSolve(), but for  HPD matrices  represented
    by their Cholesky decomposition.

    Algorithm features:
    * automatic detection of degenerate cases
    * O(N^2) complexity
    * condition number estimation
    * matrix is represented by its upper or lower triangle

    No iterative refinement is provided because such partial representation of
    matrix does not allow efficient calculation of extra-precise  matrix-vector
    products for large matrices. Use RMatrixSolve or RMatrixMixedSolve  if  you
    need iterative refinement.

    INPUT PARAMETERS
        CHA     -   array[0..N-1,0..N-1], Cholesky decomposition,
                    SPDMatrixCholesky result
        N       -   size of A
        IsUpper -   what half of CHA is provided
        B       -   array[0..N-1], right part

    OUTPUT PARAMETERS
        Info    -   same as in RMatrixSolve
        Rep     -   same as in RMatrixSolve
        X       -   same as in RMatrixSolve

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixcholeskysolve(const ap::template_2d_array< amp::campf<Precision> >& cha,
        int n,
        bool isupper,
        const ap::template_1d_array< amp::campf<Precision> >& b,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_1d_array< amp::campf<Precision> >& x)
    {
        ap::template_2d_array< amp::campf<Precision> > bm;
        ap::template_2d_array< amp::campf<Precision> > xm;
        int i_;


        if( n<=0 )
        {
            info = -1;
            return;
        }
        bm.setlength(n, 1);
        for(i_=0; i_<=n-1;i_++)
        {
            bm(i_,0) = b(i_);
        }
        hpdmatrixcholeskysolvem<Precision>(cha, n, isupper, bm, 1, info, rep, xm);
        x.setlength(n);
        for(i_=0; i_<=n-1;i_++)
        {
            x(i_) = xm(i_,0);
        }
    }


    /*************************************************************************
    Dense solver.

    This subroutine finds solution of the linear system A*X=B with non-square,
    possibly degenerate A.  System  is  solved in the least squares sense, and
    general least squares solution  X = X0 + CX*y  which  minimizes |A*X-B| is
    returned. If A is non-degenerate, solution in the  usual sense is returned

    Algorithm features:
    * automatic detection of degenerate cases
    * iterative refinement
    * O(N^3) complexity

    INPUT PARAMETERS
        A       -   array[0..NRows-1,0..NCols-1], system matrix
        NRows   -   vertical size of A
        NCols   -   horizontal size of A
        B       -   array[0..NCols-1], right part
        Threshold-  a number in [0,1]. Singular values  beyond  Threshold  are
                    considered  zero.  Set  it to 0.0, if you don't understand
                    what it means, so the solver will choose good value on its
                    own.
                    
    OUTPUT PARAMETERS
        Info    -   return code:
                    * -4    SVD subroutine failed
                    * -1    if NRows<=0 or NCols<=0 or Threshold<0 was passed
                    *  1    if task is solved
        Rep     -   solver report, see below for more info
        X       -   array[0..N-1,0..M-1], it contains:
                    * solution of A*X=B if A is non-singular (well-conditioned
                      or ill-conditioned, but not very close to singular)
                    * zeros,  if  A  is  singular  or  VERY  close to singular
                      (in this case Info=-3).

    SOLVER REPORT

    Subroutine sets following fields of the Rep structure:
    * R2        reciprocal of condition number: 1/cond(A), 2-norm.
    * N         = NCols
    * K         dim(Null(A))
    * CX        array[0..N-1,0..K-1], kernel of A.
                Columns of CX store such vectors that A*CX[i]=0.

      -- ALGLIB --
         Copyright 24.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixsolvels(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int nrows,
        int ncols,
        const ap::template_1d_array< amp::ampf<Precision> >& b,
        amp::ampf<Precision> threshold,
        int& info,
        densesolverlsreport<Precision>& rep,
        ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        ap::template_1d_array< amp::ampf<Precision> > sv;
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > vt;
        ap::template_1d_array< amp::ampf<Precision> > rp;
        ap::template_1d_array< amp::ampf<Precision> > utb;
        ap::template_1d_array< amp::ampf<Precision> > sutb;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        ap::template_1d_array< amp::ampf<Precision> > ta;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        ap::template_1d_array< amp::ampf<Precision> > buf;
        ap::template_1d_array< amp::ampf<Precision> > w;
        int i;
        int j;
        int nsv;
        int kernelidx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        bool svdfailed;
        bool zeroa;
        int rfs;
        int nrfs;
        bool terminatenexttime;
        bool smallerr;


        if( nrows<=0 || ncols<=0 || threshold<0 )
        {
            info = -1;
            return;
        }
        if( threshold==0 )
        {
            threshold = 1000*amp::ampf<Precision>::getAlgoPascalEpsilon();
        }
        
        //
        // Factorize A first
        //
        svdfailed = !svd::rmatrixsvd<Precision>(a, nrows, ncols, 1, 2, 2, sv, u, vt);
        zeroa = sv(0)==0;
        if( svdfailed || zeroa )
        {
            if( svdfailed )
            {
                info = -4;
            }
            else
            {
                info = 1;
            }
            x.setlength(ncols);
            for(i=0; i<=ncols-1; i++)
            {
                x(i) = 0;
            }
            rep.n = ncols;
            rep.k = ncols;
            rep.cx.setlength(ncols, ncols);
            for(i=0; i<=ncols-1; i++)
            {
                for(j=0; j<=ncols-1; j++)
                {
                    if( i==j )
                    {
                        rep.cx(i,j) = 1;
                    }
                    else
                    {
                        rep.cx(i,j) = 0;
                    }
                }
            }
            rep.r2 = 0;
            return;
        }
        nsv = ap::minint(ncols, nrows);
        if( nsv==ncols )
        {
            rep.r2 = sv(nsv-1)/sv(0);
        }
        else
        {
            rep.r2 = 0;
        }
        rep.n = ncols;
        info = 1;
        
        //
        // Iterative refinement of xc combined with solution:
        // 1. xc = 0
        // 2. calculate r = bc-A*xc using extra-precise dot product
        // 3. solve A*y = r
        // 4. update x:=x+r
        // 5. goto 2
        //
        // This cycle is executed until one of two things happens:
        // 1. maximum number of iterations reached
        // 2. last iteration decreased error to the lower limit
        //
        utb.setlength(nsv);
        sutb.setlength(nsv);
        x.setlength(ncols);
        tmp.setlength(ncols);
        ta.setlength(ncols+1);
        tx.setlength(ncols+1);
        buf.setlength(ncols+1);
        for(i=0; i<=ncols-1; i++)
        {
            x(i) = 0;
        }
        kernelidx = nsv;
        for(i=0; i<=nsv-1; i++)
        {
            if( sv(i)<=threshold*sv(0) )
            {
                kernelidx = i;
                break;
            }
        }
        rep.k = ncols-kernelidx;
        nrfs = densesolverrfsmaxv2<Precision>(ncols, rep.r2);
        terminatenexttime = false;
        rp.setlength(nrows);
        for(rfs=0; rfs<=nrfs; rfs++)
        {
            if( terminatenexttime )
            {
                break;
            }
            
            //
            // calculate right part
            //
            if( rfs==0 )
            {
                amp::vmove(rp.getvector(0, nrows-1), b.getvector(0, nrows-1));
            }
            else
            {
                smallerr = true;
                for(i=0; i<=nrows-1; i++)
                {
                    amp::vmove(ta.getvector(0, ncols-1), a.getrow(i, 0, ncols-1));
                    ta(ncols) = -1;
                    amp::vmove(tx.getvector(0, ncols-1), x.getvector(0, ncols-1));
                    tx(ncols) = b(i);
                    xblas::xdot<Precision>(ta, tx, ncols+1, buf, v, verr);
                    rp(i) = -v;
                    smallerr = smallerr && amp::abs<Precision>(v)<4*verr;
                }
                if( smallerr )
                {
                    terminatenexttime = true;
                }
            }
            
            //
            // solve A*dx = rp
            //
            for(i=0; i<=ncols-1; i++)
            {
                tmp(i) = 0;
            }
            for(i=0; i<=nsv-1; i++)
            {
                utb(i) = 0;
            }
            for(i=0; i<=nrows-1; i++)
            {
                v = rp(i);
                amp::vadd(utb.getvector(0, nsv-1), u.getrow(i, 0, nsv-1), v);
            }
            for(i=0; i<=nsv-1; i++)
            {
                if( i<kernelidx )
                {
                    sutb(i) = utb(i)/sv(i);
                }
                else
                {
                    sutb(i) = 0;
                }
            }
            for(i=0; i<=nsv-1; i++)
            {
                v = sutb(i);
                amp::vadd(tmp.getvector(0, ncols-1), vt.getrow(i, 0, ncols-1), v);
            }
            
            //
            // update x:  x:=x+dx
            //
            amp::vadd(x.getvector(0, ncols-1), tmp.getvector(0, ncols-1));
        }
        
        //
        // fill CX
        //
        if( rep.k>0 )
        {
            rep.cx.setlength(ncols, rep.k);
            for(i=0; i<=rep.k-1; i++)
            {
                amp::vmove(rep.cx.getcolumn(i, 0, ncols-1), vt.getrow(kernelidx+i, 0, ncols-1));
            }
        }
    }


    /*************************************************************************
    Internal LU solver

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixlusolveinternal(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        const amp::ampf<Precision>& scalea,
        int n,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        int i;
        int j;
        int k;
        int rfs;
        int nrfs;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > bc;
        ap::template_1d_array< amp::ampf<Precision> > xa;
        ap::template_1d_array< amp::ampf<Precision> > xb;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        amp::ampf<Precision> mxb;
        amp::ampf<Precision> scaleright;
        bool smallerr;
        bool terminatenexttime;


        ap::ap_error::make_assertion(scalea>0);
        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            if( p(i)>n-1 || p(i)<i )
            {
                info = -1;
                return;
            }
        }
        x.setlength(n, m);
        y.setlength(n);
        xc.setlength(n);
        bc.setlength(n);
        tx.setlength(n+1);
        xa.setlength(n+1);
        xb.setlength(n+1);
        
        //
        // estimate condition number, test for near singularity
        //
        rep.r1 = rcond::rmatrixlurcond1<Precision>(lua, n);
        rep.rinf = rcond::rmatrixlurcondinf<Precision>(lua, n);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        
        //
        // solve
        //
        for(k=0; k<=m-1; k++)
        {
            
            //
            // copy B to contiguous storage
            //
            amp::vmove(bc.getvector(0, n-1), b.getcolumn(k, 0, n-1));
            
            //
            // Scale right part:
            // * MX stores max(|Bi|)
            // * ScaleRight stores actual scaling applied to B when solving systems
            //   it is chosen to make |scaleRight*b| close to 1.
            //
            mxb = 0;
            for(i=0; i<=n-1; i++)
            {
                mxb = amp::maximum<Precision>(mxb, amp::abs<Precision>(bc(i)));
            }
            if( mxb==0 )
            {
                mxb = 1;
            }
            scaleright = 1/mxb;
            
            //
            // First, non-iterative part of solution process.
            // We use separate code for this task because
            // XDot is quite slow and we want to save time.
            //
            amp::vmove(xc.getvector(0, n-1), bc.getvector(0, n-1), scaleright);
            rbasiclusolve<Precision>(lua, p, scalea, n, xc, tx);
            
            //
            // Iterative refinement of xc:
            // * calculate r = bc-A*xc using extra-precise dot product
            // * solve A*y = r
            // * update x:=x+r
            //
            // This cycle is executed until one of two things happens:
            // 1. maximum number of iterations reached
            // 2. last iteration decreased error to the lower limit
            //
            if( havea )
            {
                nrfs = densesolverrfsmax<Precision>(n, rep.r1, rep.rinf);
                terminatenexttime = false;
                for(rfs=0; rfs<=nrfs-1; rfs++)
                {
                    if( terminatenexttime )
                    {
                        break;
                    }
                    
                    //
                    // generate right part
                    //
                    smallerr = true;
                    amp::vmove(xb.getvector(0, n-1), xc.getvector(0, n-1));
                    for(i=0; i<=n-1; i++)
                    {
                        amp::vmove(xa.getvector(0, n-1), a.getrow(i, 0, n-1), scalea);
                        xa(n) = -1;
                        xb(n) = scaleright*bc(i);
                        xblas::xdot<Precision>(xa, xb, n+1, tx, v, verr);
                        y(i) = -v;
                        smallerr = smallerr && amp::abs<Precision>(v)<4*verr;
                    }
                    if( smallerr )
                    {
                        terminatenexttime = true;
                    }
                    
                    //
                    // solve and update
                    //
                    rbasiclusolve<Precision>(lua, p, scalea, n, y, tx);
                    amp::vadd(xc.getvector(0, n-1), y.getvector(0, n-1));
                }
            }
            
            //
            // Store xc.
            // Post-scale result.
            //
            v = scalea*mxb;
            amp::vmove(x.getcolumn(k, 0, n-1), xc.getvector(0, n-1), v);
        }
    }


    /*************************************************************************
    Internal Cholesky solver

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdmatrixcholeskysolveinternal(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        const amp::ampf<Precision>& sqrtscalea,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::ampf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::ampf<Precision> >& x)
    {
        int i;
        int j;
        int k;
        int rfs;
        int nrfs;
        ap::template_1d_array< amp::ampf<Precision> > xc;
        ap::template_1d_array< amp::ampf<Precision> > y;
        ap::template_1d_array< amp::ampf<Precision> > bc;
        ap::template_1d_array< amp::ampf<Precision> > xa;
        ap::template_1d_array< amp::ampf<Precision> > xb;
        ap::template_1d_array< amp::ampf<Precision> > tx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        amp::ampf<Precision> mxb;
        amp::ampf<Precision> scaleright;
        bool smallerr;
        bool terminatenexttime;


        ap::ap_error::make_assertion(sqrtscalea>0);
        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        x.setlength(n, m);
        y.setlength(n);
        xc.setlength(n);
        bc.setlength(n);
        tx.setlength(n+1);
        xa.setlength(n+1);
        xb.setlength(n+1);
        
        //
        // estimate condition number, test for near singularity
        //
        rep.r1 = rcond::spdmatrixcholeskyrcond<Precision>(cha, n, isupper);
        rep.rinf = rep.r1;
        if( rep.r1<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        
        //
        // solve
        //
        for(k=0; k<=m-1; k++)
        {
            
            //
            // copy B to contiguous storage
            //
            amp::vmove(bc.getvector(0, n-1), b.getcolumn(k, 0, n-1));
            
            //
            // Scale right part:
            // * MX stores max(|Bi|)
            // * ScaleRight stores actual scaling applied to B when solving systems
            //   it is chosen to make |scaleRight*b| close to 1.
            //
            mxb = 0;
            for(i=0; i<=n-1; i++)
            {
                mxb = amp::maximum<Precision>(mxb, amp::abs<Precision>(bc(i)));
            }
            if( mxb==0 )
            {
                mxb = 1;
            }
            scaleright = 1/mxb;
            
            //
            // First, non-iterative part of solution process.
            // We use separate code for this task because
            // XDot is quite slow and we want to save time.
            //
            amp::vmove(xc.getvector(0, n-1), bc.getvector(0, n-1), scaleright);
            spdbasiccholeskysolve<Precision>(cha, sqrtscalea, n, isupper, xc, tx);
            
            //
            // Store xc.
            // Post-scale result.
            //
            v = amp::sqr<Precision>(sqrtscalea)*mxb;
            amp::vmove(x.getcolumn(k, 0, n-1), xc.getvector(0, n-1), v);
        }
    }


    /*************************************************************************
    Internal LU solver

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cmatrixlusolveinternal(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        const amp::ampf<Precision>& scalea,
        int n,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        int i;
        int j;
        int k;
        int rfs;
        int nrfs;
        ap::template_1d_array< amp::campf<Precision> > xc;
        ap::template_1d_array< amp::campf<Precision> > y;
        ap::template_1d_array< amp::campf<Precision> > bc;
        ap::template_1d_array< amp::campf<Precision> > xa;
        ap::template_1d_array< amp::campf<Precision> > xb;
        ap::template_1d_array< amp::campf<Precision> > tx;
        ap::template_1d_array< amp::ampf<Precision> > tmpbuf;
        amp::campf<Precision> v;
        amp::ampf<Precision> verr;
        amp::ampf<Precision> mxb;
        amp::ampf<Precision> scaleright;
        bool smallerr;
        bool terminatenexttime;
        int i_;


        ap::ap_error::make_assertion(scalea>0);
        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        for(i=0; i<=n-1; i++)
        {
            if( p(i)>n-1 || p(i)<i )
            {
                info = -1;
                return;
            }
        }
        x.setlength(n, m);
        y.setlength(n);
        xc.setlength(n);
        bc.setlength(n);
        tx.setlength(n);
        xa.setlength(n+1);
        xb.setlength(n+1);
        tmpbuf.setlength(2*n+2);
        
        //
        // estimate condition number, test for near singularity
        //
        rep.r1 = rcond::cmatrixlurcond1<Precision>(lua, n);
        rep.rinf = rcond::cmatrixlurcondinf<Precision>(lua, n);
        if( rep.r1<rcond::rcondthreshold<Precision>() || rep.rinf<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        
        //
        // solve
        //
        for(k=0; k<=m-1; k++)
        {
            
            //
            // copy B to contiguous storage
            //
            for(i_=0; i_<=n-1;i_++)
            {
                bc(i_) = b(i_,k);
            }
            
            //
            // Scale right part:
            // * MX stores max(|Bi|)
            // * ScaleRight stores actual scaling applied to B when solving systems
            //   it is chosen to make |scaleRight*b| close to 1.
            //
            mxb = 0;
            for(i=0; i<=n-1; i++)
            {
                mxb = amp::maximum<Precision>(mxb, amp::abscomplex<Precision>(bc(i)));
            }
            if( mxb==0 )
            {
                mxb = 1;
            }
            scaleright = 1/mxb;
            
            //
            // First, non-iterative part of solution process.
            // We use separate code for this task because
            // XDot is quite slow and we want to save time.
            //
            for(i_=0; i_<=n-1;i_++)
            {
                xc(i_) = scaleright*bc(i_);
            }
            cbasiclusolve<Precision>(lua, p, scalea, n, xc, tx);
            
            //
            // Iterative refinement of xc:
            // * calculate r = bc-A*xc using extra-precise dot product
            // * solve A*y = r
            // * update x:=x+r
            //
            // This cycle is executed until one of two things happens:
            // 1. maximum number of iterations reached
            // 2. last iteration decreased error to the lower limit
            //
            if( havea )
            {
                nrfs = densesolverrfsmax<Precision>(n, rep.r1, rep.rinf);
                terminatenexttime = false;
                for(rfs=0; rfs<=nrfs-1; rfs++)
                {
                    if( terminatenexttime )
                    {
                        break;
                    }
                    
                    //
                    // generate right part
                    //
                    smallerr = true;
                    for(i_=0; i_<=n-1;i_++)
                    {
                        xb(i_) = xc(i_);
                    }
                    for(i=0; i<=n-1; i++)
                    {
                        for(i_=0; i_<=n-1;i_++)
                        {
                            xa(i_) = scalea*a(i,i_);
                        }
                        xa(n) = -1;
                        xb(n) = scaleright*bc(i);
                        xblas::xcdot<Precision>(xa, xb, n+1, tmpbuf, v, verr);
                        y(i) = -v;
                        smallerr = smallerr && amp::abscomplex<Precision>(v)<4*verr;
                    }
                    if( smallerr )
                    {
                        terminatenexttime = true;
                    }
                    
                    //
                    // solve and update
                    //
                    cbasiclusolve<Precision>(lua, p, scalea, n, y, tx);
                    for(i_=0; i_<=n-1;i_++)
                    {
                        xc(i_) = xc(i_) + y(i_);
                    }
                }
            }
            
            //
            // Store xc.
            // Post-scale result.
            //
            v = scalea*mxb;
            for(i_=0; i_<=n-1;i_++)
            {
                x(i_,k) = v*xc(i_);
            }
        }
    }


    /*************************************************************************
    Internal Cholesky solver

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdmatrixcholeskysolveinternal(const ap::template_2d_array< amp::campf<Precision> >& cha,
        const amp::ampf<Precision>& sqrtscalea,
        int n,
        bool isupper,
        const ap::template_2d_array< amp::campf<Precision> >& a,
        bool havea,
        const ap::template_2d_array< amp::campf<Precision> >& b,
        int m,
        int& info,
        densesolverreport<Precision>& rep,
        ap::template_2d_array< amp::campf<Precision> >& x)
    {
        int i;
        int j;
        int k;
        int rfs;
        int nrfs;
        ap::template_1d_array< amp::campf<Precision> > xc;
        ap::template_1d_array< amp::campf<Precision> > y;
        ap::template_1d_array< amp::campf<Precision> > bc;
        ap::template_1d_array< amp::campf<Precision> > xa;
        ap::template_1d_array< amp::campf<Precision> > xb;
        ap::template_1d_array< amp::campf<Precision> > tx;
        amp::ampf<Precision> v;
        amp::ampf<Precision> verr;
        amp::ampf<Precision> mxb;
        amp::ampf<Precision> scaleright;
        bool smallerr;
        bool terminatenexttime;
        int i_;


        ap::ap_error::make_assertion(sqrtscalea>0);
        
        //
        // prepare: check inputs, allocate space...
        //
        if( n<=0 || m<=0 )
        {
            info = -1;
            return;
        }
        x.setlength(n, m);
        y.setlength(n);
        xc.setlength(n);
        bc.setlength(n);
        tx.setlength(n+1);
        xa.setlength(n+1);
        xb.setlength(n+1);
        
        //
        // estimate condition number, test for near singularity
        //
        rep.r1 = rcond::hpdmatrixcholeskyrcond<Precision>(cha, n, isupper);
        rep.rinf = rep.r1;
        if( rep.r1<rcond::rcondthreshold<Precision>() )
        {
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=m-1; j++)
                {
                    x(i,j) = 0;
                }
            }
            rep.r1 = 0;
            rep.rinf = 0;
            info = -3;
            return;
        }
        info = 1;
        
        //
        // solve
        //
        for(k=0; k<=m-1; k++)
        {
            
            //
            // copy B to contiguous storage
            //
            for(i_=0; i_<=n-1;i_++)
            {
                bc(i_) = b(i_,k);
            }
            
            //
            // Scale right part:
            // * MX stores max(|Bi|)
            // * ScaleRight stores actual scaling applied to B when solving systems
            //   it is chosen to make |scaleRight*b| close to 1.
            //
            mxb = 0;
            for(i=0; i<=n-1; i++)
            {
                mxb = amp::maximum<Precision>(mxb, amp::abscomplex<Precision>(bc(i)));
            }
            if( mxb==0 )
            {
                mxb = 1;
            }
            scaleright = 1/mxb;
            
            //
            // First, non-iterative part of solution process.
            // We use separate code for this task because
            // XDot is quite slow and we want to save time.
            //
            for(i_=0; i_<=n-1;i_++)
            {
                xc(i_) = scaleright*bc(i_);
            }
            hpdbasiccholeskysolve<Precision>(cha, sqrtscalea, n, isupper, xc, tx);
            
            //
            // Store xc.
            // Post-scale result.
            //
            v = amp::sqr<Precision>(sqrtscalea)*mxb;
            for(i_=0; i_<=n-1;i_++)
            {
                x(i_,k) = v*xc(i_);
            }
        }
    }


    /*************************************************************************
    Internal subroutine.
    Returns maximum count of RFS iterations as function of:
    1. machine epsilon
    2. task size.
    3. condition number

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int densesolverrfsmax(int n,
        amp::ampf<Precision> r1,
        amp::ampf<Precision> rinf)
    {
        int result;


        result = 5;
        return result;
    }


    /*************************************************************************
    Internal subroutine.
    Returns maximum count of RFS iterations as function of:
    1. machine epsilon
    2. task size.
    3. norm-2 condition number

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int densesolverrfsmaxv2(int n,
        amp::ampf<Precision> r2)
    {
        int result;


        result = densesolverrfsmax<Precision>(n, amp::ampf<Precision>(0), amp::ampf<Precision>(0));
        return result;
    }


    /*************************************************************************
    Basic LU solver for ScaleA*PLU*x = y.

    This subroutine assumes that:
    * L is well-scaled, and it is U which needs scaling by ScaleA.
    * A=PLU is well-conditioned, so no zero divisions or overflow may occur

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rbasiclusolve(const ap::template_2d_array< amp::ampf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        amp::ampf<Precision> scalea,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& xb,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        amp::ampf<Precision> v;


        for(i=0; i<=n-1; i++)
        {
            if( p(i)!=i )
            {
                v = xb(i);
                xb(i) = xb(p(i));
                xb(p(i)) = v;
            }
        }
        for(i=1; i<=n-1; i++)
        {
            v = amp::vdotproduct(lua.getrow(i, 0, i-1), xb.getvector(0, i-1));
            xb(i) = xb(i)-v;
        }
        xb(n-1) = xb(n-1)/(scalea*lua(n-1,n-1));
        for(i=n-2; i>=0; i--)
        {
            amp::vmove(tmp.getvector(i+1, n-1), lua.getrow(i, i+1, n-1), scalea);
            v = amp::vdotproduct(tmp.getvector(i+1, n-1), xb.getvector(i+1, n-1));
            xb(i) = (xb(i)-v)/(scalea*lua(i,i));
        }
    }


    /*************************************************************************
    Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

    This subroutine assumes that:
    * A*ScaleA is well scaled
    * A is well-conditioned, so no zero divisions or overflow may occur

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void spdbasiccholeskysolve(const ap::template_2d_array< amp::ampf<Precision> >& cha,
        amp::ampf<Precision> sqrtscalea,
        int n,
        bool isupper,
        ap::template_1d_array< amp::ampf<Precision> >& xb,
        ap::template_1d_array< amp::ampf<Precision> >& tmp)
    {
        int i;
        amp::ampf<Precision> v;


        
        //
        // A = L*L' or A=U'*U
        //
        if( isupper )
        {
            
            //
            // Solve U'*y=b first.
            //
            for(i=0; i<=n-1; i++)
            {
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
                if( i<n-1 )
                {
                    v = xb(i);
                    amp::vmove(tmp.getvector(i+1, n-1), cha.getrow(i, i+1, n-1), sqrtscalea);
                    amp::vsub(xb.getvector(i+1, n-1), tmp.getvector(i+1, n-1), v);
                }
            }
            
            //
            // Solve U*x=y then.
            //
            for(i=n-1; i>=0; i--)
            {
                if( i<n-1 )
                {
                    amp::vmove(tmp.getvector(i+1, n-1), cha.getrow(i, i+1, n-1), sqrtscalea);
                    v = amp::vdotproduct(tmp.getvector(i+1, n-1), xb.getvector(i+1, n-1));
                    xb(i) = xb(i)-v;
                }
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            }
        }
        else
        {
            
            //
            // Solve L*y=b first
            //
            for(i=0; i<=n-1; i++)
            {
                if( i>0 )
                {
                    amp::vmove(tmp.getvector(0, i-1), cha.getrow(i, 0, i-1), sqrtscalea);
                    v = amp::vdotproduct(tmp.getvector(0, i-1), xb.getvector(0, i-1));
                    xb(i) = xb(i)-v;
                }
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            }
            
            //
            // Solve L'*x=y then.
            //
            for(i=n-1; i>=0; i--)
            {
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
                if( i>0 )
                {
                    v = xb(i);
                    amp::vmove(tmp.getvector(0, i-1), cha.getrow(i, 0, i-1), sqrtscalea);
                    amp::vsub(xb.getvector(0, i-1), tmp.getvector(0, i-1), v);
                }
            }
        }
    }


    /*************************************************************************
    Basic LU solver for ScaleA*PLU*x = y.

    This subroutine assumes that:
    * L is well-scaled, and it is U which needs scaling by ScaleA.
    * A=PLU is well-conditioned, so no zero divisions or overflow may occur

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void cbasiclusolve(const ap::template_2d_array< amp::campf<Precision> >& lua,
        const ap::template_1d_array< int >& p,
        amp::ampf<Precision> scalea,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& xb,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        amp::campf<Precision> v;
        int i_;


        for(i=0; i<=n-1; i++)
        {
            if( p(i)!=i )
            {
                v = xb(i);
                xb(i) = xb(p(i));
                xb(p(i)) = v;
            }
        }
        for(i=1; i<=n-1; i++)
        {
            v = 0.0;
            for(i_=0; i_<=i-1;i_++)
            {
                v += lua(i,i_)*xb(i_);
            }
            xb(i) = xb(i)-v;
        }
        xb(n-1) = xb(n-1)/(scalea*lua(n-1,n-1));
        for(i=n-2; i>=0; i--)
        {
            for(i_=i+1; i_<=n-1;i_++)
            {
                tmp(i_) = scalea*lua(i,i_);
            }
            v = 0.0;
            for(i_=i+1; i_<=n-1;i_++)
            {
                v += tmp(i_)*xb(i_);
            }
            xb(i) = (xb(i)-v)/(scalea*lua(i,i));
        }
    }


    /*************************************************************************
    Basic Cholesky solver for ScaleA*Cholesky(A)'*x = y.

    This subroutine assumes that:
    * A*ScaleA is well scaled
    * A is well-conditioned, so no zero divisions or overflow may occur

      -- ALGLIB --
         Copyright 27.01.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void hpdbasiccholeskysolve(const ap::template_2d_array< amp::campf<Precision> >& cha,
        amp::ampf<Precision> sqrtscalea,
        int n,
        bool isupper,
        ap::template_1d_array< amp::campf<Precision> >& xb,
        ap::template_1d_array< amp::campf<Precision> >& tmp)
    {
        int i;
        amp::campf<Precision> v;
        int i_;


        
        //
        // A = L*L' or A=U'*U
        //
        if( isupper )
        {
            
            //
            // Solve U'*y=b first.
            //
            for(i=0; i<=n-1; i++)
            {
                xb(i) = xb(i)/(sqrtscalea*amp::conj<Precision>(cha(i,i)));
                if( i<n-1 )
                {
                    v = xb(i);
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        tmp(i_) = sqrtscalea*amp::conj(cha(i,i_));
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        xb(i_) = xb(i_) - v*tmp(i_);
                    }
                }
            }
            
            //
            // Solve U*x=y then.
            //
            for(i=n-1; i>=0; i--)
            {
                if( i<n-1 )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        tmp(i_) = sqrtscalea*cha(i,i_);
                    }
                    v = 0.0;
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        v += tmp(i_)*xb(i_);
                    }
                    xb(i) = xb(i)-v;
                }
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            }
        }
        else
        {
            
            //
            // Solve L*y=b first
            //
            for(i=0; i<=n-1; i++)
            {
                if( i>0 )
                {
                    for(i_=0; i_<=i-1;i_++)
                    {
                        tmp(i_) = sqrtscalea*cha(i,i_);
                    }
                    v = 0.0;
                    for(i_=0; i_<=i-1;i_++)
                    {
                        v += tmp(i_)*xb(i_);
                    }
                    xb(i) = xb(i)-v;
                }
                xb(i) = xb(i)/(sqrtscalea*cha(i,i));
            }
            
            //
            // Solve L'*x=y then.
            //
            for(i=n-1; i>=0; i--)
            {
                xb(i) = xb(i)/(sqrtscalea*amp::conj<Precision>(cha(i,i)));
                if( i>0 )
                {
                    v = xb(i);
                    for(i_=0; i_<=i-1;i_++)
                    {
                        tmp(i_) = sqrtscalea*amp::conj(cha(i,i_));
                    }
                    for(i_=0; i_<=i-1;i_++)
                    {
                        xb(i_) = xb(i_) - v*tmp(i_);
                    }
                }
            }
        }
    }
} // namespace

#endif
