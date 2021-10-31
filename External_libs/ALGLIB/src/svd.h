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

#ifndef _svd_h
#define _svd_h

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
namespace svd
{
    template<unsigned int Precision>
    bool rmatrixsvd(ap::template_2d_array< amp::ampf<Precision> > a,
        int m,
        int n,
        int uneeded,
        int vtneeded,
        int additionalmemory,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& u,
        ap::template_2d_array< amp::ampf<Precision> >& vt);


    /*************************************************************************
    Singular value decomposition of a rectangular matrix.

    The algorithm calculates the singular value decomposition of a matrix of
    size MxN: A = U * S * V^T

    The algorithm finds the singular values and, optionally, matrices U and V^T.
    The algorithm can find both first min(M,N) columns of matrix U and rows of
    matrix V^T (singular vectors), and matrices U and V^T wholly (of sizes MxM
    and NxN respectively).

    Take into account that the subroutine does not return matrix V but V^T.

    Input parameters:
        A           -   matrix to be decomposed.
                        Array whose indexes range within [0..M-1, 0..N-1].
        M           -   number of rows in matrix A.
        N           -   number of columns in matrix A.
        UNeeded     -   0, 1 or 2. See the description of the parameter U.
        VTNeeded    -   0, 1 or 2. See the description of the parameter VT.
        AdditionalMemory -
                        If the parameter:
                         * equals 0, the algorithm doesn’t use additional
                           memory (lower requirements, lower performance).
                         * equals 1, the algorithm uses additional
                           memory of size min(M,N)*min(M,N) of real numbers.
                           It often speeds up the algorithm.
                         * equals 2, the algorithm uses additional
                           memory of size M*min(M,N) of real numbers.
                           It allows to get a maximum performance.
                        The recommended value of the parameter is 2.

    Output parameters:
        W           -   contains singular values in descending order.
        U           -   if UNeeded=0, U isn't changed, the left singular vectors
                        are not calculated.
                        if Uneeded=1, U contains left singular vectors (first
                        min(M,N) columns of matrix U). Array whose indexes range
                        within [0..M-1, 0..Min(M,N)-1].
                        if UNeeded=2, U contains matrix U wholly. Array whose
                        indexes range within [0..M-1, 0..M-1].
        VT          -   if VTNeeded=0, VT isn’t changed, the right singular vectors
                        are not calculated.
                        if VTNeeded=1, VT contains right singular vectors (first
                        min(M,N) rows of matrix V^T). Array whose indexes range
                        within [0..min(M,N)-1, 0..N-1].
                        if VTNeeded=2, VT contains matrix V^T wholly. Array whose
                        indexes range within [0..N-1, 0..N-1].

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixsvd(ap::template_2d_array< amp::ampf<Precision> > a,
        int m,
        int n,
        int uneeded,
        int vtneeded,
        int additionalmemory,
        ap::template_1d_array< amp::ampf<Precision> >& w,
        ap::template_2d_array< amp::ampf<Precision> >& u,
        ap::template_2d_array< amp::ampf<Precision> >& vt)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > tauq;
        ap::template_1d_array< amp::ampf<Precision> > taup;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > e;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_2d_array< amp::ampf<Precision> > t2;
        bool isupper;
        int minmn;
        int ncu;
        int nrvt;
        int nru;
        int ncvt;
        int i;
        int j;


        result = true;
        if( m==0 || n==0 )
        {
            return result;
        }
        ap::ap_error::make_assertion(uneeded>=0 && uneeded<=2);
        ap::ap_error::make_assertion(vtneeded>=0 && vtneeded<=2);
        ap::ap_error::make_assertion(additionalmemory>=0 && additionalmemory<=2);
        
        //
        // initialize
        //
        minmn = ap::minint(m, n);
        w.setbounds(1, minmn);
        ncu = 0;
        nru = 0;
        if( uneeded==1 )
        {
            nru = m;
            ncu = minmn;
            u.setbounds(0, nru-1, 0, ncu-1);
        }
        if( uneeded==2 )
        {
            nru = m;
            ncu = m;
            u.setbounds(0, nru-1, 0, ncu-1);
        }
        nrvt = 0;
        ncvt = 0;
        if( vtneeded==1 )
        {
            nrvt = minmn;
            ncvt = n;
            vt.setbounds(0, nrvt-1, 0, ncvt-1);
        }
        if( vtneeded==2 )
        {
            nrvt = n;
            ncvt = n;
            vt.setbounds(0, nrvt-1, 0, ncvt-1);
        }
        
        //
        // M much larger than N
        // Use bidiagonal reduction with QR-decomposition
        //
        if( m>amp::ampf<Precision>("1.6")*n )
        {
            if( uneeded==0 )
            {
                
                //
                // No left singular vectors to be computed
                //
                ortfac::rmatrixqr<Precision>(a, m, n, tau);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=i-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
                ortfac::rmatrixbd<Precision>(a, n, n, tauq, taup);
                ortfac::rmatrixbdunpackpt<Precision>(a, n, n, taup, nrvt, vt);
                ortfac::rmatrixbdunpackdiagonals<Precision>(a, n, n, isupper, w, e);
                result = bdsvd::rmatrixbdsvd<Precision>(w, e, n, isupper, false, u, 0, a, 0, vt, ncvt);
                return result;
            }
            else
            {
                
                //
                // Left singular vectors (may be full matrix U) to be computed
                //
                ortfac::rmatrixqr<Precision>(a, m, n, tau);
                ortfac::rmatrixqrunpackq<Precision>(a, m, n, tau, ncu, u);
                for(i=0; i<=n-1; i++)
                {
                    for(j=0; j<=i-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
                ortfac::rmatrixbd<Precision>(a, n, n, tauq, taup);
                ortfac::rmatrixbdunpackpt<Precision>(a, n, n, taup, nrvt, vt);
                ortfac::rmatrixbdunpackdiagonals<Precision>(a, n, n, isupper, w, e);
                if( additionalmemory<1 )
                {
                    
                    //
                    // No additional memory can be used
                    //
                    ortfac::rmatrixbdmultiplybyq<Precision>(a, n, n, tauq, u, m, n, true, false);
                    result = bdsvd::rmatrixbdsvd<Precision>(w, e, n, isupper, false, u, m, a, 0, vt, ncvt);
                }
                else
                {
                    
                    //
                    // Large U. Transforming intermediate matrix T2
                    //
                    work.setbounds(1, ap::maxint(m, n));
                    ortfac::rmatrixbdunpackq<Precision>(a, n, n, tauq, n, t2);
                    blas::copymatrix<Precision>(u, 0, m-1, 0, n-1, a, 0, m-1, 0, n-1);
                    blas::inplacetranspose<Precision>(t2, 0, n-1, 0, n-1, work);
                    result = bdsvd::rmatrixbdsvd<Precision>(w, e, n, isupper, false, u, 0, t2, n, vt, ncvt);
                    blas::matrixmatrixmultiply<Precision>(a, 0, m-1, 0, n-1, false, t2, 0, n-1, 0, n-1, true, amp::ampf<Precision>("1.0"), u, 0, m-1, 0, n-1, amp::ampf<Precision>("0.0"), work);
                }
                return result;
            }
        }
        
        //
        // N much larger than M
        // Use bidiagonal reduction with LQ-decomposition
        //
        if( n>amp::ampf<Precision>("1.6")*m )
        {
            if( vtneeded==0 )
            {
                
                //
                // No right singular vectors to be computed
                //
                ortfac::rmatrixlq<Precision>(a, m, n, tau);
                for(i=0; i<=m-1; i++)
                {
                    for(j=i+1; j<=m-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
                ortfac::rmatrixbd<Precision>(a, m, m, tauq, taup);
                ortfac::rmatrixbdunpackq<Precision>(a, m, m, tauq, ncu, u);
                ortfac::rmatrixbdunpackdiagonals<Precision>(a, m, m, isupper, w, e);
                work.setbounds(1, m);
                blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
                result = bdsvd::rmatrixbdsvd<Precision>(w, e, m, isupper, false, a, 0, u, nru, vt, 0);
                blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
                return result;
            }
            else
            {
                
                //
                // Right singular vectors (may be full matrix VT) to be computed
                //
                ortfac::rmatrixlq<Precision>(a, m, n, tau);
                ortfac::rmatrixlqunpackq<Precision>(a, m, n, tau, nrvt, vt);
                for(i=0; i<=m-1; i++)
                {
                    for(j=i+1; j<=m-1; j++)
                    {
                        a(i,j) = 0;
                    }
                }
                ortfac::rmatrixbd<Precision>(a, m, m, tauq, taup);
                ortfac::rmatrixbdunpackq<Precision>(a, m, m, tauq, ncu, u);
                ortfac::rmatrixbdunpackdiagonals<Precision>(a, m, m, isupper, w, e);
                work.setbounds(1, ap::maxint(m, n));
                blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
                if( additionalmemory<1 )
                {
                    
                    //
                    // No additional memory available
                    //
                    ortfac::rmatrixbdmultiplybyp<Precision>(a, m, m, taup, vt, m, n, false, true);
                    result = bdsvd::rmatrixbdsvd<Precision>(w, e, m, isupper, false, a, 0, u, nru, vt, n);
                }
                else
                {
                    
                    //
                    // Large VT. Transforming intermediate matrix T2
                    //
                    ortfac::rmatrixbdunpackpt<Precision>(a, m, m, taup, m, t2);
                    result = bdsvd::rmatrixbdsvd<Precision>(w, e, m, isupper, false, a, 0, u, nru, t2, m);
                    blas::copymatrix<Precision>(vt, 0, m-1, 0, n-1, a, 0, m-1, 0, n-1);
                    blas::matrixmatrixmultiply<Precision>(t2, 0, m-1, 0, m-1, false, a, 0, m-1, 0, n-1, false, amp::ampf<Precision>("1.0"), vt, 0, m-1, 0, n-1, amp::ampf<Precision>("0.0"), work);
                }
                blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
                return result;
            }
        }
        
        //
        // M<=N
        // We can use inplace transposition of U to get rid of columnwise operations
        //
        if( m<=n )
        {
            ortfac::rmatrixbd<Precision>(a, m, n, tauq, taup);
            ortfac::rmatrixbdunpackq<Precision>(a, m, n, tauq, ncu, u);
            ortfac::rmatrixbdunpackpt<Precision>(a, m, n, taup, nrvt, vt);
            ortfac::rmatrixbdunpackdiagonals<Precision>(a, m, n, isupper, w, e);
            work.setbounds(1, m);
            blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
            result = bdsvd::rmatrixbdsvd<Precision>(w, e, minmn, isupper, false, a, 0, u, nru, vt, ncvt);
            blas::inplacetranspose<Precision>(u, 0, nru-1, 0, ncu-1, work);
            return result;
        }
        
        //
        // Simple bidiagonal reduction
        //
        ortfac::rmatrixbd<Precision>(a, m, n, tauq, taup);
        ortfac::rmatrixbdunpackq<Precision>(a, m, n, tauq, ncu, u);
        ortfac::rmatrixbdunpackpt<Precision>(a, m, n, taup, nrvt, vt);
        ortfac::rmatrixbdunpackdiagonals<Precision>(a, m, n, isupper, w, e);
        if( additionalmemory<2 || uneeded==0 )
        {
            
            //
            // We cant use additional memory or there is no need in such operations
            //
            result = bdsvd::rmatrixbdsvd<Precision>(w, e, minmn, isupper, false, u, nru, a, 0, vt, ncvt);
        }
        else
        {
            
            //
            // We can use additional memory
            //
            t2.setbounds(0, minmn-1, 0, m-1);
            blas::copyandtranspose<Precision>(u, 0, m-1, 0, minmn-1, t2, 0, minmn-1, 0, m-1);
            result = bdsvd::rmatrixbdsvd<Precision>(w, e, minmn, isupper, false, u, 0, t2, m, vt, ncvt);
            blas::copyandtranspose<Precision>(t2, 0, minmn-1, 0, m-1, u, 0, m-1, 0, minmn-1);
        }
        return result;
    }
} // namespace

#endif
