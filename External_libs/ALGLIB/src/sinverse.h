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

#ifndef _sinverse_h
#define _sinverse_h

#include "ap.h"
#include "amp.h"
#include "sblas.h"
#include "ldlt.h"
namespace sinverse
{
    template<unsigned int Precision>
    bool smatrixldltinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        bool isupper);
    template<unsigned int Precision>
    bool smatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper);
    template<unsigned int Precision>
    bool inverseldlt(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        bool isupper);
    template<unsigned int Precision>
    bool inversesymmetricindefinite(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper);


    /*************************************************************************
    Inversion of a symmetric indefinite matrix

    The algorithm gets an LDLT-decomposition as an input, generates matrix A^-1
    and saves the lower or upper triangle of an inverse matrix depending on the
    input (U*D*U' or L*D*L').

    Input parameters:
        A       -   LDLT-decomposition of the matrix,
                    Output of subroutine SMatrixLDLT.
        N       -   size of matrix A.
        IsUpper -   storage format. If IsUpper = True, then the symmetric matrix
                    is given as decomposition A = U*D*U' and this decomposition
                    is stored in the upper triangle of matrix A and on the main
                    diagonal, and the lower triangle of matrix A is not used.
        Pivots  -   a table of permutations, output of subroutine SMatrixLDLT.

    Output parameters:
        A       -   inverse of the matrix, whose LDLT-decomposition was stored
                    in matrix A as a subroutine input.
                    Array with elements [0..N-1, 0..N-1].
                    If IsUpper = True, then A contains the upper triangle of
                    matrix A^-1, and the elements below the main diagonal are
                    not used nor changed. The same applies if IsUpper = False.

    Result:
        True, if the matrix is not singular.
        False, if the matrix is singular and could not be inverted.

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         March 31, 1993
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixldltinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_1d_array< amp::ampf<Precision> > work2;
        int i;
        int k;
        int kp;
        int kstep;
        amp::ampf<Precision> ak;
        amp::ampf<Precision> akkp1;
        amp::ampf<Precision> akp1;
        amp::ampf<Precision> d;
        amp::ampf<Precision> t;
        amp::ampf<Precision> temp;
        int km1;
        int kp1;
        int l;
        int i1;
        int i2;
        amp::ampf<Precision> v;


        work.setbounds(1, n);
        work2.setbounds(1, n);
        result = true;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        
        //
        // Check that the diagonal matrix D is nonsingular.
        //
        for(i=0; i<=n-1; i++)
        {
            if( pivots(i)>=0 && a(i,i)==0 )
            {
                result = false;
                return result;
            }
        }
        if( isupper )
        {
            
            //
            // Compute inv(A) from the factorization A = U*D*U'.
            //
            // K+1 is the main loop index, increasing from 1 to N in steps of
            // 1 or 2, depending on the size of the diagonal blocks.
            //
            k = 0;
            while( k<=n-1 )
            {
                if( pivots(k)>=0 )
                {
                    
                    //
                    // 1 x 1 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    a(k,k) = 1/a(k,k);
                    
                    //
                    // Compute column K+1 of the inverse.
                    //
                    if( k>0 )
                    {
                        amp::vmove(work.getvector(1, k), a.getcolumn(k, 0, k-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 1-1, k+1-1-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, 0, k-1), work2.getvector(1, k));
                        v = amp::vdotproduct(work2.getvector(1, k), work.getvector(1, k));
                        a(k,k) = a(k,k)-v;
                    }
                    kstep = 1;
                }
                else
                {
                    
                    //
                    // 2 x 2 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    t = amp::abs<Precision>(a(k,k+1));
                    ak = a(k,k)/t;
                    akp1 = a(k+1,k+1)/t;
                    akkp1 = a(k,k+1)/t;
                    d = t*(ak*akp1-1);
                    a(k,k) = akp1/d;
                    a(k+1,k+1) = ak/d;
                    a(k,k+1) = -akkp1/d;
                    
                    //
                    // Compute columns K+1 and K+1+1 of the inverse.
                    //
                    if( k>0 )
                    {
                        amp::vmove(work.getvector(1, k), a.getcolumn(k, 0, k-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 0, k-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, 0, k-1), work2.getvector(1, k));
                        v = amp::vdotproduct(work.getvector(1, k), work2.getvector(1, k));
                        a(k,k) = a(k,k)-v;
                        v = amp::vdotproduct(a.getcolumn(k, 0, k-1), a.getcolumn(k+1, 0, k-1));
                        a(k,k+1) = a(k,k+1)-v;
                        amp::vmove(work.getvector(1, k), a.getcolumn(k+1, 0, k-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 0, k-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k+1, 0, k-1), work2.getvector(1, k));
                        v = amp::vdotproduct(work.getvector(1, k), work2.getvector(1, k));
                        a(k+1,k+1) = a(k+1,k+1)-v;
                    }
                    kstep = 2;
                }
                if( pivots(k)>=0 )
                {
                    kp = pivots(k);
                }
                else
                {
                    kp = n+pivots(k);
                }
                if( kp!=k )
                {
                    
                    //
                    // Interchange rows and columns K and KP in the leading
                    // submatrix
                    //
                    amp::vmove(work.getvector(1, kp), a.getcolumn(k, 0, kp-1));
                    amp::vmove(a.getcolumn(k, 0, kp-1), a.getcolumn(kp, 0, kp-1));
                    amp::vmove(a.getcolumn(kp, 0, kp-1), work.getvector(1, kp));
                    amp::vmove(work.getvector(1, k-1-kp), a.getcolumn(k, kp+1, k-1));
                    amp::vmove(a.getcolumn(k, kp+1, k-1), a.getrow(kp, kp+1, k-1));
                    amp::vmove(a.getrow(kp, kp+1, k-1), work.getvector(1, k-1-kp));
                    temp = a(k,k);
                    a(k,k) = a(kp,kp);
                    a(kp,kp) = temp;
                    if( kstep==2 )
                    {
                        temp = a(k,k+1);
                        a(k,k+1) = a(kp,k+1);
                        a(kp,k+1) = temp;
                    }
                }
                k = k+kstep;
            }
        }
        else
        {
            
            //
            // Compute inv(A) from the factorization A = L*D*L'.
            //
            // K is the main loop index, increasing from 0 to N-1 in steps of
            // 1 or 2, depending on the size of the diagonal blocks.
            //
            k = n-1;
            while( k>=0 )
            {
                if( pivots(k)>=0 )
                {
                    
                    //
                    // 1 x 1 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    a(k,k) = 1/a(k,k);
                    
                    //
                    // Compute column K+1 of the inverse.
                    //
                    if( k<n-1 )
                    {
                        amp::vmove(work.getvector(1, n-k-1), a.getcolumn(k, k+1, n-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, k+1, n-1), work2.getvector(1, n-k-1));
                        v = amp::vdotproduct(work.getvector(1, n-k-1), work2.getvector(1, n-k-1));
                        a(k,k) = a(k,k)-v;
                    }
                    kstep = 1;
                }
                else
                {
                    
                    //
                    // 2 x 2 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    t = amp::abs<Precision>(a(k,k-1));
                    ak = a(k-1,k-1)/t;
                    akp1 = a(k,k)/t;
                    akkp1 = a(k,k-1)/t;
                    d = t*(ak*akp1-1);
                    a(k-1,k-1) = akp1/d;
                    a(k,k) = ak/d;
                    a(k,k-1) = -akkp1/d;
                    
                    //
                    // Compute columns K+1-1 and K+1 of the inverse.
                    //
                    if( k<n-1 )
                    {
                        amp::vmove(work.getvector(1, n-k-1), a.getcolumn(k, k+1, n-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, k+1, n-1), work2.getvector(1, n-k-1));
                        v = amp::vdotproduct(work.getvector(1, n-k-1), work2.getvector(1, n-k-1));
                        a(k,k) = a(k,k)-v;
                        v = amp::vdotproduct(a.getcolumn(k, k+1, n-1), a.getcolumn(k-1, k+1, n-1));
                        a(k,k-1) = a(k,k-1)-v;
                        amp::vmove(work.getvector(1, n-k-1), a.getcolumn(k-1, k+1, n-1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k-1, k+1, n-1), work2.getvector(1, n-k-1));
                        v = amp::vdotproduct(work.getvector(1, n-k-1), work2.getvector(1, n-k-1));
                        a(k-1,k-1) = a(k-1,k-1)-v;
                    }
                    kstep = 2;
                }
                if( pivots(k)>=0 )
                {
                    kp = pivots(k);
                }
                else
                {
                    kp = pivots(k)+n;
                }
                if( kp!=k )
                {
                    
                    //
                    // Interchange rows and columns K and KP
                    //
                    if( kp<n-1 )
                    {
                        amp::vmove(work.getvector(1, n-kp-1), a.getcolumn(k, kp+1, n-1));
                        amp::vmove(a.getcolumn(k, kp+1, n-1), a.getcolumn(kp, kp+1, n-1));
                        amp::vmove(a.getcolumn(kp, kp+1, n-1), work.getvector(1, n-kp-1));
                    }
                    amp::vmove(work.getvector(1, kp-k-1), a.getcolumn(k, k+1, kp-1));
                    amp::vmove(a.getcolumn(k, k+1, kp-1), a.getrow(kp, k+1, kp-1));
                    amp::vmove(a.getrow(kp, k+1, kp-1), work.getvector(1, kp-k-1));
                    temp = a(k,k);
                    a(k,k) = a(kp,kp);
                    a(kp,kp) = temp;
                    if( kstep==2 )
                    {
                        temp = a(k,k-1);
                        a(k,k-1) = a(kp,k-1);
                        a(kp,k-1) = temp;
                    }
                }
                k = k-kstep;
            }
        }
        return result;
    }


    /*************************************************************************
    Inversion of a symmetric indefinite matrix

    Given a lower or upper triangle of matrix A, the algorithm generates
    matrix A^-1 and saves the lower or upper triangle depending on the input.

    Input parameters:
        A       -   matrix to be inverted (upper or lower triangle).
                    Array with elements [0..N-1, 0..N-1].
        N       -   size of matrix A.
        IsUpper -   storage format. If IsUpper = True, then the upper
                    triangle of matrix A is given, otherwise the lower
                    triangle is given.

    Output parameters:
        A       -   inverse of matrix A.
                    Array with elements [0..N-1, 0..N-1].
                    If IsUpper = True, then A contains the upper triangle of
                    matrix A^-1, and the elements below the main diagonal are
                    not used nor changed.
                    The same applies if IsUpper = False.

    Result:
        True, if the matrix is not singular.
        False, if the matrix is singular and could not be inverted.

      -- LAPACK routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         March 31, 1993
    *************************************************************************/
    template<unsigned int Precision>
    bool smatrixinverse(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< int > pivots;


        ldlt::smatrixldlt<Precision>(a, n, isupper, pivots);
        result = smatrixldltinverse<Precision>(a, pivots, n, isupper);
        return result;
    }


    template<unsigned int Precision>
    bool inverseldlt(ap::template_2d_array< amp::ampf<Precision> >& a,
        const ap::template_1d_array< int >& pivots,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< amp::ampf<Precision> > work;
        ap::template_1d_array< amp::ampf<Precision> > work2;
        int i;
        int k;
        int kp;
        int kstep;
        amp::ampf<Precision> ak;
        amp::ampf<Precision> akkp1;
        amp::ampf<Precision> akp1;
        amp::ampf<Precision> d;
        amp::ampf<Precision> t;
        amp::ampf<Precision> temp;
        int km1;
        int kp1;
        int l;
        int i1;
        int i2;
        amp::ampf<Precision> v;


        work.setbounds(1, n);
        work2.setbounds(1, n);
        result = true;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return result;
        }
        
        //
        // Check that the diagonal matrix D is nonsingular.
        //
        for(i=1; i<=n; i++)
        {
            if( pivots(i)>0 && a(i,i)==0 )
            {
                result = false;
                return result;
            }
        }
        if( isupper )
        {
            
            //
            // Compute inv(A) from the factorization A = U*D*U'.
            //
            // K is the main loop index, increasing from 1 to N in steps of
            // 1 or 2, depending on the size of the diagonal blocks.
            //
            k = 1;
            while( k<=n )
            {
                if( pivots(k)>0 )
                {
                    
                    //
                    // 1 x 1 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    a(k,k) = 1/a(k,k);
                    
                    //
                    // Compute column K of the inverse.
                    //
                    if( k>1 )
                    {
                        km1 = k-1;
                        amp::vmove(work.getvector(1, km1), a.getcolumn(k, 1, km1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 1, k-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, 1, km1), work2.getvector(1, km1));
                        v = amp::vdotproduct(work2.getvector(1, km1), work.getvector(1, km1));
                        a(k,k) = a(k,k)-v;
                    }
                    kstep = 1;
                }
                else
                {
                    
                    //
                    // 2 x 2 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    t = amp::abs<Precision>(a(k,k+1));
                    ak = a(k,k)/t;
                    akp1 = a(k+1,k+1)/t;
                    akkp1 = a(k,k+1)/t;
                    d = t*(ak*akp1-1);
                    a(k,k) = akp1/d;
                    a(k+1,k+1) = ak/d;
                    a(k,k+1) = -akkp1/d;
                    
                    //
                    // Compute columns K and K+1 of the inverse.
                    //
                    if( k>1 )
                    {
                        km1 = k-1;
                        kp1 = k+1;
                        amp::vmove(work.getvector(1, km1), a.getcolumn(k, 1, km1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 1, k-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, 1, km1), work2.getvector(1, km1));
                        v = amp::vdotproduct(work.getvector(1, km1), work2.getvector(1, km1));
                        a(k,k) = a(k,k)-v;
                        v = amp::vdotproduct(a.getcolumn(k, 1, km1), a.getcolumn(kp1, 1, km1));
                        a(k,k+1) = a(k,k+1)-v;
                        amp::vmove(work.getvector(1, km1), a.getcolumn(kp1, 1, km1));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, 1, k-1, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(kp1, 1, km1), work2.getvector(1, km1));
                        v = amp::vdotproduct(work.getvector(1, km1), work2.getvector(1, km1));
                        a(k+1,k+1) = a(k+1,k+1)-v;
                    }
                    kstep = 2;
                }
                kp = abs(pivots(k));
                if( kp!=k )
                {
                    
                    //
                    // Interchange rows and columns K and KP in the leading
                    // submatrix A(1:k+1,1:k+1)
                    //
                    l = kp-1;
                    amp::vmove(work.getvector(1, l), a.getcolumn(k, 1, l));
                    amp::vmove(a.getcolumn(k, 1, l), a.getcolumn(kp, 1, l));
                    amp::vmove(a.getcolumn(kp, 1, l), work.getvector(1, l));
                    l = k-kp-1;
                    i1 = kp+1;
                    i2 = k-1;
                    amp::vmove(work.getvector(1, l), a.getcolumn(k, i1, i2));
                    amp::vmove(a.getcolumn(k, i1, i2), a.getrow(kp, i1, i2));
                    amp::vmove(a.getrow(kp, i1, i2), work.getvector(1, l));
                    temp = a(k,k);
                    a(k,k) = a(kp,kp);
                    a(kp,kp) = temp;
                    if( kstep==2 )
                    {
                        temp = a(k,k+1);
                        a(k,k+1) = a(kp,k+1);
                        a(kp,k+1) = temp;
                    }
                }
                k = k+kstep;
            }
        }
        else
        {
            
            //
            // Compute inv(A) from the factorization A = L*D*L'.
            //
            // K is the main loop index, increasing from 1 to N in steps of
            // 1 or 2, depending on the size of the diagonal blocks.
            //
            k = n;
            while( k>=1 )
            {
                if( pivots(k)>0 )
                {
                    
                    //
                    // 1 x 1 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    a(k,k) = 1/a(k,k);
                    
                    //
                    // Compute column K of the inverse.
                    //
                    if( k<n )
                    {
                        kp1 = k+1;
                        km1 = k-1;
                        l = n-k;
                        amp::vmove(work.getvector(1, l), a.getcolumn(k, kp1, n));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, kp1, n), work2.getvector(1, l));
                        v = amp::vdotproduct(work.getvector(1, l), work2.getvector(1, l));
                        a(k,k) = a(k,k)-v;
                    }
                    kstep = 1;
                }
                else
                {
                    
                    //
                    // 2 x 2 diagonal block
                    //
                    // Invert the diagonal block.
                    //
                    t = amp::abs<Precision>(a(k,k-1));
                    ak = a(k-1,k-1)/t;
                    akp1 = a(k,k)/t;
                    akkp1 = a(k,k-1)/t;
                    d = t*(ak*akp1-1);
                    a(k-1,k-1) = akp1/d;
                    a(k,k) = ak/d;
                    a(k,k-1) = -akkp1/d;
                    
                    //
                    // Compute columns K-1 and K of the inverse.
                    //
                    if( k<n )
                    {
                        kp1 = k+1;
                        km1 = k-1;
                        l = n-k;
                        amp::vmove(work.getvector(1, l), a.getcolumn(k, kp1, n));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(k, kp1, n), work2.getvector(1, l));
                        v = amp::vdotproduct(work.getvector(1, l), work2.getvector(1, l));
                        a(k,k) = a(k,k)-v;
                        v = amp::vdotproduct(a.getcolumn(k, kp1, n), a.getcolumn(km1, kp1, n));
                        a(k,k-1) = a(k,k-1)-v;
                        amp::vmove(work.getvector(1, l), a.getcolumn(km1, kp1, n));
                        sblas::symmetricmatrixvectormultiply<Precision>(a, isupper, k+1, n, work, amp::ampf<Precision>(-1), work2);
                        amp::vmove(a.getcolumn(km1, kp1, n), work2.getvector(1, l));
                        v = amp::vdotproduct(work.getvector(1, l), work2.getvector(1, l));
                        a(k-1,k-1) = a(k-1,k-1)-v;
                    }
                    kstep = 2;
                }
                kp = abs(pivots(k));
                if( kp!=k )
                {
                    
                    //
                    // Interchange rows and columns K and KP in the trailing
                    // submatrix A(k-1:n,k-1:n)
                    //
                    if( kp<n )
                    {
                        l = n-kp;
                        kp1 = kp+1;
                        amp::vmove(work.getvector(1, l), a.getcolumn(k, kp1, n));
                        amp::vmove(a.getcolumn(k, kp1, n), a.getcolumn(kp, kp1, n));
                        amp::vmove(a.getcolumn(kp, kp1, n), work.getvector(1, l));
                    }
                    l = kp-k-1;
                    i1 = k+1;
                    i2 = kp-1;
                    amp::vmove(work.getvector(1, l), a.getcolumn(k, i1, i2));
                    amp::vmove(a.getcolumn(k, i1, i2), a.getrow(kp, i1, i2));
                    amp::vmove(a.getrow(kp, i1, i2), work.getvector(1, l));
                    temp = a(k,k);
                    a(k,k) = a(kp,kp);
                    a(kp,kp) = temp;
                    if( kstep==2 )
                    {
                        temp = a(k,k-1);
                        a(k,k-1) = a(kp,k-1);
                        a(kp,k-1) = temp;
                    }
                }
                k = k-kstep;
            }
        }
        return result;
    }


    template<unsigned int Precision>
    bool inversesymmetricindefinite(ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        bool isupper)
    {
        bool result;
        ap::template_1d_array< int > pivots;


        ldlt::ldltdecomposition<Precision>(a, n, isupper, pivots);
        result = inverseldlt<Precision>(a, pivots, n, isupper);
        return result;
    }
} // namespace

#endif
