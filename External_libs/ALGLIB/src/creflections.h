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

#ifndef _creflections_h
#define _creflections_h

#include "ap.h"
#include "amp.h"
namespace creflections
{
    template<unsigned int Precision>
    void complexgeneratereflection(ap::template_1d_array< amp::campf<Precision> >& x,
        int n,
        amp::campf<Precision>& tau);
    template<unsigned int Precision>
    void complexapplyreflectionfromtheleft(ap::template_2d_array< amp::campf<Precision> >& c,
        amp::campf<Precision> tau,
        const ap::template_1d_array< amp::campf<Precision> >& v,
        int m1,
        int m2,
        int n1,
        int n2,
        ap::template_1d_array< amp::campf<Precision> >& work);
    template<unsigned int Precision>
    void complexapplyreflectionfromtheright(ap::template_2d_array< amp::campf<Precision> >& c,
        amp::campf<Precision> tau,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int m1,
        int m2,
        int n1,
        int n2,
        ap::template_1d_array< amp::campf<Precision> >& work);


    /*************************************************************************
    Generation of an elementary complex reflection transformation

    The subroutine generates elementary complex reflection H of  order  N,  so
    that, for a given X, the following equality holds true:

         ( X(1) )   ( Beta )
    H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
         ( X(n) )   (  0   )

    where

                  ( V(1) )
    H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
                  ( V(n) )

    where the first component of vector V equals 1.

    Input parameters:
        X   -   vector. Array with elements [1..N].
        N   -   reflection order.

    Output parameters:
        X   -   components from 2 to N are replaced by vector V.
                The first component is replaced with parameter Beta.
        Tau -   scalar value Tau.

    This subroutine is the modification of CLARFG subroutines  from the LAPACK
    library. It has similar functionality except for the fact that it  doesn’t
    handle errors when intermediate results cause an overflow.

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         September 30, 1994
    *************************************************************************/
    template<unsigned int Precision>
    void complexgeneratereflection(ap::template_1d_array< amp::campf<Precision> >& x,
        int n,
        amp::campf<Precision>& tau)
    {
        int j;
        amp::campf<Precision> alpha;
        amp::ampf<Precision> alphi;
        amp::ampf<Precision> alphr;
        amp::ampf<Precision> beta;
        amp::ampf<Precision> xnorm;
        amp::ampf<Precision> mx;
        amp::campf<Precision> t;
        amp::ampf<Precision> s;
        amp::campf<Precision> v;
        int i_;


        if( n<=0 )
        {
            tau = 0;
            return;
        }
        
        //
        // Scale if needed (to avoid overflow/underflow during intermediate
        // calculations).
        //
        mx = 0;
        for(j=1; j<=n; j++)
        {
            mx = amp::maximum<Precision>(amp::abscomplex<Precision>(x(j)), mx);
        }
        s = 1;
        if( mx!=0 )
        {
            if( mx<1 )
            {
                s = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMinNumber());
                v = 1/s;
                for(i_=1; i_<=n;i_++)
                {
                    x(i_) = v*x(i_);
                }
            }
            else
            {
                s = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
                v = 1/s;
                for(i_=1; i_<=n;i_++)
                {
                    x(i_) = v*x(i_);
                }
            }
        }
        
        //
        // calculate
        //
        alpha = x(1);
        mx = 0;
        for(j=2; j<=n; j++)
        {
            mx = amp::maximum<Precision>(amp::abscomplex<Precision>(x(j)), mx);
        }
        xnorm = 0;
        if( mx!=0 )
        {
            for(j=2; j<=n; j++)
            {
                t = x(j)/mx;
                xnorm = xnorm+(t*amp::conj<Precision>(t)).x;
            }
            xnorm = amp::sqrt<Precision>(xnorm)*mx;
        }
        alphr = alpha.x;
        alphi = alpha.y;
        if( xnorm==0 && alphi==0 )
        {
            tau = 0;
            x(1) = x(1)*s;
            return;
        }
        mx = amp::maximum<Precision>(amp::abs<Precision>(alphr), amp::abs<Precision>(alphi));
        mx = amp::maximum<Precision>(mx, amp::abs<Precision>(xnorm));
        beta = -mx*amp::sqrt<Precision>(amp::sqr<Precision>(alphr/mx)+amp::sqr<Precision>(alphi/mx)+amp::sqr<Precision>(xnorm/mx));
        if( alphr<0 )
        {
            beta = -beta;
        }
        tau.x = (beta-alphr)/beta;
        tau.y = -alphi/beta;
        alpha = 1/(alpha-beta);
        if( n>1 )
        {
            for(i_=2; i_<=n;i_++)
            {
                x(i_) = alpha*x(i_);
            }
        }
        alpha = beta;
        x(1) = alpha;
        
        //
        // Scale back
        //
        x(1) = x(1)*s;
    }


    /*************************************************************************
    Application of an elementary reflection to a rectangular matrix of size MxN

    The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
    transformation  which  is  given  by  column  V  and  scalar  Tau (see the
    description of the GenerateReflection). Not the whole matrix  but  only  a
    part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
    the elements of this submatrix are changed.

    Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
    multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.

    Input parameters:
        C       -   matrix to be transformed.
        Tau     -   scalar defining transformation.
        V       -   column defining transformation.
                    Array whose index ranges within [1..M2-M1+1]
        M1, M2  -   range of rows to be transformed.
        N1, N2  -   range of columns to be transformed.
        WORK    -   working array whose index goes from N1 to N2.

    Output parameters:
        C       -   the result of multiplying the input matrix C by the
                    transformation matrix which is given by Tau and V.
                    If N1>N2 or M1>M2, C is not modified.

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         September 30, 1994
    *************************************************************************/
    template<unsigned int Precision>
    void complexapplyreflectionfromtheleft(ap::template_2d_array< amp::campf<Precision> >& c,
        amp::campf<Precision> tau,
        const ap::template_1d_array< amp::campf<Precision> >& v,
        int m1,
        int m2,
        int n1,
        int n2,
        ap::template_1d_array< amp::campf<Precision> >& work)
    {
        amp::campf<Precision> t;
        int i;
        int vm;
        int i_;


        if( tau==0 || n1>n2 || m1>m2 )
        {
            return;
        }
        
        //
        // w := C^T * conj(v)
        //
        vm = m2-m1+1;
        for(i=n1; i<=n2; i++)
        {
            work(i) = 0;
        }
        for(i=m1; i<=m2; i++)
        {
            t = amp::conj<Precision>(v(i+1-m1));
            for(i_=n1; i_<=n2;i_++)
            {
                work(i_) = work(i_) + t*c(i,i_);
            }
        }
        
        //
        // C := C - tau * v * w^T
        //
        for(i=m1; i<=m2; i++)
        {
            t = v(i-m1+1)*tau;
            for(i_=n1; i_<=n2;i_++)
            {
                c(i,i_) = c(i,i_) - t*work(i_);
            }
        }
    }


    /*************************************************************************
    Application of an elementary reflection to a rectangular matrix of size MxN

    The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
    transformation  which  is  given  by  column  V  and  scalar  Tau (see the
    description  of  the  GenerateReflection). Not the whole matrix but only a
    part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
    Only the elements of this submatrix are changed.

    Input parameters:
        C       -   matrix to be transformed.
        Tau     -   scalar defining transformation.
        V       -   column defining transformation.
                    Array whose index ranges within [1..N2-N1+1]
        M1, M2  -   range of rows to be transformed.
        N1, N2  -   range of columns to be transformed.
        WORK    -   working array whose index goes from M1 to M2.

    Output parameters:
        C       -   the result of multiplying the input matrix C by the
                    transformation matrix which is given by Tau and V.
                    If N1>N2 or M1>M2, C is not modified.

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         September 30, 1994
    *************************************************************************/
    template<unsigned int Precision>
    void complexapplyreflectionfromtheright(ap::template_2d_array< amp::campf<Precision> >& c,
        amp::campf<Precision> tau,
        ap::template_1d_array< amp::campf<Precision> >& v,
        int m1,
        int m2,
        int n1,
        int n2,
        ap::template_1d_array< amp::campf<Precision> >& work)
    {
        amp::campf<Precision> t;
        int i;
        int vm;
        int i_;
        int i1_;


        if( tau==0 || n1>n2 || m1>m2 )
        {
            return;
        }
        
        //
        // w := C * v
        //
        vm = n2-n1+1;
        for(i=m1; i<=m2; i++)
        {
            i1_ = (1)-(n1);
            t = 0.0;
            for(i_=n1; i_<=n2;i_++)
            {
                t += c(i,i_)*v(i_+i1_);
            }
            work(i) = t;
        }
        
        //
        // C := C - w * conj(v^T)
        //
        for(i_=1; i_<=vm;i_++)
        {
            v(i_) = amp::conj(v(i_));
        }
        for(i=m1; i<=m2; i++)
        {
            t = work(i)*tau;
            i1_ = (1) - (n1);
            for(i_=n1; i_<=n2;i_++)
            {
                c(i,i_) = c(i,i_) - t*v(i_+i1_);
            }
        }
        for(i_=1; i_<=vm;i_++)
        {
            v(i_) = amp::conj(v(i_));
        }
    }
} // namespace

#endif
