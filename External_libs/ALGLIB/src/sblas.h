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

#ifndef _sblas_h
#define _sblas_h

#include "ap.h"
#include "amp.h"
namespace sblas
{
    template<unsigned int Precision>
    void symmetricmatrixvectormultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool isupper,
        int i1,
        int i2,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision> alpha,
        ap::template_1d_array< amp::ampf<Precision> >& y);
    template<unsigned int Precision>
    void symmetricrank2update(ap::template_2d_array< amp::ampf<Precision> >& a,
        bool isupper,
        int i1,
        int i2,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& t,
        amp::ampf<Precision> alpha);


    template<unsigned int Precision>
    void symmetricmatrixvectormultiply(const ap::template_2d_array< amp::ampf<Precision> >& a,
        bool isupper,
        int i1,
        int i2,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision> alpha,
        ap::template_1d_array< amp::ampf<Precision> >& y)
    {
        int i;
        int ba1;
        int ba2;
        int by1;
        int by2;
        int bx1;
        int bx2;
        int n;
        amp::ampf<Precision> v;


        n = i2-i1+1;
        if( n<=0 )
        {
            return;
        }
        
        //
        // Let A = L + D + U, where
        //  L is strictly lower triangular (main diagonal is zero)
        //  D is diagonal
        //  U is strictly upper triangular (main diagonal is zero)
        //
        // A*x = L*x + D*x + U*x
        //
        // Calculate D*x first
        //
        for(i=i1; i<=i2; i++)
        {
            y(i-i1+1) = a(i,i)*x(i-i1+1);
        }
        
        //
        // Add L*x + U*x
        //
        if( isupper )
        {
            for(i=i1; i<=i2-1; i++)
            {
                
                //
                // Add L*x to the result
                //
                v = x(i-i1+1);
                by1 = i-i1+2;
                by2 = n;
                ba1 = i+1;
                ba2 = i2;
                amp::vadd(y.getvector(by1, by2), a.getrow(i, ba1, ba2), v);
                
                //
                // Add U*x to the result
                //
                bx1 = i-i1+2;
                bx2 = n;
                ba1 = i+1;
                ba2 = i2;
                v = amp::vdotproduct(x.getvector(bx1, bx2), a.getrow(i, ba1, ba2));
                y(i-i1+1) = y(i-i1+1)+v;
            }
        }
        else
        {
            for(i=i1+1; i<=i2; i++)
            {
                
                //
                // Add L*x to the result
                //
                bx1 = 1;
                bx2 = i-i1;
                ba1 = i1;
                ba2 = i-1;
                v = amp::vdotproduct(x.getvector(bx1, bx2), a.getrow(i, ba1, ba2));
                y(i-i1+1) = y(i-i1+1)+v;
                
                //
                // Add U*x to the result
                //
                v = x(i-i1+1);
                by1 = 1;
                by2 = i-i1;
                ba1 = i1;
                ba2 = i-1;
                amp::vadd(y.getvector(by1, by2), a.getrow(i, ba1, ba2), v);
            }
        }
        amp::vmul(y.getvector(1, n), alpha);
    }


    template<unsigned int Precision>
    void symmetricrank2update(ap::template_2d_array< amp::ampf<Precision> >& a,
        bool isupper,
        int i1,
        int i2,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        ap::template_1d_array< amp::ampf<Precision> >& t,
        amp::ampf<Precision> alpha)
    {
        int i;
        int tp1;
        int tp2;
        amp::ampf<Precision> v;


        if( isupper )
        {
            for(i=i1; i<=i2; i++)
            {
                tp1 = i+1-i1;
                tp2 = i2-i1+1;
                v = x(i+1-i1);
                amp::vmove(t.getvector(tp1, tp2), y.getvector(tp1, tp2), v);
                v = y(i+1-i1);
                amp::vadd(t.getvector(tp1, tp2), x.getvector(tp1, tp2), v);
                amp::vmul(t.getvector(tp1, tp2), alpha);
                amp::vadd(a.getrow(i, i, i2), t.getvector(tp1, tp2));
            }
        }
        else
        {
            for(i=i1; i<=i2; i++)
            {
                tp1 = 1;
                tp2 = i+1-i1;
                v = x(i+1-i1);
                amp::vmove(t.getvector(tp1, tp2), y.getvector(tp1, tp2), v);
                v = y(i+1-i1);
                amp::vadd(t.getvector(tp1, tp2), x.getvector(tp1, tp2), v);
                amp::vmul(t.getvector(tp1, tp2), alpha);
                amp::vadd(a.getrow(i, i1, i), t.getvector(tp1, tp2));
            }
        }
    }
} // namespace

#endif
