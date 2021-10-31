/*************************************************************************
Copyright (c) 2007, Sergey Bochkanov (ALGLIB project).

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

#ifndef _ratinterpolation_h
#define _ratinterpolation_h

#include "ap.h"
#include "amp.h"
#include "tsort.h"
namespace ratinterpolation
{
    template<unsigned int Precision>
    void buildfloaterhormannrationalinterpolant(ap::template_1d_array< amp::ampf<Precision> > x,
        int n,
        int d,
        ap::template_1d_array< amp::ampf<Precision> >& w);


    /*************************************************************************
    Rational barycentric interpolation without poles

    The subroutine constructs the rational interpolating function without real
    poles. It should be noted that the barycentric weights of the  interpolant
    constructed are independent of the values of the given function.

    Input parameters:
        X   -   interpolation nodes, array[0..N-1].
        N   -   number of nodes, N>0.
        D   -   order of the interpolation scheme, 0 <= D <= N-1.

    Output parameters:
        W   -   array of the barycentric weights which  can  be  used  in  the
                BarycentricInterpolate subroutine. Array[0..N-1]

    Note:
        this algorithm always succeeds and calculates the weights  with  close
        to machine precision.

      -- ALGLIB PROJECT --
         Copyright 17.06.2007 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void buildfloaterhormannrationalinterpolant(ap::template_1d_array< amp::ampf<Precision> > x,
        int n,
        int d,
        ap::template_1d_array< amp::ampf<Precision> >& w)
    {
        amp::ampf<Precision> s0;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;
        int i;
        int j;
        int k;
        ap::template_1d_array< int > perm;
        ap::template_1d_array< amp::ampf<Precision> > wtemp;


        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(d>=0 && d<=n);
        
        //
        // Prepare
        //
        w.setbounds(0, n-1);
        s0 = 1;
        for(k=1; k<=d; k++)
        {
            s0 = -s0;
        }
        perm.setbounds(0, n-1);
        for(i=0; i<=n-1; i++)
        {
            perm(i) = i;
        }
        tsort::tagsortfasti<Precision>(x, perm, n);
        
        //
        // Calculate Wk
        //
        for(k=0; k<=n-1; k++)
        {
            
            //
            // Wk
            //
            s = 0;
            for(i=ap::maxint(k-d, 0); i<=ap::minint(k, n-1-d); i++)
            {
                v = 1;
                for(j=i; j<=i+d; j++)
                {
                    if( j!=k )
                    {
                        v = v/amp::abs<Precision>(x(k)-x(j));
                    }
                }
                s = s+v;
            }
            w(k) = s0*s;
            
            //
            // Next S0
            //
            s0 = -s0;
        }
        
        //
        // Reorder W
        //
        wtemp.setbounds(0, n-1);
        amp::vmove(wtemp.getvector(0, n-1), w.getvector(0, n-1));
        for(i=0; i<=n-1; i++)
        {
            w(perm(i)) = wtemp(i);
        }
    }
} // namespace

#endif
