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

#ifndef _inverseupdate_h
#define _inverseupdate_h

#include "ap.h"
#include "amp.h"
namespace inverseupdate
{
    template<unsigned int Precision>
    void rmatrixinvupdatesimple(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        int updcolumn,
        amp::ampf<Precision> updval);
    template<unsigned int Precision>
    void rmatrixinvupdaterow(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        const ap::template_1d_array< amp::ampf<Precision> >& v);
    template<unsigned int Precision>
    void rmatrixinvupdatecolumn(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updcolumn,
        const ap::template_1d_array< amp::ampf<Precision> >& u);
    template<unsigned int Precision>
    void rmatrixinvupdateuv(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v);
    template<unsigned int Precision>
    void shermanmorrisonsimpleupdate(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        int updcolumn,
        amp::ampf<Precision> updval);
    template<unsigned int Precision>
    void shermanmorrisonupdaterow(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        const ap::template_1d_array< amp::ampf<Precision> >& v);
    template<unsigned int Precision>
    void shermanmorrisonupdatecolumn(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updcolumn,
        const ap::template_1d_array< amp::ampf<Precision> >& u);
    template<unsigned int Precision>
    void shermanmorrisonupdateuv(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v);


    /*************************************************************************
    Inverse matrix update by the Sherman-Morrison formula

    The algorithm updates matrix A^-1 when adding a number to an element
    of matrix A.

    Input parameters:
        InvA    -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        UpdRow  -   row where the element to be updated is stored.
        UpdColumn - column where the element to be updated is stored.
        UpdVal  -   a number to be added to the element.


    Output parameters:
        InvA    -   inverse of modified matrix A.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixinvupdatesimple(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        int updcolumn,
        amp::ampf<Precision> updval)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        ap::ap_error::make_assertion(updrow>=0 && updrow<n);
        ap::ap_error::make_assertion(updcolumn>=0 && updcolumn<n);
        t1.setbounds(0, n-1);
        t2.setbounds(0, n-1);
        
        //
        // T1 = InvA * U
        //
        amp::vmove(t1.getvector(0, n-1), inva.getcolumn(updrow, 0, n-1));
        
        //
        // T2 = v*InvA
        //
        amp::vmove(t2.getvector(0, n-1), inva.getrow(updcolumn, 0, n-1));
        
        //
        // Lambda = v * InvA * U
        //
        lambda = updval*inva(updcolumn,updrow);
        
        //
        // InvA = InvA - correction
        //
        for(i=0; i<=n-1; i++)
        {
            vt = updval*t1(i);
            vt = vt/(1+lambda);
            amp::vsub(inva.getrow(i, 0, n-1), t2.getvector(0, n-1), vt);
        }
    }


    /*************************************************************************
    Inverse matrix update by the Sherman-Morrison formula

    The algorithm updates matrix A^-1 when adding a vector to a row
    of matrix A.

    Input parameters:
        InvA    -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        UpdRow  -   the row of A whose vector V was added.
                    0 <= Row <= N-1
        V       -   the vector to be added to a row.
                    Array whose index ranges within [0..N-1].

    Output parameters:
        InvA    -   inverse of modified matrix A.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixinvupdaterow(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        const ap::template_1d_array< amp::ampf<Precision> >& v)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        int j;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(0, n-1);
        t2.setbounds(0, n-1);
        
        //
        // T1 = InvA * U
        //
        amp::vmove(t1.getvector(0, n-1), inva.getcolumn(updrow, 0, n-1));
        
        //
        // T2 = v*InvA
        // Lambda = v * InvA * U
        //
        for(j=0; j<=n-1; j++)
        {
            vt = amp::vdotproduct(v.getvector(0, n-1), inva.getcolumn(j, 0, n-1));
            t2(j) = vt;
        }
        lambda = t2(updrow);
        
        //
        // InvA = InvA - correction
        //
        for(i=0; i<=n-1; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 0, n-1), t2.getvector(0, n-1), vt);
        }
    }


    /*************************************************************************
    Inverse matrix update by the Sherman-Morrison formula

    The algorithm updates matrix A^-1 when adding a vector to a column
    of matrix A.

    Input parameters:
        InvA        -   inverse of matrix A.
                        Array whose indexes range within [0..N-1, 0..N-1].
        N           -   size of matrix A.
        UpdColumn   -   the column of A whose vector U was added.
                        0 <= UpdColumn <= N-1
        U           -   the vector to be added to a column.
                        Array whose index ranges within [0..N-1].

    Output parameters:
        InvA        -   inverse of modified matrix A.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixinvupdatecolumn(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updcolumn,
        const ap::template_1d_array< amp::ampf<Precision> >& u)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(0, n-1);
        t2.setbounds(0, n-1);
        
        //
        // T1 = InvA * U
        // Lambda = v * InvA * U
        //
        for(i=0; i<=n-1; i++)
        {
            vt = amp::vdotproduct(inva.getrow(i, 0, n-1), u.getvector(0, n-1));
            t1(i) = vt;
        }
        lambda = t1(updcolumn);
        
        //
        // T2 = v*InvA
        //
        amp::vmove(t2.getvector(0, n-1), inva.getrow(updcolumn, 0, n-1));
        
        //
        // InvA = InvA - correction
        //
        for(i=0; i<=n-1; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 0, n-1), t2.getvector(0, n-1), vt);
        }
    }


    /*************************************************************************
    Inverse matrix update by the Sherman-Morrison formula

    The algorithm computes the inverse of matrix A+u*v’ by using the given matrix
    A^-1 and the vectors u and v.

    Input parameters:
        InvA    -   inverse of matrix A.
                    Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        U       -   the vector modifying the matrix.
                    Array whose index ranges within [0..N-1].
        V       -   the vector modifying the matrix.
                    Array whose index ranges within [0..N-1].

    Output parameters:
        InvA - inverse of matrix A + u*v'.

      -- ALGLIB --
         Copyright 2005 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixinvupdateuv(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        int j;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(0, n-1);
        t2.setbounds(0, n-1);
        
        //
        // T1 = InvA * U
        // Lambda = v * T1
        //
        for(i=0; i<=n-1; i++)
        {
            vt = amp::vdotproduct(inva.getrow(i, 0, n-1), u.getvector(0, n-1));
            t1(i) = vt;
        }
        lambda = amp::vdotproduct(v.getvector(0, n-1), t1.getvector(0, n-1));
        
        //
        // T2 = v*InvA
        //
        for(j=0; j<=n-1; j++)
        {
            vt = amp::vdotproduct(v.getvector(0, n-1), inva.getcolumn(j, 0, n-1));
            t2(j) = vt;
        }
        
        //
        // InvA = InvA - correction
        //
        for(i=0; i<=n-1; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 0, n-1), t2.getvector(0, n-1), vt);
        }
    }


    template<unsigned int Precision>
    void shermanmorrisonsimpleupdate(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        int updcolumn,
        amp::ampf<Precision> updval)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(1, n);
        t2.setbounds(1, n);
        
        //
        // T1 = InvA * U
        //
        amp::vmove(t1.getvector(1, n), inva.getcolumn(updrow, 1, n));
        
        //
        // T2 = v*InvA
        //
        amp::vmove(t2.getvector(1, n), inva.getrow(updcolumn, 1, n));
        
        //
        // Lambda = v * InvA * U
        //
        lambda = updval*inva(updcolumn,updrow);
        
        //
        // InvA = InvA - correction
        //
        for(i=1; i<=n; i++)
        {
            vt = updval*t1(i);
            vt = vt/(1+lambda);
            amp::vsub(inva.getrow(i, 1, n), t2.getvector(1, n), vt);
        }
    }


    template<unsigned int Precision>
    void shermanmorrisonupdaterow(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updrow,
        const ap::template_1d_array< amp::ampf<Precision> >& v)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        int j;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(1, n);
        t2.setbounds(1, n);
        
        //
        // T1 = InvA * U
        //
        amp::vmove(t1.getvector(1, n), inva.getcolumn(updrow, 1, n));
        
        //
        // T2 = v*InvA
        // Lambda = v * InvA * U
        //
        for(j=1; j<=n; j++)
        {
            vt = amp::vdotproduct(v.getvector(1, n), inva.getcolumn(j, 1, n));
            t2(j) = vt;
        }
        lambda = t2(updrow);
        
        //
        // InvA = InvA - correction
        //
        for(i=1; i<=n; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 1, n), t2.getvector(1, n), vt);
        }
    }


    template<unsigned int Precision>
    void shermanmorrisonupdatecolumn(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        int updcolumn,
        const ap::template_1d_array< amp::ampf<Precision> >& u)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(1, n);
        t2.setbounds(1, n);
        
        //
        // T1 = InvA * U
        // Lambda = v * InvA * U
        //
        for(i=1; i<=n; i++)
        {
            vt = amp::vdotproduct(inva.getrow(i, 1, n), u.getvector(1, n));
            t1(i) = vt;
        }
        lambda = t1(updcolumn);
        
        //
        // T2 = v*InvA
        //
        amp::vmove(t2.getvector(1, n), inva.getrow(updcolumn, 1, n));
        
        //
        // InvA = InvA - correction
        //
        for(i=1; i<=n; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 1, n), t2.getvector(1, n), vt);
        }
    }


    template<unsigned int Precision>
    void shermanmorrisonupdateuv(ap::template_2d_array< amp::ampf<Precision> >& inva,
        int n,
        const ap::template_1d_array< amp::ampf<Precision> >& u,
        const ap::template_1d_array< amp::ampf<Precision> >& v)
    {
        ap::template_1d_array< amp::ampf<Precision> > t1;
        ap::template_1d_array< amp::ampf<Precision> > t2;
        int i;
        int j;
        amp::ampf<Precision> lambda;
        amp::ampf<Precision> vt;


        t1.setbounds(1, n);
        t2.setbounds(1, n);
        
        //
        // T1 = InvA * U
        // Lambda = v * T1
        //
        for(i=1; i<=n; i++)
        {
            vt = amp::vdotproduct(inva.getrow(i, 1, n), u.getvector(1, n));
            t1(i) = vt;
        }
        lambda = amp::vdotproduct(v.getvector(1, n), t1.getvector(1, n));
        
        //
        // T2 = v*InvA
        //
        for(j=1; j<=n; j++)
        {
            vt = amp::vdotproduct(v.getvector(1, n), inva.getcolumn(j, 1, n));
            t2(j) = vt;
        }
        
        //
        // InvA = InvA - correction
        //
        for(i=1; i<=n; i++)
        {
            vt = t1(i)/(1+lambda);
            amp::vsub(inva.getrow(i, 1, n), t2.getvector(1, n), vt);
        }
    }
} // namespace

#endif
