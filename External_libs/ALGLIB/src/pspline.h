/*************************************************************************
Copyright (c) 2006-2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _pspline_h
#define _pspline_h

#include "ap.h"
#include "amp.h"
#include "spline3.h"
#include "blas.h"
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
#include "matinv.h"
#include "hblas.h"
#include "sblas.h"
#include "ortfac.h"
#include "rotations.h"
#include "bdsvd.h"
#include "svd.h"
#include "xblas.h"
#include "densesolver.h"
#include "linmin.h"
#include "minlbfgs.h"
#include "minlm.h"
#include "lsfit.h"
#include "apserv.h"
#include "spline1d.h"
#include "tsort.h"
#include "hsschur.h"
#include "evd.h"
#include "gammafunc.h"
#include "gq.h"
#include "gkq.h"
#include "autogk.h"
namespace pspline
{
    /*************************************************************************
    Parametric spline inteprolant: 2-dimensional curve.

    You should not try to access its members directly - use PSpline2XXXXXXXX()
    functions instead.
    *************************************************************************/
    template<unsigned int Precision>
    class pspline2interpolant
    {
    public:
        int n;
        bool periodic;
        ap::template_1d_array< amp::ampf<Precision> > p;
        spline1d::spline1dinterpolant<Precision> x;
        spline1d::spline1dinterpolant<Precision> y;
    };


    /*************************************************************************
    Parametric spline inteprolant: 3-dimensional curve.

    You should not try to access its members directly - use PSpline3XXXXXXXX()
    functions instead.
    *************************************************************************/
    template<unsigned int Precision>
    class pspline3interpolant
    {
    public:
        int n;
        bool periodic;
        ap::template_1d_array< amp::ampf<Precision> > p;
        spline1d::spline1dinterpolant<Precision> x;
        spline1d::spline1dinterpolant<Precision> y;
        spline1d::spline1dinterpolant<Precision> z;
    };




    template<unsigned int Precision>
    void pspline2build(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline2interpolant<Precision>& p);
    template<unsigned int Precision>
    void pspline3build(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline3interpolant<Precision>& p);
    template<unsigned int Precision>
    void pspline2buildperiodic(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline2interpolant<Precision>& p);
    template<unsigned int Precision>
    void pspline3buildperiodic(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline3interpolant<Precision>& p);
    template<unsigned int Precision>
    void pspline2parametervalues(const pspline2interpolant<Precision>& p,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& t);
    template<unsigned int Precision>
    void pspline3parametervalues(const pspline3interpolant<Precision>& p,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& t);
    template<unsigned int Precision>
    void pspline2calc(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y);
    template<unsigned int Precision>
    void pspline3calc(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& z);
    template<unsigned int Precision>
    void pspline2tangent(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y);
    template<unsigned int Precision>
    void pspline3tangent(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& z);
    template<unsigned int Precision>
    void pspline2diff(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy);
    template<unsigned int Precision>
    void pspline3diff(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& z,
        amp::ampf<Precision>& dz);
    template<unsigned int Precision>
    void pspline2diff2(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& d2x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& d2y);
    template<unsigned int Precision>
    void pspline3diff2(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& d2x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& d2y,
        amp::ampf<Precision>& z,
        amp::ampf<Precision>& dz,
        amp::ampf<Precision>& d2z);
    template<unsigned int Precision>
    amp::ampf<Precision> pspline2arclength(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    amp::ampf<Precision> pspline3arclength(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b);
    template<unsigned int Precision>
    void pspline2par(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int pt,
        ap::template_1d_array< amp::ampf<Precision> >& p);
    template<unsigned int Precision>
    void pspline3par(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int pt,
        ap::template_1d_array< amp::ampf<Precision> >& p);


    /*************************************************************************
    This function  builds  non-periodic 2-dimensional parametric spline  which
    starts at (X[0],Y[0]) and ends at (X[N-1],Y[N-1]).

    INPUT PARAMETERS:
        XY  -   points, array[0..N-1,0..1].
                XY[I,0:1] corresponds to the Ith point.
                Order of points is important!
        N   -   points count, N>=5 for Akima splines, N>=2 for other types  of
                splines.
        ST  -   spline type:
                * 0     Akima spline
                * 1     parabolically terminated Catmull-Rom spline (Tension=0)
                * 2     parabolically terminated cubic spline
        PT  -   parameterization type:
                * 0     uniform
                * 1     chord length
                * 2     centripetal

    OUTPUT PARAMETERS:
        P   -   parametric spline interpolant


    NOTES:
    * this function  assumes  that  there all consequent points  are distinct.
      I.e. (x0,y0)<>(x1,y1),  (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so on.
      However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
      =(x2,y2).

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2build(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline2interpolant<Precision>& p)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(st>=0 && st<=2);
        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        if( st==0 )
        {
            ap::ap_error::make_assertion(n>=5);
        }
        else
        {
            ap::ap_error::make_assertion(n>=2);
        }
        
        //
        // Prepare
        //
        p.n = n;
        p.periodic = false;
        tmp.setlength(n);
        
        //
        // Build parameterization, check that all parameters are distinct
        //
        pspline2par<Precision>(xy, n, pt, p.p);
        ap::ap_error::make_assertion(apserv::apservaredistinct<Precision>(p.p, n));
        
        //
        // Build splines
        //
        if( st==0 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildakima<Precision>(p.p, tmp, n, p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildakima<Precision>(p.p, tmp, n, p.y);
        }
        if( st==1 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), p.y);
        }
        if( st==2 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), p.y);
        }
    }


    /*************************************************************************
    This function  builds  non-periodic 3-dimensional parametric spline  which
    starts at (X[0],Y[0],Z[0]) and ends at (X[N-1],Y[N-1],Z[N-1]).

    Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
    description here.

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3build(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline3interpolant<Precision>& p)
    {
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(st>=0 && st<=2);
        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        if( st==0 )
        {
            ap::ap_error::make_assertion(n>=5);
        }
        else
        {
            ap::ap_error::make_assertion(n>=2);
        }
        
        //
        // Prepare
        //
        p.n = n;
        p.periodic = false;
        tmp.setlength(n);
        
        //
        // Build parameterization, check that all parameters are distinct
        //
        pspline3par<Precision>(xy, n, pt, p.p);
        ap::ap_error::make_assertion(apserv::apservaredistinct<Precision>(p.p, n));
        
        //
        // Build splines
        //
        if( st==0 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildakima<Precision>(p.p, tmp, n, p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildakima<Precision>(p.p, tmp, n, p.y);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(2, 0, n-1));
            spline1d::spline1dbuildakima<Precision>(p.p, tmp, n, p.z);
        }
        if( st==1 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), p.y);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(2, 0, n-1));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), p.z);
        }
        if( st==2 )
        {
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(0, 0, n-1));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(1, 0, n-1));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), p.y);
            amp::vmove(tmp.getvector(0, n-1), xy.getcolumn(2, 0, n-1));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n, 0, amp::ampf<Precision>("0.0"), 0, amp::ampf<Precision>("0.0"), p.z);
        }
    }


    /*************************************************************************
    This  function  builds  periodic  2-dimensional  parametric  spline  which
    starts at (X[0],Y[0]), goes through all points to (X[N-1],Y[N-1]) and then
    back to (X[0],Y[0]).

    INPUT PARAMETERS:
        XY  -   points, array[0..N-1,0..1].
                XY[I,0:1] corresponds to the Ith point.
                XY[N-1,0:1] must be different from XY[0,0:1].
                Order of points is important!
        N   -   points count, N>=3 for other types of splines.
        ST  -   spline type:
                * 1     Catmull-Rom spline (Tension=0) with cyclic boundary conditions
                * 2     cubic spline with cyclic boundary conditions
        PT  -   parameterization type:
                * 0     uniform
                * 1     chord length
                * 2     centripetal

    OUTPUT PARAMETERS:
        P   -   parametric spline interpolant


    NOTES:
    * this function  assumes  that there all consequent points  are  distinct.
      I.e. (x0,y0)<>(x1,y1), (x1,y1)<>(x2,y2),  (x2,y2)<>(x3,y3)  and  so  on.
      However, non-consequent points may coincide, i.e. we can  have  (x0,y0)=
      =(x2,y2).
    * last point of sequence is NOT equal to the first  point.  You  shouldn't
      make curve "explicitly periodic" by making them equal.

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2buildperiodic(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline2interpolant<Precision>& p)
    {
        ap::template_2d_array< amp::ampf<Precision> > xyp;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(st>=1 && st<=2);
        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        ap::ap_error::make_assertion(n>=3);
        
        //
        // Prepare
        //
        p.n = n;
        p.periodic = true;
        tmp.setlength(n+1);
        xyp.setlength(n+1, 2);
        amp::vmove(xyp.getcolumn(0, 0, n-1), xy.getcolumn(0, 0, n-1));
        amp::vmove(xyp.getcolumn(1, 0, n-1), xy.getcolumn(1, 0, n-1));
        amp::vmove(xyp.getrow(n, 0, 1), xy.getrow(0, 0, 1));
        
        //
        // Build parameterization, check that all parameters are distinct
        //
        pspline2par<Precision>(xyp, n+1, pt, p.p);
        ap::ap_error::make_assertion(apserv::apservaredistinct<Precision>(p.p, n+1));
        
        //
        // Build splines
        //
        if( st==1 )
        {
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(0, 0, n));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(1, 0, n));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), p.y);
        }
        if( st==2 )
        {
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(0, 0, n));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(1, 0, n));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), p.y);
        }
    }


    /*************************************************************************
    This  function  builds  periodic  3-dimensional  parametric  spline  which
    starts at (X[0],Y[0],Z[0]), goes through all points to (X[N-1],Y[N-1],Z[N-1])
    and then back to (X[0],Y[0],Z[0]).

    Same as PSpline2Build() function, but for 3D, so we  won't  duplicate  its
    description here.

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3buildperiodic(ap::template_2d_array< amp::ampf<Precision> > xy,
        int n,
        int st,
        int pt,
        pspline3interpolant<Precision>& p)
    {
        ap::template_2d_array< amp::ampf<Precision> > xyp;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(st>=1 && st<=2);
        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        ap::ap_error::make_assertion(n>=3);
        
        //
        // Prepare
        //
        p.n = n;
        p.periodic = true;
        tmp.setlength(n+1);
        xyp.setlength(n+1, 3);
        amp::vmove(xyp.getcolumn(0, 0, n-1), xy.getcolumn(0, 0, n-1));
        amp::vmove(xyp.getcolumn(1, 0, n-1), xy.getcolumn(1, 0, n-1));
        amp::vmove(xyp.getcolumn(2, 0, n-1), xy.getcolumn(2, 0, n-1));
        amp::vmove(xyp.getrow(n, 0, 2), xy.getrow(0, 0, 2));
        
        //
        // Build parameterization, check that all parameters are distinct
        //
        pspline3par<Precision>(xyp, n+1, pt, p.p);
        ap::ap_error::make_assertion(apserv::apservaredistinct<Precision>(p.p, n+1));
        
        //
        // Build splines
        //
        if( st==1 )
        {
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(0, 0, n));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(1, 0, n));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), p.y);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(2, 0, n));
            spline1d::spline1dbuildcatmullrom<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), p.z);
        }
        if( st==2 )
        {
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(0, 0, n));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), p.x);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(1, 0, n));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), p.y);
            amp::vmove(tmp.getvector(0, n), xyp.getcolumn(2, 0, n));
            spline1d::spline1dbuildcubic<Precision>(p.p, tmp, n+1, -1, amp::ampf<Precision>("0.0"), -1, amp::ampf<Precision>("0.0"), p.z);
        }
    }


    /*************************************************************************
    This function returns vector of parameter values correspoding to points.

    I.e. for P created from (X[0],Y[0])...(X[N-1],Y[N-1]) and U=TValues(P)  we
    have
        (X[0],Y[0]) = PSpline2Calc(P,U[0]),
        (X[1],Y[1]) = PSpline2Calc(P,U[1]),
        (X[2],Y[2]) = PSpline2Calc(P,U[2]),
        ...

    INPUT PARAMETERS:
        P   -   parametric spline interpolant

    OUTPUT PARAMETERS:
        N   -   array size
        T   -   array[0..N-1]


    NOTES:
    * for non-periodic splines U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]=1
    * for periodic splines     U[0]=0, U[0]<U[1]<...<U[N-1], U[N-1]<1

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2parametervalues(const pspline2interpolant<Precision>& p,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& t)
    {
        ap::ap_error::make_assertion(p.n>=2);
        n = p.n;
        t.setlength(n);
        amp::vmove(t.getvector(0, n-1), p.p.getvector(0, n-1));
        t(0) = 0;
        if( !p.periodic )
        {
            t(n-1) = 1;
        }
    }


    /*************************************************************************
    This function returns vector of parameter values correspoding to points.

    Same as PSpline2ParameterValues(), but for 3D.

      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3parametervalues(const pspline3interpolant<Precision>& p,
        int& n,
        ap::template_1d_array< amp::ampf<Precision> >& t)
    {
        ap::ap_error::make_assertion(p.n>=2);
        n = p.n;
        t.setlength(n);
        amp::vmove(t.getvector(0, n-1), p.p.getvector(0, n-1));
        t(0) = 0;
        if( !p.periodic )
        {
            t(n-1) = 1;
        }
    }


    /*************************************************************************
    This function  calculates  the value of the parametric spline for a  given
    value of parameter T

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-position
        Y   -   Y-position


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2calc(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y)
    {
        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        x = spline1d::spline1dcalc<Precision>(p.x, t);
        y = spline1d::spline1dcalc<Precision>(p.y, t);
    }


    /*************************************************************************
    This function  calculates  the value of the parametric spline for a  given
    value of parameter T.

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-position
        Y   -   Y-position
        Z   -   Z-position


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3calc(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& z)
    {
        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        x = spline1d::spline1dcalc<Precision>(p.x, t);
        y = spline1d::spline1dcalc<Precision>(p.y, t);
        z = spline1d::spline1dcalc<Precision>(p.z, t);
    }


    /*************************************************************************
    This function  calculates  tangent vector for a given value of parameter T

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X    -   X-component of tangent vector (normalized)
        Y    -   Y-component of tangent vector (normalized)
        
    NOTE:
        X^2+Y^2 is either 1 (for non-zero tangent vector) or 0.


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2tangent(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;


        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        pspline2diff<Precision>(p, t, v0, x, v1, y);
        if( x!=0 || y!=0 )
        {
            
            //
            // this code is a bit more complex than X^2+Y^2 to avoid
            // overflow for large values of X and Y.
            //
            v = apserv::safepythag2<Precision>(x, y);
            x = x/v;
            y = y/v;
        }
    }


    /*************************************************************************
    This function  calculates  tangent vector for a given value of parameter T

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X    -   X-component of tangent vector (normalized)
        Y    -   Y-component of tangent vector (normalized)
        Z    -   Z-component of tangent vector (normalized)

    NOTE:
        X^2+Y^2+Z^2 is either 1 (for non-zero tangent vector) or 0.


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3tangent(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& z)
    {
        amp::ampf<Precision> v;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;


        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        pspline3diff<Precision>(p, t, v0, x, v1, y, v2, z);
        if( x!=0 || y!=0 || z!=0 )
        {
            v = apserv::safepythag3<Precision>(x, y, z);
            x = x/v;
            y = y/v;
            z = z/v;
        }
    }


    /*************************************************************************
    This function calculates derivative, i.e. it returns (dX/dT,dY/dT).

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-value
        DX  -   X-derivative
        Y   -   Y-value
        DY  -   Y-derivative


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2diff(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy)
    {
        amp::ampf<Precision> d2s;


        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        spline1d::spline1ddiff<Precision>(p.x, t, x, dx, d2s);
        spline1d::spline1ddiff<Precision>(p.y, t, y, dy, d2s);
    }


    /*************************************************************************
    This function calculates derivative, i.e. it returns (dX/dT,dY/dT,dZ/dT).

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-value
        DX  -   X-derivative
        Y   -   Y-value
        DY  -   Y-derivative
        Z   -   Z-value
        DZ  -   Z-derivative


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3diff(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& z,
        amp::ampf<Precision>& dz)
    {
        amp::ampf<Precision> d2s;


        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        spline1d::spline1ddiff<Precision>(p.x, t, x, dx, d2s);
        spline1d::spline1ddiff<Precision>(p.y, t, y, dy, d2s);
        spline1d::spline1ddiff<Precision>(p.z, t, z, dz, d2s);
    }


    /*************************************************************************
    This function calculates first and second derivative with respect to T.

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-value
        DX  -   derivative
        D2X -   second derivative
        Y   -   Y-value
        DY  -   derivative
        D2Y -   second derivative


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2diff2(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& d2x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& d2y)
    {
        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        spline1d::spline1ddiff<Precision>(p.x, t, x, dx, d2x);
        spline1d::spline1ddiff<Precision>(p.y, t, y, dy, d2y);
    }


    /*************************************************************************
    This function calculates first and second derivative with respect to T.

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        T   -   point:
                * T in [0,1] corresponds to interval spanned by points
                * for non-periodic splines T<0 (or T>1) correspond to parts of
                  the curve before the first (after the last) point
                * for periodic splines T<0 (or T>1) are projected  into  [0,1]
                  by making T=T-floor(T).

    OUTPUT PARAMETERS:
        X   -   X-value
        DX  -   derivative
        D2X -   second derivative
        Y   -   Y-value
        DY  -   derivative
        D2Y -   second derivative
        Z   -   Z-value
        DZ  -   derivative
        D2Z -   second derivative


      -- ALGLIB PROJECT --
         Copyright 28.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3diff2(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> t,
        amp::ampf<Precision>& x,
        amp::ampf<Precision>& dx,
        amp::ampf<Precision>& d2x,
        amp::ampf<Precision>& y,
        amp::ampf<Precision>& dy,
        amp::ampf<Precision>& d2y,
        amp::ampf<Precision>& z,
        amp::ampf<Precision>& dz,
        amp::ampf<Precision>& d2z)
    {
        if( p.periodic )
        {
            t = t-amp::floor<Precision>(t);
        }
        spline1d::spline1ddiff<Precision>(p.x, t, x, dx, d2x);
        spline1d::spline1ddiff<Precision>(p.y, t, y, dy, d2y);
        spline1d::spline1ddiff<Precision>(p.z, t, z, dz, d2z);
    }


    /*************************************************************************
    This function  calculates  arc length, i.e. length of  curve  between  t=a
    and t=b.

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        A,B -   parameter values corresponding to arc ends:
                * B>A will result in positive length returned
                * B<A will result in negative length returned

    RESULT:
        length of arc starting at T=A and ending at T=B.


      -- ALGLIB PROJECT --
         Copyright 30.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> pspline2arclength(const pspline2interpolant<Precision>& p,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;
        autogk::autogkstate<Precision> state;
        autogk::autogkreport<Precision> rep;
        amp::ampf<Precision> sx;
        amp::ampf<Precision> dsx;
        amp::ampf<Precision> d2sx;
        amp::ampf<Precision> sy;
        amp::ampf<Precision> dsy;
        amp::ampf<Precision> d2sy;


        autogk::autogksmooth<Precision>(a, b, state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            spline1d::spline1ddiff<Precision>(p.x, state.x, sx, dsx, d2sx);
            spline1d::spline1ddiff<Precision>(p.y, state.x, sy, dsy, d2sy);
            state.f = apserv::safepythag2<Precision>(dsx, dsy);
        }
        autogk::autogkresults<Precision>(state, result, rep);
        ap::ap_error::make_assertion(rep.terminationtype>0);
        return result;
    }


    /*************************************************************************
    This function  calculates  arc length, i.e. length of  curve  between  t=a
    and t=b.

    INPUT PARAMETERS:
        P   -   parametric spline interpolant
        A,B -   parameter values corresponding to arc ends:
                * B>A will result in positive length returned
                * B<A will result in negative length returned

    RESULT:
        length of arc starting at T=A and ending at T=B.


      -- ALGLIB PROJECT --
         Copyright 30.05.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> pspline3arclength(const pspline3interpolant<Precision>& p,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;
        autogk::autogkstate<Precision> state;
        autogk::autogkreport<Precision> rep;
        amp::ampf<Precision> sx;
        amp::ampf<Precision> dsx;
        amp::ampf<Precision> d2sx;
        amp::ampf<Precision> sy;
        amp::ampf<Precision> dsy;
        amp::ampf<Precision> d2sy;
        amp::ampf<Precision> sz;
        amp::ampf<Precision> dsz;
        amp::ampf<Precision> d2sz;


        autogk::autogksmooth<Precision>(a, b, state);
        while( autogk::autogkiteration<Precision>(state) )
        {
            spline1d::spline1ddiff<Precision>(p.x, state.x, sx, dsx, d2sx);
            spline1d::spline1ddiff<Precision>(p.y, state.x, sy, dsy, d2sy);
            spline1d::spline1ddiff<Precision>(p.z, state.x, sz, dsz, d2sz);
            state.f = apserv::safepythag3<Precision>(dsx, dsy, dsz);
        }
        autogk::autogkresults<Precision>(state, result, rep);
        ap::ap_error::make_assertion(rep.terminationtype>0);
        return result;
    }


    /*************************************************************************
    Builds non-periodic parameterization for 2-dimensional spline
    *************************************************************************/
    template<unsigned int Precision>
    void pspline2par(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int pt,
        ap::template_1d_array< amp::ampf<Precision> >& p)
    {
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        
        //
        // Build parameterization:
        // * fill by non-normalized values
        // * normalize them so we have P[0]=0, P[N-1]=1.
        //
        p.setlength(n);
        if( pt==0 )
        {
            for(i=0; i<=n-1; i++)
            {
                p(i) = i;
            }
        }
        if( pt==1 )
        {
            p(0) = 0;
            for(i=1; i<=n-1; i++)
            {
                p(i) = p(i-1)+apserv::safepythag2<Precision>(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1));
            }
        }
        if( pt==2 )
        {
            p(0) = 0;
            for(i=1; i<=n-1; i++)
            {
                p(i) = p(i-1)+amp::sqrt<Precision>(apserv::safepythag2<Precision>(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1)));
            }
        }
        v = 1/p(n-1);
        amp::vmul(p.getvector(0, n-1), v);
    }


    /*************************************************************************
    Builds non-periodic parameterization for 3-dimensional spline
    *************************************************************************/
    template<unsigned int Precision>
    void pspline3par(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int pt,
        ap::template_1d_array< amp::ampf<Precision> >& p)
    {
        amp::ampf<Precision> v;
        int i;


        ap::ap_error::make_assertion(pt>=0 && pt<=2);
        
        //
        // Build parameterization:
        // * fill by non-normalized values
        // * normalize them so we have P[0]=0, P[N-1]=1.
        //
        p.setlength(n);
        if( pt==0 )
        {
            for(i=0; i<=n-1; i++)
            {
                p(i) = i;
            }
        }
        if( pt==1 )
        {
            p(0) = 0;
            for(i=1; i<=n-1; i++)
            {
                p(i) = p(i-1)+apserv::safepythag3<Precision>(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1), xy(i,2)-xy(i-1,2));
            }
        }
        if( pt==2 )
        {
            p(0) = 0;
            for(i=1; i<=n-1; i++)
            {
                p(i) = p(i-1)+amp::sqrt<Precision>(apserv::safepythag3<Precision>(xy(i,0)-xy(i-1,0), xy(i,1)-xy(i-1,1), xy(i,2)-xy(i-1,2)));
            }
        }
        v = 1/p(n-1);
        amp::vmul(p.getvector(0, n-1), v);
    }
} // namespace

#endif
