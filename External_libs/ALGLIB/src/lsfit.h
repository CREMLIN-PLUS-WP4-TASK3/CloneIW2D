/*************************************************************************
Copyright (c) 2006-2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _lsfit_h
#define _lsfit_h

#include "ap.h"
#include "amp.h"
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
namespace lsfit
{
    /*************************************************************************
    Least squares fitting report:
        TaskRCond       reciprocal of task's condition number
        RMSError        RMS error
        AvgError        average error
        AvgRelError     average relative error (for non-zero Y[I])
        MaxError        maximum error
    *************************************************************************/
    template<unsigned int Precision>
    class lsfitreport
    {
    public:
        amp::ampf<Precision> taskrcond;
        amp::ampf<Precision> rmserror;
        amp::ampf<Precision> avgerror;
        amp::ampf<Precision> avgrelerror;
        amp::ampf<Precision> maxerror;
    };


    template<unsigned int Precision>
    class lsfitstate
    {
    public:
        int n;
        int m;
        int k;
        amp::ampf<Precision> epsf;
        amp::ampf<Precision> epsx;
        int maxits;
        amp::ampf<Precision> stpmax;
        ap::template_2d_array< amp::ampf<Precision> > taskx;
        ap::template_1d_array< amp::ampf<Precision> > tasky;
        ap::template_1d_array< amp::ampf<Precision> > w;
        bool cheapfg;
        bool havehess;
        bool needf;
        bool needfg;
        bool needfgh;
        int pointindex;
        ap::template_1d_array< amp::ampf<Precision> > x;
        ap::template_1d_array< amp::ampf<Precision> > c;
        amp::ampf<Precision> f;
        ap::template_1d_array< amp::ampf<Precision> > g;
        ap::template_2d_array< amp::ampf<Precision> > h;
        int repterminationtype;
        amp::ampf<Precision> reprmserror;
        amp::ampf<Precision> repavgerror;
        amp::ampf<Precision> repavgrelerror;
        amp::ampf<Precision> repmaxerror;
        minlm::minlmstate<Precision> optstate;
        minlm::minlmreport<Precision> optrep;
        amp::rcommstate<Precision> rstate;
    };




    template<unsigned int Precision>
    void lsfitlinearw(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitlinearwc(ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        ap::template_2d_array< amp::ampf<Precision> > cmatrix,
        int n,
        int m,
        int k,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitlinear(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitlinearc(ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_2d_array< amp::ampf<Precision> >& cmatrix,
        int n,
        int m,
        int k,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitnonlinearwfg(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        bool cheapfg,
        lsfitstate<Precision>& state);
    template<unsigned int Precision>
    void lsfitnonlinearfg(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        bool cheapfg,
        lsfitstate<Precision>& state);
    template<unsigned int Precision>
    void lsfitnonlinearwfgh(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        lsfitstate<Precision>& state);
    template<unsigned int Precision>
    void lsfitnonlinearfgh(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        lsfitstate<Precision>& state);
    template<unsigned int Precision>
    void lsfitnonlinearsetcond(lsfitstate<Precision>& state,
        amp::ampf<Precision> epsf,
        amp::ampf<Precision> epsx,
        int maxits);
    template<unsigned int Precision>
    void lsfitnonlinearsetstpmax(lsfitstate<Precision>& state,
        amp::ampf<Precision> stpmax);
    template<unsigned int Precision>
    bool lsfitnonlineariteration(lsfitstate<Precision>& state);
    template<unsigned int Precision>
    void lsfitnonlinearresults(const lsfitstate<Precision>& state,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitscalexy(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& xc,
        ap::template_1d_array< amp::ampf<Precision> >& yc,
        const ap::template_1d_array< int >& dc,
        int k,
        amp::ampf<Precision>& xa,
        amp::ampf<Precision>& xb,
        amp::ampf<Precision>& sa,
        amp::ampf<Precision>& sb,
        ap::template_1d_array< amp::ampf<Precision> >& xoriginal,
        ap::template_1d_array< amp::ampf<Precision> >& yoriginal);
    template<unsigned int Precision>
    void lsfitlinearinternal(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep);
    template<unsigned int Precision>
    void lsfitclearrequestfields(lsfitstate<Precision>& state);


    /*************************************************************************
    Weighted linear least squares fitting.

    QR decomposition is used to reduce task to MxM, then triangular solver  or
    SVD-based solver is used depending on condition number of the  system.  It
    allows to maximize speed and retain decent accuracy.

    INPUT PARAMETERS:
        Y       -   array[0..N-1] Function values in  N  points.
        W       -   array[0..N-1]  Weights  corresponding to function  values.
                    Each summand in square  sum  of  approximation  deviations
                    from  given  values  is  multiplied  by  the   square   of
                    corresponding weight.
        FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                    FMatrix[I, J] - value of J-th basis function in I-th point.
        N       -   number of points used. N>=1.
        M       -   number of basis functions, M>=1.

    OUTPUT PARAMETERS:
        Info    -   error code:
                    * -4    internal SVD decomposition subroutine failed (very
                            rare and for degenerate systems only)
                    * -1    incorrect N/M were specified
                    *  1    task is solved
        C       -   decomposition coefficients, array[0..M-1]
        Rep     -   fitting report. Following fields are set:
                    * Rep.TaskRCond     reciprocal of condition number
                    * RMSError          rms error on the (X,Y).
                    * AvgError          average error on the (X,Y).
                    * AvgRelError       average relative error on the non-zero Y
                    * MaxError          maximum error
                                        NON-WEIGHTED ERRORS ARE CALCULATED

    SEE ALSO
        LSFitLinear
        LSFitLinearC
        LSFitLinearWC

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitlinearw(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        lsfitlinearinternal<Precision>(y, w, fmatrix, n, m, info, c, rep);
    }


    /*************************************************************************
    Weighted constained linear least squares fitting.

    This  is  variation  of LSFitLinearW(), which searchs for min|A*x=b| given
    that  K  additional  constaints  C*x=bc are satisfied. It reduces original
    task to modified one: min|B*y-d| WITHOUT constraints,  then LSFitLinearW()
    is called.

    INPUT PARAMETERS:
        Y       -   array[0..N-1] Function values in  N  points.
        W       -   array[0..N-1]  Weights  corresponding to function  values.
                    Each summand in square  sum  of  approximation  deviations
                    from  given  values  is  multiplied  by  the   square   of
                    corresponding weight.
        FMatrix -   a table of basis functions values, array[0..N-1, 0..M-1].
                    FMatrix[I,J] - value of J-th basis function in I-th point.
        CMatrix -   a table of constaints, array[0..K-1,0..M].
                    I-th row of CMatrix corresponds to I-th linear constraint:
                    CMatrix[I,0]*C[0] + ... + CMatrix[I,M-1]*C[M-1] = CMatrix[I,M]
        N       -   number of points used. N>=1.
        M       -   number of basis functions, M>=1.
        K       -   number of constraints, 0 <= K < M
                    K=0 corresponds to absence of constraints.

    OUTPUT PARAMETERS:
        Info    -   error code:
                    * -4    internal SVD decomposition subroutine failed (very
                            rare and for degenerate systems only)
                    * -3    either   too   many  constraints  (M   or   more),
                            degenerate  constraints   (some   constraints  are
                            repetead twice) or inconsistent  constraints  were
                            specified.
                    * -1    incorrect N/M/K were specified
                    *  1    task is solved
        C       -   decomposition coefficients, array[0..M-1]
        Rep     -   fitting report. Following fields are set:
                    * RMSError          rms error on the (X,Y).
                    * AvgError          average error on the (X,Y).
                    * AvgRelError       average relative error on the non-zero Y
                    * MaxError          maximum error
                                        NON-WEIGHTED ERRORS ARE CALCULATED

    IMPORTANT:
        this subroitine doesn't calculate task's condition number for K<>0.

    SEE ALSO
        LSFitLinear
        LSFitLinearC
        LSFitLinearWC

      -- ALGLIB --
         Copyright 07.09.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitlinearwc(ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        ap::template_2d_array< amp::ampf<Precision> > cmatrix,
        int n,
        int m,
        int k,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        int i;
        int j;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > f2;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        ap::template_1d_array< amp::ampf<Precision> > c0;
        amp::ampf<Precision> v;


        if( n<1 || m<1 || k<0 )
        {
            info = -1;
            return;
        }
        if( k>=m )
        {
            info = -3;
            return;
        }
        
        //
        // Solve
        //
        if( k==0 )
        {
            
            //
            // no constraints
            //
            lsfitlinearinternal<Precision>(y, w, fmatrix, n, m, info, c, rep);
        }
        else
        {
            
            //
            // First, find general form solution of constraints system:
            // * factorize C = L*Q
            // * unpack Q
            // * fill upper part of C with zeros (for RCond)
            //
            // We got C=C0+Q2'*y where Q2 is lower M-K rows of Q.
            //
            ortfac::rmatrixlq<Precision>(cmatrix, k, m, tau);
            ortfac::rmatrixlqunpackq<Precision>(cmatrix, k, m, tau, m, q);
            for(i=0; i<=k-1; i++)
            {
                for(j=i+1; j<=m-1; j++)
                {
                    cmatrix(i,j) = amp::ampf<Precision>("0.0");
                }
            }
            if( rcond::rmatrixlurcondinf<Precision>(cmatrix, k)<1000*amp::ampf<Precision>::getAlgoPascalEpsilon() )
            {
                info = -3;
                return;
            }
            tmp.setlength(k);
            for(i=0; i<=k-1; i++)
            {
                if( i>0 )
                {
                    v = amp::vdotproduct(cmatrix.getrow(i, 0, i-1), tmp.getvector(0, i-1));
                }
                else
                {
                    v = 0;
                }
                tmp(i) = (cmatrix(i,m)-v)/cmatrix(i,i);
            }
            c0.setlength(m);
            for(i=0; i<=m-1; i++)
            {
                c0(i) = 0;
            }
            for(i=0; i<=k-1; i++)
            {
                v = tmp(i);
                amp::vadd(c0.getvector(0, m-1), q.getrow(i, 0, m-1), v);
            }
            
            //
            // Second, prepare modified matrix F2 = F*Q2' and solve modified task
            //
            tmp.setlength(ap::maxint(n, m)+1);
            f2.setlength(n, m-k);
            blas::matrixvectormultiply<Precision>(fmatrix, 0, n-1, 0, m-1, false, c0, 0, m-1, -amp::ampf<Precision>("1.0"), y, 0, n-1, amp::ampf<Precision>("1.0"));
            blas::matrixmatrixmultiply<Precision>(fmatrix, 0, n-1, 0, m-1, false, q, k, m-1, 0, m-1, true, amp::ampf<Precision>("1.0"), f2, 0, n-1, 0, m-k-1, amp::ampf<Precision>("0.0"), tmp);
            lsfitlinearinternal<Precision>(y, w, f2, n, m-k, info, tmp, rep);
            rep.taskrcond = -1;
            if( info<=0 )
            {
                return;
            }
            
            //
            // then, convert back to original answer: C = C0 + Q2'*Y0
            //
            c.setlength(m);
            amp::vmove(c.getvector(0, m-1), c0.getvector(0, m-1));
            blas::matrixvectormultiply<Precision>(q, k, m-1, 0, m-1, true, tmp, 0, m-k-1, amp::ampf<Precision>("1.0"), c, 0, m-1, amp::ampf<Precision>("1.0"));
        }
    }


    /*************************************************************************
    Linear least squares fitting, without weights.

    See LSFitLinearW for more information.

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitlinear(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        ap::template_1d_array< amp::ampf<Precision> > w;
        int i;


        if( n<1 )
        {
            info = -1;
            return;
        }
        w.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = 1;
        }
        lsfitlinearinternal<Precision>(y, w, fmatrix, n, m, info, c, rep);
    }


    /*************************************************************************
    Constained linear least squares fitting, without weights.

    See LSFitLinearWC() for more information.

      -- ALGLIB --
         Copyright 07.09.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitlinearc(ap::template_1d_array< amp::ampf<Precision> > y,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        const ap::template_2d_array< amp::ampf<Precision> >& cmatrix,
        int n,
        int m,
        int k,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        ap::template_1d_array< amp::ampf<Precision> > w;
        int i;


        if( n<1 )
        {
            info = -1;
            return;
        }
        w.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            w(i) = 1;
        }
        lsfitlinearwc<Precision>(y, w, fmatrix, cmatrix, n, m, k, info, c, rep);
    }


    /*************************************************************************
    Weighted nonlinear least squares fitting using gradient and Hessian.

    Nonlinear task min(F(c)) is solved, where

        F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,
        
        * N is a number of points,
        * M is a dimension of a space points belong to,
        * K is a dimension of a space of parameters being fitted,
        * w is an N-dimensional vector of weight coefficients,
        * x is a set of N points, each of them is an M-dimensional vector,
        * c is a K-dimensional vector of parameters being fitted
        
    This subroutine uses only f(x[i],c) and its gradient.
        
    INPUT PARAMETERS:
        X       -   array[0..N-1,0..M-1], points (one row = one point)
        Y       -   array[0..N-1], function values.
        W       -   weights, array[0..N-1]
        C       -   array[0..K-1], initial approximation to the solution,
        N       -   number of points, N>1
        M       -   dimension of space
        K       -   number of parameters being fitted
        CheapFG -   boolean flag, which is:
                    * True  if both function and gradient calculation complexity
                            are less than O(M^2).  An improved  algorithm  can
                            be  used  which corresponds  to  FGJ  scheme  from
                            MINLM unit.
                    * False otherwise.
                            Standard Jacibian-bases  Levenberg-Marquardt  algo
                            will be used (FJ scheme).

    OUTPUT PARAMETERS:
        State   -   structure which stores algorithm state between subsequent
                    calls  of   LSFitNonlinearIteration.   Used  for  reverse
                    communication.  This  structure   should   be  passed  to
                    LSFitNonlinearIteration subroutine.

    See also:
        LSFitNonlinearIteration
        LSFitNonlinearResults
        LSFitNonlinearFG (fitting without weights)
        LSFitNonlinearWFGH (fitting using Hessian)
        LSFitNonlinearFGH (fitting using Hessian, without weights)


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearwfg(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        bool cheapfg,
        lsfitstate<Precision>& state)
    {
        int i;


        state.n = n;
        state.m = m;
        state.k = k;
        lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
        lsfitnonlinearsetstpmax<Precision>(state, amp::ampf<Precision>("0.0"));
        state.cheapfg = cheapfg;
        state.havehess = false;
        if( n>=1 && m>=1 && k>=1 )
        {
            state.taskx.setlength(n, m);
            state.tasky.setlength(n);
            state.w.setlength(n);
            state.c.setlength(k);
            amp::vmove(state.c.getvector(0, k-1), c.getvector(0, k-1));
            amp::vmove(state.w.getvector(0, n-1), w.getvector(0, n-1));
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(state.taskx.getrow(i, 0, m-1), x.getrow(i, 0, m-1));
                state.tasky(i) = y(i);
            }
        }
        state.rstate.ia.setbounds(0, 4);
        state.rstate.ra.setbounds(0, 1);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Nonlinear least squares fitting, no individual weights.
    See LSFitNonlinearWFG for more information.

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearfg(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        bool cheapfg,
        lsfitstate<Precision>& state)
    {
        int i;


        state.n = n;
        state.m = m;
        state.k = k;
        lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
        lsfitnonlinearsetstpmax<Precision>(state, amp::ampf<Precision>("0.0"));
        state.cheapfg = cheapfg;
        state.havehess = false;
        if( n>=1 && m>=1 && k>=1 )
        {
            state.taskx.setlength(n, m);
            state.tasky.setlength(n);
            state.w.setlength(n);
            state.c.setlength(k);
            amp::vmove(state.c.getvector(0, k-1), c.getvector(0, k-1));
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(state.taskx.getrow(i, 0, m-1), x.getrow(i, 0, m-1));
                state.tasky(i) = y(i);
                state.w(i) = 1;
            }
        }
        state.rstate.ia.setbounds(0, 4);
        state.rstate.ra.setbounds(0, 1);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Weighted nonlinear least squares fitting using gradient/Hessian.

    Nonlinear task min(F(c)) is solved, where

        F(c) = (w[0]*(f(x[0],c)-y[0]))^2 + ... + (w[n-1]*(f(x[n-1],c)-y[n-1]))^2,

        * N is a number of points,
        * M is a dimension of a space points belong to,
        * K is a dimension of a space of parameters being fitted,
        * w is an N-dimensional vector of weight coefficients,
        * x is a set of N points, each of them is an M-dimensional vector,
        * c is a K-dimensional vector of parameters being fitted

    This subroutine uses f(x[i],c), its gradient and its Hessian.

    See LSFitNonlinearWFG() subroutine for information about function
    parameters.

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearwfgh(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        lsfitstate<Precision>& state)
    {
        int i;


        state.n = n;
        state.m = m;
        state.k = k;
        lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
        lsfitnonlinearsetstpmax<Precision>(state, amp::ampf<Precision>("0.0"));
        state.cheapfg = true;
        state.havehess = true;
        if( n>=1 && m>=1 && k>=1 )
        {
            state.taskx.setlength(n, m);
            state.tasky.setlength(n);
            state.w.setlength(n);
            state.c.setlength(k);
            amp::vmove(state.c.getvector(0, k-1), c.getvector(0, k-1));
            amp::vmove(state.w.getvector(0, n-1), w.getvector(0, n-1));
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(state.taskx.getrow(i, 0, m-1), x.getrow(i, 0, m-1));
                state.tasky(i) = y(i);
            }
        }
        state.rstate.ia.setbounds(0, 4);
        state.rstate.ra.setbounds(0, 1);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Nonlinear least squares fitting using gradient/Hessian without  individual
    weights. See LSFitNonlinearWFGH() for more information.


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearfgh(const ap::template_2d_array< amp::ampf<Precision> >& x,
        const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& c,
        int n,
        int m,
        int k,
        lsfitstate<Precision>& state)
    {
        int i;


        state.n = n;
        state.m = m;
        state.k = k;
        lsfitnonlinearsetcond<Precision>(state, amp::ampf<Precision>("0.0"), amp::ampf<Precision>("0.0"), 0);
        lsfitnonlinearsetstpmax<Precision>(state, amp::ampf<Precision>("0.0"));
        state.cheapfg = true;
        state.havehess = true;
        if( n>=1 && m>=1 && k>=1 )
        {
            state.taskx.setlength(n, m);
            state.tasky.setlength(n);
            state.w.setlength(n);
            state.c.setlength(k);
            amp::vmove(state.c.getvector(0, k-1), c.getvector(0, k-1));
            for(i=0; i<=n-1; i++)
            {
                amp::vmove(state.taskx.getrow(i, 0, m-1), x.getrow(i, 0, m-1));
                state.tasky(i) = y(i);
                state.w(i) = 1;
            }
        }
        state.rstate.ia.setbounds(0, 4);
        state.rstate.ra.setbounds(0, 1);
        state.rstate.stage = -1;
    }


    /*************************************************************************
    Stopping conditions for nonlinear least squares fitting.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be initialized
                    with LSFitNonLinearCreate???()
        EpsF    -   stopping criterion. Algorithm stops if
                    |F(k+1)-F(k)| <= EpsF*max{|F(k)|, |F(k+1)|, 1}
        EpsX    -   stopping criterion. Algorithm stops if
                    |X(k+1)-X(k)| <= EpsX*(1+|X(k)|)
        MaxIts  -   stopping criterion. Algorithm stops after MaxIts iterations.
                    MaxIts=0 means no stopping criterion.

    NOTE

    Passing EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to automatic
    stopping criterion selection (according to the scheme used by MINLM unit).


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearsetcond(lsfitstate<Precision>& state,
        amp::ampf<Precision> epsf,
        amp::ampf<Precision> epsx,
        int maxits)
    {
        ap::ap_error::make_assertion(epsf>=0);
        ap::ap_error::make_assertion(epsx>=0);
        ap::ap_error::make_assertion(maxits>=0);
        state.epsf = epsf;
        state.epsx = epsx;
        state.maxits = maxits;
    }


    /*************************************************************************
    This function sets maximum step length

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between calls and
                    which is used for reverse communication. Must be
                    initialized with LSFitNonLinearCreate???()
        StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                    want to limit step length.

    Use this subroutine when you optimize target function which contains exp()
    or  other  fast  growing  functions,  and optimization algorithm makes too
    large  steps  which  leads  to overflow. This function allows us to reject
    steps  that  are  too  large  (and  therefore  expose  us  to the possible
    overflow) without actually calculating function value at the x+stp*d.

    NOTE: non-zero StpMax leads to moderate  performance  degradation  because
    intermediate  step  of  preconditioned L-BFGS optimization is incompatible
    with limits on step size.

      -- ALGLIB --
         Copyright 02.04.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearsetstpmax(lsfitstate<Precision>& state,
        amp::ampf<Precision> stpmax)
    {
        ap::ap_error::make_assertion(stpmax>=0);
        state.stpmax = stpmax;
    }


    /*************************************************************************
    Nonlinear least squares fitting. Algorithm iteration.

    Called after inialization of the State structure with  LSFitNonlinearXXX()
    subroutine. See HTML docs for examples.

    INPUT PARAMETERS:
        State   -   structure which stores algorithm state between  subsequent
                    calls and which is used for reverse communication. Must be
                    initialized with LSFitNonlinearXXX() call first.

    RESULT
    1. If subroutine returned False, iterative algorithm has converged.
    2. If subroutine returned True, then if:
    * if State.NeedF=True,      function value F(X,C) is required
    * if State.NeedFG=True,     function value F(X,C) and gradient  dF/dC(X,C)
                                are required
    * if State.NeedFGH=True     function value F(X,C), gradient dF/dC(X,C) and
                                Hessian are required

    One and only one of this fields can be set at time.

    Function, its gradient and Hessian are calculated at  (X,C),  where  X  is
    stored in State.X[0..M-1] and C is stored in State.C[0..K-1].

    Results are stored:
    * function value            -   in State.F
    * gradient                  -   in State.G[0..K-1]
    * Hessian                   -   in State.H[0..K-1,0..K-1]

      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool lsfitnonlineariteration(lsfitstate<Precision>& state)
    {
        bool result;
        int n;
        int m;
        int k;
        int i;
        int j;
        amp::ampf<Precision> v;
        amp::ampf<Precision> relcnt;


        
        //
        // Reverse communication preparations
        // I know it looks ugly, but it works the same way
        // anywhere from C++ to Python.
        //
        // This code initializes locals by:
        // * random values determined during code
        //   generation - on first subroutine call
        // * values from previous call - on subsequent calls
        //
        if( state.rstate.stage>=0 )
        {
            n = state.rstate.ia(0);
            m = state.rstate.ia(1);
            k = state.rstate.ia(2);
            i = state.rstate.ia(3);
            j = state.rstate.ia(4);
            v = state.rstate.ra(0);
            relcnt = state.rstate.ra(1);
        }
        else
        {
            n = -983;
            m = -989;
            k = -834;
            i = 900;
            j = -287;
            v = 364;
            relcnt = 214;
        }
        if( state.rstate.stage==0 )
        {
            goto lbl_0;
        }
        if( state.rstate.stage==1 )
        {
            goto lbl_1;
        }
        if( state.rstate.stage==2 )
        {
            goto lbl_2;
        }
        if( state.rstate.stage==3 )
        {
            goto lbl_3;
        }
        if( state.rstate.stage==4 )
        {
            goto lbl_4;
        }
        
        //
        // Routine body
        //
        
        //
        // check params
        //
        if( state.n<1 || state.m<1 || state.k<1 || state.epsf<0 || state.epsx<0 || state.maxits<0 )
        {
            state.repterminationtype = -1;
            result = false;
            return result;
        }
        
        //
        // init
        //
        n = state.n;
        m = state.m;
        k = state.k;
        state.x.setlength(m);
        state.g.setlength(k);
        if( state.havehess )
        {
            state.h.setlength(k, k);
        }
        
        //
        // initialize LM optimizer
        //
        if( state.havehess )
        {
            
            //
            // use Hessian.
            // transform stopping conditions.
            //
            minlm::minlmcreatefgh<Precision>(k, state.c, state.optstate);
        }
        else
        {
            
            //
            // use one of gradient-based schemes (depending on gradient cost).
            // transform stopping conditions.
            //
            if( state.cheapfg )
            {
                minlm::minlmcreatefgj<Precision>(k, n, state.c, state.optstate);
            }
            else
            {
                minlm::minlmcreatefj<Precision>(k, n, state.c, state.optstate);
            }
        }
        minlm::minlmsetcond<Precision>(state.optstate, amp::ampf<Precision>("0.0"), state.epsf, state.epsx, state.maxits);
        minlm::minlmsetstpmax<Precision>(state.optstate, state.stpmax);
        
        //
        // Optimize
        //
    lbl_5:
        if( ! minlm::minlmiteration<Precision>(state.optstate) )
        {
            goto lbl_6;
        }
        if( ! state.optstate.needf )
        {
            goto lbl_7;
        }
        
        //
        // calculate F = sum (wi*(f(xi,c)-yi))^2
        //
        state.optstate.f = 0;
        i = 0;
    lbl_9:
        if( i>n-1 )
        {
            goto lbl_11;
        }
        amp::vmove(state.c.getvector(0, k-1), state.optstate.x.getvector(0, k-1));
        amp::vmove(state.x.getvector(0, m-1), state.taskx.getrow(i, 0, m-1));
        state.pointindex = i;
        lsfitclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 0;
        goto lbl_rcomm;
    lbl_0:
        state.optstate.f = state.optstate.f+amp::sqr<Precision>(state.w(i)*(state.f-state.tasky(i)));
        i = i+1;
        goto lbl_9;
    lbl_11:
        goto lbl_5;
    lbl_7:
        if( ! state.optstate.needfg )
        {
            goto lbl_12;
        }
        
        //
        // calculate F/gradF
        //
        state.optstate.f = 0;
        for(i=0; i<=k-1; i++)
        {
            state.optstate.g(i) = 0;
        }
        i = 0;
    lbl_14:
        if( i>n-1 )
        {
            goto lbl_16;
        }
        amp::vmove(state.c.getvector(0, k-1), state.optstate.x.getvector(0, k-1));
        amp::vmove(state.x.getvector(0, m-1), state.taskx.getrow(i, 0, m-1));
        state.pointindex = i;
        lsfitclearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 1;
        goto lbl_rcomm;
    lbl_1:
        state.optstate.f = state.optstate.f+amp::sqr<Precision>(state.w(i)*(state.f-state.tasky(i)));
        v = amp::sqr<Precision>(state.w(i))*2*(state.f-state.tasky(i));
        amp::vadd(state.optstate.g.getvector(0, k-1), state.g.getvector(0, k-1), v);
        i = i+1;
        goto lbl_14;
    lbl_16:
        goto lbl_5;
    lbl_12:
        if( ! state.optstate.needfij )
        {
            goto lbl_17;
        }
        
        //
        // calculate Fi/jac(Fi)
        //
        i = 0;
    lbl_19:
        if( i>n-1 )
        {
            goto lbl_21;
        }
        amp::vmove(state.c.getvector(0, k-1), state.optstate.x.getvector(0, k-1));
        amp::vmove(state.x.getvector(0, m-1), state.taskx.getrow(i, 0, m-1));
        state.pointindex = i;
        lsfitclearrequestfields<Precision>(state);
        state.needfg = true;
        state.rstate.stage = 2;
        goto lbl_rcomm;
    lbl_2:
        state.optstate.fi(i) = state.w(i)*(state.f-state.tasky(i));
        v = state.w(i);
        amp::vmove(state.optstate.j.getrow(i, 0, k-1), state.g.getvector(0, k-1), v);
        i = i+1;
        goto lbl_19;
    lbl_21:
        goto lbl_5;
    lbl_17:
        if( ! state.optstate.needfgh )
        {
            goto lbl_22;
        }
        
        //
        // calculate F/grad(F)/hess(F)
        //
        state.optstate.f = 0;
        for(i=0; i<=k-1; i++)
        {
            state.optstate.g(i) = 0;
        }
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                state.optstate.h(i,j) = 0;
            }
        }
        i = 0;
    lbl_24:
        if( i>n-1 )
        {
            goto lbl_26;
        }
        amp::vmove(state.c.getvector(0, k-1), state.optstate.x.getvector(0, k-1));
        amp::vmove(state.x.getvector(0, m-1), state.taskx.getrow(i, 0, m-1));
        state.pointindex = i;
        lsfitclearrequestfields<Precision>(state);
        state.needfgh = true;
        state.rstate.stage = 3;
        goto lbl_rcomm;
    lbl_3:
        state.optstate.f = state.optstate.f+amp::sqr<Precision>(state.w(i)*(state.f-state.tasky(i)));
        v = amp::sqr<Precision>(state.w(i))*2*(state.f-state.tasky(i));
        amp::vadd(state.optstate.g.getvector(0, k-1), state.g.getvector(0, k-1), v);
        for(j=0; j<=k-1; j++)
        {
            v = 2*amp::sqr<Precision>(state.w(i))*state.g(j);
            amp::vadd(state.optstate.h.getrow(j, 0, k-1), state.g.getvector(0, k-1), v);
            v = 2*amp::sqr<Precision>(state.w(i))*(state.f-state.tasky(i));
            amp::vadd(state.optstate.h.getrow(j, 0, k-1), state.h.getrow(j, 0, k-1), v);
        }
        i = i+1;
        goto lbl_24;
    lbl_26:
        goto lbl_5;
    lbl_22:
        goto lbl_5;
    lbl_6:
        minlm::minlmresults<Precision>(state.optstate, state.c, state.optrep);
        state.repterminationtype = state.optrep.terminationtype;
        
        //
        // calculate errors
        //
        if( state.repterminationtype<=0 )
        {
            goto lbl_27;
        }
        state.reprmserror = 0;
        state.repavgerror = 0;
        state.repavgrelerror = 0;
        state.repmaxerror = 0;
        relcnt = 0;
        i = 0;
    lbl_29:
        if( i>n-1 )
        {
            goto lbl_31;
        }
        amp::vmove(state.c.getvector(0, k-1), state.c.getvector(0, k-1));
        amp::vmove(state.x.getvector(0, m-1), state.taskx.getrow(i, 0, m-1));
        state.pointindex = i;
        lsfitclearrequestfields<Precision>(state);
        state.needf = true;
        state.rstate.stage = 4;
        goto lbl_rcomm;
    lbl_4:
        v = state.f;
        state.reprmserror = state.reprmserror+amp::sqr<Precision>(v-state.tasky(i));
        state.repavgerror = state.repavgerror+amp::abs<Precision>(v-state.tasky(i));
        if( state.tasky(i)!=0 )
        {
            state.repavgrelerror = state.repavgrelerror+amp::abs<Precision>(v-state.tasky(i))/amp::abs<Precision>(state.tasky(i));
            relcnt = relcnt+1;
        }
        state.repmaxerror = amp::maximum<Precision>(state.repmaxerror, amp::abs<Precision>(v-state.tasky(i)));
        i = i+1;
        goto lbl_29;
    lbl_31:
        state.reprmserror = amp::sqrt<Precision>(state.reprmserror/n);
        state.repavgerror = state.repavgerror/n;
        if( relcnt!=0 )
        {
            state.repavgrelerror = state.repavgrelerror/relcnt;
        }
    lbl_27:
        result = false;
        return result;
        
        //
        // Saving state
        //
    lbl_rcomm:
        result = true;
        state.rstate.ia(0) = n;
        state.rstate.ia(1) = m;
        state.rstate.ia(2) = k;
        state.rstate.ia(3) = i;
        state.rstate.ia(4) = j;
        state.rstate.ra(0) = v;
        state.rstate.ra(1) = relcnt;
        return result;
    }


    /*************************************************************************
    Nonlinear least squares fitting results.

    Called after LSFitNonlinearIteration() returned False.

    INPUT PARAMETERS:
        State   -   algorithm state (used by LSFitNonlinearIteration).

    OUTPUT PARAMETERS:
        Info    -   completetion code:
                        * -1    incorrect parameters were specified
                        *  1    relative function improvement is no more than
                                EpsF.
                        *  2    relative step is no more than EpsX.
                        *  4    gradient norm is no more than EpsG
                        *  5    MaxIts steps was taken
        C       -   array[0..K-1], solution
        Rep     -   optimization report. Following fields are set:
                    * Rep.TerminationType completetion code:
                    * RMSError          rms error on the (X,Y).
                    * AvgError          average error on the (X,Y).
                    * AvgRelError       average relative error on the non-zero Y
                    * MaxError          maximum error
                                        NON-WEIGHTED ERRORS ARE CALCULATED


      -- ALGLIB --
         Copyright 17.08.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitnonlinearresults(const lsfitstate<Precision>& state,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        info = state.repterminationtype;
        if( info>0 )
        {
            c.setlength(state.k);
            amp::vmove(c.getvector(0, state.k-1), state.c.getvector(0, state.k-1));
            rep.rmserror = state.reprmserror;
            rep.avgerror = state.repavgerror;
            rep.avgrelerror = state.repavgrelerror;
            rep.maxerror = state.repmaxerror;
        }
    }


    template<unsigned int Precision>
    void lsfitscalexy(ap::template_1d_array< amp::ampf<Precision> >& x,
        ap::template_1d_array< amp::ampf<Precision> >& y,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& xc,
        ap::template_1d_array< amp::ampf<Precision> >& yc,
        const ap::template_1d_array< int >& dc,
        int k,
        amp::ampf<Precision>& xa,
        amp::ampf<Precision>& xb,
        amp::ampf<Precision>& sa,
        amp::ampf<Precision>& sb,
        ap::template_1d_array< amp::ampf<Precision> >& xoriginal,
        ap::template_1d_array< amp::ampf<Precision> >& yoriginal)
    {
        amp::ampf<Precision> xmin;
        amp::ampf<Precision> xmax;
        int i;


        ap::ap_error::make_assertion(n>=1);
        ap::ap_error::make_assertion(k>=0);
        
        //
        // Calculate xmin/xmax.
        // Force xmin<>xmax.
        //
        xmin = x(0);
        xmax = x(0);
        for(i=1; i<=n-1; i++)
        {
            xmin = amp::minimum<Precision>(xmin, x(i));
            xmax = amp::maximum<Precision>(xmax, x(i));
        }
        for(i=0; i<=k-1; i++)
        {
            xmin = amp::minimum<Precision>(xmin, xc(i));
            xmax = amp::maximum<Precision>(xmax, xc(i));
        }
        if( xmin==xmax )
        {
            if( xmin==0 )
            {
                xmin = -1;
                xmax = +1;
            }
            else
            {
                xmin = amp::ampf<Precision>("0.5")*xmin;
            }
        }
        
        //
        // Transform abscissas: map [XA,XB] to [0,1]
        //
        // Store old X[] in XOriginal[] (it will be used
        // to calculate relative error).
        //
        xoriginal.setlength(n);
        amp::vmove(xoriginal.getvector(0, n-1), x.getvector(0, n-1));
        xa = xmin;
        xb = xmax;
        for(i=0; i<=n-1; i++)
        {
            x(i) = 2*(x(i)-amp::ampf<Precision>("0.5")*(xa+xb))/(xb-xa);
        }
        for(i=0; i<=k-1; i++)
        {
            ap::ap_error::make_assertion(dc(i)>=0);
            xc(i) = 2*(xc(i)-amp::ampf<Precision>("0.5")*(xa+xb))/(xb-xa);
            yc(i) = yc(i)*amp::pow<Precision>(amp::ampf<Precision>("0.5")*(xb-xa), amp::ampf<Precision>(dc(i)));
        }
        
        //
        // Transform function values: map [SA,SB] to [0,1]
        // SA = mean(Y),
        // SB = SA+stddev(Y).
        //
        // Store old Y[] in YOriginal[] (it will be used
        // to calculate relative error).
        //
        yoriginal.setlength(n);
        amp::vmove(yoriginal.getvector(0, n-1), y.getvector(0, n-1));
        sa = 0;
        for(i=0; i<=n-1; i++)
        {
            sa = sa+y(i);
        }
        sa = sa/n;
        sb = 0;
        for(i=0; i<=n-1; i++)
        {
            sb = sb+amp::sqr<Precision>(y(i)-sa);
        }
        sb = amp::sqrt<Precision>(sb/n)+sa;
        if( sb==sa )
        {
            sb = 2*sa;
        }
        if( sb==sa )
        {
            sb = sa+1;
        }
        for(i=0; i<=n-1; i++)
        {
            y(i) = (y(i)-sa)/(sb-sa);
        }
        for(i=0; i<=k-1; i++)
        {
            if( dc(i)==0 )
            {
                yc(i) = (yc(i)-sa)/(sb-sa);
            }
            else
            {
                yc(i) = yc(i)/(sb-sa);
            }
        }
    }


    /*************************************************************************
    Internal fitting subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitlinearinternal(const ap::template_1d_array< amp::ampf<Precision> >& y,
        const ap::template_1d_array< amp::ampf<Precision> >& w,
        const ap::template_2d_array< amp::ampf<Precision> >& fmatrix,
        int n,
        int m,
        int& info,
        ap::template_1d_array< amp::ampf<Precision> >& c,
        lsfitreport<Precision>& rep)
    {
        amp::ampf<Precision> threshold;
        ap::template_2d_array< amp::ampf<Precision> > ft;
        ap::template_2d_array< amp::ampf<Precision> > q;
        ap::template_2d_array< amp::ampf<Precision> > l;
        ap::template_2d_array< amp::ampf<Precision> > r;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > wmod;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        int i;
        int j;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > sv;
        ap::template_2d_array< amp::ampf<Precision> > u;
        ap::template_2d_array< amp::ampf<Precision> > vt;
        ap::template_1d_array< amp::ampf<Precision> > tmp;
        ap::template_1d_array< amp::ampf<Precision> > utb;
        ap::template_1d_array< amp::ampf<Precision> > sutb;
        int relcnt;


        if( n<1 || m<1 )
        {
            info = -1;
            return;
        }
        info = 1;
        threshold = amp::sqrt<Precision>(amp::ampf<Precision>::getAlgoPascalEpsilon());
        
        //
        // Degenerate case, needs special handling
        //
        if( n<m )
        {
            
            //
            // Create design matrix.
            //
            ft.setlength(n, m);
            b.setlength(n);
            wmod.setlength(n);
            for(j=0; j<=n-1; j++)
            {
                v = w(j);
                amp::vmove(ft.getrow(j, 0, m-1), fmatrix.getrow(j, 0, m-1), v);
                b(j) = w(j)*y(j);
                wmod(j) = 1;
            }
            
            //
            // LQ decomposition and reduction to M=N
            //
            c.setlength(m);
            for(i=0; i<=m-1; i++)
            {
                c(i) = 0;
            }
            rep.taskrcond = 0;
            ortfac::rmatrixlq<Precision>(ft, n, m, tau);
            ortfac::rmatrixlqunpackq<Precision>(ft, n, m, tau, n, q);
            ortfac::rmatrixlqunpackl<Precision>(ft, n, m, l);
            lsfitlinearinternal<Precision>(b, wmod, l, n, n, info, tmp, rep);
            if( info<=0 )
            {
                return;
            }
            for(i=0; i<=n-1; i++)
            {
                v = tmp(i);
                amp::vadd(c.getvector(0, m-1), q.getrow(i, 0, m-1), v);
            }
            return;
        }
        
        //
        // N>=M. Generate design matrix and reduce to N=M using
        // QR decomposition.
        //
        ft.setlength(n, m);
        b.setlength(n);
        for(j=0; j<=n-1; j++)
        {
            v = w(j);
            amp::vmove(ft.getrow(j, 0, m-1), fmatrix.getrow(j, 0, m-1), v);
            b(j) = w(j)*y(j);
        }
        ortfac::rmatrixqr<Precision>(ft, n, m, tau);
        ortfac::rmatrixqrunpackq<Precision>(ft, n, m, tau, m, q);
        ortfac::rmatrixqrunpackr<Precision>(ft, n, m, r);
        tmp.setlength(m);
        for(i=0; i<=m-1; i++)
        {
            tmp(i) = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            v = b(i);
            amp::vadd(tmp.getvector(0, m-1), q.getrow(i, 0, m-1), v);
        }
        b.setlength(m);
        amp::vmove(b.getvector(0, m-1), tmp.getvector(0, m-1));
        
        //
        // R contains reduced MxM design upper triangular matrix,
        // B contains reduced Mx1 right part.
        //
        // Determine system condition number and decide
        // should we use triangular solver (faster) or
        // SVD-based solver (more stable).
        //
        // We can use LU-based RCond estimator for this task.
        //
        rep.taskrcond = rcond::rmatrixlurcondinf<Precision>(r, m);
        if( rep.taskrcond>threshold )
        {
            
            //
            // use QR-based solver
            //
            c.setlength(m);
            c(m-1) = b(m-1)/r(m-1,m-1);
            for(i=m-2; i>=0; i--)
            {
                v = amp::vdotproduct(r.getrow(i, i+1, m-1), c.getvector(i+1, m-1));
                c(i) = (b(i)-v)/r(i,i);
            }
        }
        else
        {
            
            //
            // use SVD-based solver
            //
            if( !svd::rmatrixsvd<Precision>(r, m, m, 1, 1, 2, sv, u, vt) )
            {
                info = -4;
                return;
            }
            utb.setlength(m);
            sutb.setlength(m);
            for(i=0; i<=m-1; i++)
            {
                utb(i) = 0;
            }
            for(i=0; i<=m-1; i++)
            {
                v = b(i);
                amp::vadd(utb.getvector(0, m-1), u.getrow(i, 0, m-1), v);
            }
            if( sv(0)>0 )
            {
                rep.taskrcond = sv(m-1)/sv(0);
                for(i=0; i<=m-1; i++)
                {
                    if( sv(i)>threshold*sv(0) )
                    {
                        sutb(i) = utb(i)/sv(i);
                    }
                    else
                    {
                        sutb(i) = 0;
                    }
                }
            }
            else
            {
                rep.taskrcond = 0;
                for(i=0; i<=m-1; i++)
                {
                    sutb(i) = 0;
                }
            }
            c.setlength(m);
            for(i=0; i<=m-1; i++)
            {
                c(i) = 0;
            }
            for(i=0; i<=m-1; i++)
            {
                v = sutb(i);
                amp::vadd(c.getvector(0, m-1), vt.getrow(i, 0, m-1), v);
            }
        }
        
        //
        // calculate errors
        //
        rep.rmserror = 0;
        rep.avgerror = 0;
        rep.avgrelerror = 0;
        rep.maxerror = 0;
        relcnt = 0;
        for(i=0; i<=n-1; i++)
        {
            v = amp::vdotproduct(fmatrix.getrow(i, 0, m-1), c.getvector(0, m-1));
            rep.rmserror = rep.rmserror+amp::sqr<Precision>(v-y(i));
            rep.avgerror = rep.avgerror+amp::abs<Precision>(v-y(i));
            if( y(i)!=0 )
            {
                rep.avgrelerror = rep.avgrelerror+amp::abs<Precision>(v-y(i))/amp::abs<Precision>(y(i));
                relcnt = relcnt+1;
            }
            rep.maxerror = amp::maximum<Precision>(rep.maxerror, amp::abs<Precision>(v-y(i)));
        }
        rep.rmserror = amp::sqrt<Precision>(rep.rmserror/n);
        rep.avgerror = rep.avgerror/n;
        if( relcnt!=0 )
        {
            rep.avgrelerror = rep.avgrelerror/relcnt;
        }
    }


    /*************************************************************************
    Internal subroutine
    *************************************************************************/
    template<unsigned int Precision>
    void lsfitclearrequestfields(lsfitstate<Precision>& state)
    {
        state.needf = false;
        state.needfg = false;
        state.needfgh = false;
    }
} // namespace

#endif
