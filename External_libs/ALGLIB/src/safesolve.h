/*************************************************************************
This file is a part of ALGLIB project.

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

#ifndef _safesolve_h
#define _safesolve_h

#include "ap.h"
#include "amp.h"
namespace safesolve
{
    template<unsigned int Precision>
    bool rmatrixscaledtrsafesolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        amp::ampf<Precision> sa,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        bool isupper,
        int trans,
        bool isunit,
        amp::ampf<Precision> maxgrowth);
    template<unsigned int Precision>
    bool cmatrixscaledtrsafesolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        amp::ampf<Precision> sa,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& x,
        bool isupper,
        int trans,
        bool isunit,
        amp::ampf<Precision> maxgrowth);
    template<unsigned int Precision>
    bool cbasicsolveandupdate(amp::campf<Precision> alpha,
        amp::campf<Precision> beta,
        amp::ampf<Precision> lnmax,
        amp::ampf<Precision> bnorm,
        amp::ampf<Precision> maxgrowth,
        amp::ampf<Precision>& xnorm,
        amp::campf<Precision>& x);


    /*************************************************************************
    Real implementation of CMatrixScaledTRSafeSolve

      -- ALGLIB routine --
         21.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool rmatrixscaledtrsafesolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        amp::ampf<Precision> sa,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        bool isupper,
        int trans,
        bool isunit,
        amp::ampf<Precision> maxgrowth)
    {
        bool result;
        amp::ampf<Precision> lnmax;
        amp::ampf<Precision> nrmb;
        amp::ampf<Precision> nrmx;
        int i;
        amp::campf<Precision> alpha;
        amp::campf<Precision> beta;
        amp::ampf<Precision> vr;
        amp::campf<Precision> cx;
        ap::template_1d_array< amp::ampf<Precision> > tmp;


        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(trans==0 || trans==1);
        result = true;
        lnmax = amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
        
        //
        // Quick return if possible
        //
        if( n<=0 )
        {
            return result;
        }
        
        //
        // Load norms: right part and X
        //
        nrmb = 0;
        for(i=0; i<=n-1; i++)
        {
            nrmb = amp::maximum<Precision>(nrmb, amp::abs<Precision>(x(i)));
        }
        nrmx = 0;
        
        //
        // Solve
        //
        tmp.setlength(n);
        result = true;
        if( isupper && trans==0 )
        {
            
            //
            // U*x = b
            //
            for(i=n-1; i>=0; i--)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                if( i<n-1 )
                {
                    amp::vmove(tmp.getvector(i+1, n-1), a.getrow(i, i+1, n-1), sa);
                    vr = amp::vdotproduct(tmp.getvector(i+1, n-1), x.getvector(i+1, n-1));
                    beta = x(i)-vr;
                }
                else
                {
                    beta = x(i);
                }
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
                if( !result )
                {
                    return result;
                }
                x(i) = cx.x;
            }
            return result;
        }
        if( !isupper && trans==0 )
        {
            
            //
            // L*x = b
            //
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                if( i>0 )
                {
                    amp::vmove(tmp.getvector(0, i-1), a.getrow(i, 0, i-1), sa);
                    vr = amp::vdotproduct(tmp.getvector(0, i-1), x.getvector(0, i-1));
                    beta = x(i)-vr;
                }
                else
                {
                    beta = x(i);
                }
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
                if( !result )
                {
                    return result;
                }
                x(i) = cx.x;
            }
            return result;
        }
        if( isupper && trans==1 )
        {
            
            //
            // U^T*x = b
            //
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
                if( !result )
                {
                    return result;
                }
                x(i) = cx.x;
                
                //
                // update the rest of right part
                //
                if( i<n-1 )
                {
                    vr = cx.x;
                    amp::vmove(tmp.getvector(i+1, n-1), a.getrow(i, i+1, n-1), sa);
                    amp::vsub(x.getvector(i+1, n-1), tmp.getvector(i+1, n-1), vr);
                }
            }
            return result;
        }
        if( !isupper && trans==1 )
        {
            
            //
            // L^T*x = b
            //
            for(i=n-1; i>=0; i--)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, cx);
                if( !result )
                {
                    return result;
                }
                x(i) = cx.x;
                
                //
                // update the rest of right part
                //
                if( i>0 )
                {
                    vr = cx.x;
                    amp::vmove(tmp.getvector(0, i-1), a.getrow(i, 0, i-1), sa);
                    amp::vsub(x.getvector(0, i-1), tmp.getvector(0, i-1), vr);
                }
            }
            return result;
        }
        result = false;
        return result;
    }


    /*************************************************************************
    Internal subroutine for safe solution of

        SA*op(A)=b
        
    where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
    either identity transform, transposition or Hermitian transposition, SA is
    a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

    This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
    MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
    near-degenerate matrices are handled correctly (False is returned) as long
    as MaxGrowth is significantly less than MaxRealNumber/norm(b).

      -- ALGLIB routine --
         21.01.2010
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cmatrixscaledtrsafesolve(const ap::template_2d_array< amp::campf<Precision> >& a,
        amp::ampf<Precision> sa,
        int n,
        ap::template_1d_array< amp::campf<Precision> >& x,
        bool isupper,
        int trans,
        bool isunit,
        amp::ampf<Precision> maxgrowth)
    {
        bool result;
        amp::ampf<Precision> lnmax;
        amp::ampf<Precision> nrmb;
        amp::ampf<Precision> nrmx;
        int i;
        amp::campf<Precision> alpha;
        amp::campf<Precision> beta;
        amp::campf<Precision> vc;
        ap::template_1d_array< amp::campf<Precision> > tmp;
        int i_;


        ap::ap_error::make_assertion(n>0);
        ap::ap_error::make_assertion(trans==0 || trans==1 || trans==2);
        result = true;
        lnmax = amp::log<Precision>(amp::ampf<Precision>::getAlgoPascalMaxNumber());
        
        //
        // Quick return if possible
        //
        if( n<=0 )
        {
            return result;
        }
        
        //
        // Load norms: right part and X
        //
        nrmb = 0;
        for(i=0; i<=n-1; i++)
        {
            nrmb = amp::maximum<Precision>(nrmb, amp::abscomplex<Precision>(x(i)));
        }
        nrmx = 0;
        
        //
        // Solve
        //
        tmp.setlength(n);
        result = true;
        if( isupper && trans==0 )
        {
            
            //
            // U*x = b
            //
            for(i=n-1; i>=0; i--)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                if( i<n-1 )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        tmp(i_) = sa*a(i,i_);
                    }
                    vc = 0.0;
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        vc += tmp(i_)*x(i_);
                    }
                    beta = x(i)-vc;
                }
                else
                {
                    beta = x(i);
                }
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
            }
            return result;
        }
        if( !isupper && trans==0 )
        {
            
            //
            // L*x = b
            //
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                if( i>0 )
                {
                    for(i_=0; i_<=i-1;i_++)
                    {
                        tmp(i_) = sa*a(i,i_);
                    }
                    vc = 0.0;
                    for(i_=0; i_<=i-1;i_++)
                    {
                        vc += tmp(i_)*x(i_);
                    }
                    beta = x(i)-vc;
                }
                else
                {
                    beta = x(i);
                }
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
            }
            return result;
        }
        if( isupper && trans==1 )
        {
            
            //
            // U^T*x = b
            //
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
                
                //
                // update the rest of right part
                //
                if( i<n-1 )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        tmp(i_) = sa*a(i,i_);
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        x(i_) = x(i_) - vc*tmp(i_);
                    }
                }
            }
            return result;
        }
        if( !isupper && trans==1 )
        {
            
            //
            // L^T*x = b
            //
            for(i=n-1; i>=0; i--)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = a(i,i)*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
                
                //
                // update the rest of right part
                //
                if( i>0 )
                {
                    for(i_=0; i_<=i-1;i_++)
                    {
                        tmp(i_) = sa*a(i,i_);
                    }
                    for(i_=0; i_<=i-1;i_++)
                    {
                        x(i_) = x(i_) - vc*tmp(i_);
                    }
                }
            }
            return result;
        }
        if( isupper && trans==2 )
        {
            
            //
            // U^H*x = b
            //
            for(i=0; i<=n-1; i++)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = amp::conj<Precision>(a(i,i))*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
                
                //
                // update the rest of right part
                //
                if( i<n-1 )
                {
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        tmp(i_) = sa*amp::conj(a(i,i_));
                    }
                    for(i_=i+1; i_<=n-1;i_++)
                    {
                        x(i_) = x(i_) - vc*tmp(i_);
                    }
                }
            }
            return result;
        }
        if( !isupper && trans==2 )
        {
            
            //
            // L^T*x = b
            //
            for(i=n-1; i>=0; i--)
            {
                
                //
                // Task is reduced to alpha*x[i] = beta
                //
                if( isunit )
                {
                    alpha = sa;
                }
                else
                {
                    alpha = amp::conj<Precision>(a(i,i))*sa;
                }
                beta = x(i);
                
                //
                // solve alpha*x[i] = beta
                //
                result = cbasicsolveandupdate<Precision>(alpha, beta, lnmax, nrmb, maxgrowth, nrmx, vc);
                if( !result )
                {
                    return result;
                }
                x(i) = vc;
                
                //
                // update the rest of right part
                //
                if( i>0 )
                {
                    for(i_=0; i_<=i-1;i_++)
                    {
                        tmp(i_) = sa*amp::conj(a(i,i_));
                    }
                    for(i_=0; i_<=i-1;i_++)
                    {
                        x(i_) = x(i_) - vc*tmp(i_);
                    }
                }
            }
            return result;
        }
        result = false;
        return result;
    }


    /*************************************************************************
    complex basic solver-updater for reduced linear system

        alpha*x[i] = beta

    solves this equation and updates it in overlfow-safe manner (keeping track
    of relative growth of solution).

    Parameters:
        Alpha   -   alpha
        Beta    -   beta
        LnMax   -   precomputed Ln(MaxRealNumber)
        BNorm   -   inf-norm of b (right part of original system)
        MaxGrowth-  maximum growth of norm(x) relative to norm(b)
        XNorm   -   inf-norm of other components of X (which are already processed)
                    it is updated by CBasicSolveAndUpdate.
        X       -   solution

      -- ALGLIB routine --
         26.01.2009
         Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool cbasicsolveandupdate(amp::campf<Precision> alpha,
        amp::campf<Precision> beta,
        amp::ampf<Precision> lnmax,
        amp::ampf<Precision> bnorm,
        amp::ampf<Precision> maxgrowth,
        amp::ampf<Precision>& xnorm,
        amp::campf<Precision>& x)
    {
        bool result;
        amp::ampf<Precision> v;


        result = false;
        if( alpha==0 )
        {
            return result;
        }
        if( beta!=0 )
        {
            
            //
            // alpha*x[i]=beta
            //
            v = amp::log<Precision>(amp::abscomplex<Precision>(beta))-amp::log<Precision>(amp::abscomplex<Precision>(alpha));
            if( v>lnmax )
            {
                return result;
            }
            x = beta/alpha;
        }
        else
        {
            
            //
            // alpha*x[i]=0
            //
            x = 0;
        }
        
        //
        // update NrmX, test growth limit
        //
        xnorm = amp::maximum<Precision>(xnorm, amp::abscomplex<Precision>(x));
        if( xnorm>maxgrowth*bnorm )
        {
            return result;
        }
        result = true;
        return result;
    }
} // namespace

#endif
