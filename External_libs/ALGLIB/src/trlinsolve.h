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

#ifndef _trlinsolve_h
#define _trlinsolve_h

#include "ap.h"
#include "amp.h"
namespace trlinsolve
{
    template<unsigned int Precision>
    void rmatrixtrsafesolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& s,
        bool isupper,
        bool istrans,
        bool isunit);
    template<unsigned int Precision>
    void safesolvetriangular(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& s,
        bool isupper,
        bool istrans,
        bool isunit,
        bool normin,
        ap::template_1d_array< amp::ampf<Precision> >& cnorm);


    /*************************************************************************
    Utility subroutine performing the "safe" solution of system of linear
    equations with triangular coefficient matrices.

    The subroutine uses scaling and solves the scaled system A*x=s*b (where  s
    is  a  scalar  value)  instead  of  A*x=b,  choosing  s  so  that x can be
    represented by a floating-point number. The closer the system  gets  to  a
    singular, the less s is. If the system is singular, s=0 and x contains the
    non-trivial solution of equation A*x=0.

    The feature of an algorithm is that it could not cause an  overflow  or  a
    division by zero regardless of the matrix used as the input.

    The algorithm can solve systems of equations with  upper/lower  triangular
    matrices,  with/without unit diagonal, and systems of type A*x=b or A'*x=b
    (where A' is a transposed matrix A).

    Input parameters:
        A       -   system matrix. Array whose indexes range within [0..N-1, 0..N-1].
        N       -   size of matrix A.
        X       -   right-hand member of a system.
                    Array whose index ranges within [0..N-1].
        IsUpper -   matrix type. If it is True, the system matrix is the upper
                    triangular and is located in  the  corresponding  part  of
                    matrix A.
        Trans   -   problem type. If it is True, the problem to be  solved  is
                    A'*x=b, otherwise it is A*x=b.
        Isunit  -   matrix type. If it is True, the system matrix has  a  unit
                    diagonal (the elements on the main diagonal are  not  used
                    in the calculation process), otherwise the matrix is considered
                    to be a general triangular matrix.

    Output parameters:
        X       -   solution. Array whose index ranges within [0..N-1].
        S       -   scaling factor.

      -- LAPACK auxiliary routine (version 3.0) --
         Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
         Courant Institute, Argonne National Lab, and Rice University
         June 30, 1992
    *************************************************************************/
    template<unsigned int Precision>
    void rmatrixtrsafesolve(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& s,
        bool isupper,
        bool istrans,
        bool isunit)
    {
        bool normin;
        ap::template_1d_array< amp::ampf<Precision> > cnorm;
        ap::template_2d_array< amp::ampf<Precision> > a1;
        ap::template_1d_array< amp::ampf<Precision> > x1;
        int i;


        
        //
        // From 0-based to 1-based
        //
        normin = false;
        a1.setbounds(1, n, 1, n);
        x1.setbounds(1, n);
        for(i=1; i<=n; i++)
        {
            amp::vmove(a1.getrow(i, 1, n), a.getrow(i-1, 0, n-1));
        }
        amp::vmove(x1.getvector(1, n), x.getvector(0, n-1));
        
        //
        // Solve 1-based
        //
        safesolvetriangular<Precision>(a1, n, x1, s, isupper, istrans, isunit, normin, cnorm);
        
        //
        // From 1-based to 0-based
        //
        amp::vmove(x.getvector(0, n-1), x1.getvector(1, n));
    }


    /*************************************************************************
    Obsolete 1-based subroutine.
    See RMatrixTRSafeSolve for 0-based replacement.
    *************************************************************************/
    template<unsigned int Precision>
    void safesolvetriangular(const ap::template_2d_array< amp::ampf<Precision> >& a,
        int n,
        ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision>& s,
        bool isupper,
        bool istrans,
        bool isunit,
        bool normin,
        ap::template_1d_array< amp::ampf<Precision> >& cnorm)
    {
        int i;
        int imax;
        int j;
        int jfirst;
        int jinc;
        int jlast;
        int jm1;
        int jp1;
        int ip1;
        int im1;
        int k;
        int flg;
        amp::ampf<Precision> v;
        amp::ampf<Precision> vd;
        amp::ampf<Precision> bignum;
        amp::ampf<Precision> grow;
        amp::ampf<Precision> rec;
        amp::ampf<Precision> smlnum;
        amp::ampf<Precision> sumj;
        amp::ampf<Precision> tjj;
        amp::ampf<Precision> tjjs;
        amp::ampf<Precision> tmax;
        amp::ampf<Precision> tscal;
        amp::ampf<Precision> uscal;
        amp::ampf<Precision> xbnd;
        amp::ampf<Precision> xj;
        amp::ampf<Precision> xmax;
        bool notran;
        bool upper;
        bool nounit;


        upper = isupper;
        notran = !istrans;
        nounit = !isunit;
        
        //
        // Quick return if possible
        //
        if( n==0 )
        {
            return;
        }
        
        //
        // Determine machine dependent parameters to control overflow.
        //
        smlnum = amp::ampf<Precision>::getAlgoPascalMinNumber()/(amp::ampf<Precision>::getAlgoPascalEpsilon()*2);
        bignum = 1/smlnum;
        s = 1;
        if( !normin )
        {
            cnorm.setbounds(1, n);
            
            //
            // Compute the 1-norm of each column, not including the diagonal.
            //
            if( upper )
            {
                
                //
                // A is upper triangular.
                //
                for(j=1; j<=n; j++)
                {
                    v = 0;
                    for(k=1; k<=j-1; k++)
                    {
                        v = v+amp::abs<Precision>(a(k,j));
                    }
                    cnorm(j) = v;
                }
            }
            else
            {
                
                //
                // A is lower triangular.
                //
                for(j=1; j<=n-1; j++)
                {
                    v = 0;
                    for(k=j+1; k<=n; k++)
                    {
                        v = v+amp::abs<Precision>(a(k,j));
                    }
                    cnorm(j) = v;
                }
                cnorm(n) = 0;
            }
        }
        
        //
        // Scale the column norms by TSCAL if the maximum element in CNORM is
        // greater than BIGNUM.
        //
        imax = 1;
        for(k=2; k<=n; k++)
        {
            if( cnorm(k)>cnorm(imax) )
            {
                imax = k;
            }
        }
        tmax = cnorm(imax);
        if( tmax<=bignum )
        {
            tscal = 1;
        }
        else
        {
            tscal = 1/(smlnum*tmax);
            amp::vmul(cnorm.getvector(1, n), tscal);
        }
        
        //
        // Compute a bound on the computed solution vector to see if the
        // Level 2 BLAS routine DTRSV can be used.
        //
        j = 1;
        for(k=2; k<=n; k++)
        {
            if( amp::abs<Precision>(x(k))>amp::abs<Precision>(x(j)) )
            {
                j = k;
            }
        }
        xmax = amp::abs<Precision>(x(j));
        xbnd = xmax;
        if( notran )
        {
            
            //
            // Compute the growth in A * x = b.
            //
            if( upper )
            {
                jfirst = n;
                jlast = 1;
                jinc = -1;
            }
            else
            {
                jfirst = 1;
                jlast = n;
                jinc = 1;
            }
            if( tscal!=1 )
            {
                grow = 0;
            }
            else
            {
                if( nounit )
                {
                    
                    //
                    // A is non-unit triangular.
                    //
                    // Compute GROW = 1/G(j) and XBND = 1/M(j).
                    // Initially, G(0) = max{x(i), i=1,...,n}.
                    //
                    grow = 1/amp::maximum<Precision>(xbnd, smlnum);
                    xbnd = grow;
                    j = jfirst;
                    while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                    {
                        
                        //
                        // Exit the loop if the growth factor is too small.
                        //
                        if( grow<=smlnum )
                        {
                            break;
                        }
                        
                        //
                        // M(j) = G(j-1) / abs(A(j,j))
                        //
                        tjj = amp::abs<Precision>(a(j,j));
                        xbnd = amp::minimum<Precision>(xbnd, amp::minimum<Precision>(amp::ampf<Precision>(1), tjj)*grow);
                        if( tjj+cnorm(j)>=smlnum )
                        {
                            
                            //
                            // G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                            //
                            grow = grow*(tjj/(tjj+cnorm(j)));
                        }
                        else
                        {
                            
                            //
                            // G(j) could overflow, set GROW to 0.
                            //
                            grow = 0;
                        }
                        if( j==jlast )
                        {
                            grow = xbnd;
                        }
                        j = j+jinc;
                    }
                }
                else
                {
                    
                    //
                    // A is unit triangular.
                    //
                    // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                    //
                    grow = amp::minimum<Precision>(amp::ampf<Precision>(1), 1/amp::maximum<Precision>(xbnd, smlnum));
                    j = jfirst;
                    while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                    {
                        
                        //
                        // Exit the loop if the growth factor is too small.
                        //
                        if( grow<=smlnum )
                        {
                            break;
                        }
                        
                        //
                        // G(j) = G(j-1)*( 1 + CNORM(j) )
                        //
                        grow = grow*(1/(1+cnorm(j)));
                        j = j+jinc;
                    }
                }
            }
        }
        else
        {
            
            //
            // Compute the growth in A' * x = b.
            //
            if( upper )
            {
                jfirst = 1;
                jlast = n;
                jinc = 1;
            }
            else
            {
                jfirst = n;
                jlast = 1;
                jinc = -1;
            }
            if( tscal!=1 )
            {
                grow = 0;
            }
            else
            {
                if( nounit )
                {
                    
                    //
                    // A is non-unit triangular.
                    //
                    // Compute GROW = 1/G(j) and XBND = 1/M(j).
                    // Initially, M(0) = max{x(i), i=1,...,n}.
                    //
                    grow = 1/amp::maximum<Precision>(xbnd, smlnum);
                    xbnd = grow;
                    j = jfirst;
                    while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                    {
                        
                        //
                        // Exit the loop if the growth factor is too small.
                        //
                        if( grow<=smlnum )
                        {
                            break;
                        }
                        
                        //
                        // G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
                        //
                        xj = 1+cnorm(j);
                        grow = amp::minimum<Precision>(grow, xbnd/xj);
                        
                        //
                        // M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
                        //
                        tjj = amp::abs<Precision>(a(j,j));
                        if( xj>tjj )
                        {
                            xbnd = xbnd*(tjj/xj);
                        }
                        if( j==jlast )
                        {
                            grow = amp::minimum<Precision>(grow, xbnd);
                        }
                        j = j+jinc;
                    }
                }
                else
                {
                    
                    //
                    // A is unit triangular.
                    //
                    // Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                    //
                    grow = amp::minimum<Precision>(amp::ampf<Precision>(1), 1/amp::maximum<Precision>(xbnd, smlnum));
                    j = jfirst;
                    while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                    {
                        
                        //
                        // Exit the loop if the growth factor is too small.
                        //
                        if( grow<=smlnum )
                        {
                            break;
                        }
                        
                        //
                        // G(j) = ( 1 + CNORM(j) )*G(j-1)
                        //
                        xj = 1+cnorm(j);
                        grow = grow/xj;
                        j = j+jinc;
                    }
                }
            }
        }
        if( grow*tscal>smlnum )
        {
            
            //
            // Use the Level 2 BLAS solve if the reciprocal of the bound on
            // elements of X is not too small.
            //
            if( upper && notran || !upper && !notran )
            {
                if( nounit )
                {
                    vd = a(n,n);
                }
                else
                {
                    vd = 1;
                }
                x(n) = x(n)/vd;
                for(i=n-1; i>=1; i--)
                {
                    ip1 = i+1;
                    if( upper )
                    {
                        v = amp::vdotproduct(a.getrow(i, ip1, n), x.getvector(ip1, n));
                    }
                    else
                    {
                        v = amp::vdotproduct(a.getcolumn(i, ip1, n), x.getvector(ip1, n));
                    }
                    if( nounit )
                    {
                        vd = a(i,i);
                    }
                    else
                    {
                        vd = 1;
                    }
                    x(i) = (x(i)-v)/vd;
                }
            }
            else
            {
                if( nounit )
                {
                    vd = a(1,1);
                }
                else
                {
                    vd = 1;
                }
                x(1) = x(1)/vd;
                for(i=2; i<=n; i++)
                {
                    im1 = i-1;
                    if( upper )
                    {
                        v = amp::vdotproduct(a.getcolumn(i, 1, im1), x.getvector(1, im1));
                    }
                    else
                    {
                        v = amp::vdotproduct(a.getrow(i, 1, im1), x.getvector(1, im1));
                    }
                    if( nounit )
                    {
                        vd = a(i,i);
                    }
                    else
                    {
                        vd = 1;
                    }
                    x(i) = (x(i)-v)/vd;
                }
            }
        }
        else
        {
            
            //
            // Use a Level 1 BLAS solve, scaling intermediate results.
            //
            if( xmax>bignum )
            {
                
                //
                // Scale X so that its components are less than or equal to
                // BIGNUM in absolute value.
                //
                s = bignum/xmax;
                amp::vmul(x.getvector(1, n), s);
                xmax = bignum;
            }
            if( notran )
            {
                
                //
                // Solve A * x = b
                //
                j = jfirst;
                while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                {
                    
                    //
                    // Compute x(j) = b(j) / A(j,j), scaling x if necessary.
                    //
                    xj = amp::abs<Precision>(x(j));
                    flg = 0;
                    if( nounit )
                    {
                        tjjs = a(j,j)*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                        if( tscal==1 )
                        {
                            flg = 100;
                        }
                    }
                    if( flg!=100 )
                    {
                        tjj = amp::abs<Precision>(tjjs);
                        if( tjj>smlnum )
                        {
                            
                            //
                            // abs(A(j,j)) > SMLNUM:
                            //
                            if( tjj<1 )
                            {
                                if( xj>tjj*bignum )
                                {
                                    
                                    //
                                    // Scale x by 1/b(j).
                                    //
                                    rec = 1/xj;
                                    amp::vmul(x.getvector(1, n), rec);
                                    s = s*rec;
                                    xmax = xmax*rec;
                                }
                            }
                            x(j) = x(j)/tjjs;
                            xj = amp::abs<Precision>(x(j));
                        }
                        else
                        {
                            if( tjj>0 )
                            {
                                
                                //
                                // 0 < abs(A(j,j)) <= SMLNUM:
                                //
                                if( xj>tjj*bignum )
                                {
                                    
                                    //
                                    // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                                    // to avoid overflow when dividing by A(j,j).
                                    //
                                    rec = tjj*bignum/xj;
                                    if( cnorm(j)>1 )
                                    {
                                        
                                        //
                                        // Scale by 1/CNORM(j) to avoid overflow when
                                        // multiplying x(j) times column j.
                                        //
                                        rec = rec/cnorm(j);
                                    }
                                    amp::vmul(x.getvector(1, n), rec);
                                    s = s*rec;
                                    xmax = xmax*rec;
                                }
                                x(j) = x(j)/tjjs;
                                xj = amp::abs<Precision>(x(j));
                            }
                            else
                            {
                                
                                //
                                // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                                // scale = 0, and compute a solution to A*x = 0.
                                //
                                for(i=1; i<=n; i++)
                                {
                                    x(i) = 0;
                                }
                                x(j) = 1;
                                xj = 1;
                                s = 0;
                                xmax = 0;
                            }
                        }
                    }
                    
                    //
                    // Scale x if necessary to avoid overflow when adding a
                    // multiple of column j of A.
                    //
                    if( xj>1 )
                    {
                        rec = 1/xj;
                        if( cnorm(j)>(bignum-xmax)*rec )
                        {
                            
                            //
                            // Scale x by 1/(2*abs(x(j))).
                            //
                            rec = rec*amp::ampf<Precision>("0.5");
                            amp::vmul(x.getvector(1, n), rec);
                            s = s*rec;
                        }
                    }
                    else
                    {
                        if( xj*cnorm(j)>bignum-xmax )
                        {
                            
                            //
                            // Scale x by 1/2.
                            //
                            amp::vmul(x.getvector(1, n), 0.5);
                            s = s*amp::ampf<Precision>("0.5");
                        }
                    }
                    if( upper )
                    {
                        if( j>1 )
                        {
                            
                            //
                            // Compute the update
                            // x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                            //
                            v = x(j)*tscal;
                            jm1 = j-1;
                            amp::vsub(x.getvector(1, jm1), a.getcolumn(j, 1, jm1), v);
                            i = 1;
                            for(k=2; k<=j-1; k++)
                            {
                                if( amp::abs<Precision>(x(k))>amp::abs<Precision>(x(i)) )
                                {
                                    i = k;
                                }
                            }
                            xmax = amp::abs<Precision>(x(i));
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            
                            //
                            // Compute the update
                            // x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                            //
                            jp1 = j+1;
                            v = x(j)*tscal;
                            amp::vsub(x.getvector(jp1, n), a.getcolumn(j, jp1, n), v);
                            i = j+1;
                            for(k=j+2; k<=n; k++)
                            {
                                if( amp::abs<Precision>(x(k))>amp::abs<Precision>(x(i)) )
                                {
                                    i = k;
                                }
                            }
                            xmax = amp::abs<Precision>(x(i));
                        }
                    }
                    j = j+jinc;
                }
            }
            else
            {
                
                //
                // Solve A' * x = b
                //
                j = jfirst;
                while( jinc>0 && j<=jlast || jinc<0 && j>=jlast )
                {
                    
                    //
                    // Compute x(j) = b(j) - sum A(k,j)*x(k).
                    //   k<>j
                    //
                    xj = amp::abs<Precision>(x(j));
                    uscal = tscal;
                    rec = 1/amp::maximum<Precision>(xmax, amp::ampf<Precision>(1));
                    if( cnorm(j)>(bignum-xj)*rec )
                    {
                        
                        //
                        // If x(j) could overflow, scale x by 1/(2*XMAX).
                        //
                        rec = rec*amp::ampf<Precision>("0.5");
                        if( nounit )
                        {
                            tjjs = a(j,j)*tscal;
                        }
                        else
                        {
                            tjjs = tscal;
                        }
                        tjj = amp::abs<Precision>(tjjs);
                        if( tjj>1 )
                        {
                            
                            //
                            // Divide by A(j,j) when scaling x if A(j,j) > 1.
                            //
                            rec = amp::minimum<Precision>(amp::ampf<Precision>(1), rec*tjj);
                            uscal = uscal/tjjs;
                        }
                        if( rec<1 )
                        {
                            amp::vmul(x.getvector(1, n), rec);
                            s = s*rec;
                            xmax = xmax*rec;
                        }
                    }
                    sumj = 0;
                    if( uscal==1 )
                    {
                        
                        //
                        // If the scaling needed for A in the dot product is 1,
                        // call DDOT to perform the dot product.
                        //
                        if( upper )
                        {
                            if( j>1 )
                            {
                                jm1 = j-1;
                                sumj = amp::vdotproduct(a.getcolumn(j, 1, jm1), x.getvector(1, jm1));
                            }
                            else
                            {
                                sumj = 0;
                            }
                        }
                        else
                        {
                            if( j<n )
                            {
                                jp1 = j+1;
                                sumj = amp::vdotproduct(a.getcolumn(j, jp1, n), x.getvector(jp1, n));
                            }
                        }
                    }
                    else
                    {
                        
                        //
                        // Otherwise, use in-line code for the dot product.
                        //
                        if( upper )
                        {
                            for(i=1; i<=j-1; i++)
                            {
                                v = a(i,j)*uscal;
                                sumj = sumj+v*x(i);
                            }
                        }
                        else
                        {
                            if( j<n )
                            {
                                for(i=j+1; i<=n; i++)
                                {
                                    v = a(i,j)*uscal;
                                    sumj = sumj+v*x(i);
                                }
                            }
                        }
                    }
                    if( uscal==tscal )
                    {
                        
                        //
                        // Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                        // was not used to scale the dotproduct.
                        //
                        x(j) = x(j)-sumj;
                        xj = amp::abs<Precision>(x(j));
                        flg = 0;
                        if( nounit )
                        {
                            tjjs = a(j,j)*tscal;
                        }
                        else
                        {
                            tjjs = tscal;
                            if( tscal==1 )
                            {
                                flg = 150;
                            }
                        }
                        
                        //
                        // Compute x(j) = x(j) / A(j,j), scaling if necessary.
                        //
                        if( flg!=150 )
                        {
                            tjj = amp::abs<Precision>(tjjs);
                            if( tjj>smlnum )
                            {
                                
                                //
                                // abs(A(j,j)) > SMLNUM:
                                //
                                if( tjj<1 )
                                {
                                    if( xj>tjj*bignum )
                                    {
                                        
                                        //
                                        // Scale X by 1/abs(x(j)).
                                        //
                                        rec = 1/xj;
                                        amp::vmul(x.getvector(1, n), rec);
                                        s = s*rec;
                                        xmax = xmax*rec;
                                    }
                                }
                                x(j) = x(j)/tjjs;
                            }
                            else
                            {
                                if( tjj>0 )
                                {
                                    
                                    //
                                    // 0 < abs(A(j,j)) <= SMLNUM:
                                    //
                                    if( xj>tjj*bignum )
                                    {
                                        
                                        //
                                        // Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                                        //
                                        rec = tjj*bignum/xj;
                                        amp::vmul(x.getvector(1, n), rec);
                                        s = s*rec;
                                        xmax = xmax*rec;
                                    }
                                    x(j) = x(j)/tjjs;
                                }
                                else
                                {
                                    
                                    //
                                    // A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                                    // scale = 0, and compute a solution to A'*x = 0.
                                    //
                                    for(i=1; i<=n; i++)
                                    {
                                        x(i) = 0;
                                    }
                                    x(j) = 1;
                                    s = 0;
                                    xmax = 0;
                                }
                            }
                        }
                    }
                    else
                    {
                        
                        //
                        // Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                        // product has already been divided by 1/A(j,j).
                        //
                        x(j) = x(j)/tjjs-sumj;
                    }
                    xmax = amp::maximum<Precision>(xmax, amp::abs<Precision>(x(j)));
                    j = j+jinc;
                }
            }
            s = s/tscal;
        }
        
        //
        // Scale the column norms by 1/TSCAL for return.
        //
        if( tscal!=1 )
        {
            v = 1/tscal;
            amp::vmul(cnorm.getvector(1, n), v);
        }
    }
} // namespace

#endif
