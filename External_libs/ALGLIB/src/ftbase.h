/*************************************************************************
Copyright (c) 2009, Sergey Bochkanov (ALGLIB project).

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

#ifndef _ftbase_h
#define _ftbase_h

#include "ap.h"
#include "amp.h"
namespace ftbase
{
    template<unsigned int Precision>
    class ftplan
    {
    public:
        ap::template_1d_array< int > plan;
        ap::template_1d_array< amp::ampf<Precision> > precomputed;
        ap::template_1d_array< amp::ampf<Precision> > tmpbuf;
        ap::template_1d_array< amp::ampf<Precision> > stackbuf;
    };




    template<unsigned int Precision>
    void ftbasegeneratecomplexfftplan(int n,
        ftplan<Precision>& plan);
    template<unsigned int Precision>
    void ftbasegeneraterealfftplan(int n,
        ftplan<Precision>& plan);
    template<unsigned int Precision>
    void ftbasegeneraterealfhtplan(int n,
        ftplan<Precision>& plan);
    template<unsigned int Precision>
    void ftbaseexecuteplan(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        int n,
        ftplan<Precision>& plan);
    template<unsigned int Precision>
    void ftbaseexecuteplanrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        ftplan<Precision>& plan,
        int entryoffset,
        int stackptr);
    template<unsigned int Precision>
    void ftbasefactorize(int n,
        int tasktype,
        int& n1,
        int& n2);
    template<unsigned int Precision>
    bool ftbaseissmooth(int n);
    template<unsigned int Precision>
    int ftbasefindsmooth(int n);
    template<unsigned int Precision>
    int ftbasefindsmootheven(int n);
    template<unsigned int Precision>
    amp::ampf<Precision> ftbasegetflopestimate(int n);
    template<unsigned int Precision>
    void ftbasegenerateplanrec(int n,
        int tasktype,
        ftplan<Precision>& plan,
        int& plansize,
        int& precomputedsize,
        int& planarraysize,
        int& tmpmemsize,
        int& stackmemsize,
        int stackptr);
    template<unsigned int Precision>
    void ftbaseprecomputeplanrec(ftplan<Precision>& plan,
        int entryoffset,
        int stackptr);
    template<unsigned int Precision>
    void ffttwcalc(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        int n1,
        int n2);
    template<unsigned int Precision>
    void internalcomplexlintranspose(ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        int astart,
        ap::template_1d_array< amp::ampf<Precision> >& buf);
    template<unsigned int Precision>
    void internalreallintranspose(ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        int astart,
        ap::template_1d_array< amp::ampf<Precision> >& buf);
    template<unsigned int Precision>
    void ffticltrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int astart,
        int astride,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        int bstart,
        int bstride,
        int m,
        int n);
    template<unsigned int Precision>
    void fftirltrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int astart,
        int astride,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        int bstart,
        int bstride,
        int m,
        int n);
    template<unsigned int Precision>
    void ftbasefindsmoothrec(int n,
        int seed,
        int leastfactor,
        int& best);
    template<unsigned int Precision>
    void fftarrayresize(ap::template_1d_array< int >& a,
        int& asize,
        int newasize);
    template<unsigned int Precision>
    void reffht(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n,
        int offs);


    static const int ftbaseplanentrysize = 8;
    static const int ftbasecffttask = 0;
    static const int ftbaserfhttask = 1;
    static const int ftbaserffttask = 2;
    static const int fftcooleytukeyplan = 0;
    static const int fftbluesteinplan = 1;
    static const int fftcodeletplan = 2;
    static const int fhtcooleytukeyplan = 3;
    static const int fhtcodeletplan = 4;
    static const int fftrealcooleytukeyplan = 5;
    static const int fftemptyplan = 6;
    static const int fhtn2plan = 999;
    static const int ftbaseupdatetw = 4;
    static const int ftbasecodeletmax = 5;
    static const int ftbasecodeletrecommended = 5;
    template<unsigned int Precision>
    const amp::ampf<Precision>& ftbaseinefficiencyfactor()
    {
        static amp::ampf<Precision> v = amp::ampf<Precision>("1.3");
        return v;
    }
    static const int ftbasemaxsmoothfactor = 5;


    /*************************************************************************
    This subroutine generates FFT plan - a decomposition of a N-length FFT to
    the more simpler operations. Plan consists of the root entry and the child
    entries.

    Subroutine parameters:
        N               task size
        
    Output parameters:
        Plan            plan

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasegeneratecomplexfftplan(int n,
        ftplan<Precision>& plan)
    {
        int planarraysize;
        int plansize;
        int precomputedsize;
        int tmpmemsize;
        int stackmemsize;
        int stackptr;


        planarraysize = 1;
        plansize = 0;
        precomputedsize = 0;
        stackmemsize = 0;
        stackptr = 0;
        tmpmemsize = 2*n;
        plan.plan.setlength(planarraysize);
        ftbasegenerateplanrec<Precision>(n, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
        plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
        plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
        plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
        stackptr = 0;
        ftbaseprecomputeplanrec<Precision>(plan, 0, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
    }


    /*************************************************************************
    Generates real FFT plan
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasegeneraterealfftplan(int n,
        ftplan<Precision>& plan)
    {
        int planarraysize;
        int plansize;
        int precomputedsize;
        int tmpmemsize;
        int stackmemsize;
        int stackptr;


        planarraysize = 1;
        plansize = 0;
        precomputedsize = 0;
        stackmemsize = 0;
        stackptr = 0;
        tmpmemsize = 2*n;
        plan.plan.setlength(planarraysize);
        ftbasegenerateplanrec<Precision>(n, ftbaserffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
        plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
        plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
        plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
        stackptr = 0;
        ftbaseprecomputeplanrec<Precision>(plan, 0, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
    }


    /*************************************************************************
    Generates real FHT plan
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasegeneraterealfhtplan(int n,
        ftplan<Precision>& plan)
    {
        int planarraysize;
        int plansize;
        int precomputedsize;
        int tmpmemsize;
        int stackmemsize;
        int stackptr;


        planarraysize = 1;
        plansize = 0;
        precomputedsize = 0;
        stackmemsize = 0;
        stackptr = 0;
        tmpmemsize = n;
        plan.plan.setlength(planarraysize);
        ftbasegenerateplanrec<Precision>(n, ftbaserfhttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
        plan.stackbuf.setlength(ap::maxint(stackmemsize, 1));
        plan.tmpbuf.setlength(ap::maxint(tmpmemsize, 1));
        plan.precomputed.setlength(ap::maxint(precomputedsize, 1));
        stackptr = 0;
        ftbaseprecomputeplanrec<Precision>(plan, 0, stackptr);
        ap::ap_error::make_assertion(stackptr==0);
    }


    /*************************************************************************
    This subroutine executes FFT/FHT plan.

    If Plan is a:
    * complex FFT plan  -   sizeof(A)=2*N,
                            A contains interleaved real/imaginary values
    * real FFT plan     -   sizeof(A)=2*N,
                            A contains real values interleaved with zeros
    * real FHT plan     -   sizeof(A)=2*N,
                            A contains real values interleaved with zeros

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbaseexecuteplan(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        int n,
        ftplan<Precision>& plan)
    {
        int stackptr;


        stackptr = 0;
        ftbaseexecuteplanrec<Precision>(a, aoffset, plan, 0, stackptr);
    }


    /*************************************************************************
    Recurrent subroutine for the FTBaseExecutePlan

    Parameters:
        A           FFT'ed array
        AOffset     offset of the FFT'ed part (distance is measured in doubles)

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbaseexecuteplanrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        ftplan<Precision>& plan,
        int entryoffset,
        int stackptr)
    {
        int i;
        int j;
        int k;
        int n1;
        int n2;
        int n;
        int m;
        int offs;
        int offs1;
        int offs2;
        int offsa;
        int offsb;
        int offsp;
        amp::ampf<Precision> hk;
        amp::ampf<Precision> hnk;
        amp::ampf<Precision> x;
        amp::ampf<Precision> y;
        amp::ampf<Precision> bx;
        amp::ampf<Precision> by;
        ap::template_1d_array< amp::ampf<Precision> > emptyarray;
        amp::ampf<Precision> a0x;
        amp::ampf<Precision> a0y;
        amp::ampf<Precision> a1x;
        amp::ampf<Precision> a1y;
        amp::ampf<Precision> a2x;
        amp::ampf<Precision> a2y;
        amp::ampf<Precision> a3x;
        amp::ampf<Precision> a3y;
        amp::ampf<Precision> v0;
        amp::ampf<Precision> v1;
        amp::ampf<Precision> v2;
        amp::ampf<Precision> v3;
        amp::ampf<Precision> t1x;
        amp::ampf<Precision> t1y;
        amp::ampf<Precision> t2x;
        amp::ampf<Precision> t2y;
        amp::ampf<Precision> t3x;
        amp::ampf<Precision> t3y;
        amp::ampf<Precision> t4x;
        amp::ampf<Precision> t4y;
        amp::ampf<Precision> t5x;
        amp::ampf<Precision> t5y;
        amp::ampf<Precision> m1x;
        amp::ampf<Precision> m1y;
        amp::ampf<Precision> m2x;
        amp::ampf<Precision> m2y;
        amp::ampf<Precision> m3x;
        amp::ampf<Precision> m3y;
        amp::ampf<Precision> m4x;
        amp::ampf<Precision> m4y;
        amp::ampf<Precision> m5x;
        amp::ampf<Precision> m5y;
        amp::ampf<Precision> s1x;
        amp::ampf<Precision> s1y;
        amp::ampf<Precision> s2x;
        amp::ampf<Precision> s2y;
        amp::ampf<Precision> s3x;
        amp::ampf<Precision> s3y;
        amp::ampf<Precision> s4x;
        amp::ampf<Precision> s4y;
        amp::ampf<Precision> s5x;
        amp::ampf<Precision> s5y;
        amp::ampf<Precision> c1;
        amp::ampf<Precision> c2;
        amp::ampf<Precision> c3;
        amp::ampf<Precision> c4;
        amp::ampf<Precision> c5;
        ap::template_1d_array< amp::ampf<Precision> > tmp;


        if( plan.plan(entryoffset+3)==fftemptyplan )
        {
            return;
        }
        if( plan.plan(entryoffset+3)==fftcooleytukeyplan )
        {
            
            //
            // Cooley-Tukey plan
            // * transposition
            // * row-wise FFT
            // * twiddle factors:
            //   - TwBase is a basis twiddle factor for I=1, J=1
            //   - TwRow is a twiddle factor for a second element in a row (J=1)
            //   - Tw is a twiddle factor for a current element
            // * transposition again
            // * row-wise FFT again
            //
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            internalcomplexlintranspose<Precision>(a, n1, n2, aoffset, plan.tmpbuf);
            for(i=0; i<=n2-1; i++)
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+i*n1*2, plan, plan.plan(entryoffset+5), stackptr);
            }
            ffttwcalc<Precision>(a, aoffset, n1, n2);
            internalcomplexlintranspose<Precision>(a, n2, n1, aoffset, plan.tmpbuf);
            for(i=0; i<=n1-1; i++)
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+i*n2*2, plan, plan.plan(entryoffset+6), stackptr);
            }
            internalcomplexlintranspose<Precision>(a, n1, n2, aoffset, plan.tmpbuf);
            return;
        }
        if( plan.plan(entryoffset+3)==fftrealcooleytukeyplan )
        {
            
            //
            // Cooley-Tukey plan
            // * transposition
            // * row-wise FFT
            // * twiddle factors:
            //   - TwBase is a basis twiddle factor for I=1, J=1
            //   - TwRow is a twiddle factor for a second element in a row (J=1)
            //   - Tw is a twiddle factor for a current element
            // * transposition again
            // * row-wise FFT again
            //
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            internalcomplexlintranspose<Precision>(a, n2, n1, aoffset, plan.tmpbuf);
            for(i=0; i<=n1/2-1; i++)
            {
                
                //
                // pack two adjacent smaller real FFT's together,
                // make one complex FFT,
                // unpack result
                //
                offs = aoffset+2*i*n2*2;
                for(k=0; k<=n2-1; k++)
                {
                    a(offs+2*k+1) = a(offs+2*n2+2*k+0);
                }
                ftbaseexecuteplanrec<Precision>(a, offs, plan, plan.plan(entryoffset+6), stackptr);
                plan.tmpbuf(0) = a(offs+0);
                plan.tmpbuf(1) = 0;
                plan.tmpbuf(2*n2+0) = a(offs+1);
                plan.tmpbuf(2*n2+1) = 0;
                for(k=1; k<=n2-1; k++)
                {
                    offs1 = 2*k;
                    offs2 = 2*n2+2*k;
                    hk = a(offs+2*k+0);
                    hnk = a(offs+2*(n2-k)+0);
                    plan.tmpbuf(offs1+0) = +amp::ampf<Precision>("0.5")*(hk+hnk);
                    plan.tmpbuf(offs2+1) = -amp::ampf<Precision>("0.5")*(hk-hnk);
                    hk = a(offs+2*k+1);
                    hnk = a(offs+2*(n2-k)+1);
                    plan.tmpbuf(offs2+0) = +amp::ampf<Precision>("0.5")*(hk+hnk);
                    plan.tmpbuf(offs1+1) = +amp::ampf<Precision>("0.5")*(hk-hnk);
                }
                amp::vmove(a.getvector(offs, offs+2*n2*2-1), plan.tmpbuf.getvector(0, 2*n2*2-1));
            }
            if( n1%2!=0 )
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+(n1-1)*n2*2, plan, plan.plan(entryoffset+6), stackptr);
            }
            ffttwcalc<Precision>(a, aoffset, n2, n1);
            internalcomplexlintranspose<Precision>(a, n1, n2, aoffset, plan.tmpbuf);
            for(i=0; i<=n2-1; i++)
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+i*n1*2, plan, plan.plan(entryoffset+5), stackptr);
            }
            internalcomplexlintranspose<Precision>(a, n2, n1, aoffset, plan.tmpbuf);
            return;
        }
        if( plan.plan(entryoffset+3)==fhtcooleytukeyplan )
        {
            
            //
            // Cooley-Tukey FHT plan:
            // * transpose                    \
            // * smaller FHT's                |
            // * pre-process                  |
            // * multiply by twiddle factors  | corresponds to multiplication by H1
            // * post-process                 |
            // * transpose again              /
            // * multiply by H2 (smaller FHT's)
            // * final transposition
            //
            // For more details see Vitezslav Vesely, "Fast algorithms
            // of Fourier and Hartley transform and their implementation in MATLAB",
            // page 31.
            //
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            n = n1*n2;
            internalreallintranspose<Precision>(a, n1, n2, aoffset, plan.tmpbuf);
            for(i=0; i<=n2-1; i++)
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+i*n1, plan, plan.plan(entryoffset+5), stackptr);
            }
            for(i=0; i<=n2-1; i++)
            {
                for(j=0; j<=n1-1; j++)
                {
                    offsa = aoffset+i*n1;
                    hk = a(offsa+j);
                    hnk = a(offsa+(n1-j)%n1);
                    offs = 2*(i*n1+j);
                    plan.tmpbuf(offs+0) = -amp::ampf<Precision>("0.5")*(hnk-hk);
                    plan.tmpbuf(offs+1) = +amp::ampf<Precision>("0.5")*(hk+hnk);
                }
            }
            ffttwcalc<Precision>(plan.tmpbuf, 0, n1, n2);
            for(j=0; j<=n1-1; j++)
            {
                a(aoffset+j) = plan.tmpbuf(2*j+0)+plan.tmpbuf(2*j+1);
            }
            if( n2%2==0 )
            {
                offs = 2*(n2/2)*n1;
                offsa = aoffset+n2/2*n1;
                for(j=0; j<=n1-1; j++)
                {
                    a(offsa+j) = plan.tmpbuf(offs+2*j+0)+plan.tmpbuf(offs+2*j+1);
                }
            }
            for(i=1; i<=(n2+1)/2-1; i++)
            {
                offs = 2*i*n1;
                offs2 = 2*(n2-i)*n1;
                offsa = aoffset+i*n1;
                for(j=0; j<=n1-1; j++)
                {
                    a(offsa+j) = plan.tmpbuf(offs+2*j+1)+plan.tmpbuf(offs2+2*j+0);
                }
                offsa = aoffset+(n2-i)*n1;
                for(j=0; j<=n1-1; j++)
                {
                    a(offsa+j) = plan.tmpbuf(offs+2*j+0)+plan.tmpbuf(offs2+2*j+1);
                }
            }
            internalreallintranspose<Precision>(a, n2, n1, aoffset, plan.tmpbuf);
            for(i=0; i<=n1-1; i++)
            {
                ftbaseexecuteplanrec<Precision>(a, aoffset+i*n2, plan, plan.plan(entryoffset+6), stackptr);
            }
            internalreallintranspose<Precision>(a, n1, n2, aoffset, plan.tmpbuf);
            return;
        }
        if( plan.plan(entryoffset+3)==fhtn2plan )
        {
            
            //
            // Cooley-Tukey FHT plan
            //
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            n = n1*n2;
            reffht<Precision>(a, n, aoffset);
            return;
        }
        if( plan.plan(entryoffset+3)==fftcodeletplan )
        {
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            n = n1*n2;
            if( n==2 )
            {
                a0x = a(aoffset+0);
                a0y = a(aoffset+1);
                a1x = a(aoffset+2);
                a1y = a(aoffset+3);
                v0 = a0x+a1x;
                v1 = a0y+a1y;
                v2 = a0x-a1x;
                v3 = a0y-a1y;
                a(aoffset+0) = v0;
                a(aoffset+1) = v1;
                a(aoffset+2) = v2;
                a(aoffset+3) = v3;
                return;
            }
            if( n==3 )
            {
                offs = plan.plan(entryoffset+7);
                c1 = plan.precomputed(offs+0);
                c2 = plan.precomputed(offs+1);
                a0x = a(aoffset+0);
                a0y = a(aoffset+1);
                a1x = a(aoffset+2);
                a1y = a(aoffset+3);
                a2x = a(aoffset+4);
                a2y = a(aoffset+5);
                t1x = a1x+a2x;
                t1y = a1y+a2y;
                a0x = a0x+t1x;
                a0y = a0y+t1y;
                m1x = c1*t1x;
                m1y = c1*t1y;
                m2x = c2*(a1y-a2y);
                m2y = c2*(a2x-a1x);
                s1x = a0x+m1x;
                s1y = a0y+m1y;
                a1x = s1x+m2x;
                a1y = s1y+m2y;
                a2x = s1x-m2x;
                a2y = s1y-m2y;
                a(aoffset+0) = a0x;
                a(aoffset+1) = a0y;
                a(aoffset+2) = a1x;
                a(aoffset+3) = a1y;
                a(aoffset+4) = a2x;
                a(aoffset+5) = a2y;
                return;
            }
            if( n==4 )
            {
                a0x = a(aoffset+0);
                a0y = a(aoffset+1);
                a1x = a(aoffset+2);
                a1y = a(aoffset+3);
                a2x = a(aoffset+4);
                a2y = a(aoffset+5);
                a3x = a(aoffset+6);
                a3y = a(aoffset+7);
                t1x = a0x+a2x;
                t1y = a0y+a2y;
                t2x = a1x+a3x;
                t2y = a1y+a3y;
                m2x = a0x-a2x;
                m2y = a0y-a2y;
                m3x = a1y-a3y;
                m3y = a3x-a1x;
                a(aoffset+0) = t1x+t2x;
                a(aoffset+1) = t1y+t2y;
                a(aoffset+4) = t1x-t2x;
                a(aoffset+5) = t1y-t2y;
                a(aoffset+2) = m2x+m3x;
                a(aoffset+3) = m2y+m3y;
                a(aoffset+6) = m2x-m3x;
                a(aoffset+7) = m2y-m3y;
                return;
            }
            if( n==5 )
            {
                offs = plan.plan(entryoffset+7);
                c1 = plan.precomputed(offs+0);
                c2 = plan.precomputed(offs+1);
                c3 = plan.precomputed(offs+2);
                c4 = plan.precomputed(offs+3);
                c5 = plan.precomputed(offs+4);
                t1x = a(aoffset+2)+a(aoffset+8);
                t1y = a(aoffset+3)+a(aoffset+9);
                t2x = a(aoffset+4)+a(aoffset+6);
                t2y = a(aoffset+5)+a(aoffset+7);
                t3x = a(aoffset+2)-a(aoffset+8);
                t3y = a(aoffset+3)-a(aoffset+9);
                t4x = a(aoffset+6)-a(aoffset+4);
                t4y = a(aoffset+7)-a(aoffset+5);
                t5x = t1x+t2x;
                t5y = t1y+t2y;
                a(aoffset+0) = a(aoffset+0)+t5x;
                a(aoffset+1) = a(aoffset+1)+t5y;
                m1x = c1*t5x;
                m1y = c1*t5y;
                m2x = c2*(t1x-t2x);
                m2y = c2*(t1y-t2y);
                m3x = -c3*(t3y+t4y);
                m3y = c3*(t3x+t4x);
                m4x = -c4*t4y;
                m4y = c4*t4x;
                m5x = -c5*t3y;
                m5y = c5*t3x;
                s3x = m3x-m4x;
                s3y = m3y-m4y;
                s5x = m3x+m5x;
                s5y = m3y+m5y;
                s1x = a(aoffset+0)+m1x;
                s1y = a(aoffset+1)+m1y;
                s2x = s1x+m2x;
                s2y = s1y+m2y;
                s4x = s1x-m2x;
                s4y = s1y-m2y;
                a(aoffset+2) = s2x+s3x;
                a(aoffset+3) = s2y+s3y;
                a(aoffset+4) = s4x+s5x;
                a(aoffset+5) = s4y+s5y;
                a(aoffset+6) = s4x-s5x;
                a(aoffset+7) = s4y-s5y;
                a(aoffset+8) = s2x-s3x;
                a(aoffset+9) = s2y-s3y;
                return;
            }
        }
        if( plan.plan(entryoffset+3)==fhtcodeletplan )
        {
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            n = n1*n2;
            if( n==2 )
            {
                a0x = a(aoffset+0);
                a1x = a(aoffset+1);
                a(aoffset+0) = a0x+a1x;
                a(aoffset+1) = a0x-a1x;
                return;
            }
            if( n==3 )
            {
                offs = plan.plan(entryoffset+7);
                c1 = plan.precomputed(offs+0);
                c2 = plan.precomputed(offs+1);
                a0x = a(aoffset+0);
                a1x = a(aoffset+1);
                a2x = a(aoffset+2);
                t1x = a1x+a2x;
                a0x = a0x+t1x;
                m1x = c1*t1x;
                m2y = c2*(a2x-a1x);
                s1x = a0x+m1x;
                a(aoffset+0) = a0x;
                a(aoffset+1) = s1x-m2y;
                a(aoffset+2) = s1x+m2y;
                return;
            }
            if( n==4 )
            {
                a0x = a(aoffset+0);
                a1x = a(aoffset+1);
                a2x = a(aoffset+2);
                a3x = a(aoffset+3);
                t1x = a0x+a2x;
                t2x = a1x+a3x;
                m2x = a0x-a2x;
                m3y = a3x-a1x;
                a(aoffset+0) = t1x+t2x;
                a(aoffset+1) = m2x-m3y;
                a(aoffset+2) = t1x-t2x;
                a(aoffset+3) = m2x+m3y;
                return;
            }
            if( n==5 )
            {
                offs = plan.plan(entryoffset+7);
                c1 = plan.precomputed(offs+0);
                c2 = plan.precomputed(offs+1);
                c3 = plan.precomputed(offs+2);
                c4 = plan.precomputed(offs+3);
                c5 = plan.precomputed(offs+4);
                t1x = a(aoffset+1)+a(aoffset+4);
                t2x = a(aoffset+2)+a(aoffset+3);
                t3x = a(aoffset+1)-a(aoffset+4);
                t4x = a(aoffset+3)-a(aoffset+2);
                t5x = t1x+t2x;
                v0 = a(aoffset+0)+t5x;
                a(aoffset+0) = v0;
                m2x = c2*(t1x-t2x);
                m3y = c3*(t3x+t4x);
                s3y = m3y-c4*t4x;
                s5y = m3y+c5*t3x;
                s1x = v0+c1*t5x;
                s2x = s1x+m2x;
                s4x = s1x-m2x;
                a(aoffset+1) = s2x-s3y;
                a(aoffset+2) = s4x-s5y;
                a(aoffset+3) = s4x+s5y;
                a(aoffset+4) = s2x+s3y;
                return;
            }
        }
        if( plan.plan(entryoffset+3)==fftbluesteinplan )
        {
            
            //
            // Bluestein plan:
            // 1. multiply by precomputed coefficients
            // 2. make convolution: forward FFT, multiplication by precomputed FFT
            //    and backward FFT. backward FFT is represented as
            //
            //        invfft(x) = fft(x')'/M
            //
            //    for performance reasons reduction of inverse FFT to
            //    forward FFT is merged with multiplication of FFT components
            //    and last stage of Bluestein's transformation.
            // 3. post-multiplication by Bluestein factors
            //
            n = plan.plan(entryoffset+1);
            m = plan.plan(entryoffset+4);
            offs = plan.plan(entryoffset+7);
            for(i=stackptr+2*n; i<=stackptr+2*m-1; i++)
            {
                plan.stackbuf(i) = 0;
            }
            offsp = offs+2*m;
            offsa = aoffset;
            offsb = stackptr;
            for(i=0; i<=n-1; i++)
            {
                bx = plan.precomputed(offsp+0);
                by = plan.precomputed(offsp+1);
                x = a(offsa+0);
                y = a(offsa+1);
                plan.stackbuf(offsb+0) = x*bx-y*(-by);
                plan.stackbuf(offsb+1) = x*(-by)+y*bx;
                offsp = offsp+2;
                offsa = offsa+2;
                offsb = offsb+2;
            }
            ftbaseexecuteplanrec<Precision>(plan.stackbuf, stackptr, plan, plan.plan(entryoffset+5), stackptr+2*2*m);
            offsb = stackptr;
            offsp = offs;
            for(i=0; i<=m-1; i++)
            {
                x = plan.stackbuf(offsb+0);
                y = plan.stackbuf(offsb+1);
                bx = plan.precomputed(offsp+0);
                by = plan.precomputed(offsp+1);
                plan.stackbuf(offsb+0) = x*bx-y*by;
                plan.stackbuf(offsb+1) = -(x*by+y*bx);
                offsb = offsb+2;
                offsp = offsp+2;
            }
            ftbaseexecuteplanrec<Precision>(plan.stackbuf, stackptr, plan, plan.plan(entryoffset+5), stackptr+2*2*m);
            offsb = stackptr;
            offsp = offs+2*m;
            offsa = aoffset;
            for(i=0; i<=n-1; i++)
            {
                x = +plan.stackbuf(offsb+0)/m;
                y = -plan.stackbuf(offsb+1)/m;
                bx = plan.precomputed(offsp+0);
                by = plan.precomputed(offsp+1);
                a(offsa+0) = x*bx-y*(-by);
                a(offsa+1) = x*(-by)+y*bx;
                offsp = offsp+2;
                offsa = offsa+2;
                offsb = offsb+2;
            }
            return;
        }
    }


    /*************************************************************************
    Returns good factorization N=N1*N2.

    Usually N1<=N2 (but not always - small N's may be exception).
    if N1<>1 then N2<>1.

    Factorization is chosen depending on task type and codelets we have.

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasefactorize(int n,
        int tasktype,
        int& n1,
        int& n2)
    {
        int j;


        n1 = 0;
        n2 = 0;
        
        //
        // try to find good codelet
        //
        if( n1*n2!=n )
        {
            for(j=ftbasecodeletrecommended; j>=2; j--)
            {
                if( n%j==0 )
                {
                    n1 = j;
                    n2 = n/j;
                    break;
                }
            }
        }
        
        //
        // try to factorize N
        //
        if( n1*n2!=n )
        {
            for(j=ftbasecodeletrecommended+1; j<=n-1; j++)
            {
                if( n%j==0 )
                {
                    n1 = j;
                    n2 = n/j;
                    break;
                }
            }
        }
        
        //
        // looks like N is prime :(
        //
        if( n1*n2!=n )
        {
            n1 = 1;
            n2 = n;
        }
        
        //
        // normalize
        //
        if( n2==1 && n1!=1 )
        {
            n2 = n1;
            n1 = 1;
        }
    }


    /*************************************************************************
    Is number smooth?

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    bool ftbaseissmooth(int n)
    {
        bool result;
        int i;


        for(i=2; i<=ftbasemaxsmoothfactor; i++)
        {
            while( n%i==0 )
            {
                n = n/i;
            }
        }
        result = n==1;
        return result;
    }


    /*************************************************************************
    Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
    than or equal to max(N,2)

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int ftbasefindsmooth(int n)
    {
        int result;
        int best;


        best = 2;
        while( best<n )
        {
            best = 2*best;
        }
        ftbasefindsmoothrec<Precision>(n, 1, 2, best);
        result = best;
        return result;
    }


    /*************************************************************************
    Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
    greater than or equal to max(N,2)

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int ftbasefindsmootheven(int n)
    {
        int result;
        int best;


        best = 2;
        while( best<n )
        {
            best = 2*best;
        }
        ftbasefindsmoothrec<Precision>(n, 2, 2, best);
        result = best;
        return result;
    }


    /*************************************************************************
    Returns estimate of FLOP count for the FFT.

    It is only an estimate based on operations count for the PERFECT FFT
    and relative inefficiency of the algorithm actually used.

    N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> ftbasegetflopestimate(int n)
    {
        amp::ampf<Precision> result;


        result = ftbaseinefficiencyfactor<Precision>()*(4*n*amp::log<Precision>(amp::ampf<Precision>(n))/amp::log<Precision>(amp::ampf<Precision>(2))-6*n+8);
        return result;
    }


    /*************************************************************************
    Recurrent subroutine for the FFTGeneratePlan:

    PARAMETERS:
        N                   plan size
        IsReal              whether input is real or not.
                            subroutine MUST NOT ignore this flag because real
                            inputs comes with non-initialized imaginary parts,
                            so ignoring this flag will result in corrupted output
        HalfOut             whether full output or only half of it from 0 to
                            floor(N/2) is needed. This flag may be ignored if
                            doing so will simplify calculations
        Plan                plan array
        PlanSize            size of used part (in integers)
        PrecomputedSize     size of precomputed array allocated yet
        PlanArraySize       plan array size (actual)
        TmpMemSize          temporary memory required size
        BluesteinMemSize    temporary memory required size

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasegenerateplanrec(int n,
        int tasktype,
        ftplan<Precision>& plan,
        int& plansize,
        int& precomputedsize,
        int& planarraysize,
        int& tmpmemsize,
        int& stackmemsize,
        int stackptr)
    {
        int k;
        int m;
        int n1;
        int n2;
        int esize;
        int entryoffset;


        
        //
        // prepare
        //
        if( plansize+ftbaseplanentrysize>planarraysize )
        {
            fftarrayresize<Precision>(plan.plan, planarraysize, 8*planarraysize);
        }
        entryoffset = plansize;
        esize = ftbaseplanentrysize;
        plansize = plansize+esize;
        
        //
        // if N=1, generate empty plan and exit
        //
        if( n==1 )
        {
            plan.plan(entryoffset+0) = esize;
            plan.plan(entryoffset+1) = -1;
            plan.plan(entryoffset+2) = -1;
            plan.plan(entryoffset+3) = fftemptyplan;
            plan.plan(entryoffset+4) = -1;
            plan.plan(entryoffset+5) = -1;
            plan.plan(entryoffset+6) = -1;
            plan.plan(entryoffset+7) = -1;
            return;
        }
        
        //
        // generate plans
        //
        ftbasefactorize<Precision>(n, tasktype, n1, n2);
        if( tasktype==ftbasecffttask || tasktype==ftbaserffttask )
        {
            
            //
            // complex FFT plans
            //
            if( n1!=1 )
            {
                
                //
                // Cooley-Tukey plan (real or complex)
                //
                // Note that child plans are COMPLEX
                // (whether plan itself is complex or not).
                //
                tmpmemsize = ap::maxint(tmpmemsize, 2*n1*n2);
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n1;
                plan.plan(entryoffset+2) = n2;
                if( tasktype==ftbasecffttask )
                {
                    plan.plan(entryoffset+3) = fftcooleytukeyplan;
                }
                else
                {
                    plan.plan(entryoffset+3) = fftrealcooleytukeyplan;
                }
                plan.plan(entryoffset+4) = 0;
                plan.plan(entryoffset+5) = plansize;
                ftbasegenerateplanrec<Precision>(n1, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                plan.plan(entryoffset+6) = plansize;
                ftbasegenerateplanrec<Precision>(n2, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                plan.plan(entryoffset+7) = -1;
                return;
            }
            else
            {
                if( n==2 || n==3 || n==4 || n==5 )
                {
                    
                    //
                    // hard-coded plan
                    //
                    plan.plan(entryoffset+0) = esize;
                    plan.plan(entryoffset+1) = n1;
                    plan.plan(entryoffset+2) = n2;
                    plan.plan(entryoffset+3) = fftcodeletplan;
                    plan.plan(entryoffset+4) = 0;
                    plan.plan(entryoffset+5) = -1;
                    plan.plan(entryoffset+6) = -1;
                    plan.plan(entryoffset+7) = precomputedsize;
                    if( n==3 )
                    {
                        precomputedsize = precomputedsize+2;
                    }
                    if( n==5 )
                    {
                        precomputedsize = precomputedsize+5;
                    }
                    return;
                }
                else
                {
                    
                    //
                    // Bluestein's plan
                    //
                    // Select such M that M>=2*N-1, M is composite, and M's
                    // factors are 2, 3, 5
                    //
                    k = 2*n2-1;
                    m = ftbasefindsmooth<Precision>(k);
                    tmpmemsize = ap::maxint(tmpmemsize, 2*m);
                    plan.plan(entryoffset+0) = esize;
                    plan.plan(entryoffset+1) = n2;
                    plan.plan(entryoffset+2) = -1;
                    plan.plan(entryoffset+3) = fftbluesteinplan;
                    plan.plan(entryoffset+4) = m;
                    plan.plan(entryoffset+5) = plansize;
                    stackptr = stackptr+2*2*m;
                    stackmemsize = ap::maxint(stackmemsize, stackptr);
                    ftbasegenerateplanrec<Precision>(m, ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                    stackptr = stackptr-2*2*m;
                    plan.plan(entryoffset+6) = -1;
                    plan.plan(entryoffset+7) = precomputedsize;
                    precomputedsize = precomputedsize+2*m+2*n;
                    return;
                }
            }
        }
        if( tasktype==ftbaserfhttask )
        {
            
            //
            // real FHT plans
            //
            if( n1!=1 )
            {
                
                //
                // Cooley-Tukey plan
                //
                //
                tmpmemsize = ap::maxint(tmpmemsize, 2*n1*n2);
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n1;
                plan.plan(entryoffset+2) = n2;
                plan.plan(entryoffset+3) = fhtcooleytukeyplan;
                plan.plan(entryoffset+4) = 0;
                plan.plan(entryoffset+5) = plansize;
                ftbasegenerateplanrec<Precision>(n1, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                plan.plan(entryoffset+6) = plansize;
                ftbasegenerateplanrec<Precision>(n2, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr);
                plan.plan(entryoffset+7) = -1;
                return;
            }
            else
            {
                
                //
                // N2 plan
                //
                plan.plan(entryoffset+0) = esize;
                plan.plan(entryoffset+1) = n1;
                plan.plan(entryoffset+2) = n2;
                plan.plan(entryoffset+3) = fhtn2plan;
                plan.plan(entryoffset+4) = 0;
                plan.plan(entryoffset+5) = -1;
                plan.plan(entryoffset+6) = -1;
                plan.plan(entryoffset+7) = -1;
                if( n==2 || n==3 || n==4 || n==5 )
                {
                    
                    //
                    // hard-coded plan
                    //
                    plan.plan(entryoffset+0) = esize;
                    plan.plan(entryoffset+1) = n1;
                    plan.plan(entryoffset+2) = n2;
                    plan.plan(entryoffset+3) = fhtcodeletplan;
                    plan.plan(entryoffset+4) = 0;
                    plan.plan(entryoffset+5) = -1;
                    plan.plan(entryoffset+6) = -1;
                    plan.plan(entryoffset+7) = precomputedsize;
                    if( n==3 )
                    {
                        precomputedsize = precomputedsize+2;
                    }
                    if( n==5 )
                    {
                        precomputedsize = precomputedsize+5;
                    }
                    return;
                }
                return;
            }
        }
    }


    /*************************************************************************
    Recurrent subroutine for precomputing FFT plans

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbaseprecomputeplanrec(ftplan<Precision>& plan,
        int entryoffset,
        int stackptr)
    {
        int i;
        int idx;
        int n1;
        int n2;
        int n;
        int m;
        int offs;
        amp::ampf<Precision> v;
        ap::template_1d_array< amp::ampf<Precision> > emptyarray;
        amp::ampf<Precision> bx;
        amp::ampf<Precision> by;


        if( plan.plan(entryoffset+3)==fftcooleytukeyplan || plan.plan(entryoffset+3)==fftrealcooleytukeyplan || plan.plan(entryoffset+3)==fhtcooleytukeyplan )
        {
            ftbaseprecomputeplanrec<Precision>(plan, plan.plan(entryoffset+5), stackptr);
            ftbaseprecomputeplanrec<Precision>(plan, plan.plan(entryoffset+6), stackptr);
            return;
        }
        if( plan.plan(entryoffset+3)==fftcodeletplan || plan.plan(entryoffset+3)==fhtcodeletplan )
        {
            n1 = plan.plan(entryoffset+1);
            n2 = plan.plan(entryoffset+2);
            n = n1*n2;
            if( n==3 )
            {
                offs = plan.plan(entryoffset+7);
                plan.precomputed(offs+0) = amp::cos<Precision>(2*amp::pi<Precision>()/3)-1;
                plan.precomputed(offs+1) = amp::sin<Precision>(2*amp::pi<Precision>()/3);
                return;
            }
            if( n==5 )
            {
                offs = plan.plan(entryoffset+7);
                v = 2*amp::pi<Precision>()/5;
                plan.precomputed(offs+0) = (amp::cos<Precision>(v)+amp::cos<Precision>(2*v))/2-1;
                plan.precomputed(offs+1) = (amp::cos<Precision>(v)-amp::cos<Precision>(2*v))/2;
                plan.precomputed(offs+2) = -amp::sin<Precision>(v);
                plan.precomputed(offs+3) = -(amp::sin<Precision>(v)+amp::sin<Precision>(2*v));
                plan.precomputed(offs+4) = amp::sin<Precision>(v)-amp::sin<Precision>(2*v);
                return;
            }
        }
        if( plan.plan(entryoffset+3)==fftbluesteinplan )
        {
            ftbaseprecomputeplanrec<Precision>(plan, plan.plan(entryoffset+5), stackptr);
            n = plan.plan(entryoffset+1);
            m = plan.plan(entryoffset+4);
            offs = plan.plan(entryoffset+7);
            for(i=0; i<=2*m-1; i++)
            {
                plan.precomputed(offs+i) = 0;
            }
            for(i=0; i<=n-1; i++)
            {
                bx = amp::cos<Precision>(amp::pi<Precision>()*amp::sqr<Precision>(amp::ampf<Precision>(i))/n);
                by = amp::sin<Precision>(amp::pi<Precision>()*amp::sqr<Precision>(amp::ampf<Precision>(i))/n);
                plan.precomputed(offs+2*i+0) = bx;
                plan.precomputed(offs+2*i+1) = by;
                plan.precomputed(offs+2*m+2*i+0) = bx;
                plan.precomputed(offs+2*m+2*i+1) = by;
                if( i>0 )
                {
                    plan.precomputed(offs+2*(m-i)+0) = bx;
                    plan.precomputed(offs+2*(m-i)+1) = by;
                }
            }
            ftbaseexecuteplanrec<Precision>(plan.precomputed, offs, plan, plan.plan(entryoffset+5), stackptr);
            return;
        }
    }


    /*************************************************************************
    Twiddle factors calculation

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ffttwcalc(ap::template_1d_array< amp::ampf<Precision> >& a,
        int aoffset,
        int n1,
        int n2)
    {
        int i;
        int j;
        int n;
        int idx;
        int offs;
        amp::ampf<Precision> x;
        amp::ampf<Precision> y;
        amp::ampf<Precision> twxm1;
        amp::ampf<Precision> twy;
        amp::ampf<Precision> twbasexm1;
        amp::ampf<Precision> twbasey;
        amp::ampf<Precision> twrowxm1;
        amp::ampf<Precision> twrowy;
        amp::ampf<Precision> tmpx;
        amp::ampf<Precision> tmpy;
        amp::ampf<Precision> v;


        n = n1*n2;
        v = -2*amp::pi<Precision>()/n;
        twbasexm1 = -2*amp::sqr<Precision>(amp::sin<Precision>(amp::ampf<Precision>("0.5")*v));
        twbasey = amp::sin<Precision>(v);
        twrowxm1 = 0;
        twrowy = 0;
        for(i=0; i<=n2-1; i++)
        {
            twxm1 = 0;
            twy = 0;
            for(j=0; j<=n1-1; j++)
            {
                idx = i*n1+j;
                offs = aoffset+2*idx;
                x = a(offs+0);
                y = a(offs+1);
                tmpx = x*twxm1-y*twy;
                tmpy = x*twy+y*twxm1;
                a(offs+0) = x+tmpx;
                a(offs+1) = y+tmpy;
                
                //
                // update Tw: Tw(new) = Tw(old)*TwRow
                //
                if( j<n1-1 )
                {
                    if( j%ftbaseupdatetw==0 )
                    {
                        v = -2*amp::pi<Precision>()*i*(j+1)/n;
                        twxm1 = -2*amp::sqr<Precision>(amp::sin<Precision>(amp::ampf<Precision>("0.5")*v));
                        twy = amp::sin<Precision>(v);
                    }
                    else
                    {
                        tmpx = twrowxm1+twxm1*twrowxm1-twy*twrowy;
                        tmpy = twrowy+twxm1*twrowy+twy*twrowxm1;
                        twxm1 = twxm1+tmpx;
                        twy = twy+tmpy;
                    }
                }
            }
            
            //
            // update TwRow: TwRow(new) = TwRow(old)*TwBase
            //
            if( i<n2-1 )
            {
                if( j%ftbaseupdatetw==0 )
                {
                    v = -2*amp::pi<Precision>()*(i+1)/n;
                    twrowxm1 = -2*amp::sqr<Precision>(amp::sin<Precision>(amp::ampf<Precision>("0.5")*v));
                    twrowy = amp::sin<Precision>(v);
                }
                else
                {
                    tmpx = twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
                    tmpy = twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
                    twrowxm1 = twrowxm1+tmpx;
                    twrowy = twrowy+tmpy;
                }
            }
        }
    }


    /*************************************************************************
    Linear transpose: transpose complex matrix stored in 1-dimensional array

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void internalcomplexlintranspose(ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        int astart,
        ap::template_1d_array< amp::ampf<Precision> >& buf)
    {
        ffticltrec<Precision>(a, astart, n, buf, 0, m, m, n);
        amp::vmove(a.getvector(astart, astart+2*m*n-1), buf.getvector(0, 2*m*n-1));
    }


    /*************************************************************************
    Linear transpose: transpose real matrix stored in 1-dimensional array

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void internalreallintranspose(ap::template_1d_array< amp::ampf<Precision> >& a,
        int m,
        int n,
        int astart,
        ap::template_1d_array< amp::ampf<Precision> >& buf)
    {
        fftirltrec<Precision>(a, astart, n, buf, 0, m, m, n);
        amp::vmove(a.getvector(astart, astart+m*n-1), buf.getvector(0, m*n-1));
    }


    /*************************************************************************
    Recurrent subroutine for a InternalComplexLinTranspose

    Write A^T to B, where:
    * A is m*n complex matrix stored in array A as pairs of real/image values,
      beginning from AStart position, with AStride stride
    * B is n*m complex matrix stored in array B as pairs of real/image values,
      beginning from BStart position, with BStride stride
    stride is measured in complex numbers, i.e. in real/image pairs.

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ffticltrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int astart,
        int astride,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        int bstart,
        int bstride,
        int m,
        int n)
    {
        int i;
        int j;
        int idx1;
        int idx2;
        int m2;
        int m1;
        int n1;


        if( m==0 || n==0 )
        {
            return;
        }
        if( ap::maxint(m, n)<=8 )
        {
            m2 = 2*bstride;
            for(i=0; i<=m-1; i++)
            {
                idx1 = bstart+2*i;
                idx2 = astart+2*i*astride;
                for(j=0; j<=n-1; j++)
                {
                    b(idx1+0) = a(idx2+0);
                    b(idx1+1) = a(idx2+1);
                    idx1 = idx1+m2;
                    idx2 = idx2+2;
                }
            }
            return;
        }
        if( n>m )
        {
            
            //
            // New partition:
            //
            // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
            //                                  ( B2 )
            //
            n1 = n/2;
            if( n-n1>=8 && n1%8!=0 )
            {
                n1 = n1+(8-n1%8);
            }
            ap::ap_error::make_assertion(n-n1>0);
            ffticltrec<Precision>(a, astart, astride, b, bstart, bstride, m, n1);
            ffticltrec<Precision>(a, astart+2*n1, astride, b, bstart+2*n1*bstride, bstride, m, n-n1);
        }
        else
        {
            
            //
            // New partition:
            //
            // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
            //                     ( A2 )
            //
            m1 = m/2;
            if( m-m1>=8 && m1%8!=0 )
            {
                m1 = m1+(8-m1%8);
            }
            ap::ap_error::make_assertion(m-m1>0);
            ffticltrec<Precision>(a, astart, astride, b, bstart, bstride, m1, n);
            ffticltrec<Precision>(a, astart+2*m1*astride, astride, b, bstart+2*m1, bstride, m-m1, n);
        }
    }


    /*************************************************************************
    Recurrent subroutine for a InternalRealLinTranspose


      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void fftirltrec(ap::template_1d_array< amp::ampf<Precision> >& a,
        int astart,
        int astride,
        ap::template_1d_array< amp::ampf<Precision> >& b,
        int bstart,
        int bstride,
        int m,
        int n)
    {
        int i;
        int j;
        int idx1;
        int idx2;
        int m1;
        int n1;


        if( m==0 || n==0 )
        {
            return;
        }
        if( ap::maxint(m, n)<=8 )
        {
            for(i=0; i<=m-1; i++)
            {
                idx1 = bstart+i;
                idx2 = astart+i*astride;
                for(j=0; j<=n-1; j++)
                {
                    b(idx1) = a(idx2);
                    idx1 = idx1+bstride;
                    idx2 = idx2+1;
                }
            }
            return;
        }
        if( n>m )
        {
            
            //
            // New partition:
            //
            // "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
            //                                  ( B2 )
            //
            n1 = n/2;
            if( n-n1>=8 && n1%8!=0 )
            {
                n1 = n1+(8-n1%8);
            }
            ap::ap_error::make_assertion(n-n1>0);
            fftirltrec<Precision>(a, astart, astride, b, bstart, bstride, m, n1);
            fftirltrec<Precision>(a, astart+n1, astride, b, bstart+n1*bstride, bstride, m, n-n1);
        }
        else
        {
            
            //
            // New partition:
            //
            // "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
            //                     ( A2 )
            //
            m1 = m/2;
            if( m-m1>=8 && m1%8!=0 )
            {
                m1 = m1+(8-m1%8);
            }
            ap::ap_error::make_assertion(m-m1>0);
            fftirltrec<Precision>(a, astart, astride, b, bstart, bstride, m1, n);
            fftirltrec<Precision>(a, astart+m1*astride, astride, b, bstart+m1, bstride, m-m1, n);
        }
    }


    /*************************************************************************
    recurrent subroutine for FFTFindSmoothRec

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void ftbasefindsmoothrec(int n,
        int seed,
        int leastfactor,
        int& best)
    {
        ap::ap_error::make_assertion(ftbasemaxsmoothfactor<=5);
        if( seed>=n )
        {
            best = ap::minint(best, seed);
            return;
        }
        if( leastfactor<=2 )
        {
            ftbasefindsmoothrec<Precision>(n, seed*2, 2, best);
        }
        if( leastfactor<=3 )
        {
            ftbasefindsmoothrec<Precision>(n, seed*3, 3, best);
        }
        if( leastfactor<=5 )
        {
            ftbasefindsmoothrec<Precision>(n, seed*5, 5, best);
        }
    }


    /*************************************************************************
    Internal subroutine: array resize

      -- ALGLIB --
         Copyright 01.05.2009 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void fftarrayresize(ap::template_1d_array< int >& a,
        int& asize,
        int newasize)
    {
        ap::template_1d_array< int > tmp;
        int i;


        tmp.setlength(asize);
        for(i=0; i<=asize-1; i++)
        {
            tmp(i) = a(i);
        }
        a.setlength(newasize);
        for(i=0; i<=asize-1; i++)
        {
            a(i) = tmp(i);
        }
        asize = newasize;
    }


    /*************************************************************************
    Reference FHT stub
    *************************************************************************/
    template<unsigned int Precision>
    void reffht(ap::template_1d_array< amp::ampf<Precision> >& a,
        int n,
        int offs)
    {
        ap::template_1d_array< amp::ampf<Precision> > buf;
        int i;
        int j;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(n>0);
        buf.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            v = 0;
            for(j=0; j<=n-1; j++)
            {
                v = v+a(offs+j)*(amp::cos<Precision>(2*amp::pi<Precision>()*i*j/n)+amp::sin<Precision>(2*amp::pi<Precision>()*i*j/n));
            }
            buf(i) = v;
        }
        for(i=0; i<=n-1; i++)
        {
            a(offs+i) = buf(i);
        }
    }
} // namespace

#endif
