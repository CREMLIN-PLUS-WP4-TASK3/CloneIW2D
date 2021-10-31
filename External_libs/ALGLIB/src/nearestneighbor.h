/*************************************************************************
Copyright (c) 2010, Sergey Bochkanov (ALGLIB project).

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

#ifndef _nearestneighbor_h
#define _nearestneighbor_h

#include "ap.h"
#include "amp.h"
#include "tsort.h"
namespace nearestneighbor
{
    template<unsigned int Precision>
    class kdtree
    {
    public:
        int n;
        int nx;
        int ny;
        int normtype;
        int distmatrixtype;
        ap::template_2d_array< amp::ampf<Precision> > xy;
        ap::template_1d_array< int > tags;
        ap::template_1d_array< amp::ampf<Precision> > boxmin;
        ap::template_1d_array< amp::ampf<Precision> > boxmax;
        ap::template_1d_array< amp::ampf<Precision> > curboxmin;
        ap::template_1d_array< amp::ampf<Precision> > curboxmax;
        amp::ampf<Precision> curdist;
        ap::template_1d_array< int > nodes;
        ap::template_1d_array< amp::ampf<Precision> > splits;
        ap::template_1d_array< amp::ampf<Precision> > x;
        int kneeded;
        amp::ampf<Precision> rneeded;
        bool selfmatch;
        amp::ampf<Precision> approxf;
        int kcur;
        ap::template_1d_array< int > idx;
        ap::template_1d_array< amp::ampf<Precision> > r;
        ap::template_1d_array< amp::ampf<Precision> > buf;
        int debugcounter;
    };




    template<unsigned int Precision>
    void kdtreebuild(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int nx,
        int ny,
        int normtype,
        kdtree<Precision>& kdt);
    template<unsigned int Precision>
    void kdtreebuildtagged(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const ap::template_1d_array< int >& tags,
        int n,
        int nx,
        int ny,
        int normtype,
        kdtree<Precision>& kdt);
    template<unsigned int Precision>
    int kdtreequeryknn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int k,
        bool selfmatch);
    template<unsigned int Precision>
    int kdtreequeryrnn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision> r,
        bool selfmatch);
    template<unsigned int Precision>
    int kdtreequeryaknn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int k,
        bool selfmatch,
        amp::ampf<Precision> eps);
    template<unsigned int Precision>
    void kdtreequeryresultsx(const kdtree<Precision>& kdt,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int& k);
    template<unsigned int Precision>
    void kdtreequeryresultsxy(const kdtree<Precision>& kdt,
        ap::template_2d_array< amp::ampf<Precision> >& xy,
        int& k);
    template<unsigned int Precision>
    void kdtreequeryresultstags(const kdtree<Precision>& kdt,
        ap::template_1d_array< int >& tags,
        int& k);
    template<unsigned int Precision>
    void kdtreequeryresultsdistances(const kdtree<Precision>& kdt,
        ap::template_1d_array< amp::ampf<Precision> >& r,
        int& k);
    template<unsigned int Precision>
    void kdtreesplit(kdtree<Precision>& kdt,
        int i1,
        int i2,
        int d,
        amp::ampf<Precision> s,
        int& i3);
    template<unsigned int Precision>
    void kdtreegeneratetreerec(kdtree<Precision>& kdt,
        int& nodesoffs,
        int& splitsoffs,
        int i1,
        int i2,
        int maxleafsize);
    template<unsigned int Precision>
    void kdtreequerynnrec(kdtree<Precision>& kdt,
        int offs);
    template<unsigned int Precision>
    void kdtreeinitbox(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x);
    template<unsigned int Precision>
    amp::ampf<Precision> vrootfreenorm(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        int normtype);
    template<unsigned int Precision>
    amp::ampf<Precision> vrootfreecomponentnorm(amp::ampf<Precision> x,
        int normtype);
    template<unsigned int Precision>
    amp::ampf<Precision> vrangedist(amp::ampf<Precision> x,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b);


    static const int splitnodesize = 6;


    /*************************************************************************
    KD-tree creation

    This subroutine creates KD-tree from set of X-values and optional Y-values

    INPUT PARAMETERS
        XY      -   dataset, array[0..N-1,0..NX+NY-1].
                    one row corresponds to one point.
                    first NX columns contain X-values, next NY (NY may be zero)
                    columns may contain associated Y-values
        N       -   number of points, N>=1
        NX      -   space dimension, NX>=1.
        NY      -   number of optional Y-values, NY>=0.
        NormType-   norm type:
                    * 0 denotes infinity-norm
                    * 1 denotes 1-norm
                    * 2 denotes 2-norm (Euclidean norm)
                    
    OUTPUT PARAMETERS
        KDT     -   KD-tree
        
        
    NOTES

    1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
       requirements.
    2. Although KD-trees may be used with any combination of N  and  NX,  they
       are more efficient than brute-force search only when N >> 4^NX. So they
       are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
       inefficient case, because  simple  binary  search  (without  additional
       structures) is much more efficient in such tasks than KD-trees.

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreebuild(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        int n,
        int nx,
        int ny,
        int normtype,
        kdtree<Precision>& kdt)
    {
        ap::template_1d_array< int > tags;
        int i;


        ap::ap_error::make_assertion(n>=1);
        ap::ap_error::make_assertion(nx>=1);
        ap::ap_error::make_assertion(ny>=0);
        ap::ap_error::make_assertion(normtype>=0 && normtype<=2);
        tags.setlength(n);
        for(i=0; i<=n-1; i++)
        {
            tags(i) = 0;
        }
        kdtreebuildtagged<Precision>(xy, tags, n, nx, ny, normtype, kdt);
    }


    /*************************************************************************
    KD-tree creation

    This  subroutine  creates  KD-tree  from set of X-values, integer tags and
    optional Y-values

    INPUT PARAMETERS
        XY      -   dataset, array[0..N-1,0..NX+NY-1].
                    one row corresponds to one point.
                    first NX columns contain X-values, next NY (NY may be zero)
                    columns may contain associated Y-values
        Tags    -   tags, array[0..N-1], contains integer tags associated
                    with points.
        N       -   number of points, N>=1
        NX      -   space dimension, NX>=1.
        NY      -   number of optional Y-values, NY>=0.
        NormType-   norm type:
                    * 0 denotes infinity-norm
                    * 1 denotes 1-norm
                    * 2 denotes 2-norm (Euclidean norm)

    OUTPUT PARAMETERS
        KDT     -   KD-tree

    NOTES

    1. KD-tree  creation  have O(N*logN) complexity and O(N*(2*NX+NY))  memory
       requirements.
    2. Although KD-trees may be used with any combination of N  and  NX,  they
       are more efficient than brute-force search only when N >> 4^NX. So they
       are most useful in low-dimensional tasks (NX=2, NX=3). NX=1  is another
       inefficient case, because  simple  binary  search  (without  additional
       structures) is much more efficient in such tasks than KD-trees.

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreebuildtagged(const ap::template_2d_array< amp::ampf<Precision> >& xy,
        const ap::template_1d_array< int >& tags,
        int n,
        int nx,
        int ny,
        int normtype,
        kdtree<Precision>& kdt)
    {
        int i;
        int j;
        int maxnodes;
        int nodesoffs;
        int splitsoffs;


        ap::ap_error::make_assertion(n>=1);
        ap::ap_error::make_assertion(nx>=1);
        ap::ap_error::make_assertion(ny>=0);
        ap::ap_error::make_assertion(normtype>=0 && normtype<=2);
        
        //
        // initialize
        //
        kdt.n = n;
        kdt.nx = nx;
        kdt.ny = ny;
        kdt.normtype = normtype;
        kdt.distmatrixtype = 0;
        kdt.xy.setlength(n, 2*nx+ny);
        kdt.tags.setlength(n);
        kdt.idx.setlength(n);
        kdt.r.setlength(n);
        kdt.x.setlength(nx);
        kdt.buf.setlength(ap::maxint(n, nx));
        
        //
        // Initial fill
        //
        for(i=0; i<=n-1; i++)
        {
            amp::vmove(kdt.xy.getrow(i, 0, nx-1), xy.getrow(i, 0, nx-1));
            amp::vmove(kdt.xy.getrow(i, nx, 2*nx+ny-1), xy.getrow(i, 0, nx+ny-1));
            kdt.tags(i) = tags(i);
        }
        
        //
        // Determine bounding box
        //
        kdt.boxmin.setlength(nx);
        kdt.boxmax.setlength(nx);
        kdt.curboxmin.setlength(nx);
        kdt.curboxmax.setlength(nx);
        amp::vmove(kdt.boxmin.getvector(0, nx-1), kdt.xy.getrow(0, 0, nx-1));
        amp::vmove(kdt.boxmax.getvector(0, nx-1), kdt.xy.getrow(0, 0, nx-1));
        for(i=1; i<=n-1; i++)
        {
            for(j=0; j<=nx-1; j++)
            {
                kdt.boxmin(j) = amp::minimum<Precision>(kdt.boxmin(j), kdt.xy(i,j));
                kdt.boxmax(j) = amp::maximum<Precision>(kdt.boxmax(j), kdt.xy(i,j));
            }
        }
        
        //
        // prepare tree structure
        // * MaxNodes=N because we guarantee no trivial splits, i.e.
        //   every split will generate two non-empty boxes
        //
        maxnodes = n;
        kdt.nodes.setlength(splitnodesize*2*maxnodes);
        kdt.splits.setlength(2*maxnodes);
        nodesoffs = 0;
        splitsoffs = 0;
        amp::vmove(kdt.curboxmin.getvector(0, nx-1), kdt.boxmin.getvector(0, nx-1));
        amp::vmove(kdt.curboxmax.getvector(0, nx-1), kdt.boxmax.getvector(0, nx-1));
        kdtreegeneratetreerec<Precision>(kdt, nodesoffs, splitsoffs, 0, n, 8);
        
        //
        // Set current query size to 0
        //
        kdt.kcur = 0;
    }


    /*************************************************************************
    K-NN query: K nearest neighbors

    INPUT PARAMETERS
        KDT         -   KD-tree
        X           -   point, array[0..NX-1].
        K           -   number of neighbors to return, K>=1
        SelfMatch   -   whether self-matches are allowed:
                        * if True, nearest neighbor may be the point itself
                          (if it exists in original dataset)
                        * if False, then only points with non-zero distance
                          are returned

    RESULT
        number of actual neighbors found (either K or N, if K>N).

    This  subroutine  performs  query  and  stores  its result in the internal
    structures of the KD-tree. You can use  following  subroutines  to  obtain
    these results:
    * KDTreeQueryResultsX() to get X-values
    * KDTreeQueryResultsXY() to get X- and Y-values
    * KDTreeQueryResultsTags() to get tag values
    * KDTreeQueryResultsDistances() to get distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int kdtreequeryknn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int k,
        bool selfmatch)
    {
        int result;


        result = kdtreequeryaknn<Precision>(kdt, x, k, selfmatch, amp::ampf<Precision>("0.0"));
        return result;
    }


    /*************************************************************************
    R-NN query: all points within R-sphere centered at X

    INPUT PARAMETERS
        KDT         -   KD-tree
        X           -   point, array[0..NX-1].
        R           -   radius of sphere (in corresponding norm), R>0
        SelfMatch   -   whether self-matches are allowed:
                        * if True, nearest neighbor may be the point itself
                          (if it exists in original dataset)
                        * if False, then only points with non-zero distance
                          are returned

    RESULT
        number of neighbors found, >=0

    This  subroutine  performs  query  and  stores  its result in the internal
    structures of the KD-tree. You can use  following  subroutines  to  obtain
    actual results:
    * KDTreeQueryResultsX() to get X-values
    * KDTreeQueryResultsXY() to get X- and Y-values
    * KDTreeQueryResultsTags() to get tag values
    * KDTreeQueryResultsDistances() to get distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int kdtreequeryrnn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        amp::ampf<Precision> r,
        bool selfmatch)
    {
        int result;
        int i;
        int j;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vmin;
        amp::ampf<Precision> vmax;


        ap::ap_error::make_assertion(r>0);
        
        //
        // Prepare parameters
        //
        kdt.kneeded = 0;
        if( kdt.normtype!=2 )
        {
            kdt.rneeded = r;
        }
        else
        {
            kdt.rneeded = amp::sqr<Precision>(r);
        }
        kdt.selfmatch = selfmatch;
        kdt.approxf = 1;
        kdt.kcur = 0;
        
        //
        // calculate distance from point to current bounding box
        //
        kdtreeinitbox<Precision>(kdt, x);
        
        //
        // call recursive search
        // results are returned as heap
        //
        kdtreequerynnrec<Precision>(kdt, 0);
        
        //
        // pop from heap to generate ordered representation
        //
        // last element is non pop'ed because it is already in
        // its place
        //
        result = kdt.kcur;
        j = kdt.kcur;
        for(i=kdt.kcur; i>=2; i--)
        {
            tsort::tagheappopi<Precision>(kdt.r, kdt.idx, j);
        }
        return result;
    }


    /*************************************************************************
    K-NN query: approximate K nearest neighbors

    INPUT PARAMETERS
        KDT         -   KD-tree
        X           -   point, array[0..NX-1].
        K           -   number of neighbors to return, K>=1
        SelfMatch   -   whether self-matches are allowed:
                        * if True, nearest neighbor may be the point itself
                          (if it exists in original dataset)
                        * if False, then only points with non-zero distance
                          are returned
        Eps         -   approximation factor, Eps>=0. eps-approximate  nearest
                        neighbor  is  a  neighbor  whose distance from X is at
                        most (1+eps) times distance of true nearest neighbor.

    RESULT
        number of actual neighbors found (either K or N, if K>N).
        
    NOTES
        significant performance gain may be achieved only when Eps  is  is  on
        the order of magnitude of 1 or larger.

    This  subroutine  performs  query  and  stores  its result in the internal
    structures of the KD-tree. You can use  following  subroutines  to  obtain
    these results:
    * KDTreeQueryResultsX() to get X-values
    * KDTreeQueryResultsXY() to get X- and Y-values
    * KDTreeQueryResultsTags() to get tag values
    * KDTreeQueryResultsDistances() to get distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    int kdtreequeryaknn(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x,
        int k,
        bool selfmatch,
        amp::ampf<Precision> eps)
    {
        int result;
        int i;
        int j;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vmin;
        amp::ampf<Precision> vmax;


        ap::ap_error::make_assertion(k>0);
        ap::ap_error::make_assertion(eps>=0);
        
        //
        // Prepare parameters
        //
        k = ap::minint(k, kdt.n);
        kdt.kneeded = k;
        kdt.rneeded = 0;
        kdt.selfmatch = selfmatch;
        if( kdt.normtype==2 )
        {
            kdt.approxf = 1/amp::sqr<Precision>(1+eps);
        }
        else
        {
            kdt.approxf = 1/(1+eps);
        }
        kdt.kcur = 0;
        
        //
        // calculate distance from point to current bounding box
        //
        kdtreeinitbox<Precision>(kdt, x);
        
        //
        // call recursive search
        // results are returned as heap
        //
        kdtreequerynnrec<Precision>(kdt, 0);
        
        //
        // pop from heap to generate ordered representation
        //
        // last element is non pop'ed because it is already in
        // its place
        //
        result = kdt.kcur;
        j = kdt.kcur;
        for(i=kdt.kcur; i>=2; i--)
        {
            tsort::tagheappopi<Precision>(kdt.r, kdt.idx, j);
        }
        return result;
    }


    /*************************************************************************
    X-values from last query

    INPUT PARAMETERS
        KDT     -   KD-tree
        X       -   pre-allocated array, at least K rows, at least NX columns
        
    OUTPUT PARAMETERS
        X       -   K rows are filled with X-values
        K       -   number of points

    NOTE
        points are ordered by distance from the query point (first = closest)

    SEE ALSO
    * KDTreeQueryResultsXY()            X- and Y-values
    * KDTreeQueryResultsTags()          tag values
    * KDTreeQueryResultsDistances()     distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreequeryresultsx(const kdtree<Precision>& kdt,
        ap::template_2d_array< amp::ampf<Precision> >& x,
        int& k)
    {
        int i;


        k = kdt.kcur;
        for(i=0; i<=k-1; i++)
        {
            amp::vmove(x.getrow(i, 0, kdt.nx-1), kdt.xy.getrow(kdt.idx(i), kdt.nx, 2*kdt.nx-1));
        }
    }


    /*************************************************************************
    X- and Y-values from last query

    INPUT PARAMETERS
        KDT     -   KD-tree
        XY      -   pre-allocated array, at least K rows, at least NX+NY columns

    OUTPUT PARAMETERS
        X       -   K rows are filled with points: first NX columns with
                    X-values, next NY columns - with Y-values.
        K       -   number of points

    NOTE
        points are ordered by distance from the query point (first = closest)

    SEE ALSO
    * KDTreeQueryResultsX()             X-values
    * KDTreeQueryResultsTags()          tag values
    * KDTreeQueryResultsDistances()     distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreequeryresultsxy(const kdtree<Precision>& kdt,
        ap::template_2d_array< amp::ampf<Precision> >& xy,
        int& k)
    {
        int i;


        k = kdt.kcur;
        for(i=0; i<=k-1; i++)
        {
            amp::vmove(xy.getrow(i, 0, kdt.nx+kdt.ny-1), kdt.xy.getrow(kdt.idx(i), kdt.nx, 2*kdt.nx+kdt.ny-1));
        }
    }


    /*************************************************************************
    point tags from last query

    INPUT PARAMETERS
        KDT     -   KD-tree
        Tags    -   pre-allocated array, at least K elements

    OUTPUT PARAMETERS
        Tags    -   first K elements are filled with tags associated with points,
                    or, when no tags were supplied, with zeros
        K       -   number of points

    NOTE
        points are ordered by distance from the query point (first = closest)

    SEE ALSO
    * KDTreeQueryResultsX()             X-values
    * KDTreeQueryResultsXY()            X- and Y-values
    * KDTreeQueryResultsDistances()     distances

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreequeryresultstags(const kdtree<Precision>& kdt,
        ap::template_1d_array< int >& tags,
        int& k)
    {
        int i;


        k = kdt.kcur;
        for(i=0; i<=k-1; i++)
        {
            tags(i) = kdt.tags(kdt.idx(i));
        }
    }


    /*************************************************************************
    Distances from last query

    INPUT PARAMETERS
        KDT     -   KD-tree
        R       -   pre-allocated array, at least K elements

    OUTPUT PARAMETERS
        R       -   first K elements are filled with distances
                    (in corresponding norm)
        K       -   number of points

    NOTE
        points are ordered by distance from the query point (first = closest)

    SEE ALSO
    * KDTreeQueryResultsX()             X-values
    * KDTreeQueryResultsXY()            X- and Y-values
    * KDTreeQueryResultsTags()          tag values

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreequeryresultsdistances(const kdtree<Precision>& kdt,
        ap::template_1d_array< amp::ampf<Precision> >& r,
        int& k)
    {
        int i;


        k = kdt.kcur;
        
        //
        // unload norms
        //
        // Abs() call is used to handle cases with negative norms
        // (generated during KFN requests)
        //
        if( kdt.normtype==0 )
        {
            for(i=0; i<=k-1; i++)
            {
                r(i) = amp::abs<Precision>(kdt.r(i));
            }
        }
        if( kdt.normtype==1 )
        {
            for(i=0; i<=k-1; i++)
            {
                r(i) = amp::abs<Precision>(kdt.r(i));
            }
        }
        if( kdt.normtype==2 )
        {
            for(i=0; i<=k-1; i++)
            {
                r(i) = amp::sqrt<Precision>(amp::abs<Precision>(kdt.r(i)));
            }
        }
    }


    /*************************************************************************
    Rearranges nodes [I1,I2) using partition in D-th dimension with S as threshold.
    Returns split position I3: [I1,I3) and [I3,I2) are created as result.

    This subroutine doesn't create tree structures, just rearranges nodes.
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreesplit(kdtree<Precision>& kdt,
        int i1,
        int i2,
        int d,
        amp::ampf<Precision> s,
        int& i3)
    {
        int i;
        int j;
        int ileft;
        int iright;
        amp::ampf<Precision> v;


        
        //
        // split XY/Tags in two parts:
        // * [ILeft,IRight] is non-processed part of XY/Tags
        //
        // After cycle is done, we have Ileft=IRight. We deal with
        // this element separately.
        //
        // After this, [I1,ILeft) contains left part, and [ILeft,I2)
        // contains right part.
        //
        ileft = i1;
        iright = i2-1;
        while( ileft<iright )
        {
            if( kdt.xy(ileft,d)<=s )
            {
                
                //
                // XY[ILeft] is on its place.
                // Advance ILeft.
                //
                ileft = ileft+1;
            }
            else
            {
                
                //
                // XY[ILeft,..] must be at IRight.
                // Swap and advance IRight.
                //
                for(i=0; i<=2*kdt.nx+kdt.ny-1; i++)
                {
                    v = kdt.xy(ileft,i);
                    kdt.xy(ileft,i) = kdt.xy(iright,i);
                    kdt.xy(iright,i) = v;
                }
                j = kdt.tags(ileft);
                kdt.tags(ileft) = kdt.tags(iright);
                kdt.tags(iright) = j;
                iright = iright-1;
            }
        }
        if( kdt.xy(ileft,d)<=s )
        {
            ileft = ileft+1;
        }
        else
        {
            iright = iright-1;
        }
        i3 = ileft;
    }


    /*************************************************************************
    Recursive kd-tree generation subroutine.

    PARAMETERS
        KDT         tree
        NodesOffs   unused part of Nodes[] which must be filled by tree
        SplitsOffs  unused part of Splits[]
        I1, I2      points from [I1,I2) are processed
        
    NodesOffs[] and SplitsOffs[] must be large enough.

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreegeneratetreerec(kdtree<Precision>& kdt,
        int& nodesoffs,
        int& splitsoffs,
        int i1,
        int i2,
        int maxleafsize)
    {
        int n;
        int nx;
        int ny;
        int i;
        int j;
        int oldoffs;
        int i3;
        int cntless;
        int cntgreater;
        amp::ampf<Precision> minv;
        amp::ampf<Precision> maxv;
        int minidx;
        int maxidx;
        int d;
        amp::ampf<Precision> ds;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;


        ap::ap_error::make_assertion(i2>i1);
        
        //
        // Generate leaf if needed
        //
        if( i2-i1<=maxleafsize )
        {
            kdt.nodes(nodesoffs+0) = i2-i1;
            kdt.nodes(nodesoffs+1) = i1;
            nodesoffs = nodesoffs+2;
            return;
        }
        
        //
        // Load values for easier access
        //
        nx = kdt.nx;
        ny = kdt.ny;
        
        //
        // select dimension to split:
        // * D is a dimension number
        //
        d = 0;
        ds = kdt.curboxmax(0)-kdt.curboxmin(0);
        for(i=1; i<=nx-1; i++)
        {
            v = kdt.curboxmax(i)-kdt.curboxmin(i);
            if( v>ds )
            {
                ds = v;
                d = i;
            }
        }
        
        //
        // Select split position S using sliding midpoint rule,
        // rearrange points into [I1,I3) and [I3,I2)
        //
        s = kdt.curboxmin(d)+amp::ampf<Precision>("0.5")*ds;
        amp::vmove(kdt.buf.getvector(0, i2-i1-1), kdt.xy.getcolumn(d, i1, i2-1));
        n = i2-i1;
        cntless = 0;
        cntgreater = 0;
        minv = kdt.buf(0);
        maxv = kdt.buf(0);
        minidx = i1;
        maxidx = i1;
        for(i=0; i<=n-1; i++)
        {
            v = kdt.buf(i);
            if( v<minv )
            {
                minv = v;
                minidx = i1+i;
            }
            if( v>maxv )
            {
                maxv = v;
                maxidx = i1+i;
            }
            if( v<s )
            {
                cntless = cntless+1;
            }
            if( v>s )
            {
                cntgreater = cntgreater+1;
            }
        }
        if( cntless>0 && cntgreater>0 )
        {
            
            //
            // normal midpoint split
            //
            kdtreesplit<Precision>(kdt, i1, i2, d, s, i3);
        }
        else
        {
            
            //
            // sliding midpoint
            //
            if( cntless==0 )
            {
                
                //
                // 1. move split to MinV,
                // 2. place one point to the left bin (move to I1),
                //    others - to the right bin
                //
                s = minv;
                if( minidx!=i1 )
                {
                    for(i=0; i<=2*kdt.nx+kdt.ny-1; i++)
                    {
                        v = kdt.xy(minidx,i);
                        kdt.xy(minidx,i) = kdt.xy(i1,i);
                        kdt.xy(i1,i) = v;
                    }
                    j = kdt.tags(minidx);
                    kdt.tags(minidx) = kdt.tags(i1);
                    kdt.tags(i1) = j;
                }
                i3 = i1+1;
            }
            else
            {
                
                //
                // 1. move split to MaxV,
                // 2. place one point to the right bin (move to I2-1),
                //    others - to the left bin
                //
                s = maxv;
                if( maxidx!=i2-1 )
                {
                    for(i=0; i<=2*kdt.nx+kdt.ny-1; i++)
                    {
                        v = kdt.xy(maxidx,i);
                        kdt.xy(maxidx,i) = kdt.xy(i2-1,i);
                        kdt.xy(i2-1,i) = v;
                    }
                    j = kdt.tags(maxidx);
                    kdt.tags(maxidx) = kdt.tags(i2-1);
                    kdt.tags(i2-1) = j;
                }
                i3 = i2-1;
            }
        }
        
        //
        // Generate 'split' node
        //
        kdt.nodes(nodesoffs+0) = 0;
        kdt.nodes(nodesoffs+1) = d;
        kdt.nodes(nodesoffs+2) = splitsoffs;
        kdt.splits(splitsoffs+0) = s;
        oldoffs = nodesoffs;
        nodesoffs = nodesoffs+splitnodesize;
        splitsoffs = splitsoffs+1;
        
        //
        // Recirsive generation:
        // * update CurBox
        // * call subroutine
        // * restore CurBox
        //
        kdt.nodes(oldoffs+3) = nodesoffs;
        v = kdt.curboxmax(d);
        kdt.curboxmax(d) = s;
        kdtreegeneratetreerec<Precision>(kdt, nodesoffs, splitsoffs, i1, i3, maxleafsize);
        kdt.curboxmax(d) = v;
        kdt.nodes(oldoffs+4) = nodesoffs;
        v = kdt.curboxmin(d);
        kdt.curboxmin(d) = s;
        kdtreegeneratetreerec<Precision>(kdt, nodesoffs, splitsoffs, i3, i2, maxleafsize);
        kdt.curboxmin(d) = v;
    }


    /*************************************************************************
    Recursive subroutine for NN queries.

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreequerynnrec(kdtree<Precision>& kdt,
        int offs)
    {
        amp::ampf<Precision> ptdist;
        int i;
        int j;
        int k;
        int ti;
        int nx;
        int i1;
        int i2;
        int k1;
        int k2;
        amp::ampf<Precision> r1;
        amp::ampf<Precision> r2;
        int d;
        amp::ampf<Precision> s;
        amp::ampf<Precision> v;
        amp::ampf<Precision> t1;
        int childbestoffs;
        int childworstoffs;
        int childoffs;
        amp::ampf<Precision> prevdist;
        bool todive;
        bool bestisleft;
        bool updatemin;


        
        //
        // Leaf node.
        // Process points.
        //
        if( kdt.nodes(offs)>0 )
        {
            i1 = kdt.nodes(offs+1);
            i2 = i1+kdt.nodes(offs);
            for(i=i1; i<=i2-1; i++)
            {
                
                //
                // Calculate distance
                //
                ptdist = 0;
                nx = kdt.nx;
                if( kdt.normtype==0 )
                {
                    for(j=0; j<=nx-1; j++)
                    {
                        ptdist = amp::maximum<Precision>(ptdist, amp::abs<Precision>(kdt.xy(i,j)-kdt.x(j)));
                    }
                }
                if( kdt.normtype==1 )
                {
                    for(j=0; j<=nx-1; j++)
                    {
                        ptdist = ptdist+amp::abs<Precision>(kdt.xy(i,j)-kdt.x(j));
                    }
                }
                if( kdt.normtype==2 )
                {
                    for(j=0; j<=nx-1; j++)
                    {
                        ptdist = ptdist+amp::sqr<Precision>(kdt.xy(i,j)-kdt.x(j));
                    }
                }
                
                //
                // Skip points with zero distance if self-matches are turned off
                //
                if( ptdist==0 && !kdt.selfmatch )
                {
                    continue;
                }
                
                //
                // We CAN'T process point if R-criterion isn't satisfied,
                // i.e. (RNeeded<>0) AND (PtDist>R).
                //
                if( kdt.rneeded==0 || ptdist<=kdt.rneeded )
                {
                    
                    //
                    // R-criterion is satisfied, we must either:
                    // * replace worst point, if (KNeeded<>0) AND (KCur=KNeeded)
                    //   (or skip, if worst point is better)
                    // * add point without replacement otherwise
                    //
                    if( kdt.kcur<kdt.kneeded || kdt.kneeded==0 )
                    {
                        
                        //
                        // add current point to heap without replacement
                        //
                        tsort::tagheappushi<Precision>(kdt.r, kdt.idx, kdt.kcur, ptdist, i);
                    }
                    else
                    {
                        
                        //
                        // New points are added or not, depending on their distance.
                        // If added, they replace element at the top of the heap
                        //
                        if( ptdist<kdt.r(0) )
                        {
                            if( kdt.kneeded==1 )
                            {
                                kdt.idx(0) = i;
                                kdt.r(0) = ptdist;
                            }
                            else
                            {
                                tsort::tagheapreplacetopi<Precision>(kdt.r, kdt.idx, kdt.kneeded, ptdist, i);
                            }
                        }
                    }
                }
            }
            return;
        }
        
        //
        // Simple split
        //
        if( kdt.nodes(offs)==0 )
        {
            
            //
            // Load:
            // * D  dimension to split
            // * S  split position
            //
            d = kdt.nodes(offs+1);
            s = kdt.splits(kdt.nodes(offs+2));
            
            //
            // Calculate:
            // * ChildBestOffs      child box with best chances
            // * ChildWorstOffs     child box with worst chances
            //
            if( kdt.x(d)<=s )
            {
                childbestoffs = kdt.nodes(offs+3);
                childworstoffs = kdt.nodes(offs+4);
                bestisleft = true;
            }
            else
            {
                childbestoffs = kdt.nodes(offs+4);
                childworstoffs = kdt.nodes(offs+3);
                bestisleft = false;
            }
            
            //
            // Navigate through childs
            //
            for(i=0; i<=1; i++)
            {
                
                //
                // Select child to process:
                // * ChildOffs      current child offset in Nodes[]
                // * UpdateMin      whether minimum or maximum value
                //                  of bounding box is changed on update
                //
                if( i==0 )
                {
                    childoffs = childbestoffs;
                    updatemin = !bestisleft;
                }
                else
                {
                    updatemin = bestisleft;
                    childoffs = childworstoffs;
                }
                
                //
                // Update bounding box and current distance
                //
                if( updatemin )
                {
                    prevdist = kdt.curdist;
                    t1 = kdt.x(d);
                    v = kdt.curboxmin(d);
                    if( t1<=s )
                    {
                        if( kdt.normtype==0 )
                        {
                            kdt.curdist = amp::maximum<Precision>(kdt.curdist, s-t1);
                        }
                        if( kdt.normtype==1 )
                        {
                            kdt.curdist = kdt.curdist-amp::maximum<Precision>(v-t1, amp::ampf<Precision>(0))+s-t1;
                        }
                        if( kdt.normtype==2 )
                        {
                            kdt.curdist = kdt.curdist-amp::sqr<Precision>(amp::maximum<Precision>(v-t1, amp::ampf<Precision>(0)))+amp::sqr<Precision>(s-t1);
                        }
                    }
                    kdt.curboxmin(d) = s;
                }
                else
                {
                    prevdist = kdt.curdist;
                    t1 = kdt.x(d);
                    v = kdt.curboxmax(d);
                    if( t1>=s )
                    {
                        if( kdt.normtype==0 )
                        {
                            kdt.curdist = amp::maximum<Precision>(kdt.curdist, t1-s);
                        }
                        if( kdt.normtype==1 )
                        {
                            kdt.curdist = kdt.curdist-amp::maximum<Precision>(t1-v, amp::ampf<Precision>(0))+t1-s;
                        }
                        if( kdt.normtype==2 )
                        {
                            kdt.curdist = kdt.curdist-amp::sqr<Precision>(amp::maximum<Precision>(t1-v, amp::ampf<Precision>(0)))+amp::sqr<Precision>(t1-s);
                        }
                    }
                    kdt.curboxmax(d) = s;
                }
                
                //
                // Decide: to dive into cell or not to dive
                //
                if( kdt.rneeded!=0 && kdt.curdist>kdt.rneeded )
                {
                    todive = false;
                }
                else
                {
                    if( kdt.kcur<kdt.kneeded || kdt.kneeded==0 )
                    {
                        
                        //
                        // KCur<KNeeded (i.e. not all points are found)
                        //
                        todive = true;
                    }
                    else
                    {
                        
                        //
                        // KCur=KNeeded, decide to dive or not to dive
                        // using point position relative to bounding box.
                        //
                        todive = kdt.curdist<=kdt.r(0)*kdt.approxf;
                    }
                }
                if( todive )
                {
                    kdtreequerynnrec<Precision>(kdt, childoffs);
                }
                
                //
                // Restore bounding box and distance
                //
                if( updatemin )
                {
                    kdt.curboxmin(d) = v;
                }
                else
                {
                    kdt.curboxmax(d) = v;
                }
                kdt.curdist = prevdist;
            }
            return;
        }
    }


    /*************************************************************************
    Copies X[] to KDT.X[]
    Loads distance from X[] to bounding box.
    Initializes CurBox[].

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    void kdtreeinitbox(kdtree<Precision>& kdt,
        const ap::template_1d_array< amp::ampf<Precision> >& x)
    {
        int i;
        amp::ampf<Precision> vx;
        amp::ampf<Precision> vmin;
        amp::ampf<Precision> vmax;


        
        //
        // calculate distance from point to current bounding box
        //
        kdt.curdist = 0;
        if( kdt.normtype==0 )
        {
            for(i=0; i<=kdt.nx-1; i++)
            {
                vx = x(i);
                vmin = kdt.boxmin(i);
                vmax = kdt.boxmax(i);
                kdt.x(i) = vx;
                kdt.curboxmin(i) = vmin;
                kdt.curboxmax(i) = vmax;
                if( vx<vmin )
                {
                    kdt.curdist = amp::maximum<Precision>(kdt.curdist, vmin-vx);
                }
                else
                {
                    if( vx>vmax )
                    {
                        kdt.curdist = amp::maximum<Precision>(kdt.curdist, vx-vmax);
                    }
                }
            }
        }
        if( kdt.normtype==1 )
        {
            for(i=0; i<=kdt.nx-1; i++)
            {
                vx = x(i);
                vmin = kdt.boxmin(i);
                vmax = kdt.boxmax(i);
                kdt.x(i) = vx;
                kdt.curboxmin(i) = vmin;
                kdt.curboxmax(i) = vmax;
                if( vx<vmin )
                {
                    kdt.curdist = kdt.curdist+vmin-vx;
                }
                else
                {
                    if( vx>vmax )
                    {
                        kdt.curdist = kdt.curdist+vx-vmax;
                    }
                }
            }
        }
        if( kdt.normtype==2 )
        {
            for(i=0; i<=kdt.nx-1; i++)
            {
                vx = x(i);
                vmin = kdt.boxmin(i);
                vmax = kdt.boxmax(i);
                kdt.x(i) = vx;
                kdt.curboxmin(i) = vmin;
                kdt.curboxmax(i) = vmax;
                if( vx<vmin )
                {
                    kdt.curdist = kdt.curdist+amp::sqr<Precision>(vmin-vx);
                }
                else
                {
                    if( vx>vmax )
                    {
                        kdt.curdist = kdt.curdist+amp::sqr<Precision>(vx-vmax);
                    }
                }
            }
        }
    }


    /*************************************************************************
    Returns norm_k(x)^k (root-free = faster, but preserves ordering)

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> vrootfreenorm(const ap::template_1d_array< amp::ampf<Precision> >& x,
        int n,
        int normtype)
    {
        amp::ampf<Precision> result;
        int i;


        result = 0;
        if( normtype==0 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = amp::maximum<Precision>(result, amp::abs<Precision>(x(i)));
            }
            return result;
        }
        if( normtype==1 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = result+amp::abs<Precision>(x(i));
            }
            return result;
        }
        if( normtype==2 )
        {
            result = 0;
            for(i=0; i<=n-1; i++)
            {
                result = result+amp::sqr<Precision>(x(i));
            }
            return result;
        }
        return result;
    }


    /*************************************************************************
    Returns norm_k(x)^k (root-free = faster, but preserves ordering)

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> vrootfreecomponentnorm(amp::ampf<Precision> x,
        int normtype)
    {
        amp::ampf<Precision> result;


        result = 0;
        if( normtype==0 )
        {
            result = amp::abs<Precision>(x);
        }
        if( normtype==1 )
        {
            result = amp::abs<Precision>(x);
        }
        if( normtype==2 )
        {
            result = amp::sqr<Precision>(x);
        }
        return result;
    }


    /*************************************************************************
    Returns range distance: distance from X to [A,B]

      -- ALGLIB --
         Copyright 28.02.2010 by Bochkanov Sergey
    *************************************************************************/
    template<unsigned int Precision>
    amp::ampf<Precision> vrangedist(amp::ampf<Precision> x,
        amp::ampf<Precision> a,
        amp::ampf<Precision> b)
    {
        amp::ampf<Precision> result;


        if( x<a )
        {
            result = a-x;
        }
        else
        {
            if( x>b )
            {
                result = x-b;
            }
            else
            {
                result = 0;
            }
        }
        return result;
    }
} // namespace

#endif
