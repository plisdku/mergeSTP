/*
 *  ShellMapper.h
 *  NefLab
 *
 *  Created by Paul C Hansen on 2/24/14.
 *  Copyright 2014 Stanford University. All rights reserved.
 *
 */

#ifndef SHELLMAPPER_H
#define SHELLMAPPER_H

#include "utility/VectorMatrix2.h"
#include "CGALUtilities.h"

bool operator<(Point_3 p1, Point_3 p2)
{
    for (int nn = 0; nn < 3; nn++)
    {
        if (p1[nn] < p2[nn])
            return true;
        else if (p1[nn] > p2[nn])
            return false;
    }
    return false;
}

bool operator<(const NefPolyhedron::Vertex_const_handle & v1,
    const NefPolyhedron::Vertex_const_handle & v2)
{
    return v1->point() < v2->point();
}

bool operator<(const NefPolyhedron::Plane_3 & p1,
    const NefPolyhedron::Plane_3 & p2)
{
    if (p1.a() < p2.a())
        return true;
    else if (p1.a() > p2.a())
        return false;
    else if (p1.b() < p2.b())
        return true;
    else if (p1.b() > p2.b())
        return false;
    else if (p1.c() < p2.c())
        return true;
    else if (p1.c() > p2.c())
        return false;
    else if (p1.d() < p2.d())
        return true;
    return false;
}

// A Segment is a sorted pair of points with p1 < p2.  (So, it's not
// oriented.  The sorting is to wipe out orientation.)
typedef Vector2<Point_3> Segment;
Segment toSegment(Point_3 p1, Point_3 p2)
{
    if (p1 < p2)
        return Segment(p1,p2);
    else
        return Segment(p2,p1);
}

// A Contour should of course be sorted and oriented and all that,
// but if I sort the points alphabetically then I can compare contours and
// use them as keys in data structures.  Hee.  :-D

typedef vector<Point_3> Contour;
typedef vector<Point_3> ContourKey;

//vector<Point_3> toSortedContour(
bool operator<(const vector<Point_3> & contour1,
    const vector<Point_3> & contour2)
{
    if (contour1.size() < contour2.size())
        return true;
    else if (contour1.size() > contour2.size())
        return false;
    
    for (int nn = 0; nn < contour1.size(); nn++)
    {
        if (contour1[nn] < contour2[nn])
            return true;
        else if (contour2[nn] < contour1[nn])
            return false;
    }
    return false;
}

// Faces are made of one or more contours.  To make them sortable I will
// define the same ordering as I did for ContourKey, above.  So that's a
// FaceKey.  Woo.

typedef vector<Contour> Face;
typedef vector<ContourKey> FaceKey;

bool operator<(const FaceKey & f1, const FaceKey & f2)
{
    if (f1.size() < f2.size())
        return true;
    else if (f1.size() > f2.size())
        return false;
    
    for (int nn = 0; nn < f1.size(); nn++)
    {
        if (f1[nn] < f2[nn])
            return true;
        else if (f2[nn] < f1[nn])
            return false;
    }
    return false;
}




class ShellMapper
{
public:
    typedef NefPolyhedron::Halffacet_cycle_const_iterator
        HalffacetCycleConstIterator;
    typedef NefPolyhedron::Nef_polyhedron_S2 NefPolyhedronS2;
    typedef NefPolyhedronS2::SHalfedge_const_handle SHalfedgeConstHandle;
    
    ShellMapper()
    {
        
    }
    
    void visit(NefPolyhedron::Halffacet_const_handle facetItr)
    {
        std::vector<Rect3d> boundingBoxes;
        NefPolyhedron::Plane_3 plane(facetItr->plane().opposite());
        // I flip the plane because shells point into volumes, not out of them.
        // This seems ridiculous.
        
        // This vector will become the FaceKey once I've sorted it.
        std::vector<ContourKey> faceKey;
        std::vector<Contour> contours;
        
        HalffacetCycleConstIterator itr;
        for (itr = facetItr->facet_cycles_begin();
            itr != facetItr->facet_cycles_end();
            itr++)
        {
            
            if (itr.is_shalfedge())
            {
                SHalfedgeConstHandle handle(itr);
                SHalfedgeConstHandle edgeHandle = handle;
                
                std::vector<Point_3> contour;
                
                Rect3d boundingBox(std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max(),
                    std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max(),
                    -std::numeric_limits<double>::max());
                
                do {
                    edgeHandle = edgeHandle->prev(); // to make it face outwards
                    NefPolyhedron::Point_3 p3(
                        edgeHandle->source()->center_vertex()->point());
                    contour.push_back(p3);
                    
                    Vector3d vert(toVector3d(p3));
                    boundingBox.p1 = vec_min(boundingBox.p1, vert);
                    boundingBox.p2 = vec_max(boundingBox.p2, vert);
                    
                    mCGALPoints.insert(p3); // *** Save the point!
                } while (edgeHandle != handle);
                
                // *** Save the segments!
                for (int nn = 1; nn < contour.size(); nn++)
                {
                    Segment seg = toSegment(contour[nn-1], contour[nn]);
                    mCGALSegments.insert(seg);
                }
                mCGALSegments.insert(toSegment(contour[0], contour[contour.size()-1]));
                
                for (int nn = 0; nn < contour.size()-1; nn++)
                {
                    cout << "Contour " << contour[nn] << " to " << contour[nn+1] << "\n";
                }
                
                // *** Save the contours!
                boundingBoxes.push_back(boundingBox);
                contours.push_back(contour);
                
                ContourKey contourKey = contour;
                std::sort(contourKey.begin(), contourKey.end());
                mContourKeys.insert(contourKey);
                faceKey.push_back(contourKey);
                mContourMap[contourKey] = contour;
            }
            else if (itr.is_shalfloop())
            {
                //cerr << "Found isolated vertex; ignoring it!\n";
            }
            else
            {
                //cerr << "\tFound cycle, it's a mys it's a mys it's a mystery\n";
            }
        }
        
        // Find the largest contour.  Its bounding box is biggest...
        int biggestBoxIndex = 0;
        for (int bb = 1; bb < boundingBoxes.size(); bb++)
        {
            if (sumSquares(boundingBoxes[bb].size()) >
                sumSquares(boundingBoxes[biggestBoxIndex].size()))
            {
                biggestBoxIndex = bb;
            }
        }
        swap(contours[0], contours[biggestBoxIndex]);
        // Now the contours are ordered with the outside box first, and all
        // the "hole" contours are after that.
        
        // *** Save the face!
        std::sort(faceKey.begin(), faceKey.end());
        mFaceKeys.insert(faceKey);
        mFaceMap[faceKey] = contours;
        
        // *** Save the face for this shell, too
        
        assert(mShellFaceKeys.size());
        mShellFaceKeys[mShellFaceKeys.size()-1].push_back(faceKey);
    }
    
    void newShell()
    {
        mShellFaceKeys.push_back(std::vector<FaceKey>());
    }
    
    void visit(NefPolyhedron::Vertex_const_handle h) {}
    void visit(NefPolyhedron::Halfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfloop_const_handle h) {}
    void visit(NefPolyhedron::SFace_const_handle) {}
    
    const std::set<Point_3> & CGALPoints() const { return mCGALPoints; }
    const std::set<Segment> & CGALSegments() const { return mCGALSegments; }
    const std::set<ContourKey> & contourKeys() const { return mContourKeys; }
    const std::set<FaceKey> & faceKeys() const { return mFaceKeys; }
    const std::map<ContourKey, Contour> & contourMap() const { return mContourMap; }
    const std::map<FaceKey, Face> & faceMap() const { return mFaceMap; }
    const std::vector<std::vector<FaceKey> > shellFaceKeys() const { return mShellFaceKeys; }
    
private:
    std::set<Point_3> mCGALPoints;
    std::set<Segment> mCGALSegments;
    std::set<ContourKey> mContourKeys;
    std::set<FaceKey> mFaceKeys;
    
    std::map<ContourKey, Contour> mContourMap;
    std::map<FaceKey, Face> mFaceMap;
    
    std::vector<std::vector<FaceKey> > mShellFaceKeys;
};


#endif



