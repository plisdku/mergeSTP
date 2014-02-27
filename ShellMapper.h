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

// This is essentially the unique key of a plane AND ITS OPPOSITE.
Vector3<Point_3> planeProjections(const NefPolyhedron::Plane_3 & plane)
{
    return Vector3<Point_3>(plane.projection(Point_3(1,0,0)),
        plane.projection(Point_3(0,1,0)),
        plane.projection(Point_3(0,0,1)));
}


struct PlaneLT
{
    bool operator()(const Plane_3 & p1, const Plane_3 & p2) const
    {
        return planeProjections(p1) < planeProjections(p2);
    }
};

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

struct FaceAndPlane
{
    Face face;
    Plane_3 plane;
};

bool operator<(const FaceAndPlane & lhs, const FaceAndPlane & rhs)
{
    return lhs.face < rhs.face;
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
        
        // *** Save the plane!  Actually this will essentially discard orientation.
        mCGALPlanes.insert(plane);
        
        FaceAndPlane faceRecord;
        faceRecord.plane = plane;
        
        HalffacetCycleConstIterator itr;
        for (itr = facetItr->facet_cycles_begin();
            itr != facetItr->facet_cycles_end();
            itr++)
        {
            if (itr.is_shalfedge())
            {
                SHalfedgeConstHandle handle(itr);
                SHalfedgeConstHandle edgeHandle = handle;
                
                Contour contour;
                
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
                    Segment seg = Segment(contour[nn-1], contour[nn]);
                    mCGALSegments.insert(seg);
                }
                mCGALSegments.insert(Segment(contour[contour.size()-1], contour[0]));
                
//                for (int nn = 0; nn < contour.size()-1; nn++)
//                {
//                    cout << "Contour " << contour[nn] << " to " << contour[nn+1] << "\n";
//                }
                
                // *** Save the contours!
                boundingBoxes.push_back(boundingBox);
                faceRecord.face.push_back(contour);
                
                mContours.insert(contour);
                
                //faceKey.push_back(contourKey);
                //mContourMap[contourKey] = contour;
                //mContourMap[contour] = contour;
                
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
        swap(faceRecord.face[0], faceRecord.face[biggestBoxIndex]);
        // Now the contours are ordered with the outside box first, and all
        // the "hole" contours are after that.
        
        // *** Save the face!
//        faceKey = contours; // now all the contours are indeed oriented.
//        std::sort(faceKey.begin(), faceKey.end());
//        mFaceKeys.insert(faceKey);
//        //mFaceMap[faceKey] = contours;
//        mFaceMap[contours] = contours;
        
        mFaces.insert(faceRecord); // TODO: fix this.
        // The face key is, still, a disoriented face.  Hmm.
        
        // *** Save the face for this shell, too
        
        assert(mShellFaces.size());
        mShellFaces[mShellFaces.size()-1].push_back(faceRecord.face);
    }
    
    void newShell()
    {
        mShellFaces.push_back(std::vector<Face>());
    }
    
    void visit(NefPolyhedron::Vertex_const_handle h) {}
    void visit(NefPolyhedron::Halfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfloop_const_handle h) {}
    void visit(NefPolyhedron::SFace_const_handle) {}
    
    const std::set<Point_3> & CGALPoints() const { return mCGALPoints; }
    const std::set<Plane_3, PlaneLT> & planes() const { return mCGALPlanes; }
    
    const std::set<Segment> & CGALSegments() const { return mCGALSegments; }
    const std::set<Contour> & contours() const { return mContours; }
    const std::set<FaceAndPlane> & faces() const { return mFaces; }
    
    const std::vector<std::vector<Face> > shellFaces() const { return mShellFaces; }
    
private:
    std::set<Point_3> mCGALPoints;
    std::set<Segment> mCGALSegments;
    std::set<Contour> mContours;
    std::set<FaceAndPlane> mFaces;
    std::set<Plane_3, PlaneLT> mCGALPlanes;
    
    std::vector<std::vector<Face> > mShellFaces;
};


#endif



