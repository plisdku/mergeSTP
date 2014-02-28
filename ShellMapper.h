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

Vector2<Point_3> lineProjections(const Line_3 & line)
{
    return Vector2<Point_3>(line.projection(Point_3(1,0,0)),
        line.projection(Point_3(0,1,0)));
}

struct LineLT
{
    bool operator()(const Line_3 & l1, const Line_3 & l2) const
    {
        return lineProjections(l1) < lineProjections(l2);
    }
};


class Contour
{
public:
    Contour()
    {
    }
    
    Contour(const std::vector<Point_3> & points) :
        mPoints(points),
        mLines(points.size())
    {
        for (int vv = 0; vv < points.size(); vv++)
        {
            int ww = (vv+1) % points.size();
            mLines[vv] = Line_3(points[vv], points[ww]);
        }
    }
    
    Contour(const std::vector<Point_3> & points,
        const std::vector<Line_3> & lines) :
        mPoints(points), mLines(lines)
    {
    }
    
    const std::vector<Point_3> & points() const { return mPoints; }
    const std::vector<Line_3> & lines() const { return mLines; }
    
private:
    std::vector<Point_3> mPoints;
    std::vector<Line_3> mLines;
};

class Face
{
public:
    Face();
    Face(const Plane_3 & p) :
        mPlane(p)
    {
    }
    
    Face(const Plane_3 & p, const std::vector<Contour> & contours) :
        mPlane(p),
        mContours(contours)
    {
    }
    
    void plane(const Plane_3 & p) { mPlane = p; }
    void addContour(const Contour & c) { mContours.push_back(c); }
    
    const Plane_3 & plane() const { return mPlane; }
    const std::vector<Contour> & contours() const { return mContours; }
    
private:
    std::vector<Contour> mContours;
    Plane_3 mPlane;
};

typedef vector<Face> Shell;

class ShellMapper
{
public:
    typedef NefPolyhedron::Halffacet_cycle_const_iterator
        HalffacetCycleConstIterator;
    typedef NefPolyhedron::Nef_polyhedron_S2 NefPolyhedronS2;
    typedef NefPolyhedronS2::SHalfedge_const_handle SHalfedgeConstHandle;
    
    ShellMapper(std::set<Point_3> & inSharedPoints,
        std::set<Line_3, LineLT> & inSharedLines,
        std::set<Plane_3, PlaneLT> & inSharedPlanes) :
        mUniquePoints(inSharedPoints),
        mUniqueLines(inSharedLines),
        mUniquePlanes(inSharedPlanes)
    {
    }
    
    void visit(NefPolyhedron::Halffacet_const_handle facetItr)
    {
        std::vector<Rect3d> boundingBoxes;
        NefPolyhedron::Plane_3 plane(facetItr->plane().opposite());
        // I flip the plane because shells point into volumes, not out of them.
        // This seems ridiculous.
        
        // *** Save the plane!  This will discard plane orientations.
        
        mUniquePlanes.insert(plane);
        
        vector<Contour> faceContours;
        
        HalffacetCycleConstIterator itr;
        for (itr = facetItr->facet_cycles_begin();
            itr != facetItr->facet_cycles_end();
            itr++)
        {
            if (itr.is_shalfedge())
            {
                SHalfedgeConstHandle firstEdgeHandle(itr);
                SHalfedgeConstHandle edgeHandle = firstEdgeHandle;
                
                Contour newContour;
                
                vector<Point_3> contourPoints;
                vector<Line_3> contourLines;
                
                Rect3d boundingBox(sMaxBox());
                
                do {
                    edgeHandle = edgeHandle->prev(); // to make it face outwards
                    Point_3 pointOnContour(edgeHandle->source()->center_vertex()->point());
                    //newContour.push_back(pointOnContour);
                    contourPoints.push_back(pointOnContour);
                    
                    Vector3d vert(toVector3d(pointOnContour));
                    boundingBox.p1 = vec_min(boundingBox.p1, vert);
                    boundingBox.p2 = vec_max(boundingBox.p2, vert);
                    
//                    if (0 == mUniquePoints.count(pointOnContour))
//                        std::cout << "Inserting point " << pointOnContour << ".\n";
                    
                    mUniquePoints.insert(pointOnContour); // *** Save the point!
                } while (edgeHandle != firstEdgeHandle);
                
                // *** Save the contour!
                boundingBoxes.push_back(boundingBox);
                
                // *** Save the unique lines!
                for (int nn = 0; nn < contourPoints.size(); nn++)
                {
                    int mm = (nn+1) % contourPoints.size();
                    
                    Line_3 ln(contourPoints[nn], contourPoints[mm]);
                    contourLines.push_back(ln);
                    mUniqueLines.insert(ln);
                }
                
                faceContours.push_back(Contour(contourPoints, contourLines));
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
        swap(faceContours[0], faceContours[biggestBoxIndex]);
        // Now the contours are ordered with the outside box first, and all
        // the "hole" contours are after that.
        
        assert(mShells.size() > 0);
        Shell & currentShell = mShells[mShells.size() - 1];
        currentShell.push_back(Face(plane, faceContours));
    }
    
    void newShell()
    {
        mShells.push_back(std::vector<Face>());
    }
    
    void visit(NefPolyhedron::Vertex_const_handle h) {}
    void visit(NefPolyhedron::Halfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfedge_const_handle h) {}
    void visit(NefPolyhedron::SHalfloop_const_handle h) {}
    void visit(NefPolyhedron::SFace_const_handle) {}
    
    const std::vector<Shell> shells() const { return mShells; }
    
    static Rect3d sMaxBox()
    {
        return Rect3d(std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max(),
            -std::numeric_limits<double>::max());
    }
private:
    std::set<Point_3> & mUniquePoints;
    std::set<Line_3, LineLT> & mUniqueLines;
    std::set<Plane_3, PlaneLT> & mUniquePlanes;
        
    std::vector<Shell> mShells;
};


#endif



