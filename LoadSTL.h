/*
 *  LoadSTL.h
 *  mergeSTP
 *
 *  Created by Paul C Hansen on 2/28/14.
 *  Copyright 2014 Stanford University. All rights reserved.
 *
 */

#ifndef LOADSTL_H
#define LOADSTL_H

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Polyhedron_3.h>
#include <sstream>
#include <cstdio>
#include <string>

namespace LoadSTL
{

// A modifier creating a triangle with the incremental builder.
template <class HDS>
class BuildSTL : public CGAL::Modifier_base<HDS>
{
public:
    typedef typename HDS::Vertex Vertex;
    typedef typename Vertex::Point Point;
    
    class ComparePoint
    {
    public:
        bool operator()(const Point & lhs, const Point & rhs) const
        {
            for (int nn = 0; nn < 3; nn++)
            if (lhs[nn] < rhs[nn])
                return true;
            else if (rhs[nn] < lhs[nn])
                return false;
            
            return false;
        }
    };
    
    BuildSTL(std::istream & str) :
        mStream(str)
    {
    }
    
    void operator()( HDS& hds)
    {
        // Postcondition: hds is a valid polyhedral surface.
        CGAL::Polyhedron_incremental_builder_3<HDS> builder( hds, true);
        
        std::string line, token;
        
        std::map<Point, int, ComparePoint> vertexIndices;
        std::vector<Point> verticesInOrder;
        
        std::vector<int> triVerts[3]; // vertex indices for triangles!
        
        while (mStream)
        {
            getline(mStream, line);
            istringstream stream(line);
            
            cout << "Line is " << line << "\n";
            
            stream >> token;
            
            int triPoints[3];
            int curTriangleVertex;
            
            if (token == "solid")
            {
                // yay.
            }
            else if (token == "facet")
            {   
//                builder.begin_facet();
            }
            else if (token == "outer")
            {
                // outer loop
                curTriangleVertex = 0;
            }
            else if (token == "vertex")
            {
                float vx, vy, vz;
                
                stream >> vx >> vy >> vz;
                
                Point p3(vx,vy,vz);
                
                int vIndex; // the global index of this vertex!
                
                if (0 == vertexIndices.count(p3))
                {
                    int numVertices = verticesInOrder.size();
                    verticesInOrder.push_back(p3);
                    vertexIndices[p3] = numVertices;
                    vIndex = numVertices;
                    
                    cout << "Added new vertex.\n";
                }
                else
                {
                    vIndex = vertexIndices[p3];
                    cout << "Noted duplicate vertex.\n";
                }
                
                if (curTriangleVertex >= 0 && curTriangleVertex < 3)
                    triVerts[curTriangleVertex++].push_back(vIndex);
                else
                    throw(std::runtime_error("Too many vertices in loop."));
            }
            else if (token == "endloop")
            {
                if (curTriangleVertex == 3)
                {
//                    builder.add_vertex_to_facet(triPoints[0]);
//                    builder.add_vertex_to_facet(triPoints[1]);
//                    builder.add_vertex_to_facet(triPoints[2]);
                }
                else
                    throw(std::runtime_error("Too few vertices to close loop."));
            }
            else if (token == "endfacet")
            {
//                builder.end_facet();
            }
        }
        assert(triVerts[0].size() == triVerts[1].size());
        assert(triVerts[1].size() == triVerts[2].size());
        
        // Use the builder!!  Here I dump the points and faces into the
        // polyhedral data structure.
        builder.begin_surface(verticesInOrder.size(), triVerts[0].size());
        
        for (int vv = 0; vv < verticesInOrder.size(); vv++)
        {
            builder.add_vertex(verticesInOrder[vv]);
            cout << "Vertex " << vv << " is " << verticesInOrder[vv]
                << " map " << vertexIndices[verticesInOrder[vv]] << "\n";
        }
        
        for (int ff = 0; ff < triVerts[0].size(); ff++)
        {
            cout << "Face " << ff << " is vertices ";
            builder.begin_facet();
            for (int vv = 0; vv < 3; vv++)
            {
                builder.add_vertex_to_facet(triVerts[vv][ff]);
                cout << triVerts[vv][ff] << " ";
            }
            cout << "\n";
            builder.end_facet();
        }
        
        builder.end_surface();
    }
private:
    std::istream & mStream;
};

} // namespace LoadSTL


#endif
