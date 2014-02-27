#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> NefPolyhedron;
typedef NefPolyhedron::Point_3 Point_3;
typedef NefPolyhedron::Plane_3 Plane_3;
typedef NefPolyhedron::Traits Traits;
typedef CGAL::Polyhedron_3<Traits> Polyhedron;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Direction_3 Direction_3;

using namespace std;

//#include "ShellVisitor.h"
#include "InteriorVolumes-inl.h"

#include <gp_Pln.hxx>
#include <TopoDS.hxx> // for casting one shape to another
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_CompSolid.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeShell.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepTools.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <BRepCheck_Face.hxx>
#include <BRepCheck_Wire.hxx>
#include <BRepCheck_Shell.hxx>
#include <ShapeUpgrade_ShellSewing.hxx>

#include <StlAPI.hxx>
#include <StlAPI_Reader.hxx>

#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>

#include "ShellMapper.h"

void printHelp();


void translatePoints(const set<Point_3> & inPoints,
    map<Point_3, TopoDS_Vertex> & outMap);
void translateLines(const set<Line_3, LineLT> & inLines,
    map<Line_3, gp_Lin, LineLT> & outMap);
void translatePlanes(const set<Plane_3, PlaneLT> & inPlanes,
    map<Plane_3, gp_Pln, PlaneLT> & outMap);
TopoDS_Solid translateSolid(const vector<Shell> & shells,
    const map<Point_3, TopoDS_Vertex> & pointMap,
    const map<Line_3, gp_Lin, LineLT> & lineMap,
    const map<Plane_3, gp_Pln, PlaneLT> & planeMap);


std::vector<Shell> shellNefPolyhedron(NefPolyhedron p,
    set<Point_3> & accumPoints,
    set<Line_3, LineLT> & accumLines,
    set<Plane_3, PlaneLT> accumPlanes);
//    , 
//    map<Point_3, gp_Pnt> & convertPoints,
//    map<Line_3, gp_Lin, LineLT> & convertLines,
//    map<Plane_3, gp_Pln, PlaneLT> & convertPlanes);

TopoDS_Compound mergeSolids(const vector<vector<Shell> > & solids,
    const set<Point_3> & uniquePoints,
    const set<Line_3, LineLT> & uniqueLines,
    const set<Plane_3, PlaneLT> & uniquePlanes);

NefPolyhedron loadMultiOFF(istream & instr);
NefPolyhedron loadSTL(ifstream & instr);
NefPolyhedron loadOFF(ifstream & instr);

//void nefToBRep(const NefPolyhedron & poly);

string sFileExtension(const string & fname)
{
    int dotPosition = fname.find_last_of(".");
    
    if (dotPosition == std::string::npos)
        throw(std::runtime_error("File has no extension"));
    
    return fname.substr(dotPosition+1);
}

int main(int argc, char const** argv)
{
    if (argc == 1)
    {
        printHelp();
        return 0;
    }
    
    set<Point_3> uniquePoints;
    set<Line_3, LineLT> uniqueLines;
    set<Plane_3, PlaneLT> uniquePlanes;
    vector<vector<Shell> > solids;
    
    int numInputFiles = argc - 1;
    for (int inputFileNum = 0; inputFileNum < numInputFiles; inputFileNum++)
    {
        string fileName(argv[inputFileNum+1]);
        string fileExt;
        
        NefPolyhedron nef;
        vector<Shell> shells;
        
        bool loaded = true;
        
        try
        {
            ifstream inStream(fileName.c_str());
            
            if (inStream.bad())
                throw(std::runtime_error("Stream is bad."));
            else if (inStream.fail())
                throw(std::runtime_error("Stream failed."));
            else if (inStream.eof())
                throw(std::runtime_error("Stream empty."));
            else if (false == inStream.good())
                throw(std::runtime_error("Stream is not good."));
            
            fileExt = sFileExtension(fileName);
            
            if (fileExt == "txt")
            {
                cout << "Loading " << fileName << " as multi-OFF file.\n";
                
                nef = loadMultiOFF(inStream);
            }
            else if (fileExt == "stl")
            {
                cout << "Loading " << fileName << " as STL file.\n";
                
                nef = loadSTL(inStream);
            }
            else if (fileExt == "off")
            {
                cout << "Loading " << fileName << " as OFF file.\n";
                
                nef = loadOFF(inStream);
            }
            else
            {
                cout << "Unsure how to handle " << fileExt << " file; skipping "
                    << fileName << ".\n";
            }
        }
        catch (const std::exception & exc)
        {
            cerr << "Can't load " << fileName << ": ";
            cerr << exc.what() << "\n";
            loaded = false;
        }
        
        if (loaded)
        {
            shells = shellNefPolyhedron(nef,
                uniquePoints, uniqueLines, uniquePlanes);
            solids.push_back(shells);
        }
    }
    
    TopoDS_Compound allSolids = mergeSolids(solids,
        uniquePoints, uniqueLines, uniquePlanes);
    
    STEPControl_Writer stepWriter;
    
    // 0: writes STEP files without assemblies
    // 1: writes all shapes in the form of assemblies
    // 2: writes shapes having a structure of (possibly nested) TopoDS_Compounds in the form of
    //    STEP assemblies; single shapes are written without assembly structures.
    assert(Interface_Static::SetIVal("write.step.assembly", 2));
    
//    assert(Interface_Static::SetCVal("write.step.schema","1"));
    
    stepWriter.Transfer(allSolids, STEPControl_AsIs);
    
    return 0;
}

void printHelp()
{
    cout << "Usage:\n";
    cout << "\tNefCascade\n";
}


NefPolyhedron loadMultiOFF(istream & instr)
{
    int numPositiveShells, numNegativeShells;
    int numIgnoredFacets;
    
    instr >> numPositiveShells;
    instr >> numNegativeShells;
    
    NefPolyhedron nef;
    
    for (int nn = 0; nn < numPositiveShells; nn++)
    {
        cerr << "Adding " << nn << " of " << numPositiveShells << ".\n";
        NefPolyhedron addend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(instr, addend);
        if (numIgnoredFacets > 0)
            cerr << "Ignored " << numIgnoredFacets << " facets!\n";
        nef = nef + addend;
    }
    
    for (int nn = 0; nn < numNegativeShells; nn++)
    {
        cerr << "Subtracting " << nn << " of " << numNegativeShells << ".\n";
        NefPolyhedron subtrahend;
        
        numIgnoredFacets = CGAL::OFF_to_nef_3(instr, subtrahend);
        if (numIgnoredFacets > 0)
            cerr << "Ignored " << numIgnoredFacets << " facets!\n";
        nef = nef - subtrahend;
    }
    
    return nef;
}

NefPolyhedron loadSTL(ifstream & instr)
{
    throw(std::logic_error("I can't load STL yet."));
}

NefPolyhedron loadOFF(ifstream & instr)
{
    NefPolyhedron poly;
    int numIgnoredFacets = CGAL::OFF_to_nef_3(instr, poly);
    
    if (numIgnoredFacets > 0)
        cerr << "Ignored " << numIgnoredFacets << " facets!\n";
    
    return poly;
}


// This function converts one Nef polyhedron to a vector of shells.
// Ideally solids with coplanar, touching faces should be identifiable as
// such to programs that use the STEP file we're making here.  I try to make
// this as likely as possible by not creating any redundant TopoDS_Vertex or
// gp_Pln objects (vertices and surfaces can be shared).
//
// Therefore:
// accumPoints builds up the list of UNIQUE vertices for the STEP file
// accumPlanes builds up the list of UNIQUE surfaces for the STEP file
// convertPoints stores the unique TopoDS_Vertex located at each Point_3
// convertPlanes stores the unique gp_Pln located at each Plane_3
//
// Planes are orientable (there's an up-plane and a down-plane) but I only save
// one of each orientation, i.e. to me a plane has no up or down orientation.
// I don't believe that OpenCASCADE needs to know the surface orientation for
// anything that I'm doing with it.
std::vector<Shell> shellNefPolyhedron(NefPolyhedron poly,
    set<Point_3> & accumPoints,
    set<Line_3, LineLT> & accumLines,
    set<Plane_3, PlaneLT> accumPlanes)
{
    typedef NefPolyhedron::Halffacet_cycle_const_iterator
        HalffacetCycleConstIterator;
    typedef NefPolyhedron::Volume_const_handle VolumeConstHandle;
    typedef NefPolyhedron::Shell_entry_const_iterator ShellEntryConstIterator;
    typedef NefPolyhedron::SFace_const_handle SFaceConstHandle;
    
    ShellMapper shellMapper(accumPoints, accumLines, accumPlanes);
    
    vector<VolumeConstHandle> insideVolumes =
        InteriorVolumes::interiorVolumes(poly);
    
    if (insideVolumes.size() > 1)
    {
        cerr << "Error: I can't handle more than one volume yet.\n";
        assert(false);
    }
    
    for (int volNum = 0; volNum < insideVolumes.size(); volNum++)
    {
        cerr << "Volume " << volNum << "\n";
        
        int shell = 0;
        
        ShellEntryConstIterator itr;
        for (itr = insideVolumes[volNum]->shells_begin();
            itr != insideVolumes[volNum]->shells_end(); itr++)
        {
            shellMapper.newShell();
            
            cout << "\tShell " << shell++ << "\n";
            poly.visit_shell_objects(SFaceConstHandle(itr),
                shellMapper);
        }
    }
    
    // Now the NefPolyhedron has been decomposed into its shells by the
    // shell mapper.
    
    return shellMapper.shells();
}

// This function just tosses a bunch of TopoDS_Shapes into one TopoDS_Compound.
TopoDS_Compound mergeSolids(const vector<vector<Shell> > & solids,
    const set<Point_3> & uniquePoints,
    const set<Line_3, LineLT> & uniqueLines,
    const set<Plane_3, PlaneLT> & uniquePlanes)
{
    map<Point_3, TopoDS_Vertex> convertPoints;
    map<Line_3, gp_Lin, LineLT> convertLines;
    map<Plane_3, gp_Pln, PlaneLT> convertPlanes;
    
    translatePoints(uniquePoints, convertPoints);
    translateLines(uniqueLines, convertLines);
    translatePlanes(uniquePlanes, convertPlanes);
    
    // Time to handle the topology for each solid!!
    
    TopoDS_Compound compound;
    BRep_Builder buildCompound;
    buildCompound.MakeCompound(compound);
    
    for (int iSolid = 0; iSolid < solids.size(); iSolid++)
    {
        TopoDS_Solid solid = translateSolid(solids[iSolid], 
            convertPoints, convertLines, convertPlanes);
        buildCompound.Add(compound, solid);
    }
    
    return compound;
}

void translatePoints(const set<Point_3> & inPoints,
    map<Point_3, TopoDS_Vertex> & outMap)
{
    for (set<Point_3>::const_iterator itr = inPoints.begin();
        itr != inPoints.end(); itr++)
    {
        outMap[*itr] = BRepBuilderAPI_MakeVertex(gp_Pnt(
            CGAL::to_double(itr->x()),
            CGAL::to_double(itr->y()),
            CGAL::to_double(itr->z())));
    }
}

void translateLines(const set<Line_3, LineLT> & inLines,
    map<Line_3, gp_Lin, LineLT> & outMap)
{
    for (set<Line_3, LineLT>::const_iterator itr = inLines.begin();
        itr != inLines.end(); itr++)
    {
        Point_3 inPt = itr->point(0); // the 0 is sort of a seed value.
        gp_Pnt outPt(CGAL::to_double(inPt[0]),
            CGAL::to_double(inPt[1]),
            CGAL::to_double(inPt[2]));
            
        Direction_3 inDir = itr->direction();
        gp_Dir outDir(CGAL::to_double(inDir.dx()),
            CGAL::to_double(inDir.dy()),
            CGAL::to_double(inDir.dz()));
        
        outMap[*itr] = gp_Lin(outPt, outDir);
    }
}

void translatePlanes(const set<Plane_3, PlaneLT> & inPlanes,
    map<Plane_3, gp_Pln, PlaneLT> & outMap)
{
    for (set<Plane_3, PlaneLT>::const_iterator itr = inPlanes.begin();
        itr != inPlanes.end(); itr++)
    {
        outMap[*itr] = gp_Pln(CGAL::to_double(itr->a()),
            CGAL::to_double(itr->b()),
            CGAL::to_double(itr->c()),
            CGAL::to_double(itr->d()));
    }
}

TopoDS_Solid translateSolid(const vector<Shell> & shells,
    const map<Point_3, TopoDS_Vertex> & pointMap,
    const map<Line_3, gp_Lin, LineLT> & lineMap,
    const map<Plane_3, gp_Pln, PlaneLT> & planeMap)
{
    vector<TopoDS_Shell> newShells(shells.size());
    
    for (int iShell = 0; iShell < shells.size(); iShell++)
    {
        const Shell & shell = shells.at(iShell);
        
        TopoDS_Shell newShell;
        BRep_Builder buildShell;
        buildShell.MakeShell(newShell);
        
        for (int iFace = 0; iFace < shell.size(); iFace++)
        {
            const Face & face = shell.at(iFace);
            
            vector<TopoDS_Wire> newWires(face.contours().size());
            
            for (int iContour = 0; iContour < face.contours().size(); iContour++)
            {
                const Contour & contour = face.contours()[iContour];
                
                BRep_Builder buildWire;
                buildWire.MakeWire(newWires[iContour]);
                
                for (int iVertex = 0; iVertex < contour.points().size(); iVertex++)
                {
                    int iNext = (iVertex+1) % contour.points().size();
                    
                    assert(pointMap.count(contour.points()[iVertex]));
                    assert(pointMap.count(contour.points()[iNext]));
                    buildWire.Add(newWires[iContour], BRepBuilderAPI_MakeEdge(
                        lineMap.find(contour.lines()[iVertex])->second,
                        pointMap.find(contour.points()[iVertex])->second,
                        pointMap.find(contour.points()[iNext])->second));
                }
                
                BRepCheck_Wire checkWire(newWires[iContour]);
                assert(checkWire.Closed() == BRepCheck_NoError);
                assert(checkWire.Orientation(TopoDS_Face()) == BRepCheck_NoError);
            }
            
            TopoDS_Face newFace;
            BRepBuilderAPI_MakeFace makeFace(planeMap.find(face.plane())->second,
                newWires[0]);
            
            for (int nn = 1; nn < newWires.size(); nn++)
                makeFace.Add(newWires[nn]);
            
            BRepCheck_Face checkFace(newFace);
            assert(checkFace.IntersectWires() == BRepCheck_NoError);
            assert(checkFace.ClassifyWires() == BRepCheck_NoError);
            assert(checkFace.OrientationOfWires() == BRepCheck_NoError);
            
            buildShell.Add(newShell, newFace);
        }
        
        ShapeUpgrade_ShellSewing sewing;
        TopoDS_Shape sewedShell = sewing.ApplySewing(newShell);
        assert(sewedShell.ShapeType() == TopAbs_SHELL);
        TopoDS_Shell & finalShell = newShells[iShell];
        finalShell = TopoDS::Shell(sewedShell);
        
        BRepCheck_Shell checkShell(finalShell);
        assert(checkShell.Closed() == BRepCheck_NoError);
        assert(checkShell.Orientation() == BRepCheck_NoError);
    }
    
    TopoDS_Solid newSolid;
    BRep_Builder buildSolid;
    buildSolid.MakeSolid(newSolid);
    
    for (int iShell = 0; iShell < newShells.size(); iShell++)
    {
        buildSolid.Add(newSolid, newShells.at(iShell));
    }
    
    return newSolid;
}

/*

void nefToBRep(const NefPolyhedron & poly)
{
    map<Point_3, TopoDS_Vertex> vertices;
    map<Segment, TopoDS_Edge> edges;
    map<Plane_3, gp_Pln, PlaneLT> planes;
    map<Contour, TopoDS_Wire> wires;
    map<Face, TopoDS_Face> faces;
    vector<TopoDS_Shell> shells;
    
    const set<Point_3> & points = shellVisitor.CGALPoints();
    const set<Segment> & segs = shellVisitor.CGALSegments();
    const set<Contour> & contours = shellVisitor.contours();
    const set<FaceAndPlane> & inFaces = shellVisitor.faces();
    const set<Plane_3, PlaneLT> & inPlanes = shellVisitor.planes();
    
    const vector<vector<Face> > & shellFaces = shellVisitor.shellFaces();
    
    cout << "\n---- Point list:\n";
    for (set<Point_3>::const_iterator itr = points.begin();
        itr != points.end(); itr++)
    {
        Point_3 p = *itr;
        //cout << p << "\n";
        
        vertices[p] = BRepBuilderAPI_MakeVertex(gp_Pnt(
            CGAL::to_double(p[0]),
            CGAL::to_double(p[1]),
            CGAL::to_double(p[2])));
    }
    
    cout << "\n---- Segment list:\n";
    for (set<Segment>::const_iterator itr = segs.begin();
        itr != segs.end(); itr++)
    {
        Segment s = *itr;
        cout << s[0] << " to " << s[1] << "\n";
        
        edges[s] = BRepBuilderAPI_MakeEdge(
            vertices[s[0]], vertices[s[1]]);
    }
    
    cout << "\n---- Contour list:\n";
    for (set<Contour>::const_iterator itr = contours.begin();
        itr != contours.end(); itr++)
    {
        const Contour & contour = *itr;
        cout << "Contour of length " << contour.size() << "\n";
        
        TopoDS_Wire newWire;
        BRep_Builder builder;
        builder.MakeWire(newWire);
        
        for (int nn = 0; nn < contour.size(); nn++)
        {
            int mm = (nn+1) % contour.size();
            // used to use toSegment, with the disoriented ShellMapper.
            Segment mySegment = Segment(contour[nn], contour[mm]);
            
            assert(edges.count(mySegment));
            builder.Add(newWire, edges[mySegment]);
        }
        
        //newWire.Closed(true);
        wires[contour] = newWire;
        
        BRepCheck_Wire checkWire(newWire);
        BRepCheck_Status status = checkWire.Closed();
        assert(status == BRepCheck_NoError);
        status = checkWire.Orientation(TopoDS_Face()); // using a null face since no face exists
        assert(status == BRepCheck_NoError);
        
        //for (int nn = 0; nn < contour.size(); nn++)
        //    cout << contour[nn] << " ";
        //cout << "\n";
    }
    
    cout << "\n---- Surface list:\n";
    for (set<Plane_3>::const_iterator itr = inPlanes.begin();
        itr != inPlanes.end(); itr++)
    {
        const Plane_3 & plane = *itr;
        
        planes[plane] = gp_Pln(CGAL::to_double(plane.a()),
            CGAL::to_double(plane.b()),
            CGAL::to_double(plane.c()),
            CGAL::to_double(plane.d()));
    }
    
    cout << "\n---- Face list:\n";
    for (set<FaceAndPlane>::const_iterator itr = inFaces.begin();
        itr != inFaces.end(); itr++)
    {
        const Face & face = itr->face;
        cout << "Face of length " << face.size() << "\n";
        
        assert(planes.count(itr->plane));
        
        BRepBuilderAPI_MakeFace makeFace(planes[itr->plane], wires[face[0]]);
        for (int nn = 1; nn < face.size(); nn++)
        {
            assert(wires.count(face[nn]));
            makeFace.Add(wires[face[nn]]);
        }
        
        TopoDS_Face newFace = makeFace;
        faces[face] = newFace;
        
        BRepCheck_Face checkFace(newFace);
        BRepCheck_Status status = checkFace.IntersectWires();
        assert(status == BRepCheck_NoError);
        status = checkFace.ClassifyWires();
        assert(status == BRepCheck_NoError);
        status = checkFace.OrientationOfWires();
        assert(status == BRepCheck_NoError);
    }
    
    cout << "\n---- Shell list:\n";
    for (int nn = 0; nn < shellFaces.size(); nn++)
    {
        
        TopoDS_Shell newShell;
        BRep_Builder builder;
        builder.MakeShell(newShell);
        
        for (int mm = 0; mm < shellFaces[nn].size(); mm++)
        {
            const Face & face = shellFaces[nn][mm];
            assert(faces.count(face));
            builder.Add(newShell, faces[face]);
        }
        
        ShapeUpgrade_ShellSewing sewing;
        TopoDS_Shape sewedShell = sewing.ApplySewing(newShell);
        assert(sewedShell.ShapeType() == TopAbs_SHELL);
        TopoDS_Shell newerShell = TopoDS::Shell(sewedShell);
        
        BRepCheck_Shell checkShell(newerShell);
        BRepCheck_Status status = checkShell.Closed();
        assert(status == BRepCheck_NoError);
        status = checkShell.Orientation();
        assert(status == BRepCheck_NoError);
        
        shells.push_back(newerShell);
        
        cout << "Shell " << nn << " has " << shellFaces[nn].size() << " faces.\n";
    }
    
    cout << "\n---- Build solid:\n";
    
    TopoDS_Solid newSolid;
    BRep_Builder builder;
    builder.MakeSolid(newSolid);
    
    for (int nn = 0; nn < shells.size(); nn++)
    {
        builder.Add(newSolid, shells[nn]);
    }
    
    BRepTools::Dump(newSolid, cout);

    cout << "\n---- Write solid to STP file:\n";

    STEPControl_Writer stepWriter;
    
    // 0: writes STEP files without assemblies
    // 1: writes all shapes in the form of assemblies
    // 2: writes shapes having a structure of (possibly nested) TopoDS_Compounds in the form of
    //    STEP assemblies; single shapes are written without assembly structures.
    assert(Interface_Static::SetIVal("write.step.assembly", 2));
    
//    assert(Interface_Static::SetCVal("write.step.schema","1"));
    
    stepWriter.Transfer(newSolid, STEPControl_AsIs);
}
*/

