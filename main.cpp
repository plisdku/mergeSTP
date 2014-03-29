#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <map>
#include <algorithm>

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/OFF_to_nef_3.h>
#include <CGAL/Segment_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel Kernel;
typedef CGAL::Nef_polyhedron_3<Kernel> NefPolyhedron;
typedef NefPolyhedron::Point_3 Point_3;
typedef NefPolyhedron::Plane_3 Plane_3;
typedef NefPolyhedron::Traits Traits;
typedef NefPolyhedron::Aff_transformation_3 Aff_transformation_3;
typedef CGAL::Polyhedron_3<Traits> Polyhedron;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
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
#include <ShapeFix_Shape.hxx>

// things to check orientations
#include <BRepClass3d_SolidClassifier.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>

#include <StlAPI.hxx>
#include <StlAPI_Reader.hxx>

#include <STEPControl_Writer.hxx>
#include <Interface_Static.hxx>

#include "ShellMapper.h"
#include "LoadSTL.h"

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
    set<Plane_3, PlaneLT> & accumPlanes);

TopoDS_Compound mergeSolids(const vector<vector<Shell> > & solids,
    const set<Point_3> & uniquePoints,
    const set<Line_3, LineLT> & uniqueLines,
    const set<Plane_3, PlaneLT> & uniquePlanes);

NefPolyhedron loadMultiOFF(istream & instr);
NefPolyhedron loadSTL(ifstream & instr);
NefPolyhedron loadOFF(ifstream & instr);

void parseInput(int argc, char const** argv,
    std::vector<std::string> & outFileNames, double & outScaleFactor, 
    std::string & outUnitString);

string sFileExtension(const string & fname)
{
    int dotPosition = fname.find_last_of(".");
    
    if (dotPosition == std::string::npos)
        throw(std::runtime_error("File has no extension"));
    
    return fname.substr(dotPosition+1);
}

int main(int argc, char const** argv)
{
    std::vector<std::string> inputFileNames;
    
    double scaleFactor = 1.0;
    std::string unitString = "MM";
    
    try
    {
        parseInput(argc, argv, inputFileNames, scaleFactor, unitString);
    }
    catch (const std::exception & ex)
    {
        cerr << "Error parsing input: " << ex.what() << "\n";
        printHelp();
        return 1;
    }
    
    int numInputFiles = inputFileNames.size();
    if (numInputFiles == 0)
    {
        printHelp();
        return 0;
    }
    
    set<Point_3> uniquePoints;
    set<Line_3, LineLT> uniqueLines;
    set<Plane_3, PlaneLT> uniquePlanes;
    vector<vector<Shell> > solids;
    
    Aff_transformation_3 scalingTx(CGAL::SCALING, scaleFactor);
    
    for (int inputFileNum = 0; inputFileNum < numInputFiles; inputFileNum++)
    {
        const string & fileName = inputFileNames[inputFileNum];
        string fileExt;
        
        NefPolyhedron nef;
        vector<Shell> shells;
        
        bool loaded = true;
        
        try
        {
            ifstream inStream(fileName.c_str());
            
            if (!inStream)
                throw(std::runtime_error("Stream is not 'true'."));
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
                loaded = false;
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
            nef.transform(scalingTx); // as-needed to convert units!!
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
    assert(Interface_Static::SetCVal("write.step.unit", unitString.c_str()));
    
    stepWriter.Transfer(allSolids, STEPControl_AsIs);
    stepWriter.Write("outStep.step");
    
    return 0;
}

void parseInput(int argc, char const** argv,
    std::vector<std::string> & outFileNames, double & outScaleFactor,
    std::string & outUnitString)
{
    for (int nn = 1; nn < argc; nn++)
    {
        string token(argv[nn]);
        
        if (token == "-unit")
        {
            if (nn+1 >= argc)
                throw(std::runtime_error("No unit given after -unit"));
            
            token = argv[++nn];
            
            if (token == "nm")
            {
                outScaleFactor = 1e-6;
                outUnitString = "MM";
            }
            else if (token == "um")
            {
                outScaleFactor = 1e-3;
                outUnitString = "MM";
            }
            else if (token == "mm")
            {
                outScaleFactor = 1.0;
                outUnitString = "MM";
            }
            else if (token == "cm")
            {
                outScaleFactor = 1.0;
                outUnitString = "CM";
            }
            else if (token == "m")
            {
                outScaleFactor = 1.0;
                outUnitString = "M";
            }
            else if (token == "km")
            {
                outScaleFactor = 1.0;
                outUnitString = "KM";
            }
            else
                throw(std::runtime_error("Unrecognized unit string"));
        }
        else
        {
            outFileNames.push_back(token);
        }
    }
}

void printHelp()
{
    cout << "Usage:\n";
    cout << "\tmergeSTP [file1] [file2] ...\n";
    cout << "\tmergeSTP -unit nm [file1] [file2] ...\n";
    cout << "-unit determines how lengths in input files are interpreted.\n";
    cout << "Valid units are nm, um, mm, cm, m, km.\n";
    cout << "Default unit is mm (standard for STEP files).\n";
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
    Polyhedron polyhedron;
    
    LoadSTL::BuildSTL<Polyhedron::HalfedgeDS> builder(instr);
    polyhedron.delegate(builder);
    
    return NefPolyhedron(polyhedron);
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
    set<Plane_3, PlaneLT> & accumPlanes)
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
//        cerr << "Volume " << volNum << "\n";
        
        int shell = 0;
        
        ShellEntryConstIterator itr;
        for (itr = insideVolumes[volNum]->shells_begin();
            itr != insideVolumes[volNum]->shells_end(); itr++)
        {
            shellMapper.newShell();
            
//            cout << "\tShell " << shell++ << "\n";
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

    BRepCheck_Status status;
    
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
                    
                    assert(lineMap.count(contour.lines()[iVertex]));
                    assert(pointMap.count(contour.points()[iVertex]));
                    assert(pointMap.count(contour.points()[iNext]));
                    
                    const TopoDS_Vertex & p1 = pointMap.find(contour.points()[iVertex])->second;
                    const TopoDS_Vertex & p2 = pointMap.find(contour.points()[iNext])->second;
                    
//                    cout << "\tadding " << contour.points()[iVertex]
//                        << "[" << p1.HashCode(99999) << "]" << " to "
//                        << contour.points()[iNext]
//                        << "[" << p2.HashCode(99999) << "]" << "\n";
                        
                    const gp_Lin & line = lineMap.find(contour.lines()[iVertex])->second;
                    
                    // An edge's orientation depends on the orientation of its
                    // underlying curve.  It seems important to construct the
                    // edge so that it points "forward" parallel to the line
                    // it lies on.
                    //
                    // I didn't see a straightforward way to find the gp_Pnt
                    // of a TopoDS_Vertex, so I'm building the p1p2 vector from
                    // the CGAL points and translating that to an OpenCASCADE
                    // vector, then comparing to the OpenCASCADE line direction.
                    
                    const Point_3 & p1cgal = contour.points()[iVertex];
                    const Point_3 & p2cgal = contour.points()[iNext];
                    Kernel::Vector_3 p1p2cgal = Kernel::Segment_3(p1cgal, p2cgal).to_vector();
                    
                    gp_Vec p1p2(CGAL::to_double(p1p2cgal[0]),
                        CGAL::to_double(p1p2cgal[1]), 
                        CGAL::to_double(p1p2cgal[2]));
                    gp_Vec lineVec(line.Direction());
                    
                    TopoDS_Edge edge;
                    
                    // here's the step where I care about the orientations.
                    // there are two ways to build the edge.  choose wisely. :-)
                    if (p1p2.Dot(lineVec) > 0)
                        edge = BRepBuilderAPI_MakeEdge(line, p1, p2);
                    else
                    {
                        gp_Lin floppedLine = line.Reversed();
                        edge = BRepBuilderAPI_MakeEdge(floppedLine, p1, p2);
                    }
                    
//                    if (iContour > 0)
//                    {
//                        cout << "v" << iVertex << ": " << p1p2.Dot(lineVec)
//                            << " is " << p1cgal << " to " << p2cgal << "\n";
//                        BRepTools::Dump(edge, cout);
//                    }
                    
                    buildWire.Add(newWires[iContour], edge);
                }
                
                BRepCheck_Wire checkWire(newWires[iContour]);
                assert( (status = checkWire.Closed()) == BRepCheck_NoError);
                assert( (status = checkWire.Orientation(TopoDS_Face())) == BRepCheck_NoError);
                
//                if (iContour > 0)
//                {
//                    cout << "Inner wire dump:\n";
//                    BRepTools::Dump(newWires[iContour], cout);
//                    cout << "Reversed wire dump:\n";
//                    BRepTools::Dump(newWires[iContour].Reversed(), cout);
//                }
            }
            
            TopoDS_Face newFace;
            assert(planeMap.count(face.plane()));
            BRepBuilderAPI_MakeFace makeFace(planeMap.find(face.plane())->second,
                newWires[0]);
//            cout << "Plane is " << face.plane() << "\n";
            
            
            // ok look, i don't understand why an inner wire that goes clockwise 
            // is not automatically seen as a hole within an outer wire that goes 
            // counterclockwise, but reversing the inner ones manually seems necessary,
            // so here goes.  sigh.
            for (int nn = 1; nn < newWires.size(); nn++)
            {
//                makeFace.Add(TopoDS::Wire(newWires[nn].Reversed()));
                makeFace.Add(newWires[nn]);
                //BRepCheck_Wire checkWire(newWires[nn]);
                //assert( (status = checkWire.Orientation(makeFace)) == BRepCheck_NoError);
            }
            
            
            newFace = makeFace;
            BRepCheck_Face checkFace(newFace);

//            if (newWires.size() > 1)
//            {
//                BRepTools::Dump(newFace, cout);
//            }
            
            assert( (status = checkFace.IntersectWires()) == BRepCheck_NoError);
            assert( (status = checkFace.ClassifyWires()) == BRepCheck_NoError);
            
            status = checkFace.OrientationOfWires();
            if (status != BRepCheck_NoError)
            {
                cout << "A face needs fixing.\n";
//                cout << "Face that needs fixing:\n";
//                BRepTools::Dump(newFace, cout);
                
                BRepCheck_Face chk(newFace);
                status = chk.OrientationOfWires();
                ShapeFix_Face fixFace;
                fixFace.Init(newFace);
                fixFace.FixOrientation();
                
                newFace = fixFace.Face();
                status = BRepCheck_Face(newFace).OrientationOfWires();
                
//                cout << "Face that is fixed:\n";
//                BRepTools::Dump(newFace, cout);
            }
            assert( status == BRepCheck_NoError);
            
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
    
    // NOW for some diagnostics!  Copied from Cauchy Ding's forum post.
    BRepClass3d_SolidClassifier classify;
    classify.Load(newSolid);
    classify.PerformInfinitePoint(1.0e-4); // classify a distant point, in/out
    TopAbs_State state = classify.State();
    
    if (state == TopAbs_IN)
        cerr << "Warning: a solid is inside-out.\n";
//    else
//        cout << "Solid is oriented correctly.\n";
    
    // Next the volume diagnostic:
    
    GProp_GProps properties;
    BRepGProp::VolumeProperties(newSolid, properties);
    
//    cout << "Solid mass is " << properties.Mass() << ".\n";
    
    
    return newSolid;
}
