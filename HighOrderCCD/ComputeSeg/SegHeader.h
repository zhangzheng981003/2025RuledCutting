#ifndef SEGHEADER_H
#define SEGHEADER_H

#include "HighOrderCCD\CrossField\crossField.h"
#include "HighOrderCCD\CrossField\covering_space.h"
#include <random>
#include "HighOrderCCD\CrossField/LoopGen.h"
#include<math.h>
#include"HighOrderCCD/Readtxt.h"
#include"HighOrderCCD/GraphAlgorithm.h"
#include "HighOrderCCD/SetCover/SetCoverRun.h"
#include "HighOrderCCD/SplineFitting/SplineSurface.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <list>
#include <tuple> 
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>

#include <CGAL/Polyhedron_3.h>
#include <Eigen/Dense>
#include "HighOrderCCD/SplineFitting/Mesh/MeshDefinition.h"
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include <CGAL/algorithm.h>

#include <TopoDS_Face.hxx>
#include <TopoDS_Edge.hxx>
#include <BRep_Tool.hxx>
#include <Geom_Plane.hxx>
#include <Geom_CylindricalSurface.hxx>
#include <Geom_ConicalSurface.hxx>
#include <Geom_SphericalSurface.hxx>
#include <Geom_ToroidalSurface.hxx>
#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_SurfaceOfRevolution.hxx>
#include <Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_OffsetSurface.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <BRep_Builder.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS.hxx>
#include <XBRepMesh.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Poly_Triangulation.hxx>
#include <Poly_Array1OfTriangle.hxx>
#include <Standard.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <Poly_Polygon3D.hxx>
#include <IntPolyh_Intersection.hxx>

using namespace std;
using namespace Eigen;


using K = typename CGAL::Simple_cartesian<double>;
using FT = typename K::FT;
using Ray = typename K::Ray_3;
using Line = typename K::Line_3;
using Segment = typename K::Segment_3;
using Point = typename K::Point_3;
using Triangle = typename K::Triangle_3;
typedef K::Plane_3 Plane;
typedef CGAL::Surface_mesh<Point> PolygonMesh;

typedef Eigen::MatrixXd Data;
typedef std::tuple< int, std::pair<double, double>, Eigen::MatrixXd >  Node;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> inner_derivative_t;//3*(order_num+1)
typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
typedef Eigen::Matrix<inner_scalar_t, Eigen::Dynamic, 1> derivative_t;
typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Vec12;
typedef Eigen::Matrix<scalar_t, 1, 3> Vec3;

/*
struct MyTriangle : public Triangle
{
public:
	MyTriangle(Point& p0, Point& p1, Point& p2, int _index, Vec3& _face_normal)
		: Triangle(p0, p1, p2)
	{
		index = _index;
		my_face_normal = _face_normal;
	}
	int index;
	Vec3 my_face_normal;
};

typedef std::vector<MyTriangle>::iterator MyIterator;
typedef CGAL::AABB_triangle_primitive<K, MyIterator> MyPrimitive;//1.
typedef CGAL::AABB_traits<K, MyPrimitive> MyAABB_triangle_traits;
typedef CGAL::AABB_tree<MyAABB_triangle_traits> MyTree;

typedef std::vector<Triangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;//1.
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;

*/

#endif