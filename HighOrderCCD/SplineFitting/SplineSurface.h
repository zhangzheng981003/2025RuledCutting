#ifndef SPLINESURFACE_H
#define SPLINESURFACE_H

#include <iostream>
#include <list>
#include <vector>
#include <tuple> 
#include <cmath>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>
#include "Mesh/MeshDefinition.h"
//#include "Gradient.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/squared_distance_3.h>
#include "BezierSubvision.h"

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

#include <CGAL/algorithm.h>



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

typedef std::vector<MyTriangle>::iterator Iterator;
typedef CGAL::AABB_triangle_primitive<K, Iterator> Primitive;//1.
typedef CGAL::AABB_traits<K, Primitive> AABB_triangle_traits;
typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;


typedef std::vector<Triangle>::iterator Iterator_;
typedef CGAL::AABB_triangle_primitive<K, Iterator_> Primitive_;//1.
typedef CGAL::AABB_traits<K, Primitive_> AABB_triangle_traits_;
typedef CGAL::AABB_tree<AABB_triangle_traits_> Tree_;
typedef boost::optional<Tree_::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef Eigen::MatrixXd Data;
/*
typedef boost::optional<Tree::Intersection_and_primitive_id<Ray>::Type> Ray_intersection;
typedef boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
*/

class SplineSurface
{
public:

	struct Bezier_info
	{
		double u_cord;
		double v_cord;
		Eigen::Vector3d bezier_point;
		Vec3 mesh_point;
		Vec3 face_normal;
		int u_idx;
		int v_idx;
		bool is_intersection = 1;
	};

	struct Forward_Projection
	{
		double u_cord;
		double v_cord;
		Vec3 input_point;
		int input_idx;
		Vec3 output_point;
		Vec3 normal_on_output_surface;
		int is_orth;
		double avrdis;
		double mindis;
		double maxdis;
		double initdis;
	};

	SplineSurface(Mesh& mesh);
	void init_ready_data(Mesh& mesh);
	void uv_sample(int row_num);
	int combination(int n, int m);

	std::vector<int> is_or;
	double calc_bernstein_basis(int n, int m, double value);
	void calc_data_term_energy_derivative(VectorXd& spline, scalar_t& data_value, scalar_t& smooth_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, Eigen::VectorXd& decrease_direction,Tree& tree,Tree_& tree_ );
	pair<double,double> calc_data_term_energy(VectorXd& spline);
	pair<double, double> gcl_calc_data_term_energy(const Data& spline);


	double line_search(Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction,  Eigen::MatrixXd& spline1, double c_parm, double beta, double data_term, double smooth_term);
	void optimization(Mesh& mesh, Eigen::MatrixXd& spline);// vector<vector<Eigen::Vector3d>>& control_points);
	std::vector<double> generate_equally_spaced_vector(int num_elements);
	void bezier_visualization(int row_num, vector<vector<Eigen::Vector3d>>& control_points, int iter_num, int bezier_num);
	Eigen::Vector3d calc_bezier_surface_value(int u_degree, int v_degree, double u_cord, double v_cord, vector<vector<Eigen::Vector3d>>& control_points);
	double spline_ray_intersection(double x_cord, double y_cord, Tree_& tree_);
	void test();

	//***********IPC OPT
	double ipc_line_search(Mesh& mesh, vector<vector<Eigen::MatrixXd>>& all_bezier, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction, Eigen::VectorXd& barrier_grad, scalar_t& barrier_function, VectorXd& spline,
		double c_parm, double beta, double data_term, Tree_& tree_,
		vector<Eigen::RowVector3d>& P_input_set, vector<vector<Eigen::RowVector3d>>& Edge_input_set, vector<vector<Eigen::RowVector3d>>& face_input_set);
	void ipc_total_derivative(Mesh& mesh, VectorXd& spline, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, scalar_t& barrier_value, Eigen::VectorXd& barrier_grad, Eigen::MatrixXd& barrier_hessian, Eigen::VectorXd& decrease_direction, Tree& tree, Tree_& tree_,
		vector<vector<Eigen::MatrixXd>> all_bezier,
		vector<Eigen::RowVector3d>& P_input_set, vector<vector<Eigen::RowVector3d>>& Edge_input_set, vector<vector<Eigen::RowVector3d>>& face_input_set);
	void ipc_optimization(Mesh& mesh, vector<vector<Eigen::Vector3d>>& control_points);
	double step_CCD_check(VectorXd& spline, VectorXd& new_spline, vector<vector<Eigen::MatrixXd>>& all_bezier);
	double calc_dist_convex_hull_and_input_mesh(vector<Point>& points);
	void bezier_surface_subvision(Mesh& mesh, VectorXd& spline, vector<vector<Eigen::MatrixXd>>& all_bezier);

	//************正向投影
	void knot_vectors();
	pair<Eigen::MatrixXd, Eigen::MatrixXd> calc_uv_derivative(double u_, double v_);
	void calc_smooth_term_derivative(VectorXd& spline, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian);
	Eigen::Vector3d calc_bezier_surface_normal(int u_degree, int v_degree, double u_cord, double v_cord, vector<vector<Eigen::Vector3d>>& control_points);
	void calc_data_smooth_term_derivative(VectorXd& spline,scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian);
	void gcl_calc_data_smooth_term_derivative(const Data& spline, scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian);
	//void calc_smooth_term(VectorXd& spline, scalar_t& smooth_value);
	void gcl_bezier_visualization(int row_num, Data& spline, int iter_num, int bezier_num);


	std::set<int> search_projection_area(Data& spline, int sample_num);
	int spline_ray_intersection_(Eigen::Vector3d start_point, Eigen::Vector3d end_point);



	void fitting_init(int sample_num);



public:
	Mesh mesh;
	vector<Bezier_info> bi_vec;
	int U_Degree = 3;
	int V_Degree = 1;

	vector<Point> input_mesh_point_init;
	vector<Segment> input_mesh_line_init;
	vector<Triangle> input_mesh_triangle_init;

	double d0 = 1e-3;
	double x0 = 1e-3;
	double convex_hull_diameter_threshold = 0.2;
	double lamda = 1;
	double diagonal_length = 1.24;
	VectorXd A_vec;
	int ipc_search_step = 0;
	int line_search_num = 0;

	std::vector<MyTriangle> cgal_triangles;
	std::vector<Triangle> cgal_triangles_;
	double smooth_weight = 0;

	vector<Forward_Projection> all_data_forward_projection;
	vector<double> u_knot;
	vector<double> v_knot;
	Eigen::MatrixXd knot_weight;

	int calc_iter_num = 0;
	
	Tree tree;
	Tree_ tree_;

	Tree_ ordinary_tree;

	int smooth_model =2;
	double u_smooth_weight = 0.005;
	double v_smooth_weight = 0.2;

	double data_mesh2surface_weight = 0.2;
	double data_surface2mesh_weight = 0.005;
	

	vector<vector<Vec3>> face_point_temp;
	vector<int> cover_area_temp;

};

#endif