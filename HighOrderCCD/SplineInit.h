#ifndef SPLINEINIT_H
#define SPLINEINIT_H

#include <HighOrderCCD/SplineFitting/SplineSurface.h>
#include <Eigen/Eigenvalues>
#include "HighOrderCCD/ComputeSeg/SegSurface.h"

#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <Geom_BezierCurve.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <Eigen/Eigenvalues>

#include <CGAL/Projection_traits_yz_3.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/algorithm.h>

#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <gp_Pnt.hxx>
#include <GeomFill_BSplineCurves.hxx> 
#include <Geom_BSplineCurve.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Pnt.hxx>

using namespace std;
using namespace Eigen;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> T;
typedef K::Vector_3 Vector_3;
typedef CGAL::Projection_traits_yz_3<K> Projection_traits;
typedef Projection_traits::Point_2 Projection_point;

class SplineInit
{
public:
	struct PathInfo
	{
		int rank;
		int idx;
		double edge_length;
		double arc_length_param;
		Vec3 input_point;
		Vector3d e_input_point;
		Standard_Real u;
		vector<double> u_basis_temp;
		int u_interval_order;
	};
	struct CtrInfo
	{
		gp_Pnt position;
		Standard_Real weight;
	};

	struct BSplinePoint
	{
		Standard_Real u;
		Standard_Real v;
		gp_Pnt position;
		Eigen::Vector3d bsp_position;
	};

	struct VInfo
	{
		int rank;
		int idx;
		Vector3d p;
		Vector3d edge_derection;
		Vector3d rulling_derection;
		Standard_Real u;
	};

	struct SampleCurve
	{
		Vector3d p;
		Vector3d tangent_vector;
		Vector3d normal_vector;
		Vector3d rulling_vector;
		double u;
		int u_interval_order;
		vector<double> u_basis_temp;
	};

	SplineInit(Mesh& mesh);
	void initialization(SegSurface& ss);

	double spline_basis(int spline_order, int cur_order, double arc_param);

	int u_k_position(double u, vector<double>& u_knot);
	void calc_basis_fuction(double u, int k, int poly_degree, std::vector<double>& basis_func);
	void add_ctr_point(Eigen::VectorXd& BSpline, double u);
	void calc_control_point(int spline_order, int pc_num);
	int combination(int n, int m);
	void sample_points();
	std::vector<double> generate_equally_spaced_vector(int num_elements);
	void SplineInit::create_BSpline_surface(vector<CtrInfo>& ctr_point1, vector<CtrInfo>& ctr_point2, Handle(Geom_BSplineSurface)& surface);
	void SplineInit::create_BSpline_curve(vector<CtrInfo>& ctr_point, Handle(Geom_BSplineCurve)& curve);
	void sample_init_spline(SegSurface::Path& path);
	void calc_projection_distance2_gradient(VectorXd& spline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, Eigen::VectorXd& decrease_direction);
	double calc_energy(VectorXd& spline);
	double line_search(VectorXd& spline, double data_value, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction);
	void energy_opt(VectorXd& spline, SegSurface::Path& path);
	std::vector<double> uni_para(int samp_num);
	void draw_point();
	void SplineInit::BSpline_surface_viewer_2(const Data& spline, int sample_num, int out_number);
	void draw_surface();
	MatrixXd calc_up_degree_matrix(int spline_order);

	void calc_pca_energy(VectorXd& curve, int sample_num, SegSurface& ss, SegSurface::Path& path);
	void calc_point_projection(Point& point, Plane& plane);
	void calc_spline_init_info(Vector3d& center, Vector3d& direction, Vector3d& ruling_direction, SegSurface& ss, SegSurface::Path& path, double min_translation_length);

	void calc_soomth_surface(VectorXd& curve, int sample_num);
	void calc_soomth_surface_pca_plane(int sample_num, SegSurface& ss, SegSurface::Path& path);
public:
	Mesh mesh;
	SpMat m_A,m_L, m_P;
	Eigen::SimplicialLDLT<SpMat> solver;
	Eigen::VectorXd m_b,m_x, m_B;
	Eigen::MatrixXd cross_field;
	
	vector<int> point_cloud;
	vector<PathInfo> loop;
	
	int spline_order = 3;
	int poly_degree = 3;//B样条阶数
	int init_ctr_num = 6;//初始控制点个数
	int insert1_idx;//插入点位置
	Eigen::MatrixXd control_points;
	VectorXd spline_init;
	double c_parm = 1e-6;
	double beta = 0.5;
	double stop_threshold = 1e-6;
	double data_vvv = 0.0;

	double rulings_length = 0.0;
	double average_edge_length;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	double bbox_mix_length = 0.0;
	vector<VInfo> vifo_temp;
	vector<vector<Vector3d>> cp_;
	vector<SampleCurve> sample_curve_temp;
	std::vector<MyTriangle> cgal_triangles;
	vector<Point> curve_pro_point;
	vector<pair<Point,Point>> rulings_pro_point;
	vector<double> u_knot_temp;
	vector<SampleCurve> sample_surface_temp;
};
#endif