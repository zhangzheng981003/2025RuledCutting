#ifndef BSPLINEFITTING_H
#define BSPLINEFITTING_H
#include "HighOrderCCD/Config/Config.h"
#include <iostream>
#include <list>
#include <vector>
#include <tuple> 
#include <cmath>

#include "SplineFitting/Mesh/MeshDefinition.h"
#include "HighOrderCCD/SplineFitting/SplineSurface.h"

// Open CASCADE
#include <Geom_BSplineSurface.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array2OfReal.hxx>
#include <gp_Pnt.hxx>
#include <GeomFill_BSplineCurves.hxx> 
#include <Geom_BSplineCurve.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <gp_Pnt.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>

using namespace std;
using namespace Eigen;
extern string result_path_fitting1;
class BSplineFitting
{
public:
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

	//储存mesh到曲面的拟合点(投影点)
	struct Mesh2Surface
	{
		double u;
		double v;
		Vec3 point_on_surface;
		Vec3 point_on_mesh;
		bool is_found = true;
		vector<double> u_basis_temp;
		int u_interval_order;
		OpenMesh::Vec3d normal;
		double energy;
	};

	//储存曲面到mesh的拟合点(最近点)
	struct Surface2Mesh
	{
		double u;
		double v;
		Vec3 point_on_surface;
		Vec3 point_on_mesh;
	};

	BSplineFitting(Mesh& mesh);

	//初始化
	void sample_fitting_point_on_face(int sample_num);
	void BSpline_initialization();

	//计算BSpline值
	//void calc_BSpline_function_value();
	//void calc_BSpline_curve();
	//void calc_BSpline_surface();

	//构造BSpline曲线和曲面
	std::vector<Standard_Real> generate_equally_spaced_vector(int num_elements, double up_bound);
	void create_BSpline_curve(vector<CtrInfo>& ctr_point, Handle(Geom_BSplineCurve)& curve);
	void create_BSpline_surface(vector<CtrInfo>& ctr_point1, vector<CtrInfo>& ctr_point2, Handle(Geom_BSplineSurface)& surface);

	//BSpline曲面可视化
	void BSpline_surface_viewer(Handle(Geom_BSplineSurface)& surface, int sample_num);
	void BSpline_surface_viewer_2(const Data& spline, int u_sample_num, int v_sample_num, int out_number, int it_num);
	//测试
	void run();

	//计算BSpline的基函数值
	int u_k_position(double u, vector<double>& u_knot);
	void calc_basis_fuction(double u, int k, int poly_degree, std::vector<double>& basis_func);

	//优化 
	//1.计算导数
	void calc_data_term_grad_hessien(Eigen::VectorXd& BSpline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian);
	void calc_smooth_term_grad_hessien(VectorXd& BSpline, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian);
	void calc_objective_grad_hessien(const Data& spline, scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian);

	//2.计算能量
	double calc_data_term_energy(const Data& spline);
	double calc_smooth_term_energyn(const Data& spline);
	void calc_objective_energy(const Data& spline);

	//3.线搜索
	double line_search(VectorXd& BSpline, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction, double c_parm, double beta, double data_term);

	//4.优化
	void optimization(Mesh& mesh, Eigen::VectorXd& BSpline);

	void test();
	//5.细分
	void divide();
	//6.增加控制点
	void add_ctr_point(Eigen::VectorXd& BSpline, double u);

	void CCD_initialization(Data& spline);
	std::vector<double> uni_para(int samp_num);

public:
	Mesh mesh;
	Standard_Integer poly_degree = 3;
	int control_number = 0;

	Tree my_tree;//带面索引：只能用最近点，不能用射线求交
	Tree_ tree;//不带面索引
	std::vector<MyTriangle> my_cgal_triangles;
	std::vector<Triangle> cgal_triangles;

	vector<vector<Vec3>> face_point_temp;
	vector<OpenMesh::Vec3d> face_normal_temp;
	vector<int> cover_area_temp;

	vector<Mesh2Surface> mesh2surface_temp;
	vector<Surface2Mesh> surface2mesh_temp;

	//vector<double> u_knot_temp ;
	int is_add = 0;
	vector<int> is_change;
	Eigen::MatrixXd control_points_;
	int p_ = 3;
	double interval_ = 1;
	int flag123 = 0;
	int n_, m_;
	int is_good = 1;
	Eigen::VectorXd u_;
	vector<double> u_knot_temp, v_knot_temp;
	vector<double> e_intervals;
	vector<double> e_len;
	vector<double> avr_e_intervals;
	vector<int> e_num;
	double interval_up_bound = 1;
	int insert1_idx;
	vector<int> count_divide;
	double avr_energy = 10.0;
	double max_energy = 0.0;
	double min_energy = 10.0;
	double data_mesh2surface_weight = 1;
	double data_surface2mesh_weight = 0;
	int is_neg = 1;
	int zero_num = 0;
	int divide_num = 0;
	vector<vector<Eigen::VectorXd>> divide_ctrpoint;
	vector<vector<vector<double>>> divide_u_knot;
	std::vector<std::vector<Eigen::MatrixXd>> divide_basis;
	int change_idx;
	int iter_num = 0;
	vector<int> interested_area;
	double u_smooth_weight = 1;
	double v_smooth_weight = 1;

	int weight_order = 0;
	vector<double> m2s_weight;
	vector<double> s2m_weight;

};

#endif