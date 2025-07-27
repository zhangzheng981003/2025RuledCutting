#ifndef CURVEFITTING_H
#define CURVEFITTING_H


#include "SplineInit.h"
#include <iostream>
#include <Standard.hxx>
#include <gp.hxx>
#include <Geom_BSplineSurface.hxx>
#include <GeomLProp_SLProps.hxx>


class CurveFitting
{
public:
	CurveFitting(Mesh& mesh);

	struct PointCloud
	{
		int index = 0;
		Vector3d position;
		Vec3 autodiff_position;
		Standard_Real u = 0;
	};

	struct BezierPoint
	{
		int index = 0;
		double u = 0;
		Vector3d position;
		Vector3d rulling_direction;
		Vector3d tangent_direction;
		Vector3d normal_direction;
	};

	int combination(int n, int m);
	void initialization();
	void curve_init(int u_degree, VectorXd& curve, vector<Vector3d>& point_cloud);
	void curve_fitting(int degree, VectorXd& curve, vector<Vector3d>& point_cloud);
	void calc_projection_distance2_gradient(int u_degree, VectorXd& spline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, Eigen::VectorXd& decrease_direction);
	double calc_energy(int u_degree, VectorXd& spline);
	double line_search(int u_degree, VectorXd& spline, double data_value, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction);
	void energy_opt(int u_degree, VectorXd& spline);
	std::vector<double> generate_equally_spaced_vector(int num_elements);
	void calc_soomth_surface(VectorXd& curve, int sample_num);
public:
	Mesh mesh;
	vector<PointCloud> point_cloud_temp;
	vector<BezierPoint> bezier_point_temp;
	double c_parm = 1e-4;
	double beta = 0.5;
	double stop_threshold = 1e-6;
	double data_energy = 0.0;
	vector<Triangle> cgal_triangles;
	vector<MyTriangle> my_cgal_triangles;
	Tree my_tree;
	Tree_ tree;
	vector<SplineInit::SampleCurve> sample_curve_temp;


	SpMat m_A, m_L, m_P;
	Eigen::SimplicialLDLT<SpMat> solver;
	Eigen::VectorXd m_b, m_x, m_B;
	double rulings_length = 0.0;
	double average_edge_length;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
};


#endif