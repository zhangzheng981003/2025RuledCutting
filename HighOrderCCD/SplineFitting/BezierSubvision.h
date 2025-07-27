#pragma once
#include <iostream>
#include <list>
#include <vector>
#include <tuple> 
#include <cmath>
#include "Mesh/MeshDefinition.h"
#include "CommonFunction.h"
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>

using namespace std;
using namespace Eigen;
class BezierSubvision
{
public:
	BezierSubvision(Mesh& mesh);
	typedef Eigen::Matrix<double, 4, 4> Mat_4_4;
	typedef Eigen::MatrixXd Data;
	typedef std::vector< std::tuple< int, std::pair<double, double>, Eigen::MatrixXd > > Tree;
	typedef std::tuple< int, std::pair<double, double>, Eigen::MatrixXd >  Node;
	typedef Eigen::Matrix<double, Eigen::Dynamic, 1> inner_derivative_t;//3*(order_num+1)
	typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
	typedef Eigen::Matrix<inner_scalar_t, Eigen::Dynamic, 1> derivative_t;
	typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
	typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Vec12;
	typedef Eigen::Matrix<scalar_t, 1, 3> Vec3;
	
	Eigen::MatrixXd get_ori_two_sub_bezier_curve_mat(int sub_num, int n_order_num);
	vector<Eigen::MatrixXd> get_all_sub_bezier_curve_mat(int sub_num, int n_order_num);

public:
	Mesh mesh;


};