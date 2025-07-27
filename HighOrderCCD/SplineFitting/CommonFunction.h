#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <Eigen/Dense>

using namespace std;
typedef Eigen::Matrix<double, 4, 4> Mat_4_4;
class CommonFunction
{
public:
	CommonFunction();
	
	vector<Eigen::MatrixXd> total_subvision(int num, int n_order_num);
	vector<Eigen::MatrixXd> sub0(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub1(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub2(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub3(int var_numf, int n_order_num);
	vector<Eigen::MatrixXd> sub4(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub5(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub6(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub7(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub8(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub9(int var_num, int n_order_num);
	vector<Eigen::MatrixXd> sub10(int var_num, int n_order_num);


	vector<Eigen::MatrixXd> n_degree_bezier_curve_subdivide(int n_order_num);


public:

private:

};