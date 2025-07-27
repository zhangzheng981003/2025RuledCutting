#include "BezierSubvision.h"


BezierSubvision::BezierSubvision(Mesh& mesh):mesh(mesh)
{

}


Eigen::MatrixXd BezierSubvision::get_ori_two_sub_bezier_curve_mat(int sub_num,int n_order_num)
{
	CommonFunction cf;
	auto mat_blocks = cf.total_subvision(sub_num, n_order_num);
	Eigen::MatrixXd big_mat;
	big_mat.resize((n_order_num + 1) * mat_blocks.size(), n_order_num + 1);
	big_mat.setZero();
	for (int i = 0; i < mat_blocks.size(); i++)
	{
		big_mat.block((n_order_num + 1) * i, 0, (n_order_num + 1), (n_order_num + 1)) = mat_blocks[i];
	}

	return big_mat;
}

// (1-t)A+tB
vector<Eigen::MatrixXd> BezierSubvision::get_all_sub_bezier_curve_mat(int sub_num, int n_order_num)
{
	Eigen::MatrixXd big_mat = get_ori_two_sub_bezier_curve_mat(sub_num, n_order_num);
	double gap = 1.0 / double(pow(2, sub_num));
	vector<double> t_value;
	t_value.clear();
	t_value.resize(pow(2, sub_num));
	for (int i=0; i<t_value.size();i++)
	{
		t_value[i] = i * gap;
	}
	t_value.push_back(1);

	Eigen::MatrixXd Zeros_mat;
	Zeros_mat.resize(big_mat.rows(), (n_order_num + 1));
	Zeros_mat.setZero();

	Eigen::MatrixXd BindCols_A(big_mat.rows(), big_mat.cols() + Zeros_mat.cols());
	BindCols_A << big_mat, Zeros_mat;

	Eigen::MatrixXd BindCols_B(big_mat.rows(), big_mat.cols() + Zeros_mat.cols());
	BindCols_B << Zeros_mat, big_mat;
	
	vector<Eigen::MatrixXd> all_mat_;
	all_mat_.clear();
	
	for (int i = 0; i < t_value.size(); i++)
	{
		Eigen::MatrixXd C_ = BindCols_A * (1 - t_value[i]) + BindCols_B * t_value[i];
		all_mat_.push_back(C_);
	}

	return all_mat_;
}



