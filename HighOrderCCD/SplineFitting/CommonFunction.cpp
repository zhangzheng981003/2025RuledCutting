#include "CommonFunction.h"

CommonFunction::CommonFunction()
{

}

vector<Eigen::MatrixXd> CommonFunction::total_subvision(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> mat_blocks;
	switch (var_num)
	{
		case 0:mat_blocks = sub0(var_num,  n_order_num); break;
		case 1:mat_blocks = sub1(var_num,  n_order_num); break;
		case 2:mat_blocks = sub2(var_num,  n_order_num); break;
		case 3:mat_blocks = sub3(var_num,  n_order_num); break;
		case 4:mat_blocks = sub4(var_num,  n_order_num); break;
		case 5:mat_blocks = sub5(var_num,  n_order_num); break;
		case 6:mat_blocks = sub6(var_num,  n_order_num); break;
		case 7:mat_blocks = sub7(var_num,  n_order_num); break;
		case 8:mat_blocks = sub8(var_num,  n_order_num); break;
		case 9:mat_blocks = sub9(var_num,  n_order_num); break;
		case 10:mat_blocks = sub10(var_num,n_order_num); break;
		default:cout << endl;
				cout << "The number of variables exceeds 10." << endl;
				abort();
				exit(EXIT_FAILURE);;
	}
	return mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub0(int var_num, int n_order_num)
{
	
	Eigen::MatrixXd E;
	E.setIdentity(n_order_num + 1, n_order_num + 1);
	vector<Eigen::MatrixXd> ori = { E };
	return ori;
}


vector<Eigen::MatrixXd> CommonFunction::sub1(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	return ori;
}

vector<Eigen::MatrixXd> CommonFunction::sub2(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			Eigen::MatrixXd C_ = ori[i] * ori[j];
			new_mat_blocks.push_back(C_);
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub3(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k];
				new_mat_blocks.push_back(C_);
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub4(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4];
					new_mat_blocks.push_back(C_);
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub5(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5];
						new_mat_blocks.push_back(C_);
					}
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub6(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						for (int k6 = 0; k6 < ori.size(); k6++)
						{
							Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5] * ori[k6];
							new_mat_blocks.push_back(C_);
						}
					}
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub7(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						for (int k6 = 0; k6 < ori.size(); k6++)
						{
							for (int k7 = 0; k7 < ori.size(); k7++)
							{
								Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5] * ori[k6] * ori[k7];
								new_mat_blocks.push_back(C_);
							}
						}
					}
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub8(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						for (int k6 = 0; k6 < ori.size(); k6++)
						{
							for (int k7 = 0; k7 < ori.size(); k7++)
							{
								for (int k8 = 0; k8 < ori.size(); k8++)
								{
									Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5] * ori[k6] * ori[k7] * ori[k8];
									new_mat_blocks.push_back(C_);
								}
							}
						}
					}
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::sub9(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						for (int k6 = 0; k6 < ori.size(); k6++)
						{
							for (int k7 = 0; k7 < ori.size(); k7++)
							{
								for (int k8 = 0; k8 < ori.size(); k8++)
								{
									for (int k9 = 0; k9 < ori.size(); k9++)
									{
										Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5] * ori[k6] * ori[k7] * ori[k8] * ori[k9];
										new_mat_blocks.push_back(C_);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return new_mat_blocks;
}


vector<Eigen::MatrixXd> CommonFunction::sub10(int var_num, int n_order_num)
{
	vector<Eigen::MatrixXd> ori = n_degree_bezier_curve_subdivide(n_order_num);
	vector<Eigen::MatrixXd> new_mat_blocks;
	new_mat_blocks.clear();
	for (int i = 0; i < ori.size(); i++)
	{
		for (int j = 0; j < ori.size(); j++)
		{
			for (int k = 0; k < ori.size(); k++)
			{
				for (int k4 = 0; k4 < ori.size(); k4++)
				{
					for (int k5 = 0; k5 < ori.size(); k5++)
					{
						for (int k6 = 0; k6 < ori.size(); k6++)
						{
							for (int k7 = 0; k7 < ori.size(); k7++)
							{
								for (int k8 = 0; k8 < ori.size(); k8++)
								{
									for (int k9 = 0; k9 < ori.size(); k9++)
									{
										for (int k10 = 0; k10 < ori.size(); k10++)
										{
											Eigen::MatrixXd C_ = ori[i] * ori[j] * ori[k] * ori[k4] * ori[k5] * ori[k6] * ori[k7] * ori[k8] * ori[k9] * ori[k10];
											new_mat_blocks.push_back(C_);
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return new_mat_blocks;
}

vector<Eigen::MatrixXd> CommonFunction::n_degree_bezier_curve_subdivide(int n_order_num)
{
	Eigen::MatrixXd basis1, basis2;
	basis1.setZero(n_order_num + 1, n_order_num + 1);
	basis2.setZero(n_order_num + 1, n_order_num + 1);

	for (int i = 0; i < n_order_num + 1; i++)
	{
		for (int j = 0; j < i + 1; j++)
		{
			if (j == 0 || j == i)
			{
				basis1(i, j) = std::pow(0.5, i);
				basis2(n_order_num - i, n_order_num - j) = std::pow(0.5, i);
			}
			else
			{
				basis1(i, j) = 0.5 * basis1(i - 1, j - 1) + 0.5 * basis1(i - 1, j);
				basis2(n_order_num - i, n_order_num - j) = basis1(i, j);
			}	
		}
	}

	vector<Eigen::MatrixXd> mat;
	mat.resize(2);
	mat[0] = basis1; mat[1] = basis2;

	return mat;
}




















