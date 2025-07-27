#include "BSplineNoCollision.h"

BSplineNoCollision::BSplineNoCollision(Mesh& mesh) :mesh(mesh)
{
	//
}

BSplineNoCollision::~BSplineNoCollision()
{

}

void BSplineNoCollision::bspline_fitting_init()
{
	//
}

void BSplineNoCollision::bspline_fitting_init_plane(SegSurface::Path& one_path)
{
	SplineInit* si = new SplineInit(mesh);
	SegSurface* ss = new SegSurface(mesh);
	si->initialization(*ss);
	VectorXd spline;
	si->energy_opt(spline, one_path);
	ctr_num = spline.size() / 3;
	si->calc_soomth_surface(spline, 20);
	si->calc_soomth_surface_pca_plane(200, *ss, one_path);

	delete si, ss;
	si = nullptr;
	ss = nullptr;
}

BSplineNoCollision::ResInfo BSplineNoCollision::bspline_fitting_no_collision(const std::string& input_mesh_file, Mesh& input_mesh, int patch_num, SegSurface::Path& one_path)
{
	MeshTools::WriteMesh(input_mesh, "D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/inputmesh.obj");
	int zero_num = 0;
	int init_num = 0;
	vector<double> res_avre(5, 100.0);
	vector<double> res_mine(5, 100.0);
	vector<double> res_maxe(5, 100.0);
	vector<Data> res_spline;
	res_spline.resize(5);
	vector<vector<double>> res_knot;
	res_knot.resize(5);
	zero_num = 0;
	iscollision.clear();
	iscollision.resize(5);

	auto ruling_direction = one_path.ruling_direction;
	auto center_point = one_path.center_point;
	auto translation_direction = one_path.translation_direction;

	flag2:Data bspline_mat;
	bspline_mat.resize(ctr_num * 2, 3);
	Eigen::Vector3d center_trans = (one_path.center_point + one_path.translation_length * one_path.translation_direction);
	Eigen::Vector3d v0 = center_trans + one_path.ruling_direction + one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v1 = center_trans + one_path.ruling_direction - one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v2 = center_trans - one_path.ruling_direction + one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v3 = center_trans - one_path.ruling_direction - one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;

	bspline_mat.row(0) << v0(0), v0(1), v0(2);
	bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
	bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
	bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
	for (int i = 1; i < ctr_num - 1; i++)
	{
		bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
		bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
	}

	flag:VectorXd bspline;
	bspline.resize(bspline_mat.rows() * 3);
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		bspline[3 * i] = bspline_mat.coeffRef(i, 0);
		bspline[3 * i + 1] = bspline_mat.coeffRef(i, 1);
		bspline[3 * i + 2] = bspline_mat.coeffRef(i, 2);
	}

	// 2.normalize axis and aabb_axis
	int dim = axis.size();
	for (int k = 0; k < dim; k++)
	{
		axis[k].normalize();
	}

	dim = aabb_axis.size();
	for (int k = 0; k < dim; k++)
	{
		aabb_axis[k].normalize();
	}

	//3.设置节点向量
	BSplineFitting bsf(input_mesh);
	bsf.interested_area.resize(input_mesh.n_faces());
	for (int i = 0; i < input_mesh.n_faces(); i++)
	{
		bsf.interested_area[i] = 0;
	}

	for (int j = 0; j < one_path.cover_face.size(); j++)
	{
		bsf.interested_area[one_path.cover_face[j]] = 1;
	}


	bsf.CCD_initialization(bspline_mat);

	//4.初始化 细分的控制点及节点向量
	int max_sub_num = 13;
	divide_ctrpoint_.resize(max_sub_num);
	divide_u_knot_.resize(max_sub_num);
	divide_basis_.resize(max_sub_num);
	for (int i = 0; i < max_sub_num; i++)
	{
		divide_ctrpoint_[i].resize(pow(2, i));
		divide_u_knot_[i].resize(pow(2, i));
		divide_basis_[i].resize(pow(2, i));
	}

	divide_ctrpoint_[0][0] = bspline;
	divide_u_knot_[0][0] = bsf.u_knot_temp;
	divide_basis_[0][0].resize(bspline.size() / 6, bspline.size() / 6);
	divide_basis_[0][0].setIdentity();

	subdivide_tree.resize(1);
	std::pair<double, double> range(1, 1);
	Eigen::MatrixXd basis = divide_basis_[0][0];
	Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
	combined.setZero();
	combined.block(0, 0, basis.rows(), basis.cols()) = basis;
	combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
	subdivide_tree[0] = std::make_tuple(0, range, combined);

	stop_condition = 1;
	data_energy0 = 1;
	epsilon2 = 1e-10;
	max_step = 1.0;
	offset = 1e-10;
	margin = 1e-10;

	data_weight = 1;
	smooth_weight = 1;
	fast_barrier_weight = 1;

	subdivision_number = 0;
	barrier_lambda = 0.8;

	Eigen::MatrixXd V;// vn*3
	Eigen::MatrixXi F;// fn*3
	igl::read_triangle_mesh(input_mesh_file, V, F);
	std::cout << "before bvh init\n";
	BVH bvh;
	bvh.InitObstacle(V, F);

	int opt_num = 0;
	bsf.count_divide.resize(bspline_mat.size() / 6 - 4 + 1); //统计每个初始参数区间段的细分次数
	for (int i = 0; i < bsf.count_divide.size(); i++)
	{
		bsf.count_divide[i] = 0;
	}
	bsf.BSpline_surface_viewer_2(bspline_mat, 5,100, -1,patch_num);

	//储存初始B样条曲面
	ofstream ioszz("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/spline_0.txt");
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		ioszz << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz.close();
	ofstream ioszz1("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/knot_0.txt");
	for (int i = 0; i < bsf.u_knot_temp.size(); i++)
	{
		ioszz1 << bsf.u_knot_temp[i] << endl;
	}
	ioszz1.close();

	iter_ = 0;
	iscollision[zero_num] = 0;
	is_good_init = 1;

	while (iscollision[zero_num] == 0&&(is_good_init)&& bsf.iter_num < 10 && bsf.avr_energy>2e-4 && (bsf.iter_num < 4 || bsf.avr_energy < 0.05))
	{

		cout << "*** 迭代次数： " << opt_num << " , MAX_step: " << max_step << " ***" << endl;
		clock_t opt_time1 = clock();
		Optimization3D_time::gcl_optimization(bspline_mat, V, F, bvh, bsf);
		bsf.BSpline_surface_viewer_2(bspline_mat,5, 100, opt_num,patch_num);
		clock_t opt_time2 = clock();
		opt_num++;
		iter_++;
		bsf.iter_num++;
	}

	if (is_good_init!=1&& init_num<4&&zero_num==0)
	{
		init_num++;
		bspline_mat.resize(ctr_num * 2, 3);
		Matrix3d Rota;
		Rota.setZero();
	    ruling_direction = one_path.ruling_direction;
		//center_point = one_path.center_point;
		translation_direction = one_path.translation_direction;
		double c;
		double s;
		if (init_num == 1)//绕rullings顺时针旋转
		{
			c = cos(-20.0 / 180.0 * 3.1415926);
			s = sin(-20.0 / 180.0 * 3.1415926);
		}
		if (init_num == 2)
		{
			c = cos(20.0 / 180.0 * 3.1415926);
			s = sin(20.0 / 180.0 * 3.1415926);
		}
		if (init_num == 3)
		{
			c = cos(-40.0 / 180.0 * 3.1415926);
			s = sin(-40.0 / 180.0 * 3.1415926);
		}
		if (init_num == 4)
		{
			c = cos(40.0 / 180.0 * 3.1415926);
			s = sin(40.0 / 180.0 * 3.1415926);
		}
		Rota(0, 0) = c + (1 - c) * pow(ruling_direction[0], 2);
		Rota(0, 1) = (1 - c) * ruling_direction[0] * ruling_direction[1] - s * ruling_direction[2];
		Rota(0, 2) = (1 - c) * ruling_direction[0] * ruling_direction[2] + s * ruling_direction[1];
		Rota(1, 0) = (1 - c) * ruling_direction[0] * ruling_direction[1] + s * ruling_direction[2];
		Rota(1, 1) = c + (1 - c) * pow(ruling_direction[1], 2);
		Rota(1, 2) = (1 - c) * ruling_direction[1] * ruling_direction[2] - s * ruling_direction[0];
		Rota(2, 0) = (1 - c) * ruling_direction[0] * ruling_direction[2] - s * ruling_direction[1];
		Rota(2, 1) = (1 - c) * ruling_direction[1] * ruling_direction[2] + s * ruling_direction[0];
		Rota(2, 2) = c + (1 - c) * pow(ruling_direction[2], 2);
		center_point = Rota * center_point;
		translation_direction = Rota * translation_direction;


		if (init_num ==12321)
		{
			c = cos(-70.0 / 180.0 * 3.1415926);
			s = sin(-70.0 / 180.0 * 3.1415926);
			Rota(0, 0) = c + (1 - c) * pow(translation_direction[0], 2);
			Rota(0, 1) = (1 - c) * translation_direction[0] * translation_direction[1] - s * translation_direction[2];
			Rota(0, 2) = (1 - c) * translation_direction[0] * translation_direction[2] + s * translation_direction[1];
			Rota(1, 0) = (1 - c) * translation_direction[0] * translation_direction[1] + s * translation_direction[2];
			Rota(1, 1) = c + (1 - c) * pow(translation_direction[1], 2);
			Rota(1, 2) = (1 - c) * translation_direction[1] * translation_direction[2] - s * translation_direction[0];
			Rota(2, 0) = (1 - c) * translation_direction[0] * translation_direction[2] - s * translation_direction[1];
			Rota(2, 1) = (1 - c) * translation_direction[1] * translation_direction[2] + s * translation_direction[0];
			Rota(2, 2) = c + (1 - c) * pow(translation_direction[2], 2);
			ruling_direction = Rota * ruling_direction;
			auto dir = ruling_direction.cross(translation_direction);
			c = cos(30.0 / 180.0 * 3.1415926);
			s = sin(30.0 / 180.0 * 3.1415926);
			Rota(0, 0) = c + (1 - c) * pow(dir[0], 2);
			Rota(0, 1) = (1 - c) * dir[0] * dir[1] - s * dir[2];
			Rota(0, 2) = (1 - c) * dir[0] * dir[2] + s * dir[1];
			Rota(1, 0) = (1 - c) * dir[0] * dir[1] + s * dir[2];
			Rota(1, 1) = c + (1 - c) * pow(dir[1], 2);
			Rota(1, 2) = (1 - c) * dir[1] * dir[2] - s * dir[0];
			Rota(2, 0) = (1 - c) * dir[0] * dir[2] - s * dir[1];
			Rota(2, 1) = (1 - c) * dir[1] * dir[2] + s * dir[0];
			Rota(2, 2) = c + (1 - c) * pow(dir[2], 2);
			translation_direction = Rota * translation_direction;
			ruling_direction = dir.cross(translation_direction);
			c = cos(-10.0 / 180.0 * 3.1415926);
			s = sin(-10.0 / 180.0 * 3.1415926);
			Rota(0, 0) = c + (1 - c) * pow(ruling_direction[0], 2);
			Rota(0, 1) = (1 - c) * ruling_direction[0] * ruling_direction[1] - s * ruling_direction[2];
			Rota(0, 2) = (1 - c) * ruling_direction[0] * ruling_direction[2] + s * ruling_direction[1];
			Rota(1, 0) = (1 - c) * ruling_direction[0] * ruling_direction[1] + s * ruling_direction[2];
			Rota(1, 1) = c + (1 - c) * pow(ruling_direction[1], 2);
			Rota(1, 2) = (1 - c) * ruling_direction[1] * ruling_direction[2] - s * ruling_direction[0];
			Rota(2, 0) = (1 - c) * ruling_direction[0] * ruling_direction[2] - s * ruling_direction[1];
			Rota(2, 1) = (1 - c) * ruling_direction[1] * ruling_direction[2] + s * ruling_direction[0];
			Rota(2, 2) = c + (1 - c) * pow(ruling_direction[2], 2);
			//center_point = Rota * center_point;
			translation_direction = Rota * translation_direction;
		}


		cout << "init_translation_direction=" << one_path.translation_direction.norm() << endl;
		cout <<"初始导线长度：" <<ruling_direction.cross(one_path.translation_direction).norm() << endl;
		cout <<"translation_direction="<< translation_direction.norm() << endl;
		cout << "导线长度：" << ruling_direction.cross(translation_direction).norm() << endl;
		Eigen::Vector3d center_trans = (center_point + one_path.translation_length * translation_direction);
		Eigen::Vector3d v0 = center_trans + ruling_direction + ruling_direction.cross(translation_direction)* 0.5 * one_path.path_length;
		Eigen::Vector3d v1 = center_trans + ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
		Eigen::Vector3d v2 = center_trans - ruling_direction + ruling_direction.cross(translation_direction) * 0.5 * one_path.path_length;
		Eigen::Vector3d v3 = center_trans - ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;

		bspline_mat.row(0) << v0(0), v0(1), v0(2);
		bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
		bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
		bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
		for (int i = 1; i < ctr_num - 1; i++)
		{
			bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
			bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
		}
		goto flag;
	}
	if (is_good_init && init_num > 0&&zero_num==0)
	{
		one_path.ruling_direction  =  ruling_direction;
		one_path.center_point = center_point;
		one_path.translation_direction = translation_direction;
		init_num = 0;
	}

	if (iscollision[0] == 1)
	{
		iscollision[0] = 0;
		one_path.translation_length += 0.2;
		goto flag2;
	}

	if (bsf.avr_energy > 1.5 * 1e-3 && zero_num <= -1)
	{
		res_avre[zero_num] = bsf.avr_energy;
		res_maxe[zero_num] = bsf.max_energy;
		res_mine[zero_num] = bsf.min_energy;
		res_knot[zero_num] = bsf.u_knot_temp;
		res_spline[zero_num] = bspline_mat;
		zero_num++;
		if (zero_num <= 4)
		{
			bspline_mat.resize(ctr_num * 2, 3);
			Matrix3d Rota;
			Rota.setZero();
			ruling_direction = one_path.ruling_direction;
			//center_point = one_path.center_point;
			translation_direction = one_path.translation_direction;
			double c;
			double s;
			if (zero_num == 1)//绕rullings顺时针旋转
			{
				c = cos(-32.5 / 180.0 * 3.1415926);
				s = sin(-32.5 / 180.0 * 3.1415926);
			}
			if (zero_num == 2)
			{
				c = cos(22.5 / 180.0 * 3.1415926);
				s = sin(22.5 / 180.0 * 3.1415926);
			}
			if (zero_num == 3)
			{
				c = cos(-45.0 / 180.0 * 3.1415926);
				s = sin(-45.0 / 180.0 * 3.1415926);
			}
			if (zero_num == 4)
			{
				c = cos(45.0 / 180.0 * 3.1415926);
				s = sin(45.0 / 180.0 * 3.1415926);
			}
			Rota(0, 0) = c + (1 - c) * pow(ruling_direction[0], 2);
			Rota(0, 1) = (1 - c) * ruling_direction[0] * ruling_direction[1] - s * ruling_direction[2];
			Rota(0, 2) = (1 - c) * ruling_direction[0] * ruling_direction[2] + s * ruling_direction[1];
			Rota(1, 0) = (1 - c) * ruling_direction[0] * ruling_direction[1] + s * ruling_direction[2];
			Rota(1, 1) = c + (1 - c) * pow(ruling_direction[1], 2);
			Rota(1, 2) = (1 - c) * ruling_direction[1] * ruling_direction[2] - s * ruling_direction[0];
			Rota(2, 0) = (1 - c) * ruling_direction[0] * ruling_direction[2] - s * ruling_direction[1];
			Rota(2, 1) = (1 - c) * ruling_direction[1] * ruling_direction[2] + s * ruling_direction[0];
			Rota(2, 2) = c + (1 - c) * pow(ruling_direction[2], 2);
			center_point = Rota * center_point;
			translation_direction = Rota * translation_direction;
			translation_direction.normalize();
			Eigen::Vector3d center_trans = (center_point + one_path.translation_length * translation_direction);
			Eigen::Vector3d v0 = center_trans + ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v1 = center_trans + ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v2 = center_trans - ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v3 = center_trans - ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;

			bspline_mat.row(0) << v0(0), v0(1), v0(2);
			bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
			bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
			bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
			for (int i = 1; i < ctr_num - 1; i++)
			{
				bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
				bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
			}
			goto flag;
		}
	}
	double min_avre = 1e6;
	double min_mine = 1e6;
	double min_maxe = 1e6;
	if (zero_num >0)
	{
		
		for (int i = 0; i < zero_num; i++)
		{
			if (res_avre[i] < min_avre && iscollision[i] == 0)
			{
				min_avre = res_avre[i];
				min_mine = res_mine[i];
				min_maxe = res_maxe[i];
				bspline_mat = res_spline[i];
				bsf.u_knot_temp = res_knot[i];
			}
		}

	}



	//储存最终B样条曲面
	ofstream ioszz_("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/spline_1.txt");
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		ioszz_ << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz_.close();
	ofstream ioszz_1("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/knot_1.txt");
	for (int i = 0; i < bsf.u_knot_temp.size(); i++)
	{
		ioszz_1 << bsf.u_knot_temp[i] << endl;
	}
	ioszz_1.close();


	ofstream ioszz_11("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/spline_11.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_11 << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz_11.close();

	ofstream ioszz_12("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/spline_12.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_12 << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 0) << " " << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 1) << " " << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 2) << endl;
	}
	ioszz_12.close();
	ofstream ioszz_33("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/knot_2.txt");
	for (int i = 1; i < bsf.u_knot_temp.size() - 1; i++)
	{
		ioszz_33 << bsf.u_knot_temp[i] << endl;
	}
	ioszz_33.close();

	ofstream ioszz_34("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/weight.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_34 << 1 << endl;
	}
	ioszz_34.close();
	bsf.BSpline_surface_viewer_2(bspline_mat,5, 100, -2,patch_num);



	ResInfo res_temp;
	res_temp.spline = bspline_mat;
	res_temp.u_knot = bsf.u_knot_temp;
	if (zero_num >0)
	{
		res_temp.arv_e = min_avre;
		res_temp.max_e = min_maxe;
		res_temp.min_e = min_mine;
	}
	else
	{
		res_temp.arv_e = bsf.avr_energy;
		res_temp.min_e = bsf.min_energy;
		res_temp.max_e = bsf.max_energy;
	}
	res_temp.data_e = res_temp.arv_e * bsf.mesh2surface_temp.size();
	return res_temp;
}

BSplineNoCollision::ResInfo BSplineNoCollision::bspline_fitting_no_collision1(const std::string& input_mesh_file, Mesh& input_mesh, int patch_num,int seg_num, SegSurface::Path& one_path)
{
	int zero_num = 0;
	int init_num = 0;
	vector<double> res_avre(5, 100.0);
	vector<double> res_mine(5, 100.0);
	vector<double> res_maxe(5, 100.0);
	vector<Data> res_spline;
	res_spline.resize(5);
	vector<vector<double>> res_knot;
	res_knot.resize(5);
	zero_num = 0;
	iscollision.clear();
	iscollision.resize(5);
	auto ruling_direction = one_path.ruling_direction;
	auto center_point = one_path.center_point;
	auto translation_direction = one_path.translation_direction;

flag2:Data bspline_mat;
	bspline_mat.resize(ctr_num * 2, 3);
	Eigen::Vector3d center_trans = (one_path.center_point + one_path.translation_length * one_path.translation_direction);
	Eigen::Vector3d v0 = center_trans + one_path.ruling_direction + one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v1 = center_trans + one_path.ruling_direction - one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v2 = center_trans - one_path.ruling_direction + one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;
	Eigen::Vector3d v3 = center_trans - one_path.ruling_direction - one_path.ruling_direction.cross(one_path.translation_direction) * 0.7 * one_path.path_length;

	bspline_mat.row(0) << v0(0), v0(1), v0(2);
	bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
	bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
	bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
	for (int i = 1; i < ctr_num - 1; i++)
	{
		bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
		bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
	}

flag:VectorXd bspline;
	bspline.resize(bspline_mat.rows() * 3);
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		bspline[3 * i] = bspline_mat.coeffRef(i, 0);
		bspline[3 * i + 1] = bspline_mat.coeffRef(i, 1);
		bspline[3 * i + 2] = bspline_mat.coeffRef(i, 2);
	}

	// 2.normalize axis and aabb_axis
	int dim = axis.size();
	for (int k = 0; k < dim; k++)
	{
		axis[k].normalize();
	}

	dim = aabb_axis.size();
	for (int k = 0; k < dim; k++)
	{
		aabb_axis[k].normalize();
	}

	//3.设置节点向量
	BSplineFitting bsf(input_mesh);
	bsf.interested_area.resize(input_mesh.n_faces());
	for (int i = 0; i < input_mesh.n_faces(); i++)
	{
		bsf.interested_area[i] = 0;
	}

	for (int j = 0; j < one_path.cover_face.size(); j++)
	{
		bsf.interested_area[one_path.cover_face[j]] = 1;
	}

	bsf.CCD_initialization(bspline_mat);

	//4.初始化 细分的控制点及节点向量
	int max_sub_num = 13;
	divide_ctrpoint_.resize(max_sub_num);
	divide_u_knot_.resize(max_sub_num);
	divide_basis_.resize(max_sub_num);
	for (int i = 0; i < max_sub_num; i++)
	{
		divide_ctrpoint_[i].resize(pow(2, i));
		divide_u_knot_[i].resize(pow(2, i));
		divide_basis_[i].resize(pow(2, i));
	}

	divide_ctrpoint_[0][0] = bspline;
	divide_u_knot_[0][0] = bsf.u_knot_temp;
	divide_basis_[0][0].resize(bspline.size() / 6, bspline.size() / 6);
	divide_basis_[0][0].setIdentity();

	subdivide_tree.resize(1);
	std::pair<double, double> range(1, 1);
	Eigen::MatrixXd basis = divide_basis_[0][0];
	Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
	combined.setZero();
	combined.block(0, 0, basis.rows(), basis.cols()) = basis;
	combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
	subdivide_tree[0] = std::make_tuple(0, range, combined);

	stop_condition = 1;
	data_energy0 = 1;
	epsilon2 = 1e-10;
	max_step = 1.0;
	offset = 1e-10;
	margin = 1e-10;

	data_weight = 1;
	smooth_weight = 1;
	fast_barrier_weight = 1;

	subdivision_number = 0;
	barrier_lambda = 0.8;

	Eigen::MatrixXd V;// vn*3
	Eigen::MatrixXi F;// fn*3
	igl::read_triangle_mesh(input_mesh_file, V, F);
	std::cout << "before bvh init\n";
	BVH bvh;
	bvh.InitObstacle(V, F);

	int opt_num = 0;
	bsf.count_divide.resize(bspline_mat.size() / 6 - 4 + 1); //统计每个初始参数区间段的细分次数
	for (int i = 0; i < bsf.count_divide.size(); i++)
	{
		bsf.count_divide[i] = 0;
	}
	bsf.BSpline_surface_viewer_2(bspline_mat, 5, 100, -50-seg_num, patch_num);

	//储存初始B样条曲面
	ofstream ioszz("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_"+to_string(seg_num)+"spline_0.txt");
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		ioszz << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz.close();
	ofstream ioszz1("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "knot_0.txt");
	for (int i = 0; i < bsf.u_knot_temp.size(); i++)
	{
		ioszz1 << bsf.u_knot_temp[i] << endl;
	}
	ioszz1.close();

	iter_ = 0;
	iscollision[zero_num] = 0;
	is_good_init = 1;

	while (iscollision[zero_num] == 0 && is_good_init && bsf.iter_num < 10 && bsf.avr_energy>2e-4 && (bsf.iter_num < 6 || bsf.avr_energy < 0.05))
	{

		cout << "*** 迭代次数： " << opt_num << " , MAX_step: " << max_step << " ***" << endl;
		clock_t opt_time1 = clock();
		Optimization3D_time::gcl_optimization(bspline_mat, V, F, bvh, bsf);
		bsf.BSpline_surface_viewer_2(bspline_mat, 5, 100, opt_num, patch_num);
		clock_t opt_time2 = clock();
		opt_num++;
		iter_++;
		bsf.iter_num++;
	}

	if (is_good_init != 1 && init_num < 4 && zero_num == 0)
	{
		init_num++;
		bspline_mat.resize(ctr_num * 2, 3);
		Matrix3d Rota;
		Rota.setZero();
		ruling_direction = one_path.ruling_direction;
		//center_point = one_path.center_point;
		translation_direction = one_path.translation_direction;
		double c;
		double s;
		if (init_num == 1)//绕rullings顺时针旋转
		{
			c = cos(-20.0 / 180.0 * 3.1415926);
			s = sin(-20.0 / 180.0 * 3.1415926);
		}
		if (init_num == 2)
		{
			c = cos(20.0 / 180.0 * 3.1415926);
			s = sin(20.0 / 180.0 * 3.1415926);
		}
		if (init_num == 3)
		{
			c = cos(-40.0 / 180.0 * 3.1415926);
			s = sin(-40.0 / 180.0 * 3.1415926);
		}
		if (init_num == 4)
		{
			c = cos(40.0 / 180.0 * 3.1415926);
			s = sin(40.0 / 180.0 * 3.1415926);
		}
		Rota(0, 0) = c + (1 - c) * pow(ruling_direction[0], 2);
		Rota(0, 1) = (1 - c) * ruling_direction[0] * ruling_direction[1] - s * ruling_direction[2];
		Rota(0, 2) = (1 - c) * ruling_direction[0] * ruling_direction[2] + s * ruling_direction[1];
		Rota(1, 0) = (1 - c) * ruling_direction[0] * ruling_direction[1] + s * ruling_direction[2];
		Rota(1, 1) = c + (1 - c) * pow(ruling_direction[1], 2);
		Rota(1, 2) = (1 - c) * ruling_direction[1] * ruling_direction[2] - s * ruling_direction[0];
		Rota(2, 0) = (1 - c) * ruling_direction[0] * ruling_direction[2] - s * ruling_direction[1];
		Rota(2, 1) = (1 - c) * ruling_direction[1] * ruling_direction[2] + s * ruling_direction[0];
		Rota(2, 2) = c + (1 - c) * pow(ruling_direction[2], 2);
		center_point = Rota * center_point;
		translation_direction = Rota * translation_direction;
		cout << "init_translation_direction=" << one_path.translation_direction.norm() << endl;
		cout << "初始导线长度：" << ruling_direction.cross(one_path.translation_direction).norm() << endl;
		cout << "translation_direction=" << translation_direction.norm() << endl;
		cout << "导线长度：" << ruling_direction.cross(translation_direction).norm() << endl;
		Eigen::Vector3d center_trans = (center_point + one_path.translation_length * translation_direction);
		Eigen::Vector3d v0 = center_trans + ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
		Eigen::Vector3d v1 = center_trans + ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
		Eigen::Vector3d v2 = center_trans - ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
		Eigen::Vector3d v3 = center_trans - ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;

		bspline_mat.row(0) << v0(0), v0(1), v0(2);
		bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
		bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
		bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
		for (int i = 1; i < ctr_num - 1; i++)
		{
			bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
			bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
		}
		goto flag;
	}
	if (is_good_init && init_num > 0 && zero_num == 0)
	{
		one_path.ruling_direction = ruling_direction;
		one_path.center_point = center_point;
		one_path.translation_direction = translation_direction;
		init_num = 0;
	}

	if (iscollision[0] == 1)
	{
		iscollision[0] = 1;
		one_path.translation_length += 0.2;
		goto flag2;
	}

	if (bsf.avr_energy > 1.5 * 1e-3 && zero_num <= -1)
	{
		res_avre[zero_num] = bsf.avr_energy;
		res_maxe[zero_num] = bsf.max_energy;
		res_mine[zero_num] = bsf.min_energy;
		res_knot[zero_num] = bsf.u_knot_temp;
		res_spline[zero_num] = bspline_mat;
		zero_num++;
		if (zero_num <= 4)
		{
			bspline_mat.resize(ctr_num * 2, 3);
			Matrix3d Rota;
			Rota.setZero();
			ruling_direction = one_path.ruling_direction;
			//center_point = one_path.center_point;
			translation_direction = one_path.translation_direction;
			double c;
			double s;
			if (zero_num == 1)//绕rullings顺时针旋转
			{
				c = cos(-22.5 / 180.0 * 3.1415926);
				s = sin(-22.5 / 180.0 * 3.1415926);
			}
			if (zero_num == 2)
			{
				c = cos(22.5 / 180.0 * 3.1415926);
				s = sin(22.5 / 180.0 * 3.1415926);
			}
			if (zero_num == 3)
			{
				c = cos(-45.0 / 180.0 * 3.1415926);
				s = sin(-45.0 / 180.0 * 3.1415926);
			}
			if (zero_num == 4)
			{
				c = cos(45.0 / 180.0 * 3.1415926);
				s = sin(45.0 / 180.0 * 3.1415926);
			}
			Rota(0, 0) = c + (1 - c) * pow(ruling_direction[0], 2);
			Rota(0, 1) = (1 - c) * ruling_direction[0] * ruling_direction[1] - s * ruling_direction[2];
			Rota(0, 2) = (1 - c) * ruling_direction[0] * ruling_direction[2] + s * ruling_direction[1];
			Rota(1, 0) = (1 - c) * ruling_direction[0] * ruling_direction[1] + s * ruling_direction[2];
			Rota(1, 1) = c + (1 - c) * pow(ruling_direction[1], 2);
			Rota(1, 2) = (1 - c) * ruling_direction[1] * ruling_direction[2] - s * ruling_direction[0];
			Rota(2, 0) = (1 - c) * ruling_direction[0] * ruling_direction[2] - s * ruling_direction[1];
			Rota(2, 1) = (1 - c) * ruling_direction[1] * ruling_direction[2] + s * ruling_direction[0];
			Rota(2, 2) = c + (1 - c) * pow(ruling_direction[2], 2);
			center_point = Rota * center_point;
			translation_direction = Rota * translation_direction;
			translation_direction.normalize();
			Eigen::Vector3d center_trans = (center_point + one_path.translation_length * translation_direction);
			Eigen::Vector3d v0 = center_trans + ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v1 = center_trans + ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v2 = center_trans - ruling_direction + ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;
			Eigen::Vector3d v3 = center_trans - ruling_direction - ruling_direction.cross(translation_direction) * 0.7 * one_path.path_length;

			bspline_mat.row(0) << v0(0), v0(1), v0(2);
			bspline_mat.row(ctr_num * 2 - 1) << v3(0), v3(1), v3(2);
			bspline_mat.row(ctr_num - 1) << v1(0), v1(1), v1(2);
			bspline_mat.row(ctr_num) << v2(0), v2(1), v2(2);
			for (int i = 1; i < ctr_num - 1; i++)
			{
				bspline_mat.row(i) = (bspline_mat.row(0) * (ctr_num - 1 - i) + bspline_mat.row(ctr_num - 1) * (i)) / (ctr_num - 1);
				bspline_mat.row(i + ctr_num) = (bspline_mat.row(ctr_num) * (ctr_num - 1 - i) + bspline_mat.row(2 * ctr_num - 1) * (i)) / (ctr_num - 1);
			}
			goto flag;
		}
	}
	double min_avre = 1e6;
	double min_mine = 1e6;
	double min_maxe = 1e6;
	if (zero_num > 0)
	{

		for (int i = 0; i < zero_num; i++)
		{
			if (res_avre[i] < min_avre && iscollision[i] == 0)
			{
				min_avre = res_avre[i];
				min_mine = res_mine[i];
				min_maxe = res_maxe[i];
				bspline_mat = res_spline[i];
				bsf.u_knot_temp = res_knot[i];
			}
		}

	}



	//储存最终B样条曲面
	ofstream ioszz_("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "spline_1.txt");
	for (int i = 0; i < bspline_mat.rows(); i++)
	{
		ioszz_ << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz_.close();
	ofstream ioszz_1("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "knot_1.txt");
	for (int i = 0; i < bsf.u_knot_temp.size(); i++)
	{
		ioszz_1 << bsf.u_knot_temp[i] << endl;
	}
	ioszz_1.close();

	ofstream ioszz_11("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "spline_11.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_11 << bspline_mat.coeffRef(i, 0) << " " << bspline_mat.coeffRef(i, 1) << " " << bspline_mat.coeffRef(i, 2) << endl;
	}
	ioszz_11.close();

	ofstream ioszz_12("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "spline_12.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_12 << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 0) << " " << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 1) << " " << bspline_mat.coeffRef(i + bspline_mat.rows() / 2, 2) << endl;
	}
	ioszz_12.close();
	ofstream ioszz_33("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "knot_2.txt");
	for (int i = 1; i < bsf.u_knot_temp.size() - 1; i++)
	{
		ioszz_33 << bsf.u_knot_temp[i] << endl;
	}
	ioszz_33.close();

	ofstream ioszz_34("D:/programfiles/3/11/build/result/" + result_path_seg + "/123/" + to_string(patch_num) + "/seg_" + to_string(seg_num) + "weight.txt");
	for (int i = 0; i < bspline_mat.rows() / 2; i++)
	{
		ioszz_34 << 1 << endl;
	}
	ioszz_34.close();
	bsf.BSpline_surface_viewer_2(bspline_mat, 5, 100, -2-1-seg_num, patch_num);



	ResInfo res_temp;
	res_temp.spline = bspline_mat;
	res_temp.u_knot = bsf.u_knot_temp;
	if (zero_num > 0)
	{
		res_temp.arv_e = min_avre;
		res_temp.max_e = min_maxe;
		res_temp.min_e = min_mine;
	}
	else
	{
		res_temp.arv_e = bsf.avr_energy;
		res_temp.min_e = bsf.min_energy;
		res_temp.max_e = bsf.max_energy;
	}
	res_temp.data_e = res_temp.arv_e * bsf.mesh2surface_temp.size();
	return res_temp;
}