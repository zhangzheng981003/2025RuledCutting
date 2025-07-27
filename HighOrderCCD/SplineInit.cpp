#include "SplineInit.h"

SplineInit::SplineInit(Mesh& mesh) :mesh(mesh)
{
	average_edge_length = 0.0;
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		average_edge_length += mesh.calc_edge_length(mesh.edge_handle(i));
	}
	average_edge_length = average_edge_length / mesh.n_edges();

	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}
	rulings_length = (ptMax - ptMin).norm();
}

void SplineInit::initialization(SegSurface& ss)
{
	average_edge_length = 0.0;
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		average_edge_length += mesh.calc_edge_length(mesh.edge_handle(i));
	}
	average_edge_length = average_edge_length / mesh.n_edges();

	ss.collect_init_informarion();

	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}

	rulings_length = (ptMax - ptMin).norm();

}


double SplineInit::spline_basis(int spline_order, int cur_order, double arc_param)
{
	return combination(spline_order, cur_order) * pow(arc_param, cur_order) * pow(1 - arc_param, spline_order - cur_order);
}

int SplineInit::combination(int n, int m)
{
	std::vector<std::vector<int>> dp(n + 1, std::vector<int>(m + 1, 0));

	for (int i = 0; i <= n; ++i) {
		for (int j = 0; j <= std::min(i, m); ++j) {
			if (j == 0 || j == i) {
				dp[i][j] = 1;
			}
			else {
				dp[i][j] = dp[i - 1][j - 1] + dp[i - 1][j];
			}
		}
	}

	return dp[n][m];
}
int SplineInit::u_k_position(double u, vector<double>& u_knot)
{
	for (int i = 0; i < u_knot.size(); i++)
	{
		if (u == u_knot[0])
		{
			return 0;
		}
		else if (u < u_knot[i])
		{
			return i - 1;
		}
		else if (u == u_knot.back())
		{
			return u_knot.size() - 1;
		}
	}
}


void SplineInit::calc_basis_fuction(double u, int k, int poly_degree, std::vector<double>& basis_func)
{
	/*!
	*\brief b样条曲线基函数计算
	*\ param double u 参数值
	*\ param int k 区间 k。可以通过std::upper_bound()-1获得
	*\ param std::vector<double> & basis_func basis_func基函数值，对应N_{k-p},...,N_k
	*\ Returns:   void
	*/

	const int& p = poly_degree;
	const std::vector<double>& knots = u_knot_temp;
	basis_func.clear();
	basis_func.resize(p + 1);

	//2p+1个 N_{i,0}
	int n = 2 * p + 1;
	vector<double> temp(n, 0);
	temp[p] = 1;

	//迭代p次
	for (int j = 1; j <= p; ++j)
	{
		//区间 [k-p,k+p+1)
		for (int i = k - p, h = 0; h < (n - j); ++h, ++i)
		{
			//递推公式
			double a = (u - knots[i]);
			double dev = (knots[i + j] - knots[i]);
			a = (dev != 0) ? a / dev : 0;

			double b = (knots[i + j + 1] - u);
			dev = (knots[i + j + 1] - knots[i + 1]);
			b = (dev != 0) ? b / dev : 0;

			temp[h] = a * temp[h] + b * temp[h + 1];
		}
	}

	//拷贝前 p+1个值到basis_func
	std::copy(temp.begin(), temp.begin() + p + 1, basis_func.begin());
}

void SplineInit::calc_control_point(int spline_order,int pc_num)
{
	std::vector<T> trip;
	for (int i = 0; i < pc_num; i++)
	{
		for (int j = 0; j < (spline_order + 1); j++)
		{
			trip.emplace_back(3 * i, 3 * j, combination(spline_order, j) * pow(loop[i].arc_length_param, j) * pow(1 - loop[i].arc_length_param, spline_order - j));
			trip.emplace_back(3 * i + 1, 3 * j + 1, combination(spline_order, j) * pow(loop[i].arc_length_param, j) * pow(1 - loop[i].arc_length_param, spline_order - j));
			trip.emplace_back(3 * i + 2, 3 * j + 2, combination(spline_order, j) * pow(loop[i].arc_length_param, j) * pow(1 - loop[i].arc_length_param, spline_order - j));
		}
	}
	
	m_A.resize(3 * pc_num, 3 * (spline_order + 1));
	m_A.setFromTriplets(trip.begin(), trip.end());

	solver.compute(m_A.transpose() * m_A);
	m_x = solver.solve(m_A.transpose() * m_b);

	std::cout << "能量值： " << (m_A * m_x - m_b).norm() << endl;
	std::cout << "残差： " << (m_A.transpose() * m_A * m_x - m_A.transpose() * m_b).norm() << endl;
	sample_points();
}


std::vector<double> SplineInit::generate_equally_spaced_vector(int num_elements)//yes
{
	if (num_elements <= 1) {
		std::cerr << "Invalid input. Number of elements should be a positive integer." << std::endl;
		return std::vector<double>();
	}

	std::vector<double> result;
	double step = 1.0 / (num_elements - 1); // 计算等分的步长

	for (int i = 0; i < num_elements; ++i) {
		double value = i * step; // 计算等分的值
		result.push_back(value);
	}

	return result;
}



void SplineInit::sample_points()
{
	control_points.resize(spline_order + 1, 3);
	for (int i = 0; i < control_points.rows(); i++)
	{
		control_points(i, 0) = m_b[3 * i];
		control_points(i, 1) = m_b[3 * i + 1];
		control_points(i, 2) = m_b[3 * i + 2];
	}

	vector<double> res = generate_equally_spaced_vector(1000);
	ofstream iosss("res0.txt");
	for (int i = 0; i < res.size(); i++)
	{
		VectorXd p0; p0.setZero(3);
		for (int j = 0; j < control_points.rows(); j++)
		{
			p0[0] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 0);
			p0[1] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 1);
			p0[2] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 2);
		}
		iosss << p0[0] << " " << p0[1] << " " << p0[2] << endl;
		///cout << res[i] << endl;
	}
	iosss.close();

	/*
	ofstream iosss1("bezier.txt");
	for (int i = 0; i < m_x.size(); i++)
	{
		iosss1 << m_x[i] << endl;
		//cout << m_x[i] << endl;
	}
	iosss1.close();*/
	
	//cout << control_points << endl;
}

void SplineInit::sample_init_spline(SegSurface::Path& path)
{
	loop.clear();
	for (int i = 0; i < path.path.size(); i++)
	{
		PathInfo pi;
		pi.rank = i;
		pi.input_point[0] = path.path[i][0];
		pi.input_point[1] = path.path[i][1]; 
		pi.input_point[2] = path.path[i][2]; 
		pi.e_input_point = path.path[i];
		loop.push_back(pi);
	}

	int gap = std::round(loop.size() / init_ctr_num);
	spline_init.resize(init_ctr_num * 3);
	spline_init.setZero();
	for (int i = 0; i < init_ctr_num; i++)
	{
		if (i != init_ctr_num-1)
		{
			spline_init[3 * i + 0] = loop[i * gap].e_input_point[0];
			spline_init[3 * i + 1] = loop[i * gap].e_input_point[1];
			spline_init[3 * i + 2] = loop[i * gap].e_input_point[2];
		}
		else
		{
			spline_init[3 * i + 0] = loop[loop.size() - 1].e_input_point[0];
			spline_init[3 * i + 1] = loop[loop.size() - 1].e_input_point[1];
			spline_init[3 * i + 2] = loop[loop.size() - 1].e_input_point[2];
		}
	}
	cout << spline_init[0] << " " << spline_init[1] << " " << spline_init[2] << " " << path.path[0][0] << " " << path.path[0][1] << " " << path.path[0][2] << endl;
	cout << spline_init[init_ctr_num * 3-3] << " " << spline_init[init_ctr_num * 3-2] << " " << spline_init[init_ctr_num * 3-1] << " " << path.path[path.path.size()-1][0] << " " << path.path[path.path.size() - 1][1] << " " << path.path[path.path.size() - 1][2] << endl;
}

void SplineInit::calc_projection_distance2_gradient(VectorXd& spline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,Eigen::VectorXd& decrease_direction)
{
	data_value = 0;
	grad.resize(spline.size());
	grad.setZero();
	hessian.resize(spline.size(), spline.size());
	hessian.setZero();
	


	Handle(Geom_BSplineCurve) curve = nullptr;
	int ctr_num = spline.size()/3;
	cout << ctr_num << endl;
	TColgp_Array1OfPnt Poles(1, ctr_num);
	TColStd_Array1OfReal PolesWeight(1, ctr_num);
	for (int i = 0; i < ctr_num; i++)
	{
		Poles.SetValue(i + 1, gp_Pnt(spline[3 * i], spline[3 * i + 1], spline[3 * i + 2]));
		PolesWeight.SetValue(i + 1, 1);
	}
	Standard_Integer PNum = spline.size()/3;
	Standard_Integer KNum = PNum - 3 + 1;
	TColStd_Array1OfReal knots(1, KNum);
	TColStd_Array1OfInteger mults(1, KNum);


	for (int i = 0; i < KNum; ++i)
	{
		//cout << i <<" " << u_knot_temp[i + poly_degree] << "success" << endl;
		knots.SetValue(i + 1, u_knot_temp[i + 3]);

		if (i == 0 || i == KNum - 1)
		{
			mults.SetValue(i + 1, 3 + 1);
		}
		else
		{
			mults.SetValue(i + 1, 1);
		}
	}
	curve = new Geom_BSplineCurve(Poles, PolesWeight, knots, mults, 3);




	/*
	TColgp_Array1OfPnt controlPoints(1, spline_order + 1);
	for (int i = 0; i < spline_order + 1; i++)
	{
		controlPoints.SetValue(i + 1, gp_Pnt(spline[3 * i], spline[3 * i + 1], spline[3 * i + 2]));
	}
	//Standard_Real tola = 1e-6;
	Handle(Geom_Curve) bezierCurve = new Geom_BezierCurve(controlPoints);
	*/
	for (int k = 0; k < loop.size(); k++)
	{
		// 定义要计算投影距离的点
		gp_Pnt pointToProject(loop[k].input_point[0].value().value(), loop[k].input_point[1].value().value(), loop[k].input_point[2].value().value());
		// 计算点到贝塞尔曲线的投影
		GeomAPI_ProjectPointOnCurve projectPointOnCurve(pointToProject, curve);

		if (projectPointOnCurve.NbPoints() > 0)
		{
			// 获取投影点
			Standard_Real u = projectPointOnCurve.LowerDistanceParameter();
			//cout << projectPointOnCurve.LowerDistance() << endl;
			loop[k].u = u;
			//cout <<"u="<< u << endl;
			Vec12 X;
			X.resize(3 * (ctr_num));
			for (int r = 0; r < ctr_num ; r++)
			{
				X(3 * r).value() = spline(3 * r);
				X(3 * r + 1).value() = spline(3 * r + 1);
				X(3 * r + 2).value() = spline(3 * r + 2);
			}

			// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
			for (int id = 0; id < 3 * ctr_num ; id++)
			{
				X(id).derivatives().resize(3 * (ctr_num));
				X(id).derivatives().setZero();
				X(id).derivatives()(id) = 1;
				X(id).value().derivatives() = inner_derivative_t::Unit(3 * (ctr_num), id);
			}

			for (int idx = 0; idx < 3 * (ctr_num); idx++) {
				for (int id = 0; id < 3 * (ctr_num); id++)
				{
					X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (ctr_num));
				}
			}


			std::vector<Vec3> P;
			P.resize((ctr_num));

			for (int i = 0; i < P.size(); i++)
			{
				P[i][0] = X[3 * i];
				P[i][1] = X[3 * i + 1];
				P[i][2] = X[3 * i + 2];
			}

			int u_interval_order = u_k_position(loop[k].u, u_knot_temp);
			vector<double> u_basis_temp;
			if (u_interval_order == 0)
			{
				vector<double> basis_temp(3+ 1, 0);
				basis_temp[0] = 1;
				u_basis_temp = basis_temp;
			}
			else if (u_interval_order == u_knot_temp.size() - 1)
			{
				vector<double> basis_temp(3 + 1, 0);
				basis_temp[3] = 1;
				u_basis_temp = basis_temp;
			}
			else
			{
				calc_basis_fuction(loop[k].u, u_interval_order, 3, u_basis_temp);
			}

			loop[k].u_basis_temp = u_basis_temp;
			loop[k].u_interval_order = u_interval_order;

			
			Vec3 bspline_surface_value;
			bspline_surface_value.setZero();
			if (u_interval_order == 0)
			{
				Vec3 pi0;
				pi0 = P[0];
				bspline_surface_value[0] = pi0[0];
				bspline_surface_value[1] = pi0[1];
				bspline_surface_value[2] = pi0[2];
			}
			else if (u_interval_order == u_knot_temp.size() - 1)
			{
				Vec3 pi0;
				pi0 = P[ctr_num - 1];
				bspline_surface_value[0] = pi0[0];
				bspline_surface_value[1] = pi0[1];
				bspline_surface_value[2] = pi0[2];
			}
			else
			{
				for (int i = 0; i < u_basis_temp.size(); i++)
				{
					Vec3 pi0;
					pi0 = P[u_interval_order - 3 + i];
					bspline_surface_value[0] +=  u_basis_temp[i] * pi0[0];
					bspline_surface_value[1] +=  u_basis_temp[i] * pi0[1];
					bspline_surface_value[2] +=  u_basis_temp[i] * pi0[2];
				}
			}
			
			scalar_t e = (bspline_surface_value - loop[k].input_point) * (bspline_surface_value - loop[k].input_point).transpose();
			data_value += e;
			grad += e.value().derivatives();
			
			Eigen::MatrixXd B;
			B.resize(3 * ctr_num, 3 * ctr_num);
			for (int r = 0; r < 3 * ctr_num; r++)
			{
				B.row(r) = e.derivatives()(r).derivatives().transpose();
				//cout << "r="<<r<<" " << e.derivatives()(r).derivatives().transpose().norm() << endl;
			}
			//cout << "u1=" << u << endl;
			hessian += B;
			
		}
		
		else 
		{
			continue;
			std::cout << "无法计算投影点" << std::endl;
		}
	}
	
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
	Eigen::MatrixXd U = svd.matrixU();
	Eigen::MatrixXd V = svd.matrixV();
	Eigen::VectorXd singularValues = svd.singularValues();
	Eigen::VectorXd inverse_singularValues;
	inverse_singularValues.resize(singularValues.size());
	inverse_singularValues.setZero();
	for (int i = 0; i < singularValues.size(); i++)//4.62223e-33数值不稳定
	{
		if (singularValues[i] < 1e-6)
		{
			singularValues[i] = 1e-6;
			inverse_singularValues[i] = 1.0 / singularValues[i];
		}
		else
		{
			inverse_singularValues[i] = 1.0 / singularValues[i];
		}
	}

	auto pseudo_inverse_matrix = V * inverse_singularValues.asDiagonal() * U.transpose();
	decrease_direction = -pseudo_inverse_matrix * grad;
	
	//decrease_direction.normalize();
}


double SplineInit::calc_energy(VectorXd& spline)
{
	double data_value = 0.0;
	int ctr_num = spline.size() / 3;
	vector<double> e_intervals(u_knot_temp.size(), 0.0);
	
	for (int k = 0; k < loop.size(); k++)
	{
		if (loop[k].u_basis_temp.size() > 0)
		{
			std::vector<Vector3d> P;
			P.resize(ctr_num);
			for (int i = 0; i < P.size(); i++)
			{
				P[i][0] = spline[3 * i];
				P[i][1] = spline[3 * i + 1];
				P[i][2] = spline[3 * i + 2];
			}

			Vector3d bspline_surface_value;
			bspline_surface_value.setZero();
			if (loop[k].u_interval_order == 0)
			{
				Vector3d pi0;
				pi0 = P[0];
				bspline_surface_value[0] = pi0[0];
				bspline_surface_value[1] = pi0[1];
				bspline_surface_value[2] = pi0[2];
			}
			else if (loop[k].u_interval_order == u_knot_temp.size() - 1)
			{
				Vector3d pi0;
				pi0 = P[ctr_num - 1];
				bspline_surface_value[0] = pi0[0];
				bspline_surface_value[1] = pi0[1];
				bspline_surface_value[2] = pi0[2];
			}
			else
			{
				for (int i = 0; i < loop[k].u_basis_temp.size(); i++)
				{
					Vector3d pi0;
					pi0 = P[loop[k].u_interval_order - 3 + i];
					bspline_surface_value[0] += loop[k].u_basis_temp[i] * pi0[0];
					bspline_surface_value[1] += loop[k].u_basis_temp[i] * pi0[1];
					bspline_surface_value[2] += loop[k].u_basis_temp[i] * pi0[2];
				}
			}

			data_value += (bspline_surface_value - loop[k].e_input_point).transpose() * (bspline_surface_value - loop[k].e_input_point);

			e_intervals[loop[k].u_interval_order] += (bspline_surface_value - loop[k].e_input_point).transpose() * (bspline_surface_value - loop[k].e_input_point);
		}
	}
	
	double max_temp_e = 0;
	for (int i = 0; i < e_intervals.size() - 1; i++)
	{
		int init_interval = floor(u_knot_temp[i]);
		if (max_temp_e < e_intervals[i] )
		{

			max_temp_e = e_intervals[i];
			if (u_knot_temp[i + 1] == 0)
			{
				int i_ = i;
				while (u_knot_temp[i_ + 1] == 0)
				{
					i_++;
				}
				insert1_idx = i_;
			}
			else if (u_knot_temp[i] == u_knot_temp.back())
			{
				int i_ = i;
				while (u_knot_temp[i_] == u_knot_temp.back())
				{
					i_--;
				}
				insert1_idx = i_;
			}
			else
			{
				insert1_idx = i;
			}

		}
	}

	return data_value;
}


double SplineInit::line_search(VectorXd& spline,double data_value, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction)
{
	
	double step_size = 1.0;// Initial step size
	int iter_num = 0;
	while (step_size>0.01)
	{
		Eigen::VectorXd new_spline = spline + step_size * decrease_direction;
		double new_data_value = calc_energy(new_spline);
		std::cout << endl;
		std::cout << "step: " << step_size << endl;
		std::cout << "data_value_new: " << new_data_value << endl;
		std::cout << "old_total: " << data_value << endl; //data_value.value() <<  endl;

		if (new_data_value > data_value+ c_parm * step_size)
		{
			step_size *= beta;
		}
		else
		{
			return step_size;
			break;
		}

		iter_num++;
		if (step_size <= 0.01)
		{
			step_size = 0;
			return step_size;
		}
	}
}

void SplineInit::energy_opt(VectorXd& spline, SegSurface::Path& path)
{
	
	sample_init_spline(path);
	int iter_num = 0;
	spline = spline_init;
	cout << loop.size() << endl;
	u_knot_temp.clear();
	int KnotNum = spline.size() / 3 + 3 + 1;
	u_knot_temp.resize(KnotNum);
	int no_repeat_interval = spline.size() / 3 - 3;
	for (int i = 0; i < KnotNum; i++)
	{
		if (i < 3 + 1)
		{
			u_knot_temp[i] = 0;
		}
		else if (i > KnotNum - 3 - 1)
		{
			u_knot_temp[i] = no_repeat_interval;
		}
		else
		{
			u_knot_temp[i] = i - 3;
		}
	}
	while (iter_num==0||(data_vvv>loop.size()*1e-5&&spline.size()<40))
	{
		std::vector<Vector3d> P;
		P.resize(spline.size() / 3);
		for (int i = 0; i < P.size(); i++)
		{
			P[i][0] = spline[3 * i];
			P[i][1] = spline[3 * i + 1];
			P[i][2] = spline[3 * i + 2];
		}
		vector<double> res = uni_para(20);
		vector<Vector3d> curve_res;

		for (int i = 0; i < res.size(); i++)
		{
			int k = u_k_position(res[i], u_knot_temp);
			vector<double> u_basis_temp;
			calc_basis_fuction(res[i], k, poly_degree, u_basis_temp);
			Vector3d bspline_surface_value;
			bspline_surface_value.setZero();
			if (k == 0)
			{
				Vector3d pi0;
				pi0 = P[0];
				bspline_surface_value = pi0;
			}
			else if (k == u_knot_temp.size() - 1)
			{
				Vector3d pi0;
				pi0 = P[spline.size() / 3 - 1];
				bspline_surface_value = pi0;
			}
			else
			{
				for (int i = 0; i < u_basis_temp.size(); i++)
				{
					Vector3d pi0;
					pi0 = P[k - 3 + i];
					bspline_surface_value += u_basis_temp[i] * pi0;
				}
			}
			curve_res.push_back(bspline_surface_value);
		}
		ofstream iosss("curve" + to_string(iter_num) + ".txt");
		for (int i = curve_res.size() - 2; i >= 0; i--)
		{
			iosss << curve_res[i + 1][0] << " " << curve_res[i + 1][1] << " " << curve_res[i + 1][2] << " " << curve_res[i][0] << " " << curve_res[i][1] << " " << curve_res[i][2] << endl;
		}
		iosss.close();



		Eigen::VectorXd grad, new_spline,derection;
		scalar_t data_value;
		Data hessien;
		
		calc_projection_distance2_gradient(spline, data_value, grad, hessien, derection);
		//grad.normalize();
		//cout << "梯度："<<grad << endl;
		derection[0] = 0;
		derection[1] = 0;
		derection[2] = 0;
		derection[derection.size() - 3] = 0;
		derection[derection.size() - 2] = 0;
		derection[derection.size() - 1] = 0;
		double data_term = calc_energy(spline);
		
		double step = line_search(spline, data_term, grad, derection);
		//Eigen::VectorXd direction;
		new_spline = spline + step * derection;
		
		
		spline = new_spline;
		iter_num++;
		data_vvv = calc_energy(spline);
		std::cout << "第 " << iter_num << " 次迭代：" << data_vvv << endl;
		double termination_condition = 1.0-data_vvv/data_term;
		std::cout << "termination_condition= "<<termination_condition << endl;
		if (termination_condition < 0.4 )
		{
			double u_loc = u_knot_temp[insert1_idx] * 0.5 + u_knot_temp[insert1_idx + 1] *0.5;
			add_ctr_point(spline, u_loc);
			loop.clear();
			for (int i = 0; i < path.path.size(); i++)
			{
				PathInfo pi;
				pi.rank = i;
				pi.input_point[0] = path.path[i][0];
				pi.input_point[1] = path.path[i][1];
				pi.input_point[2] = path.path[i][2];
				pi.e_input_point = path.path[i];
				loop.push_back(pi);
			}
		}
	}

	cout << "迭代结束！" << endl;
	control_points.resize(spline.size()/3 ,3);
	for (int i = 0; i < control_points.rows(); i++)
	{
		control_points(i, 0) = spline[3 * i];
		control_points(i, 1) = spline[3 * i + 1];
		control_points(i, 2) = spline[3 * i + 2];
	}

	/*
	vector<double> res = generate_equally_spaced_vector(1000);
	for (int i = 0; i < res.size(); i++)
	{
		VectorXd p0; p0.setZero(3);
		for (int j = 0; j < control_points.rows(); j++)
		{
			p0[0] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 0);
			p0[1] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 1);
			p0[2] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 2);
		}
	}*/
}

void SplineInit::draw_point()
{
	std::string pathname("spline.txt");
	ifstream infile;
	infile.open(pathname.data());
	assert(infile.is_open());
	vector<double> suanz;
	string s;
	vector<vector<double>> res0;
	while (getline(infile, s)) {
		istringstream is(s);
		double d;
		while (!is.eof()) {
			is >> d;
			suanz.push_back(d);
		}
		res0.push_back(suanz);
		suanz.clear();
		s.clear();
	}
	infile.close();

	//cout << res0.size() << endl;
	control_points.resize(spline_order + 1, 3);
	for (int i = 0; i < control_points.rows(); i++)
	{
		control_points(i, 0) = res0[3 * i][0];
		control_points(i, 1) = res0[3 * i + 1][0];
		control_points(i, 2) = res0[3 * i + 2][0];
	}

	vector<double> res = generate_equally_spaced_vector(1000);
	//ofstream iosss("res0.txt");
	for (int i = 0; i < res.size(); i++)
	{
		VectorXd p0; p0.setZero(3);
		for (int j = 0; j < control_points.rows(); j++)
		{
			p0[0] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 0);
			p0[1] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 1);
			p0[2] += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * control_points.coeff(j, 2);
		}
		//iosss << p0[0] << " " << p0[1] << " " << p0[2] << endl;
	}
	//iosss.close();
}


void SplineInit::draw_surface()
{
	//vector<double> u_vec = generate_equally_spaced_vector(100);
	vector<double> res = uni_para(20);
	ofstream iosss("ccc.txt");
	for (int i = 0; i < res.size(); i++)
	{
		for (int k = 0; k < res.size(); k++)
		{
			VectorXd p0; p0.setZero(3);
			for (int j = 0; j < spline_order+1; j++)
			{
				for (int j2 = 0; j2 < 2; j2++)
				{
					p0 += combination(spline_order, j) * pow(res[i], j) * pow(1 - res[i], spline_order - j) * combination(1, j2) * pow(res[k], j2) * pow(1 - res[k], 1 - j2) * cp_[j][j2];
				}
			
			}
			iosss << p0[0] << " " << p0[1] << " " << p0[2] << endl;
		}
	}
	iosss.close();
}

MatrixXd SplineInit::calc_up_degree_matrix(int spline_order)
{
	MatrixXd up_mat(spline_order + 2, spline_order + 1); up_mat.setZero();
	up_mat(0, 0) = 1; up_mat(spline_order + 1, spline_order) = 1;
	for (int row = 1; row < spline_order + 1; row++)
	{
		up_mat(row, row - 1) = double(row) / double((spline_order + 1));
		up_mat(row, row) = 1 - (double(row) / double((spline_order + 1)));
	}
	return up_mat;
}

void SplineInit::calc_pca_energy(VectorXd& curve, int sample_num, SegSurface& ss, SegSurface::Path& path)
{
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		vector<Point> v_p_3;
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}

		auto f_normal = mesh.calc_face_normal(*f_it);
		Vec3 f_n;
		f_n[0] = f_normal[0];
		f_n[1] = f_normal[1];
		f_n[2] = f_normal[2];
		cgal_triangles.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()), f_it->idx(), f_n);
		v_p_3.clear();
	}

	Tree tree;
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();

	vector<double> res = generate_equally_spaced_vector(sample_num);//yes
	sample_curve_temp.clear();
	int order = (curve.size() / 3) - 1;

	Vector3d p_mean; p_mean.setZero();
	for (int k = 0; k < res.size(); k++)
	{
		Vector3d curve_p; curve_p.setZero();
		for (int i = 0; i < order + 1; i++)
		{
			Vector3d pi = { curve[3 * i],curve[3 * i + 1] ,curve[3 * i + 2] };
			curve_p += combination(order, i) * (pow(1 - res[k], order - i) * pow(res[k], i)) * pi;
		}

		SampleCurve sc;
		sc.p = curve_p;
		sc.u = res[k];

		Vector3d tangent; tangent.setZero();
		for (int i = 0; i < order; i++)
		{
			Vector3d pi = { curve[3 * i],curve[3 * i + 1] ,curve[3 * i + 2] };
			Vector3d pi1 = { curve[3 * (i + 1)],curve[3 * (i + 1) + 1] ,curve[3 * (i + 1) + 2] };
			tangent += combination(order - 1, i) * (pow(1 - res[k], order - 1 - i) * pow(res[k], i)) * (pi1 - pi);
		}

		sc.tangent_vector = tangent; sc.tangent_vector.normalize();
		Point p_close(curve_p[0], curve_p[1], curve_p[2]);
		auto result = tree.closest_point_and_primitive(p_close);

		sc.normal_vector = { result.second->my_face_normal[0].value().value(), result.second->my_face_normal[1].value().value() , result.second->my_face_normal[2].value().value() };
		sc.normal_vector.normalize();
		sc.rulling_vector = sc.normal_vector.cross(sc.tangent_vector);
		sc.rulling_vector.normalize();
		sample_curve_temp.push_back(sc);
		p_mean += sc.p;
	}

	p_mean /= res.size(); std::cout << "中心点： " << p_mean << endl;
	Matrix3d mat1; mat1.setZero();
	Matrix3d mat2; mat2.setZero();
	Matrix3d mat; mat.setZero();
	for (int i = 0; i < sample_curve_temp.size(); i++)
	{
		mat1 += ((p_mean - sample_curve_temp[i].p).normalized()) * ((p_mean - sample_curve_temp[i].p).normalized()).transpose();
		mat2 += sample_curve_temp[i].rulling_vector * sample_curve_temp[i].rulling_vector.transpose();
	}

	double weight = 0.1;
	mat = weight * mat1 + mat2;
	EigenSolver<Matrix3d> es(mat);
	Matrix3d D = es.pseudoEigenvalueMatrix();
	Matrix3d V = es.pseudoEigenvectors();

	vector<double> evv = { D.coeffRef(0,0),D.coeffRef(1,1), D.coeffRef(2,2) };
	auto min_it = std::min_element(evv.begin(), evv.end());
	int min_index = 0;
	if (min_it != evv.end())
	{
		min_index = std::distance(evv.begin(), min_it);
		std::cout << "最小值为: " << *min_it << " 在索引位置: " << min_index << std::endl;
	}
	else 
	{
		std::cout << "向量为空" << std::endl;
	}
	
	std::cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
	std::cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
	std::cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
	std::cout << "最小特征值： " << *min_it << endl;
	std::cout << "最小特征向量： " << V.col(min_index).transpose() << endl;
	Vector3d plane_normal = V.col(min_index);// .transpose();
	double term1 = plane_normal.transpose() * mat1 * plane_normal;
	double term2 = plane_normal.transpose() * mat2 * plane_normal;
	//VectorXd plane_n; plane_.re
	std::cout << "第一项能量： " << weight * term1 << " 第二项能量： " << term2 << endl;

	//计算平移的距离
	Mesh::Point p0(p_mean.x(), p_mean.y(), p_mean.z());
	double pv = -p_mean.transpose() * plane_normal; Mesh::Point pn(plane_normal[0], plane_normal[1], plane_normal[2]);
	double p_z = (-pv - plane_normal[0] * (p_mean.x() + 1) - plane_normal[1] * (p_mean.y() + 1)) / plane_normal.z();
	Mesh::Point d0(1, 1, p_z - p_mean.z()); d0.normalize();
	Mesh::Point p1 = p0 + d0 * rulings_length * 10;
	Mesh::Point p2 = p0 - d0 * rulings_length * 10;
	Mesh::Point p3 = p0 + d0.cross(pn) * rulings_length * 10;
	Mesh::Point p4 = p0 - d0.cross(pn) * rulings_length * 10;

	curve_pro_point.clear();
	rulings_pro_point.clear();
	Plane plane(plane_normal.x(), plane_normal.y(), plane_normal.z(), pv);
	Vector_3 normal = plane.orthogonal_vector();
	for (int i = 0; i < sample_curve_temp.size(); i++)
	{
		Point point_(sample_curve_temp[i].p.x(), sample_curve_temp[i].p.y(), sample_curve_temp[i].p.z());
		Vector_3 pointToPlaneVec(point_ - plane.projection(point_));
		Point projectedPoint = point_ - pointToPlaneVec;
		curve_pro_point.push_back(projectedPoint);


		Vector3d point1 = sample_curve_temp[i].p + sample_curve_temp[i].rulling_vector * 0.5 * rulings_length;
		Vector3d point2 = sample_curve_temp[i].p - sample_curve_temp[i].rulling_vector * 0.5 * rulings_length;
		
		Point p1(point1.x(), point1.y(), point1.z());
		Point p2(point2.x(), point2.y(), point2.z());
		Vector_3 pointToPlaneVec1(p1 - plane.projection(p1));
		Point projectedPoint1 = p1 - pointToPlaneVec1;
		Vector_3 pointToPlaneVec2(p2 - plane.projection(p2));
		Point projectedPoint2 = p2 - pointToPlaneVec2;
		rulings_pro_point.push_back(make_pair(projectedPoint1, projectedPoint2));
	}

	Vector3d dir_init; dir_init.setZero();
	for (int i = 0; i < rulings_pro_point.size(); i++)
	{
		Vector3d p1 = { rulings_pro_point[i].first.x(),rulings_pro_point[i].first.y(), rulings_pro_point[i].first.z() };
		Vector3d p2 = { rulings_pro_point[i].second.x(),rulings_pro_point[i].second.y(), rulings_pro_point[i].second.z() };
		dir_init += (p1 - p2).normalized();	
	}
	dir_init.normalize();

	double min_translation_length = 2 * average_edge_length;
	calc_spline_init_info(p_mean, plane_normal, dir_init, ss,  path, min_translation_length);

}

void SplineInit::calc_point_projection(Point& point, Plane& plane)
{
	///Point point(1.0, 9.0, -3.0); // Replace with your point coordinates
	//Plane plane(0.0, 0.0, 1.0, -100.0); // Replace with plane equation coefficients (ax + by + cz + d = 0)

	// Calculate the vector from the point to the plane
	Vector_3 normal = plane.orthogonal_vector();
	Vector_3 pointToPlaneVec(point - plane.projection(point));

	// Calculate the projected point
	Point projectedPoint = point - pointToPlaneVec;

	// Output the coordinates of the projected point
	std::cout << "Coordinates of the projected point: " << projectedPoint << std::endl;
}

void SplineInit::calc_spline_init_info(Vector3d& center, Vector3d& direction, Vector3d& ruling_direction, SegSurface& ss, SegSurface::Path& path, double min_translation_length)
{
	path.center_point = center;
	
	Vector3d weight_face_normal; weight_face_normal.setZero();
	for (int i = 0; i < path.cover_face.size(); i++)
	{
		weight_face_normal += ss.face_basis_info_temp[path.cover_face[i]].face_normal.normalized();
	}
	weight_face_normal.normalize();

	if (weight_face_normal.dot(direction) > 0)
	{
		//
	}
	else
	{
		direction = -direction;
	}

	path.translation_direction = direction;

	double path_length = 0.0;
	Vector3d path_direction; path_direction.setZero();
	for (int i = 1; i < path.path.size(); i++)
	{
		
		Vector3d dir_ = path.path[i] - path.path[i - 1];
		path_length += dir_.norm();
		path_direction += dir_;
	}
	
	path_direction.normalize();
	std::cout << path_direction[0] << " " << path_direction[1] << " " << path_direction[2] << " " << endl;

	path.path_direction = path_direction;
	path.path_length = path_length;
	std::cout << "lent" << path_length << endl;
	path.ruling_direction = ruling_direction;

	int translation_num = 0;
	while (true)
	{
		Eigen::Vector3d center_trans = (center + translation_num * min_translation_length * direction);
		Eigen::Vector3d v0 = center_trans + ruling_direction * rulings_length + ruling_direction.cross(direction) * 0.7 * path_length;
		Eigen::Vector3d v1 = center_trans + ruling_direction * rulings_length - ruling_direction.cross(direction) * 0.7 * path_length;
		Eigen::Vector3d v2 = center_trans - ruling_direction * rulings_length + ruling_direction.cross(direction) * 0.7 * path_length;
		Eigen::Vector3d v3 = center_trans - ruling_direction * rulings_length - ruling_direction.cross(direction) * 0.7 * path_length;
		
		vector<Triangle> plane_triangle; plane_triangle.clear();
		plane_triangle.push_back(Triangle(Point(v0[0], v0[1], v0[2]), Point(v1[0], v1[1], v1[2]), Point(v3[0], v3[1], v3[2])));
		plane_triangle.push_back(Triangle(Point(v0[0], v0[1], v0[2]), Point(v2[0], v2[1], v2[2]), Point(v1[0], v1[1], v1[2])));

		int intersetion_num = 0;
		for (auto fit : plane_triangle)
		{
			if (ss.tree.do_intersect(fit))
			{
				intersetion_num++;
				translation_num++;
				break;
			}
		}

		Mesh::Point p1(v0.x(), v0.y(), v0.z());
		Mesh::Point p2(v1.x(), v1.y(), v1.z());
		Mesh::Point p3(v2.x(), v2.y(), v2.z());
		Mesh::Point p4(v3.x(), v3.y(), v3.z());

		Mesh mesh_res;
		auto mv0 = mesh_res.add_vertex(p1);
		auto mv1 = mesh_res.add_vertex(p2);
		auto mv2 = mesh_res.add_vertex(p3);
		auto mv3 = mesh_res.add_vertex(p4);
		mesh_res.add_face(mv0, mv1, mv3);
		mesh_res.add_face(mv0, mv3, mv2);

		MeshTools::WriteMesh(mesh_res, "./init_mesh_res.obj");

		if (intersetion_num == 0)
		{
			break;
		}
	}
	path.translation_length = translation_num * min_translation_length+0.1;
}


void SplineInit::calc_soomth_surface(VectorXd& curve, int sample_num)
{
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		vector<Point> v_p_3;
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}

		auto f_normal = mesh.calc_face_normal(*f_it);
		Vec3 f_n;
		f_n[0] = f_normal[0];
		f_n[1] = f_normal[1];
		f_n[2] = f_normal[2];
		cgal_triangles.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()), f_it->idx(), f_n);
		v_p_3.clear();
	}

	Tree tree;
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();

	vector<double> res = uni_para(sample_num);//yes
	sample_curve_temp.clear();
	int order = (curve.size() / 3) - 1;
	ofstream ioss("1122.txt");
	int ctr_num = curve.size() / 3;
	ofstream ioss1("tangent.txt");
	for (int k = 0; k < res.size(); k++)
	{
		//cout << res[k] << endl;
		int u_interval_order = u_k_position(res[k], u_knot_temp);
		vector<double> u_basis_temp;
		calc_basis_fuction(res[k], u_interval_order, poly_degree, u_basis_temp);
		std::vector<Vector3d> P;
		P.resize(ctr_num);
		for (int i = 0; i < P.size(); i++)
		{
			P[i][0] = curve[3 * i];
			P[i][1] = curve[3 * i + 1];
			P[i][2] = curve[3 * i + 2];
		}

		Vector3d bspline_surface_value;
		bspline_surface_value.setZero();
		Vector3d tangent; 
		tangent.setZero();
		if (u_interval_order == 0)
		{
			Vector3d pi0;
			pi0 = P[0];
			bspline_surface_value = pi0;
			tangent = -3.0 / (u_knot_temp[poly_degree + 1]-u_knot_temp[0]) * pi0;
		}
		else if (u_interval_order == u_knot_temp.size() - 1)
		{
			Vector3d pi0;
			pi0 = P[ctr_num - 1];
			bspline_surface_value = pi0;
			tangent = 3.0 / (u_knot_temp[u_knot_temp.size() - 1] - u_knot_temp[u_knot_temp.size()-5]) * pi0;
		}
		else
		{
			for (int i = 0; i < u_basis_temp.size(); i++)
			{
				Vector3d pi0;
				pi0 = P[u_interval_order - 3 + i];
				bspline_surface_value += u_basis_temp[i] * pi0;
				if (i == 0)
				{
					tangent += -3.0 * pow(u_knot_temp[u_interval_order + 1] - res[k], 2) / ((u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 2]))*pi0;
				}
				if (i == 1)
				{
					tangent += ((pow(u_knot_temp[u_interval_order + 1] - res[k], 2) - 2.0 * (u_knot_temp[u_interval_order + 1] - res[k]) * (res[k] - u_knot_temp[u_interval_order - 2])) / ((u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 2]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order])) + ((u_knot_temp[u_interval_order + 2]-res[k])*(u_knot_temp[u_interval_order + 1]-res[k])-(res[k]- u_knot_temp[u_interval_order -1])*(u_knot_temp[u_interval_order + 1]-res[k])-(res[k]- u_knot_temp[u_interval_order -1])*(u_knot_temp[u_interval_order + 2]-res[k])) / ((u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order])) + (pow(u_knot_temp[u_interval_order + 2] - res[k], 2) - 2.0 * (u_knot_temp[u_interval_order + 2] - res[k]) * (res[k] - u_knot_temp[u_interval_order])) / ((u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order]))) * pi0;
				}
				if (i == 2)
				{
					tangent += ((2.0 * (res[k] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - res[k]) - pow(res[k] - u_knot_temp[u_interval_order - 1], 2)) / ((u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order])) + ((u_knot_temp[u_interval_order + 2] - res[k]) * (res[k] - u_knot_temp[u_interval_order]) + (res[k] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 2] - res[k]) - (res[k] - u_knot_temp[u_interval_order - 1]) * (res[k] - u_knot_temp[u_interval_order])) / ((u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order - 1]) * (u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order])) + (2.0 * (res[k] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 3] - res[k]) - pow(res[k] - u_knot_temp[u_interval_order], 2)) / ((u_knot_temp[u_interval_order + 3] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order]))) * pi0;
				}
				if (i == 3)
				{
					tangent += (3.0 * pow(res[k] - u_knot_temp[u_interval_order], 2) / ((u_knot_temp[u_interval_order + 3] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 2] - u_knot_temp[u_interval_order]) * (u_knot_temp[u_interval_order + 1] - u_knot_temp[u_interval_order]))) * pi0;
				}
			}
			
		}

		SampleCurve sc;
		sc.p = bspline_surface_value;
		sc.u = res[k];
		sc.tangent_vector = tangent; 
		sc.tangent_vector.normalize();
		sc.u_interval_order = u_interval_order;
		sc.u_basis_temp = u_basis_temp;
		auto mini_tan = bspline_surface_value+tangent * 0.01;
		ioss1 << bspline_surface_value[0] << " " << bspline_surface_value[1] << " " << bspline_surface_value[2] << " " << mini_tan[0] << " " << mini_tan[1] << " " << mini_tan[2] << endl;
		Point p_close(bspline_surface_value[0], bspline_surface_value[1], bspline_surface_value[2]);
		auto result = tree.closest_point_and_primitive(p_close);

		sc.normal_vector = { result.second->my_face_normal[0].value().value(), result.second->my_face_normal[1].value().value() , result.second->my_face_normal[2].value().value() };
		sc.normal_vector.normalize();
		sc.rulling_vector = sc.normal_vector.cross(sc.tangent_vector);
		sc.rulling_vector.normalize();
		sample_curve_temp.push_back(sc);

		ioss << sc.p[0] << " " << sc.p[1] << " " << sc.p[2] << endl;
		//ioss << sc.tangent_vector[0] << " " << sc.tangent_vector[1] << " " << sc.tangent_vector[2] << endl;
		//ioss << sc.normal_vector[0] << " " << sc.normal_vector[1] << " " << sc.normal_vector[2] << endl;
		ioss << sc.rulling_vector[0] << " " << sc.rulling_vector[1] << " " << sc.rulling_vector[2] << endl;
	}
	ioss.close();
	ioss1.close();
	//sample_curve_temp.pop_back();

	/*
	ioss << sample_curve_temp[0].p[0] << " " << sample_curve_temp[0].p[1] << " " << sample_curve_temp[0].p[2] << endl;
	//ioss << sample_curve_temp[0].tangent_vector[0] << " " << sample_curve_temp[0].tangent_vector[1] << " " << sample_curve_temp[0].tangent_vector[2] << endl;
	//ioss << sample_curve_temp[0].normal_vector[0] << " " << sample_curve_temp[0].normal_vector[1] << " " << sample_curve_temp[0].normal_vector[2] << endl;
	ioss << sample_curve_temp[0].rulling_vector[0] << " " << sample_curve_temp[0].rulling_vector[1] << " " << sample_curve_temp[0].rulling_vector[2] << endl;
	ioss.close();
	*/

	//curve2surface_init(curve);
	std::vector<T> trip;
	m_b.resize(3 * sample_curve_temp.size());// m_b.setConstant(rulings_length);
	m_B.resize(curve.size()/3);
	for (int i = 0; i < sample_curve_temp.size(); i++)
	{
		m_b[3 * i] = sample_curve_temp[i].rulling_vector[0];
		m_b[3 * i + 1] = sample_curve_temp[i].rulling_vector[1];
		m_b[3 * i + 2] = sample_curve_temp[i].rulling_vector[2];
		if (sample_curve_temp[i].u_interval_order == 0)
		{
			trip.emplace_back(3 * i, 0, 1);
			trip.emplace_back(3 * i+1, 1, 1);
			trip.emplace_back(3 * i + 2, 2, 1);
		}
		else if (sample_curve_temp[i].u_interval_order == u_knot_temp.size() - 1)
		{
			trip.emplace_back(3 * i, curve.size() - 3, 1);
			trip.emplace_back(3 * i + 1, curve.size() - 2, 1);
			trip.emplace_back(3 * i + 2, curve.size() - 1, 1);
		}
		else
		{
			for (int j = 0; j < sample_curve_temp[i].u_basis_temp.size(); j++)
			{
				trip.emplace_back(3 * i, 3 * (sample_curve_temp[i].u_interval_order - 3 + j), sample_curve_temp[i].u_basis_temp[j]);
				trip.emplace_back(3 * i + 1, 1 + 3 * (sample_curve_temp[i].u_interval_order - 3 + j), sample_curve_temp[i].u_basis_temp[j]);
				trip.emplace_back(3 * i + 2, 2 + 3 * (sample_curve_temp[i].u_interval_order - 3 + j), sample_curve_temp[i].u_basis_temp[j]);
			}
		}
	}

	std::vector<T> trip2; trip2.clear();
	for (int i = 0; i < curve.size()/3-1; i++)
	{
		trip2.emplace_back(3 * i, 3 * i, -1);
		trip2.emplace_back(3 * i + 1, 3 * i + 1, -1);
		trip2.emplace_back(3 * i + 2, 3 * i + 2, -1);
		trip2.emplace_back(3 * i, 3 * (i + 1), 1);
		trip2.emplace_back(3 * i + 1, 3 * (i + 1) + 1, 1);
		trip2.emplace_back(3 * i + 2, 3 * (i + 1) + 2, 1);
	}

	m_P.resize(3 * (curve.size()/3-1), 3 * curve.size()/3);
	m_P.setFromTriplets(trip2.begin(), trip2.end());

	//cout << rulings_length << endl;
	//rulings_length = 5;
	m_b = m_b * rulings_length;
	m_L.resize(3 * sample_curve_temp.size(), 3 * curve.size()/3);
	m_L.setFromTriplets(trip.begin(), trip.end());


	double lamda = 50;
	solver.compute(m_L.transpose() * m_L + lamda * m_P.transpose() * m_P);
	m_x = solver.solve(m_L.transpose() * m_b);

	std::cout << "能量值： " << (m_L * m_x - m_b).norm() + lamda * (m_P * m_x).norm() << endl;
	std::cout << "残差： " << ((m_L.transpose() * m_L + lamda * m_P.transpose() * m_P) * m_x - m_L.transpose() * m_b).norm() << endl;

	vector<Vector3d> rulings_d; rulings_d.resize(curve.size()/3);
	for (int i = 0; i < (curve.size() / 3); i++)
	{
		Vector3d p = { m_x[3 * i],m_x[3 * i + 1] ,m_x[3 * i + 2] }; //p.normalize();
		rulings_d[i] = p; //cout << p << endl;
		//cout << p << endl;
	}


	Data spline_; spline_.resize(2 * (curve.size() / 3), 3);
	for (int i = 0; i < curve.size() / 3; i++)
	{
		spline_(i, 0) = curve[3 * i] + rulings_d[i][0] * 0.5 * rulings_length;
		spline_(i, 1) = curve[3 * i + 1] + rulings_d[i][1] * 0.5 * rulings_length;
		spline_(i, 2) = curve[3 * i + 2] + rulings_d[i][2] * 0.5 * rulings_length;
		spline_(curve.size() / 3 + i, 0) = curve[3 * i] - rulings_d[i][0] * 0.5 * rulings_length;
		spline_(curve.size() / 3 + i, 1) = curve[3 * i + 1] - rulings_d[i][1] * 0.5 * rulings_length;
		spline_(curve.size() / 3 + i, 2) = curve[3 * i + 2] - rulings_d[i][2] * 0.5 * rulings_length;
	}
	//BSpline_surface_viewer_2(spline_, 100, 0);
	int cp_row = spline_.rows() / 2;
	cp_.resize(cp_row);
	for (int i = 0; i < cp_row; i++)
	{
		cp_[i].resize(2);
		cp_[i][0] = spline_.row(i);
		cp_[i][1] = spline_.row(i + cp_row);
	}

	//draw_surface();
	std::vector<vector<Vector3d>> P1;
	P1.resize(spline_.rows()/2);
	for (int i = 0; i < P1.size(); i++)
	{
		P1[i].resize(2);
		P1[i][0][0] = spline_.coeffRef(i, 0);
		P1[i][0][1] = spline_.coeffRef(i, 1);
		P1[i][0][2] = spline_.coeffRef(i, 2);
		P1[i][1][0] = spline_.coeffRef(i + P1.size(), 0);
		P1[i][1][1] = spline_.coeffRef(i + P1.size(), 1);
		P1[i][1][2] = spline_.coeffRef(i + P1.size(), 2);
	}
	sample_surface_temp.clear();
	for (int k = 0; k < res.size(); k++)
	{
		SampleCurve sc;
		Vector3d p_mid; p_mid.setZero();
		Vector3d rule_v; rule_v.setZero();
		sc.u = sample_curve_temp[k].u;
		sc.u_interval_order = sample_curve_temp[k].u_interval_order;
		sc.u_basis_temp = sample_curve_temp[k].u_basis_temp;
		if (sc.u_interval_order == 0)
		{
			p_mid = (P1[0][0] + P1[0][1]) / 2.0;
			rule_v = (P1[0][0] - P1[0][1]);
		}
		else if (sc.u_interval_order == u_knot_temp.size() - 1)
		{
			p_mid = (P1[P1.size() - 1][0] + P1[P1.size() - 1][1]) / 2.0;
			rule_v = (P1[P1.size() - 1][0] - P1[P1.size() - 1][1]);
		}
		else
		{
			for (int it = 0; it < sc.u_basis_temp.size(); it++)
			{
				p_mid += (P1[sc.u_interval_order - 3 + it][0] + P1[sc.u_interval_order - 3 + it][1]) * 0.5 * sc.u_basis_temp[it];
				rule_v += (P1[sc.u_interval_order - 3 + it][0] - P1[sc.u_interval_order - 3 + it][1]) * sc.u_basis_temp[it];
			}
		}

		sc.p = p_mid;
		sc.u = res[k];
		sc.rulling_vector = rule_v.normalized();
		sample_surface_temp.push_back(sc);
	}

}

void SplineInit::calc_soomth_surface_pca_plane(int sample_num, SegSurface& ss, SegSurface::Path& path)
{

	Vector3d p_mean; p_mean.setZero();
	for (int i = 0; i < sample_surface_temp.size(); i++)
	{
		p_mean += sample_surface_temp[i].p;
	}

	p_mean /= sample_surface_temp.size(); std::cout << "中心点： " << p_mean << endl;
	Matrix3d mat1; mat1.setZero();
	Matrix3d mat2; mat2.setZero();
	Matrix3d mat; mat.setZero();
	for (int i = 0; i < sample_surface_temp.size(); i++)
	{
		mat1 += ((p_mean - sample_surface_temp[i].p).normalized()) * ((p_mean - sample_surface_temp[i].p).normalized()).transpose();
		mat2 += sample_surface_temp[i].rulling_vector * sample_surface_temp[i].rulling_vector.transpose();
	}
	double weight = 0.1;
	mat = weight * mat1 + mat2;
	EigenSolver<Matrix3d> es(mat);
	Matrix3d D = es.pseudoEigenvalueMatrix();
	Matrix3d V = es.pseudoEigenvectors();

	vector<double> evv = { D.coeffRef(0,0),D.coeffRef(1,1), D.coeffRef(2,2) };
	auto min_it = std::min_element(evv.begin(), evv.end());
	int min_index = 0;
	if (min_it != evv.end())
	{
		min_index = std::distance(evv.begin(), min_it);
		std::cout << "最小值为: " << *min_it << " 在索引位置: " << min_index << std::endl;
	}
	else
	{
		std::cout << "向量为空" << std::endl;
	}

	std::cout << "The pseudo-eigenvalue matrix D is:" << endl << D << endl;
	std::cout << "The pseudo-eigenvector matrix V is:" << endl << V << endl;
	std::cout << "Finally, V * D * V^(-1) = " << endl << V * D * V.inverse() << endl;
	std::cout << "最小特征值： " << *min_it << endl;
	std::cout << "最小特征向量： " << V.col(min_index).transpose() << endl;
	Vector3d plane_normal = V.col(min_index);// .transpose();
	double term1 = plane_normal.transpose() * mat1 * plane_normal;
	double term2 = plane_normal.transpose() * mat2 * plane_normal;
	//VectorXd plane_n; plane_.re
	std::cout << "第一项能量： " << weight * term1 << " 第二项能量： " << term2 << endl;


	//构建平面网格
	Mesh::Point p0(p_mean.x(), p_mean.y(), p_mean.z());
	double pv = -p_mean.transpose() * plane_normal;
	Mesh::Point pn(plane_normal[0], plane_normal[1], plane_normal[2]);
	double p_z = (-pv - plane_normal[0] * (p_mean.x() + 1) - plane_normal[1] * (p_mean.y() + 1)) / plane_normal.z();
	Mesh::Point d0(1, 1, p_z - p_mean.z()); 
	d0.normalize();
	Mesh::Point p1 = p0 + d0 * rulings_length * 2;
	Mesh::Point p2 = p0 - d0 * rulings_length * 2;
	Mesh::Point p3 = p0 + d0.cross(pn) * rulings_length * 1;
	Mesh::Point p4 = p0 - d0.cross(pn) * rulings_length * 1;

	Mesh mesh_res;
	auto v0 = mesh_res.add_vertex(p1);
	auto v1 = mesh_res.add_vertex(p2);
	auto v2 = mesh_res.add_vertex(p3);
	auto v3 = mesh_res.add_vertex(p4);
	mesh_res.add_face(v0, v2, v3);
	mesh_res.add_face(v1, v3, v2);

	MeshTools::WriteMesh(mesh_res, "./mesh_res.obj");


	
	curve_pro_point.clear();
	rulings_pro_point.clear();
	Plane plane(plane_normal.x(), plane_normal.y(), plane_normal.z(), pv);
	Vector_3 normal = plane.orthogonal_vector();
	ofstream iospro("pro_point.txt"); ofstream ioslinepro("pro_line.txt");
	for (int i = 0; i < sample_surface_temp.size(); i++)
	{
		Point point_(sample_surface_temp[i].p.x(), sample_surface_temp[i].p.y(), sample_surface_temp[i].p.z());
		Vector_3 pointToPlaneVec(point_ - plane.projection(point_));
		Point projectedPoint = point_ - pointToPlaneVec;
		curve_pro_point.push_back(projectedPoint);
		iospro << projectedPoint.x() << " " << projectedPoint.y() << " " << projectedPoint.z() << endl;

		Vector3d point1 = sample_surface_temp[i].p + sample_surface_temp[i].rulling_vector * 0.5 * rulings_length;
		Vector3d point2 = sample_surface_temp[i].p - sample_surface_temp[i].rulling_vector * 0.5 * rulings_length;

		Point p1(point1.x(), point1.y(), point1.z());
		Point p2(point2.x(), point2.y(), point2.z());
		Vector_3 pointToPlaneVec1(p1 - plane.projection(p1));
		Point projectedPoint1 = p1 - pointToPlaneVec1;
		Vector_3 pointToPlaneVec2(p2 - plane.projection(p2));
		Point projectedPoint2 = p2 - pointToPlaneVec2;
		ioslinepro << projectedPoint1.x() << " " << projectedPoint1.y() << " " << projectedPoint1.z() << endl;
		ioslinepro << projectedPoint2.x() << " " << projectedPoint2.y() << " " << projectedPoint2.z() << endl;
		rulings_pro_point.push_back(make_pair(projectedPoint1, projectedPoint2));
	}
	ioslinepro.close();
	iospro.close();

	Vector3d dir_init; dir_init.setZero();
	for (int i = 0; i < rulings_pro_point.size(); i++)
	{
		Vector3d p1 = { rulings_pro_point[i].first.x(),rulings_pro_point[i].first.y(), rulings_pro_point[i].first.z() };
		Vector3d p2 = { rulings_pro_point[i].second.x(),rulings_pro_point[i].second.y(), rulings_pro_point[i].second.z() };
		dir_init += (p1 - p2).normalized();
	}
	dir_init.normalize();
	std::cout << dir_init << endl;

	double min_translation_length = 2 * average_edge_length;
	calc_spline_init_info(p_mean, plane_normal, dir_init, ss, path, min_translation_length);

}

void SplineInit::add_ctr_point(Eigen::VectorXd& BSpline, double u)
{
	int k = u_k_position(u, u_knot_temp);
	vector<double> a(poly_degree, 0);
	for (int i = 0; i < poly_degree; i++)
	{
		a[i] = (u - u_knot_temp[i + k - poly_degree + 1]) / (u_knot_temp[i + k + 1] - u_knot_temp[i + k - poly_degree + 1]);
	}

	auto ctr_temp = BSpline;
	BSpline.resize(BSpline.size() + 3);
	for (int i = 0; i <= k - poly_degree; i++)
	{
		BSpline[3 * i] = ctr_temp[3 * i];
		BSpline[3 * i + 1] = ctr_temp[3 * i + 1];
		BSpline[3 * i + 2] = ctr_temp[3 * i + 2];
	}

	for (int i = k - poly_degree + 1; i <= k; i++)
	{
		BSpline[3 * i] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1)] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i];
		BSpline[3 * i + 1] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 1] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 1];
		BSpline[3 * i + 2] = (1 - a[i - (k - poly_degree + 1)]) * ctr_temp[3 * (i - 1) + 2] + a[i - (k - poly_degree + 1)] * ctr_temp[3 * i + 2];

	}

	for (int i = k + 1; i < BSpline.size() / 3; i++)
	{
		BSpline[3 * i] = ctr_temp[3 * (i - 1)];
		BSpline[3 * i + 1] = ctr_temp[3 * (i - 1) + 1];
		BSpline[3 * i + 2] = ctr_temp[3 * (i - 1) + 2];

	}

	auto knot_temp = u_knot_temp;
	u_knot_temp.clear();
	u_knot_temp.resize(knot_temp.size() + 1);
	for (int i = 0; i <= k; i++)
	{
		u_knot_temp[i] = knot_temp[i];
	}
	u_knot_temp[k + 1] = u;
	for (int i = k + 2; i < u_knot_temp.size(); i++)
	{
		u_knot_temp[i] = knot_temp[i - 1];
	}
}

std::vector<double> SplineInit::uni_para(int samp_num)
{
	vector<double> res;
	res.resize((u_knot_temp.size() - poly_degree * 2 - 1) * samp_num+1);
	for (int i = 0; i < (u_knot_temp.size() - poly_degree * 2 - 1); i++)
	{
		for (int j = 0; j < samp_num; j++)
		{
			res[i * samp_num + j] = (u_knot_temp[i +1+ poly_degree] * (j)+u_knot_temp[i  + poly_degree] * (samp_num - j)) / samp_num;
		}
	}
	res[(u_knot_temp.size() - poly_degree * 2 - 1)* samp_num] = u_knot_temp.back();
	return res;
}


void SplineInit::BSpline_surface_viewer_2(const Data& spline, int sample_num, int out_number)
{
	VectorXd BSpline;
	BSpline.resize(spline.rows() * 3);
	for (int i = 0; i < spline.rows(); i++)
	{

		BSpline[3 * i] = spline.coeff(i, 0);
		BSpline[3 * i + 1] = spline.coeff(i, 1);
		BSpline[3 * i + 2] = spline.coeff(i, 2);
	}


	vector<vector<CtrInfo>> ctr_point_temp;
	int curve_ctr_num = BSpline.size() / 6;

	//for (int i = 0; i < 2; i++)
	//{
	vector<CtrInfo> ctr_point1;
	for (int j = 0; j < curve_ctr_num; j++)
	{
		CtrInfo ci;
		Standard_Real px = BSpline[3 * j];
		Standard_Real py = BSpline[3 * j + 1];
		Standard_Real pz = BSpline[3 * j + 2];
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point1.push_back(ci);
	}
	ctr_point_temp.push_back(ctr_point1);

	vector<CtrInfo> ctr_point2;
	for (int j = 0; j < curve_ctr_num; j++)
	{
		CtrInfo ci;
		Standard_Real px = BSpline[3 * (j + curve_ctr_num)];
		Standard_Real py = BSpline[3 * (j + curve_ctr_num) + 1];
		Standard_Real pz = BSpline[3 * (j + curve_ctr_num) + 2];
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point2.push_back(ci);
	}
	ctr_point_temp.push_back(ctr_point2);
	//}

	Handle(Geom_BSplineSurface) surface = nullptr;
	create_BSpline_surface(ctr_point_temp[0], ctr_point_temp[1], surface);
	//interval_up_bound = divide_u_knot[2][0].back();
	//cout << interval_up_bound << endl;
	vector<Standard_Real> u_temp;// = generate_equally_spaced_vector(sample_num, interval_up_bound);
	for (int i = 0; i < u_knot_temp.size() - 1; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			Standard_Real val = u_knot_temp[i] + j / 5.0 * (u_knot_temp[i + 1] - u_knot_temp[i]);
			u_temp.push_back(val);
		}
	}
	u_temp.push_back(u_knot_temp.back());
	vector<Standard_Real> v_temp = generate_equally_spaced_vector(100);

	//2.Bezier_info构造
	Matrix<SplineInit::BSplinePoint, Eigen::Dynamic, Eigen::Dynamic> bsp_temp;
	bsp_temp.resize(u_temp.size(), v_temp.size());
	for (int i = 0; i < u_temp.size(); i++)
	{
		for (int j = 0; j < v_temp.size(); j++)
		{
			BSplinePoint bsp;
			bsp.u = u_temp[i];
			bsp.v = v_temp[j];
			surface->D0(bsp.u, bsp.v, bsp.position);
			bsp.bsp_position[0] = bsp.position.X();
			bsp.bsp_position[1] = bsp.position.Y();
			bsp.bsp_position[2] = bsp.position.Z();
			bsp_temp(i, j) = bsp;
		}
	}

	//3.
	Mesh outmesh;
	Matrix<int, Eigen::Dynamic, Eigen::Dynamic> v_idx_mat;
	v_idx_mat.resize(u_temp.size(), v_temp.size());
	int init_v_idx = 0;
	for (int i = 0; i < bsp_temp.rows(); i++)
	{
		for (int j = 0; j < bsp_temp.cols(); j++)
		{
			Mesh::VertexHandle v0 = outmesh.add_vertex(Mesh::Point(bsp_temp(i, j).bsp_position.data()));
			v_idx_mat(i, j) = init_v_idx;
			init_v_idx++;
		}
	}

	for (int i = 0; i < u_temp.size() - 1; ++i)
	{
		for (int j = 0; j < v_temp.size() - 1; ++j)
		{
			Mesh::VertexHandle v0 = outmesh.vertex_handle(v_idx_mat(i, j));
			Mesh::VertexHandle v1 = outmesh.vertex_handle(v_idx_mat(i + 1, j));
			Mesh::VertexHandle v2 = outmesh.vertex_handle(v_idx_mat(i + 1, j + 1));
			Mesh::VertexHandle v3 = outmesh.vertex_handle(v_idx_mat(i, j + 1));

			std::vector<Mesh::VertexHandle> face_vhandles0;
			face_vhandles0.push_back(v0);
			face_vhandles0.push_back(v1);
			face_vhandles0.push_back(v2);
			outmesh.add_face(face_vhandles0);

			std::vector<Mesh::VertexHandle> face_vhandles1;
			face_vhandles1.push_back(v0);
			face_vhandles1.push_back(v2);
			face_vhandles1.push_back(v3);
			outmesh.add_face(face_vhandles1);
		}
	}

	//4.构造并输出网格
	MeshTools::WriteMesh(outmesh, "./out_" + to_string(out_number) + ".obj");
	std::cout << "success" << endl << endl;
}

void SplineInit::create_BSpline_surface(vector<CtrInfo>& ctr_point1, vector<CtrInfo>& ctr_point2, Handle(Geom_BSplineSurface)& surface)
{
	Handle(Geom_BSplineCurve) curve1 = nullptr;
	create_BSpline_curve(ctr_point1, curve1);
	Handle(Geom_BSplineCurve) curve2 = nullptr;
	create_BSpline_curve(ctr_point2, curve2);

	GeomFill_FillingStyle Type = GeomFill_StretchStyle;
	GeomFill_BSplineCurves aGeomFill(curve1, curve2, Type);
	surface = aGeomFill.Surface();
}

void SplineInit::create_BSpline_curve(vector<CtrInfo>& ctr_point, Handle(Geom_BSplineCurve)& curve)
{
	int ctr_num = ctr_point.size();
	TColgp_Array1OfPnt Poles(1, ctr_num);
	TColStd_Array1OfReal PolesWeight(1, ctr_num);
	for (int i = 0; i < ctr_num; i++)
	{
		Poles.SetValue(i + 1, ctr_point[i].position);
		PolesWeight.SetValue(i + 1, ctr_point[i].weight);
	}

	Standard_Integer PNum = ctr_point.size();

	Standard_Integer KNum = PNum - poly_degree + 1;

	TColStd_Array1OfReal knots(1, KNum);
	TColStd_Array1OfInteger mults(1, KNum);
	for (int i = 0; i < KNum; ++i)
	{
		//cout << i <<" " << u_knot_temp[i + poly_degree] << "success" << endl;
		knots.SetValue(i + 1, u_knot_temp[i + poly_degree]);

		if (i == 0 || i == KNum - 1)
		{
			mults.SetValue(i + 1, poly_degree + 1);
		}
		else
		{
			mults.SetValue(i + 1, 1);
		}

	}

	curve = new Geom_BSplineCurve(Poles, PolesWeight, knots, mults, poly_degree);

}