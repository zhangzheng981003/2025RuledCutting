#include "CurveFitting.h"

CurveFitting::CurveFitting(Mesh& mesh) :mesh(mesh)
{

}

int CurveFitting::combination(int n, int m)
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

void CurveFitting::initialization()
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
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()));

		my_cgal_triangles.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()), f_it->idx(), f_n);
		v_p_3.clear();
	}

	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();

	my_tree.insert(my_cgal_triangles.begin(), my_cgal_triangles.end());
	my_tree.build();
	my_tree.accelerate_distance_queries();

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

void CurveFitting::curve_init(int u_degree, VectorXd& curve, vector<Vector3d>& point_cloud)
{
	int gap = std::round(point_cloud.size() / (u_degree + 1));
	curve.resize((u_degree + 1) * 3);
	curve.setZero();
	for (int i = 0; i < (u_degree + 1); i++)
	{
		if (i != u_degree)
		{
			curve[3 * i + 0] = point_cloud[i * gap].x();
			curve[3 * i + 1] = point_cloud[i * gap].y();
			curve[3 * i + 2] = point_cloud[i * gap].x();
		}
		else
		{
			curve[3 * i + 0] = point_cloud[point_cloud.size() - 1].x();
			curve[3 * i + 1] = point_cloud[point_cloud.size() - 1].y();
			curve[3 * i + 2] = point_cloud[point_cloud.size() - 1].x();
		}
	}

	point_cloud_temp.clear();
	for (int i = 0; i < point_cloud.size(); i++)
	{
		PointCloud pc;
		pc.index = i;
		pc.position = point_cloud[i];
		pc.autodiff_position[0] = pc.position[0];
		pc.autodiff_position[1] = pc.position[1];
		pc.autodiff_position[2] = pc.position[2];
		point_cloud_temp.push_back(pc);
	}

}

void CurveFitting::calc_projection_distance2_gradient(int u_degree, VectorXd& spline, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, Eigen::VectorXd& decrease_direction)
{
	data_value = 0;
	grad.resize(3 * (u_degree + 1));
	grad.setZero();
	hessian.resize(3 * (u_degree + 1), 3 * (u_degree + 1));
	hessian.setZero();

	TColgp_Array1OfPnt controlPoints(1, u_degree + 1);
	for (int i = 0; i < u_degree + 1; i++)
	{
		controlPoints.SetValue(i + 1, gp_Pnt(spline[3 * i], spline[3 * i + 1], spline[3 * i + 2]));
	}
	//Standard_Real tola = 1e-6;
	Handle(Geom_Curve) bezierCurve = new Geom_BezierCurve(controlPoints);

	for (int k = 0; k < point_cloud_temp.size(); k++)
	{
		// 定义要计算投影距离的点
		gp_Pnt pointToProject(point_cloud_temp[k].position.x(), point_cloud_temp[k].position.y(), point_cloud_temp[k].position.z());
		// 计算点到贝塞尔曲线的投影
		GeomAPI_ProjectPointOnCurve projectPointOnCurve(pointToProject, bezierCurve);

		if (projectPointOnCurve.NbPoints() > 0)
		{
			// 获取投影点
			Standard_Real u = projectPointOnCurve.LowerDistanceParameter();
			point_cloud_temp[k].u = u;
			Vec12 X;
			X.resize(3 * (u_degree + 1));
			for (int r = 0; r < (u_degree + 1); r++)
			{
				X(3 * r).value() = spline(3 * r);
				X(3 * r + 1).value() = spline(3 * r + 1);
				X(3 * r + 2).value() = spline(3 * r + 2);
			}

			// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
			for (int id = 0; id < 3 * (u_degree + 1); id++)
			{
				X(id).derivatives().resize(3 * (u_degree + 1));
				X(id).derivatives().setZero();
				X(id).derivatives()(id) = 1;
				X(id).value().derivatives() = inner_derivative_t::Unit(3 * (u_degree + 1), id);
			}

			for (int idx = 0; idx < 3 * (u_degree + 1); idx++) {
				for (int id = 0; id < 3 * (u_degree + 1); id++)
				{
					X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (u_degree + 1));
				}
			}


			std::vector<Vec3> P;
			P.resize((u_degree + 1));

			for (int i = 0; i < P.size(); i++)
			{
				P[i][0] = X[3 * i];
				P[i][1] = X[3 * i + 1];
				P[i][2] = X[3 * i + 2];
			}

			Vec3 bezier_surface_value;
			bezier_surface_value.setZero();
			for (int i = 0; i <= u_degree; i++)
			{
				Vec3 pi;
				pi = P[i];
				bezier_surface_value[0] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[0];
				bezier_surface_value[1] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[1];
				bezier_surface_value[2] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[2];
			}
			scalar_t e = (bezier_surface_value - point_cloud_temp[k].autodiff_position) * (bezier_surface_value - point_cloud_temp[k].autodiff_position).transpose();
			data_value += e;
			grad += e.value().derivatives();

			Eigen::MatrixXd B;
			B.resize(3 * (u_degree + 1), 3 * (u_degree + 1));
			for (int r = 0; r < 3 * (u_degree + 1); r++)
			{
				B.row(r) = e.derivatives()(r).derivatives().transpose();
			}
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
	decrease_direction.normalize();
}


double CurveFitting::calc_energy(int u_degree, VectorXd& spline)
{
	double data_value = 0.0;
	for (int k = 0; k < point_cloud_temp.size(); k++)
	{
		std::vector<Vector3d> P;
		P.resize((u_degree + 1));
		for (int i = 0; i < P.size(); i++)
		{
			P[i][0] = spline[3 * i];
			P[i][1] = spline[3 * i + 1];
			P[i][2] = spline[3 * i + 2];
		}

		Vector3d bezier_surface_value;
		bezier_surface_value.setZero();
		for (int i = 0; i <= u_degree; i++)
		{
			Vector3d pi;
			pi = P[i];
			bezier_surface_value[0] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[0];
			bezier_surface_value[1] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[1];
			bezier_surface_value[2] += combination(u_degree, i) * pow((1 - point_cloud_temp[k].u), u_degree - i) * pow(point_cloud_temp[k].u, i) * pi[2];
		}
		double e = (bezier_surface_value - point_cloud_temp[k].position).transpose() * (bezier_surface_value - point_cloud_temp[k].position);
		data_value += e;
	}
	return data_value;
}


double CurveFitting::line_search(int u_degree,VectorXd& spline, double data_value, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction)
{

	double step_size = 1.0;// Initial step size
	int iter_num = 0;
	while (iter_num < 50)
	{
		Eigen::VectorXd new_spline = spline - step_size * grad;
		double new_data_value = calc_energy(u_degree, new_spline);
		cout << endl;
		cout << "setp: " << step_size << endl;
		cout << "data_value_new: " << new_data_value << endl;
		cout << "old_total: " << data_value << endl; //data_value.value() <<  endl;

		if (new_data_value > data_value + c_parm * step_size * decrease_direction.transpose() * grad)
		{
			step_size *= beta;
		}
		else
		{
			return step_size;
			break;
		}

		iter_num++;
		if (iter_num == 50)
		{
			step_size = 0;
			return step_size;
		}
	}
}

void CurveFitting::energy_opt(int u_degree, VectorXd& spline)
{
	int iter_num = 0;
	while (true)
	{
		Eigen::VectorXd grad, new_spline, derection;
		scalar_t data_value;
		Data hessien;
		//cout << spline.size() << endl;
		calc_projection_distance2_gradient(u_degree, spline, data_value, grad, hessien, derection);
		grad.normalize();
		//cout << "梯度："<<grad << endl;
		double data_term = data_value.value().value();
		double step = line_search(u_degree, spline, data_term, grad, derection);
		//Eigen::VectorXd direction;
		new_spline = spline + step * derection;
		double termination_condition = (new_spline - spline).squaredNorm();

		spline = new_spline;
		iter_num++;
		data_energy = calc_energy(u_degree, spline);
		cout << "第 " << iter_num << " 迭代：" << data_energy << endl;

		cout << termination_condition << endl;
		if (termination_condition < stop_threshold || iter_num>500)
		{
			break;
		}
	}

	/*
	control_points.resize(u_degree + 1, 3);
	for (int i = 0; i < control_points.rows(); i++)
	{
		control_points(i, 0) = spline[3 * i];
		control_points(i, 1) = spline[3 * i + 1];
		control_points(i, 2) = spline[3 * i + 2];
	}

	vector<double> res = generate_equally_spaced_vector(1000);
	for (int i = 0; i < res.size(); i++)
	{
		VectorXd p0; p0.setZero(3);
		for (int j = 0; j < control_points.rows(); j++)
		{
			p0[0] += combination(u_degree, j) * pow(res[i], j) * pow(1 - res[i], u_degree - j) * control_points.coeff(j, 0);
			p0[1] += combination(u_degree, j) * pow(res[i], j) * pow(1 - res[i], u_degree - j) * control_points.coeff(j, 1);
			p0[2] += combination(u_degree, j) * pow(res[i], j) * pow(1 - res[i], u_degree - j) * control_points.coeff(j, 2);
		}
	}*/
}

void CurveFitting::curve_fitting(int u_degree, VectorXd& curve, vector<Vector3d>& point_cloud)
{
	curve_init(u_degree, curve, point_cloud);
	energy_opt(u_degree, curve);
}

std::vector<double> CurveFitting::generate_equally_spaced_vector(int num_elements)//yes
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

void CurveFitting::calc_soomth_surface(VectorXd& curve, int sample_num)
{
	vector<double> res = generate_equally_spaced_vector(sample_num);//yes
	sample_curve_temp.clear();
	int order = (curve.size() / 3) - 1;
	for (int k = 0; k < res.size(); k++)
	{
		Vector3d curve_p; curve_p.setZero();
		for (int i = 0; i < order + 1; i++)
		{
			Vector3d pi = { curve[3 * i],curve[3 * i + 1] ,curve[3 * i + 2] };
			curve_p += combination(order, i) * (pow(1 - res[k], order - i) * pow(res[k], i)) * pi;
		}

		SplineInit::SampleCurve sc;
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
		auto result = my_tree.closest_point_and_primitive(p_close);

		sc.normal_vector = { result.second->my_face_normal[0].value().value(), result.second->my_face_normal[1].value().value() , result.second->my_face_normal[2].value().value() };
		sc.normal_vector.normalize();
		sc.rulling_vector = sc.normal_vector.cross(sc.tangent_vector);
		sc.rulling_vector.normalize();
		sample_curve_temp.push_back(sc);
	}

	std::vector<T> trip;
	m_b.resize(3 * sample_curve_temp.size());
	m_B.resize(order + 1);
	for (int i = 0; i < sample_curve_temp.size(); i++)
	{
		m_b[3 * i] = sample_curve_temp[i].rulling_vector[0];
		m_b[3 * i + 1] = sample_curve_temp[i].rulling_vector[1];
		m_b[3 * i + 2] = sample_curve_temp[i].rulling_vector[2];
		for (int j = 0; j < order + 1; j++)
		{
			trip.emplace_back(3 * i, 3 * j, combination(order, j) * pow(sample_curve_temp[i].u, j) * pow(1 - sample_curve_temp[i].u, order - j));
			trip.emplace_back(3 * i + 1, 3 * j + 1, combination(order, j) * pow(sample_curve_temp[i].u, j) * pow(1 - sample_curve_temp[i].u, order - j));
			trip.emplace_back(3 * i + 2, 3 * j + 2, combination(order, j) * pow(sample_curve_temp[i].u, j) * pow(1 - sample_curve_temp[i].u, order - j));
		}
	}

	std::vector<T> trip2; trip2.clear();
	for (int i = 0; i < order; i++)
	{
		trip2.emplace_back(3 * i, 3 * i, -1);
		trip2.emplace_back(3 * i + 1, 3 * i + 1, -1);
		trip2.emplace_back(3 * i + 2, 3 * i + 2, -1);
		trip2.emplace_back(3 * i, 3 * (i + 1), 1);
		trip2.emplace_back(3 * i + 1, 3 * (i + 1) + 1, 1);
		trip2.emplace_back(3 * i + 2, 3 * (i + 1) + 2, 1);
	}

	m_P.resize(3 * order, 3 * (order + 1));
	m_P.setFromTriplets(trip2.begin(), trip2.end());

	m_b = m_b * rulings_length;
	m_L.resize(3 * sample_curve_temp.size(), 3 * (order + 1));
	m_L.setFromTriplets(trip.begin(), trip.end());


	double lamda = 50;
	solver.compute(m_L.transpose() * m_L + lamda * m_P.transpose() * m_P);
	m_x = solver.solve(m_L.transpose() * m_b);

	cout << "能量值： " << (m_L * m_x - m_b).norm() + lamda * (m_P * m_x).norm() << endl;
	cout << "残差： " << ((m_L.transpose() * m_L + lamda * m_P.transpose() * m_P) * m_x - m_L.transpose() * m_b).norm() << endl;

	vector<Vector3d> rulings_d; rulings_d.resize((order + 1));
	for (int i = 0; i < (order + 1); i++)
	{
		Vector3d p = { m_x[3 * i],m_x[3 * i + 1] ,m_x[3 * i + 2] }; 
		p.normalize();
		rulings_d[i] = p; 
	}


	Data spline_; spline_.resize(2 * (order + 1), 3);
	for (int i = 0; i < order + 1; i++)
	{
		spline_(i, 0) = curve[3 * i] + rulings_d[i][0] * 0.5 * rulings_length;
		spline_(i, 1) = curve[3 * i + 1] + rulings_d[i][1] * 0.5 * rulings_length;
		spline_(i, 2) = curve[3 * i + 2] + rulings_d[i][2] * 0.5 * rulings_length;
		spline_(order + 1 + i, 0) = curve[3 * i] - rulings_d[i][0] * 0.5 * rulings_length;
		spline_(order + 1 + i, 1) = curve[3 * i + 1] - rulings_d[i][1] * 0.5 * rulings_length;
		spline_(order + 1 + i, 2) = curve[3 * i + 2] - rulings_d[i][2] * 0.5 * rulings_length;
	}

	vector<vector<Vector3d>> control_points; control_points.clear();
	control_points.resize((order + 1));
	for (int i = 0; i < (order + 1); i++)
	{
		control_points[i].resize(2);
		control_points[i][0] = spline_.row(i);
		control_points[i][1] = spline_.row(i + (order + 1));
	}

	TColgp_Array2OfPnt controlPoints(1, order + 1, 1, 1 + 1);
	for (int i = 0; i < order + 1; i++)
	{
		controlPoints.SetValue(i + 1, 1, gp_Pnt(spline_(i, 0), spline_(i, 1), spline_(i, 2)));
		controlPoints.SetValue(i + 1, 2, gp_Pnt(spline_(order + 1 + i, 0), spline_(order + 1 + i, 1), spline_(order + 1 + i, 2)));
	}
	Standard_Real tola = 1e-6;
	Handle(Geom_BezierSurface) bezierSurface = new Geom_BezierSurface(controlPoints);

	bezier_point_temp.clear();
	for (int k = 0; k < res.size(); k++)
	{
		BezierPoint bp;
		bp.index = k;
		bp.u = res[k];

		Standard_Real u = res[k]; // 曲面参数 u
		Standard_Real v = 0.5; // 曲面参数 v

		// 获取曲面上参数点的切向量和法向量
		gp_Pnt point;
		gp_Dir tangentU, tangentV;
		gp_Dir normal;

		// 计算曲面上点的位置
		bezierSurface->D0(u, v, point);
		
		// 使用 LProp_SLProps 计算切向量和法向量
		GeomLProp_SLProps properties(bezierSurface, u, v, 1, 0.01);
		if (properties.IsNormalDefined() && properties.IsTangentUDefined() && properties.IsTangentVDefined()) 
		{
			properties.TangentU(tangentU);
			properties.TangentV(tangentV);
			normal = properties.Normal();

			bp.position = { point.X(), point.Y(),  point.Z() };
			bp.tangent_direction = { tangentV.X(),tangentV.Y() ,tangentV.Z() };
			bp.normal_direction = { normal.X(),normal.Y() ,normal.Z() };
			bp.rulling_direction = bp.tangent_direction.cross(bp.normal_direction);
			bp.rulling_direction.normalize();

			// 移动到表面 offset
			bezier_point_temp.push_back(bp);
		}

		// 打印结果
		//std::cout << "点位置: " << point.X() << ", " << point.Y() << ", " << point.Z() << std::endl;
		//std::cout << "切向量 U: " << tangentU.X() << ", " << tangentU.Y() << ", " << tangentU.Z() << std::endl;
		//std::cout << "切向量 V: " << tangentV.X() << ", " << tangentV.Y() << ", " << tangentV.Z() << std::endl;
		//std::cout << "法向量: " << normal.X() << ", " << normal.Y() << ", " << normal.Z() << std::endl;
	}

	//输出贝塞尔曲面，以及采样点的信息

	/*
	ofstream iossss("spline_surface.txt");
	for (int i = 0; i < 2 * (order + 1); i++)
	{
		iossss << spline_(i, 0) << " " << spline_(i, 1) << " " << spline_(i, 2) << endl;
	}
	iossss.close();

	int cp_row = spline_.rows() / 2;
	cp_.resize(cp_row);
	for (int i = 0; i < cp_row; i++)
	{
		cp_[i].resize(2);
		cp_[i][0] = spline_.row(i);
		cp_[i][1] = spline_.row(i + cp_row);
	}

	//draw_surface();

	sample_surface_temp.clear();
	for (int k = 0; k < res.size(); k++)
	{
		SampleCurve sc;
		Vector3d p_mid; p_mid.setZero();
		Vector3d rule_v; rule_v.setZero();

		for (int i = 0; i < order + 1; i++)
		{
			for (int j = 0; j < 2; j++)
			{
				p_mid += combination(order, i) * pow((1 - res[k]), order - i) * pow(res[k], i) * combination(1, j) * pow(0.5, 1 - j) * pow(0.5, j) * cp_[i][j];
			}
			rule_v = combination(order, i) * pow((1 - res[k]), order - i) * pow(res[k], i) * (cp_[i][0] - cp_[i][1]);
		}

		sc.p = p_mid;
		sc.u = res[k];
		sc.rulling_vector = rule_v.normalized();
		sample_surface_temp.push_back(sc);
	}*/
}


