#include "SplineSurface.h"

using namespace std;

SplineSurface::SplineSurface(Mesh& mesh) :mesh(mesh)
{
	
}

void SplineSurface::init_ready_data(Mesh& mesh)
{
	//1.p2p
	fitting_init(100);
	//cover_area_temp.clear();
	for (int i = 0; i < cover_area_temp.size(); i++)
	{
		for (int j = 0; j < face_point_temp[cover_area_temp[i]].size(); j++)
		{
			Forward_Projection fp;
			fp.input_point[0] = face_point_temp[cover_area_temp[i]][j].x();
			fp.input_point[1] = face_point_temp[cover_area_temp[i]][j].y();
			fp.input_point[2] = face_point_temp[cover_area_temp[i]][j].z();
			all_data_forward_projection.push_back(fp);
		}
	}
	/*
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		MeshTraits::Point v = mesh.point(mesh.vertex_handle(i));
		Point v_cgal(v[0], v[1], v[2]);
		input_mesh_point_init.push_back(v_cgal);
		Forward_Projection fp;
		fp.input_idx = i;
		fp.input_point[0] = v[0];
		fp.input_point[1] = v[1];
		fp.input_point[2] = v[2];
		all_data_forward_projection.push_back(fp);
	}*/

	//1.line2line
	for (Mesh::EdgeIter edgeIter = mesh.edges_begin(); edgeIter != mesh.edges_end(); ++edgeIter) {
		Mesh::HalfedgeHandle heh = mesh.halfedge_handle(*edgeIter, 0);

		Mesh::VertexHandle fromVertex = mesh.from_vertex_handle(heh);
		Mesh::VertexHandle toVertex = mesh.to_vertex_handle(heh);

		MeshTraits::Point v_f = mesh.point(fromVertex);
		Point v_cgal_1(v_f[0], v_f[1], v_f[2]);
		MeshTraits::Point v_t = mesh.point(toVertex);
		Point v_cgal_2(v_t[0], v_t[1], v_t[2]);

		Segment s1(v_cgal_1, v_cgal_2);
		input_mesh_line_init.push_back(s1);
	}

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
		cgal_triangles_.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()));
		input_mesh_triangle_init.push_back(Triangle(v_p_3[0], v_p_3[1], v_p_3[2]));
		v_p_3.clear();
	}

	ordinary_tree.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	ordinary_tree.build();
	ordinary_tree.accelerate_distance_queries();

	u_knot = generate_equally_spaced_vector(U_Degree + 1);
	v_knot = generate_equally_spaced_vector(V_Degree + 1);

	knot_weight.resize((U_Degree + 1), (V_Degree + 1));
	knot_weight.setConstant(0.5);
	for (int i = 1; i < U_Degree; i++)
	{
		for (int j = 0; j < V_Degree + 1; j++)
		{
			knot_weight(i, j) = 1;
		}
	}
}

void SplineSurface::uv_sample(int row_num)
{
	double gap = 1.0 / (double(row_num - 1));
	for (int i = 0; i < row_num; i++)
	{
		for (int j = 0; j < row_num; j++)
		{
			Bezier_info bi;
			double u_cord = gap * double(i);
			double v_cord = gap * double(j);

			if (gap * i > 1)
			{
				u_cord = 1;
			}
			if (gap * j > 1)
			{
				v_cord = 1;
			}
			bi.u_cord = u_cord;
			//cout << u_cord << endl;
			bi.v_cord = v_cord;
			bi.u_idx = i;
			bi.v_idx = j;
			bi_vec.push_back(bi);
		}
	}

}

int SplineSurface::combination(int n, int m)
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

double SplineSurface::calc_bernstein_basis(int n, int m, double value)
{
	return combination(n, m) * pow((1 - value), n - m) * pow(value, m);
}

void SplineSurface::calc_data_term_energy_derivative(Eigen::VectorXd& spline, scalar_t& data_value, scalar_t& smooth_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, Eigen::VectorXd& decrease_direction, Tree& tree, Tree_& tree_)
{
	ofstream iosout("pro_points_" + to_string(calc_iter_num) + "_.txt");
	data_value = 0;
	int num = 3 * (U_Degree + 1) * (V_Degree + 1);//X维数
	grad.resize(num);
	grad.setZero();
	hessian.resize(num, num);
	hessian.setZero();
	
	vector<vector<Eigen::Vector3d>> control_points;
	control_points.resize((U_Degree + 1));
	for (int i = 0; i < control_points.size(); i++)
	{
		control_points[i].resize(2);
	}

	for (int i = 0; i < control_points.size(); i++)
	{
		control_points[i][0][0] = spline[3 * i];
		control_points[i][0][1] = spline[3 * i + 1];
		control_points[i][0][2] = spline[3 * i + 2];
		control_points[i][1][0] = spline[(U_Degree + 1) * 3 + 3 * i];
		control_points[i][1][1] = spline[(U_Degree + 1) * 3 + 3 * i + 1];
		control_points[i][1][2] = spline[(U_Degree + 1) * 3 + 3 * i + 2];
	}

	//
	TColgp_Array2OfPnt controlPoints(1, U_Degree + 1, 1, V_Degree + 1);
	for (int i = 0; i < U_Degree + 1; i++)
	{
		controlPoints.SetValue(i + 1, 1, gp_Pnt(spline[3 * i], spline[3 * i + 1], spline[3 * i + 2]));
		controlPoints.SetValue(i + 1, 2, gp_Pnt(spline[3 * i + 3 * (U_Degree + 1)], spline[3 * i + 1 + 3 * (U_Degree + 1)], spline[3 * i + 2 + 3 * (U_Degree + 1)]));
	}


	
	Standard_Real tola = 1e-6;
	//定义贝塞尔曲面
	Handle(Geom_BezierSurface) bezierSurface = new Geom_BezierSurface(controlPoints);

	//ofstream ioss("search_points.txt");

	int nv = mesh.n_vertices();
	int ornum = 0;
	double maxdis = 0;
	
	scalar_t data_mesh2surface = 0;
	MatrixXd data_mesh2surface_hessian;
	VectorXd data_mesh2surface_grad;
	data_mesh2surface_grad.resize(num);
	data_mesh2surface_grad.setZero();
	data_mesh2surface_hessian.resize(num, num);
	data_mesh2surface_hessian.setZero();

	for (int k = 0; k < all_data_forward_projection.size(); k++)
	{
		
		 //if(SplineSurface::is_or[k]==1) 
		 {

			gp_Pnt pointToProject(all_data_forward_projection[k].input_point[0].value().value(), all_data_forward_projection[k].input_point[1].value().value(), all_data_forward_projection[k].input_point[2].value().value());
			GeomAPI_ProjectPointOnSurf projection(pointToProject, bezierSurface, tola);
			
			projection.Perform(pointToProject);
			
			Standard_Real u, v;

			projection.LowerDistanceParameters(u, v);
			
			all_data_forward_projection[k].u_cord = u;
			all_data_forward_projection[k].v_cord = v;

			Vec12 X;
			X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
			for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
			{
				X(3 * r).value() = spline(3 * r);
				X(3 * r + 1).value() = spline(3 * r + 1);
				X(3 * r + 2).value() = spline(3 * r + 2);
			}

			// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
				X(id).derivatives().setZero();
				X(id).derivatives()(id) = 1;
				X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
			}

			//  3.初始化海塞矩阵;set the hessian matrix to zero
			for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++) {
				for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
				{
					X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
				}
			}

			std::vector<vector<Vec3>> P;
			P.resize((U_Degree + 1));
			for (int i = 0; i < P.size(); i++)
			{
				P[i].resize((V_Degree + 1));
			}

			for (int i = 0; i < P.size(); i++)
			{
				P[i][0][0] = X[3 * i];
				P[i][0][1] = X[3 * i + 1];
				P[i][0][2] = X[3 * i + 2];
				P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
				P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
				P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
			}

			Vec3 bezier_surface_value;
			bezier_surface_value.setZero();
			for (int i = 0; i <= U_Degree; i++)
			{
				for (int j = 0; j <= V_Degree; j++)
				{

					Vec3 pij;
					pij = P[i][j];
					bezier_surface_value[0] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[0];
					bezier_surface_value[1] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[1];
					bezier_surface_value[2] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[2];
				}
			}
			
			scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * (bezier_surface_value - all_data_forward_projection[k].input_point).transpose();
			iosout << bezier_surface_value[0].value().value() << " " << bezier_surface_value[1].value().value() << " " << bezier_surface_value[2].value().value() << endl;
			iosout << all_data_forward_projection[k].input_point[0].value().value() << " " << all_data_forward_projection[k].input_point[1].value().value() << " " << all_data_forward_projection[k].input_point[2].value().value() << endl; //<< " " << all_data_forward_projection[k].input_point.transpose().value().value() << endl;
			/*double a1 = bezier_surface_value[0].value().value();
			double b1 = all_data_forward_projection[k].input_point[0].value().value();
			double a2 = bezier_surface_value[1].value().value();
			double b2 = all_data_forward_projection[k].input_point[1].value().value();
			double a3 = bezier_surface_value[2].value().value();
			double b3 = all_data_forward_projection[k].input_point[2].value().value();
			if (calc_iter_num == 0) {
				Eigen::Vector3d surface_normal = calc_bezier_surface_normal(3, 1, u, v, control_points);
				Eigen::Vector3d line1;

				line1[0] = a1 - b1;
				line1[1] = a2 - b2;
				line1[2] = a3 - b3;
				auto line2 = line1;
				line1.normalize();
				surface_normal.normalize();



				double theta = abs(line1.dot(surface_normal));
				if (theta > 0.99)
				{
					if (line2.norm() > all_data_forward_projection[0].maxdis) all_data_forward_projection[0].maxdis = line2.norm();
					if (line2.norm() < all_data_forward_projection[0].mindis) all_data_forward_projection[0].mindis = line2.norm();
					all_data_forward_projection[0].avrdis += line2.norm();
					all_data_forward_projection[k].is_orth = 1;
					all_data_forward_projection[k].initdis = line2.norm();

				}
				else {
					all_data_forward_projection[k].is_orth = 0;
				}
			}
			if (calc_iter_num == 1)
			{
				if (all_data_forward_projection[k].is_orth == 1)
				{
					if (all_data_forward_projection[k].initdis > (all_data_forward_projection[0].avrdis * 2 / 3))
					{
						all_data_forward_projection[k].is_orth = 0;
					}
				}
			}*/
			{
				
				//cout << "initdis=" << all_data_forward_projection[k].initdis << endl << "all_data_forward_projection[0].avrdis / 3.0=" << all_data_forward_projection[0].avrdis / 3.0 << endl;
				ornum++;
				data_mesh2surface += e;
				data_mesh2surface_grad += e.value().derivatives();
				Eigen::MatrixXd B;
				B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));
				
				for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
				{
					B.row(r) = e.derivatives()(r).derivatives().transpose();
				}
				data_mesh2surface_hessian += B;
			}
			
		}
	}

	scalar_t data_surface2mesh = 0;
	MatrixXd data_surface2mesh_hessian;
	VectorXd data_surface2mesh_grad;
	data_surface2mesh_grad.resize(num);
	data_surface2mesh_grad.setZero();
	data_surface2mesh_hessian.resize(num, num);
	data_surface2mesh_hessian.setZero();
	for (int k = 0; k < bi_vec.size(); k++)
	{
		Vec12 X;
		X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
		{
			X(3 * r).value() = spline(3 * r);
			X(3 * r + 1).value() = spline(3 * r + 1);
			X(3 * r + 2).value() = spline(3 * r + 2);
		}

		// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
		for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
		{
			X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
			X(id).derivatives().setZero();
			X(id).derivatives()(id) = 1;
			X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
		}

		//  3.初始化海塞矩阵;set the hessian matrix to zero
		for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++) {
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
			}
		}

		std::vector<vector<Vec3>> P;
		P.resize((U_Degree + 1));
		for (int i = 0; i < P.size(); i++)
		{
			P[i].resize((V_Degree + 1));
		}

		for (int i = 0; i < P.size(); i++)
		{
			P[i][0][0] = X[3 * i];
			P[i][0][1] = X[3 * i + 1];
			P[i][0][2] = X[3 * i + 2];
			P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
			P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
			P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
		}

		Vec3 bezier_surface_value;
		bezier_surface_value.setZero();
		for (int i = 0; i <= U_Degree; i++)
		{
			for (int j = 0; j <= V_Degree; j++)
			{
				Vec3 pij;
				pij = P[i][j];
				bezier_surface_value[0] += combination(U_Degree, i) * pow((1 - bi_vec[k].u_cord), U_Degree - i) * pow(bi_vec[k].u_cord, i) * combination(V_Degree, j) * pow((1 - bi_vec[k].v_cord), V_Degree - j) * pow(bi_vec[k].v_cord , j) * pij[0];
				bezier_surface_value[1] += combination(U_Degree, i) * pow((1 - bi_vec[k].u_cord), U_Degree - i) * pow(bi_vec[k].u_cord, i) * combination(V_Degree, j) * pow((1 - bi_vec[k].v_cord), V_Degree - j) * pow(bi_vec[k].v_cord , j) * pij[1];
				bezier_surface_value[2] += combination(U_Degree, i) * pow((1 - bi_vec[k].u_cord), U_Degree - i) * pow(bi_vec[k].u_cord, i) * combination(V_Degree, j) * pow((1 - bi_vec[k].v_cord), V_Degree - j) * pow(bi_vec[k].v_cord , j) * pij[2];
			}
		}

		bi_vec[k].bezier_point[0] = bezier_surface_value.x().value().value();
		bi_vec[k].bezier_point[1] = bezier_surface_value.y().value().value();
		bi_vec[k].bezier_point[2] = bezier_surface_value.z().value().value();

		auto result_closest = ordinary_tree.closest_point_and_primitive(Point(bi_vec[k].bezier_point[0], bi_vec[k].bezier_point[1], bi_vec[k].bezier_point[2]));
		bi_vec[k].mesh_point[0] = result_closest.first.x();
		bi_vec[k].mesh_point[1] = result_closest.first.y();
		bi_vec[k].mesh_point[2] = result_closest.first.z();

		scalar_t e = (bezier_surface_value - bi_vec[k].mesh_point) * (bezier_surface_value - bi_vec[k].mesh_point).transpose();
		data_surface2mesh += e;
		data_surface2mesh_grad += e.value().derivatives();
		Eigen::MatrixXd B;
		B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));

		for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
		{
			B.row(r) = e.derivatives()(r).derivatives().transpose();
		}
		data_surface2mesh_hessian += B;
	}

	grad = data_mesh2surface_weight * data_mesh2surface_grad + data_surface2mesh_weight * data_surface2mesh_grad;
	hessian = data_mesh2surface_weight * data_mesh2surface_hessian + data_surface2mesh_weight * data_surface2mesh_hessian;
	//ioss.close();
	std::cout << "ornum:" << ornum << endl;
	smooth_value = 0;
}


pair<double, double> SplineSurface::calc_data_term_energy(VectorXd& spline)
{

	double data_value = 0.0;
	std::vector<vector<Vec3>> P;
	P.resize((U_Degree + 1));
	for (int i = 0; i < P.size(); i++)
	{
		P[i].resize((V_Degree + 1));
	}

	for (int i = 0; i < P.size(); i++)
	{
		P[i][0][0] = spline[3 * i];
		P[i][0][1] = spline[3 * i + 1];
		P[i][0][2] = spline[3 * i + 2];
		P[i][1][0] = spline[3 * (U_Degree + 1) + 3 * i];
		P[i][1][1] = spline[3 * (U_Degree + 1) + 3 * i + 1];
		P[i][1][2] = spline[3 * (U_Degree + 1) + 3 * i + 2];
	}

	for (int k = 0; k < all_data_forward_projection.size(); k++)
	{
		Vec3 bezier_surface_value;
		bezier_surface_value.setZero();
		for (int i = 0; i <= U_Degree; i++)
		{
			for (int j = 0; j <= V_Degree; j++)
			{

				Vec3 pij;
				pij = P[i][j];

				bezier_surface_value[0] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[0];
				bezier_surface_value[1] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[1];
				bezier_surface_value[2] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[2];
			}
		}

		//scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * all_data_forward_projection[k].normal_on_output_surface.transpose();
		scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * (bezier_surface_value - all_data_forward_projection[k].input_point).transpose();
		//e = e * e;
		data_value += e.value().value();
		/*
		scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * (bezier_surface_value - all_data_forward_projection[k].input_point).transpose();
		data_value += e.value().value();*/
		/*
		if (bi_vec[k].is_intersection)
		{
			scalar_t e = ((bezier_surface_value - bi_vec[k].mesh_point) * (bezier_surface_value - bi_vec[k].mesh_point).transpose());
			//e = e * e;
			data_value += e.value().value();
		}*/
	}

	
	double smooth_value = 0;
	/*
	for (int ku = 0; ku < u_knot.size(); ku++)
	{
		for (int kv = 0; kv < v_knot.size(); kv++)
		{
			Vec3 s_value;
			s_value.setZero();
			pair<Eigen::MatrixXd, Eigen::MatrixXd> smooth_coff = calc_uv_derivative(u_knot[ku], v_knot[kv]);
			for (int dim_row = 0; dim_row < U_Degree + 1; dim_row++)
			{
				for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
				{
					Vec3 pij; pij = P[dim_row][dim_col];
					s_value[0] += smooth_coff.first.coeff(dim_row, dim_col) * pij[0] + smooth_coff.second.coeff(dim_row, dim_col) * pij[0];
					s_value[1] += smooth_coff.first.coeff(dim_row, dim_col) * pij[1] + smooth_coff.second.coeff(dim_row, dim_col) * pij[1];
					s_value[2] += smooth_coff.first.coeff(dim_row, dim_col) * pij[2] + smooth_coff.second.coeff(dim_row, dim_col) * pij[2];
				}
			}
			scalar_t s = s_value * s_value.transpose();
			smooth_value += s.value().value();
		}
	}*/


	/*
	for (int i = 0; i < P.size() - 2; i++)
	{
		for (int j = 0; j < P[i].size(); j++)
		{
			Vec3 pi2_j = P[i + 2][j];
			Vec3 pi1_j = P[i + 1][j];
			Vec3 pi0_j = P[i][j];
			Vec3 sub_ennery;
			sub_ennery(0) = pi2_j(0) - 2 * pi1_j(0) + pi0_j(0);
			sub_ennery(1) = pi2_j(1) - 2 * pi1_j(1) + pi0_j(1);
			sub_ennery(2) = pi2_j(2) - 2 * pi1_j(2) + pi0_j(2);

			scalar_t s = smooth_weight * sub_ennery * sub_ennery.transpose();
			smooth_value += s.value().value();
		}
	}


	for (int i = 0; i < P.size(); i++)
	{
		Vec3 pi1_j = P[i][0];
		Vec3 pi0_j = P[i][1];
		Vec3 sub_ennery;
		sub_ennery(0) = pi1_j(0) - pi0_j(0)-1;
		sub_ennery(1) = pi1_j(1) - pi0_j(1)-1;
		sub_ennery(2) = pi1_j(2) - pi0_j(2)-1;

		scalar_t s = smooth_weight * sub_ennery * sub_ennery.transpose();
		smooth_value += s.value().value();
	}*/

	pair<double, double> pair_value = make_pair(data_value, smooth_value);
	return pair_value;
}

pair<double, double> SplineSurface::gcl_calc_data_term_energy(const Data& spline)
{
	VectorXd spline_v;
	double d1, d2, d3;
	d1 = 1000;
	d2 = 0;
	d3 = 0;
	spline_v.resize(spline.rows() * 3);
	for (int i = 0; i < spline.rows(); i++)
	{
		spline_v[3 * i] = spline.coeff(i, 0);
		spline_v[3 * i + 1] = spline.coeff(i, 1);
		spline_v[3 * i + 2] = spline.coeff(i, 2);
	}

	double data_value = 0.0;
	std::vector<vector<Vec3>> P;
	P.resize((U_Degree + 1));
	for (int i = 0; i < P.size(); i++)
	{
		P[i].resize((V_Degree + 1));
	}

	for (int i = 0; i < P.size(); i++)
	{
		P[i][0][0] = spline_v[3 * i];
		P[i][0][1] = spline_v[3 * i + 1];
		P[i][0][2] = spline_v[3 * i + 2];
		P[i][1][0] = spline_v[3 * (U_Degree + 1) + 3 * i];
		P[i][1][1] = spline_v[3 * (U_Degree + 1) + 3 * i + 1];
		P[i][1][2] = spline_v[3 * (U_Degree + 1) + 3 * i + 2];
	}


	double data_mesh2surface_value = 0.0;
	for (int k = 0; k < all_data_forward_projection.size(); k++)
	{
		//if (SplineSurface::is_or[k] == 1)
		{
			Vec3 bezier_surface_value;
			bezier_surface_value.setZero();
			for (int i = 0; i <= U_Degree; i++)
			{
				for (int j = 0; j <= V_Degree; j++)
				{
					Vec3 pij;
					pij = P[i][j];
					bezier_surface_value[0] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[0];
					bezier_surface_value[1] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[1];
					bezier_surface_value[2] += combination(U_Degree, i) * pow((1 - all_data_forward_projection[k].u_cord), U_Degree - i) * pow(all_data_forward_projection[k].u_cord, i) * combination(V_Degree, j) * pow((1 - all_data_forward_projection[k].v_cord), V_Degree - j) * pow(all_data_forward_projection[k].v_cord, j) * pij[2];
				}
			}

			//scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * all_data_forward_projection[k].normal_on_output_surface.transpose();
			scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * (bezier_surface_value - all_data_forward_projection[k].input_point).transpose();
			{
				data_mesh2surface_value += e.value().value();
				if (e.value().value() < d1) d1 = e.value().value();
				if (e.value().value() > d2) d2 = e.value().value();
				//d3 += e.value().value() / all_data_forward_projection.size();
			}
			/*
			scalar_t e = (bezier_surface_value - all_data_forward_projection[k].input_point) * (bezier_surface_value - all_data_forward_projection[k].input_point).transpose();
			data_value += e.value().value();*/
			/*
			if (bi_vec[k].is_intersection)
			{
				scalar_t e = ((bezier_surface_value - bi_vec[k].mesh_point) * (bezier_surface_value - bi_vec[k].mesh_point).transpose());
				//e = e * e;
				data_value += e.value().value();
			}*/
		}
	}

	double data_surface2mesh = 0.0;
	for (int k = 0; k < bi_vec.size(); k++)
	{
		auto result_closest = ordinary_tree.closest_point_and_primitive(Point(bi_vec[k].bezier_point[0], bi_vec[k].bezier_point[1], bi_vec[k].bezier_point[2]));
		Vector3d close_point;
		close_point[0] = result_closest.first.x();
		close_point[1] = result_closest.first.y();
		close_point[2] = result_closest.first.z();

		double e = (bi_vec[k].bezier_point - close_point).transpose() *(bi_vec[k].bezier_point - close_point);
		data_surface2mesh += e;
	}

	data_value = data_mesh2surface_weight * data_mesh2surface_value + data_surface2mesh_weight * data_surface2mesh;


	std::cout << " 最小距离" << d1 << endl;
	std::cout << " 最大距离" << d2 << endl;
	//cout << " 平均距离" << d3 << endl;
	double smooth_value = 0.0;
	if (smooth_model == 0)
	{
		scalar_t s = 0;
		for (int dim_row = 0; dim_row < U_Degree; dim_row++)
		{
			Vec3 pi0; pi0 = P[dim_row][0];
			Vec3 pi1; pi1 = P[dim_row][1];

			Vec3 pii0; pii0 = P[dim_row + 1][0];
			Vec3 pii1; pii1 = P[dim_row + 1][1];

			Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
			s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
		}
		smooth_value = s.value().value();
	}
	else if (smooth_model == 1)
	{
		scalar_t s = 0;
		for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
		{
			for (int dim_row = 0; dim_row < U_Degree - 1; dim_row++)
			{
				Vec3 pij; pij = P[dim_row][dim_col];
				Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
				Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
				Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
				s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
			}
		}
		smooth_value = s.value().value();
	}
	else if (smooth_model == 2)
	{
		scalar_t s = 0;
		for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
		{
			for (int dim_row = 0; dim_row < U_Degree - 1; dim_row++)
			{
				Vec3 pij; pij = P[dim_row][dim_col];
				Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
				Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
				Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
				s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
			}
		}

		scalar_t s1 = 0;
		for (int dim_row = 0; dim_row < U_Degree; dim_row++)
		{
			Vec3 pi0; pi0 = P[dim_row][0];
			Vec3 pi1; pi1 = P[dim_row][1];

			Vec3 pii0; pii0 = P[dim_row + 1][0];
			Vec3 pii1; pii1 = P[dim_row + 1][1];

			Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
			s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
		}
		scalar_t s_total = u_smooth_weight * s + v_smooth_weight * s1;
		//cout<<"u光滑项："<<
		smooth_value = s_total.value().value();
	}

	/*
	for (int ku = 0; ku < u_knot.size(); ku++)
	{
		for (int kv = 0; kv < v_knot.size(); kv++)
		{
			Vec3 s_value;
			s_value.setZero();
			pair<Eigen::MatrixXd, Eigen::MatrixXd> smooth_coff = calc_uv_derivative(u_knot[ku], v_knot[kv]);
			for (int dim_row = 0; dim_row < U_Degree + 1; dim_row++)
			{
				for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
				{
					Vec3 pij; pij = P[dim_row][dim_col];
					s_value[0] += smooth_coff.first.coeff(dim_row, dim_col) * pij[0] + smooth_coff.second.coeff(dim_row, dim_col) * pij[0];
					s_value[1] += smooth_coff.first.coeff(dim_row, dim_col) * pij[1] + smooth_coff.second.coeff(dim_row, dim_col) * pij[1];
					s_value[2] += smooth_coff.first.coeff(dim_row, dim_col) * pij[2] + smooth_coff.second.coeff(dim_row, dim_col) * pij[2];
				}
			}
			scalar_t s = s_value * s_value.transpose();
			smooth_value += s.value().value();
		}
	}*/


	/*
	for (int i = 0; i < P.size() - 2; i++)
	{
		for (int j = 0; j < P[i].size(); j++)
		{
			Vec3 pi2_j = P[i + 2][j];
			Vec3 pi1_j = P[i + 1][j];
			Vec3 pi0_j = P[i][j];
			Vec3 sub_ennery;
			sub_ennery(0) = pi2_j(0) - 2 * pi1_j(0) + pi0_j(0);
			sub_ennery(1) = pi2_j(1) - 2 * pi1_j(1) + pi0_j(1);
			sub_ennery(2) = pi2_j(2) - 2 * pi1_j(2) + pi0_j(2);

			scalar_t s = smooth_weight * sub_ennery * sub_ennery.transpose();
			smooth_value += s.value().value();
		}
	}


	for (int i = 0; i < P.size(); i++)
	{
		Vec3 pi1_j = P[i][0];
		Vec3 pi0_j = P[i][1];
		Vec3 sub_ennery;
		sub_ennery(0) = pi1_j(0) - pi0_j(0)-1;
		sub_ennery(1) = pi1_j(1) - pi0_j(1)-1;
		sub_ennery(2) = pi1_j(2) - pi0_j(2)-1;

		scalar_t s = smooth_weight * sub_ennery * sub_ennery.transpose();
		smooth_value += s.value().value();
	}*/

	pair<double, double> pair_value = make_pair(data_value, smooth_value);
	return pair_value;
}



double SplineSurface::line_search(Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction, Eigen::MatrixXd& spline1, double c_parm, double beta, double data_term, double smooth_term)
{
	VectorXd spline;
	spline.resize(spline1.rows() * 3);
	for (int i = 0; i < spline1.rows(); i++)
	{
		spline[3 * i] = spline1.coeff(i, 0);
		spline[3 * i + 1] = spline1.coeff(i, 1);
		spline[3 * i + 2] = spline1.coeff(i, 2);
	}

	double step_size = 1.0;// Initial step size
	int iter_num = 0;
	while (iter_num < 10)
	{
		Eigen::VectorXd new_spline = spline + step_size * decrease_direction;
		auto pair_value = calc_data_term_energy(new_spline);
		std::cout << endl;
		std::cout << "setp: " << step_size << endl;
		std::cout << "new_total: " << pair_value.first + smooth_weight * pair_value.second << ", data_value_new: " << pair_value.first << ", smooth_value_new: " << smooth_weight * pair_value.second << endl;
		std::cout << "old_total: " << data_term + smooth_weight * smooth_term << ", data_value: " << data_term << ", smooth_value: " << smooth_weight * smooth_term << endl;

		if (pair_value.first + smooth_weight * pair_value.second > data_term + smooth_weight * smooth_term + c_parm * step_size * decrease_direction.transpose() * grad)
		{
			step_size *= beta;
		}
		else
		{
			return step_size;
			break;
		}

		iter_num++;
		if (iter_num == 10)
		{
			step_size = 0;
			return step_size;
		}
	}
}

void SplineSurface::optimization(Mesh& mesh, Eigen::MatrixXd& spline)//vector<vector<Eigen::Vector3d>>& control_points)
{
	init_ready_data(mesh);
	Eigen::MatrixXd new_spline;
	new_spline.resize((U_Degree + 1) * (V_Degree + 1), 3); new_spline.setZero();

	/*
	for (int i = 0; i < control_points.size(); i++)
	{
		spline[3 * i] = control_points[i][0][0];
		spline[3 * i + 1] = control_points[i][0][1];
		spline[3 * i + 2] = control_points[i][0][2];
		spline[12 + 3 * i] = control_points[i][1][0];
		spline[12 + 3 * i + 1] = control_points[i][1][1];
		spline[12 + 3 * i + 2] = control_points[i][1][2];
	}*/

	double stop_threshold = 1e-6;
	uv_sample(100);
	int iter_num = 0;


	Tree tree;
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();


	Tree_ tree_;
	tree_.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	tree_.build();
	tree_.accelerate_distance_queries();

	gcl_bezier_visualization(100, spline, -1, -1);
	while (true)
	{
		Eigen::VectorXd grad, decrease_direction;
		Eigen::MatrixXd hassien;
		Eigen::VectorXd smooth_grad;
		Eigen::MatrixXd smooth_hassien;
		Eigen::VectorXd total_grad;
		Eigen::MatrixXd total_hassien;
		scalar_t y_value; scalar_t smooth_value;
		//calc_data_term_energy_derivative(spline, y_value, smooth_value, grad, hassien, decrease_direction,tree,tree_);
		gcl_calc_data_smooth_term_derivative(spline, y_value, grad, hassien, decrease_direction, smooth_value, smooth_grad, smooth_hassien, total_grad, total_hassien);
		double data_term = y_value.value().value();
		double smooth_term = smooth_value.value().value();
		double step = line_search(total_grad, decrease_direction, spline, 1e-4, 0.5, data_term, smooth_term);

		Eigen::MatrixXd direction;
		direction.resize((U_Degree + 1) * (V_Degree + 1), 3);
		for (int i = 0; i < direction.rows(); i++)
		{
			direction(i, 0) = decrease_direction[3 * i];
			direction(i, 1) = decrease_direction[3 * i + 1];
			direction(i, 2) = decrease_direction[3 * i + 2];
		}

		new_spline = spline + step * direction;
		double termination_condition = (new_spline - spline).squaredNorm();
		spline = new_spline;

		gcl_bezier_visualization(100, spline, iter_num, 0);
		iter_num++;
		if (termination_condition < stop_threshold)
		{
			break;
		}
	}
}

std::vector<double> SplineSurface::generate_equally_spaced_vector(int num_elements)//yes
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

Eigen::Vector3d SplineSurface::calc_bezier_surface_value(int u_degree, int v_degree, double u_cord, double v_cord, vector<vector<Eigen::Vector3d>>& control_points)//控制点数据类型 Matrix<Vector3d>
{
	Eigen::Vector3d bezier_surface_value(0, 0, 0);
	for (int i = 0; i <= u_degree; i++)
	{
		for (int j = 0; j <= v_degree; j++)
		{
			Eigen::Vector3d pij;
			pij = control_points[i][j];
			bezier_surface_value += combination(u_degree, i) * pow((1 - u_cord), u_degree - i) * pow(u_cord, i) * combination(v_degree, j) * pow((1 - v_cord), v_degree - j) * pow(v_cord, j) * pij;
		}
	}
	return bezier_surface_value;
}

void SplineSurface::bezier_visualization(int row_num, vector<vector<Eigen::Vector3d>>& control_points, int iter_num, int bezier_num)
{
	//1.基于采样
	vector<double> u_vec = generate_equally_spaced_vector(row_num);
	vector<double> v_vec = generate_equally_spaced_vector(row_num);

	//2.Bezier_info构造
	Matrix<SplineSurface::Bezier_info, Eigen::Dynamic, Eigen::Dynamic> bezier_info_mat;
	bezier_info_mat.resize(row_num, row_num);
	for (int i = 0; i < u_vec.size(); i++)
	{
		for (int j = 0; j < v_vec.size(); j++)
		{
			Bezier_info bi;
			bi.u_cord = u_vec[i];
			bi.v_cord = v_vec[j];
			bi.u_idx = i;
			bi.v_idx = j;
			bi.bezier_point = calc_bezier_surface_value(U_Degree, V_Degree, bi.u_cord, bi.v_cord, control_points);
			bezier_info_mat(i, j) = bi;
		}
	}

	//3.
	Mesh outmesh;
	Matrix<int, Eigen::Dynamic, Eigen::Dynamic> v_idx_mat;
	v_idx_mat.resize(row_num, row_num);
	int init_v_idx = 0;
	for (int i = 0; i < bezier_info_mat.rows(); i++)
	{
		for (int j = 0; j < bezier_info_mat.cols(); j++)
		{
			Mesh::VertexHandle v0 = outmesh.add_vertex(Mesh::Point(bezier_info_mat(i, j).bezier_point.data()));
			v_idx_mat(i, j) = init_v_idx;
			init_v_idx++;
		}
	}

	for (int i = 0; i < row_num - 1; ++i)
	{
		for (int j = 0; j < row_num - 1; ++j)
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
	MeshTools::WriteMesh(mesh, "./inputmesh.obj");
	MeshTools::WriteMesh(outmesh, "./outputmesh_" + to_string(iter_num) + "_" + to_string(bezier_num) + ".obj");
}


double SplineSurface::spline_ray_intersection(double x_cord, double y_cord, Tree_& tree_)
{


	Point ray_source(x_cord, y_cord, 0); // 射线源点坐标
	Point ray_direction(x_cord, y_cord, 1); // 射线的方向向量
	Ray ray_query(ray_source, ray_direction);

	/*
	if (tree.do_intersect(ray_query))
		std::cout << "intersection(s)" << std::endl;
	else
		std::cout << "no intersection" << std::endl;
		// computes #intersections with segment query
		//std::cout << tree.number_of_intersected_primitives(ray_query) << " intersection(s)" << std::endl;*/


	Ray_intersection intersection = tree_.first_intersection(ray_query);
	if (intersection)
	{
		if (boost::get<Point>(&(intersection->first)))
		{
			const Point* p = boost::get<Point>(&(intersection->first));
			return p->z();
		}
	}
	else
	{
		return -1003;
	}
}

void SplineSurface::test()
{
	Bezier_info bi;
	std::cout << bi.is_intersection << endl;
	bi.is_intersection = 0;
	std::cout << bi.is_intersection << endl;
}

//******************************  IPC_OPTIMIZATION  ******************************
double SplineSurface::ipc_line_search(Mesh& mesh, vector<vector<Eigen::MatrixXd>>& all_bezier, Eigen::VectorXd& grad, Eigen::VectorXd& decrease_direction, Eigen::VectorXd& barrier_grad, scalar_t& barrier_function, VectorXd& spline,
	double c_parm, double beta, double data_term, Tree_& tree_,
	vector<Eigen::RowVector3d>& P_input_set, vector<vector<Eigen::RowVector3d>>& Edge_input_set, vector<vector<Eigen::RowVector3d>>& face_input_set)
{
	ofstream ioss("line_search_" + to_string(line_search_num) + ".txt");
	double step_size = 1.0;// Initial step size
	int iter_num = 0;
	while (iter_num < 20)
	{
		bool is_safe = true;
		Eigen::VectorXd new_spline = spline + step_size * decrease_direction;

		double ccd_dist = step_CCD_check(spline, new_spline, all_bezier);
		if (ccd_dist < d0)
		{
			is_safe = false;
		}

		pair<double, double> values = calc_data_term_energy(new_spline);
		double data_term_new = values.first;
		scalar_t y_new_barrier;
		//Gradient::calc_barrier_function(mesh, new_spline, y_new_barrier, all_bezier, P_input_set, Edge_input_set, face_input_set);
		double barrier_value_new = y_new_barrier.value().value();
		double wlf_new = data_term_new + lamda * barrier_value_new;
		double wlf_old = data_term + lamda * barrier_function.value().value() + c_parm * step_size * decrease_direction.transpose() * (grad + lamda * barrier_grad);

		ioss << "step: " << step_size << endl;
		ioss << "new_total: " << "new_data_value: " << data_term_new << " ,new_barrier_value: " << barrier_value_new << endl;
		ioss << "old_total: " << "old_data_value: " << data_term << " ,old_barrier_value: " << barrier_function.value().value() << endl;
		ioss << "wlf_new " << wlf_new << " ,wlf_old: " << wlf_old << endl;
		ioss << endl;

		std::cout << endl;
		std::cout << "step: " << step_size << endl;
		std::cout << "new_total: " << "new_data_value: " << data_term_new << " ,new_barrier_value: " << barrier_value_new << endl;
		std::cout << "old_total: " << "old_data_value: " << data_term << " ,old_barrier_value: " << barrier_function.value().value() << endl;
		std::cout << "wlf_new " << wlf_new << " ,wlf_old: " << wlf_old << endl;

		if ((wlf_new <= wlf_old) && is_safe)
		{
			return step_size;
		}
		else
		{
			step_size *= beta;
		}

		iter_num++;
		if (iter_num == 20)
		{
			step_size = 0;
			return step_size;
		}
	}
	line_search_num++;
	ioss.close();

}

void SplineSurface::ipc_total_derivative(Mesh& mesh, VectorXd& spline, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian, scalar_t& data_value, Eigen::VectorXd& grad, Eigen::MatrixXd& hessian, scalar_t& barrier_value, Eigen::VectorXd& barrier_grad, Eigen::MatrixXd& barrier_hessian, Eigen::VectorXd& decrease_direction, Tree& tree, Tree_& tree_,
	vector<vector<Eigen::MatrixXd>> all_bezier,
	vector<Eigen::RowVector3d>& P_input_set, vector<vector<Eigen::RowVector3d>>& Edge_input_set, vector<vector<Eigen::RowVector3d>>& face_input_set)
{
	Eigen::VectorXd decrease_direction_data_term; scalar_t smooth_value;
	calc_data_term_energy_derivative(spline, data_value, smooth_value, grad, hessian, decrease_direction, tree, tree_);
	//Gradient::ipc_barrier_gradient(mesh, spline, barrier_value, barrier_grad, barrier_hessian, all_bezier,
		//P_input_set, Edge_input_set, face_input_set);

	total_grad = grad + lamda * barrier_grad;
	total_hessian = hessian + lamda * barrier_hessian;
	//cout << total_grad << endl;
	//cout << total_hessian << endl;
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


void SplineSurface::ipc_optimization(Mesh& mesh, vector<vector<Eigen::Vector3d>>& control_points)
{
	init_ready_data(mesh);
	Eigen::VectorXd spline, new_spline;
	spline.resize((U_Degree + 1) * (V_Degree + 1) * 3); spline.setZero(); new_spline.resize((U_Degree + 1) * (V_Degree + 1) * 3); new_spline.setZero();
	for (int i = 0; i < control_points.size(); i++)
	{
		spline[3 * i] = control_points[i][0][0];
		spline[3 * i + 1] = control_points[i][0][1];
		spline[3 * i + 2] = control_points[i][0][2];
		spline[(U_Degree + 1) * 3 + 3 * i] = control_points[i][1][0];
		spline[(U_Degree + 1) * 3 + 3 * i + 1] = control_points[i][1][1];
		spline[(U_Degree + 1) * 3 + 3 * i + 2] = control_points[i][1][2];
	}

	double stop_threshold = 1e-10;
	uv_sample(100);
	int iter_num = 0;
	Tree tree;

	Tree_ tree_;
	tree_.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	tree_.build();
	tree_.accelerate_distance_queries();

	bezier_visualization(100, control_points, -1, -1);
	vector<Eigen::RowVector3d> P_input_set; vector<vector<Eigen::RowVector3d>> Edge_input_set; vector<vector<Eigen::RowVector3d>> face_input_set;
	//BezierSubvision::get_input_point_edge_face(mesh, P_input_set, Edge_input_set, face_input_set);

	while (true)
	{
		vector<vector<Eigen::MatrixXd>> all_bezier;
		bezier_surface_subvision(mesh, spline, all_bezier);

		Eigen::VectorXd data_grad, decrease_direction;
		Eigen::MatrixXd hassien;
		scalar_t data_value; scalar_t barrier_value;
		Eigen::VectorXd total_grad; Eigen::MatrixXd total_hessian; Eigen::VectorXd barrier_grad; Eigen::MatrixXd barrier_hessian;
		ipc_total_derivative(mesh, spline, total_grad, total_hessian, data_value, data_grad, hassien, barrier_value, barrier_grad, barrier_hessian, decrease_direction, tree, tree_,
			all_bezier, P_input_set, Edge_input_set, face_input_set);


		double step = ipc_line_search(mesh, all_bezier, data_grad, decrease_direction, barrier_grad, barrier_value, spline,
			1e-4, 0.6, data_value.value().value(), tree_,
			P_input_set, Edge_input_set, face_input_set);

		new_spline = spline + step * decrease_direction;
		ofstream ios("spline_" + to_string(iter_num) + "_.txt");
		ios << new_spline << endl;
		ios.close();

		Eigen::VectorXd data_grad_coff, decrease_direction_coff;
		Eigen::MatrixXd hassien_coff;
		scalar_t data_value_coff; scalar_t barrier_value_coff;
		Eigen::VectorXd total_grad_coff; Eigen::MatrixXd total_hessian_coff; Eigen::VectorXd barrier_grad_coff; Eigen::MatrixXd barrier_hessian_coff;
		//lamda=update_barrier_coffient(mesh, new_spline, total_grad_coff, total_hessian_coff, data_value_coff, data_grad_coff, hassien_coff, barrier_value_coff, barrier_grad_coff, barrier_hessian_coff, decrease_direction_coff, tree,
		//	all_bezier, P_input_set, Edge_input_set, face_input_set);
		double termination_condition = (new_spline - spline).squaredNorm();
		spline = new_spline;

		for (int i = 0; i < control_points.size(); i++)
		{
			control_points[i][0][0] = spline[3 * i];
			control_points[i][0][1] = spline[3 * i + 1];
			control_points[i][0][2] = spline[3 * i + 2];
			control_points[i][1][0] = spline[(U_Degree + 1) * 3 + 3 * i];
			control_points[i][1][1] = spline[(U_Degree + 1) * 3 + 3 * i + 1];
			control_points[i][1][2] = spline[(U_Degree + 1) * 3 + 3 * i + 2];
		}

		bezier_visualization(100, control_points, iter_num, 0);
		iter_num++;
		if (termination_condition < stop_threshold)
		{
			break;
		}
	}
}

double SplineSurface::step_CCD_check(VectorXd& spline, VectorXd& new_spline, vector<vector<Eigen::MatrixXd>>& all_bezier)
{
	vector<vector<vector<Point>>> two_bezier_points;
	two_bezier_points.resize(all_bezier.size());
	for (int i = 0; i < two_bezier_points.size(); i++)
	{
		two_bezier_points[i].resize(all_bezier.size());
	}

	for (int i = 0; i < all_bezier.size(); i++)
	{
		for (int j = 0; j < all_bezier[i].size(); j++)
		{
			std::vector<Point> P;
			P.resize(2 * (U_Degree + 1) * (V_Degree + 1));
			for (int dim1 = 0; dim1 < (U_Degree + 1) * (V_Degree + 1); dim1++)
			{
				Eigen::Vector3d pi, pj;
				pi.setZero(); pj.setZero();
				for (int dim2 = 0; dim2 < (U_Degree + 1) * (V_Degree + 1); dim2++)
				{
					pi[0] += all_bezier[i][j](dim1, dim2) * spline[3 * dim2];//基函数
					pi[1] += all_bezier[i][j](dim1, dim2) * spline[3 * dim2 + 1];
					pi[2] += all_bezier[i][j](dim1, dim2) * spline[3 * dim2 + 2];

					pj[0] += all_bezier[i][j](dim1, dim2) * new_spline[3 * dim2];//基函数
					pj[1] += all_bezier[i][j](dim1, dim2) * new_spline[3 * dim2 + 1];
					pj[2] += all_bezier[i][j](dim1, dim2) * new_spline[3 * dim2 + 2];

				}

				P[dim1] = { pi[0] ,pi[1], pi[2] };
				P[(U_Degree + 1) * (V_Degree + 1) + dim1] = { pj[0] ,pj[1], pj[2] };
			}
			two_bezier_points[i][j] = P;
		}
	}

	double all_cmp_dist = 1e10;
	for (int i = 0; i < two_bezier_points.size(); i++)
	{
		for (int j = 0; j < two_bezier_points[i].size(); j++)
		{
			double dist = calc_dist_convex_hull_and_input_mesh(two_bezier_points[i][j]);
			all_cmp_dist = std::min(all_cmp_dist, dist);
		}
	}

	return all_cmp_dist;
}

double SplineSurface::calc_dist_convex_hull_and_input_mesh(vector<Point>& points) //求双向的距离才行
{

	vector<Point> convex_hull_point = points;
	//1.line2line
	vector<Segment> convex_hull_line;
	for (int i = 0; i < convex_hull_point.size() - 1; i++)
	{
		for (int j = i + 1; j < convex_hull_point.size(); j++)
		{
			Segment s1(convex_hull_point[i], convex_hull_point[j]);
			convex_hull_line.push_back(s1);
		}
	}

	double sum_clog_line2line = 1e10;
	for (int i = 0; i < input_mesh_line_init.size(); i++)
	{
		for (int j = 0; j < convex_hull_line.size(); j++)
		{
			double dist_line2line = std::sqrtf(CGAL::squared_distance(input_mesh_line_init[i], convex_hull_line[j]));
			sum_clog_line2line = std::min(dist_line2line, sum_clog_line2line);
		}
	}

	vector<Triangle> convex_hull_triangle;
	std::vector<bool> mask(3, true);
	mask.resize(convex_hull_point.size(), false);

	do {
		std::vector<Point> combination;
		for (int i = 0; i < convex_hull_point.size(); ++i) {
			if (mask[i]) {
				combination.push_back(convex_hull_point[i]);
			}
		}
		Triangle tri(combination[0], combination[1], combination[2]);
		convex_hull_triangle.push_back(tri);
	} while (std::prev_permutation(mask.begin(), mask.end()));

	double sum_clog_p2tri = 1e10;
	for (int i = 0; i < input_mesh_point_init.size(); i++)
	{
		for (int j = 0; j < convex_hull_triangle.size(); j++)
		{
			double dist_p2tri = std::sqrtf(CGAL::squared_distance(input_mesh_point_init[i], convex_hull_triangle[j]));
			sum_clog_p2tri = std::min(sum_clog_p2tri, dist_p2tri);
		}
	}

	double sum_clog_tri2p = 1e10;
	for (int i = 0; i < convex_hull_point.size(); i++)
	{
		for (int j = 0; j < input_mesh_triangle_init.size(); j++)
		{
			double dist_tri2p = std::sqrtf(CGAL::squared_distance(convex_hull_point[i], input_mesh_triangle_init[j]));
			sum_clog_tri2p = std::min(sum_clog_tri2p, dist_tri2p);
		}
	}

	double value = min(sum_clog_line2line, sum_clog_tri2p);
	return min(value, sum_clog_p2tri);
}

void SplineSurface::bezier_surface_subvision(Mesh& mesh, VectorXd& spline, vector<vector<Eigen::MatrixXd>>& all_bezier)
{
	BezierSubvision bs(mesh);
	int sub_num = 2;
	bool is_sub = false;
	vector<Eigen::MatrixXd> result = bs.get_all_sub_bezier_curve_mat(sub_num, U_Degree);
	all_bezier.resize(pow(2, sub_num));
	for (int i = 0; i < all_bezier.size(); i++)
	{
		all_bezier[i].resize(pow(2, sub_num));
	}
	for (int i = 0; i < pow(2, sub_num); i++)
	{
		for (int j = 0; j < pow(2, sub_num); j++)
		{
			auto m_1 = result[j].block((U_Degree + 1) * i, 0, (U_Degree + 1), (U_Degree + 1) * (V_Degree + 1));
			auto m_2 = result[j + 1].block((U_Degree + 1) * i, 0, (U_Degree + 1), (U_Degree + 1) * (V_Degree + 1));
			MatrixXd BindRows(m_1.rows() + m_2.rows(), (U_Degree + 1) * (V_Degree + 1));
			BindRows << m_1, m_2;
			all_bezier[i][j] = BindRows;
		}
	}
}

Eigen::Vector3d SplineSurface::calc_bezier_surface_normal(int u_degree, int v_degree, double u_cord, double v_cord, vector<vector<Eigen::Vector3d>>& control_points)
{
	Eigen::Vector3d bezier_surface_u_derivative(0, 0, 0);
	Eigen::Vector3d bezier_surface_v_derivative(0, 0, 0);

	bezier_surface_u_derivative = (-3 * pow((1 - u_cord), 2)) * (1 - v_cord) * control_points[0][0] + (-3 * pow((1 - u_cord), 2)) * v_cord * control_points[0][1]
		+ 3 * (3 * u_cord * u_cord - 4 * u_cord + 1) * (1 - v_cord) * control_points[1][0] + 3 * (3 * u_cord * u_cord - 4 * u_cord + 1) * v_cord * control_points[1][1]
		+ 3 * (2 * u_cord - 3 * u_cord * u_cord) * (1 - v_cord) * control_points[2][0] + 3 * (2 * u_cord - 3 * u_cord * u_cord) * v_cord * control_points[2][1]
		+ (3 * u_cord * u_cord) * (1 - v_cord) * control_points[3][0] + (3 * u_cord * u_cord) * v_cord * control_points[3][1];

	vector<double> dv = { -1,1 };
	for (int i = 0; i <= u_degree; i++)
	{
		for (int j = 0; j <= v_degree; j++)
		{
			Eigen::Vector3d pij;
			pij = control_points[i][j];
			bezier_surface_v_derivative += combination(u_degree, i) * pow((1 - u_cord), u_degree - i) * pow(u_cord, i) * dv[j] * pij;
		}
	}

	Eigen::Vector3d bezier_surface_normal = bezier_surface_u_derivative.cross(bezier_surface_v_derivative);
	return bezier_surface_normal;
}

void SplineSurface::knot_vectors()
{

}

pair<Eigen::MatrixXd, Eigen::MatrixXd> SplineSurface::calc_uv_derivative(double u_, double v_)
{
	Eigen::MatrixXd Fuu; Fuu.setZero((U_Degree + 1), (V_Degree + 1));
	Eigen::MatrixXd Fuv; Fuv.setZero((U_Degree + 1), (V_Degree + 1));

	Fuu(0, 0) = 6 * (1 - u_) * (1 - v_); Fuu(0, 1) = 6 * (1 - u_) * v_;
	Fuu(1, 0) = 3 * (6 * u_ - 4) * (1 - v_); Fuu(1, 1) = 3 * (6 * u_ - 4) * v_;
	Fuu(2, 0) = 3 * (2 - 6 * u_) * (1 - v_); Fuu(2, 1) = 3 * (2 - 6 * u_) * v_;
	Fuu(3, 0) = 6 * u_ * (1 - v_); Fuu(3, 1) = 6 * u_ * v_;

	Fuv(0, 0) = 3 * (1 - u_) * (1 - u_); Fuv(0, 1) = -3 * (1 - u_) * (1 - u_);
	Fuv(1, 0) = -3 * (3 * u_ * u_ - 4 * u_ + 1); Fuv(1, 1) = 3 * (3 * u_ * u_ - 4 * u_ + 1);
	Fuv(2, 0) = -3 * (2 * u_ - 3 * u_ * u_); Fuv(2, 1) = 3 * (2 * u_ - 3 * u_ * u_);
	Fuv(3, 0) = -3 * u_ * u_; Fuv(3, 1) = 3 * u_ * u_;

	pair<Eigen::MatrixXd, Eigen::MatrixXd> result = make_pair(Fuu, Fuv);
	return result;
}

void SplineSurface::calc_smooth_term_derivative(VectorXd& spline, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian)
{
	smooth_value = 0;
	int num = 3 * (U_Degree + 1) * (V_Degree + 1);//X维数
	smooth_grad.resize(num);
	smooth_grad.setZero();
	smooth_hessian.resize(num, num);
	smooth_hessian.setZero();
	
	if (smooth_model == 0)
	{
		Vec12 X;
		X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
		{
			X(3 * r).value() = spline(3 * r);
			X(3 * r + 1).value() = spline(3 * r + 1);
			X(3 * r + 2).value() = spline(3 * r + 2);
		}

		// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
		for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
		{
			X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
			X(id).derivatives().setZero();
			X(id).derivatives()(id) = 1;
			X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
		}

		//  3.初始化海塞矩阵;set the hessian matrix to zero
		for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++)
		{
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
			}
		}

		std::vector<vector<Vec3>> P;
		P.resize((U_Degree + 1));
		for (int i = 0; i < P.size(); i++)
		{
			P[i].resize((V_Degree + 1));
		}

		for (int i = 0; i < P.size(); i++)
		{
			P[i][0][0] = X[3 * i];
			P[i][0][1] = X[3 * i + 1];
			P[i][0][2] = X[3 * i + 2];
			P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
			P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
			P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
		}

		//Vec3 s_value;
		//s_value.setZero();
		scalar_t s = 0;
		for (int dim_row = 0; dim_row < U_Degree; dim_row++)
		{
			Vec3 pi0; pi0 = P[dim_row][0];
			Vec3 pi1; pi1 = P[dim_row][1];

			Vec3 pii0; pii0 = P[dim_row + 1][0];
			Vec3 pii1; pii1 = P[dim_row + 1][1];

			Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
			s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
		}
		//scalar_t s = s_value * s_value.transpose();
		smooth_value = s;
		smooth_grad = s.value().derivatives();
		Eigen::MatrixXd B;
		B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
		{
			B.row(r) = s.derivatives()(r).derivatives().transpose();
		}
		smooth_hessian += B;
	}
	else if (smooth_model == 1)
	{
		Vec12 X;
		X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
		{
			X(3 * r).value() = spline(3 * r);
			X(3 * r + 1).value() = spline(3 * r + 1);
			X(3 * r + 2).value() = spline(3 * r + 2);
		}

		// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
		for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
		{
			X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
			X(id).derivatives().setZero();
			X(id).derivatives()(id) = 1;
			X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
		}

		//  3.初始化海塞矩阵;set the hessian matrix to zero
		for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++)
		{
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
			}
		}

		std::vector<vector<Vec3>> P;
		P.resize((U_Degree + 1));
		for (int i = 0; i < P.size(); i++)
		{
			P[i].resize((V_Degree + 1));
		}

		for (int i = 0; i < P.size(); i++)
		{
			P[i][0][0] = X[3 * i];
			P[i][0][1] = X[3 * i + 1];
			P[i][0][2] = X[3 * i + 2];
			P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
			P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
			P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
		}

		//Vec3 s_value;
		//s_value.setZero();
		scalar_t s = 0;
		for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
		{
			for (int dim_row = 0; dim_row < U_Degree - 1; dim_row++)
			{
				Vec3 pij; pij = P[dim_row][dim_col];
				Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
				Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
				Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
				s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
			}
		}

		//scalar_t s = s_value * s_value.transpose();
		smooth_value = s;
		smooth_grad = s.value().derivatives();
		Eigen::MatrixXd B;
		B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
		{
			B.row(r) = s.derivatives()(r).derivatives().transpose();
		}
		smooth_hessian += B;
	}
	else if (smooth_model == 2)
	{
		Vec12 X;
		X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
		{
			X(3 * r).value() = spline(3 * r);
			X(3 * r + 1).value() = spline(3 * r + 1);
			X(3 * r + 2).value() = spline(3 * r + 2);
		}

		// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
		for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
		{
			X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
			X(id).derivatives().setZero();
			X(id).derivatives()(id) = 1;
			X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
		}

		//  3.初始化海塞矩阵;set the hessian matrix to zero
		for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++)
		{
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
			}
		}

		std::vector<vector<Vec3>> P;
		P.resize((U_Degree + 1));
		for (int i = 0; i < P.size(); i++)
		{
			P[i].resize((V_Degree + 1));
		}

		for (int i = 0; i < P.size(); i++)
		{
			P[i][0][0] = X[3 * i];
			P[i][0][1] = X[3 * i + 1];
			P[i][0][2] = X[3 * i + 2];
			P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
			P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
			P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
		}


		scalar_t s = 0;
		for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
		{
			for (int dim_row = 0; dim_row < U_Degree - 1; dim_row++)
			{
				Vec3 pij; pij = P[dim_row][dim_col];
				Vec3 pi1j; pi1j = P[dim_row + 1][dim_col];
				Vec3 pi2j; pi2j = P[dim_row + 2][dim_col];
				Vec3 ruls = (pi2j - pi1j) - (pi1j - pij);
				s += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
			}
		}

		scalar_t s1 = 0;
		for (int dim_row = 0; dim_row < U_Degree; dim_row++)
		{
			Vec3 pi0; pi0 = P[dim_row][0];
			Vec3 pi1; pi1 = P[dim_row][1];

			Vec3 pii0; pii0 = P[dim_row + 1][0];
			Vec3 pii1; pii1 = P[dim_row + 1][1];

			Vec3 ruls = (pii1 - pii0) - (pi1 - pi0);
			s1 += ruls[0] * ruls[0] + ruls[1] * ruls[1] + ruls[2] * ruls[2];
		}


		scalar_t s_total = u_smooth_weight * s + v_smooth_weight * s1;
		smooth_value = s_total;
		smooth_grad = s_total.value().derivatives();
		Eigen::MatrixXd B;
		B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));
		for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
		{
			B.row(r) = s_total.derivatives()(r).derivatives().transpose();
		}
		smooth_hessian += B;
	}
	/*
	for (int ku = 0; ku < u_knot.size(); ku++)
	{
		for (int kv = 0; kv < v_knot.size(); kv++)
		{
			Vec12 X;
			X.resize(3 * (U_Degree + 1) * (V_Degree + 1));
			for (int r = 0; r < (U_Degree + 1) * (V_Degree + 1); r++)
			{
				X(3 * r).value() = spline(3 * r);
				X(3 * r + 1).value() = spline(3 * r + 1);
				X(3 * r + 2).value() = spline(3 * r + 2);
			}

			// 2.初始化一阶导数; repeat partial derivatives for the inner AutoDiffScalar
			for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
			{
				X(id).derivatives().resize(3 * (U_Degree + 1) * (V_Degree + 1));
				X(id).derivatives().setZero();
				X(id).derivatives()(id) = 1;
				X(id).value().derivatives() = inner_derivative_t::Unit(3 * (U_Degree + 1) * (V_Degree + 1), id);
			}

			//  3.初始化海塞矩阵;set the hessian matrix to zero
			for (int idx = 0; idx < 3 * (U_Degree + 1) * (V_Degree + 1); idx++) {
				for (int id = 0; id < 3 * (U_Degree + 1) * (V_Degree + 1); id++)
				{
					X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * (U_Degree + 1) * (V_Degree + 1));
				}
			}

			std::vector<vector<Vec3>> P;
			P.resize((U_Degree + 1));
			for (int i = 0; i < P.size(); i++)
			{
				P[i].resize((V_Degree + 1));
			}

			for (int i = 0; i < P.size(); i++)
			{
				P[i][0][0] = X[3 * i];
				P[i][0][1] = X[3 * i + 1];
				P[i][0][2] = X[3 * i + 2];
				P[i][1][0] = X[3 * (U_Degree + 1) + 3 * i];
				P[i][1][1] = X[3 * (U_Degree + 1) + 3 * i + 1];
				P[i][1][2] = X[3 * (U_Degree + 1) + 3 * i + 2];
			}

			Vec3 s_value;
			s_value.setZero();
			pair<Eigen::MatrixXd, Eigen::MatrixXd> smooth_coff = calc_uv_derivative(u_knot[ku], v_knot[kv]);
			for (int dim_row = 0; dim_row < U_Degree + 1; dim_row++)
			{
				for (int dim_col = 0; dim_col < V_Degree + 1; dim_col++)
				{
					Vec3 pij;pij = P[dim_row][dim_col];
					s_value[0] += smooth_coff.first.coeff(dim_row, dim_col) * pij[0] + smooth_coff.second.coeff(dim_row, dim_col) * pij[0] * knot_weight.coeff(dim_row, dim_col);
					s_value[1] += smooth_coff.first.coeff(dim_row, dim_col) * pij[1] + smooth_coff.second.coeff(dim_row, dim_col) * pij[1] * knot_weight.coeff(dim_row, dim_col);
					s_value[2] += smooth_coff.first.coeff(dim_row, dim_col) * pij[2] + smooth_coff.second.coeff(dim_row, dim_col) * pij[2] * knot_weight.coeff(dim_row, dim_col);
				}
			}
			scalar_t s = s_value * s_value.transpose();
			smooth_value += s;
			smooth_grad += s.value().derivatives();
			Eigen::MatrixXd B;
			B.resize(3 * (U_Degree + 1) * (V_Degree + 1), 3 * (U_Degree + 1) * (V_Degree + 1));
			for (int r = 0; r < 3 * (U_Degree + 1) * (V_Degree + 1); r++)
			{
				B.row(r) = s.derivatives()(r).derivatives().transpose();
			}
			smooth_hessian += B;
		}
	}
	*/
}





void SplineSurface::calc_data_smooth_term_derivative(VectorXd& spline, scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian)
{
	Tree tree; Tree_ tree_;
	/*
	Tree tree;
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();


	Tree_ tree_;
	tree_.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	tree_.build();
	tree_.accelerate_distance_queries();*/


	calc_data_term_energy_derivative(spline, data_value, smooth_value, data_grad, data_hessian, decrease_direction, tree, tree_);
	calc_smooth_term_derivative(spline, smooth_value, smooth_grad, smooth_hessian);
	total_grad = data_grad + smooth_weight * smooth_grad;
	total_hessian = data_hessian + smooth_weight * smooth_hessian;

	Eigen::JacobiSVD<Eigen::MatrixXd> svd(total_hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
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
	decrease_direction = -pseudo_inverse_matrix * total_grad;
	decrease_direction.normalize();
}



void SplineSurface::gcl_calc_data_smooth_term_derivative(const Data& spline, scalar_t& data_value, Eigen::VectorXd& data_grad, Eigen::MatrixXd& data_hessian, Eigen::VectorXd& decrease_direction, scalar_t& smooth_value, Eigen::VectorXd& smooth_grad, Eigen::MatrixXd& smooth_hessian, Eigen::VectorXd& total_grad, Eigen::MatrixXd& total_hessian)
{
	
	VectorXd spline_v;
	spline_v.resize(spline.rows() * 3);
	for (int i = 0; i < spline.rows(); i++)
	{
		
		spline_v[3 * i] = spline.coeff(i, 0);
		spline_v[3 * i + 1] = spline.coeff(i, 1);
		spline_v[3 * i + 2] = spline.coeff(i, 2);
	}
	
	Tree tree; Tree_ tree_;
	/*
	Tree tree;
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();


	Tree_ tree_;
	tree_.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	tree_.build();
	tree_.accelerate_distance_queries();*/

	
	calc_data_term_energy_derivative(spline_v, data_value, smooth_value, data_grad, data_hessian, decrease_direction, tree, tree_);
	
	calc_smooth_term_derivative(spline_v, smooth_value, smooth_grad, smooth_hessian);

	/*
	total_grad = data_grad + smooth_weight * smooth_grad;
	total_hessian = data_hessian + smooth_weight * smooth_hessian;


	Eigen::JacobiSVD<Eigen::MatrixXd> svd(total_hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
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
	decrease_direction = -pseudo_inverse_matrix * total_grad;
	decrease_direction.normalize();*/
}


void SplineSurface::gcl_bezier_visualization(int row_num, Data& spline, int iter_num, int bezier_num)
{
	vector<vector<Eigen::Vector3d>> control_points;
	int cp_row = spline.rows() / 2;
	control_points.resize(cp_row);
	for (int i = 0; i < cp_row; i++)
	{
		control_points[i].resize(2);
		control_points[i][0] = spline.row(i);
		control_points[i][1] = spline.row(i + cp_row);
	}

	//1.基于采样
	vector<double> u_vec = generate_equally_spaced_vector(row_num);
	vector<double> v_vec = generate_equally_spaced_vector(row_num);

	//2.Bezier_info构造
	Matrix<SplineSurface::Bezier_info, Eigen::Dynamic, Eigen::Dynamic> bezier_info_mat;
	bezier_info_mat.resize(row_num, row_num);
	for (int i = 0; i < u_vec.size(); i++)
	{
		for (int j = 0; j < v_vec.size(); j++)
		{
			Bezier_info bi;
			bi.u_cord = u_vec[i];
			bi.v_cord = v_vec[j];
			bi.u_idx = i;
			bi.v_idx = j;
			bi.bezier_point = calc_bezier_surface_value(U_Degree, V_Degree, bi.u_cord, bi.v_cord, control_points);
			bezier_info_mat(i, j) = bi;
		}
	}

	//3.
	Mesh outmesh;
	Matrix<int, Eigen::Dynamic, Eigen::Dynamic> v_idx_mat;
	v_idx_mat.resize(row_num, row_num);
	int init_v_idx = 0;
	for (int i = 0; i < bezier_info_mat.rows(); i++)
	{
		for (int j = 0; j < bezier_info_mat.cols(); j++)
		{
			Mesh::VertexHandle v0 = outmesh.add_vertex(Mesh::Point(bezier_info_mat(i, j).bezier_point.data()));
			v_idx_mat(i, j) = init_v_idx;
			init_v_idx++;
		}
	}

	for (int i = 0; i < row_num - 1; ++i)
	{
		for (int j = 0; j < row_num - 1; ++j)
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
	MeshTools::WriteMesh(mesh, "./inputmesh.obj");
	MeshTools::WriteMesh(outmesh, "./outputmesh_" + to_string(iter_num) + "_" + to_string(bezier_num) + ".obj");
}


std::set<int> SplineSurface::search_projection_area(Data& spline, int sample_num)
{
	vector<vector<Eigen::Vector3d>> control_points;
	vector<int> mesh_area_face_id;
	mesh_area_face_id.clear();
	int cp_row = spline.rows() / 2;
	control_points.resize(cp_row);
	for (int i = 0; i < cp_row; i++)
	{
		control_points[i].resize(2);
		control_points[i][0] = spline.row(i);
		control_points[i][1] = spline.row(i + cp_row);
	}

	//1.按照法相向网格发射线
	//1.1 确定法相
	tree.insert(cgal_triangles.begin(), cgal_triangles.end());
	tree.build();
	tree.accelerate_distance_queries();

	tree_.insert(cgal_triangles_.begin(), cgal_triangles_.end());
	tree_.build();
	tree_.accelerate_distance_queries();

	Point p1(control_points[0][0].x(), control_points[0][0].y(), control_points[0][0].z());
	Point p2(control_points[0][1].x(), control_points[0][1].y(), control_points[0][1].z());
	Point p3(control_points[1][0].x(), control_points[1][0].y(), control_points[1][0].z());
	K::Vector_3 normal = CGAL::normal(p1, p2, p3);
	Eigen::Vector3d plane_normal(normal.x(), normal.y(), normal.z());

	plane_normal.normalize();
	Eigen::Vector3d start_point = calc_bezier_surface_value(U_Degree, V_Degree, 0.5, 0.5, control_points);
	Eigen::Vector3d end_point = start_point + plane_normal;
	Point ray_source(start_point.x(), start_point.y(), start_point.z()); // 射线源点坐标
	Point ray_direction(end_point.x(), end_point.y(), end_point.z()); // 射线的方向向量
	Ray ray_query(ray_source, ray_direction);

	Ray_intersection intersection_test = tree_.first_intersection(ray_query);
	if (intersection_test)
	{

	}
	else
	{
		plane_normal = -1 * plane_normal;
	}

	//2.对贝塞尔曲面均匀采样
	vector<double> u_vec = generate_equally_spaced_vector(sample_num);
	vector<double> v_vec = generate_equally_spaced_vector(sample_num);
	//2.1.Bezier_info构造
	Matrix<SplineSurface::Bezier_info, Eigen::Dynamic, Eigen::Dynamic> bezier_info_mat;
	bezier_info_mat.resize(sample_num, sample_num);
	for (int i = 0; i < u_vec.size(); i++)
	{
		for (int j = 0; j < v_vec.size(); j++)
		{
			Bezier_info bi;
			bi.u_cord = u_vec[i];
			bi.v_cord = v_vec[j];
			bi.u_idx = i;
			bi.v_idx = j;
			bi.bezier_point = calc_bezier_surface_value(U_Degree, V_Degree, bi.u_cord, bi.v_cord, control_points);
			Eigen::Vector3d end_ = bi.bezier_point + plane_normal;
			int face_id = spline_ray_intersection_(bi.bezier_point, end_);
			bezier_info_mat(i, j) = bi;
			mesh_area_face_id.push_back(face_id);

		}
	}

	mesh_area_face_id.erase(std::remove(mesh_area_face_id.begin(), mesh_area_face_id.end(), -1), mesh_area_face_id.end());
	std::set<int> uniqueNumbers(mesh_area_face_id.begin(), mesh_area_face_id.end());
	return uniqueNumbers;
	//3.可能不连通，搜索1或2邻域，直至全覆盖

}

int SplineSurface::spline_ray_intersection_(Eigen::Vector3d start_point, Eigen::Vector3d end_point)
{
	Point ray_source(start_point.x(), start_point.y(), start_point.z());
	Point ray_direction(end_point.x(), end_point.y(), end_point.z());
	Ray ray_query(ray_source, ray_direction);
	Ray_intersection intersection = tree_.first_intersection(ray_query);

	auto result_closest = tree.closest_point_and_primitive(ray_source);
	Point p_closest_point = result_closest.first;
	double distance_closest = std::sqrt(CGAL::squared_distance(p_closest_point, ray_source));

	if (intersection)
	{
		if (boost::get<Point>(&(intersection->first)))
		{
			const Point* p = boost::get<Point>(&(intersection->first));
			Point p_close(p->x(), p->y(), p->z());

			double distance_pro = std::sqrt(CGAL::squared_distance(p_close, ray_source));

			if (distance_pro <=1.2* distance_closest)
			{
				auto result = tree.closest_point_and_primitive(p_close);
				return result.second->index;
			}
			else
			{
				return -1;
			}

		}
	}
	return -1;
}

void SplineSurface::fitting_init(int sample_num)
{
	face_point_temp.clear();
	face_point_temp.resize(mesh.n_faces());
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		vector<Point> v_p_3; v_p_3.clear();
		FH fh = mesh.face_handle(i);
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}
		Point A(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()), B(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()), C(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z());
		Triangle triangle(A, B, C);
		std::vector<Point> sample_points;
		CGAL::Random_points_in_triangle_3<Point> rand_points_in_triangle(triangle);
		for (unsigned int j = 0; j < sample_num; ++j)
		{
			sample_points.push_back(*rand_points_in_triangle++);
		}

		std::vector<Vec3> sample_points_vec3; sample_points_vec3.clear();
		for (unsigned int j = 0; j < sample_points.size(); ++j)
		{
			Vec3 p;
			p[0] = sample_points[j].x();
			p[1] = sample_points[j].y();
			p[2] = sample_points[j].z();
			sample_points_vec3.push_back(p);
		}
		face_point_temp[i] = sample_points_vec3;
	}
}

