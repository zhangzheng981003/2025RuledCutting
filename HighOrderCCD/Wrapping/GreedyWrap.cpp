#include "GreedyWrap.h"
#include "GreedyWrap.h"
#include <algorithm >
using namespace std;

GreedyWrap::GreedyWrap(Mesh& mesh):mesh(mesh)
{
	//
}

GreedyWrap::~GreedyWrap()
{
}

void GreedyWrap::initialization()
{
	//构建aabb树

	my_cgal_triangles.clear();
	cgal_triangles.clear();
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		vector<Point> v_p_3;
		v_p_3.clear();
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(mesh.point(*fv_it)[0], mesh.point(*fv_it)[1], mesh.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}

		cgal_triangles.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()));

		auto f_normal = mesh.calc_face_normal(*f_it);
		Vec3 f_n;
		f_n[0] = f_normal[0];
		f_n[1] = f_normal[1];
		f_n[2] = f_normal[2];
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

	//储存每个面的基本信息
	face_basis_info_temp.clear();
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		SegSurface::FaceBasisInfo fhi;
		FH fh = mesh.face_handle(i);
		fhi.face_area = mesh.calc_face_area(fh);
		Mesh::Point fc = mesh.calc_face_centroid(fh);
		Mesh::Point fn = mesh.calc_face_normal(fh);
		fhi.face_normal = { fn[0], fn[1], fn[2] }; fhi.face_normal.normalize();
		fhi.centeriod = { fc[0], fc[1], fc[2] };
		fhi.face_idx = i;
		face_basis_info_temp.push_back(fhi);
	}

	// 平均长度
	average_edge_length = 0.0;
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		average_edge_length += mesh.calc_edge_length(mesh.edge_handle(i));
	}
	average_edge_length = average_edge_length / mesh.n_edges();

	// bbx长度===采样的直母线长度
	Mesh::Point ptMin; Mesh::Point ptMax;
	ptMin[0] = ptMin[1] = ptMin[2] = DBL_MAX;
	ptMax[0] = ptMax[1] = ptMax[2] = -DBL_MAX;
	for (const auto& vh : mesh.vertices())
	{
		ptMin.minimize(mesh.point(vh));
		ptMax.maximize(mesh.point(vh));
	}
	rulings_length = (ptMax - ptMin).norm();

	//clear 向量
	vertex_fitting_temp.clear();
	inaccessible_area.clear();
	sucess_cover_point.clear();
	failed_cover_point.clear();
	inaccessible_point.clear();
	inaccessible_point_status.clear();
	inaccessible_point_status.resize(mesh.n_vertices(), true);
	inaccessible_face_status.clear();
	inaccessible_face_status.resize(mesh.n_faces(), true);
	
	// 分块采用的是1.5* average_edge_length
	// 寻找不可到达的点
	search_inaccessible_area(2.0 * average_edge_length, 300);

	/*
	vector<vector<double>> id1;
	read_txt("vaseid1_.txt", id1);

	for (int i = 0; i < id1.size(); i++)
	{
		if (inaccessible_face_status[id1[i][0]] == 1)
		{
			inaccessible_face_status[id1[i][0]] = 0;
			inaccessible_area.push_back(id1[i][0]);
		}
	}*/
	vector<vector<double>> id2;
	read_txt(result_path_seg +".txt", id2);
	
	for (int i = 0; i < id2.size(); i++)
	{
		if (inaccessible_face_status[id2[i][0]] == 1)
		{
			inaccessible_face_status[id2[i][0]] = 0;
			inaccessible_area.push_back(id2[i][0]);
		}
	}
	
	write_txt("inaccessible_area.txt", inaccessible_area);

	// 每个顶点的邻域面
	vertex_temp.resize(mesh.n_vertices());
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		vertex_temp[i].v_idx = i;
		vertex_temp[i].f_idx_temp.clear();
		for (auto vf : mesh.vf_range(mesh.vertex_handle(i)))
		{
			vertex_temp[i].f_idx_temp.push_back(vf.idx());
		}
	}

}

//点集转成面集
vector<int> GreedyWrap::vertex2face(vector<int>& vertexId)
{
	vector<int> faceId; faceId.clear();
	for (auto ele : vertexId)
	{
		for (int i = 0; i < vertex_temp[ele].f_idx_temp.size(); i++)
		{
			faceId.push_back(vertex_temp[ele].f_idx_temp[i]);
		}
	}

	std::sort(faceId.begin(), faceId.end());
	auto last = std::unique(faceId.begin(), faceId.end());
	faceId.erase(last, faceId.end());

	return faceId;
}

int GreedyWrap::select_big_patches(vector<SegSurface::Path>& pathes_temp, double coverage_rate_threshold, int cover_face_nums)
{
	std::sort(pathes_temp.begin(), pathes_temp.end(), cmp_smaller_path_temp);
	double coverage_rate = 0.0;
	vector<int> is_cover(mesh.n_faces(), 0);
	vector<int> is_cover_temp(mesh.n_faces(), 0);
	set<int> cover_face_set; cover_face_set.clear();
	int first_big_patches_nums = 0;
	auto path1 = pathes_temp;
	path1.clear();
	//cout<<"dsadaad" << pathes_temp.size() << endl;
	for (int i = 0; i < pathes_temp.size(); i++)
	{


		int chosen_num = 0;
		is_cover_temp.clear();
		is_cover_temp.shrink_to_fit();
		is_cover_temp = is_cover;
		for (auto ele : pathes_temp[i].cover_face)
		{
			//if (i == 1) cout << is_cover[ele] << " "<<ele<<endl;
			if (is_cover[ele] == 0)
			{
				is_cover[ele] = 1;
				chosen_num++;
			}
		}
		if (chosen_num > pathes_temp[i].cover_face.size() / 10&& pathes_temp[i].info_every_face.size()>4)
		{
			cout << "chosen_num=" << chosen_num << ";pathes_temp[i].cover_face.size()= " << pathes_temp[i].cover_face.size() << endl;
			for (auto ele : pathes_temp[i].cover_face)
			{
				cover_face_set.insert(ele);
			}
			
			path1.push_back(pathes_temp[i]);
		}
		else
		{
			pathes_temp[i].cover_face.clear();
			is_cover = is_cover_temp;
		}
		coverage_rate = (double(cover_face_set.size()) / double(cover_face_nums));
		cout << cover_face_set.size() << endl;
		if (coverage_rate >= coverage_rate_threshold)
		{
			first_big_patches_nums = i;
			break;
		}


	}

	if (first_big_patches_nums == 0)
	{
		cout << "采样较少！！！" << endl;
		
	}
	else {  }
	pathes_temp = path1;
	first_big_patches_nums = pathes_temp.size();

	ofstream ios("cover_id.txt");
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		if (is_cover[i] == 1)
		{
			ios << i << endl;
		}
	}
	ios.close();

	return first_big_patches_nums;
}

void GreedyWrap::rotate_vector(Vector3d& n, double angle, Vector3d& v, Vector3d& result)
{
	double cosTheta = cos(angle);
	double sinTheta = sin(angle);
	n.normalize();
	v.normalize();
	result[0] = (cosTheta + (1 - cosTheta) * n[0] * n[0]) * v[0] +
		((1 - cosTheta) * n[0] * n[1] - sinTheta * n[2]) * v[1] +
		((1 - cosTheta) * n[0] * n[2] + sinTheta * n[1]) * v[2];

	result[1] = ((1 - cosTheta) * n[0] * n[1] + sinTheta * n[2]) * v[0] +
		(cosTheta + (1 - cosTheta) * n[1] * n[1]) * v[1] +
		((1 - cosTheta) * n[1] * n[2] - sinTheta * n[0]) * v[2];

	result[2] = ((1 - cosTheta) * n[0] * n[2] - sinTheta * n[1]) * v[0] +
		((1 - cosTheta) * n[1] * n[2] + sinTheta * n[0]) * v[1] +
		(cosTheta + (1 - cosTheta) * n[2] * n[2]) * v[2];

	result.normalize();
}

void GreedyWrap::write_txt(const string& file_path, vector<int>& input_vec)
{
	ofstream out(file_path);
	for (auto ele : input_vec)
	{
		out << ele << endl;
	}
	out.close();
}

void GreedyWrap::read_txt(const string& pathname, vector<vector<double>>& res)
{
	res.clear();
	ifstream infile;
	infile.open(pathname.data());
	assert(infile.is_open());
	vector<double> suanz;
	string s;
	while (getline(infile, s))
	{
		istringstream is(s);
		double d;
		while (!is.eof())
		{
			is >> d;
			suanz.push_back(d);
		}
		res.push_back(suanz);
		suanz.clear();
		s.clear();
	}
	infile.close();
}

bool GreedyWrap::is_segment_intersection(Eigen::Vector3d start_point, Eigen::Vector3d end_point)
{
	Point p0(start_point.x(), start_point.y(), start_point.z());
	Point p1(end_point.x(), end_point.y(), end_point.z());
	Segment seg_line(p0, p1);
	if (tree.do_intersect(seg_line))
	{
		return true;
	}
	else
	{
		return false;
	}
}

void GreedyWrap::search_inaccessible_area(double high_value, int sample_rulings_num_one_face)
{
	double rotate_angle = double(double(M_PI) / double(sample_rulings_num_one_face));
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		//1.旋转向量
		Vector3d cp = face_basis_info_temp[i].centeriod;
		Vector3d fn = face_basis_info_temp[i].face_normal;
		VH fv = mesh.fv_begin(mesh.face_handle(i));
		Vec3d p = mesh.point(fv);
		Vector3d p0 = { p[0],p[1],p[2] };
		Vector3d start_vec = (p0 - cp).normalized();
		
		//2.判断rulings求交
		int intersection_num = 0;
		cp = cp + high_value * fn;
		for (int rn = 0; rn < sample_rulings_num_one_face; rn++)
		{
			Vector3d ro_vec;
			rotate_vector(fn, rotate_angle * rn, start_vec, ro_vec);
			Eigen::Vector3d start_point = cp + 0.5 * rulings_length * ro_vec;
			Eigen::Vector3d end_point = cp - 0.5 * rulings_length * ro_vec;
			if (is_segment_intersection(start_point, end_point))
			{
				intersection_num++;
			}
			else
			{
				break;
			}
		}

		if (intersection_num == sample_rulings_num_one_face)
		{
			inaccessible_area.push_back(i);
		}
	}

	for (int i = 0; i < inaccessible_area.size(); i++)
	{
		inaccessible_face_status[inaccessible_area[i]] = false;
		FH fh = mesh.face_handle(inaccessible_area[i]);
		for (auto fv : mesh.fv_range(fh))
		{
			inaccessible_point.push_back(fv.idx());
		}
	}

	std::sort(inaccessible_point.begin(), inaccessible_point.end());
	auto last = std::unique(inaccessible_point.begin(), inaccessible_point.end());
	inaccessible_point.erase(last, inaccessible_point.end());

	for (auto ele : inaccessible_point)
	{
		inaccessible_point_status[ele] = false;
	}
}

void GreedyWrap::cut_mesh(Mesh& mesh, Mesh& seg_mesh, SegSurface::Path& path)
{
	//2.CutMesh：ori_face_id --> cutMesh_face_id:以及每一个面的拟合区域
	vector<int> new_patch(mesh.n_faces(), 0);
	ofstream ioszhzh("mesh_data/final_seg.txt");
	for (int i = 0; i < path.info_every_face.size(); i++)
	{
		new_patch[path.info_every_face[i].face_idx] = 1;
	}
	for (int i = 0; i < new_patch.size(); i++)
	{
		ioszhzh << new_patch[i] << endl;
	}
	ioszhzh.close();

	std::string  seg_face_path = "mesh_data/final_seg.txt";
	Mesh output_mesh = mesh;
	int nf = output_mesh.n_faces();
	std::vector<int> face_status(nf, -1);
	std::ifstream face_status_(seg_face_path);
	int value;
	int line_number = 0;
	while (face_status_ >> value)
	{
		face_status[line_number] = value;
		line_number++;
	}

	int ne = output_mesh.n_edges();
	std::vector<int> seam_status;
	for (auto& e : output_mesh.edges())
	{
		if (!e.is_boundary())
		{
			OpenMesh::HalfedgeHandle heh = output_mesh.halfedge_handle(e, 0);
			OpenMesh::FaceHandle fh1 = output_mesh.face_handle(heh);

			OpenMesh::HalfedgeHandle opp_heh = output_mesh.opposite_halfedge_handle(heh);
			OpenMesh::FaceHandle fh2 = output_mesh.face_handle(opp_heh);

			if (face_status[fh1.idx()] != face_status[fh2.idx()]) seam_status.push_back(e.idx());
		}
	}

	CutMesh cut_mesh("mesh_data/", output_mesh, seam_status);
	Mesh cut_output_mesh = cut_mesh.Compute();

	for (int i = 0; i < cut_mesh.res_two_patches.size(); i++)
	{
		if (cut_mesh.res_two_patches[i].n_faces() == path.info_every_face.size())
		{
			seg_mesh = cut_mesh.res_two_patches[i];
			break;
		}
	}
}

void GreedyWrap::build_bspline_surface_handle(Mesh& mesh, Handle(Geom_BSplineSurface)& bspline_handle, Eigen::MatrixXd& ctr_points, vector<double>& u_knots_temp)
{
	BSplineFitting* sf = new BSplineFitting(mesh);
	int curve_ctr_num = ctr_points.rows() / 2;
	vector<BSplineFitting::CtrInfo> ctr_point1; vector<BSplineFitting::CtrInfo> ctr_point2;
	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = ctr_points.coeffRef(j, 0);
		Standard_Real py = ctr_points.coeffRef(j, 1);
		Standard_Real pz = ctr_points.coeffRef(j, 2);
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point1.push_back(ci);
	}

	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = ctr_points.coeffRef(j + curve_ctr_num, 0);
		Standard_Real py = ctr_points.coeffRef(j + curve_ctr_num, 1);
		Standard_Real pz = ctr_points.coeffRef(j + curve_ctr_num, 2);
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point2.push_back(ci);
	}

	sf->u_knot_temp.clear(); sf->u_knot_temp = u_knots_temp;
	bspline_handle = nullptr;
	sf->create_BSpline_surface(ctr_point1, ctr_point2, bspline_handle);
	delete sf;
	sf = nullptr;
}

void GreedyWrap::calc_seg_energy(Mesh& mesh, Mesh& seg_mesh, SegSurface::Path path, Handle(Geom_BSplineSurface)& bspline_handle
	, Eigen::MatrixXd& ctr_points, vector<double>& u_knots_temp, vector<SegEnergy::InnerEdge>& ie_temp, vector<double>& rulings_energy
	, vector<int>& real2img, vector<int>& img2real)
{
	BSplineFitting* sf = new BSplineFitting(mesh);
	int curve_ctr_num = ctr_points.rows() / 2;
	vector<BSplineFitting::CtrInfo> ctr_point1; vector<BSplineFitting::CtrInfo> ctr_point2;
	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = ctr_points.coeffRef(j, 0);
		Standard_Real py = ctr_points.coeffRef(j, 1);
		Standard_Real pz = ctr_points.coeffRef(j, 2);
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point1.push_back(ci);
	}

	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = ctr_points.coeffRef(j + curve_ctr_num, 0);
		Standard_Real py = ctr_points.coeffRef(j + curve_ctr_num, 1);
		Standard_Real pz = ctr_points.coeffRef(j + curve_ctr_num, 2);
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point2.push_back(ci);
	}
	sf->u_knot_temp.clear(); sf->u_knot_temp = u_knots_temp;
	bspline_handle = nullptr;
	sf->create_BSpline_surface(ctr_point1, ctr_point2, bspline_handle);
	
	//int curve_ctr_num = ctr_points.rows() / 2;
	std::vector<vector<Vector3d>> P;
	P.resize(curve_ctr_num);
	for (int i = 0; i < P.size(); i++)
	{
		P[i].resize(2);
	}

	for (int i = 0; i < P.size(); i++)
	{
		P[i][0][0] = ctr_points.coeffRef(i, 0);
		P[i][0][1] = ctr_points.coeffRef(i, 1);
		P[i][0][2] = ctr_points.coeffRef(i, 2);
		P[i][1][0] = ctr_points.coeffRef(i + curve_ctr_num, 0);
		P[i][1][1] = ctr_points.coeffRef(i + curve_ctr_num, 1);
		P[i][1][2] = ctr_points.coeffRef(i + curve_ctr_num, 2);
	}

	Tree seg_tree;
	std::vector<MyTriangle> segMesh_cgal_triangles;
	my_cgal_triangles.clear();
	for (Mesh::FaceIter f_it = seg_mesh.faces_begin(); f_it != seg_mesh.faces_end(); ++f_it)
	{
		vector<Point> v_p_3;
		v_p_3.clear();
		for (Mesh::FaceVertexIter fv_it = seg_mesh.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(seg_mesh.point(*fv_it)[0], seg_mesh.point(*fv_it)[1], seg_mesh.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}

		auto f_normal = seg_mesh.calc_face_normal(*f_it);
		Vec3 f_n;
		f_n[0] = f_normal[0];
		f_n[1] = f_normal[1];
		f_n[2] = f_normal[2];
		segMesh_cgal_triangles.emplace_back(
			Point(v_p_3[0].x(), v_p_3[0].y(), v_p_3[0].z()),
			Point(v_p_3[1].x(), v_p_3[1].y(), v_p_3[1].z()),
			Point(v_p_3[2].x(), v_p_3[2].y(), v_p_3[2].z()), f_it->idx(), f_n);
		v_p_3.clear();
	}
	seg_tree.insert(segMesh_cgal_triangles.begin(), segMesh_cgal_triangles.end());
	seg_tree.build();
	seg_tree.accelerate_distance_queries();

	//vector<int> real2img; 
	real2img.resize(mesh.n_faces(), 0);
	//vector<int> img2real; 
	img2real.resize(path.info_every_face.size(), 0);
	for (int i = 0; i < path.info_every_face.size(); i++)
	{
		Point p0 = { path.info_every_face[i].centeriod[0]  ,path.info_every_face[i].centeriod[1]  ,path.info_every_face[i].centeriod[2] };
		auto result = seg_tree.closest_point_and_primitive(p0);
		real2img[path.info_every_face[i].face_idx] = result.second->index;
		img2real[result.second->index] = path.info_every_face[i].face_idx;
	}


	vector<int> id_in_info_every_face; id_in_info_every_face.clear(); id_in_info_every_face.resize(mesh.n_faces(), 0);
	for (int i = 0; i < path.info_every_face.size(); i++)
	{
		id_in_info_every_face[path.info_every_face[i].face_idx] = i;
	}

	vector<SegEnergy::AdjFace> adj_face_temp; adj_face_temp.clear(); adj_face_temp.resize(mesh.n_faces());
	for (int i = 0; i < seg_mesh.n_faces(); i++)
	{
		SegEnergy::AdjFace af;
		af.face_idx = i;
		af.rulling_f0_ori = { path.info_every_face[id_in_info_every_face[img2real[i]]].ruling_direction[0],path.info_every_face[id_in_info_every_face[img2real[i]]].ruling_direction[1],path.info_every_face[id_in_info_every_face[img2real[i]]].ruling_direction[2] };
		FH fh = seg_mesh.face_handle(i);
		Vec3d p_ = seg_mesh.calc_face_centroid(fh);
		gp_Pnt pointToProject(p_[0], p_[1], p_[2]);
		GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, bspline_handle);
		if (projectPoint.NbPoints() > 0)
		{
			Standard_Real u, v;
			projectPoint.LowerDistanceParameters(u, v);
			Standard_Real dist;
			dist = projectPoint.LowerDistance();
			int u_interval_order = sf->u_k_position(u, sf->u_knot_temp);

			vector<double> u_basis_temp;
			if (u_interval_order == 0)
			{
				vector<double> basis_temp(sf->poly_degree + 1, 0);
				basis_temp[0] = 1;
				u_basis_temp = basis_temp;
			}
			else if (u_interval_order == sf->u_knot_temp.size() - 1)
			{
				vector<double> basis_temp(sf->poly_degree + 1, 0);
				basis_temp[sf->poly_degree] = 1;
				u_basis_temp = basis_temp;
			}
			else
			{
				sf->calc_basis_fuction(u, u_interval_order, sf->poly_degree, u_basis_temp);
			}
			Vector3d rulling;
			rulling.setZero();
			if (u_interval_order == 0)
			{
				Vector3d pi0, pi1;
				pi0 = P[0][0];
				pi1 = P[0][1];
				rulling = pi1 - pi0;
			}
			else if (u_interval_order == sf->u_knot_temp.size() - 1)
			{
				Vector3d pi0, pi1;
				pi0 = P[curve_ctr_num - 1][0];
				pi1 = P[curve_ctr_num - 1][1];
				rulling = pi1 - pi0;
			}
			else
			{
				for (int i = 0; i < u_basis_temp.size(); i++)
				{
					Vector3d pi0, pi1;
					pi0 = P[u_interval_order - sf->poly_degree + i][0];
					pi1 = P[u_interval_order - sf->poly_degree + i][1];
					rulling += u_basis_temp[i] * (pi1 - pi0);
				}
			}
			af.rulling_f0_bspline = rulling;
		}
		adj_face_temp[i] = af;
	}

	cout << "flag" << 123 << endl;

	ie_temp.clear();
	for (int k0 = 0; k0 < seg_mesh.n_edges(); k0++)
	{
		EH eh = seg_mesh.edge_handle(k0);
		if (!seg_mesh.is_boundary(eh))
		{
			SegEnergy::InnerEdge ie;
			ie.inner_edge_idx = eh.idx();
			HEH heh0 = seg_mesh.halfedge_handle(eh, 0);
			HEH heh1 = seg_mesh.halfedge_handle(eh, 1);

			FH f0 = seg_mesh.face_handle(heh0);
			FH f1 = seg_mesh.face_handle(heh1);

			ie.f0_id = f0.idx();
			ie.rulling_f0_ori = adj_face_temp[ie.f0_id].rulling_f0_ori;
			ie.rulling_f0_bspline = adj_face_temp[ie.f0_id].rulling_f0_bspline;

			ie.f1_id = f1.idx();
			ie.rulling_f1_ori = adj_face_temp[ie.f1_id].rulling_f0_ori;
			ie.rulling_f1_bspline = adj_face_temp[ie.f1_id].rulling_f0_bspline;

			set<int> v_id_temp; v_id_temp.clear();
			for (auto vv : seg_mesh.fv_range(f0))
			{
				v_id_temp.insert(vv.idx());
			}
			for (auto vv : seg_mesh.fv_range(f1))
			{
				v_id_temp.insert(vv.idx());
			}
			ie.Adj_v = v_id_temp;

			double total_dist = 0.0;
			for (auto ele : v_id_temp)
			{
				Vec3d p0 = seg_mesh.point(seg_mesh.vertex_handle(ele));
				gp_Pnt pointToProject(p0[0], p0[1], p0[2]);
				GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, bspline_handle);
				if (projectPoint.NbPoints() > 0)
				{
					Standard_Real u, v;
					projectPoint.LowerDistanceParameters(u, v);
					Standard_Real dist;
					dist = projectPoint.LowerDistance();
					total_dist += dist;
				}
			}
			ie.vertices_dist = total_dist;
			ie_temp.push_back(ie);
		}
	}
}

void GreedyWrap::calc_vertex_dist_for_bspline_surfaces(vector<vector<vector<BSplineSurfaceInfo>>>& bspline_res_3temp, vector<VertexFittingStatus>& vertex_fitting_temp, int start_pos)
{
	for (int wrap_idx = start_pos; wrap_idx < bspline_res_3temp.size(); wrap_idx++)
	{
		for (int branch_idx = 0; branch_idx < bspline_res_3temp[wrap_idx].size(); branch_idx++)
		{
			for (int surface_idx = 0; surface_idx < bspline_res_3temp[wrap_idx][branch_idx].size(); surface_idx++)
			{
				for (int v_id = 0; v_id < mesh.n_vertices(); v_id++)
				{
					if (inaccessible_point_status[v_id])
					{
						VertexFittingStatus vfs;
						vfs.vertex_idx = v_id;
						Vec3d p = mesh.point(mesh.vertex_handle(v_id));
						gp_Pnt point(p[0], p[1], p[2]);
						GeomAPI_ProjectPointOnSurf projectPoint(point, bspline_res_3temp[wrap_idx][branch_idx][surface_idx].bspline_handle);
						if (projectPoint.NbPoints() > 0)
						{
							gp_Pnt closestPoint = projectPoint.NearestPoint();
							Standard_Real distance = point.Distance(closestPoint);
							vfs.vertex_bspline_dist_temp.push_back(make_pair(bspline_res_3temp[wrap_idx][branch_idx][surface_idx], distance));
							vertex_fitting_temp.push_back(vfs);
						}
						else
						{
							std::cout << "计算点 (" << v_id << ") 距离失败" << endl;
						}
					}
				}
			}
		}
	}

	for (auto ele : vertex_fitting_temp)
	{
		std::sort(ele.vertex_bspline_dist_temp.begin(), ele.vertex_bspline_dist_temp.end(), cmp_bigger);
	}
}

void GreedyWrap::greedy_wrap()
{
	//1.初始化
	//1.1 连通分支
	vector<vector<int>> connect_branches;
	connect_branches.clear(); connect_branches.resize(1);
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		connect_branches[0].push_back(i);
	}
	
	//1.2 
	initialization();
	int wrapping_idx = 0;
	vector<vector<vector<GreedyWrap::BSplineSurfaceInfo>>> bspline_res_3temp; bspline_res_3temp.clear();
	// wrapping_idx ----> branch_idx ----> surface_idx

	//2.
	int cut_num = 0;
	while (wrapping_idx < 1)
	{
		vector<vector<BSplineSurfaceInfo>>  bspline_res_2temp; bspline_res_2temp.clear();
		for (int i = 0; i < connect_branches.size(); i++)
		{
			vector<GreedyWrap::BSplineSurfaceInfo> bspline_res_1temp; bspline_res_1temp.clear();
			//1.分块
			SegSurface ss(mesh);
			vector<SegSurface::Path> pathes_temp;
			vector<int> branch_accessible = connect_branches[i];
			ss.seg_surface_no_iter(connect_branches[i], pathes_temp, inaccessible_area);
			//debug_read_pathes_temp("D:/programfiles/3/11/build/result/"+ result_path_fitting+"/pathes_temp/", pathes_temp);
			for (int i = 0; i < pathes_temp[0].info_every_face.size(); i++)
			{
				cout << pathes_temp[0].info_every_face.size() << " " << i << " " << pathes_temp[0].info_every_face[i].face_idx << endl;
			}
			//2.拟合
			//2.1选取大块拟合,计算块数
			branch_accessible.erase(std::remove_if(branch_accessible.begin(), branch_accessible.end(), [&](int num) {
				return std::find(inaccessible_area.begin(), inaccessible_area.end(), num) != inaccessible_area.end();
				}), branch_accessible.end());

			cout << "flag" <<pathes_temp[0].cover_face.size()<< endl;
			int fitting_nums = select_big_patches(pathes_temp, 0.95, branch_accessible.size());

		    //fitting_nums = 1;
			cout << "zz: " << fitting_nums << endl;
			vector<SegSurface::Path> pathes_temp_copy(pathes_temp.begin(), pathes_temp.begin() + fitting_nums);
			debug_write_pathes_temp("D:/programfiles/3/11/build/result/" + result_path_fitting + "/pathes_temp/", pathes_temp_copy);
			pathes_temp.clear(); pathes_temp.shrink_to_fit();
			pathes_temp = pathes_temp_copy;
			ss.print_seg(pathes_temp);
			
			//2.2开始拟合
			for (int patch_num = 0; patch_num < fitting_nums && pathes_temp[patch_num].cover_face.size()>0; patch_num++)
			{
				BSplineNoCollision bnc(mesh);
				bnc.bspline_fitting_init_plane(pathes_temp[patch_num]);

				//返回拟合结果：数据能量与距离
				auto bspline_res = bnc.bspline_fitting_no_collision(result_path_fitting + ".obj", mesh, patch_num, pathes_temp[patch_num]);

				BSplineSurfaceInfo bsi;
				bsi.ctr_points = bspline_res.spline;
				bsi.u_knots_temp = bspline_res.u_knot;
				bsi.dist_energy = bspline_res.data_e;
				bsi.max_dist = bspline_res.max_e;
				bsi.mean_dist = bspline_res.arv_e;
				bsi.path = pathes_temp[patch_num];
				bsi.wrapping_idx = wrapping_idx;
				bsi.surface_idx = patch_num;
				bsi.branch_idx = i;

				bspline_res_1temp.push_back(bsi);


				if (bsi.mean_dist > 8e-4)
				{
					// 1.切割网格，构造bspline曲面handle
					Mesh seg_mesh; vector<int> real2img; vector<int> img2real;
					cut_mesh(mesh, seg_mesh, bsi.path);

					vector<SegEnergy::InnerEdge> ie_temp; vector<double> rulings_energy;
					calc_seg_energy(mesh, seg_mesh, bsi.path, bsi.bspline_handle, bsi.ctr_points, bsi.u_knots_temp,
						ie_temp, rulings_energy, real2img, img2real);

					SegEnergy::calc_edge_energy0(seg_mesh, 1.0, rulings_energy, ie_temp);

					//2.Seg: 谱分割
					SegMesh seg_mesh_(seg_mesh);
					Optimization opt(seg_mesh_);
					std::vector<double> sub_fiedler;
					opt.dev_info_.dev_energy_ = rulings_energy;
					opt.DevelopAndSegmentation(sub_fiedler, real2img);

					//3.判断从哪里断开
					int seg_bounday_f0_id = 0; int seg_bounday_f1_id = 0;
					int ed_num = 0;
					for (auto eh : seg_mesh.edges())
					{
						if (!seg_mesh.is_boundary(eh))
						{
							OpenMesh::HalfedgeHandle heh = seg_mesh.halfedge_handle(eh, 0);
							OpenMesh::FaceHandle fh1 = seg_mesh.face_handle(heh);
							OpenMesh::HalfedgeHandle opp_heh = seg_mesh.opposite_halfedge_handle(heh);
							OpenMesh::FaceHandle fh2 = seg_mesh.face_handle(opp_heh);

							if (seg_mesh_.seg_id_[fh1.idx()] != seg_mesh_.seg_id_[fh2.idx()])
							{
								seg_bounday_f0_id = fh1.idx();
								seg_bounday_f1_id = fh2.idx();
								ed_num++;
								//break;
							}
						}
					}

					if (ed_num >= 2)
					{
						cout << "不联通！！！" << endl;
						continue;
					}
					int real_seg_bounday_f0_id = img2real[seg_bounday_f0_id]; int real_seg_bounday_f1_id = img2real[seg_bounday_f1_id];

					//
					int break_face_id = 0;
					for (int break_num = 0; break_num < bsi.path.info_every_face.size(); break_num++)
					{
						if ((bsi.path.info_every_face[break_num].face_idx == real_seg_bounday_f0_id)
							|| (bsi.path.info_every_face[break_num].face_idx == real_seg_bounday_f1_id))
						{
							break_face_id = break_num;
							break;
						}
					}

					vector<SegSurface::Path> seg_patches_temp;
					seg_one_path_to_two_path(bsi.path, break_face_id, seg_patches_temp);

					//4.拟合
					vector<GreedyWrap::BSplineSurfaceInfo> bspline_res_seg; bspline_res_seg.clear();
					for (int seg_num = 0; seg_num < seg_patches_temp.size(); seg_num++)
					{
						cout << seg_num << endl;
						BSplineNoCollision bnc(mesh);

						bnc.bspline_fitting_init_plane(seg_patches_temp[seg_num]);
						//返回拟合结果：数据能量与距离
						auto bspline_res1 = bnc.bspline_fitting_no_collision1(result_path_fitting + ".obj", mesh, patch_num, seg_num, seg_patches_temp[seg_num]);

						BSplineSurfaceInfo bsi1;
						bsi1.ctr_points = bspline_res1.spline;
						bsi1.u_knots_temp = bspline_res1.u_knot;
						bsi1.dist_energy = bspline_res1.data_e;
						bsi1.max_dist = bspline_res1.max_e;
						bsi1.mean_dist = bspline_res1.arv_e;
						bsi1.path = seg_patches_temp[seg_num];
						bsi1.wrapping_idx = wrapping_idx;
						bsi1.surface_idx = patch_num;
						bsi1.branch_idx = i;
						bspline_res_seg.push_back(bsi1);
					}

					//5.检验分割结果 
					if (bsi.dist_energy > bspline_res_seg[0].dist_energy + bspline_res_seg[1].dist_energy)
					{
						//更新拟合结果


						ofstream ioszz_("D:/programfiles/3/11/build/result/" + result_path_fitting + "/123/" + to_string(patch_num) + "/spline_seg_0.txt");
						for (int it_ = 0; it_ < bspline_res_seg[0].ctr_points.rows(); it_++)
						{
							ioszz_ << bspline_res_seg[0].ctr_points.coeffRef(it_, 0) << " " << bspline_res_seg[0].ctr_points.coeffRef(it_, 1) << " " << bspline_res_seg[0].ctr_points.coeffRef(it_, 2) << endl;
						}
						ioszz_.close();
						ofstream ioszz_1("D:/programfiles/3/11/build/result/" + result_path_fitting + "/123/" + to_string(patch_num) + "/spline_seg_1.txt");
						for (int it_ = 0; it_ < bspline_res_seg[1].ctr_points.rows(); it_++)
						{
							ioszz_1 << bspline_res_seg[1].ctr_points.coeffRef(it_, 0) << " " << bspline_res_seg[1].ctr_points.coeffRef(it_, 1) << " " << bspline_res_seg[1].ctr_points.coeffRef(it_, 2) << endl;
						}
						ioszz_1.close();
						ofstream ioszz_2("D:/programfiles/3/11/build/result/" + result_path_fitting + "/123/" + to_string(patch_num) + "/knot_seg_0.txt");
						for (int it_ = 0; it_ < bspline_res_seg[0].u_knots_temp.size(); it_++)
						{
							ioszz_2 << bspline_res_seg[0].u_knots_temp[it_] << endl;
						}
						ioszz_2.close();
						ofstream ioszz_3("D:/programfiles/3/11/build/result/" + result_path_fitting + "/123/" + to_string(patch_num) + "/knot_seg_1.txt");
						for (int it_ = 0; it_ < bspline_res_seg[1].u_knots_temp.size(); it_++)
						{
							ioszz_3 << bspline_res_seg[1].u_knots_temp[it_] << endl;
						}
						ioszz_3.close();


					}
					bspline_res_seg.clear();
				}
			}

			

			//3.评估每一块的拟合情况决定是否分割
			/*for (int patch_num = 0; patch_num < fitting_nums; patch_num++)
			{
				if (bspline_res_1temp[patch_num].mean_dist > 9e-4)
				{
					// 1.切割网格，构造bspline曲面handle
					Mesh seg_mesh; vector<int> real2img; vector<int> img2real;
					cut_mesh(mesh, seg_mesh, bspline_res_1temp[patch_num].path);
					
					vector<SegEnergy::InnerEdge> ie_temp; vector<double> rulings_energy;
					calc_seg_energy(mesh, seg_mesh, bspline_res_1temp[patch_num].path, bspline_res_1temp[patch_num].bspline_handle, bspline_res_1temp[patch_num].ctr_points, bspline_res_1temp[patch_num].u_knots_temp,
						ie_temp, rulings_energy, real2img, img2real);
					
					SegEnergy::calc_edge_energy0(seg_mesh, 1.0, rulings_energy, ie_temp);
				
					//2.Seg: 谱分割
					SegMesh seg_mesh_(seg_mesh);
					Optimization opt(seg_mesh_);
					std::vector<double> sub_fiedler;
					opt.dev_info_.dev_energy_ = rulings_energy;
					opt.DevelopAndSegmentation(sub_fiedler, real2img);

					//3.判断从哪里断开
					int seg_bounday_f0_id = 0; int seg_bounday_f1_id = 0;
					int ed_num = 0;
					for (auto eh : seg_mesh.edges())
					{
						if (!seg_mesh.is_boundary(eh))
						{
							OpenMesh::HalfedgeHandle heh = seg_mesh.halfedge_handle(eh, 0);
							OpenMesh::FaceHandle fh1 = seg_mesh.face_handle(heh);
							OpenMesh::HalfedgeHandle opp_heh = seg_mesh.opposite_halfedge_handle(heh);
							OpenMesh::FaceHandle fh2 = seg_mesh.face_handle(opp_heh);

							if (seg_mesh_.seg_id_[fh1.idx()] != seg_mesh_.seg_id_[fh2.idx()])
							{
								seg_bounday_f0_id = fh1.idx();
								seg_bounday_f1_id = fh2.idx();
								ed_num++;
								//break;
							}
						}
					}

					if (ed_num >= 2)
					{
						cout << "不联通！！！" << endl;
						continue;
					}
					int real_seg_bounday_f0_id = img2real[seg_bounday_f0_id]; int real_seg_bounday_f1_id = img2real[seg_bounday_f1_id];

					//
					int break_face_id = 0;
					for (int break_num = 0; break_num < bspline_res_1temp[patch_num].path.info_every_face.size(); break_num++)
					{
						if ((bspline_res_1temp[patch_num].path.info_every_face[break_num].face_idx == real_seg_bounday_f0_id)
							|| (bspline_res_1temp[patch_num].path.info_every_face[break_num].face_idx == real_seg_bounday_f1_id))
						{
							break_face_id = break_num;
							break;
						}
					}

					vector<SegSurface::Path> seg_patches_temp;
					seg_one_path_to_two_path(bspline_res_1temp[patch_num].path, break_face_id, seg_patches_temp);

					//4.拟合
					vector<GreedyWrap::BSplineSurfaceInfo> bspline_res_seg; bspline_res_seg.clear();
					for (int seg_num = 0; seg_num < seg_patches_temp.size(); seg_num++)
					{
						BSplineNoCollision bnc(mesh);
						bnc.bspline_fitting_init_plane(seg_patches_temp[patch_num]);

						//返回拟合结果：数据能量与距离
						auto bspline_res = bnc.bspline_fitting_no_collision("heart.obj", mesh, seg_num, seg_patches_temp[seg_num]);

						BSplineSurfaceInfo bsi;
						bsi.ctr_points = bspline_res.spline;
						bsi.u_knots_temp = bspline_res.u_knot;
						bsi.dist_energy = bspline_res.data_e;
						bsi.max_dist = bspline_res.max_e;
						bsi.mean_dist = bspline_res.arv_e;
						bsi.path = pathes_temp[patch_num];
						bsi.wrapping_idx = wrapping_idx;
						bsi.surface_idx = patch_num;
						bsi.branch_idx = i;
						bspline_res_seg.push_back(bsi);
					}

					//5.检验分割结果 
					if (bspline_res_1temp[patch_num].dist_energy > bspline_res_seg[0].dist_energy + bspline_res_seg[1].dist_energy)
					{
						//更新拟合结果


						ofstream ioszz_("D:/programfiles/3/5/build/123/" + to_string(patch_num) + "/spline_seg_0.txt");
						for (int it_ = 0; it_ < bspline_res_seg[0].ctr_points.rows(); it_++)
						{
							ioszz_ << bspline_res_seg[0].ctr_points.coeffRef(it_, 0) << " " << bspline_res_seg[0].ctr_points.coeffRef(it_, 1) << " " << bspline_res_seg[0].ctr_points.coeffRef(it_, 2) << endl;
						}
						ioszz_.close();
						ofstream ioszz_1("D:/programfiles/3/5/build/123/" + to_string(patch_num) + "/spline_seg_1.txt");
						for (int it_ = 0; it_ < bspline_res_seg[1].ctr_points.rows(); it_++)
						{
							ioszz_1 << bspline_res_seg[1].ctr_points.coeffRef(it_, 0) << " " << bspline_res_seg[1].ctr_points.coeffRef(it_, 1) << " " << bspline_res_seg[1].ctr_points.coeffRef(it_, 2) << endl;
						}
						ioszz_1.close();
						ofstream ioszz_2("D:/programfiles/3/5/build/123/" + to_string(patch_num) + "/knot_seg_0.txt");
						for (int it_ = 0; it_ < bspline_res_seg[0].u_knots_temp.size(); it_++)
						{
							ioszz_2 << bspline_res_seg[0].u_knots_temp[it_] << endl;
						}
						ioszz_2.close();
						ofstream ioszz_3("D:/programfiles/3/5/build/123/" + to_string(patch_num) + "/knot_seg_1.txt");
						for (int it_ = 0; it_ < bspline_res_seg[1].u_knots_temp.size(); it_++)
						{
							ioszz_3 << bspline_res_seg[1].u_knots_temp[it_] << endl;
						}
						ioszz_3.close();
					}
					bspline_res_seg.clear();
				}
			}*/
			bspline_res_2temp.push_back(bspline_res_1temp);
		}
		bspline_res_3temp.push_back(bspline_res_2temp);
		wrapping_idx++;
		//2.计算点距
		/*
		calc_vertex_dist_for_bspline_surfaces(bspline_res_3temp, vertex_fitting_temp, wrapping_idx);

		//3.确定覆盖区域和未覆盖区域,更新分支
		for (int v_id = 0; v_id < vertex_fitting_temp.size(); v_id++)
		{
			if (vertex_fitting_temp[v_id].vertex_bspline_dist_temp[0].second < vertex_dist_threshold)
			{
				sucess_cover_point.push_back(vertex_fitting_temp[v_id].vertex_idx);
			}
			else
			{
				failed_cover_point.push_back(vertex_fitting_temp[v_id].vertex_idx);
			}
		}

		std::sort(sucess_cover_point.begin(), sucess_cover_point.end());
		auto last0 = std::unique(sucess_cover_point.begin(), sucess_cover_point.end());
		sucess_cover_point.erase(last0, sucess_cover_point.end());

		std::sort(failed_cover_point.begin(), failed_cover_point.end());
		auto last1 = std::unique(failed_cover_point.begin(), failed_cover_point.end());
		failed_cover_point.erase(last1, failed_cover_point.end());

		connect_branches.clear();
		GraphAlgorithm gam(mesh);
		auto connect_branches_point = gam.search_connect_branchs(failed_cover_point);
		for (auto ele: connect_branches_point)
		{
			vector<int> faceId = vertex2face(ele);
			PatchSeg ps(mesh);
			ps.connect_one_branch(faceId, inaccessible_face_status);
			connect_branches.push_back(faceId);
		}
		wrapping_idx++;
	}*/
	}
}

void GreedyWrap::seg_one_path_to_two_path(SegSurface::Path& ori_path, int break_face_id, vector<SegSurface::Path>& patches_temp)
{
	patches_temp.clear();
	vector<SegSurface::OneFaceInfoInLine> patch_one(ori_path.info_every_face.begin(), ori_path.info_every_face.begin() + break_face_id + 1);
	vector<SegSurface::OneFaceInfoInLine> patch_two(ori_path.info_every_face.begin() + break_face_id + 1, ori_path.info_every_face.end());
	vector<vector<SegSurface::OneFaceInfoInLine>> seg_patches_temp = { patch_one ,patch_two };
	
	for (int i = 0; i < seg_patches_temp.size(); i++)
	{
		SegSurface::Path pa;
		pa.info_every_face = seg_patches_temp[i];
		for (int fid = 0; fid < pa.info_every_face.size(); fid++)
		{
			pa.path.push_back(pa.info_every_face[fid].v_pos);
			pa.rullings.push_back(pa.info_every_face[fid].rulings);
		}
		for (int k1 = 0; k1 < pa.rullings.size(); k1++)
		{
			for (int k2 = 0; k2 < pa.rullings[k1].size(); k2++)
			{
				for (int k3 = 0; k3 < pa.rullings[k1][k2].one_ruling_cover_face.size(); k3++)
				{
					// 排除不联通的部分
					int cover_fid = pa.rullings[k1][k2].one_ruling_cover_face[k3];
					auto it = std::find(ori_path.cover_face.begin(), ori_path.cover_face.end(), cover_fid);
					if (it != ori_path.cover_face.end())
					{
						pa.cover_face.push_back(cover_fid);
					}
				}	
			}
		}
		patches_temp.push_back(pa);
	}

	for (int i = 0; i < patches_temp.size(); i++)
	{
		std::sort(patches_temp[i].cover_face.begin(), patches_temp[i].cover_face.end());
		auto last = std::unique(patches_temp[i].cover_face.begin(), patches_temp[i].cover_face.end());
		patches_temp[i].cover_face.erase(last, patches_temp[i].cover_face.end());
		PatchSeg ps(mesh);
		ps.connect_one_branch(patches_temp[i].cover_face, inaccessible_face_status);
	}
}

void GreedyWrap::debug_write_pathes_temp(const string& file_name, vector<SegSurface::Path>& pathes_temp)
{
	/*
	1.Path
	储存在 file_name + to_string(i) + "/info_every_face/
	vector<OneFaceInfoInLine> info_every_face;
	储存在 file_name + to_string(i) + "/path.txt
	vector<Vector3d> path;
	储存在 file_name + to_string(i) + "/cover_face.txt
	vector<int> cover_face;
	//
	vector<vector<OneRuling>> rullings;
	//init.txt
	0.Vector3d center_point;
	1.Vector3d translation_direction;
	2.Vector3d path_direction;
	3.Vector3d ruling_direction;
	4.double path_length = 0;
	5.double translation_length = 0;

	2.OneFaceInfoInLine
	储存在 file_name + to_string(i) + "/info_every_face/info_every_face_" + to_string(j) + ".txt"
		0.int rank;
		1.int face_idx;
		2.int field_id;
		3.Eigen::Vector3d face_normal;
		4.Eigen::Vector3d v_pos;
		5.Eigen::Vector3d direction;
		6.Eigen::Vector3d real_dir;
		7.Eigen::Vector4d plane;
		8.Eigen::Vector3d ruling_direction;
		9.Eigen::Vector3d centeriod;
		10.set<int> cover_area_one_face;
		11.int rulings_num = 0;
		12.bool status = true;
		vector<OneRuling> rulings;
		rulings[0]
		13. p0 p1
		14. rulings[0].one_ruling_cover_face
		rulings[1]
		15. p0 p1
		16. rulings[1].one_ruling_cover_face
		.......
		.......

	struct OneRuling
	{
		// 1.p0与p1为直母线的两个端点
		Vector3d p0;Vector3d p1;
		// 2.直母线上采取的点集 计算拟合区域
		vector<Vector3d> one_rulling_sample_point_temp;
		// 3.拟合区域
		vector<int> one_ruling_cover_face;
	};

	*/
	ofstream out_info(file_name + "/info.txt");
	out_info << pathes_temp.size() << endl;
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		out_info << pathes_temp[i].info_every_face.size() << endl;
	}
	out_info.close();

	for (int i = 0; i < pathes_temp.size(); i++)
	{
		ofstream out_path(file_name + to_string(i) + "/path.txt");
		for (int j = 0; j < pathes_temp[i].path.size(); j++)
		{
			out_path << pathes_temp[i].path[j].x() << " " << pathes_temp[i].path[j].y() << " " << pathes_temp[i].path[j].z() << endl;
		}
		out_path.close();

		ofstream out_cover_face(file_name + to_string(i) + "/cover_face.txt");
		for (int j = 0; j < pathes_temp[i].cover_face.size(); j++)
		{
			out_cover_face << pathes_temp[i].cover_face[j] << endl;
		}
		out_cover_face.close();

		/*
		ofstream out_init(file_name + to_string(i) + "/init.txt");
		out_init << pathes_temp[i].center_point.x() << " " << pathes_temp[i].center_point.y() << " " << pathes_temp[i].center_point.z() << endl;
		out_init << pathes_temp[i].translation_direction.x() << " " << pathes_temp[i].translation_direction.y() << " " << pathes_temp[i].translation_direction.z() << endl;
		out_init << pathes_temp[i].path_direction.x() << " " << pathes_temp[i].path_direction.y() << " " << pathes_temp[i].path_direction.z() << endl;
		out_init << pathes_temp[i].ruling_direction.x() << " " << pathes_temp[i].ruling_direction.y() << " " << pathes_temp[i].ruling_direction.z() << endl;
		out_init << pathes_temp[i].path_length << endl;
		out_init << pathes_temp[i].translation_length << endl;
		out_init.close();
		*/

		for (int j = 0; j < pathes_temp[i].info_every_face.size(); j++)
		{
			ofstream out_info_every_face(file_name + to_string(i) + "/info_every_face/info_every_face_" + to_string(j) + ".txt");
			out_info_every_face << pathes_temp[i].info_every_face[j].rank << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].face_idx << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].field_id << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].face_normal.x() << " " << pathes_temp[i].info_every_face[j].face_normal.y() << " " << pathes_temp[i].info_every_face[j].face_normal.z() << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].v_pos.x() << " " << pathes_temp[i].info_every_face[j].v_pos.y() << " " << pathes_temp[i].info_every_face[j].v_pos.z() << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].direction.x() << " " << pathes_temp[i].info_every_face[j].direction.y() << " " << pathes_temp[i].info_every_face[j].direction.z() << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].real_dir.x() << " " << pathes_temp[i].info_every_face[j].real_dir.y() << " " << pathes_temp[i].info_every_face[j].real_dir.z() << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].plane[0] << " " << pathes_temp[i].info_every_face[j].plane[1] << " " << pathes_temp[i].info_every_face[j].plane[2] << " " << pathes_temp[i].info_every_face[j].plane[3] << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].ruling_direction.x() << " " << pathes_temp[i].info_every_face[j].ruling_direction.y() << " " << pathes_temp[i].info_every_face[j].ruling_direction.z() << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].centeriod.x() << " " << pathes_temp[i].info_every_face[j].centeriod.y() << " " << pathes_temp[i].info_every_face[j].centeriod.z() << endl;
			for (auto ele : pathes_temp[i].info_every_face[j].cover_area_one_face)
			{
				out_info_every_face << ele << " ";
			}
			out_info_every_face << endl;
			out_info_every_face << pathes_temp[i].info_every_face[j].rulings_num << endl;
			out_info_every_face << int(pathes_temp[i].info_every_face[j].status) << endl;
			for (auto ele : pathes_temp[i].info_every_face[j].rulings)
			{
				out_info_every_face << ele.p0.x() << " " << ele.p0.y() << " " << ele.p0.z() << " " << ele.p1.x() << " " << ele.p1.y() << " " << ele.p1.z() << endl;
				for (auto ele_ : ele.one_ruling_cover_face)
				{
					out_info_every_face << ele_ << " ";
				}
				out_info_every_face << endl;
			}
			out_info_every_face.close();
		}
	}
	
}


void GreedyWrap::debug_read_pathes_temp(const string& file_name, vector<SegSurface::Path>& pathes_temp)
{
	vector<vector<double>> res_first;
	read_txt(file_name + "info.txt", res_first);

	pathes_temp.resize(res_first[0][0]);
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		SegSurface::Path pa;
		vector<vector<double>> res_path;
		read_txt(file_name + to_string(i) + "/path.txt", res_path);
		for (auto ele : res_path)
		{
			Vector3d p(ele[0], ele[1], ele[2]);
			pa.path.push_back(p);
		}

		vector<vector<double>> res_cover_face;
		read_txt(file_name + to_string(i) + "/cover_face.txt", res_cover_face);
		for (auto ele : res_cover_face)
		{
			pa.cover_face.push_back(ele[0]);
		}

		vector<SegSurface::OneFaceInfoInLine> info_every_face;
		info_every_face.clear();
		//info_every_face.resize(res_first[i + 1][0]);

		//cout << res_first[i + 1][0] << endl;
		for (int j = 0; j < res_first[i + 1][0]; j++)
		{
			vector<vector<double>> res_info_every_face;
			//cout << info_every_face.size() << endl;
			read_txt(file_name + to_string(i) + "/info_every_face/info_every_face_" + to_string(j) + ".txt", res_info_every_face);
			SegSurface::OneFaceInfoInLine ofiil;
			ofiil.rank = res_info_every_face[0][0];
			ofiil.face_idx = res_info_every_face[1][0];
			ofiil.field_id = res_info_every_face[2][0];
			ofiil.face_normal = { res_info_every_face[3][0],res_info_every_face[3][1],res_info_every_face[3][2] };
			ofiil.v_pos = { res_info_every_face[4][0],res_info_every_face[4][1],res_info_every_face[4][2] };
			ofiil.direction = { res_info_every_face[5][0],res_info_every_face[5][1],res_info_every_face[5][2] };
			ofiil.real_dir = { res_info_every_face[6][0],res_info_every_face[6][1],res_info_every_face[6][2] };
			ofiil.plane = { res_info_every_face[7][0],res_info_every_face[7][1],res_info_every_face[7][2],res_info_every_face[7][3] };
			ofiil.ruling_direction = { res_info_every_face[8][0],res_info_every_face[8][1],res_info_every_face[8][2] };
			ofiil.centeriod = { res_info_every_face[9][0],res_info_every_face[9][1],res_info_every_face[9][2] };
			set<int> cover_area_one_face; cover_area_one_face.clear();
			for (auto ele : res_info_every_face[10])
			{
				cover_area_one_face.insert(ele);
			}
			ofiil.rulings_num = res_info_every_face[11][0];
			ofiil.status = res_info_every_face[12][0];

			vector<SegSurface::OneRuling> rulings; rulings.clear();
			for (int k0 = 0; k0 < (res_info_every_face.size() - 13) / 2; k0++)
			{
				SegSurface::OneRuling oru;
				oru.p0 = { res_info_every_face[13 + 2 * k0][0],res_info_every_face[13 + 2 * k0][1], res_info_every_face[13 + 2 * k0][2] };
				oru.p1 = { res_info_every_face[13 + 2 * k0][3],res_info_every_face[13 + 2 * k0][4], res_info_every_face[13 + 2 * k0][5] };
				vector<int> one_ruling_cover_face; one_ruling_cover_face.clear();
				for (auto ele : res_info_every_face[13 + 2 * k0 + 1])
				{
					one_ruling_cover_face.push_back(ele);
				}
				oru.one_ruling_cover_face = one_ruling_cover_face;
				rulings.push_back(oru);
			}
			ofiil.rulings = rulings;
			info_every_face.push_back(ofiil);
		}
		pa.info_every_face = info_every_face;
		pathes_temp[i] = pa;
		
	}
}



void GreedyWrap::greedy_wrap3()
{
	//1.初始化
	//1.1 连通分支
	vector<vector<int>> connect_branches;
	connect_branches.clear(); connect_branches.resize(1);
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		connect_branches[0].push_back(i);
	}

	//1.2 
	initialization();
	int wrapping_idx = 0;
	vector<vector<vector<GreedyWrap::BSplineSurfaceInfo>>> bspline_res_3temp; bspline_res_3temp.clear();
	// wrapping_idx ----> branch_idx ----> surface_idx

	//2.
	int iter_branchs = 0;
	while (iter_branchs < 1)
	{
		vector<vector<BSplineSurfaceInfo>>  bspline_res_2temp; bspline_res_2temp.clear();
		for (int i = 0; i < connect_branches.size(); i++)
		{
			vector<GreedyWrap::BSplineSurfaceInfo> bspline_res_1temp; bspline_res_1temp.clear();
			//1.分块
			SegSurface ss(mesh);
			vector<SegSurface::Path> pathes_temp;
			vector<int> branch_accessible = connect_branches[i];
			ss.seg_surface_no_iter(connect_branches[i], pathes_temp, inaccessible_area);

			//2.拟合
			//2.1选取大块拟合,计算块数
			branch_accessible.erase(std::remove_if(branch_accessible.begin(), branch_accessible.end(), [&](int num) {
				return std::find(inaccessible_area.begin(), inaccessible_area.end(), num) != inaccessible_area.end();
				}), branch_accessible.end());
			int fitting_nums = select_big_patches(pathes_temp, 0.95, branch_accessible.size());

			vector<SegSurface::Path> pathes_temp_copy(pathes_temp.begin(), pathes_temp.begin() + fitting_nums);
			debug_write_pathes_temp("D:/programfiles/3/11/build/result/" + result_path_fitting + "/pathes_temp/", pathes_temp_copy);

			//2.2开始拟合
			for (int patch_num = 0; patch_num < fitting_nums; patch_num++)
			{
				BSplineNoCollision bnc(mesh);
				bnc.bspline_fitting_init_plane(pathes_temp[patch_num]);

				//返回拟合结果：数据能量与距离
				auto bspline_res = bnc.bspline_fitting_no_collision("heart.obj", mesh, patch_num, pathes_temp[patch_num]);

				BSplineSurfaceInfo bsi;
				bsi.ctr_points = bspline_res.spline;
				bsi.u_knots_temp = bspline_res.u_knot;
				bsi.dist_energy = bspline_res.data_e;
				bsi.max_dist = bspline_res.max_e;
				bsi.mean_dist = bspline_res.arv_e;
				bsi.path = pathes_temp[patch_num];
				bsi.wrapping_idx = wrapping_idx;
				bsi.surface_idx = patch_num;
				bsi.branch_idx = i;
				bspline_res_1temp.push_back(bsi);
			}
		}
		iter_branchs++;
	}
}

void GreedyWrap::greedy_segment()
{
	vector<vector<int>> connect_branches;
	connect_branches.clear(); connect_branches.resize(1);
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		connect_branches[0].push_back(i);
	}

	//1.2 
	initialization();
	int wrapping_idx = 0;
	vector<vector<vector<GreedyWrap::BSplineSurfaceInfo>>> bspline_res_3temp; bspline_res_3temp.clear();
	// wrapping_idx ----> branch_idx ----> surface_idx

	//2.
	int cut_num = 0;
	vector<vector<BSplineSurfaceInfo>>  bspline_res_2temp; bspline_res_2temp.clear();
	vector<GreedyWrap::BSplineSurfaceInfo> bspline_res_1temp; bspline_res_1temp.clear();
	//1.分块
	SegSurface ss(mesh);
	vector<SegSurface::Path> pathes_temp;
	vector<int> branch_accessible = connect_branches[0];
	ss.seg_surface_no_iter(connect_branches[0], pathes_temp, inaccessible_area);
	//debug_read_pathes_temp("D:/programfiles/3/11/build/result/"+ result_path_fitting+"/pathes_temp/", pathes_temp);
	for (int i = 0; i < pathes_temp[0].info_every_face.size(); i++)
	{
		cout << pathes_temp[0].info_every_face.size() << " " << i << " " << pathes_temp[0].info_every_face[i].face_idx << endl;
	}
	//2.拟合
	//2.1选取大块拟合,计算块数
	branch_accessible.erase(std::remove_if(branch_accessible.begin(), branch_accessible.end(), [&](int num) {
		return std::find(inaccessible_area.begin(), inaccessible_area.end(), num) != inaccessible_area.end();
		}), branch_accessible.end());

	cout << "flag" << pathes_temp[0].cover_face.size() << endl;
	int fitting_nums = select_big_patches(pathes_temp, 0.95, branch_accessible.size());

	//fitting_nums = 1;
	cout << "zz: " << fitting_nums << endl;
	vector<SegSurface::Path> pathes_temp_copy(pathes_temp.begin(), pathes_temp.begin() + fitting_nums);
	debug_write_pathes_temp("D:/programfiles/3/11/build/result/" + result_path_fitting + "/pathes_temp/", pathes_temp_copy);
	pathes_temp.clear(); pathes_temp.shrink_to_fit();
	pathes_temp = pathes_temp_copy;
	ss.print_seg(pathes_temp);
	
}
