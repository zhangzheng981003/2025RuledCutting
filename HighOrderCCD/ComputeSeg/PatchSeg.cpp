#include "PatchSeg.h"

PatchSeg::PatchSeg(Mesh& mesh) :mesh(mesh)
{

}

PatchSeg::~PatchSeg()
{
}

void PatchSeg::read_txt(const string& pathname, vector<vector<double>>& res)
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

void PatchSeg::write_txt(const string& path, vector<int>& res)
{
	ofstream out(path);
	for (auto ele : res)
	{
		out << ele << endl;
	}
	out.close();
}

void PatchSeg::generate_connect_branch(int numId)
{
	//string filename0 = "D:/Zheng2024/0Code/BSpline/0314/build/123/1/filter.txt";
	//string filename1 = "D:/Zheng2024/0Code/BSpline/0314/build/123/1/faceid.txt";
	//seg_great_distortion_rulings_patch(filename0, filename1);
	string filename0 = "mesh_data/area"+to_string(numId)+".txt";
	vector<vector<double>> res_faceid;
	read_txt(filename0, res_faceid);
	vector<int> faceId; faceId.clear();
	for (auto ele : res_faceid)
	{
		faceId.push_back(ele[0]);
	}
	vector<FaceConnect> threeAdj_non_in_area(mesh.n_faces());
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		auto it_0 = std::find(faceId.begin(), faceId.end(), i);
		if (it_0 == faceId.end())
		{
			int Adj_num = 0;
			FaceConnect fc;
			fc.f_idx = i;
			for (auto ff : mesh.ff_range(mesh.face_handle(i)))
			{
				fc.Adj_faces.push_back(ff.idx());
				auto it_1 = std::find(faceId.begin(), faceId.end(), ff.idx());
				if (it_1 != faceId.end())
				{
					Adj_num++;
					fc.Adj_connect_faces.push_back(ff.idx());
				}
				else
				{
					fc.Adj_unconnect_faces.push_back(ff.idx());
					//fc.Adj_noconnect_faces
				}
			}

			fc.connnect_num = Adj_num;
			threeAdj_non_in_area[i] = fc;
		}
	}

	for (int i = 0; i < threeAdj_non_in_area.size(); i++)
	{
		if (threeAdj_non_in_area[i].connnect_num == 3)
		{
			faceId.push_back(i);
		}
	}

	//write_txt("mesh_data/area1_3.txt", faceId);
	for (int i = 0; i < threeAdj_non_in_area.size(); i++)
	{
		if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 2))
		{
			faceId.push_back(i);
		}
		else if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 1))
		{
			faceId.push_back(i);
			faceId.push_back(threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].f_idx);
		}
		else if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 0))
		{
			faceId.push_back(i);
		}
	}

	//write_txt("mesh_data/area" + to_string(numId) + ".txt", faceId);
	//检测圆环：设置为圆环则无法拟合
}

void PatchSeg::seg_great_distortion_rulings_patch(const string& filename0, const string& filename1)
{
	//1.挑取分块结果 
	vector<vector<double>> res; res.clear();
	read_txt(filename0, res);
	vector<vector<double>> res_faceid; res_faceid.clear();
	read_txt(filename1, res_faceid);
	vector<int> faceId; faceId.clear();
	for (auto ele: res_faceid)
	{
		faceId.push_back(ele[0]);
	}

	vector<vector<int>> all_face_area(mesh.n_faces());
	for (int i = 0; i < res.size(); i++)
	{
		vector<int> f_area; f_area.clear();
		for (int j = 7; j < res[i].size(); j++)
		{
			f_area.push_back(res[i][j]);
		}
		all_face_area[res[i][0]] = f_area;
	}

	//2.CutMesh：ori_face_id --> cutMesh_face_id:以及每一个面的拟合区域
	vector<int> new_patch(mesh.n_faces(), 0);
	ofstream ioszhzh("mesh_data/final_seg.txt");
	for (int i = 0; i < res.size(); i++)
	{
		new_patch[res[i][0]] = 1;
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

	vector<double> dist_temp;
	Mesh mesh_k;
	for (int i = 0; i < cut_mesh.res_two_patches.size(); i++)
	{
		if (cut_mesh.res_two_patches[i].n_faces() == res.size())
		{
			mesh_k = cut_mesh.res_two_patches[i];
			break;
		}
	}

	Tree my_tree;
	std::vector<MyTriangle> my_cgal_triangles;
	my_cgal_triangles.clear();
	for (Mesh::FaceIter f_it = mesh_k.faces_begin(); f_it != mesh_k.faces_end(); ++f_it)
	{
		vector<Point> v_p_3;
		v_p_3.clear();
		for (Mesh::FaceVertexIter fv_it = mesh_k.fv_iter(*f_it); fv_it.is_valid(); ++fv_it)
		{
			Point vertex(mesh_k.point(*fv_it)[0], mesh_k.point(*fv_it)[1], mesh_k.point(*fv_it)[2]);
			v_p_3.push_back(vertex);
		}

		auto f_normal = mesh_k.calc_face_normal(*f_it);
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

	my_tree.insert(my_cgal_triangles.begin(), my_cgal_triangles.end());
	my_tree.build();
	my_tree.accelerate_distance_queries();
	
	
	//2.1
	vector<Point> point_temp; point_temp.clear();
	vector<int> real2img; real2img.resize(mesh.n_faces(), 0);
	vector<int> img2real; img2real.resize(res.size(), 0);

	for (int i = 0; i < res.size(); i++)
	{
		Point p0 = { res[i][4]  ,res[i][5] ,res[i][6] };
		auto result = my_tree.closest_point_and_primitive(p0);
		real2img[res[i][0]] = result.second->index;
		img2real[result.second->index] = res[i][0];
	}


	vector<SegEnergy::InnerEdge> ie_temp;
	vector<vector<double>> res_u; res_u.clear();
	read_txt("E:/0RuledCutting/code/smooth/build/result/" + result_path_fitting1 + "/123/3/knot11.txt", res_u);
	vector<double> u_temp; u_temp.clear();
	for (auto ele : res_u)
	{
		u_temp.push_back(ele[0]);
	}

	vector<vector<double>> res_spline; res_spline.clear();
	read_txt("D:/Zheng2024/0Code/BSpline/0314/build/123/3/spline11.txt", res_spline);
	VectorXd BSpline; BSpline.resize(3*res_spline.size());
	for (int i = 0; i < res_spline.size(); i++)
	{
		BSpline[3 * i] = res_spline[i][0];
		BSpline[3 * i + 1] = res_spline[i][1];
		BSpline[3 * i + 2] = res_spline[i][2];
	}

	calc_edge_seg_energy00(mesh_k,BSpline, u_temp, ie_temp, img2real);
	vector<double> rullings_energy; rullings_energy.clear();
	SegEnergy::calc_edge_energy0(mesh_k, 1.0, rullings_energy, ie_temp);

	//3.Seg：谱分割
	SegMesh seg_mesh_(mesh_k);
	Optimization opt(seg_mesh_);
	std::vector<double> sub_fiedler;
	opt.dev_info_.dev_energy_ = rullings_energy;
	//opt.DevelopAndSegmentation(sub_fiedler, res, real2img);
	
	//4.获取每一块的拟合区域：cutMesh_face_id --> ori_face_id:以及每一个patch的拟合区域
	vector<int> area1;
	vector<int> area2;
	for (int i = 0; i < seg_mesh_.seg_id_.size(); i++)
	{
		if (seg_mesh_.seg_id_[i] == 0)
		{
			for (auto ele : all_face_area[img2real[i]])
			{
				auto it = std::find(faceId.begin(), faceId.end(), ele);
				if (it != faceId.end())
				{
					area1.push_back(ele);
				}
			}
		}
		else
		{
			for (auto ele : all_face_area[img2real[i]])
			{
				auto it = std::find(faceId.begin(), faceId.end(), ele);
				if (it != faceId.end())
				{
					area2.push_back(ele);
				}
			}
		}
	}

	write_txt("mesh_data/area1.txt",area1);
	write_txt("mesh_data/area2.txt", area2);
}

void PatchSeg::search_inner_hole()
{
	string filename0 = "C:/Users/ustc-gcl/Desktop/ressss.txt";
	vector<vector<double>> res_faceid; res_faceid.clear();
	read_txt(filename0, res_faceid);
	vector<int> faceId; faceId.clear();
	vector<bool> face_area_stutas(mesh.n_faces(), false);
	for (auto ele : res_faceid)
	{
		faceId.push_back(ele[0]);
		face_area_stutas[ele[0]] = true;
	}

	std::deque<int> v_q;
	vector<bool> half_edge_stutas(mesh.n_halfedges(), false);
	vector<bool> edge_stutas(mesh.n_edges(), false);
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		EH eh = mesh.edge_handle(i);
		HEH heh0 = mesh.halfedge_handle(eh, 0);
		HEH heh1 = mesh.halfedge_handle(eh, 1);
		if (face_area_stutas[mesh.face_handle(heh0).idx()] != face_area_stutas[mesh.face_handle(heh1).idx()])
		{
			v_q.push_back(i);
		}
		else if ((face_area_stutas[mesh.face_handle(heh0).idx()] == true) && (face_area_stutas[mesh.face_handle(heh1).idx()] == true))
		{
			edge_stutas[i] = true;
		}
	}

	std::deque<HEH> seam_seque;
	int e_idx = v_q.front(); v_q.pop_front();
	EH start_eh = mesh.edge_handle(e_idx);

	vector<int> add_area;
	while (true)
	{
		vector<int> Area_0; Area_0.clear();
		vector<int> seam_0; seam_0.clear();
		while (true)
		{
			HEH heh0 = mesh.halfedge_handle(start_eh, 0);
			int f0_id = mesh.face_handle(heh0).idx();
			HEH heh1 = mesh.halfedge_handle(start_eh, 1);
			int f1_id = mesh.face_handle(heh1).idx();
			if (face_area_stutas[f0_id] != face_area_stutas[f1_id])
			{
				if (!face_area_stutas[f0_id])
				{
					Area_0.push_back(f0_id);
					edge_stutas[start_eh.idx()] = true;
					for (auto eh : mesh.fe_range(mesh.face_handle(f0_id)))
					{
						if (!edge_stutas[eh.idx()])
						{
							seam_0.push_back(eh.idx());
						}
					}

				}
				else if (!face_area_stutas[f1_id])
				{
					Area_0.push_back(f1_id);
					edge_stutas[start_eh.idx()] = true;
					for (auto eh : mesh.fe_range(mesh.face_handle(f1_id)))
					{
						if (!edge_stutas[eh.idx()])
						{
							seam_0.push_back(eh.idx());
						}
					}
				}
			}

			VH vh = mesh.to_vertex_handle(heh0);
			int seam_num = 0;
			for (auto vh_to_he : mesh.voh_range(vh))
			{
				if (face_area_stutas[mesh.face_handle(vh_to_he).idx()] != face_area_stutas[mesh.face_handle(mesh.opposite_halfedge_handle(vh_to_he)).idx()])
				{
					EH e0 = mesh.edge_handle(vh_to_he);
					if (!edge_stutas[e0.idx()])
					{
						start_eh = e0;
						break;
					}
				}
				else
				{
					seam_num++;
				}

			}

			if (seam_num == mesh.valence(vh))
			{
				break;
			}
		}

		for (int i = 0; i < Area_0.size(); i++)
		{
			face_area_stutas[Area_0[i]] = true;
			add_area.push_back(Area_0[i]);
		}

		if (seam_0.size() != 0)
		{
			start_eh = mesh.edge_handle(seam_0[0]);
		}
		else
		{
			break;
		}
	}

	write_txt("add_area.txt", add_area);
}


void PatchSeg::calc_edge_seg_energy00(Mesh& seg_mesh, VectorXd& BSpline, vector<double>& u_temp, vector<SegEnergy::InnerEdge>& ie_temp, vector<int>& img2real)
{
	
	BSplineFitting sf(seg_mesh);
	int curve_ctr_num = BSpline.size()/6;
	vector<BSplineFitting::CtrInfo> ctr_point1; vector<BSplineFitting::CtrInfo> ctr_point2;
	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = BSpline[3 * j];
		Standard_Real py = BSpline[3 * j + 1];
		Standard_Real pz = BSpline[3 * j + 2];
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point1.push_back(ci);
	}

	for (int j = 0; j < curve_ctr_num; j++)
	{
		BSplineFitting::CtrInfo ci;
		Standard_Real px = BSpline[3 * (j + curve_ctr_num)];
		Standard_Real py = BSpline[3 * (j + curve_ctr_num) + 1];
		Standard_Real pz = BSpline[3 * (j + curve_ctr_num) + 2];
		ci.position = { px,py,pz };
		ci.weight = 1;
		ctr_point2.push_back(ci);
	}

	sf.u_knot_temp.clear(); sf.u_knot_temp = u_temp;
	Handle(Geom_BSplineSurface) BSpline_surface = nullptr;
	sf.create_BSpline_surface(ctr_point1, ctr_point2, BSpline_surface);

	std::vector<vector<Vector3d>> P;
	//int curve_ctr_num = BSpline.size() / 6;
	P.resize(curve_ctr_num);
	for (int i = 0; i < P.size(); i++)
	{
		P[i].resize(2);
	}

	for (int i = 0; i < P.size(); i++)
	{
		P[i][0][0] = BSpline[3 * i];
		P[i][0][1] = BSpline[3 * i + 1];
		P[i][0][2] = BSpline[3 * i + 2];
		P[i][1][0] = BSpline[3 * curve_ctr_num + 3 * i];
		P[i][1][1] = BSpline[3 * curve_ctr_num + 3 * i + 1];
		P[i][1][2] = BSpline[3 * curve_ctr_num + 3 * i + 2];
	}

	vector<vector<double>> res; res.clear();
	read_txt("D:/Zheng2024/0Code/BSpline/0314/build/123/3/filter.txt", res);
	vector<int> real_in_txt; real_in_txt.clear(); real_in_txt.resize(mesh.n_faces());
	for (int i = 0; i < res.size(); i++)
	{
		real_in_txt[res[i][0]] = i;
	}

	vector<SegEnergy::AdjFace> adj_face_temp; adj_face_temp.clear(); adj_face_temp.resize(seg_mesh.n_faces());
	for (int i = 0; i < seg_mesh.n_faces(); i++)
	{
		SegEnergy::AdjFace af;
		af.face_idx = i;
		af.rulling_f0_ori = { res[real_in_txt[img2real[i]]][1],res[real_in_txt[img2real[i]]][2],res[real_in_txt[img2real[i]]][3] };
		FH fh = seg_mesh.face_handle(i);
		Vec3d p_ = seg_mesh.calc_face_centroid(fh);
		gp_Pnt pointToProject(p_[0], p_[1], p_[2]);
		GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, BSpline_surface);
		if (projectPoint.NbPoints() > 0)
		{
			Standard_Real u, v;
			projectPoint.LowerDistanceParameters(u, v);
			Standard_Real dist;
			dist = projectPoint.LowerDistance();
			int u_interval_order = sf.u_k_position(u, sf.u_knot_temp);

			vector<double> u_basis_temp;
			if (u_interval_order == 0)
			{
				vector<double> basis_temp(sf.poly_degree + 1, 0);
				basis_temp[0] = 1;
				u_basis_temp = basis_temp;
			}
			else if (u_interval_order == sf.u_knot_temp.size() - 1)
			{
				vector<double> basis_temp(sf.poly_degree + 1, 0);
				basis_temp[sf.poly_degree] = 1;
				u_basis_temp = basis_temp;
			}
			else
			{
				sf.calc_basis_fuction(u, u_interval_order, sf.poly_degree, u_basis_temp);
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
			else if (u_interval_order == sf.u_knot_temp.size() - 1)
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
					pi0 = P[u_interval_order - sf.poly_degree + i][0];
					pi1 = P[u_interval_order - sf.poly_degree + i][1];
					rulling += u_basis_temp[i] * (pi1 - pi0);
				}
			}
			af.rulling_f0_bspline = rulling;
		}
		adj_face_temp[i] = af;
	}

	ie_temp.clear(); //ie_temp.resize(mesh.n_edges());
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
				GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, BSpline_surface);
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
		else
		{
			//SegEnergy::InnerEdge ie;
			//ie_temp.push_back(ie);
		}
	}
}

void PatchSeg::run()
{
	string filename0 = "d:/zheng2024/0code/bspline/0314/build/123/3/filter.txt";
	string filename1 = "d:/zheng2024/0code/bspline/0314/build/123/3/faceid.txt";
	seg_great_distortion_rulings_patch(filename0, filename1);
	generate_connect_branch(1);
	generate_connect_branch(2);
	
}


void PatchSeg::connect_one_branch(vector<int>& faceId, vector<bool>& inaccessible_face_status)
{
	vector<FaceConnect> threeAdj_non_in_area(mesh.n_faces());
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		auto it_0 = std::find(faceId.begin(), faceId.end(), i);
		if (it_0 == faceId.end())
		{
			int Adj_num = 0;
			FaceConnect fc;
			fc.f_idx = i;
			for (auto ff : mesh.ff_range(mesh.face_handle(i)))
			{
				fc.Adj_faces.push_back(ff.idx());
				auto it_1 = std::find(faceId.begin(), faceId.end(), ff.idx());
				if (it_1 != faceId.end())
				{
					Adj_num++;
					fc.Adj_connect_faces.push_back(ff.idx());
				}
				else
				{
					fc.Adj_unconnect_faces.push_back(ff.idx());
					//fc.Adj_noconnect_faces
				}
			}

			fc.connnect_num = Adj_num;
			threeAdj_non_in_area[i] = fc;
		}
	}

	for (int i = 0; i < threeAdj_non_in_area.size(); i++)
	{
		if ((threeAdj_non_in_area[i].connnect_num == 3) && (inaccessible_face_status[i]))
		{
			faceId.push_back(i);
		}
	}

	for (int i = 0; i < threeAdj_non_in_area.size(); i++)
	{
		if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 2) && (inaccessible_face_status[i]))
		{
			faceId.push_back(i);
		}
		else if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 1) && (inaccessible_face_status[i]))
		{
			faceId.push_back(i);
			faceId.push_back(threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].f_idx);
		}
		else if ((threeAdj_non_in_area[i].connnect_num == 2) && (threeAdj_non_in_area[threeAdj_non_in_area[i].Adj_unconnect_faces[0]].connnect_num == 0) && (inaccessible_face_status[i]))
		{
			faceId.push_back(i);
		}
	}

}