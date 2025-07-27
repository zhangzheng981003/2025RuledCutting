#include "SegSurface.h"

SegSurface::SegSurface(Mesh& mesh) :mesh(mesh)
{

}

SegSurface::~SegSurface()
{
	if (lg)
	{
		delete lg;
		lg = NULL;
	}
}

void SegSurface::collect_init_informarion()
{
	//1.生成主曲率场
	lg = new LoopGen::LoopGen(mesh);
	std::string name = "123";
	lg->SetModelName(name);
	lg->InitializeField();
	crossfield = lg->cf->getCrossField();
	matching = lg->cf->getMatching();
	point_singularity = lg->cf->getSingularity();

	//2.对每个面生成基础信息
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		FaceBasisInfo fhi;
		FH fh = mesh.face_handle(i);
		fhi.face_area = mesh.calc_face_area(fh);
		Mesh::Point fc = mesh.calc_face_centroid(fh);
		Mesh::Point fn = mesh.calc_face_normal(fh);
		fhi.face_normal = { fn[0], fn[1], fn[2] }; fhi.face_normal.normalize();
		fhi.centeriod = { fc[0], fc[1], fc[2] };
		fhi.face_idx = i;
		face_basis_info_temp.push_back(fhi);
	}

	//3.初始奇异面
	face_singlirity.clear();
	for (int i = 0; i < point_singularity.size(); i++)
	{
		for (int j = 0; j < point_singularity[i].size(); j++)
		{
			int v_id = point_singularity[i][j];
			auto vh = mesh.vertex_handle(v_id);
			for (Mesh::VertexFaceIter vf_it = mesh.vf_begin(vh); vf_it != mesh.vf_end(vh); ++vf_it)
			{
				face_singlirity.push_back(vf_it->idx());
			}
		}
	}
	std::sort(face_singlirity.begin(), face_singlirity.end());
	auto last = std::unique(face_singlirity.begin(), face_singlirity.end());
	face_singlirity.erase(last, face_singlirity.end());
	
	//4.构造AABB tree
	for (Mesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<Point> v_p_3;
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
	ruling_length = (ptMax - ptMin).norm();

	//5.


}

vector<int> SegSurface::sample_target_face(int target_num, vector<int>& available_areas)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<int> intDistribution(0, available_areas.size() - 1); // 定义整数范围
	std::vector<int> randomIntegers;
	for (int i = 0; i < target_num; ++i)
	{
		int randomInt = intDistribution(gen); // 生成随机整数
		randomIntegers.push_back(randomInt);
	}
	return randomIntegers;
}


vector<vector<SegSurface::OneFaceInfoInLine>> SegSurface::two_principal_curvatures_lines_on_one_face(int one_face_idx, int line_length,vector<int>& face_stop)
{
	vector<vector<OneFaceInfoInLine>> two_path_one_face_temp;
	two_path_one_face_temp.resize(4);
	FH fh = mesh.face_handle(one_face_idx);
	Vector3d cp = face_basis_info_temp[one_face_idx].centeriod;
	Vector3d fn = face_basis_info_temp[one_face_idx].face_normal;

	//1.计算起点四个点和面的信息
	for (int loop_num = 0; loop_num < 2; loop_num++)
	{
		//1.从一个面的中心点开始搜索
		Vector3d start_point = cp + crossfield.col(4 * one_face_idx + loop_num);
		
		//2.根据面法向与曲率方向计算这个平面的法相
		Vector3d plane_normal = crossfield.col(4 * one_face_idx + loop_num).cross(fn); 
		plane_normal.normalize();
		
		//3.平面的点法(a,b,c)式方程: ax+by+cz+d=0 
		double d_value = -plane_normal.transpose() * cp;
		//[a,b,c,d]*[x,y,z,1]T=0
		Vector4d plane_expression = { plane_normal[0],plane_normal[1],plane_normal[2],d_value };
		
		//4.计算平面在那条边上，有两个边
		//4.1
		vector<VH> vh_temp; vh_temp.clear();
		for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it.is_valid(); ++fv_it)
		{
			VH vhh = mesh.vertex_handle(fv_it->idx());
			vh_temp.push_back(vhh);
		}
		vector<vector<VH>> try_temp = { {vh_temp[0],vh_temp[1]},{vh_temp[0],vh_temp[2]}, {vh_temp[1],vh_temp[2]} };
		vector<vector<VH>> valid_temp; valid_temp.clear();
		//三条边种选两条
		for (int i = 0; i < 3; i++)
		{
			Mesh::Point p0 = mesh.point(try_temp[i][0]);
			Mesh::Point p1 = mesh.point(try_temp[i][1]);

			Vector4d vh0(p0[0], p0[1], p0[2], 1);
			Vector4d vh1(p1[0], p1[1], p1[2], 1);

			double res1 = (vh0.transpose() * plane_expression);
			double res2 = (vh1.transpose() * plane_expression);
			if (res1 * res2 < 0)
			{
				valid_temp.push_back(try_temp[i]);
			}
		}

		//4.2 两条边代表两个方向
		for (int i = 0; i < 2; i++)
		{
			Mesh::Point v0(mesh.point(valid_temp[i][0]));
			Mesh::Point v1(mesh.point(valid_temp[i][1]));
			Vector3d segmentStart = { v0[0], v0[1] ,v0[2] }; Vector4d segmentStart4d = { v0[0], v0[1] ,v0[2],1 };
			Vector3d segmentEnd = { v1[0], v1[1] ,v1[2] };
			Vector3d segmentDirection = { (v1 - v0)[0],(v1 - v0)[1] ,(v1 - v0)[2] };
			double numerator = plane_expression.transpose() * segmentStart4d;
			double denominator = plane_expression[0] * segmentDirection[0] + plane_expression[1] * segmentDirection[1] + plane_expression[2] * segmentDirection[2];
			double intersectionParam = -numerator / denominator;
			HEH eh = mesh.find_halfedge(valid_temp[i][0], valid_temp[i][1]);
			if (mesh.face_handle(eh).idx() == one_face_idx)
			{
				//
			}
			else
			{
				eh = mesh.opposite_halfedge_handle(eh);
			}

			int field_id = 0;
			OneFaceInfoInLine sp;
			sp.rank = 0;
			sp.face_idx = one_face_idx;
			sp.face_normal = fn;
			sp.v_pos = segmentStart + intersectionParam * segmentDirection;
			sp.heh = eh;
			double value = ((sp.v_pos - cp).normalized()).dot(crossfield.col(4 * one_face_idx + loop_num));
			value = std::round(value);
			if (value == 1)
			{
				field_id = loop_num;
			}
			else
			{
				field_id = loop_num + 2;
			}

			sp.field_id = field_id;
			sp.plane = plane_expression;
			sp.real_dir = (sp.v_pos - cp).normalized();
			sp.direction = crossfield.col(4 * sp.face_idx + field_id);

			sp.ruling_direction = sp.face_normal.cross(sp.direction);
			sp.centeriod = face_basis_info_temp[sp.face_idx].centeriod;

			//前两个为一条，后两个为一条
			two_path_one_face_temp[2 * loop_num + i].push_back(sp);
		}
	}

	////2.计算路途的点和面的信息
	int target_num = std::round(line_length / 2);
	int search_num = 0;
	vector<int> dir_num_temp = { 0,1,2,3 };
	while ((search_num < target_num - 1) && (dir_num_temp.size() > 0))
	{
		vector<int> delete_ele; delete_ele.clear();
		for (int k = 0; k < dir_num_temp.size(); k++)
		{
			int i = dir_num_temp[k];
			OneFaceInfoInLine last_sp = two_path_one_face_temp[i][two_path_one_face_temp[i].size() - 1];
			HEH next_eh = mesh.opposite_halfedge_handle(last_sp.heh);
			int next_face_id = mesh.face_handle(next_eh).idx();
			Vector3d field_dir = crossfield.col(4 * next_face_id + ((matching[last_sp.heh.idx()] + last_sp.field_id) % 4));
			Vector3d face_normal = face_basis_info_temp[next_face_id].face_normal;

			Vector3d pn = face_normal.cross(field_dir); pn.normalize();
			double dv = -pn.transpose() * last_sp.v_pos;
			Vector4d plane_ = { pn[0],pn[1],pn[2],dv };

			Mesh::VertexHandle check_v;
			Mesh::VertexHandle vf = mesh.from_vertex_handle(last_sp.heh);
			int vf_idx = vf.idx();
			Mesh::VertexHandle vt = mesh.to_vertex_handle(last_sp.heh);
			int vt_idx = vt.idx();
			for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(mesh.face_handle(next_eh)); fv_it.is_valid(); ++fv_it)
			{
				if (fv_it->idx() != vf_idx && fv_it->idx() != vt_idx)
				{
					check_v = mesh.vertex_handle(fv_it->idx());
				}
			}
			vector<vector<Mesh::VertexHandle>> check_edge = { {vf,check_v}, {vt,check_v} };
			vector<Mesh::VertexHandle> interection_edge; interection_edge.clear();
			for (int check_num = 0; check_num < check_edge.size(); check_num++)
			{
				Mesh::Point p0 = mesh.point(check_edge[check_num][0]);
				Mesh::Point p1 = mesh.point(check_edge[check_num][1]);

				Vector4d vh0(p0[0], p0[1], p0[2], 1);
				Vector4d vh1(p1[0], p1[1], p1[2], 1);

				double res1 = (vh0.transpose() * plane_);
				double res2 = (vh1.transpose() * plane_);
				if (res1 * res2 < 0)
				{
					interection_edge = check_edge[check_num];
				}
			}

			Mesh::Point vc0(mesh.point(interection_edge[0]));
			Mesh::Point vc1(mesh.point(interection_edge[1]));
			Vector3d segmentStart = { vc0[0], vc0[1] ,vc0[2] }; Vector4d segmentStart4d = { vc0[0], vc0[1] ,vc0[2],1 };
			Vector3d segmentEnd = { vc1[0], vc1[1] ,vc1[2] };
			Vector3d segmentDirection = { (vc1 - vc0)[0],(vc1 - vc0)[1] ,(vc1 - vc0)[2] };
			double numerator = plane_.transpose() * segmentStart4d;
			double denominator = plane_[0] * segmentDirection[0] + plane_[1] * segmentDirection[1] + plane_[2] * segmentDirection[2];
			double intersectionParam = -numerator / denominator;

			HEH eh = mesh.find_halfedge(interection_edge[0], interection_edge[1]);
			if (mesh.face_handle(eh).idx() == next_face_id)
			{
				//
			}
			else
			{
				eh = mesh.opposite_halfedge_handle(eh);
			}

			OneFaceInfoInLine next_sp;
			next_sp.rank = two_path_one_face_temp[i].size();
			next_sp.face_idx = next_face_id;
			next_sp.face_normal = face_normal;
			next_sp.v_pos = segmentStart + intersectionParam * segmentDirection;
			next_sp.heh = eh;
			next_sp.field_id = (matching[last_sp.heh.idx()] + last_sp.field_id) % 4;
			next_sp.plane = plane_;
			next_sp.direction = field_dir;
			next_sp.real_dir = (next_sp.v_pos - last_sp.v_pos).normalized();
			
			next_sp.ruling_direction = next_sp.face_normal.cross(next_sp.direction);
			next_sp.centeriod = face_basis_info_temp[next_sp.face_idx].centeriod;

			double vector_angle = next_sp.real_dir.dot(last_sp.real_dir);
			auto it = std::find(face_stop.begin(), face_stop.end(), next_sp.face_idx);
			if (it == face_stop.end() && vector_angle > line_angle_threshold)
			{
				two_path_one_face_temp[i].push_back(next_sp);
			}
			else
			{
				delete_ele.push_back(i);
			}
		}

		for (int k = 0; k < delete_ele.size(); k++)
		{
			int valueToRemove = delete_ele[k];
			auto newEnd = std::remove(dir_num_temp.begin(), dir_num_temp.end(), valueToRemove);
			dir_num_temp.erase(newEnd, dir_num_temp.end());
		}
		delete_ele.clear();
		search_num++;
	}

	return two_path_one_face_temp;
}

void SegSurface::sample_principal_curvatures_lines_on_many_faces(int target_num, vector<int>& many_faces_idx_temp, int line_length, vector<SampleTwoPathesOnOneFace>& sample_lines_temp)
{
	//1.只计算一次
	
	vector<int> face_stop = calc_face_stop_temp(many_faces_idx_temp);

	/*
	if (target_num> many_faces_idx_temp.size())
	{
		target_num = many_faces_idx_temp.size();
	}
	*/

	vector<int> sample_res = sample_target_face(target_num, many_faces_idx_temp);

	// 每一个被采到的面有两条
	for (int i = 0; i < sample_res.size(); i++)
	{
		// 四条合成有序的两条
		SampleTwoPathesOnOneFace spath;
		spath.start_face_idx = many_faces_idx_temp[sample_res[i]];
		//cout << " sss: " << spath.start_face_idx << endl;
		spath.four_lines_temp = two_principal_curvatures_lines_on_one_face(spath.start_face_idx, line_length, face_stop);
		vector<vector<OneFaceInfoInLine>> spvv; spvv.clear();
		for (int k = 0; k < spath.four_lines_temp.size() / 2; k++)
		{
			vector<OneFaceInfoInLine> spv; spv.clear();
			for (int k2 = 0; k2 < spath.four_lines_temp[2 * k].size(); k2++)
			{
				//spath.four_lines_temp[2 * k][spath.four_lines_temp[2 * k].size() - 1 - k2].rank = k2;
				spv.push_back(spath.four_lines_temp[2 * k][spath.four_lines_temp[2 * k].size() - 1 - k2]);
			}

			auto mid_res = spath.four_lines_temp[2 * k + 1];
			for (int k2 = 0; k2 < spath.four_lines_temp[2 * k + 1].size() - 1; k2++)
			{
				spath.four_lines_temp[2 * k + 1][k2 + 1].v_pos = mid_res[k2].v_pos;
			}

			for (int k2 = 1; k2 < spath.four_lines_temp[2 * k + 1].size(); k2++)
			{
				spv.push_back(spath.four_lines_temp[2 * k + 1][k2]);
			}
			spvv.push_back(spv);
		}
		spath.two_lines_temp = spvv;
		sample_lines_temp.push_back(spath);
	}
}

void SegSurface::generate_rullings_on_one_face(int ruling_gap, double high_value, int rulings_sample_num, OneFaceInfoInLine& ofiil)
{
	
}

void SegSurface::generate_rullings_on_many_faces(int ruling_gap, double high_value, int rulings_sample_num, vector<SampleTwoPathesOnOneFace>& sample_lines_temp)
{
	double gap = average_edge_length / double(ruling_gap);
	for (int i = 0; i < sample_lines_temp.size(); i++)
	{
		vector<vector<OneFaceInfoInLine>> two_lines_temp_on_one_face = sample_lines_temp[i].two_lines_temp;
		for (int k = 0; k < two_lines_temp_on_one_face.size(); k++)
		{
			for (int k2 = 0; k2 < two_lines_temp_on_one_face[k].size() - 1; k2++)
			{
				//1.计算一个面上的主曲率线的方向
				Vector3d line_direction = two_lines_temp_on_one_face[k][k2 + 1].v_pos - two_lines_temp_on_one_face[k][k2].v_pos;
				double dir_length = line_direction.norm();
				int rulling_num = double(dir_length) / double(gap);
				line_direction.normalize();

				//2.初始位置
				Vector3d vp = two_lines_temp_on_one_face[k][k2].v_pos + two_lines_temp_on_one_face[k][k2].face_normal * high_value;
				Vector3d rl_dir = two_lines_temp_on_one_face[k][k2].face_normal.cross(two_lines_temp_on_one_face[k][k2].direction); 
				rl_dir.normalize();

				//3.只能采一条
				if (rulling_num == 0)
				{
					two_lines_temp_on_one_face[k][k2].rulings_num = 1;
					two_lines_temp_on_one_face[k][k2].rulings.clear();
					OneRuling rl;
					rl.p0 = vp + rl_dir * 0.5 * ruling_length;
					rl.p1 = vp - rl_dir * 0.5 * ruling_length;
					double s_gap = double(ruling_length) / double(rulings_sample_num);
					rl.one_rulling_sample_point_temp.clear();
					for (int k3 = 0; k3 < rulings_sample_num; k3++)
					{
						Vector3d v0_sample = rl.p0 - k3 * s_gap * rl_dir;
						rl.one_rulling_sample_point_temp.push_back(v0_sample);
					}
					two_lines_temp_on_one_face[k][k2].rulings.push_back(rl);
				}
				else
				{
					two_lines_temp_on_one_face[k][k2].rulings_num = rulling_num + 1;
					two_lines_temp_on_one_face[k][k2].rulings.clear();
					double s_gap = double(ruling_length) / double(rulings_sample_num);
					for (int k3 = 0; k3 < rulling_num + 1; k3++)
					{
						OneRuling rl;
						Vector3d vp_ = vp + k3 * gap * line_direction;
						rl.p0 = vp_ + rl_dir * 0.5 * ruling_length;
						rl.p1 = vp_ - rl_dir * 0.5 * ruling_length;
						rl.one_rulling_sample_point_temp.clear();
						for (int k4 = 0; k4 < rulings_sample_num; k4++)
						{
							Vector3d v0_sample = rl.p0 - k4 * s_gap * rl_dir;
							rl.one_rulling_sample_point_temp.push_back(v0_sample);
						}
						two_lines_temp_on_one_face[k][k2].rulings.push_back(rl);
					}
				}
			}

			if (two_lines_temp_on_one_face[k].size() != 0)
			{
				two_lines_temp_on_one_face[k].pop_back();
			}

		}
		sample_lines_temp[i].two_lines_temp = two_lines_temp_on_one_face;
	}
}

void SegSurface::calc_fitting_area_on_many_faces(double epsilon, vector<SampleTwoPathesOnOneFace>& sample_lines_temp, vector<Path>& pathes_temp)
{
	for (int i = 0; i < sample_lines_temp.size(); i++)
	{
		for (int k = 0; k < sample_lines_temp[i].two_lines_temp.size(); k++)
		{
			for (int k1 = 0; k1 < sample_lines_temp[i].two_lines_temp[k].size(); k1++)
			{
				//1.判断相交
				vector<OneRuling> rulings_one_face_candicate = sample_lines_temp[i].two_lines_temp[k][k1].rulings;
				Vector3d ray_dir = -sample_lines_temp[i].two_lines_temp[k][k1].face_normal;
				int intersection_num = 0;
				for (int k2 = 0; k2 < rulings_one_face_candicate.size(); k2++)
				{
					Point rl_p0(rulings_one_face_candicate[k2].p0.x(), rulings_one_face_candicate[k2].p0.y(), rulings_one_face_candicate[k2].p0.z());
					Point rl_p1(rulings_one_face_candicate[k2].p1.x(), rulings_one_face_candicate[k2].p1.y(), rulings_one_face_candicate[k2].p1.z());
					Segment seg_line(rl_p0, rl_p1);
					if (tree.do_intersect(seg_line))
					{
						intersection_num++;
						break;
					}
				}

				if (intersection_num == 0)
				{
					//2.不相交:计算面积
					for (int k2 = 0; k2 < rulings_one_face_candicate.size(); k2++)
					{
						for (int k3 = 0; k3 < rulings_one_face_candicate[k2].one_rulling_sample_point_temp.size(); k3++)
						{
							//判断点与网格相交
							Vector3d p1 = rulings_one_face_candicate[k2].one_rulling_sample_point_temp[k3];
							Vector3d p2 = p1 + ray_dir;
							int face_id = spline_ray_intersection_(p1, p2, epsilon);
							if (face_id != -1)
							{
								rulings_one_face_candicate[k2].one_ruling_cover_face.push_back(face_id);
							}
						}
					}
				}
				else
				{
					sample_lines_temp[i].two_lines_temp[k][k1].status = false;
				}
				sample_lines_temp[i].two_lines_temp[k][k1].rulings.clear();
				sample_lines_temp[i].two_lines_temp[k][k1].rulings = rulings_one_face_candicate;
			}
		}
	}
	
	pathes_temp.clear();
	for (int i = 0; i < sample_lines_temp.size(); i++)
	{
		for (int k = 0; k < sample_lines_temp[i].two_lines_temp.size(); k++)
		{
			int flag_gap = 0;
			int false_num = 0;
		flag:
			Path pa; pa.cover_face.clear(); pa.path.clear();
			for (int k1 = flag_gap; k1 < sample_lines_temp[i].two_lines_temp[k].size(); k1++)
			{
				if (sample_lines_temp[i].two_lines_temp[k][k1].status)
				{
					pa.path.push_back(sample_lines_temp[i].two_lines_temp[k][k1].v_pos);
					for (int zhzh = 0; zhzh < sample_lines_temp[i].two_lines_temp[k][k1].rulings.size(); zhzh++)
					{
						for (auto zhzh0 : sample_lines_temp[i].two_lines_temp[k][k1].rulings[zhzh].one_ruling_cover_face)
						{
							sample_lines_temp[i].two_lines_temp[k][k1].cover_area_one_face.insert(zhzh0);
						}
					}
					pa.info_every_face.push_back(sample_lines_temp[i].two_lines_temp[k][k1]);
					
					//----0315 新的信息
					/*
					FilterInfo fii;
					fii.face_idx = sample_path_temp[i].two_smooth_path_temp[k][k1].face_idx;
					fii.rule_line_dir = sample_path_temp[i].two_smooth_path_temp[k][k1].face_normal.cross(sample_path_temp[i].two_smooth_path_temp[k][k1].direction);
					fii.rule_line_dir.normalize();
					fii.face_center = mesh.calc_face_centroid(mesh.face_handle(fii.face_idx));
					for (int zhzh = 0; zhzh < sample_path_temp[i].two_smooth_path_temp[k][k1].rulings.size(); zhzh++)
					{
						for (auto zhzh0 : sample_path_temp[i].two_smooth_path_temp[k][k1].rulings[zhzh].cover_face)
						{
							fii.cover_area_one_face.insert(zhzh0);
						}
					}
					pa.filter_info_temp.push_back(fii);
					*/
					//----0315 新的信息

					pa.rullings.push_back(sample_lines_temp[i].two_lines_temp[k][k1].rulings);
					for (const auto& ele : sample_lines_temp[i].two_lines_temp[k][k1].rulings)
					{
						for (int k2 = 0; k2 < ele.one_ruling_cover_face.size(); k2++)
						{
							pa.cover_face.push_back(ele.one_ruling_cover_face[k2]);
						}
					}
				}
				else
				{
					false_num++;
					if (pa.path.size() != 0)
					{
						pathes_temp.push_back(pa);
					}

					flag_gap = flag_gap + 1;
					if (flag_gap >= sample_lines_temp[i].two_lines_temp[k].size())
					{
						break;
					}
					else
					{
						goto flag;
					}

				}
			}

			if (false_num == 0)
			{
				pathes_temp.push_back(pa);
			}
		}
	}

	for (int i = 0; i < pathes_temp.size(); i++)
	{
		std::sort(pathes_temp[i].cover_face.begin(), pathes_temp[i].cover_face.end());
		auto last = std::unique(pathes_temp[i].cover_face.begin(), pathes_temp[i].cover_face.end());
		pathes_temp[i].cover_face.erase(last, pathes_temp[i].cover_face.end());
	}
}

//!!!!!
int SegSurface::spline_ray_intersection_(Eigen::Vector3d start_point, Eigen::Vector3d end_point, double epsilon)
{
	Point ray_source(start_point.x(), start_point.y(), start_point.z());
	Point ray_direction(end_point.x(), end_point.y(), end_point.z());
	Ray ray_query(ray_source, ray_direction);
	Ray_intersection intersection = tree.first_intersection(ray_query);

	if (intersection)
	{
		if (boost::get<Point>(&(intersection->first)))
		{
			const Point* p = boost::get<Point>(&(intersection->first));
			Point p_close(p->x(), p->y(), p->z());
			double distance_pro = std::sqrt(CGAL::squared_distance(p_close, ray_source));
			if (distance_pro < epsilon)
			{
				auto result = my_tree.closest_point_and_primitive(p_close);
				return result.second->index;
			}
			else
			{
				return -1;
			}

		}
		else
		{
			return -1;
		}
	}
	else
	{
		return -1;
	}
}

void SegSurface::generate_connect_branch(vector<int>& faceId)
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
		if (threeAdj_non_in_area[i].connnect_num == 3)
		{
			faceId.push_back(i);
		}
	}

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
}


void SegSurface::set_cover_algorithm(vector<Path>& pathes_temp)
{
	/*
	ofstream out("path_file.txt");
	GraphAlgorithm gam(mesh);
	gam.delete_small_patches(pathes_temp[2].cover_face);
	generate_connect_branch(pathes_temp[2].cover_face);
	for (int i = 0; i < pathes_temp[2].cover_face.size(); i++)
	{
		out << pathes_temp[2].cover_face[i] << endl;
	}
	out.close();*/
	
	/*
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		string path_file = "D:/Zheng2024/0Code/BSpline/0314/build3/123/" + to_string(i) + "//faceid.txt";
		ofstream out(path_file);
		for (int j = 0; j < pathes_temp[i].cover_face.size(); j++)
		{
			out << pathes_temp[i].cover_face[j] << endl;
		}
		out.close();
	}*/

	/*
	vector<Path> pathes_temp_copy; pathes_temp_copy.clear();
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		if ((pathes_temp[i].cover_face.size() >= candidate_set_cover_min_num) && (pathes_temp[i].rullings.size() > candidate_u_face_min_num))
		{
			GraphAlgorithm gam(mesh);
			gam.delete_small_patches(pathes_temp[i].cover_face);
			generate_connect_branch(pathes_temp[i].cover_face);
			if (pathes_temp[i].cover_face.size() >= candidate_set_cover_min_num)
			{
				pathes_temp_copy.push_back(pathes_temp[i]);
			}
		}
	}
	pathes_temp = pathes_temp_copy;
	*/

	vector<Path> pathes_temp_copy; pathes_temp_copy.clear();
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		if (pathes_temp[i].cover_face.size() > 0)
		{
			pathes_temp_copy.push_back(pathes_temp[i]);
		}
	}
	pathes_temp.clear(); pathes_temp.shrink_to_fit();
	pathes_temp = pathes_temp_copy;
	pathes_temp_copy.clear(); pathes_temp_copy.shrink_to_fit();



	cout << "1" << endl;
	ofstream ios_("set_cover.txt");
	ios_ << mesh.n_faces() << " " << pathes_temp.size() << endl;
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		GraphAlgorithm gam(mesh);
		gam.delete_small_patches(pathes_temp[i].cover_face);
		//generate_connect_branch(pathes_temp[i].cover_face);
		for (int j = 0; j < pathes_temp[i].cover_face.size(); j++)
		{
			ios_ << pathes_temp[i].cover_face[j] + 1 << " ";
		}
		ios_ << std::endl;
	}
	ios_.close();
	SetCoverRun::set_cover_run();


	vector<vector<double>> res;
	read_txt("set_cover_result.txt", res);

	vector<Path> pathes_temp_select; pathes_temp_select.clear();
	for (int i = 0; i < res.size(); i++)
	{
		pathes_temp_select.push_back(pathes_temp[res[i][0]]);
	}

	pathes_temp.clear();
	
	pathes_temp = pathes_temp_select;
}

void SegSurface::read_txt(const string& pathname, vector<vector<double>>& res)
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

void SegSurface::print_seg(vector<Path>& pathes_temp)
{
	string RES_FILE_PATH = "D:/programfiles/3/11/build/result/" + result_path_seg + "/123/";
	for (int i = 0; i < pathes_temp.size(); i++)
	{
		string path_file0 = RES_FILE_PATH + to_string(i) + "//faceid.txt";
		ofstream out0(path_file0);
		for (int j = 0; j < pathes_temp[i].cover_face.size(); j++)
		{
			out0 << pathes_temp[i].cover_face[j] << endl;
		}
		out0.close();
		//cout << pathes_temp[i].cover_face.size() << endl;

		string path_file1 = RES_FILE_PATH + to_string(i) + "//pri.txt";
		ofstream out1(path_file1);
		for (int j = 0; j < pathes_temp[i].path.size() - 1; j++)
		{
			out1 << pathes_temp[i].path[j].x() << " " << pathes_temp[i].path[j].y() << " " << pathes_temp[i].path[j].z() << " " <<
				pathes_temp[i].path[j + 1].x() << " " << pathes_temp[i].path[j + 1].y() << " " << pathes_temp[i].path[j + 1].z() <<
				endl;
		}
		out1.close();
		//cout << pathes_temp[i].path.size() << endl;

		string path_file2 = RES_FILE_PATH + to_string(i) + "//rullings.txt";
		ofstream out2(path_file2);
		for (int j = 0; j < pathes_temp[i].rullings.size(); j++)
		{
			for (int j1 = 0; j1 < pathes_temp[i].rullings[j].size(); j1++)
			{
				out2 << pathes_temp[i].rullings[j][j1].p0.x() << " " << pathes_temp[i].rullings[j][j1].p0.y() << " " << pathes_temp[i].rullings[j][j1].p0.z() << " " <<
					pathes_temp[i].rullings[j][j1].p1.x() << " " << pathes_temp[i].rullings[j][j1].p1.y() << " " << pathes_temp[i].rullings[j][j1].p1.z() <<
					endl;
			}
		}
		out2.close();
		///cout << pathes_temp[i].rullings.size() << endl;

		string path_file3 = RES_FILE_PATH + to_string(i) + "//filter.txt";
		ofstream out3(path_file3);
		for (int j = 0; j < pathes_temp[i].info_every_face.size(); j++)
		{
			out3 << pathes_temp[i].info_every_face[j].face_idx << " " << pathes_temp[i].info_every_face[j].ruling_direction.x() << " " << pathes_temp[i].info_every_face[j].ruling_direction.y() << " " <<
				pathes_temp[i].info_every_face[j].ruling_direction.z() << " " << pathes_temp[i].info_every_face[j].centeriod[0] << " " << pathes_temp[i].info_every_face[j].centeriod[1] << " " <<
				pathes_temp[i].info_every_face[j].centeriod[2] << " ";
			for (auto cf : pathes_temp[i].info_every_face[j].cover_area_one_face)
			{
				out3 << cf << " ";
			}
			out3 << endl;
		}
		out3.close();
	}
}

vector<int> SegSurface::calc_face_stop_temp(vector<int>& many_faces_idx_temp)
{
	vector<int> face_stop; face_stop.clear();
	std::sort(many_faces_idx_temp.begin(), many_faces_idx_temp.end());
	int iters = 0;
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		if (i == many_faces_idx_temp[iters])
		{
			iters++;
		}
		else
		{
			face_stop.push_back(i);
		}
	}
	/*
	for (auto ele : face_singlirity)
	{
		face_stop.push_back(ele);
	}*/

	return face_stop;
}

void SegSurface::seg_surface_no_iter(vector<int>& branch, vector<Path>& pathes_temp, vector<int>& inaccessible_area)
{
	pathes_temp.clear();
	collect_init_informarion();
	
	for (auto ele : inaccessible_area)
	{
		face_singlirity.push_back(ele);
	}
	vector<SampleTwoPathesOnOneFace> sample_lines_temp; sample_lines_temp.clear();
	cout << branch.size() << endl;

	branch.erase(std::remove_if(branch.begin(), branch.end(), [&](int num) {
		return std::find(face_singlirity.begin(), face_singlirity.end(), num) != face_singlirity.end();
		}), branch.end());

	int target_num = sample_lines_num(branch);
	cout << branch.size() << " " << target_num << endl;
	
	sample_principal_curvatures_lines_on_many_faces(seg_param1, branch, seg_param2, sample_lines_temp);
	cout << "sample_principal_curvatures_lines_on_many_faces" << endl;
	generate_rullings_on_many_faces(seg_param3, seg_param4 * average_edge_length, seg_param5, sample_lines_temp);
	cout << "generate_rullings_on_many_faces" << endl;
	calc_fitting_area_on_many_faces(seg_param6 * average_edge_length, sample_lines_temp, pathes_temp);
	cout << "calc_fitting_area_on_many_faces" << endl;
	set_cover_algorithm(pathes_temp);
}

void SegSurface::data_clear()
{
	sample_path_temp.clear();
}

void SegSurface::test()
{
	cout << seg_param1 << endl;
	/*
	vector<int> many_faces_idx_temp = { 0,1,2,3,4,5,6 };
	vector<SampleTwoPathesOnOneFace> sample_lines_temp;
	sample_principal_curvatures_lines_on_many_faces(1, many_faces_idx_temp, 100, sample_lines_temp);
	
	ofstream out("point.txt");
	for (int i = 0; i < sample_lines_temp[0].two_lines_temp[0].size(); i++)
	{
		out << sample_lines_temp[0].two_lines_temp[0][i].v_pos[0] << " " << sample_lines_temp[0].two_lines_temp[0][i].v_pos[1] << " " << sample_lines_temp[0].two_lines_temp[0][i].v_pos[2] << endl;
	}
	out.close();*/
}

int SegSurface::sample_lines_num(vector<int>& branch)
{
	int sample_num = 0;
	if (branch.size() > 10000)
	{
		sample_num = 200;
	}
	else if (branch.size() <= 10000 && branch.size() > 5000)
	{
		sample_num = 100;
	}
	else if (branch.size() <= 5000 && branch.size() > 1000)
	{
		sample_num = 50;
	}
	else if (branch.size() <= 1000 && branch.size() > 100)
	{
		sample_num = 30;
	}
	else if (branch.size() <= 100 && branch.size() >= 1)
	{
		int nums = branch.size() / 2;
		sample_num = std::max(1, nums);
	}
	else
	{
		cout << "branch size is zero!!!" << endl;
	}
	
	return sample_num;
}



