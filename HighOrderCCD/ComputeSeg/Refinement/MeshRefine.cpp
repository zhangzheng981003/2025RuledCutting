#include "MeshRefine.h"

MeshRefine::MeshRefine(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh)
{
	Mesh& mesh = seg_mesh_.GetMesh();

	e_min_ = MeshTools::AverageEdgeLength(mesh) * e_min_rate_;
	e_max_ = MeshTools::AverageEdgeLength(mesh) * e_max_rate_;
	a_min_ = MeshTools::Area(mesh) / mesh.n_faces() * a_min_rate_;
}

void MeshRefine::CollapseShortEdges()
{
	std::cout << "------- Collapse ------" << std::endl;
	Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	OpenMesh::EPropHandleT<bool> seam;
	mesh.add_property(seam);
	for (const EH& e_h : mesh.edges())
	{
		mesh.property(seam, e_h) = seam_status[e_h.idx()];
	}

	//std::cout << 1 << std::endl;
	
	bool ok;
	int i;
	for (ok = false, i = 0; !ok && i < iter_num_; ++i) {
		ok = true;

		mesh.update_face_normals();

		for (auto ite = mesh.edges_begin(); ite != mesh.edges_end(); ite++)
		{
			//ite.enable_skipping();

			double e_len = mesh.calc_edge_length(*ite);

			HEH h10 = mesh.halfedge_handle(*ite, 0);
			HEH h01 = mesh.halfedge_handle(*ite, 1);

			double area_0 = mesh.calc_sector_area(h10);
			double area_1 = mesh.calc_sector_area(h01);

			if (e_len < e_min_ || area_0 < a_min_ || area_1 < a_min_)
			{
				//std::cout << 2 << std::endl;

				bool hcol01 = mesh.is_collapse_ok(h01);
				bool hcol10 = mesh.is_collapse_ok(h10);

				VH v_0 = mesh.to_vertex_handle(h10);
				VH v_1 = mesh.to_vertex_handle(h01);

				// boundary point
				if (mesh.is_boundary(v_0))
				{
					hcol01 = false;
				}

				if (mesh.is_boundary(v_1))
				{
					hcol10 = false;
				}

				//std::cout << 3 << std::endl;

				// property
				EH prev10, next10, prev01, next01;
				bool p10 = false, p01 = false;
				if (mesh.is_valid_handle(mesh.face_handle(h10)))
				{
					prev10 = mesh.edge_handle(mesh.prev_halfedge_handle(h10));
					next10 = mesh.edge_handle(mesh.next_halfedge_handle(h10));

					// boundary
					if (mesh.is_boundary(prev10)) hcol10 = false;
					if (mesh.is_boundary(next10)) hcol01 = false;

					// property
					p10 = mesh.property(seam, prev10) || mesh.property(seam, next10);
				}

				//std::cout << 4 << std::endl;

				if (mesh.is_valid_handle(mesh.face_handle(h01)))
				{
					prev01 = mesh.edge_handle(mesh.prev_halfedge_handle(h01));
					next01 = mesh.edge_handle(mesh.next_halfedge_handle(h01));


					// boundary
					if (mesh.is_boundary(prev01)) hcol01 = false;
					if (mesh.is_boundary(next01)) hcol10 = false;

					// property
					p01 = mesh.property(seam, prev01) || mesh.property(seam, next01);
				}

				//std::cout << 5 << std::endl;

				// collapse to feature edge
				//double cov_0 = 0, cov_1 = 0;
				double cov_0 = DBL_MAX, cov_1 = DBL_MAX;
				int cov_num_0 = 0, cov_num_1 = 0;
				OpenMesh::Vec3d adj_n_0, adj_n_1, he_vec;
				FH f_0, f_1;
				for (const HEH& adj_he : mesh.voh_range(v_0))
				{
					f_0 = mesh.face_handle(adj_he);
					f_1 = mesh.opposite_face_handle(adj_he);
					if (!mesh.is_valid_handle(f_0) || !mesh.is_valid_handle(f_1)) continue;

					adj_n_0 = mesh.normal(f_0);
					adj_n_1 = mesh.normal(f_1);
					double temp_cov = adj_n_0 | adj_n_1;

					//cov_0 += dot(adj_n_0, adj_n_1);
					//++cov_num_0;

					he_vec = mesh.calc_edge_vector(adj_he).normalized();
					double sign = he_vec | (adj_n_0 % adj_n_1);

					if (sign < 0 && temp_cov < collpase_cos_bound_)
					{
						cov_0 = DBL_MAX;
						break;
					}

					cov_0 = std::min(cov_0, dot(adj_n_0, adj_n_1));
				}

				for (const HEH& adj_he : mesh.voh_range(v_1))
				{
					f_0 = mesh.face_handle(adj_he);
					f_1 = mesh.opposite_face_handle(adj_he);
					if (!mesh.is_valid_handle(f_0) || !mesh.is_valid_handle(f_1)) continue;

					adj_n_0 = mesh.normal(f_0);
					adj_n_1 = mesh.normal(f_1);
					double temp_cov = adj_n_0 | adj_n_1;

					//cov_1 += dot(adj_n_0, adj_n_1);
					//++cov_num_1;

					he_vec = mesh.calc_edge_vector(adj_he).normalized();
					double sign = he_vec | (adj_n_0 % adj_n_1);

					if (sign < 0 && temp_cov < collpase_cos_bound_)
					{
						cov_1 = DBL_MAX;
						break;
					}

					cov_1 = std::min(cov_1, dot(adj_n_0, adj_n_1));
				}

				//cov_0 /= cov_num_0;
				//cov_1 /= cov_num_1;

				if (cov_0 < 0)
				{
					hcol01 = false;
				}

				if (cov_1 < 0)
				{
					hcol10 = false;
				}

				// avoid flip triangle
				OpenMesh::Vec3d normal_ori, normal_new;
				for (const HEH& adj_he : mesh.voh_range(v_1))
				{
					HEH next_he = mesh.next_halfedge_handle(adj_he);
					if (!next_he.is_valid()) continue;

					if (mesh.to_vertex_handle(adj_he) == v_0
						|| mesh.to_vertex_handle(next_he) == v_0) continue;

					normal_ori = (mesh.calc_edge_vector(adj_he)
						% mesh.calc_edge_vector(next_he)).normalized();

					normal_new = ((mesh.point(mesh.to_vertex_handle(adj_he)) - mesh.point(v_0))
						% mesh.calc_edge_vector(next_he)).normalized();

					if ((normal_new | normal_ori) < 0)
					{
						hcol10 = false;
						break;
					}
				}

				for (const HEH& adj_he : mesh.voh_range(v_0))
				{
					HEH next_he = mesh.next_halfedge_handle(adj_he);
					if (!next_he.is_valid()) continue;

					if (mesh.to_vertex_handle(adj_he) == v_1
						|| mesh.to_vertex_handle(next_he) == v_1) continue;

					normal_ori = (mesh.calc_edge_vector(adj_he)
						% mesh.calc_edge_vector(next_he)).normalized();

					normal_new = ((mesh.point(mesh.to_vertex_handle(adj_he)) - mesh.point(v_1))
						% mesh.calc_edge_vector(next_he)).normalized();

					if ((normal_new | normal_ori) < 0)
					{
						hcol01 = false;
						break;
					}
				}

				bool collaplse_status = false;
				if (cov_0 < cov_1 && hcol10)
				{
					bool skip = false;
					for (const VH& adj_v : mesh.vv_range(v_1))
					{
						if (skip) break;
						
						if (IsTooLong(mesh, v_0, adj_v))
						{
							skip = true;
						}
					}
						
					if (!skip)
					{
						mesh.collapse(h10);
						collaplse_status = true;
					}
				}
				else if (cov_0 >= cov_1 && hcol01)
				{
					bool skip = false;
					for (const VH& adj_v : mesh.vv_range(v_0))
					{
						if (skip) break;

						if (IsTooLong(mesh, v_1, adj_v))
						{
							skip = true;
						}
					}

					if (!skip)
					{
						mesh.collapse(h01);
						collaplse_status = true;
					}
				}

				if (collaplse_status)
				{
					if (mesh.is_valid_handle(mesh.face_handle(h10)))
					{
						if (!mesh.status(prev10).deleted())
						{
							mesh.property(seam, prev10) = p10;
						}

						if (!mesh.status(next10).deleted())
						{
							mesh.property(seam, next10) = p10;
						}
					}

					if (mesh.is_valid_handle(mesh.face_handle(h01)))
					{
						if (!mesh.status(prev01).deleted())
						{
							mesh.property(seam, prev01) = p01;
						}

						if (!mesh.status(next01).deleted())
						{
							mesh.property(seam, next01) = p01;
						}
					}

					ok = false;
				}

				//std::cout << 7 << std::endl;
			}
		}


		mesh.garbage_collection();
	}

	seam_status.assign(mesh.n_edges(), false);
	for (const EH& e_h : mesh.edges())
	{
		seam_status[e_h.idx()] = mesh.property(seam, e_h);
	}
	seg_mesh_.BoundToIdx();
	seg_mesh_.IdxToBound();
	seg_mesh_.BoundToIdx();

	//std::cout << 8 << std::endl;

	mesh.remove_property(seam);
}

void MeshRefine::SplitLongEdges()
{
	std::cout << "------- Split Long ------" << std::endl;
	Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	OpenMesh::EPropHandleT<bool> seam;
	mesh.add_property(seam);
	for (const EH& e_h : mesh.edges())
	{
		mesh.property(seam, e_h) = seam_status[e_h.idx()];
		mesh.status(e_h).set_locked(false);
	}

	bool ok;
	int i;
	for (ok = false, i = 0; !ok && i < iter_num_; ++i) {
		ok = true;

		for (const EH& e_h : mesh.edges())
		{
			if (mesh.status(e_h).locked()) continue;

			VH v0 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 0));
			VH v1 = mesh.to_vertex_handle(mesh.halfedge_handle(e_h, 1));

			if (mesh.calc_edge_length(e_h) > e_max_) {
				const typename Mesh::Point& p0 = mesh.point(v0);
				const typename Mesh::Point& p1 = mesh.point(v1);

				bool e_status = mesh.property(seam, e_h);

				auto vh = mesh.add_vertex((p0 + p1) * 0.5);
				mesh.split(e_h, vh);

				if (e_status)
				{
					for (const HEH& he_h : mesh.voh_range(vh))
					{
						if (mesh.to_vertex_handle(he_h) == v0 || mesh.to_vertex_handle(he_h) == v1)
						{
							mesh.property(seam, mesh.edge_handle(he_h)) = true;
						}
						else
						{
							mesh.property(seam, mesh.edge_handle(he_h)) = false;
						}
					}
				}
				else
				{
					for (const HEH& he_h : mesh.voh_range(vh))
					{
						mesh.property(seam, mesh.edge_handle(he_h)) = false;
					}
				}

				ok = false;
				continue;
			}

			//HEH he_h = mesh.halfedge_handle(e_h, 0);
			//
			//OpenMesh::Vec3d p0 = mesh.point(v0);
			//OpenMesh::Vec3d p1 = mesh.point(v1);
			//OpenMesh::Vec3d p2 = mesh.point(mesh.opposite_vh(he_h));
			//OpenMesh::Vec3d p3 = mesh.point(mesh.opposite_he_opposite_vh(he_h));

			//double angle0 = ((p2 - p0).normalized() | (p2 - p1).normalized());
			//double angle1 = ((p3 - p0).normalized() | (p3 - p1).normalized());

			//if (angle0 < large_split_cos_ || angle1 < large_split_cos_)
			//{
			//	bool e_status = mesh.property(seam, e_h);

			//	auto vh = mesh.add_vertex((p0 + p1) * 0.5);
			//	mesh.split(e_h, vh);

			//	if (e_status)
			//	{
			//		for (const HEH& adj_he : mesh.voh_range(vh))
			//		{
			//			if (mesh.to_vertex_handle(adj_he) == v0
			//				|| mesh.to_vertex_handle(adj_he) == v1)
			//			{
			//				mesh.property(seam, mesh.edge_handle(adj_he)) = true;
			//			}
			//			else
			//			{
			//				mesh.property(seam, mesh.edge_handle(adj_he)) = false;
			//			}
			//		}
			//	}
			//	else
			//	{
			//		for (const HEH& adj_he : mesh.voh_range(vh))
			//		{
			//			mesh.property(seam, mesh.edge_handle(adj_he)) = false;
			//		}
			//	}

			//	for (const HEH& adj_he : mesh.voh_range(vh))
			//	{
			//		mesh.status(mesh.edge_handle(adj_he)).set_locked(true);
			//		mesh.status(mesh.edge_handle(mesh.next_halfedge_handle(adj_he))).set_locked(true);
			//	}

			//	ok = false;
			//}
		
		}

	}

	seam_status.assign(mesh.n_edges(), false);
	for (const EH& e_h : mesh.edges())
	{
		seam_status[e_h.idx()] = mesh.property(seam, e_h);
	}

	seg_mesh_.BoundToIdx();
	seg_mesh_.IdxToBound();

	mesh.remove_property(seam);
}

void MeshRefine::SplitObtuseEdges()
{
	//std::cout << "------- Split Obtuse ------" << std::endl;

	//Mesh& mesh = seg_mesh_.GetMesh();
	//mesh.update_face_normals();
	//std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	//OpenMesh::EPropHandleT<bool> seam;
	//mesh.add_property(seam);
	//for (const EH& e_h : mesh.edges())
	//{
	//	mesh.property(seam, e_h) = seam_status[e_h.idx()];
	//	mesh.status(e_h).set_locked(false);
	//}


	//std::vector<std::pair<double, EH>> split_e_stack(0);
	//for (EH e_h : mesh.edges())
	//{
	//	if (mesh.is_boundary(e_h)) continue;

	//	HEH he_h = mesh.halfedge_handle(e_h, 0);

	//	OpenMesh::Vec3d normal_0 = mesh.normal(mesh.face_handle(he_h));
	//	OpenMesh::Vec3d normal_1 = mesh.normal(mesh.opposite_face_handle(he_h));

	//	double normal_cor = (normal_0 | normal_1);

	//	if (normal_cor < large_split_cos_)
	//	{
	//		split_e_stack.push_back({ normal_cor, e_h });
	//	}
	//}


	//std::vector<EH> split_e(0);
	//for (EH e_h:mesh.edges())
	//{
	//	if (mesh.is_boundary(e_h)) continue;
	//	
	//	HEH he_h = mesh.halfedge_handle(e_h, 0);

	//	OpenMesh::Vec3d normal_0 = mesh.normal(mesh.face_handle(he_h));
	//	OpenMesh::Vec3d normal_1 = mesh.normal(mesh.opposite_face_handle(he_h));
	//	
	//	if ((normal_0 | normal_1) < large_split_cos_)
	//	{
	//		split_e.push_back(e_h);
	//	}
	//}



	//for (EH e_h:split_e)
	//{
	//	bool ori_status = mesh.property(seam, e_h);

	//	OpenMesh::Vec3d mid_point = mesh.calc_edge_midpoint(e_h);

	//	HEH he_0 = mesh.halfedge_handle(e_h, 0);

	//	VH to_v = mesh.to_vertex_handle(he_0);
	//	VH from_v = mesh.from_vertex_handle(he_0);

	//	VH new_v = mesh.split_copy(e_h, mid_point);

	//	for (const HEH& adj_he:mesh.voh_range(new_v))
	//	{
	//		if (mesh.to_vertex_handle(adj_he) == to_v
	//			|| mesh.to_vertex_handle(adj_he) == from_v)
	//		{
	//			mesh.property(seam, mesh.edge_handle(adj_he)) = ori_status;
	//		}
	//		else
	//		{
	//			mesh.property(seam, mesh.edge_handle(adj_he)) = false;
	//		}
	//	}
	//}


	//seam_status.assign(mesh.n_edges(), false);
	//for (const EH& e_h : mesh.edges())
	//{
	//	seam_status[e_h.idx()] = mesh.property(seam, e_h);
	//}

	//seg_mesh_.BoundToIdx();
	//seg_mesh_.IdxToBound();

	//mesh.remove_property(seam);
}

void MeshRefine::SeamFlip()
{
	Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();

	OpenMesh::EPropHandleT<bool> seam;
	mesh.add_property(seam);
	for (const EH& e_h : mesh.edges())
	{
		mesh.property(seam, e_h) = seam_status[e_h.idx()];
	}

	std::vector<int> v_cout;
	seg_mesh_.ComputeVertexCount(v_cout);

	bool ok;
	int i;
	HEH prev_he_0, prev_he_1, next_he_0, next_he_1, he_0, he_1;
	EH prev_e_0, prev_e_1, next_e_0, next_e_1;
	OpenMesh::Vec3d vec_prev_0, vec_prev_1, vec_next_0, vec_next_1;
	for (ok = false, i = 0; !ok && i < iter_num_; ++i) {
		ok = true;

		for (const EH& e_h : mesh.edges())
		{
			if (mesh.property(seam, e_h) || mesh.is_boundary(e_h)) continue;

			if (!mesh.is_flip_ok(e_h)) continue;

			he_0 = mesh.halfedge_handle(e_h, 0);
			he_1 = mesh.halfedge_handle(e_h, 1);

			prev_he_0 = mesh.prev_halfedge_handle(he_0);
			next_he_0 = mesh.next_halfedge_handle(he_0);
			prev_he_1 = mesh.prev_halfedge_handle(he_1);
			next_he_1 = mesh.next_halfedge_handle(he_1);

			prev_e_0 = mesh.edge_handle(prev_he_0);
			next_e_0 = mesh.edge_handle(next_he_0);
			prev_e_1 = mesh.edge_handle(prev_he_1);
			next_e_1 = mesh.edge_handle(next_he_1);

			//is arrow?
			// arrow or triangle
			//    /|\
			//   / | \
			//  /  |  \
			// /___|___\

			double sign_0 = (mesh.calc_edge_vector(next_he_0) % mesh.calc_edge_vector(he_0))
				| (mesh.calc_edge_vector(next_he_0) % mesh.calc_edge_vector(prev_he_1));
			double sign_1 = (mesh.calc_edge_vector(next_he_1) % mesh.calc_edge_vector(he_1))
				| (mesh.calc_edge_vector(next_he_1) % mesh.calc_edge_vector(prev_he_0));

			if (sign_0 < 0 || sign_1 < 0) continue;


			double cos_0 = mesh.calc_edge_vector(next_he_0).normalized() 
				| mesh.calc_edge_vector(prev_he_1).normalized();
			//double cos_1 = mesh.calc_edge_vector(next_he_1).normalized()
			//	| mesh.calc_edge_vector(prev_he_0).normalized();
			//double sign_1 = (mesh.calc_edge_vector(next_he_1) % mesh.calc_edge_vector(he_1))
			//	| (mesh.calc_edge_vector(next_he_1) % mesh.calc_edge_vector(prev_he_0));


			//if (mesh.property(seam, prev_e_1) 
			//	&& mesh.property(seam, next_e_0) 
			//	&& cos_0 > flip_cos_
			//	&& sign_0 > 0
			//	&& v_cout[mesh.to_vertex_handle(he_0).idx()] == 2)
			//{
			//	// diamond
			//	//   /|\
			//	//  / | \
			//	//  \ | /
			//	//   \|/
			//	//
			//	if (cos_1 > flip_cos_ && sign_1 > 0)
			//	{
			//		mesh.flip(e_h);

			//		mesh.property(seam, prev_e_1) = false;
			//		mesh.property(seam, next_e_0) = false;
			//		mesh.property(seam, e_h) = true;

			//		v_cout[mesh.to_vertex_handle(he_0).idx()] = 0;
			//	}

			//	// arrow or triangle
			//	//    /|\
			//	//   / | \
			//	//  /  |  \
			//	// /___|___\
			//	//
			//	else if (cos_0 < cos_1){
			//		mesh.property(seam, prev_e_1) = false;
			//		mesh.property(seam, next_e_0) = false;

			//		mesh.property(seam, prev_e_0) = true;
			//		mesh.property(seam, next_e_1) = true;

			//		v_cout[mesh.to_vertex_handle(he_0).idx()] = 0;
			//		v_cout[mesh.to_vertex_handle(he_1).idx()] = 2;
			//	}
			//	
			//	ok = false;
			//}
			//else if (mesh.property(seam, prev_e_0)
			//		&& mesh.property(seam, next_e_1)
			//		&& cos_1 > flip_cos_
			//		&& sign_1 > 0
			//		&& v_cout[mesh.to_vertex_handle(he_1).idx()] == 2)
			//{
			//	// diamond
			//	//   /|\
			//	//  / | \
			//	//  \ | /
			//	//   \|/
			//	//
			//	if (cos_0 > flip_cos_ && sign_0 > 0)
			//	{
			//		mesh.flip(e_h);

			//		mesh.property(seam, prev_e_0) = false;
			//		mesh.property(seam, next_e_1) = false;
			//		mesh.property(seam, e_h) = true;

			//		v_cout[mesh.to_vertex_handle(he_1).idx()] = 0;
			//	}

			//	// arrow or triangle
			//	//    /|\
			//	//   / | \
			//	//  /  |  \
			//	// /___|___\
			//	//
			//	else if (cos_1 < cos_0){
			//		mesh.property(seam, prev_e_0) = false;
			//		mesh.property(seam, next_e_1) = false;

			//		mesh.property(seam, prev_e_1) = true;
			//		mesh.property(seam, next_e_0) = true;

			//		v_cout[mesh.to_vertex_handle(he_1).idx()] = 0;
			//		v_cout[mesh.to_vertex_handle(he_0).idx()] = 2;
			//	}

			//	ok = false;
			//}
		
		}

	}

	seam_status.assign(mesh.n_edges(), false);
	for (const EH& e_h : mesh.edges())
	{
		seam_status[e_h.idx()] = mesh.property(seam, e_h);
	}

	seg_mesh_.BoundToIdx();
	seg_mesh_.IdxToBound();

	mesh.remove_property(seam);
}

void MeshRefine::ValenceFlip()
{
	Mesh& mesh = seg_mesh_.GetMesh();
	const auto& seam_status = seg_mesh_.GetSeam();

	mesh.update_face_normals();

	bool ok;
	int i;
	HEH he_h;
	VH v_0, v_1, v_2, v_3;
	double val_0, val_1, val_2, val_3;
	double val_opt_0, val_opt_1, val_opt_2, val_opt_3;
	double ve_before, ve_after;
	OpenMesh::Vec3d normal_0, normal_1;
	for (ok = false, i = 0; !ok && i < iter_num_; ++i) {
		ok = true;

		for (const EH& e_h : mesh.edges())
		{
			if (mesh.is_boundary(e_h) || seam_status[e_h.idx()]) continue;

			he_h = mesh.halfedge_handle(e_h, 0);

			normal_0 = mesh.normal(mesh.face_handle(he_h));
			normal_1 = mesh.normal(mesh.opposite_face_handle(he_h));

			if ((normal_0 | normal_1) < flip_dih_cos_bound_) continue;

			v_0 = mesh.to_vertex_handle(he_h);
			v_2 = mesh.to_vertex_handle(mesh.next_halfedge_handle(he_h));

			he_h = mesh.halfedge_handle(e_h, 1);
			v_1 = mesh.to_vertex_handle(he_h);
			v_3 = mesh.to_vertex_handle(mesh.next_halfedge_handle(he_h));

			val_0 = mesh.valence(v_0);
			val_1 = mesh.valence(v_1);
			val_2 = mesh.valence(v_2);
			val_3 = mesh.valence(v_3);

			val_opt_0 = (mesh.is_boundary(v_0) ? 4 : 6);
			val_opt_1 = (mesh.is_boundary(v_1) ? 4 : 6);
			val_opt_2 = (mesh.is_boundary(v_2) ? 4 : 6);
			val_opt_3 = (mesh.is_boundary(v_3) ? 4 : 6);

			ve_before = (val_0 - val_opt_0) * (val_0 - val_opt_0)
				+ (val_1 - val_opt_1) * (val_1 - val_opt_1)
				+ (val_2 - val_opt_2) * (val_2 - val_opt_2)
				+ (val_3 - val_opt_3) * (val_3 - val_opt_3);

			--val_0;
			--val_1;
			++val_2;
			++val_3;

			ve_after = (val_0 - val_opt_0) * (val_0 - val_opt_0)
				+ (val_1 - val_opt_1) * (val_1 - val_opt_1)
				+ (val_2 - val_opt_2) * (val_2 - val_opt_2)
				+ (val_3 - val_opt_3) * (val_3 - val_opt_3);

			if (ve_before > ve_after && mesh.is_flip_ok(e_h) && !edge_flip_flips_normal(mesh, e_h)) {
				
				mesh.flip(e_h);
				ok = false;
			}
		}
	}

	seg_mesh_.BoundToIdx();
}

void MeshRefine::TangentialSmoothing(bool _use_projection)
{
	Mesh& mesh = seg_mesh_.GetMesh();

	typename Mesh::VIter     v_it, v_end(mesh.vertices_end());
	typename Mesh::CVVIter   vv_it;
	typename Mesh::Point     u, n;

	OpenMesh::VPropHandleT<Mesh::Point> update;
	mesh.add_property(update);

	// tag active vertices
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		mesh.status(*v_it).set_tagged(!mesh.status(*v_it).locked() &&
			!mesh.status(*v_it).feature() &&
			!mesh.is_boundary(*v_it));

	mesh.update_face_normals();

	double damping = 1.0;

	// smooth
	for (int iters = 0; iters < 10; ++iters)
	{
		for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		{
			if (mesh.status(*v_it).tagged())
			{
				u.vectorize(0.0);
				int valence = 0;
				int feature_valence = 0;

				for (auto heh : mesh.voh_range(*v_it))
				{
					u += mesh.point(mesh.to_vertex_handle(heh));
					if (mesh.status(mesh.edge_handle(heh)).feature())
						++feature_valence;
					++valence;
				}

				if (valence)
				{
					u *= (1.0 / valence);
					u -= mesh.point(*v_it);
					n = mesh.normal(*v_it);
					n *= (u | n);
					u -= n;
				}

				if (feature_valence == 2)
				{
					// project update onto feature directions
					typename Mesh::Point feature_dir(0.0, 0.0, 0.0);
					for (auto heh : mesh.voh_range(*v_it))
					{
						if (mesh.status(mesh.edge_handle(heh)).feature())
						{
							auto dir = mesh.point(mesh.to_vertex_handle(heh)) - mesh.point(mesh.from_vertex_handle(heh));
							if ((dir | feature_dir) > 0)
								feature_dir += dir;
							else
								feature_dir -= dir;
						}
					}
					if (feature_dir.sqrnorm() > 0)
						feature_dir.normalize();
					u = (feature_dir | u) * feature_dir;
				}
				else if (feature_valence != 0)
				{
					// more than two or exactly one incident feature edge means vertex should be preserved
					u.vectorize(0.0);
				}

				mesh.property(update, *v_it) = u;
			}
		}

		for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
			if (mesh.status(*v_it).tagged())
				mesh.point(*v_it) += damping * mesh.property(update, *v_it);

		// check if normals changed
		bool normals_changed = false;
		for (auto fh : mesh.faces())
		{
			if ((mesh.calc_face_normal(fh) | mesh.normal(fh)) < 0)
			{
				normals_changed = true;
				break;
			}
		}

		if (normals_changed)
		{
			// revert update and try again with smaller step
			for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
				if (mesh.status(*v_it).tagged())
					mesh.point(*v_it) -= damping * mesh.property(update, *v_it);
			damping *= 0.5;
		}
	}


	// reset tagged status
	for (v_it = mesh.vertices_begin(); v_it != v_end; ++v_it)
		mesh.status(*v_it).set_tagged(false);


	//// project
	//if (_use_projection)
	//	for (v_it = mesh_.vertices_begin(); v_it != v_end; ++v_it)
	//		if (!mesh_.status(*v_it).locked() && !mesh_.status(*v_it).feature())
	//			project_to_reference(*v_it);


	mesh.remove_property(update);
}

bool MeshRefine::edge_flip_flips_normal(const Mesh& mesh, EH e_h)
{
	if (mesh.is_boundary(e_h))
		return true;

	auto heh = mesh.halfedge_handle(e_h, 0);

	// get the four points of the two incident faces in ccw order
	auto p0 = mesh.point(mesh.to_vertex_handle(heh));
	auto p1 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(heh)));
	auto p2 = mesh.point(mesh.from_vertex_handle(heh));
	auto p3 = mesh.point(mesh.to_vertex_handle(mesh.next_halfedge_handle(mesh.opposite_halfedge_handle(heh))));

	// compute normals before flip
	auto n0_before = (p1 - p0) % (p2 - p0);
	auto n1_before = (p2 - p0) % (p3 - p0);

	// compute normals after flip
	auto n0_after = (p1 - p0) % (p3 - p0);
	auto n1_after = (p3 - p2) % (p1 - p2);

	// compare all pairs of before and after normals
	if ((n0_before | n0_after) < 0 || (n1_before | n1_after) < 0 ||
		(n0_before | n1_after) < 0 || (n1_before | n0_after) < 0)
		return true;

	return false;
}
