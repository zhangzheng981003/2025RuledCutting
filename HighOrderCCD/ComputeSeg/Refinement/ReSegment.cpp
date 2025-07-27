#include "ReSegment.h"

ReSegment::ReSegment(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh), ori_mesh_(seg_mesh)
{
	seg_class_ = new Segmentation(ori_mesh_, Segmentation::SegMode::CUR);
}

ReSegment::~ReSegment()
{
	if (seg_class_)
	{
		delete seg_class_;
		seg_class_ = NULL;
	}
}

bool ReSegment::Run(int v_ring_radius, double coarse_angle_bound)
{
	std::cout << "------ coarse segment ------" << std::endl;

	std::vector<int> v_seg;

	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();

	std::vector<bool>& ori_seam_status = ori_mesh_.GetSeam();

	std::vector<double> abs_v_gauss;

	int max_id;
	double max_cur;

	//do
	//{
		abs_v_gauss = seg_mesh_.GetAbsGauss(true);

		std::cout << "seg_num: " << seg_num << std::endl;

		seg_mesh_.ComputeVertexSeg(v_seg);

		ComputeLargeCurvature(abs_v_gauss, v_seg, v_ring_radius, coarse_angle_bound);

		ComputeSegCurvature(abs_v_gauss, max_id, max_cur);

		std::cout << max_cur << std::endl;
		std::cout << coarse_angle_bound << std::endl;

		//system("pause");

		if (max_cur < coarse_angle_bound)
		{
			std::cout << "Min curvature!" << std::endl;
			return false;
		}

		std::vector<bool> f_status(mesh.n_faces(), false);
		for (const FH& f_h : mesh.faces())
		{
			if (seg_id[f_h.idx()] == max_id)
			{
				f_status[f_h.idx()] = true;
			}
		}

		ori_seam_status = seam_status;
		ori_mesh_.BoundToIdx();

		std::cout << "old num: " << ori_mesh_.GetSegNum() << std::endl;

		seg_class_->Init(f_status);
		bool run_status = seg_class_->Run();

		std::cout << "run status: " << run_status << std::endl;

		if (run_status)
		{
			seam_status = ori_seam_status;
			seg_mesh_.BoundToIdx();

			std::cout << "new num: " << seg_mesh_.GetSegNum() << std::endl;

			return true;
		}
		else
		{
			return false;
		}


	//} while (true);
}

void ReSegment::ComputeLargeCurvature(
	std::vector<double>& abs_v_gauss,
	const std::vector<int>& v_seg, 
	const int& v_ring_radius, 
	const double& seg_angle_bound)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();

	std::vector<double> neibor_max_cur(abs_v_gauss.size(), -DBL_MAX);
	for (int i = 0; i < mesh.n_vertices(); i++)
	{
		std::vector<int> v_stack = {i};
		int cout = 0, ori_size;
		while (v_stack.size() != 0 && cout < v_ring_radius)
		{
			ori_size = v_stack.size();
			for (size_t j = 0; j < ori_size; j++)
			{
				int temp_idx = v_stack[j];
				neibor_max_cur[temp_idx] =
					std::max(abs_v_gauss[temp_idx], neibor_max_cur[temp_idx]);

				for (const VH& adj_v : mesh.vv_range(VH(temp_idx)))
				{
					if (v_seg[adj_v.idx()] == -1) continue;

					v_stack.push_back(adj_v.idx());
				}
			}

			v_stack.erase(v_stack.begin(), v_stack.begin() + ori_size);
			++cout;
		}
	}
	
	std::cout << 1 << std::endl;

	abs_v_gauss = neibor_max_cur;

	//// vertices on seam
	//std::vector<int> v_stack(0);
	//for (const VH& v_h:mesh.vertices())
	//{
	//	if (v_seg[v_h.idx()] != -1) continue;
	//	
	//	v_stack.push_back(v_h.idx());
	//}
	//
	//// 5-ring of seam
	//int cout = 0, ori_size;
	//while (v_stack.size() != 0 && cout < v_ring_radius)
	//{
	//	ori_size = v_stack.size();
	//	for (size_t i = 0; i < ori_size; i++)
	//	{
	//		VH temp_v = VH(v_stack[i]);

	//		abs_v_gauss[temp_v.idx()] = 0;

	//		for (const VH& adj_v:mesh.vv_range(temp_v))
	//		{
	//			if (v_seg[adj_v.idx()] == -1) continue;
	//			
	//			v_stack.push_back(adj_v.idx());
	//		}
	//	}

	//	v_stack.erase(v_stack.begin(), v_stack.begin() + ori_size);

	//	cout++;
	//}
	//

	//// 5-ring of inner vertices
	//if (v_ring_radius != 0)
	//{
	//	for (const VH& v_h : mesh.vertices())
	//	{
	//		if (abs_v_gauss[v_h.idx()] < seg_angle_bound) continue;

	//		v_stack = { v_h.idx() };
	//		bool high_status = false;

	//		cout = 0;
	//		while (v_stack.size() != 0 && cout < v_ring_radius)
	//		{
	//			ori_size = v_stack.size();
	//			for (size_t i = 0; i < ori_size; i++)
	//			{
	//				VH temp_v = VH(v_stack[i]);

	//				for (const VH& adj_v : mesh.vv_range(temp_v))
	//				{
	//					if (abs_v_gauss[adj_v.idx()] < seg_angle_bound)
	//					{
	//						v_stack.push_back(adj_v.idx());
	//						continue;
	//					}
	//					else
	//					{
	//						high_status = true;
	//						break;
	//					}
	//				}

	//				if (high_status) break;
	//			}

	//			if (high_status) break;

	//			v_stack.erase(v_stack.begin(), v_stack.begin() + ori_size);

	//			cout++;
	//		}

	//		if (!high_status)
	//		{
	//			abs_v_gauss[v_h.idx()] = 0;
	//		}
	//	}
	//}

	std::cout << *std::max_element(abs_v_gauss.begin(), abs_v_gauss.end()) << std::endl;
	std::cout << seg_angle_bound << std::endl;

	for (const VH& v_h : mesh.vertices())
	{
		if (abs_v_gauss[v_h.idx()] > seg_angle_bound) continue;

		abs_v_gauss[v_h.idx()] = 0;
	}

	std::cout << *std::max_element(abs_v_gauss.begin(), abs_v_gauss.end()) << std::endl;

	std::cout << 2 << std::endl;

	ori_mesh_.SetModifiedGauss(abs_v_gauss);
}

void ReSegment::ComputeSegCurvature(
	const std::vector<double>& abs_v_gauss, 
	int& max_id,
	double& max_cur)
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	int& seg_num = seg_mesh_.GetSegNum();
	std::vector<bool>& seam_status = seg_mesh_.GetSeam();
	std::vector<int>& seg_id = seg_mesh_.GetSegId();

	max_cur = -DBL_MAX;
	for (EH e_h : mesh.edges())
	{
		if (seam_status[e_h.idx()] || mesh.is_boundary(e_h)) continue;

		HEH he_h = mesh.halfedge_handle(e_h, 0);

		VH v_0 = mesh.to_vertex_handle(he_h);
		VH v_1 = mesh.from_vertex_handle(he_h);

		int temp_f_seg = seg_id[mesh.face_handle(he_h).idx()];

		double temp_cur = std::max(abs_v_gauss[v_0.idx()], abs_v_gauss[v_1.idx()]);

		if (temp_cur > max_cur)
		{
			max_cur = temp_cur;
			max_id = temp_f_seg;
		}
	}
}
