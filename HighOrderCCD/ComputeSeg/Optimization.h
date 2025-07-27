#pragma once
#include <queue>
#include <fstream>
#include "SegMesh.h"
#include "Development\MeshDevelop.h"
#include "Segmentation\Segmentation.h"
#include "Refinement\ReSegment.h"
#include "Refinement\Merge.h"
#include "Refinement\MeshRefine.h"
#include "Refinement\GaussMin.h"

const std::string FILE_PATH = "mesh_data/";
using namespace std;
class DeforInfo
{
public:
	void UpdateCurvatureEnergy(SegMesh& seg_mesh, bool remove_neibor = false)
	{
		const Mesh& mesh  = seg_mesh.GetMesh();
		const auto& seam_status = seg_mesh.GetSeam();
		const auto& seg_id = seg_mesh.GetSegId();
		const int& seg_num = seg_mesh.GetSegNum();

		high_cur_.assign(seg_num, 0);

		std::vector<double> abs_v_gauss = seg_mesh.GetAbsGauss(true);

		if (remove_neibor)
		{
			for (EH e_h : mesh.edges())
			{
				if (!mesh.is_boundary(e_h) && !seam_status[e_h.idx()]) continue;

				HEH he_0 = mesh.halfedge_handle(e_h, 0);

				VH to_v = mesh.to_vertex_handle(he_0);
				VH from_v = mesh.from_vertex_handle(he_0);

				for (const VH& adj_v:mesh.vv_range(to_v))
				{
					abs_v_gauss[adj_v.idx()] = 0;
				}

				for (const VH& adj_v : mesh.vv_range(from_v))
				{
					abs_v_gauss[adj_v.idx()] = 0;
				}
			}
		}

		for (EH e_h : mesh.edges())
		{
			if (mesh.is_boundary(e_h)) continue;

			if (seam_status[e_h.idx()]) continue;

			HEH he_0 = mesh.halfedge_handle(e_h, 0);

			int to_idx = mesh.to_vertex_handle(he_0).idx();
			int from_idx = mesh.from_vertex_handle(he_0).idx();

			int temp_seg_id = seg_id[mesh.face_handle(he_0).idx()];

			high_cur_[temp_seg_id] = std::max(high_cur_[temp_seg_id], abs_v_gauss[to_idx]);
			high_cur_[temp_seg_id] = std::max(high_cur_[temp_seg_id], abs_v_gauss[from_idx]);
		}

		max_cur_ = *std::max_element(high_cur_.begin(), high_cur_.end());
	}

	bool IsSmallCurvature()
	{
		return max_cur_ < M_PI / 180 / 2;
	}

	void UpdateSegEnergy(SegMesh& seg_mesh)
	{
		const Mesh& mesh = seg_mesh.GetMesh();
		const auto& seam_status = seg_mesh.GetSeam();
		const auto& seg_id = seg_mesh.GetSegId();
		const int& seg_num = seg_mesh.GetSegNum();

		seg_energy_.assign(seg_num, 0);

		// seg total energy
		for (EH e_h : mesh.edges())
		{
			if (mesh.is_boundary(e_h)) continue;

			if (seam_status[e_h.idx()]) continue;

			HEH he_0 = mesh.halfedge_handle(e_h, 0);
			HEH he_1 = mesh.halfedge_handle(e_h, 1);

			EH prev_0 = mesh.edge_handle(mesh.prev_halfedge_handle(he_0));
			EH next_0 = mesh.edge_handle(mesh.next_halfedge_handle(he_0));
			EH prev_1 = mesh.edge_handle(mesh.prev_halfedge_handle(he_1));
			EH next_1 = mesh.edge_handle(mesh.next_halfedge_handle(he_1));

			if (seam_status[prev_0.idx()] ||
				seam_status[next_0.idx()] ||
				seam_status[prev_1.idx()] ||
				seam_status[next_1.idx()]) continue;

			FH f_0 = mesh.face_handle(he_0);
			FH f_1 = mesh.face_handle(he_1);

			if (seg_id[f_0.idx()] != seg_id[f_1.idx()]) continue;

			//
			seg_energy_[seg_id[f_0.idx()]] = std::max(dev_energy_[e_h.idx()], seg_energy_[seg_id[f_0.idx()]]);
		}
	}

	void ResetDiffBound(int n_e)
	{
		energy_diff_bound_ *= 4 * n_e;
	}

	bool IsSmallDiff()
	{
		return abs(prev_energy_ - cur_energy_) < energy_diff_bound_;
	}

	void UpdateCurEnergy(double cur_energy)
	{
		prev_energy_ = cur_energy_;
		cur_energy_ = cur_energy;
	}

	void OutputInfo()
	{
		std::ofstream file_infor(FILE_PATH + "deform_and_seg_information.txt");

		file_infor << "Iter Num: " << cout_ << std::endl;
		file_infor << "Prev energy: " << prev_energy_ << std::endl;
		file_infor << "Cur energy: " << cur_energy_ << std::endl;
		file_infor << "Energy diff: " << prev_energy_ - cur_energy_ << std::endl;
		file_infor << "Max Energy: " << max_energy_ << std::endl;

		file_infor << "Max Curvature: " << max_cur_ << std::endl;

		file_infor << "Seg Num: " << high_cur_.size() << std::endl;

		file_infor << "Seg Energy and Curvature: " << std::endl;
		for (size_t i = 0; i < high_cur_.size(); i++)
		{
			file_infor << i << "th energy: " << seg_energy_[i] << "  high curvature: " << high_cur_[i] << std::endl;
		}

		file_infor.close();
	}

public:
	int cout_ = 0;
	double prev_energy_ = DBL_MAX;
	double cur_energy_ = DBL_MAX / 2;
	double max_energy_ = DBL_MAX;
	double max_cur_ = DBL_MAX;

	std::vector<double> dev_energy_;

	std::vector<double> high_cur_;
	std::vector<double> seg_energy_;

	double energy_diff_bound_ = 1e-8;
};

class Optimization
{
public:
	//Optimization(Mesh& mesh);

	Optimization(SegMesh& seg_mesh);

	~Optimization();

	void DevelopAndSegmentation(std::vector<double>& sub_fiedler, vector<int>& real2img);

	//void Collapse();
	
	//void ReSegment();

	void Run();

	void Refinement(/*std::vector<OpenMesh::Vec3d>& v_pos = {}*/);

	//void Merge();

	//void Smooth();

	void calc_rullings_energy_based_edge(Mesh& mesh, vector<double>& rullings_energy, vector<vector<double>>& res, vector<int>& real2img);
	void pre_processing();
	void read_txt(const string& pathname, vector<vector<double>>& res);

private:
	int GetLargeSegFaces(
		const std::vector<int>& false_seg_idx,
		std::vector<bool>& f_status,
		int& max_idx);

	bool IsSingleHigh();

	bool RemoveHighTriangle(const double& high_cur)
	{
		Mesh& mesh = seg_mesh_.GetMesh();
		std::vector<bool>& seam_status = seg_mesh_.GetSeam();
		std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss(true);

		bool change_status = false;
		for (const EH& e_h:mesh.edges())
		{
			if (!seam_status[e_h.idx()]) continue;
			
			HEH he_0 = mesh.halfedge_handle(e_h, 0);
			HEH he_1 = mesh.halfedge_handle(e_h, 1);

			VH oppo_v_0 = mesh.opposite_vh(he_0);
			VH oppo_v_1 = mesh.opposite_vh(he_1);

			if (oppo_v_0.is_valid() && abs_v_gauss[oppo_v_0.idx()] > high_cur - DBL_EPSILON)
			{
				seam_status[e_h.idx()] = false;
				seam_status[mesh.edge_handle(
					mesh.next_halfedge_handle(he_0)).idx()] = true;
				seam_status[mesh.edge_handle(
					mesh.prev_halfedge_handle(he_0)).idx()] = true;

				change_status = true;
			}

			if (oppo_v_1.is_valid() && abs_v_gauss[oppo_v_1.idx()] > high_cur - DBL_EPSILON)
			{
				seam_status[e_h.idx()] = false;
				seam_status[mesh.edge_handle(
					mesh.next_halfedge_handle(he_1)).idx()] = true;
				seam_status[mesh.edge_handle(
					mesh.prev_halfedge_handle(he_1)).idx()] = true;

				change_status = true;
			}
		}

		seg_mesh_.BoundToIdx();

		return change_status;
	}

	bool IsInnerHigh(const double& max_cur)
	{
		const Mesh& mesh = seg_mesh_.GetMesh();
		const std::vector<bool>& seam_status = seg_mesh_.GetSeam();
		std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss(true);

		std::vector<bool> seam_v(mesh.n_vertices());
		seg_mesh_.ComputeVertexStatus(seam_v);
		
		std::vector<int> v_stack(0);
		for (size_t i = 0; i < seam_v.size(); i++)
		{
			if (!seam_v[i]) continue;
			
			v_stack.push_back(i);
		}

		std::cout << *std::max_element(abs_v_gauss.begin(), abs_v_gauss.end()) << std::endl;

		int cout = 0, ori_size;
		while (v_stack.size() != 0 && cout < 5)
		{
			ori_size = v_stack.size();
			for (size_t j = 0; j < ori_size; j++)
			{
				int temp_idx = v_stack[j];
				abs_v_gauss[temp_idx] = 0;
				seam_v[temp_idx] = true;
				
				for (const VH& adj_v : mesh.vv_range(VH(temp_idx)))
				{
					if (seam_v[adj_v.idx()]) continue;

					v_stack.push_back(adj_v.idx());
				}
			}

			v_stack.erase(v_stack.begin(), v_stack.begin() + ori_size);
			++cout;
		}

		std::cout << *std::max_element(abs_v_gauss.begin(), abs_v_gauss.end()) << std::endl;
		double new_max_cur = *std::max_element(abs_v_gauss.begin(), abs_v_gauss.end());

		std::cout << max_cur << std::endl;

		//system("pause");

		if (new_max_cur < max_cur - DBL_EPSILON)
		{
			return false;
		}
		else
		{
			return true;
		}
	}

	void Smoothing()
	{
		Mesh& mesh = seg_mesh_.GetMesh();
		std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss(true);

		std::vector<std::pair<double, int>> gauss_idx;
		for (size_t i = 0; i < abs_v_gauss.size(); i++)
		{
			gauss_idx.push_back({ abs_v_gauss[i], i });
		}

		std::sort(gauss_idx.begin(), gauss_idx.end(), std::greater<std::pair<double, int>>());

		for (size_t i = 0; i < gauss_idx.size(); i++)
		{
			VH v_h = VH(i);
			
			if (mesh.is_boundary(v_h)) continue;

			OpenMesh::Vec3d temp_direction(0, 0, 0);
			int cont = 0;
			for (const VH& adj_v : mesh.vv_range(v_h))
			{
				temp_direction += mesh.point(adj_v) - mesh.point(v_h);
				cont++;
			}
			temp_direction /= cont;

			mesh.set_point(v_h, mesh.point(v_h) + temp_direction);
		}
	}


	
private:
	SegMesh& seg_mesh_;

	//Mesh mesh_;
	
	const int sub_iter_ = 5000;
	const double target_dist_ = 0.005;
	const int max_iter_ = 20;//20
	const double min_seg_energy_ = 1e-5;
	const double min_cur_bound_ = M_PI / 180;

	bool seg_status_ = false;


public:
	DeforInfo dev_info_;

};


