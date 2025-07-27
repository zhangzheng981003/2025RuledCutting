#pragma once
#include "../SplineFitting/Mesh/MeshDefinition.h"

typedef OpenMesh::VertexHandle VH;
typedef OpenMesh::FaceHandle FH;
typedef OpenMesh::EdgeHandle EH;
typedef OpenMesh::HalfedgeHandle HEH;

class CutMesh
{
public:
	CutMesh(std::string file_path, Mesh mesh, std::vector<int> init_seam = {})
	{
		ori_mesh_ = mesh;

		seam_status_.resize(mesh.n_edges(), 0);
		for (int& e_idx:init_seam)
		{
			seam_status_[e_idx] = 1;
		}
		FILE_PATH = file_path;
	}

	void AddNewSeam(std::vector<int> new_ids)
	{
		if (new_ids.size() != 0)
		{
			add_record_.push_back(new_ids);

			for (int& new_id : new_ids)
			{
				seam_status_[new_id]++;
			}
		}
	}

	void UndoAdd()
	{
		if (add_record_.size() > 0)
		{
			std::vector<int> last_add = add_record_.back(); add_record_.pop_back();

			for (int& undo_id : last_add)
			{
				seam_status_[undo_id]--;
			}
		}
	}

	Mesh Compute()
	{
		Mesh cut_mesh;
		res_two_patches.clear();
		//Mesh ori_mesh = cut_mesh_;
		std::vector<VH> halfedge_split_v_idx(ori_mesh_.n_halfedges());
		for (const VH& v_h : ori_mesh_.vertices())
		{
			if (ori_mesh_.is_boundary(v_h))
			{
				int v_split_num = 1;
				for (const EH& e_h : ori_mesh_.ve_range(v_h))
				{
					if (seam_status_[e_h.idx()] > 0) v_split_num++;
				}

				std::vector<VH> add_v(0);
				for (size_t i = 0; i < v_split_num; i++)
				{
					add_v.push_back(cut_mesh.add_vertex(ori_mesh_.point(v_h)));
				}

				int cout = 0;
				for (const HEH& adj_he : ori_mesh_.voh_range(v_h))
				{
					if (!ori_mesh_.face_handle(adj_he).is_valid()) continue;

					halfedge_split_v_idx[
						ori_mesh_.next_halfedge_handle(adj_he).idx()] = add_v[cout % v_split_num];

						if (seam_status_[ori_mesh_.edge_handle(adj_he).idx()] > 0) cout++;
				}
			}
			else
			{
				int v_split_num = 0;
				for (const EH& e_h : ori_mesh_.ve_range(v_h))
				{
					if (seam_status_[e_h.idx()] > 0) v_split_num++;
				}

				if (v_split_num == 0) v_split_num++;

				std::vector<VH> add_v(0);
				for (size_t i = 0; i < v_split_num; i++)
				{
					add_v.push_back(cut_mesh.add_vertex(ori_mesh_.point(v_h)));
				}

				int cout = 0;
				for (const HEH& adj_he : ori_mesh_.voh_range(v_h))
				{
					halfedge_split_v_idx[
						ori_mesh_.next_halfedge_handle(adj_he).idx()] = add_v[cout % v_split_num];

						if (seam_status_[ori_mesh_.edge_handle(adj_he).idx()] > 0) cout++;
				}
			}
		}

		for (const FH& f_h : ori_mesh_.faces())
		{
			std::vector<VH> f_v_stack(0);
			for (const HEH& adj_he : ori_mesh_.fh_range(f_h))
			{
				f_v_stack.push_back(halfedge_split_v_idx[adj_he.idx()]);
			}

			cut_mesh.add_face(f_v_stack);
		}
		Mesh cut_mesh_ = cut_mesh;
		MeshTools::WriteOBJ(cut_mesh, FILE_PATH + "/bound_mesh.obj");
		//OpenMesh::IO::write_mesh(cut_mesh, FILE_PATH+"/bound_mesh.obj", 0x0000, 12);

		OpenMesh::FPropHandleT<int> component_id;
		cut_mesh.add_property(component_id);
		for (const FH& f_h : cut_mesh.faces())
		{
			cut_mesh.property(component_id, f_h) = -1;
		}

		int component_num = 0;
		do
		{
			FH start_f = FH(-1);
			for (const FH& f_h : cut_mesh.faces())
			{
				if (cut_mesh.property(component_id, f_h) != -1) continue;

				start_f = f_h; break;
			}

			if (!start_f.is_valid()) break;

			std::vector<FH> f_stack;
			f_stack.push_back(start_f);
			do
			{
				FH top_f = f_stack.back(); f_stack.pop_back();
				cut_mesh.property(component_id, top_f) = component_num;

				for (HEH fh_h : cut_mesh.fh_range(top_f))
				{
					FH oppo_f_h = cut_mesh.opposite_face_handle(fh_h);

					if (!cut_mesh.is_valid_handle(oppo_f_h)) continue;
					if (cut_mesh.property(component_id, oppo_f_h) != -1) continue;

					f_stack.push_back(oppo_f_h);
				}
			} while (f_stack.size() != 0);

			++component_num;

		} while (true);

		for (int i = 0; i < component_num; i++)
		{
			Mesh temp_mesh = cut_mesh;

			for (FH f_h : temp_mesh.faces())
			{
				if (cut_mesh.property(component_id, f_h) != i)
				{
					temp_mesh.delete_face(f_h, true);
				}
			}

			temp_mesh.garbage_collection();

			const std::string seg_mesh_name = FILE_PATH + "/seg_" + std::to_string(i) + ".obj";
			MeshTools::WriteMesh(temp_mesh, seg_mesh_name);
			res_two_patches.push_back(temp_mesh);
		}
		cut_mesh.remove_property(component_id);
		return cut_mesh_;
	}

	std::vector<std::vector<int>>& GetAddRecord()
	{
		return add_record_;
	}

	std::vector<Mesh> res_two_patches;

private:
	Mesh ori_mesh_;
	//Mesh cut_mesh_;
	std::vector<int> seam_status_;

	std::vector<std::vector<int>> add_record_;

	std::string FILE_PATH;
};

