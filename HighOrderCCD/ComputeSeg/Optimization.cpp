#include "Optimization.h"
#include "CutMesh.h"

Optimization::Optimization(SegMesh& seg_mesh)
	:seg_mesh_(seg_mesh)
{
	Mesh& mesh = seg_mesh.GetMesh();

	dev_info_.ResetDiffBound(mesh.n_edges());
}

Optimization::~Optimization()
{
}

void Optimization::DevelopAndSegmentation(std::vector<double>& sub_fiedler, vector<int>& real2img)
{
	std::cout << "======== Dev and Seg ========" << std::endl;

	Mesh& mesh = seg_mesh_.GetMesh();
	const std::vector<bool>& seam_status = seg_mesh_.GetSeam();
	
	//calc_rullings_energy_based_edge(mesh, dev_info_.dev_energy_, res,real2img);
	auto& dev_energy = dev_info_.dev_energy_;
	int& cout = dev_info_.cout_;

	Segmentation seg(seg_mesh_, Segmentation::SegMode::DIFF);
	//
	//MeshDevelop deve;
	////deve.SmoothStatus(true);
	//deve.oriMesh(mesh);
	//deve.tarMesh(mesh, target_dist_);

	Merge merge(seg_mesh_);

	std::vector<int> false_seg_idx(0);

	//cout = 0;
	//do
	//{
	//	if (!seg_status_)
	//	{
	//		cout++;
	//		if (cout >= max_iter_)
	//		{
	//			std::cout << "------ Max Iter ------" << std::endl;
	//			break;
	//		}

	//		std::cout << "------ Seg Reset ------" << std::endl;
	//		seg_mesh_.Init();
	//	}

		std::cout << "------ Get Large ------" << std::endl;
		std::vector<bool> f_status;
		int max_idx;
		if (!GetLargeSegFaces(false_seg_idx, f_status, max_idx)) 
		{
			seg_status_ = false;
			false_seg_idx.clear();
			//continue;
			//break;
		}

		std::cout << "------ Segmentation ------" << std::endl;
		seg.Init(f_status);
		seg_status_ = seg.Run();
		seg.GetFiedler(sub_fiedler);

	//	if (!seg_status_)
	//	{
	//		false_seg_idx.push_back(max_idx);
	//		seg_status_ = true;
	//	}

	//} while (true);

	std::cout << "------ Output Dev Info ------ " << std::endl;
	seg_mesh_.BoundToIdx();
	dev_info_.UpdateCurvatureEnergy(seg_mesh_);
	dev_info_.UpdateSegEnergy(seg_mesh_);
	
	std::cout << "energy: " << dev_info_.seg_energy_[0] << endl;
	dev_info_.OutputInfo();

	std::cout << "------ Output Seg and Mesh ------ " << std::endl;
	seg_mesh_.WriteMesh(FILE_PATH + "seg_mesh.obj");
	seg_mesh_.WriteSeg(FILE_PATH + "seg_face_id.txt");
}

void Optimization::Run()
{
	//pre_processing();
	
	//std::vector<double> sub_fiedler;
	//DevelopAndSegmentation(sub_fiedler);
	//Refinement();
}

void Optimization::Refinement(/*std::vector<OpenMesh::Vec3d>& v_pos*/)
{
	
}

//void Optimization::FineRefinement()
//{
//	ReSegment re_seg(seg_mesh_);
//	re_seg.Run(ReSegment::RunStatus::FINE);
//	seg_mesh_.WriteSeg(FILE_PATH + "fine_seg_seg.txt");
//
//	Merge mer(seg_mesh_);
//	mer.SegMerge(Merge::MergeStatus::FINE);
//	seg_mesh_.WriteSeg(FILE_PATH + "fine_merge.txt");
//
//	std::cout << "------ Output Seg and Mesh ------ " << std::endl;
//	seg_mesh_.WriteMesh(FILE_PATH + "fine_deform.obj");
//	seg_mesh_.WriteSeg(FILE_PATH + "fine_deform_seg.txt");
//}

//void Optimization::Merge()
//{
//	//Refine::Merge(seg_mesh_);
//	//seg_mesh_.WriteSeg(FILE_PATH + "merge_seg.txt");
//}

//void Optimization::Smooth()
//{
//	for (size_t i = 0; i < 100; i++)
//	{
//		Refine::MoveSeam(seg_mesh_);
//	}
//	Refine::SmoothSeam(seg_mesh_);
//
//	Refine::SplitLongEdges(seg_mesh_);
//	Refine::SplitObtuseEdges(seg_mesh_);
//
//	Refine::CollapseShortEdges(seg_mesh_);
//
//	seg_mesh_.WriteMesh(FILE_PATH + "smooth_mesh.obj");
//	seg_mesh_.WriteSeg(FILE_PATH + "smooth_seg.txt");
//}





void Optimization::calc_rullings_energy_based_edge(Mesh& mesh, vector<double>& rullings_energy, vector<vector<double>>& res, vector<int>& real2img)
{
	vector<Vector3d> face_rulings; face_rulings.resize(mesh.n_faces());
	
	//测试的网格一共12个面，前6个面一个方向，后六个面反向
	//理想的输出前六个面是一块，后6个面是一块
	//vector<vector<double>> res;
	//read_txt("D://Zheng2024//0Code//BSpline//0226//build//new_filter.txt", res);
	for (int i = 0; i < mesh.n_faces(); i++)
	{
		face_rulings[real2img[res[i][0]]] = { res[i][1],res[i][2],res[i][3] };
		//face_rulings[i + mesh.n_faces() / 2] = { 0,1,0 };
	}

	rullings_energy.clear();
	rullings_energy.resize(mesh.n_edges());
	for (auto eh : mesh.edges())
	{
		if (mesh.is_boundary(eh))
		{
			rullings_energy[eh.idx()] = 0;
		}
		else
		{
			FH f0 = mesh.face_handle(mesh.halfedge_handle(eh, 0));
			FH f1 = mesh.face_handle(mesh.halfedge_handle(eh, 1));
			double dist = (face_rulings[f0.idx()] - face_rulings[f1.idx()]).norm();
			rullings_energy[eh.idx()] = dist;
		}
	}
}

void Optimization::read_txt(const string& pathname, vector<vector<double>>& res)
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


void Optimization::pre_processing()
{
	vector<vector<double>> res;
	read_txt("E:/0RuledCutting/code/smooth/build/new_filter.txt", res);
	std::string filename = "C:/Users/ustc-gcl/Desktop/inputmesh.obj";
	std::cout << "file: " << filename << std::endl;

	std::cout << "read mesh " << std::endl;
	Mesh mesh;
	MeshTools::ReadMesh(mesh, filename);
	
	vector<int> new_patch(mesh.n_faces(), 0);
	ofstream ioszhzh("final_seg.txt");
	for (int i = 0; i < res.size(); i++)
	{
		new_patch[res[i][0]] = 1;
	}

	for (int i = 0; i < new_patch.size(); i++)
	{
		ioszhzh << new_patch[i] << endl;
	}
	ioszhzh.close();

	/*
	vector<int> seam_stutas(mesh.n_edges(), 0);
	int iter = 0;
	for (auto a : mesh.edges())
	{
		auto f1 = mesh.face_handle(mesh.halfedge_handle(a, 0));
		auto f2 = mesh.face_handle(mesh.halfedge_handle(a, 1));

		if (new_patch[f1.idx()] != new_patch[f2.idx()])
		{
			seam_stutas[a.idx()] = 1;
			iter++;
		}
	}*/
}

int Optimization::GetLargeSegFaces(
	const std::vector<int>& false_seg_idx,
	std::vector<bool>& f_status,
	int& max_idx)
{
	//std::cout << 1 << std::endl;
	// seg idx
	std::vector<int>& seg_id = seg_mesh_.GetSegId();
	Mesh& mesh = seg_mesh_.GetMesh();

	//std::cout << 2 << std::endl;

	// def info
	dev_info_.UpdateCurvatureEnergy(seg_mesh_, true);
	auto high_cur = dev_info_.high_cur_;

	//std::cout << 3 << std::endl;

	dev_info_.UpdateSegEnergy(seg_mesh_);
	auto seg_energy = dev_info_.seg_energy_;

	for (int i : false_seg_idx)
	{
		seg_energy[i] = 0;
	}

	//std::cout << 4 << std::endl;

	double& max_energy = dev_info_.max_energy_;

	//std::cout << 5 << std::endl;

	// seg with largest energy
	max_idx = -1;
	do
	{
		auto max_iter = std::max_element(seg_energy.begin(), seg_energy.end());

		if ((*max_iter) < min_seg_energy_) break;

		int temp_idx = std::distance(seg_energy.begin(), max_iter);
		if (true)//if (high_cur[temp_idx] > min_cur_bound_)
		{
			//std::cout << "max energy: " << *max_iter << std::endl;
			max_energy = *max_iter;
			max_idx = temp_idx;
			break;
		}
		*max_iter = 0;

	} while (true);

	//std::cout << 6 << std::endl;

	if (max_idx == -1)
	{
		return false;
	}
	else
	{
		f_status.assign(mesh.n_faces(), false);
		for (const FH& f_h : mesh.faces())
		{
			if (seg_id[f_h.idx()] == max_idx)
			{
				f_status[f_h.idx()] = true;
			}
		}

		return true;
	}

	//std::cout << 7 << std::endl;
}

bool Optimization::IsSingleHigh()
{
	const Mesh& mesh = seg_mesh_.GetMesh();
	std::vector<double> abs_v_gauss = seg_mesh_.GetAbsGauss(true);

	int max_id = std::max_element(abs_v_gauss.begin(), abs_v_gauss.end()) - abs_v_gauss.begin();

	for (const VH& adj_v:mesh.vv_range(VH(max_id)))
	{
		if (abs_v_gauss[adj_v.idx()] > abs_v_gauss[max_id] * 0.1)
		{
			return false;
		}
	}

	return true;

	//return false;
}
