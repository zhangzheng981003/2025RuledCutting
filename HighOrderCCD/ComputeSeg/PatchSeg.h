#ifndef PATCHSEG_H
#define PATCHSEG_H

#include "Optimization.h"
#include "CutMesh.h"
#include "HighOrderCCD/ComputeSeg/SegEnergy.h"

class PatchSeg
{
public:
	struct FaceConnect
	{
		int f_idx;
		int connnect_num = -1;
		vector<int> Adj_faces;
		vector<int> Adj_connect_faces;
		vector<int> Adj_unconnect_faces;
	};


	PatchSeg(Mesh& mesh);
	~PatchSeg();
	void read_txt(const string& path, vector<vector<double>>& res);
	void write_txt(const string& path, vector<int>& res);
	void seg_great_distortion_rulings_patch(const string& filename0, const string& filename1);
	void generate_connect_branch(int numId);
	void connect_one_branch(vector<int>& faceId, vector<bool>& inaccessible_face_status);
	void search_inner_hole();
	void calc_edge_seg_energy00(Mesh& seg_mesh, VectorXd& BSpline, vector<double>& u_temp, vector<SegEnergy::InnerEdge>& ie_temp, vector<int>& img2real);
	
	void run();

public:
	Mesh mesh;

};
#endif
