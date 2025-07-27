#ifndef GREEDYWRAP_H
#define GREEDYWRAP_H

#include "HighOrderCCD/SplineFitting/Mesh/MeshDefinition.h"
#include "HighOrderCCD/SplineFitting/BSplineNoCollision.h"
#include "HighOrderCCD/ComputeSeg/PatchSeg.h"

using namespace std;
//typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;
extern string result_path_fitting;
class GreedyWrap
{
public:
	struct BSplineSurfaceInfo
	{
		//λ��
		int wrapping_idx;
		int branch_idx;
		int surface_idx;
		//������Ϣ
		Handle(Geom_BSplineSurface) bspline_handle = nullptr;
		Eigen::MatrixXd ctr_points;
		vector<double> u_knots_temp;
		//�������
		SegSurface::Path path;
		//������
		double mean_dist = 0.0;
		double max_dist = 0.0;
		double dist_energy = 0.0;


		~BSplineSurfaceInfo()
		{
			if (bspline_handle)
			{
				bspline_handle = nullptr;
			}
		}
	};
	
	struct VertexFittingStatus
	{
		int vertex_idx;
		vector<pair<BSplineSurfaceInfo, double>>  vertex_bspline_dist_temp;
	};
	
	struct FaceAdjVertex
	{
		int face_idx;
		int p0;
		int p1;
		int p2;
	};





	static bool cmp_bigger(const pair<BSplineSurfaceInfo, double>& a, const pair<BSplineSurfaceInfo, double>& b)
	{
		return a.second < b.second;
	}

	static bool cmp_smaller_path_temp(const SegSurface::Path& path0, const SegSurface::Path& path1)
	{
		return path0.cover_face.size() > path1.cover_face.size();
	}

	GreedyWrap(Mesh& mesh);
	~GreedyWrap();
	void initialization();
	void write_txt(const string& file_path, vector<int>&input_vec);
	void read_txt(const string& pathname, vector<vector<double>>& res);
	vector<int> vertex2face(vector<int>& vertexId);
	int select_big_patches(vector<SegSurface::Path>& pathes_temp,double coverage_rate_threshold, int cover_face_nums);
	void rotate_vector(Vector3d& n, double angle, Vector3d& v, Vector3d& result);
	bool is_segment_intersection(Eigen::Vector3d start_point, Eigen::Vector3d end_point);
	void search_inaccessible_area(double high_value, int sample_rulings_num_one_face);
	void build_bspline_surface_handle(Mesh& mesh, Handle(Geom_BSplineSurface)& bspline_handle, Eigen::MatrixXd& ctr_points, vector<double>& u_knots_temp);
	void cut_mesh(Mesh& mesh, Mesh& seg_mesh, SegSurface::Path& path);
	void calc_seg_energy(Mesh& mesh, Mesh& seg_mesh, SegSurface::Path path, Handle(Geom_BSplineSurface)& bspline_handle, Eigen::MatrixXd& ctr_points, 
		vector<double>& u_knots_temp, vector<SegEnergy::InnerEdge>& ie_temp, vector<double>& rulings_energy, vector<int>& real2img, vector<int>& img2real);
	void calc_vertex_dist_for_bspline_surfaces(vector<vector<vector<BSplineSurfaceInfo>>>& bspline_res_3temp, vector<VertexFittingStatus>& vertex_fitting_temp, int start_pos);
	void seg_one_path_to_two_path(SegSurface::Path& ori_path, int break_face_id, vector<SegSurface::Path>& patches_temp);
	
	void greedy_wrap();
	void greedy_wrap3();
	void greedy_segment();
	//debug
	void debug_write_pathes_temp(const string& file_name, vector<SegSurface::Path>& pathes_temp);
	void debug_read_pathes_temp(const string& file_name, vector<SegSurface::Path>& pathes_temp);


private:
	Mesh mesh;
	Tree my_tree;
	Tree_ tree;
	std::vector<MyTriangle> my_cgal_triangles;

	std::vector<Triangle> cgal_triangles;
	vector<SegSurface::FaceBasisInfo> face_basis_info_temp;
	vector<GraphAlgorithm::VertexInfo> vertex_temp;

	double average_edge_length = 0.0;
	double rulings_length = 0.0;
	double dist_threshold = 0.0;
	int small_patch_threshold = 20;
	double mean_dist_threshold = 0.0;
	double max_dist_threshold = 0.0;
	double dist_energy_threshold = 0.0;
	double vertex_dist_threshold = 0.0;

	vector<VertexFittingStatus> vertex_fitting_temp;
	vector<int> inaccessible_area;//���������ߵ�face_stop
	vector<int> inaccessible_point;//��������BSpline�ľ���
	vector<bool> inaccessible_point_status;
	vector<bool> inaccessible_face_status;
	vector<int> sucess_cover_point;//�޳���inaccessible_point
	vector<int> failed_cover_point;//�޳���inaccessible_point

};
#endif