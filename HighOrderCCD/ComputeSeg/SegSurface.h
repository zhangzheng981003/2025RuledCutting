#ifndef SEGSURFACE_H
#define SEGSURFACE_H

#include "HighOrderCCD\ComputeSeg\SegHeader.h"
//#include "HighOrderCCD/Config/Config.h"

using namespace std;
using namespace Eigen;
extern double seg_param1, seg_param2, seg_param3, seg_param4, seg_param5, seg_param6;
extern string result_path_seg;
class SegSurface
{
public:
	struct FaceBasisInfo
	{
		int face_idx;
		double face_area;
		Eigen::Vector3d face_normal;
		Eigen::Vector3d centeriod;
	};

	struct OneRuling
	{
		// 1.p0与p1为直母线的两个端点
		Vector3d p0;Vector3d p1;
		// 2.直母线上采取的点集 计算拟合区域
		vector<Vector3d> one_rulling_sample_point_temp;
		// 3.拟合区域
		vector<int> one_ruling_cover_face;
	};

	struct OneFaceInfoInLine
	{
		int rank;
		int face_idx;
		int field_id;
		HEH heh;// = nullptr;
		Eigen::Vector3d face_normal;
		Eigen::Vector3d v_pos;
		Eigen::Vector3d direction;
		Eigen::Vector3d real_dir;
		Eigen::Vector4d plane;
		vector<OneRuling> rulings;
		Eigen::Vector3d ruling_direction;
		Eigen::Vector3d centeriod;
		set<int> cover_area_one_face;

		int rulings_num = 0;
		bool status = true;
	};

	struct SampleTwoPathesOnOneFace
	{
		int start_face_idx;
		vector<vector<OneFaceInfoInLine>> four_lines_temp;
		vector<vector<OneFaceInfoInLine>> two_lines_temp;
	};

	struct Path
	{
		vector<OneFaceInfoInLine> info_every_face;
		vector<Vector3d> path;
		vector<int> cover_face;
		vector<vector<OneRuling>> rullings;
		//vector<FilterInfo> filter_info_temp;
		
		//初始化信息
		Vector3d center_point;
		Vector3d translation_direction;
		Vector3d path_direction;
		Vector3d ruling_direction;
		double path_length = 0;
		double translation_length = 0;
		int ctr_num;
	};

	struct FaceConnect
	{
		int f_idx;
		int connnect_num = -1;
		vector<int> Adj_faces;
		vector<int> Adj_connect_faces;
		vector<int> Adj_unconnect_faces;
	};

	SegSurface(Mesh& mesh);
	~SegSurface();
	void collect_init_informarion();
	vector<int> sample_target_face(int target_num, vector<int>& available_areas);
	vector<vector<OneFaceInfoInLine>> two_principal_curvatures_lines_on_one_face(int one_face_idx, int line_length, vector<int>& face_stop);
	void sample_principal_curvatures_lines_on_many_faces(int target_num, vector<int>& many_faces_idx_temp, int line_length, vector<SampleTwoPathesOnOneFace>& sample_lines_temp);
	void generate_rullings_on_one_face(int ruling_gap, double high_value, int rulings_sample_num, OneFaceInfoInLine& ofiil);
	void generate_rullings_on_many_faces(int ruling_gap, double high_value, int rulings_sample_num, vector<SampleTwoPathesOnOneFace>& sample_lines_temp);
	void calc_fitting_area_on_many_faces(double epsilon, vector<SampleTwoPathesOnOneFace>& sample_lines_temp, vector<Path>& pathes_temp);
	//void calc_fitting_area_on_many_faces(vector<SampleTwoPathesOnOneFace>& sample_lines_temp);
	int spline_ray_intersection_(Eigen::Vector3d start_point, Eigen::Vector3d end_point, double epsilon);
	void generate_connect_branch(vector<int>& faceId);
	void set_cover_algorithm(vector<Path>& pathes_temp);
	void read_txt(const string& pathname, vector<vector<double>>& res);
	void print_seg(vector<Path>& pathes_temp);

	vector<int> calc_face_stop_temp(vector<int>& many_faces_idx_temp);
	void seg_surface_no_iter(vector<int>& branch, vector<Path>& pathes_temp, vector<int>& inaccessible_area);
	void data_clear();
	void test();
	int sample_lines_num(vector<int>& branch);


public:
	//Mesh mesh;

	Mesh mesh;
	Tree my_tree;
	Tree_ tree;
	LoopGen::LoopGen* lg = nullptr;

	Eigen::Matrix3Xd crossfield;
	std::vector<int> matching;
	vector<int> face_singlirity;
	std::vector<std::vector<int>> point_singularity;
	vector<FaceBasisInfo> face_basis_info_temp;
	double ruling_gap;
	double line_angle_threshold =0;
	std::vector<Triangle> cgal_triangles; std::vector<MyTriangle> my_cgal_triangles;
	double average_edge_length;
	Mesh::Point ptMin;
	Mesh::Point ptMax;
	double ruling_length = 0.0;
	int candidate_set_cover_min_num = 50;
	int candidate_u_face_min_num = 10;
	
	vector<OneFaceInfoInLine> sample_path_temp;
	
};

#endif