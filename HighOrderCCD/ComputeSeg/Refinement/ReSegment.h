#pragma once
#include "../SegMesh.h"
#include "../Segmentation/Segmentation.h"

class ReSegment
{
public:
	enum RunStatus {COARSE, FINE};
public:
	ReSegment(SegMesh& seg_mesh);

	~ReSegment();

	bool Run(int v_ring_radius, double coarse_angle_bound);

private:
	void ComputeLargeCurvature(
		std::vector<double>& abs_v_gauss,
		const std::vector<int>& v_seg, 
		const int& v_ring_radius,
		const double& seg_angle_bound);

	void ComputeSegCurvature(
		const std::vector<double>& abs_v_gauss,
		int& max_id,
		double& max_cur);

private:
	SegMesh& seg_mesh_;
	SegMesh ori_mesh_;

	Segmentation* seg_class_;
};

