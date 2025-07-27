#pragma once
#include "../SegMesh.h"

class MeshRefine
{
public:
	MeshRefine(SegMesh& seg_mesh);

	void CollapseShortEdges();

	void SplitLongEdges();

	void SplitObtuseEdges();

	void SeamFlip();

	void ValenceFlip();

	void TangentialSmoothing(bool _use_projection);


private:
	bool edge_flip_flips_normal(const Mesh& mesh, EH e_h);

	bool IsTooLong(const Mesh& mesh, const VH& v, const VH& adj_v)
	{
		return (mesh.point(v) - mesh.point(adj_v)).norm() > e_max_;
	}

private:
	SegMesh& seg_mesh_;

	int iter_num_ = 100;

	double e_min_;
	double a_min_;
	double e_max_;

	double e_min_rate_ = 0.4;
	double e_max_rate_ = 1.5;
	double a_min_rate_ = 1e-5;

	double collpase_cos_bound_ = cos(M_PI / 180 * 170);
	double large_split_cos_ = cos(M_PI / 180 * 150);
	double flip_cos_ = cos(M_PI / 180 * 120);
	double flip_dih_cos_bound_ = cos(M_PI / 180 * 5);
};

