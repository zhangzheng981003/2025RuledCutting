#ifndef BSPLINENOCOLLISION_H
#define BSPLINENOCOLLISION_H

#include "HighOrderCCD/Config/Config.h"
#include "HighOrderCCD/SplineInit.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include "HighOrderCCD/Optimization/Optimization3D_time.h"

USE_PRJ_NAMESPACE
#define M_PI 3.14159265358979323846
using namespace std;
using namespace Eigen;

PRJ_BEGIN

class BSplineNoCollision
{
public:
	struct ResInfo
	{
		Data spline;
		vector<double> u_knot;
		double arv_e;
		double max_e;
		double min_e;
		double data_e;
	};


	BSplineNoCollision(Mesh& mesh);
	~BSplineNoCollision();
	void bspline_fitting_init();
	void bspline_fitting_init_plane(SegSurface::Path& one_path);
	BSplineNoCollision::ResInfo bspline_fitting_no_collision(const std::string& input_mesh_file, Mesh& input_mesh, int patch_num, SegSurface::Path& one_path);
	BSplineNoCollision::ResInfo bspline_fitting_no_collision1(const std::string& input_mesh_file, Mesh& input_mesh, int patch_num,int seg_num, SegSurface::Path& one_path);
public:
	Mesh mesh;
	//SurfaceSegment* ss = nullptr;

	//SegSurface* ss = nullptr;
	//SplineInit* si = nullptr;
	int ctr_num;
};

PRJ_END
#endif