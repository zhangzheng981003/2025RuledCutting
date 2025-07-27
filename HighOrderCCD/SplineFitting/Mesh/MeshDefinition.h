#pragma once
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>
//#include "interface.h"
#include <Eigen/Dense>
#include "HighOrderCCD/Config/Config.h"
struct MeshTraits : public OpenMesh::DefaultTraits
{
	typedef OpenMesh::Vec3d Point;
	typedef OpenMesh::Vec3d Normal;
	typedef OpenMesh::Vec2d TexCoord2D;
	VertexAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	FaceAttributes(OpenMesh::Attributes::Status | OpenMesh::Attributes::Normal);
	EdgeAttributes(OpenMesh::Attributes::Status);
	HalfedgeAttributes(OpenMesh::Attributes::Status);

	FaceTraits
	{
		FaceT()
		{};
	};

	EdgeTraits
	{
		EdgeT()
		{
		};
	public:
		double weight;
	};

	HalfedgeTraits
	{
		HalfedgeT()
		{
		};
	};

	VertexTraits
	{
		VertexT() : new_pos_fixed(false)
		{
		};
	private:
		OpenMesh::Vec3d new_pos;//can be used for deformation and parameterization
		bool new_pos_fixed;
	public:
		void set_New_Pos(const OpenMesh::Vec3d& n_p) { new_pos = n_p; };
		OpenMesh::Vec3d& get_New_Pos() { return new_pos; };
		void set_new_pos_fixed(bool f) { new_pos_fixed = f; };
		bool get_new_pos_fixed() { return new_pos_fixed; };
	};
};
typedef OpenMesh::TriMesh_ArrayKernelT<MeshTraits> Mesh;
typedef OpenMesh::EdgeHandle EH;
typedef OpenMesh::FaceHandle FH;
typedef OpenMesh::HalfedgeHandle HEH;
typedef OpenMesh::VertexHandle VH;
typedef OpenMesh::Vec3d Vec3d;
typedef Mesh::HalfedgeHandle HalfedgeHandle;

double calc_mesh_ave_edge_length(Mesh* mesh_);
bool is_flip_ok_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);//just copy the code from openmesh
bool flip_openmesh(Mesh::EdgeHandle& eh, Mesh& mesh_);
bool check_in_triangle_face(const std::vector<OpenMesh::Vec3d>& tri, const OpenMesh::Vec3d& p);
bool baryCoord(const OpenMesh::Vec3d& _p, const OpenMesh::Vec3d& _u, const OpenMesh::Vec3d& _v, const OpenMesh::Vec3d& _w, OpenMesh::Vec3d& _result);
void compute_point_area(Mesh* mesh_, std::vector<std::map<int, double>>& cornerArea, std::vector<double>& pointArea, bool use_np = false);
void rot_coord_sys(const OpenMesh::Vec3d& old_u, const OpenMesh::Vec3d& old_v,
	const OpenMesh::Vec3d& new_norm,
	OpenMesh::Vec3d& new_u, OpenMesh::Vec3d& new_v);
void proj_curv(const OpenMesh::Vec3d& old_u, const OpenMesh::Vec3d& old_v,
	double old_ku, double old_kuv, double old_kv,
	const OpenMesh::Vec3d& new_u, const OpenMesh::Vec3d& new_v,
	double& new_ku, double& new_kuv, double& new_kv);
//// Given a curvature tensor, find principal directions and curvatures
// Makes sure that pdir1 and pdir2 are perpendicular to normal
void diagonalize_curv(const OpenMesh::Vec3d& old_u, const OpenMesh::Vec3d& old_v,
	double ku, double kuv, double kv,
	const OpenMesh::Vec3d& new_norm,
	OpenMesh::Vec3d& pdir1, OpenMesh::Vec3d& pdir2, double& vk1, double& vk2);
void compute_principal_curvature(Mesh* mesh_,
	std::vector<double>& K1, std::vector<double>& K2,
	std::vector<OpenMesh::Vec3d>& dir1, std::vector<OpenMesh::Vec3d>& dir2);
template <typename T>
void initMeshStatusAndNormal(T& m)
{
	m.request_vertex_status();
	m.request_edge_status();
	m.request_face_status();

	m.request_face_normals();
	m.request_vertex_normals();

	m.update_face_normals();
	m.update_vertex_normals();
}

class MeshTools
{
public:
	static bool ReadMesh(Mesh & mesh, const std::string & filename);
	static bool ReadOBJ(Mesh & mesh, const std::string & filename);
	//static bool ReadOFF(Mesh & mesh, const std::string & filename);
	static bool WriteMesh(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static bool WriteOBJ(const Mesh & mesh, const std::string & filename, const std::streamsize & precision = 6);
	static bool WriteVertexTextutedOBJ(const Mesh& mesh, const std::string& filename, const std::streamsize& precision = 6);
	static bool WriteFaceTextutedOBJ(const Mesh& mesh, const std::string& filename, const std::streamsize& precision = 6);
	static double Area(const Mesh & mesh);
	static double AverageEdgeLength(const Mesh & mesh);
	static bool HasBoundary(const Mesh & mesh);
	static bool HasOneComponent(const Mesh & mesh);
	static int Genus(const Mesh & mesh);
	static void BoundingBox(const Mesh & mesh, Mesh::Point & bmax, Mesh::Point & bmin);
	static void Reassign(const Mesh & mesh1, Mesh & mesh2);
	static void ComputeGaussianCurvature(const Mesh& mesh, std::vector<double>& v_gauss, const bool& bound_remove = true);
};

//void convertMeshTo(Mesh& src, GCLF::Simplification::SMeshT& dst);
//void convertMeshFrom(GCLF::Simplification::SMeshT& src, Mesh& dst);