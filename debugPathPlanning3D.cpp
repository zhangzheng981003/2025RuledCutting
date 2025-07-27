#include "HighOrderCCD/Config/Config.h"
#include "HighOrderCCD/Optimization/Optimization3D_time.h"
#include "HighOrderCCD/Optimization/Optimization3D_point.h"
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <vector>
#include <ctime>
#include <fstream>          
#include <io.h>
#include <process.h>

#include "HighOrderCCD/Wrapping/GreedyWrap.h"

USE_PRJ_NAMESPACE
#define M_PI 3.14159265358979323846
using namespace std;

int mesh_normlized(Mesh& mesh)
{
	Mesh::Point minCoord(std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max());
	Mesh::Point maxCoord(-std::numeric_limits<float>::max(), -std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Mesh::Point coord = mesh.point(*v_it);
		minCoord.minimize(coord);
		maxCoord.maximize(coord);
	}

	Mesh::Point center = (minCoord + maxCoord) * 0.5f;

	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Mesh::Point& coord = mesh.point(*v_it);
		coord -= center;
	}

	float scaleFactor = 1.0f / (maxCoord - minCoord).max();


	for (Mesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		Mesh::Point& coord = mesh.point(*v_it);
		coord *= scaleFactor;
	}

	
	if (!MeshTools::WriteMesh(mesh, "normalized_mesh.obj"))
	{
		std::cerr << "error" << std::endl;
		return 1;
	}
}

void read_txt(const string& pathname, vector<vector<double>>& res)
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

void read_txt_string(const string& pathname, vector<vector<string>>& res)
{
	res.clear();
	ifstream infile;
	infile.open(pathname.data());
	assert(infile.is_open());
	vector<string> suanz;
	string s;
	while (getline(infile, s))
	{
		istringstream is(s);
		string d;
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

double seg_param1, seg_param2, seg_param3, seg_param4, seg_param5, seg_param6;
string result_path_seg, result_path_fitting,result_path_fitting1;
int main()
{
	vector<vector<string>> res;
	read_txt_string("input_info.txt", res);
	for (int i = 0; i < res.size(); i++)
	{
		const std::string input_mesh_file = res[i][0];
		result_path_seg = res[i][0];
		result_path_fitting = res[i][0];
		result_path_fitting1= res[i][0];
		Mesh input_mesh;
		MeshTools::ReadMesh(input_mesh, input_mesh_file);
		seg_param1 = std::stold(res[i][1]);
		seg_param2 = std::stold(res[i][2]);
		seg_param3 = std::stold(res[i][3]);
		seg_param4 = std::stold(res[i][4]); 
		seg_param5 = std::stold(res[i][5]);
		seg_param6 = std::stold(res[i][6]);
		
		GreedyWrap gw(input_mesh);
		gw.greedy_wrap();
	}

	return 0;
}
