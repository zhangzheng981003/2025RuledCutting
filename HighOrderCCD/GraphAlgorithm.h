#ifndef GRAPHALGORITHM_H
#define GRAPHALGORITHM_H

#include "HighOrderCCD\SplineFitting/Mesh/MeshDefinition.h"
#include "HighOrderCCD\graph.hpp"
#include <fstream>
#include <iostream>
#include <deque>
#include <queue>

using namespace std;
class GraphAlgorithm
{
public:
	struct VertexInfo
	{
		int v_idx;
		vector<int> f_idx_temp;
	};
	
	struct EdgeInfo 
	{
		int e_idx;
		int v_idx_i;
		int v_idx_j;
	};

	GraphAlgorithm(Mesh& mesh);
	void mesh2graph();
	void fast_search_connected_components(std::vector<int>& interest_area);
	vector<vector<int>> search_connect_branchs(std::vector<int>& interest_area_point);
	std::vector<int> face2vertex(std::vector<int>& interest_area_face);
	std::vector<int> vertex2face(std::vector<int>& interest_area_point);
	void delete_small_patches(std::vector<int>& interest_area_face);



public:
	Mesh mesh;
	vector<VertexInfo> vertex_temp;

};
#endif