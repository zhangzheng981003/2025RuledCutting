#include "GraphAlgorithm.h"

GraphAlgorithm::GraphAlgorithm(Mesh& mesh):mesh(mesh)
{

}

void GraphAlgorithm::fast_search_connected_components(std::vector<int>& interest_area)
{
	int fnum = mesh.n_faces();
	std::vector<int> visit_status(fnum, -1);
	std::vector<bool> interest_status(fnum, false);
	std::deque<int> v_q;

	for (int i = 0; i < interest_area.size(); i++)
	{
		interest_status[interest_area[i]] = true;
	}

	visit_status[interest_area[0]] = 1;
	for (auto ff : mesh.ff_range(mesh.face_handle(interest_area[0])))
	{
		if (interest_status[ff.idx()] == true)
		{
			v_q.push_front(ff.idx());
		}
		else
		{
			v_q.push_back(ff.idx());
		}
	}

	std::vector<std::vector<int>> branchs;
	std::vector<int> branch; branch.clear();
	vector <bool> lables(fnum, false);// vector <bool> lablesss;
	while (!v_q.empty())
	{
		int front = v_q.front();
		//cout << front << " "<< v_q.size() << endl;
		visit_status[front] = 1;
		v_q.pop_front();

		if (visit_status[front] == 1)
		{
			continue;
		}

		if (lables[front] == false)
		{
			branchs.push_back(branch);
			branch.clear();
		}

		if (interest_status[front] == 1)
		{
			branch.push_back(front);
			for (auto ff : mesh.ff_range(mesh.face_handle(front)))
			{
				if ((visit_status[ff.idx()] = -1))
				{
					if (interest_status[ff.idx()] == true)
					{
						lables[ff.idx()] = true;
						v_q.push_front(ff.idx());
					}
					else
					{
						lables[ff.idx()] = false;
						v_q.push_back(ff.idx());
					}
				}
			}
		}
		else
		{
			for (auto ff : mesh.ff_range(mesh.face_handle(front)))
			{
				if ((visit_status[ff.idx()] = -1))
				{
					if (interest_status[ff.idx()] == true)
					{
						lables[ff.idx()] = true;
						v_q.push_front(ff.idx());
					}
					else
					{
						lables[ff.idx()] = false;
						v_q.push_back(ff.idx());
					}
				}
			}
		}
	}

	for (int i = 0; i < branchs.size(); i++)
	{
		std::cout << "第 " << i << "分支：" << std::endl;
		for (int j = 0; j < branchs[i].size(); j++)
		{
			std::cout << branchs[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

vector<vector<int>> GraphAlgorithm::search_connect_branchs(std::vector<int>& interest_area_point)
{
	//1.计算拟合区域的信息：即子图
	vector<EdgeInfo> edgeinfo;
	edgeinfo.clear();
	for (int i = 0; i < mesh.n_edges(); i++)
	{
		Mesh::EdgeHandle e_h = mesh.edge_handle(i);
		int id1 = mesh.from_vertex_handle(mesh.halfedge_handle(mesh.edge_handle(i), 0)).idx();
		int id2 = mesh.to_vertex_handle(mesh.halfedge_handle(mesh.edge_handle(i), 0)).idx();

		vector<int>::iterator iter1 = std::find(interest_area_point.begin(), interest_area_point.end(), id1);
		vector<int>::iterator iter2 = std::find(interest_area_point.begin(), interest_area_point.end(), id2);
		vector<int>::iterator iter_end = interest_area_point.end();

		if (iter1 != iter_end && iter2 != iter_end)
		{
			EdgeInfo einfo;
			einfo.e_idx = i;
			einfo.v_idx_i = id1;
			einfo.v_idx_j = id2;
			edgeinfo.push_back(einfo);
		}
	}

	// 2. 构建图：已验证正确
	Graph<int> graph;
	//2.1添加顶点
	for (int i = 0; i < interest_area_point.size(); i++)
	{
		graph.add_vertex(interest_area_point[i]);
	}

	//2.2添加顶点之间的连接关系即边（默认权重1）  graph.add_edge(v1.idx, v2.idx, weight)
	for (int i = 0; i < edgeinfo.size(); i++)
	{
		graph.add_edge(edgeinfo[i].v_idx_i, edgeinfo[i].v_idx_j,1);
	}

	// 3.计算连通分支
	vector<vector<int>> connected_components = graph.get_connected_components();
	return connected_components;
	//graph.print_connected_components(connected_components);
}

std::vector<int> GraphAlgorithm::face2vertex(std::vector<int>& interest_area_face)
{
	vector<int> point_area; point_area.clear();
	vertex_temp.resize(mesh.n_vertices());
	for (int i = 0; i < interest_area_face.size(); i++)
	{
		for (auto fv : mesh.fv_range(mesh.face_handle(interest_area_face[i])))
		{
			vertex_temp[fv.idx()].v_idx = fv.idx();
			vertex_temp[fv.idx()].f_idx_temp.push_back(interest_area_face[i]);
			point_area.push_back(fv.idx());
		}
	}

	std::sort(point_area.begin(), point_area.end());
	auto last = std::unique(point_area.begin(), point_area.end());
	point_area.erase(last, point_area.end());
	return point_area;
}

std::vector<int> GraphAlgorithm::vertex2face(std::vector<int>& interest_area_point)
{
	vector<int> face_area; face_area.clear();
	for (int i=0; i < interest_area_point.size(); i++)
	{
		for (int j = 0; j < vertex_temp[interest_area_point[i]].f_idx_temp.size(); j++)
		{
			face_area.push_back(vertex_temp[interest_area_point[i]].f_idx_temp[j]);
		}
	}
	std::sort(face_area.begin(), face_area.end());
	auto last = std::unique(face_area.begin(), face_area.end());
	face_area.erase(last, face_area.end());
	return face_area;
}

void GraphAlgorithm::delete_small_patches(std::vector<int>& interest_area_face)
{
	vector<int> point_area = face2vertex(interest_area_face);
	vector<vector<int>> connected_components = search_connect_branchs(point_area);
	
	size_t max_size = 0;
	size_t max_index = 0;
	// 遍历 vector<vector<int>> 中的每个元素
	for (size_t i = 0; i < connected_components.size(); ++i) {
		if (connected_components[i].size() > max_size) {
			// 如果找到更大的尺寸，更新最大尺寸和对应的索引值
			max_size = connected_components[i].size();
			max_index = i;
		}
	}

	vector<int> max_branch = connected_components[max_index];
	interest_area_face = vertex2face(max_branch);
	//std::vector<int> face_area = vertex2face(max_branch);

	/*
	for (auto aa: face_area)
	{
		cout << aa << endl;
	}*/
}
