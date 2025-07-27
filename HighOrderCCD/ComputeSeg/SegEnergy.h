#ifndef SEGENERGY_H
#define SEGENERGY_H

#include "Optimization.h"
#include "CutMesh.h"
#include "HighOrderCCD/BSplineFitting.h"

class SegEnergy
{
public:
	struct InnerEdge
	{
		int inner_edge_idx;
		int f0_id;
		int f1_id;
		set<int> Adj_v;
		double vertices_dist;
		Vector3d rulling_f0_ori;
		Vector3d rulling_f1_ori;
		Vector3d rulling_f0_bspline;
		Vector3d rulling_f1_bspline;
	};

    struct AdjFace
    {
        int face_idx;
        Vector3d rulling_f0_ori;
        Vector3d rulling_f0_bspline;
    };
    
    
    static void calc_inner_edge_info(Mesh& mesh0, Mesh& seg_mesh, vector<InnerEdge>& inner_edge_temp, vector<AdjFace>& adj_face_temp, vector<vector<double>>& res, vector<int>& real_in_txt,
		Handle(Geom_BSplineSurface) BSpline_surface, BSplineFitting& sf,vector<int>& img2real, int curve_ctr_num, std::vector<vector<Vector3d>>& P)
    {
		real_in_txt.clear(); real_in_txt.resize(mesh0.n_faces());
		for (int i = 0; i < res.size(); i++)
		{
			real_in_txt[res[i][0]] = i;
		}

		adj_face_temp.clear(); adj_face_temp.resize(seg_mesh.n_faces());
		for (int i = 0; i < seg_mesh.n_faces(); i++)
		{
			SegEnergy::AdjFace af;
			af.face_idx = i;
			af.rulling_f0_ori = {res[real_in_txt[img2real[i]]][1],res[real_in_txt[img2real[i]]][2],res[real_in_txt[img2real[i]]][3] };
			FH fh = seg_mesh.face_handle(i);
			Vec3d p_ = seg_mesh.calc_face_centroid(fh);
			gp_Pnt pointToProject(p_[0], p_[1], p_[2]);
			GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, BSpline_surface);
			if (projectPoint.NbPoints() > 0)
			{
				Standard_Real u, v;
				projectPoint.LowerDistanceParameters(u, v);
				Standard_Real dist;
				dist = projectPoint.LowerDistance();
				int u_interval_order = sf.u_k_position(u, sf.u_knot_temp);

				vector<double> u_basis_temp;
				if (u_interval_order == 0)
				{
					vector<double> basis_temp(sf.poly_degree + 1, 0);
					basis_temp[0] = 1;
					u_basis_temp = basis_temp;
				}
				else if (u_interval_order == sf.u_knot_temp.size() - 1)
				{
					vector<double> basis_temp(sf.poly_degree + 1, 0);
					basis_temp[sf.poly_degree] = 1;
					u_basis_temp = basis_temp;
				}
				else
				{
					sf.calc_basis_fuction(u, u_interval_order, sf.poly_degree, u_basis_temp);
				}
				Vector3d rulling;
				rulling.setZero();
				if (u_interval_order == 0)
				{
					Vector3d pi0, pi1;
					pi0 = P[0][0];
					pi1 = P[0][1];
					rulling = pi1 - pi0;
				}
				else if (u_interval_order == sf.u_knot_temp.size() - 1)
				{
					Vector3d pi0, pi1;
					pi0 = P[curve_ctr_num - 1][0];
					pi1 = P[curve_ctr_num - 1][1];
					rulling = pi1 - pi0;
				}
				else
				{
					for (int i = 0; i < u_basis_temp.size(); i++)
					{
						Vector3d pi0, pi1;
						pi0 = P[u_interval_order - sf.poly_degree + i][0];
						pi1 = P[u_interval_order - sf.poly_degree + i][1];
						rulling += u_basis_temp[i] * (pi1 - pi0);
					}
				}
				af.rulling_f0_bspline = rulling;
			}
			adj_face_temp[i] = af;
		}

		inner_edge_temp.clear(); //ie_temp.resize(mesh.n_edges());
		for (int k0 = 0; k0 < seg_mesh.n_edges(); k0++)
		{
			EH eh = seg_mesh.edge_handle(k0);
			if (!seg_mesh.is_boundary(eh))
			{
				SegEnergy::InnerEdge ie;
				ie.inner_edge_idx = eh.idx();
				HEH heh0 = seg_mesh.halfedge_handle(eh, 0);
				HEH heh1 = seg_mesh.halfedge_handle(eh, 1);

				FH f0 = seg_mesh.face_handle(heh0);
				FH f1 = seg_mesh.face_handle(heh1);

				ie.f0_id = f0.idx();
				ie.rulling_f0_ori = adj_face_temp[ie.f0_id].rulling_f0_ori;
				ie.rulling_f0_bspline = adj_face_temp[ie.f0_id].rulling_f0_bspline;

				ie.f1_id = f1.idx();
				ie.rulling_f1_ori = adj_face_temp[ie.f1_id].rulling_f0_ori;
				ie.rulling_f1_bspline = adj_face_temp[ie.f1_id].rulling_f0_bspline;

				set<int> v_id_temp; v_id_temp.clear();
				for (auto vv : seg_mesh.fv_range(f0))
				{
					v_id_temp.insert(vv.idx());
				}
				for (auto vv : seg_mesh.fv_range(f1))
				{
					v_id_temp.insert(vv.idx());
				}
				ie.Adj_v = v_id_temp;

				double total_dist = 0.0;
				for (auto ele : v_id_temp)
				{
					Vec3d p0 = seg_mesh.point(seg_mesh.vertex_handle(ele));
					gp_Pnt pointToProject(p0[0], p0[1], p0[2]);
					GeomAPI_ProjectPointOnSurf projectPoint(pointToProject, BSpline_surface);
					if (projectPoint.NbPoints() > 0)
					{
						Standard_Real u, v;
						projectPoint.LowerDistanceParameters(u, v);
						Standard_Real dist;
						dist = projectPoint.LowerDistance();
						total_dist += dist;
					}
				}
				ie.vertices_dist = total_dist;
				inner_edge_temp.push_back(ie);
			}
			else
			{
				SegEnergy::InnerEdge ie;
				inner_edge_temp.push_back(ie);
			}
		}
    }

    static Vector4d calc_quaternion(Vector3d& v0, Vector3d& v1)
    {
      
        Vector4d quaternion = { 0,0,0,0 };
        Vector3d nv0 = v0.normalized();
        Vector3d nv1 = v1.normalized();
      
        if ((nv0 + nv1).norm() == 0)
        {
            quaternion = { 0,1,0,0 };
        }
        else
        {
            Vector3d half = (nv0 + nv1).normalized();
            quaternion[0] = nv0.dot(half);
            quaternion[1] = nv0.cross(half)[0];
            quaternion[2] = nv0.cross(half)[1];
            quaternion[3] = nv0.cross(half)[2];
    
        }
        quaternion.normalize();
        return quaternion;
    }
    
    
    static void calc_edge_energy0(Mesh& mesh, double weight, vector<double>& rullings_energy,vector<InnerEdge>& inner_edge_temp)
    {
        
        rullings_energy.clear();
		rullings_energy.resize(mesh.n_edges(), 0);

		for (auto ie : inner_edge_temp)
		{
            //1.rulings项
            Vector4d quaternion_ori = calc_quaternion(ie.rulling_f0_ori, ie.rulling_f1_ori);
            Vector4d quaternion_bspline = calc_quaternion(ie.rulling_f0_bspline, ie.rulling_f1_bspline);
            
            if (quaternion_ori[0] < 0)
            {
                quaternion_ori = -quaternion_ori;
            }
            if (quaternion_bspline[0] < 0)
            {
                quaternion_bspline = -quaternion_bspline;
            }

            //2.距离项
            double dist_term = ie.vertices_dist;

            //3.总能量
            rullings_energy[ie.inner_edge_idx] = (quaternion_ori - quaternion_bspline).norm() + weight * dist_term;
		}
    }

    static double calc_edge_energy1(Mesh& mesh, vector<double>& rullings_energy, vector<vector<double>>& res, vector<int>& real2img)
    {
        vector<Vector3d> face_rulings; face_rulings.resize(mesh.n_faces());
        for (int i = 0; i < mesh.n_faces(); i++)
        {
            face_rulings[real2img[res[i][0]]] = { res[i][1],res[i][2],res[i][3] };
        }

        rullings_energy.clear();
        rullings_energy.resize(mesh.n_edges());
        for (auto eh : mesh.edges())
        {
            if (mesh.is_boundary(eh))
            {
                rullings_energy[eh.idx()] = 0;
            }
            else
            {
                FH f0 = mesh.face_handle(mesh.halfedge_handle(eh, 0));
                FH f1 = mesh.face_handle(mesh.halfedge_handle(eh, 1));
                double dist = (face_rulings[f0.idx()] - face_rulings[f1.idx()]).norm();
                rullings_energy[eh.idx()] = dist;
            }
        }
    }
};
#endif