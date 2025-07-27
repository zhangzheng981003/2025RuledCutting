#ifndef ENERGY_H
#define ENERGY_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/Distance_der.h"

#include <vector>

#include <HighOrderCCD/BSplineFitting.h>

PRJ_BEGIN

class Energy
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;
    
    //@brief 0927 start.
    /*
    static double gcl_whole_energy(const Data& spline, const double& piece_time,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
        return data_weight * gcl_data_energy(spline, piece_time) +
               smooth_weight * gcl_smooth_energy(spline, piece_time) +
               barrier_weight * gcl_barrier_energy(spline,V, F,bvh);
    }*/

    static double gcl_fast_whole_energy(const Data& spline, 
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh, BSplineFitting& ssf)
    {

       double data_value = ssf.calc_data_term_energy(spline);
       double smooth_value = ssf.calc_smooth_term_energyn(spline);
       return data_weight * data_value +
           smooth_weight * smooth_value;
          //+ fast_barrier_weight * gcl_fast_barrier_energy(spline, V, F, bvh);
    }

    //@brief 0927 end.
    //@brief 0927 start.
    static double gcl_data_energy(const Data& spline, int udegree, int vdegree)
    {
        double energy = 0;
        return energy;
    }

    static double gcl_smooth_energy(const Data& spline, int udegree, int vdegree)
    {
        double energy = 0;
        return energy;
    }

    static double gcl_fast_barrier_energy(const Data& spline,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
        double energy = 0;
        bvh.gcl_InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);
        
        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CheckCollision(collision_pair, offset + margin);
        //cout << "collisonnum" << collision_pair.size() << endl;
        for (unsigned int i = 0; i < collision_pair.size(); i++)
        {
            int tr_id = collision_pair[i].first;
            int ob_id = collision_pair[i].second;

            int sp_id = std::get<0>(subdivide_tree[tr_id]);
            //double weight = std::get<1>(subdivide_tree[tr_id]).second - std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[tr_id]);

            Eigen::MatrixXd bz;
            bz = spline;

            std::vector<Eigen::RowVector3d> P(basis.rows());
            for (int j = 0; j < basis.rows(); j++)
            {
                P[j].setZero();
                for (int j0 = 0; j0 < spline.rows(); j0++)
                {
                    P[j] += basis(j, j0) * bz.row(j0);
                }

            }

            std::vector<Eigen::RowVector3d> S(3);
            for (int j = 0; j < 3; j++)
            {
                S[j] = V.row(F(ob_id, j));

            }

            Eigen::RowVector3d C, C0, C1;
            double d;

            //point face

            for (int j = 0; j < basis.rows(); j++)
            {
                Distance<double, 3>::point_face(P[j], S[0], S[1], S[2], d, C);
                d -= offset;
                if (d <= 0)
                    return INFINITY;

                if (d < margin)
                {
                    energy += -(d - margin) * (d - margin) / d * log(d / margin); ///d

                    //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                }
            }

            //segment segment
            //for(int j=0;j<=order_num;j++)
            {
                int j = basis.rows() -1;
                for (int k = 0; k < 3; k++)
                {
                    int k1 = (k + 1) % 3;
                    Distance<double, 3>::segment_segment(P[j], P[(j + 1) % (basis.rows())],
                        S[k], S[k1], d, C0, C1);
                    d -= offset;
                    if (d <= 0)
                        return INFINITY;
                    if (d < margin)
                    {
                        energy += - (d - margin) * (d - margin) / d * log(d / margin);///d

                        //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                    }
                }
            }

        }

        return energy;
    }
    //@brief 0927 end.
    static double gcl_barrier_energy(const Data& spline,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
        double energy = 0;
        bvh.gcl_InitTrajectory(spline);

        //bvh.UpdateTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        //bvh.CheckCollision(collision_pair,margin);
        bvh.CheckCollision(collision_pair, offset + margin);
        
        for (unsigned int i = 0; i < collision_pair.size(); i++)
        {
            int tr_id = collision_pair[i].first;
            int ob_id = collision_pair[i].second;

            int sp_id = std::get<0>(subdivide_tree[tr_id]);
            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[tr_id]);

            Eigen::MatrixXd bz;
            bz = spline;

            std::vector<Eigen::RowVector3d> P(basis.rows());
            for (int j = 0; j < basis.rows(); j++)
            {
                P[j].setZero();
                for (int j0 = 0; j0 < spline.rows(); j0++)
                {
                    P[j] += basis(j, j0) * bz.row(j0);
                }
            }

            std::vector<Eigen::RowVector3d> S(3);
            for (int j = 0; j < 3; j++)
            {
                S[j] = V.row(F(ob_id, j));

            }

            Eigen::RowVector3d C, C0, C1;
            double d;

            //point face
            /*
            for (int j = 0; j < (udegree + 1) * (vdegree + 1); j++)
            {
                Distance<double, 3>::point_face(P[j], S[0], S[1], S[2], d, C);
                d -= offset;
                if (d <= 0)
                {
                    return INFINITY;
                }


                if (d < margin)
                {
                    energy += -(d - margin) * (d - margin) / d * log(d / margin); ///d

                    //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                }
            }
            */
            //face point

            for (int j0 = 0; j0 < basis.rows() - 2; j0++)
            {
                for (int j1 = j0 + 1; j1 < basis.rows() -1; j1++)
                {
                    for (int j2 = j1 + 1; j2 < basis.rows(); j2++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            Distance<double, 3>::face_point(P[j0], P[j1], P[j2],
                                S[k], d, C);
                            d -= offset;
                            if (d <= 0)
                                return INFINITY;
                            if (d < margin)
                            {
                                energy += -(d - margin) * (d - margin) / d * log(d / margin); ///d

                                //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                            }
                        }
                    }
                }
            }


            //segment segment
            /*
            for (int j0 = 0; j0 < (udegree + 1) * (vdegree + 1)-1; j0++)
            {
                for (int j1 = j0 + 1; j1 < (udegree + 1) * (vdegree + 1); j1++)
                {
                    for (int k = 0; k < 3; k++)
                    {
                        int k1 = (k + 1) % 3;
                        Distance<double, 3>::segment_segment(P[j0], P[j1],
                            S[k], S[k1], d, C0, C1);
                        d -= offset;
                        if (d <= 0)
                        {
                            return INFINITY;
                        }

                        if (d < margin)
                        {
                            energy += -(d - margin) * (d - margin) / d * log(d / margin); ///d
                            //cout << "energy:" << energy << endl << endl << endl;
                            //energy+=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);   
                        }
                    }
                }
            }
            */
        }

        return energy;
    }

};

PRJ_END

#endif