#ifndef SUBDIVIDE_H
#define SUBDIVIDE_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
//#include <HighOrderCCD/SplineFitting/BezierSubvision.h>
//#include <HighOrderCCD/BSplineFitting.h>
#include "HighOrderCCD/BSplineSubdivide.h"
#include <vector>

PRJ_BEGIN

class Subdivide
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;

    static double gcl_find_farthest_dist2(const Data& position)
    {
        double maxDistance = 0.0;
        for (size_t i = 0; i < position.rows(); ++i)
        {
            for (size_t j = i + 1; j < position.rows(); ++j) 
            {
                double dist = (position.coeff(i, 0) - position.coeff(j, 0)) * (position.coeff(i, 0) - position.coeff(j, 0))
                    + (position.coeff(i, 1) - position.coeff(j, 1)) * (position.coeff(i, 1) - position.coeff(j, 1)) +
                    (position.coeff(i, 2) - position.coeff(j, 2)) * (position.coeff(i, 2) - position.coeff(j, 2));
                
                maxDistance = std::max(dist, maxDistance);
            }
        }
        return maxDistance;
    }

    static void gcl_update_tree(const Data& spline,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
       // flag_subdivision:  
        bvh.gcl_InitTrajectory(spline);
        std::vector<std::pair<unsigned int, unsigned int>> collision_pair;
        bvh.CheckCollision(collision_pair, offset + margin);

        //bool is_subdivision = false;
        for (auto it = collision_pair.begin(); it != collision_pair.end(); ++it)
        {
            int tr_id = (*it).first;
            int ob_id = (*it).second;

            int sp_id = std::get<0>(subdivide_tree[tr_id]);
            int f0 = F(ob_id, 0); int f1 = F(ob_id, 1); int f2 = F(ob_id, 2);
            Eigen::Matrix3d _position;
            _position << V.row(f0), V.row(f1), V.row(f2);

            Eigen::MatrixXd bz;
            bz = spline;
            Eigen::MatrixXd x;
            x = bz;

            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[tr_id]);;//线性变换
            Eigen::MatrixXd temp = basis * x;//新控制点
            
            //double diam2 = (temp.row(0) - temp.row(temp.rows()-1)).squaredNorm();//计算直径
            //double diam2 = gcl_find_farthest_dist2(temp);
            //if (CCD::gcl_KDOPDCD(temp, _position, offset + margin))
            {
                cout << "divide_success" << endl;
                if (subdivision_number<13)
                {
                    BSplineSubdivide::divide(3, divide_ctrpoint_, divide_u_knot_, divide_basis_);
                    subdivide_tree.clear();
                    subdivide_tree.resize(divide_basis_[subdivision_number].size());
                    //int stn = 0;
                    for (int i = 0; i < divide_basis_[subdivision_number].size(); i++)
                    {
                        Eigen::MatrixXd basis = divide_basis_[subdivision_number][i];
                        Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
                        combined.setZero();
                        combined.block(0, 0, basis.rows(), basis.cols()) = basis;
                        combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
                        std::pair<double, double> range(1, 1);
                        subdivide_tree[i] = std::make_tuple(i, range, combined);
                    }

                    break;
                }
            }
        }

        //flag_end:
    }

    static void gcl_update_tree2(const Data& spline,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
        if (subdivision_number < 13)
        {
            BSplineSubdivide::divide(3, divide_ctrpoint_, divide_u_knot_, divide_basis_);
            subdivide_tree.clear();
            subdivide_tree.resize(divide_basis_[subdivision_number].size());
            //int stn = 0;
            for (int i = 0; i < divide_basis_[subdivision_number].size(); i++)
            {
                Eigen::MatrixXd basis = divide_basis_[subdivision_number][i];
                Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
                combined.setZero();
                combined.block(0, 0, basis.rows(), basis.cols()) = basis;
                combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
                std::pair<double, double> range(1, 1);
                subdivide_tree[i] = std::make_tuple(i, range, combined);
            }
        }
    }


};

PRJ_END

#endif