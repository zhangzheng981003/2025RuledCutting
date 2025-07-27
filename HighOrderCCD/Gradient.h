#ifndef GRADIENT_H
#define GRADIENT_H

#include "Config/Config.h"

#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/Distance_der.h"

#include <vector>

#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>

PRJ_BEGIN

class Gradient
{
  public:
    typedef Eigen::MatrixXd Data;
    typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;

    typedef Eigen::Matrix<double,Eigen::Dynamic,1> inner_derivative_t;//3*(order_num+1)
    typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
    typedef Eigen::Matrix<inner_scalar_t,Eigen::Dynamic,1> derivative_t;
    typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
    typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1> Vec12;
    typedef Eigen::Matrix<scalar_t,1,3> Vec3;

    static void gcl_fast_barrier_gradient(const Data& spline,
        Eigen::VectorXd& grad, Eigen::MatrixXd& hessian,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {

        int num = 3 * spline.rows();

        grad.resize(num);
        grad.setZero();

        hessian.resize(num, num);
        hessian.setZero();

        //Eigen::VectorXd auto_grad=grad;
        //Eigen::MatrixXd auto_hessian=hessian;


        bvh.gcl_InitTrajectory(spline);

        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;

        //clock_t time1 = clock();
        bvh.CheckCollision(collision_pair, offset + margin);
        //clock_t time2 = clock();
        //std::cout << "time_bvh:" << (time2 - time1) / (CLOCKS_PER_SEC / 1000) << std::endl << std::endl;
        //std::cout << "bvh_size:" << collision_pair.size() << std::endl;
        
        double dmin = 1.0;

        std::vector<std::vector<std::pair< int, int>>> segment_lists;

        std::vector<std::vector<int>> segment_ob_lists;

        segment_lists.resize(subdivide_tree.size());

        segment_ob_lists.resize(subdivide_tree.size());

        std::vector<std::vector<std::vector<int>>> face_lists;


        std::vector<std::vector<int>> face_ob_lists;


        face_lists.resize(subdivide_tree.size());

        face_ob_lists.resize(subdivide_tree.size());
        
        for (unsigned int i = 0; i < collision_pair.size(); i++)
        {
            
            int tr_id = collision_pair[i].first;
            int ob_id = collision_pair[i].second;

            int sp_id = std::get<0>(subdivide_tree[tr_id]);
            //double weight = std::get<1>(subdivide_tree[tr_id]).second - std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;
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


            for (int j0 = 0; j0 < basis.rows() - 2; j0++)
            {
        
                for (int j1 = j0 + 1; j1 < basis.rows() - 1; j1++)
                {
                    for (int j2 = j1 + 1; j2 < basis.rows(); j2++)
                    {
                        for (int k = 0; k < 3; k++)
                        {
                            Distance<double, 3>::face_point(P[j0], P[j1], P[j2],
                                S[k], d, C);
                            //d -= offset;
                            
                            if (d < margin)
                            {
                                if (d < dmin)
                                    dmin = d;
                                std::vector<int> face_list; face_list.clear();
                                face_list.push_back(j0);
                                face_list.push_back(j1);
                                face_list.push_back(j2);
                                face_list.push_back(k);
                                
                                face_lists[tr_id].push_back(face_list);

                                face_ob_lists[tr_id].push_back(ob_id);
                            }
                        }
                    }
                }
            }
            
        }

        for (unsigned int tr_id = 0; tr_id < subdivide_tree.size(); tr_id++)
        {
            int sp_id = std::get<0>(subdivide_tree[tr_id]);
            //double weight = std::get<1>(subdivide_tree[tr_id]).second - std::get<1>(subdivide_tree[tr_id]).first;
            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[tr_id]);

            //std::cout<<"sp_id:"<<sp_id<<std::endl;

            Eigen::MatrixXd x = spline;
            Vec12 X;
            X.resize(3 * spline.rows());
            for (int r = 0; r < spline.rows(); r++)
            {
                X(3 * r).value() = x(r, 0);
                X(3 * r + 1).value() = x(r, 1);
                X(3 * r + 2).value() = x(r, 2);
            }
            //repeat partial derivatives for the inner AutoDiffScalar
            for (int id = 0; id < 3 * spline.rows(); id++)
            {
                X(id).derivatives().resize(3 * spline.rows());
                X(id).derivatives().setZero();
                X(id).derivatives()(id) = 1;
                X(id).value().derivatives() = inner_derivative_t::Unit(3 * spline.rows(), id);
            }
            //set the hessian matrix to zero
            
            for (int idx = 0; idx < 3 * spline.rows(); idx++) 
            {
                for (int id = 0; id < 3 * spline.rows(); id++)
                {
                    X(id).derivatives()(idx).derivatives() = inner_derivative_t::Zero(3 * spline.rows());
                }
            }

            std::vector<Vec3> P(basis.rows());

            for (int j = 0; j < basis.rows(); j++)
            {
                P[j].setZero();
                for (int j0 = 0; j0 < spline.rows(); j0++)
                {
                    P[j](0) += basis(j, j0) * X(3 * j0);
                    P[j](1) += basis(j, j0) * X(3 * j0 + 1);
                    P[j](2) += basis(j, j0) * X(3 * j0 + 2);
                }
            }

            Vec3 C, C0, C1;
            scalar_t d;
            /*
            std::vector<std::pair< int, int>> segment_list = segment_lists[tr_id];
            
            //segment segment
            for (unsigned int j = 0; j < segment_list.size(); j++)
            {cout << tr_id << endl;
                int ob_id = segment_ob_lists[tr_id][j];
                std::vector<Eigen::RowVector3d> S(3);
                for (int i = 0; i < 3; i++)
                {
                    S[i] = V.row(F(ob_id, i));
                }
                int j0 = segment_list[j].first;
                int k = segment_list[j].second;
                Distance<scalar_t, 3>::segment_segment(P[j0], P[(j0 + 1) % ((udegree + 1) * (vdegree + 1))],
                    S[k], S[(k + 1) % 3], d, C0, C1);
                d = d - offset;

                scalar_t e = -(d - margin) * (d - margin) / d * log(d / margin);///d

                grad+= e.value().derivatives();
                cout <<"grad        :"<< grad << endl;
                Eigen::MatrixXd B;
                B.resize(3* (udegree + 1) * (vdegree + 1), 3 * (udegree + 1) * (vdegree + 1));
                for (int r = 0; r < 3 * (udegree + 1) * (vdegree + 1); r++)
                {
                    B.row(r) = e.derivatives()(r).derivatives().transpose();
                }
                hessian+= B;
            }
            */

            std::vector<std::vector<int>> face_list = face_lists[tr_id];
           // cout << "listsize2:" << face_list.size() << endl;
            for (unsigned int j = 0; j < face_list.size(); j++)
            {
                int ob_id = face_ob_lists[tr_id][j];
                std::vector<Eigen::RowVector3d> S(3);
                for (int i = 0; i < 3; i++)
                {
                    S[i] = V.row(F(ob_id, i));
                }
                int j0 = face_list[j][0];
                int j1 = face_list[j][1];
                int j2 = face_list[j][2];
                int k = face_list[j][3];

                Distance<scalar_t, 3>::face_point(P[j0], P[j1], P[j2],
                    S[k], d, C);
               // d = d - offset;
                //scalar_t e=-weight*(d-margin)*(d-margin)*log(d/margin);

                scalar_t e = - (d - margin) * (d - margin) / d * log(d / margin);///d

                //scalar_t e=weight*(1-d/margin*d/margin)*(1-d/margin*d/margin);
                grad += e.value().derivatives();

                Eigen::MatrixXd B;
                B.resize(3 * spline.rows(), 3 * spline.rows());
                for (int r = 0; r < 3 * spline.rows(); r++)
                {
                    B.row(r) = e.derivatives()(r).derivatives().transpose();
                }
                hessian += B;
            }

        }

        std::cout << "dmin:" << dmin << std::endl;
    }
    
};

PRJ_END

#endif