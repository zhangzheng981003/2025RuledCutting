#ifndef OPTIMIZATION3D_TIME_H
#define OPTIMIZATION3D_TIME_H

#include "HighOrderCCD/Config/Config.h"

#include "HighOrderCCD/Distance.h"
#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"

#include "HighOrderCCD/Subdivide.h"
#include "HighOrderCCD/Energy.h"
#include "HighOrderCCD/Gradient.h"
#include "HighOrderCCD/Step.h"

#include <vector>
#include <ctime>
#include <Eigen/SparseCholesky>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/KroneckerProduct>

#include <HighOrderCCD/BSplineFitting.h>

//#include "gurobi_c++.h"
PRJ_BEGIN

class Optimization3D_time
{
public:

    typedef Eigen::MatrixXd Data;
    //Eigen::Dynamic
    typedef Eigen::Matrix<double, Eigen::Dynamic, 1> inner_derivative_t;//3*(order_num+1)
    typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
    typedef Eigen::Matrix<inner_scalar_t, Eigen::Dynamic, 1> derivative_t;
    typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
    typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> Vec12;
    typedef Eigen::Matrix<scalar_t, 1, 3> Vec3;

    typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

    typedef std::vector< std::tuple< int, std::pair<double, double>, Eigen::MatrixXd > > Tree;
    typedef std::tuple< int, std::pair<double, double>, Eigen::MatrixXd >  Node;


    static void gcl_optimization(Data& spline,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh, BSplineFitting& ssf)
    {
        cout << "stop_condition=" << stop_condition << endl;
        if (max_step < 0)
        {
            ssf.zero_num++;
        }
        else if (max_step > 1e-4)
        {
            ssf.zero_num = 0;
        }


        if (subdivision_number < 13)
        {



            Data direction;
            if (ssf.iter_num > 0 && (ssf.insert1_idx!=-1))//|| stop_condition < 1e-3)//max_step < 1e-6 || 
            {
                cout << "增加控制点" << endl;

                Eigen::VectorXd bspline;
                bspline.resize(3 * spline.rows());

                for (int i = 0; i < spline.rows(); i++)
                {
                    bspline[3 * i] = spline.coeffRef(i, 0);
                    bspline[3 * i + 1] = spline.coeffRef(i, 1);
                    bspline[3 * i + 2] = spline.coeffRef(i, 2);
                }




                if (data_e > 1 && (bspline.size() / 6) < 0)
                {
                    double u_loc = ssf.u_knot_temp[ssf.insert1_idx] * 1 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 5 / 6;
                    double u_loc1 = ssf.u_knot_temp[ssf.insert1_idx] * 2 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 4 / 6;
                    double u_loc2 = ssf.u_knot_temp[ssf.insert1_idx] * 3 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 3 / 6;
                    double u_loc3 = ssf.u_knot_temp[ssf.insert1_idx] * 4 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 2 / 6;
                    double u_loc4 = ssf.u_knot_temp[ssf.insert1_idx] * 5 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 1 / 6;

                    cout << "插入点位置:" << u_loc << " " << u_loc1 << " " << u_loc2 << " " << u_loc3 << " " << u_loc4 << endl;
                    if (ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])] < 8)
                    {
                        ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])]++;
                        ssf.add_ctr_point(bspline, u_loc);
                        ssf.add_ctr_point(bspline, u_loc1);
                        ssf.add_ctr_point(bspline, u_loc2);
                        ssf.add_ctr_point(bspline, u_loc3);
                        ssf.add_ctr_point(bspline, u_loc4);
                    }
                    //ssf.insert1_idx++;
                    //.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])]++;
                    ssf.is_change.clear();
                    ssf.is_change.resize(bspline.size() / 3);
                    for (int i = 0; i < ssf.is_change.size(); i++)
                    {
                        ssf.is_change[i] = 1;
                    }
                    ssf.is_add = 1;
                }
                else if (bspline.size() / 6 < 300)
                {
                    double u_loc = ssf.u_knot_temp[ssf.insert1_idx] * 1 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 5 / 6;
                    double u_loc1 = ssf.u_knot_temp[ssf.insert1_idx] * 2 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 4 / 6;
                    double u_loc2 = ssf.u_knot_temp[ssf.insert1_idx] * 3 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 3 / 6;
                    double u_loc3 = ssf.u_knot_temp[ssf.insert1_idx] * 4 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 2 / 6;
                    double u_loc4 = ssf.u_knot_temp[ssf.insert1_idx] * 5 / 6 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 1 / 6;

                    cout << "avr_e" << ssf.avr_e_intervals[ssf.insert1_idx] << endl;
                    if ((ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])] < 1 ||(ssf.e_intervals[ssf.insert1_idx]>0.5||ssf.avr_e_intervals[ssf.insert1_idx]>1e-3))&&ssf.e_len[ssf.insert1_idx]>0.1)
                    {
                        cout << "插入点位置:" << u_loc << " " << u_loc1 << " " << u_loc2 << " " << u_loc3 << " " << u_loc4 << endl;
                       // ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])]++;
                        ssf.add_ctr_point(bspline, u_loc);
                        ssf.add_ctr_point(bspline, u_loc1);
                        ssf.add_ctr_point(bspline, u_loc2);
                        ssf.add_ctr_point(bspline, u_loc3);
                        ssf.add_ctr_point(bspline, u_loc4);
                    }
                    
                    else if (ssf.e_intervals[ssf.insert1_idx] <= 0.5 && (ssf.e_intervals[ssf.insert1_idx] > 0.1 || ssf.avr_e_intervals[ssf.insert1_idx] > 1e-3)&& ssf.e_len[ssf.insert1_idx] > 0.05)
                    {
                        u_loc1 = ssf.u_knot_temp[ssf.insert1_idx] * 1 / 4 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 3 / 4;
                        u_loc2 = ssf.u_knot_temp[ssf.insert1_idx] * 2 / 4 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 2 / 4;
                        u_loc3 = ssf.u_knot_temp[ssf.insert1_idx] * 3 / 4 + ssf.u_knot_temp[ssf.insert1_idx + 1] * 1 / 4;
                        cout << "插入点位置:" << u_loc1 << " " << u_loc2 << " " << u_loc3 << endl;
                        //ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])]++;
                        ssf.add_ctr_point(bspline, u_loc1);
                        ssf.add_ctr_point(bspline, u_loc2);
                        ssf.add_ctr_point(bspline, u_loc3);
                    }
                    else if(ssf.e_len[ssf.insert1_idx] > 0.02)
                    {
                        cout << "插入点位置:" << u_loc2 << endl;
                        
                        ssf.add_ctr_point(bspline, u_loc2);
                    }
                    ssf.count_divide[(int)floor(ssf.u_knot_temp[ssf.insert1_idx])]++;
                    ssf.insert1_idx += 1;
                    ssf.is_add = 2;
                }
                ssf.is_change.clear();
                ssf.is_change.resize(bspline.size() / 3);
                for (int i = 0; i < ssf.is_change.size(); i++)
                {
                    ssf.is_change[i] = 0;
                }
                ssf.is_change[ssf.insert1_idx] = 1;
                ssf.is_change[ssf.insert1_idx + bspline.size() / 6] = 1;

                spline.resize(bspline.size() / 3, 3);
                for (int i = 0; i < spline.rows(); i++)
                {
                    spline(i, 0) = bspline[3 * i];
                    spline(i, 1) = bspline[3 * i + 1];
                    spline(i, 2) = bspline[3 * i + 2];
                }

                divide_ctrpoint_[0][0] = bspline;
                divide_u_knot_[0][0] = ssf.u_knot_temp;
                divide_basis_[0][0].resize(bspline.size() / 6, bspline.size() / 6);
                divide_basis_[0][0].setIdentity();

                subdivide_tree.resize(1);
                std::pair<double, double> range(1, 1);
                Eigen::MatrixXd basis = divide_basis_[0][0];
                Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
                combined.setZero();
                combined.block(0, 0, basis.rows(), basis.cols()) = basis;
                combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
                subdivide_tree[0] = std::make_tuple(0, range, combined);
                subdivision_number = 0;


                ssf.BSpline_initialization();
                //ssf.BSpline_surface_viewer_2(spline, 50, -2);
            }
            //cout <<"size:=" << ssf.mesh2surface_temp[0].point_on_mesh[0].value().derivatives().size() << endl;
            cout << "细分次数： " << subdivision_number << endl;
            while (subdivision_number < 8)
            {
                cout << "准备细分" << endl;
                //Subdivide::gcl_update_tree(spline, V, F, bvh);
                Subdivide::gcl_update_tree2(spline, V, F, bvh);
                cout << "细分次数： " << subdivision_number << endl;
            }






            clock_t time0 = clock();
            
            gcl_descent_direction(spline, direction, V, F, bvh, ssf);

            is_good_init = ssf.is_good;


            clock_t time1 = clock();


            gcl_line_search(spline, direction, V, F, bvh, ssf, false);//vtxs

            //cout << direction << endl;

            //update_time(spline, bvh);
            clock_t time2 = clock();
            //Subdivide::gcl_update_tree(spline, V, F, bvh, udegree, vdegree);
            cout << "stop_condition" << stop_condition << endl;
            std::cout << std::endl << "time10:" << (time1 - time0) / (CLOCKS_PER_SEC / 1000) << std::endl << std::endl;
            std::cout << "time21:" << (time2 - time1) / (CLOCKS_PER_SEC / 1000) << std::endl << std::endl;
            Mesh out1;
            for (int i = 0; i < divide_ctrpoint_[8].size(); i++)
            {
                //cout << "basis:=" << divide_basis_[10][i].rows() << ";" << divide_basis_[10][i].cols() << endl;

                Eigen::MatrixXd basis = divide_basis_[8][i];

                Eigen::MatrixXd combined(2 * basis.rows(), 2 * basis.cols());
                combined.setZero();
                combined.block(0, 0, basis.rows(), basis.cols()) = basis;
                combined.block(basis.rows(), basis.cols(), basis.rows(), basis.cols()) = basis;
                //auto ctr = divide_ctrpoint_[10][i];

                Eigen::MatrixXd bctr = combined * spline;
                Eigen::VectorXd ctr;
                ctr.resize(3 * bctr.rows());

                for (int l = 0; l < bctr.rows(); l++)
                {
                    ctr[3 * l] = bctr.coeffRef(l, 0);
                    ctr[3 * l + 1] = bctr.coeffRef(l, 1);
                    ctr[3 * l + 2] = bctr.coeffRef(l, 2);
                }

                //cout << "ctr:" << ctr.size() << endl;
                for (int j = 0; j < ctr.size() / 6 - 1; j++)
                {
                    auto v0 = out1.add_vertex(Mesh::Point(ctr[3 * j], ctr[3 * j + 1], ctr[3 * j + 2]));
                    auto v1 = out1.add_vertex(Mesh::Point(ctr[3 * (j + 1)], ctr[3 * (j + 1) + 1], ctr[3 * (j + 1) + 2]));
                    auto v2 = out1.add_vertex(Mesh::Point(ctr[3 * j + ctr.size() / 2], ctr[3 * j + 1 + ctr.size() / 2], ctr[3 * j + 2 + ctr.size() / 2]));
                    auto v3 = out1.add_vertex(Mesh::Point(ctr[3 * (j + 1) + ctr.size() / 2], ctr[3 * (j + 1) + 1 + ctr.size() / 2], ctr[3 * (j + 1) + 2 + ctr.size() / 2]));
                    std::vector<Mesh::VertexHandle> face_vhandles0;
                    face_vhandles0.push_back(v0);
                    face_vhandles0.push_back(v1);
                    face_vhandles0.push_back(v3);
                    out1.add_face(face_vhandles0);

                    std::vector<Mesh::VertexHandle> face_vhandles1;
                    face_vhandles1.push_back(v0);
                    face_vhandles1.push_back(v3);
                    face_vhandles1.push_back(v2);
                    out1.add_face(face_vhandles1);
                }
            }

            MeshTools::WriteMesh(out1, "./ctr_point" + to_string(ssf.iter_num) + ".obj");
        }
        //else subdivision_number++;
    }

    static void gcl_descent_direction(const Data& spline, Data& direction,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh, BSplineFitting& ssf)
    {
        Eigen::VectorXd grad, g1, g2;
        Eigen::MatrixXd hessian, h1, h2;

        // data and smooth
        scalar_t data_value; Eigen::VectorXd data_grad; Eigen::MatrixXd data_hessian;
        scalar_t smooth_value; Eigen::VectorXd smooth_grad; Eigen::MatrixXd smooth_hessian;
        Eigen::VectorXd total_obgrad; Eigen::MatrixXd total_obhessian; Eigen::VectorXd de_direction;
        ssf.calc_objective_grad_hessien(spline, data_value, data_grad, data_hessian, de_direction, smooth_value, smooth_grad, smooth_hessian, total_obgrad, total_obhessian);

        if (subdivision_number < 13)
        {
            grad = data_weight * data_grad + smooth_weight * smooth_grad;
            hessian = data_weight * data_hessian + smooth_weight * smooth_hessian;
        }


        int non_0_nb = 0;
        vector<int> non_0_idx;

        if (grad.size() / 6 < 12 || 1)
        {
            for (int i = 0; i < grad.size() / 3; i++)
            {
                double temp_ = abs(grad[3 * i]) + abs(grad[3 * i + 1]) + abs(grad[3 * i + 2]);
                if (temp_ > 1.0 / (spline.rows() / 2.0))
                {
                    non_0_nb++;
                    non_0_idx.push_back(i);
                }
            }
        }
        else
        {
            double grad_mean = 0;
            int count1 = 0;
            for (int i = 0; i < grad.size() / 6; i++)
            {
                double temp_ = abs(grad[3 * i] * grad[3 * i] + grad[3 * i + 1] * grad[3 * i + 1] + grad[3 * i + 2] * grad[3 * i + 2]);
                if (ssf.is_change[i] == 1)
                {
                    non_0_nb++;
                    non_0_idx.push_back(i);
                }
            }
        }
        cout << non_0_nb << endl;
        Eigen::MatrixXd sub_hessian;

        sub_hessian.resize(non_0_nb * 3, non_0_nb * 3);
        sub_hessian.setZero();

        Eigen::VectorXd sub_grad;
        sub_grad.resize(non_0_nb * 3);
        sub_grad.setZero();
        for (int i = 0; i < non_0_nb; i++)
        {
            sub_grad[3 * i] = grad[3 * non_0_idx[i]];
            sub_grad[3 * i + 1] = grad[3 * non_0_idx[i] + 1];
            sub_grad[3 * i + 2] = grad[3 * non_0_idx[i] + 2];
            for (int j = 0; j < non_0_nb; j++)
            {
                sub_hessian(3 * i, 3 * j) = hessian(3 * non_0_idx[i], 3 * non_0_idx[j]);
                sub_hessian(3 * i, 3 * j + 1) = hessian(3 * non_0_idx[i], 3 * non_0_idx[j] + 1);
                sub_hessian(3 * i, 3 * j + 2) = hessian(3 * non_0_idx[i], 3 * non_0_idx[j] + 2);
                sub_hessian(3 * i + 1, 3 * j) = hessian(3 * non_0_idx[i] + 1, 3 * non_0_idx[j]);
                sub_hessian(3 * i + 1, 3 * j + 1) = hessian(3 * non_0_idx[i] + 1, 3 * non_0_idx[j] + 1);
                sub_hessian(3 * i + 1, 3 * j + 2) = hessian(3 * non_0_idx[i] + 1, 3 * non_0_idx[j] + 2);
                sub_hessian(3 * i + 2, 3 * j) = hessian(3 * non_0_idx[i] + 2, 3 * non_0_idx[j]);
                sub_hessian(3 * i + 2, 3 * j + 1) = hessian(3 * non_0_idx[i] + 2, 3 * non_0_idx[j] + 1);
                sub_hessian(3 * i + 2, 3 * j + 2) = hessian(3 * non_0_idx[i] + 2, 3 * non_0_idx[j] + 2);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> sub_svd(sub_hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::MatrixXd sub_Umat = sub_svd.matrixU();
        Eigen::MatrixXd sub_Vmat = sub_svd.matrixV();
        Eigen::VectorXd sub_singularValues = sub_svd.singularValues();
        Eigen::VectorXd sub_inverse_singularValues;
        sub_inverse_singularValues.resize(sub_singularValues.size());
        sub_inverse_singularValues.setZero();

        // cout << hessian.determinant() << endl;
        for (int i = 0; i < sub_singularValues.size(); i++)
        {
            if (sub_singularValues[i] < 1e-6)
            {
                sub_singularValues[i] = 1e-6;
                sub_inverse_singularValues[i] = 1.0 / sub_singularValues[i];
            }
            else
            {
                sub_inverse_singularValues[i] = 1.0 / sub_singularValues[i];
            }
        }
        auto sub_pseudo_inverse_matrix = sub_Vmat * sub_inverse_singularValues.asDiagonal() * sub_Umat.transpose();
        Eigen::VectorXd sub_direction_v = -sub_pseudo_inverse_matrix * sub_grad;
        /*
        Eigen::MatrixXd sub_hessian1;

        sub_hessian1.resize(non_0_nb * 3, non_0_nb * 3);
        sub_hessian1.setZero();

        Eigen::VectorXd sub_grad1;
        sub_grad1.resize(non_0_nb * 3);
        sub_grad1.setZero();
        for (int i = 0; i < non_0_nb; i++)
        {
            sub_grad1[3 * i] = grad[3 * non_0_idx[i]+grad.size()/2];
            sub_grad1[3 * i + 1] = grad[3 * non_0_idx[i] + 1 + grad.size() / 2];
            sub_grad1[3 * i + 2] = grad[3 * non_0_idx[i] + 2 + grad.size() / 2];
            for (int j = 0; j < non_0_nb; j++)
            {
                sub_hessian1(3 * i, 3 * j) = hessian(3 * non_0_idx[i] + grad.size() / 2, 3 * non_0_idx[j] + grad.size() / 2);
                sub_hessian1(3 * i, 3 * j + 1) = hessian(3 * non_0_idx[i] + grad.size() / 2, 3 * non_0_idx[j] + 1 + grad.size() / 2);
                sub_hessian1(3 * i, 3 * j + 2) = hessian(3 * non_0_idx[i] + grad.size() / 2, 3 * non_0_idx[j] + 2 + grad.size() / 2);
                sub_hessian1(3 * i + 1, 3 * j) = hessian(3 * non_0_idx[i] + 1 + grad.size() / 2, 3 * non_0_idx[j] + grad.size() / 2);
                sub_hessian1(3 * i + 1, 3 * j + 1) = hessian(3 * non_0_idx[i] + 1 + grad.size() / 2, 3 * non_0_idx[j] + 1 + grad.size() / 2);
                sub_hessian1(3 * i + 1, 3 * j + 2) = hessian(3 * non_0_idx[i] + 1 + grad.size() / 2, 3 * non_0_idx[j] + 2 + grad.size() / 2);
                sub_hessian1(3 * i + 2, 3 * j) = hessian(3 * non_0_idx[i] + 2 + grad.size() / 2, 3 * non_0_idx[j] + grad.size() / 2);
                sub_hessian1(3 * i + 2, 3 * j + 1) = hessian(3 * non_0_idx[i] + 2 + grad.size() / 2, 3 * non_0_idx[j] + 1 + grad.size() / 2);
                sub_hessian1(3 * i + 2, 3 * j + 2) = hessian(3 * non_0_idx[i] + 2 + grad.size() / 2, 3 * non_0_idx[j] + 2 + grad.size() / 2);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> sub_svd1(sub_hessian1, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::MatrixXd sub_Umat1 = sub_svd1.matrixU();
        Eigen::MatrixXd sub_Vmat1 = sub_svd1.matrixV();
        Eigen::VectorXd sub_singularValues1 = sub_svd1.singularValues();
        Eigen::VectorXd sub_inverse_singularValues1;
        sub_inverse_singularValues1.resize(sub_singularValues1.size());
        sub_inverse_singularValues1.setZero();

        // cout << hessian.determinant() << endl;
        for (int i = 0; i < sub_singularValues1.size(); i++)
        {
            if (sub_singularValues1[i] < 1e-6)
            {
                sub_singularValues1[i] = 1e-6;
                sub_inverse_singularValues1[i] = 1.0 / sub_singularValues1[i];
            }
            else
            {
                sub_inverse_singularValues1[i] = 1.0 / sub_singularValues1[i];
            }
        }
        auto sub_pseudo_inverse_matrix1 = sub_Vmat1 * sub_inverse_singularValues1.asDiagonal() * sub_Umat1.transpose();
        Eigen::VectorXd sub_direction_v1 = -sub_pseudo_inverse_matrix1 * sub_grad1;
        */

        Eigen::VectorXd direction_v;
        direction_v.resize(grad.size());
        direction_v.setZero();
        for (int i = 0; i < non_0_nb; i++)
        {
            direction_v[3 * non_0_idx[i]] = sub_direction_v[3 * i];
            direction_v[3 * non_0_idx[i] + 1] = sub_direction_v[3 * i + 1];
            direction_v[3 * non_0_idx[i] + 2] = sub_direction_v[3 * i + 2];
            //direction_v[3 * non_0_idx[i]+direction_v.size()/2] = sub_direction_v[3 * i+sub_direction_v.size()/2];
            //direction_v[3 * non_0_idx[i] + 1 + direction_v.size() / 2] = sub_direction_v[3 * i + 1 + sub_direction_v.size() / 2];
            //direction_v[3 * non_0_idx[i] + 2 + direction_v.size() / 2] = sub_direction_v[3 * i + 2 + sub_direction_v.size() / 2];
        }
        /*
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(hessian, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::MatrixXd Umat = svd.matrixU();
        Eigen::MatrixXd Vmat = svd.matrixV();
        Eigen::VectorXd singularValues = svd.singularValues();
        Eigen::VectorXd inverse_singularValues;
        inverse_singularValues.resize(singularValues.size());
        inverse_singularValues.setZero();

       // cout << hessian.determinant() << endl;
        for (int i = 0; i < singularValues.size(); i++)
        {
            if (singularValues[i] < 1e-6)
            {
                singularValues[i] = 1e-6;
                inverse_singularValues[i] = 1.0 / singularValues[i];
            }
            else
            {
                inverse_singularValues[i] = 1.0 / singularValues[i];
            }
        }

        auto pseudo_inverse_matrix = Vmat * inverse_singularValues.asDiagonal() * Umat.transpose();
        Eigen::VectorXd direction_v = -pseudo_inverse_matrix * grad;
        */
        /*
        for (int i = 0; i < grad.size()/3; i++)
        {
            if (abs(grad[3 * i] * grad[3 * i] + grad[3 * i + 1] * grad[3 * i + 1] + grad[3 * i + 2] * grad[3 * i + 2]) < 1e-3)
            {
                direction_v[3 * i] = 0;
                direction_v[3 * i + 1] = 0;
                direction_v[3 * i + 2] = 0;
            }
        }*/
        //direction_v = -grad;
        double max1 = 0;



        direction.resize(spline.rows(), 3);
        for (int i = 0; i < direction.rows(); i++)
        {
            direction(i, 0) = direction_v[3 * i];
            direction(i, 1) = direction_v[3 * i + 1];
            direction(i, 2) = direction_v[3 * i + 2];

            if (sqrt(direction(i, 0) * direction(i, 0) + direction(i, 1) * direction(i, 1) + direction(i, 2) * direction(i, 2)) > max1)
            {
                max1 = sqrt(direction(i, 0) * direction(i, 0) + direction(i, 1) * direction(i, 1) + direction(i, 2) * direction(i, 2));
            }
        }
        //max1 = 1;
        cout << "max=" << max1 << endl;
        //grad = grad / max1;
        //direction_v = direction_v / max1;
        wolfe = direction_v.dot(grad);
        //direction = direction / max1;


        for (int i = 0; i < grad.size() / 3; i++)
        {
            cout << "第" << i << "号点grad: " << grad[3 * i] << " " << grad[3 * i + 1] << " " << grad[3 * i + 2] << endl;
        }

        for (int i = 0; i < direction.rows(); i++)
        {
            cout << "第" << i << "号点direct: " << direction.coeffRef(i, 0) << " " << direction.coeffRef(i, 1) << " " << direction.coeffRef(i, 2) << endl;
        }

        //cout << "第一个端点" << grad[0] << "zuihouyige" << grad[20*3] << endl << endl;
    }

    static void gcl_line_search(Data& spline, const Data& direction,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh, BSplineFitting& ssf, bool exact)
    {
        int flag1 = 1;
        double data_value_ori;
        double data_value = 0;
        double data_value_0 = 0;
        double smooth_value = 0;

        for (int l = 0; l < spline.rows() / 2 && iscollision[zero_num] == 0; l++)
        {
            double tmp1 = direction.coeffRef(l, 0) * direction.coeffRef(l, 0) + direction.coeffRef(l, 1) * direction.coeffRef(l, 1) + direction.coeffRef(l, 2) * direction.coeffRef(l, 2) + direction.coeffRef(l + spline.rows() / 2, 0) * direction.coeffRef(l + spline.rows() / 2, 0) + direction.coeffRef(l + spline.rows() / 2, 1) * direction.coeffRef(l + spline.rows() / 2, 1) + direction.coeffRef(l + spline.rows() / 2, 2) * direction.coeffRef(l + spline.rows() / 2, 2);
            if (tmp1 > 0)
            {
                ssf.change_idx = l;
                bvh.change_idx = l;
                auto dir_temp = direction;
                dir_temp.setZero();
                if (l == spline.rows() / 2 - 1 || 1)
                {
                    dir_temp(l, 0) = direction.coeff(l, 0);
                    dir_temp(l, 1) = direction.coeff(l, 1);
                    dir_temp(l, 2) = direction.coeff(l, 2);
                    dir_temp(l + spline.rows() / 2, 0) = direction.coeff(l + spline.rows() / 2, 0);
                    dir_temp(l + spline.rows() / 2, 1) = direction.coeff(l + spline.rows() / 2, 1);
                    dir_temp(l + spline.rows() / 2, 2) = direction.coeff(l + spline.rows() / 2, 2);
                }
                double step = Step::gcl_position_step(spline, dir_temp, V, F, bvh);
                if (ssf.iter_num == 0)
                {
                    step /= 0.3;
                    step *= 0.5;
                }
                else
                {
                    step = min(pow(1.1, ssf.iter_num) * step, step * 2);
                }
                std::cout.precision(10);
                std::cout << "wolfe:" << wolfe << std::endl;
                // backtracking
                //time0 = clock();
               // double e = Energy::gcl_fast_whole_energy(spline, V, F, bvh, ssf);//fast_
                data_value = ssf.calc_data_term_energy(spline);
                smooth_value = ssf.calc_smooth_term_energyn(spline);
                double e = data_weight * data_value + smooth_weight * smooth_value;
                //cout << "下降方向： " << direction << endl;
                cout << "碰撞检测安全最大步数：" << step << endl;
                double cmp_old = e + 1e-10 * wolfe * step;
                double data_value1 = ssf.calc_data_term_energy(spline + step * dir_temp);
                double smooth_value1 = ssf.calc_smooth_term_energyn(spline + step * dir_temp);
                double e1 = data_weight * data_value1 + smooth_weight * smooth_value1;
                double cmp_new = e1;
                while (cmp_old < cmp_new)//fast_
                {
                    cout << "点编号：" << l << ";step" << step << ";wolfe： " << cmp_old << " " << cmp_new << endl;
                    //cout << "变量： " << spline + step * direction << endl;
                    step *= 0.8;

                    if (step < 1e-3)
                    {
                        step = 0;
                        break;
                    }
                    cmp_old = e + 1e-10 * wolfe * step;
                    data_value1 = ssf.calc_data_term_energy(spline + step * dir_temp);
                    smooth_value1 = ssf.calc_smooth_term_energyn(spline + step * dir_temp);
                    e1 = data_weight * data_value1 + smooth_weight * smooth_value1;
                    cmp_new = e1;

                }

                max_step = step;

                data_value_ori = ssf.calc_data_term_energy(spline);
                data_value = ssf.calc_data_term_energy(spline + step * dir_temp);

                if (flag1 == 1)
                {
                    data_value_0 = data_value_ori;
                }
                data_e = data_weight * data_value;
                double smooth_value_ori = ssf.calc_smooth_term_energyn(spline);
                smooth_value = ssf.calc_smooth_term_energyn(spline + step * dir_temp);

                //stop_condition = abs(data_weight * data_value + smooth_weight * smooth_value - data_weight * data_value_ori - smooth_weight * smooth_value_ori) / (data_weight * data_value_ori + smooth_weight * smooth_value_ori);

                std::cout << "移动最大步长:" << step << std::endl;
                std::cout << "数据能量:" << data_weight * data_value <<
                    " 光滑能量: " << smooth_weight * smooth_value << endl;
                std::cout << "平均能量:" << ssf.avr_energy << endl;
                // " 障碍函数: " << fast_barrier_weight * Energy::gcl_barrier_energy(spline + step * direction, V, F, bvh) << endl;

             //<<" "<<limit_energy(spline+step*direction)
                spline = spline + step * dir_temp;
                flag1 = 0;
            }
        }

        //cout << "线搜索步长：" << step << endl;
        stop_condition = 1 - data_value / data_value_0;
        if (ssf.iter_num == 17)
        {
            ssf.flag123 = 1;
            smooth_value = ssf.calc_smooth_term_energyn(spline);
            ssf.u_smooth_weight = min(1.0 * data_value / (smooth_weight * smooth_value), 1.0);
            //ssf.v_smooth_weight = min(0.5 * data_value / (smooth_weight * smooth_value), 1.0);

        }
        /*
        if (step<1e-8)
        {
            Subdivide::gcl_update_tree(spline, V, F, bvh, udegree, vdegree);
            //bvh.gcl_InitTrajectory(spline, udegree, vdegree);
            cout << " question" << endl;
        }
        */
    }

};


PRJ_END

#endif
