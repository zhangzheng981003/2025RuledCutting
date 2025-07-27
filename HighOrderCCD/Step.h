#ifndef STEP_H
#define STEP_H

#include "Config/Config.h"
#include "HighOrderCCD/BVH/BVH.h"
#include "HighOrderCCD/CCD/CCD.h"
#include <vector>

PRJ_BEGIN

class Step
{
public:
    typedef Eigen::MatrixXd Data;

    static double gcl_position_step(const Data& spline, const Data& direction,
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F,
        BVH& bvh)
    {
        cout << subdivide_tree.size() << endl;
        clock_t time1 = clock();
        double step = 1;
        vector<int> colnum(subdivide_tree.size() * 128, 0);
        std::vector<std::pair<unsigned int, unsigned int>> collision_pair, temp_pair;
        bvh.gcl_InitTrajectory(spline);
        bvh.CheckCollision(collision_pair, offset);
        int num_collision_pair = collision_pair.size();
        if (iter_ == 0 && num_collision_pair > 0)
        {
            step = 0;
            iscollision[zero_num] = 1;
        }

        cout << "coll_size0000:====" << collision_pair.size() << endl;
        ofstream ioszz6("coverarea2.txt");
        vector<int> color1(16446, 0);

        for (int i = 0; i < collision_pair.size(); i++)
        {
            int tr_id = collision_pair[i].first;
            int ob_id = collision_pair[i].second;
            //cout << "树的编号为:" << tr_id << ";" ;
            //cout << "面的编号为：" << ob_id << ";" << endl << endl;
            //int sp_id = std::get<0>(subdivide_tree[tr_id]);
            if (color1[ob_id] == 0)
            {
                ioszz6 << ob_id << endl;
                color1[ob_id] = 1;
            }
            colnum[tr_id]++;
        }
        //clock_t time0 = clock();
        bvh.gcl_ccdInitTrajectory(spline, step * direction);
        bvh.CheckCollision(collision_pair, offset);
        cout << "coll_size1111:====" << collision_pair.size() << endl;
        //std::cout << "bvh_init: " << collision_pair.size() << std::endl;
        double init_step = step;
        if ((int)collision_pair.size() > std::max(1000, 2 * num_collision_pair))
        {
            while (true)
            {

                init_step *= 0.5;
                temp_pair = collision_pair;
                bvh.gcl_ccdInitTrajectory(spline, init_step * direction);
                bvh.CheckCollision(collision_pair, offset);
                std::cout << std::endl << "bvhsize:" << collision_pair.size() << std::endl;
                if ((int)collision_pair.size() <= std::max(1000, 2 * num_collision_pair) && (int)collision_pair.size() > num_collision_pair)
                {
                    break;
                }
                else if ((int)collision_pair.size() <= std::max(1000, 2 * num_collision_pair) && (int)collision_pair.size() <= num_collision_pair)
                {
                    init_step *= 5;
                    collision_pair = temp_pair;
                    break;
                }
                cout << "zz" << endl;
            }
        }
        step = init_step;

        //clock_t time1 = clock();
        //std::cout<<std::endl<<"bvhtime:"<<(time1-time0)/(CLOCKS_PER_SEC/1000)<<std::endl;
        bool is_collided = true;
        //std::cout << "bvh: " << collision_pair.size() << std::endl;

        double temp_step = step;
        double step0 = 0.0;
        //std::cout<<"start:"<<collision_pair.size()<<std::endl;
        //time0 = clock();

        for (int iter0 = 0; iter0 < 1; iter0++)//
        {
            if (collision_pair.size() > 0)
            {
                is_collided = true;
            }
            else
            {
                is_collided = false;
            }

            ///std::vector<std::pair<unsigned int, unsigned int>> temp_pair = collision_pair;

            temp_step = step - step0;

            for (unsigned int it = 0; it < collision_pair.size(); it++)
            {
                int tr_id = collision_pair[it].first;
                int ob_id = collision_pair[it].second;
                //cout << "树的编号为:" << tr_id << ";" ;
                //cout << "面的编号为：" << ob_id << ";" << endl << endl;
                //int sp_id = std::get<0>(subdivide_tree[tr_id]);

                Eigen::MatrixXd basis = std::get<2>(subdivide_tree[bvh.chosen_subtree_idx[tr_id / 128]]);
                Eigen::MatrixXd bz, bz_d;
                bz = spline;
                bz_d = direction;
                auto bz1 = bz;
                int l = tr_id % 128;
                for (int k = 0; k < bz.rows() / 2; k++)
                {
                    bz1(k, 0) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 0) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 0);
                    bz1(k, 1) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 1) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 1);
                    bz1(k, 2) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 2) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 2);
                    bz1(k + bz.rows() / 2, 0) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 0) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 0);
                    bz1(k + bz.rows() / 2, 1) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 1) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 1);
                    bz1(k + bz.rows() / 2, 2) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 2) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 2);
                }

                auto bz2 = bz_d;

                for (int k = 0; k < bz_d.rows() / 2; k++)
                {
                    bz2(k, 0) = (128.0 - (double)l) / 128.0 * bz_d.coeff(k, 0) + (double)l / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 0);
                    bz2(k, 1) = (128.0 - (double)l) / 128.0 * bz_d.coeff(k, 1) + (double)l / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 1);
                    bz2(k, 2) = (128.0 - (double)l) / 128.0 * bz_d.coeff(k, 2) + (double)l / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 2);
                    bz2(k + bz_d.rows() / 2, 0) = (127.0 - (double)l) / 128.0 * bz_d.coeff(k, 0) + (1 + (double)l) / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 0);
                    bz2(k + bz_d.rows() / 2, 1) = (127.0 - (double)l) / 128.0 * bz_d.coeff(k, 1) + (1 + (double)l) / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 1);
                    bz2(k + bz_d.rows() / 2, 2) = (127.0 - (double)l) / 128.0 * bz_d.coeff(k, 2) + (1 + (double)l) / 128.0 * bz_d.coeff(k + bz_d.rows() / 2, 2);
                }
                int sub_ctr_num = basis.rows();
                Eigen::MatrixXd P(sub_ctr_num, 3), D(sub_ctr_num, 3);
                P = basis * bz1;
                D = basis * bz2;

                int f0 = F(ob_id, 0); int f1 = F(ob_id, 1); int f2 = F(ob_id, 2);

                Eigen::Matrix3d _position;
                _position << V.row(f0), V.row(f1), V.row(f2);

                is_collided = CCD::gcl_KDOPCCD(P, D, _position, offset, step0, step0 + temp_step);
                //std::cout<<"is_collided:"<<is_collided<<" \n";

                //is_collided=CCD::ExactCCD(P,D,_position);
                while (is_collided)
                {
                    //is_collided= CCD::KDOPCCD(P,D,_position,offset,step0,step0+temp_step); 
                    is_collided = CCD::gcl_KDOPCCD(P, D, _position, offset, step0, step0 + temp_step);  //cgal
                    if (is_collided)
                    {
                        temp_step *= 0.8;
                        // cout << temp_step << endl;
                        clock_t time2 = clock();
                        double consume_time = (time2 - time1) / CLOCKS_PER_SEC;

                        if (consume_time > 5)
                        {
                            temp_step = 0;
                            goto flag;
                        }
                    }
                }
            }

            if (temp_step / (step0 + temp_step) < 1e-8)
            {
                step0 = step0 + temp_step;
                temp_step = step0;
                std::cout << "temp_step:" << step0 << std::endl;
                break;
            }

            step0 = step0 + temp_step;
            temp_step = step0;
            std::cout << "temp_step:" << step0 << std::endl;
        }


    flag:
        step = temp_step * 0.3;
        for (int i = 0; i < 1024 * 8; i++)
        {
            if (colnum[i] > 0)
            {
                cout << "第" << i << "个树" << " 碰撞面的个数为" << colnum[i] << endl;
            }
        }
    flag1: return step;
    }

};

PRJ_END

#endif