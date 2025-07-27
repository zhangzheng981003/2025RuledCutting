
#include "BVH.h"

PRJ_BEGIN

BVH::BVH()
{

}
BVH::~BVH()
{

}
void BVH::InitObstacle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    unsigned int n_ob = F.rows();
    int dim = aabb_axis.size();
    /*
    for(int k=0;k<dim;k++)
    {
      aabb_axis[k].normalize();
    }
    */
    ob_tree = aabb::Tree(dim, 0.0, n_ob, true);

    std::cout << "\nInserting ob particles into AABB tree ...\n";
    for (unsigned int i = 0; i < n_ob; i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for (int k = 0; k < dim; k++)
        {
            upperBound[k] = -INFINITY;
            lowerBound[k] = INFINITY;
        }
        for (int k = 0; k < dim; k++)
        {
            for (int j = 0; j < 3; j++)
            {
                double level = aabb_axis[k].dot(V.row(F(i, j)));
                if (level < lowerBound[k])
                    lowerBound[k] = level;
                if (level > upperBound[k])
                    upperBound[k] = level;
                if (abs(upperBound[k] - lowerBound[k]) > 0.05)
                {
                    std::cout << "第" << i << "个面:" << lowerBound[k] << " " << upperBound[k] << " " << upperBound[k] - lowerBound[k] << std::endl;
                }
            }

        }


        ob_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::InitPointcloud(const Eigen::MatrixXd& V)
{
    unsigned int n_pc = V.rows();
    int dim = aabb_axis.size();
    /*
    for(int k=0;k<dim;k++)
    {
      aabb_axis[k].normalize();
    }
    */
    pc_tree = aabb::Tree(dim, 0.0, n_pc, true);

    std::cout << "\nInserting pc particles into AABB tree ...\n";
    for (unsigned int i = 0; i < n_pc; i++)
    {
        std::vector<double> lowerBound(dim);
        std::vector<double> upperBound(dim);

        for (int k = 0; k < dim; k++)
        {
            upperBound[k] = -INFINITY;
            lowerBound[k] = INFINITY;
        }


        for (int k = 0; k < dim; k++)
        {
            double level = aabb_axis[k].dot(V.row(i));
            if (level < lowerBound[k])
                lowerBound[k] = level;
            if (level > upperBound[k])
                upperBound[k] = level;
        }
        pc_tree.insertParticle(i, lowerBound, upperBound);
    }

}

void BVH::gcl_InitTrajectory(const Data& spline)
{
    chosen_subtree_idx.clear();

    for (int it = 0; it < subdivide_tree.size(); it++)
    {
        Eigen::MatrixXd basis = std::get<2>(subdivide_tree[it]);
        if (basis.col(change_idx).norm() > 0)
        {
            chosen_subtree_idx.push_back(it);
        }
    }


    unsigned int n_tr = chosen_subtree_idx.size() * 128;
    std::cout << "树的大小为：" << n_tr << std::endl;
    int dim = aabb_axis.size();
    tr_tree = aabb::Tree(dim, 0.0, n_tr, true);

    int fitting_ctr_num = spline.rows();
    Eigen::MatrixXd bz;
    bz = spline;

    //std::cout << n_tr<<"\nInserting tr particles into AABB tree ...\n";
    for (unsigned int it = 0; it < chosen_subtree_idx.size(); it++)
    {
        unsigned int i = chosen_subtree_idx[it];
        for (int l = 0; l < 128; l++)
        {

            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[i]);
            auto bz1 = bz;
            for (int k = 0; k < bz.rows() / 2; k++)
            {
                bz1(k, 0) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 0) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 0);
                bz1(k, 1) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 1) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 1);
                bz1(k, 2) = (128.0 - (double)l) / 128.0 * bz.coeff(k, 2) + (double)l / 128.0 * bz.coeff(k + bz.rows() / 2, 2);
                bz1(k + bz.rows() / 2, 0) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 0) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 0);
                bz1(k + bz.rows() / 2, 1) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 1) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 1);
                bz1(k + bz.rows() / 2, 2) = (127.0 - (double)l) / 128.0 * bz.coeff(k, 2) + (1 + (double)l) / 128.0 * bz.coeff(k + bz.rows() / 2, 2);
            }
            std::vector<double> lowerBound(dim);
            std::vector<double> upperBound(dim);

            //int sp_id = std::get<0>(subdivide_tree[i]);

            int sub_ctr_num = basis.rows();
            //std::cout << basis<< std::endl;
            std::vector<Eigen::RowVector3d> P(sub_ctr_num);
            for (int k = 0; k < dim; k++)
            {
                upperBound[k] = -INFINITY;
                lowerBound[k] = INFINITY;
            }

            for (int j = 0; j < sub_ctr_num; j++)
            {
                P[j].setZero();
                for (int j0 = 0; j0 < fitting_ctr_num; j0++)
                {
                    P[j] += basis(j, j0) * bz1.row(j0);
                }
                //std::cout << P[j][0] << " " << P[j][1] << " " << P[j][2] << std::endl;
                for (int k = 0; k < dim; k++)
                {
                    double level = aabb_axis[k].dot(P[j]);
                    if (level < lowerBound[k])
                        lowerBound[k] = level;
                    if (level > upperBound[k])
                        upperBound[k] = level;
                }
            }
            //std::cout << std::endl;
            tr_tree.insertParticle(128 * it + l, lowerBound, upperBound);

        }
    }

}

void BVH::gcl_ccdInitTrajectory(const Data& spline, const Data& direction)
{
    chosen_subtree_idx.clear();
    for (int it = 0; it < subdivide_tree.size(); it++)
    {
        Eigen::MatrixXd basis = std::get<2>(subdivide_tree[it]);
        if (basis.col(change_idx).norm() != 0)
        {
            chosen_subtree_idx.push_back(it);
        }
    }


    unsigned int n_tr = chosen_subtree_idx.size() * 128;
    int dim = aabb_axis.size();
    tr_tree = aabb::Tree(dim, 0.0, n_tr, true);

    int fitting_ctr_num = spline.rows();

    Eigen::MatrixXd bz, bz_d;
    bz = spline;
    bz_d = direction;

    //std::cout << "\nUpdate tr particles in AABB tree ...\n";
    for (unsigned int it = 0; it < chosen_subtree_idx.size(); it++)
    {
        unsigned int i = chosen_subtree_idx[it];
        for (int l = 0; l < 128; l++)
        {
            auto bz1 = bz;

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

            std::vector<double> lowerBound(dim);
            std::vector<double> upperBound(dim);

            //int sp_id = std::get<0>(subdivide_tree[i]);
            Eigen::MatrixXd basis = std::get<2>(subdivide_tree[i]);
            int sub_ctr_num = basis.rows();
            std::vector<Eigen::RowVector3d> P(sub_ctr_num), D(sub_ctr_num);

            for (int k = 0; k < dim; k++)
            {
                upperBound[k] = -INFINITY;
                lowerBound[k] = INFINITY;
            }
            for (int j = 0; j < sub_ctr_num; j++)
            {
                P[j].setZero();
                D[j].setZero();

                for (int j0 = 0; j0 < fitting_ctr_num; j0++)
                {
                    P[j] += basis(j, j0) * bz1.row(j0);
                    D[j] += basis(j, j0) * bz2.row(j0);
                }

                for (int k = 0; k < dim; k++)
                {

                    double level = aabb_axis[k].dot(P[j]);
                    if (level < lowerBound[k])
                        lowerBound[k] = level;
                    if (level > upperBound[k])
                        upperBound[k] = level;

                    level = aabb_axis[k].dot(P[j] + D[j]);
                    if (level < lowerBound[k])
                        lowerBound[k] = level;
                    if (level > upperBound[k])
                        upperBound[k] = level;
                }
            }

            tr_tree.insertParticle(128 * it + l, lowerBound, upperBound);
        }
    }

}


// 返回环境小于margin的面的索引值和分片样条的当中的第几条样条的索引值。
void BVH::CheckCollision(std::vector<id_pair>& collision_pair, double margin)
{
    collision_pair = tr_tree.query(ob_tree, margin);
}

void BVH::pcCheckCollision(std::vector<id_pair>& collision_pair, double margin)
{
    collision_pair = tr_tree.query(pc_tree, margin);
}





PRJ_END

