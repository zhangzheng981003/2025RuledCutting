#ifndef BVH_H
#define BVH_H

#include "HighOrderCCD/Config/Config.h"
//#include "HighOrderCCD/Element.h"
//#include "HighOrderCCD/Subdivide.h"
//#include "HighOrderCCD/ElementCD.h"
//#include "HighOrderCCD/Distance.h"
#include "src/AABB.h"

PRJ_BEGIN

class BVH
{
public:
    typedef std::vector< std::tuple< int, std::pair<double, double>, Eigen::MatrixXd > > SubdivideTree;
    typedef Eigen::MatrixXd Data;
    typedef std::pair<unsigned int, unsigned int> id_pair;
    int change_idx;
    int flag123 = 0;
    std::vector<int> chosen_subtree_idx;
    aabb::Tree tr_tree;
    aabb::Tree ob_tree;
    aabb::Tree pc_tree;

    BVH();
    ~BVH();

    void InitObstacle(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
    void InitPointcloud(const Eigen::MatrixXd& V);

    void gcl_InitTrajectory(const Data& spline);
    void gcl_ccdInitTrajectory(const Data& spline, const Data& direction);

    void CheckCollision(std::vector<id_pair>& collision_pair, double margin);
    void pcCheckCollision(std::vector<id_pair>& collision_pair, double margin);

};

PRJ_END

#endif