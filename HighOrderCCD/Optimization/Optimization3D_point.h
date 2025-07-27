#ifndef OPTIMIZATION3D_POINT_H
#define OPTIMIZATION3D_POINT_H

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
//#include "gurobi_c++.h"
PRJ_BEGIN

class Optimization3D_point 
{
public:

  typedef Eigen::MatrixXd Data;
 //Eigen::Dynamic
  typedef Eigen::Matrix<double,Eigen::Dynamic,1> inner_derivative_t;//3*(order_num+1)
  typedef Eigen::AutoDiffScalar<inner_derivative_t> inner_scalar_t;
  typedef Eigen::Matrix<inner_scalar_t,Eigen::Dynamic,1> derivative_t;
  typedef Eigen::AutoDiffScalar<derivative_t> scalar_t;
  typedef Eigen::Matrix<scalar_t,Eigen::Dynamic,1> Vec12;
  typedef Eigen::Matrix<scalar_t,1,3> Vec3;

  typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double

  typedef std::vector< std::tuple< int, std::pair<double,double>, Eigen::MatrixXd > > Tree;
  typedef std::tuple< int, std::pair<double,double>, Eigen::MatrixXd >  Node;
};


PRJ_END

#endif
