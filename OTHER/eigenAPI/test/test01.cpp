#include <iostream>
#include "../EigenAPI.h"
 
using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
// using Eigen::Matrix5d;
// using Eigen::Matrix6d
 
int main()
{
  // MatrixXd m(2,2);
  Matrix2d m(2,2);
  m(0,0) = 3;
  m(1,0) = 2.5;
  m(0,1) = -1;
  m(1,1) = m(1,0) + m(0,1);
  std::cout << m << std::endl;
}