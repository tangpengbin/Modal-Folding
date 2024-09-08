#ifndef CODEGEN_COMPUTEDTHETA_H
#define CODEGEN_COMPUTEDTHETA_H

#include <Eigen/Core>

namespace Codegen { 
void computedTheta(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, double& theta);
void computedThetaGradient(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, Eigen::Matrix<double, 12, 1>& thetagradient);
void computedThetaHessian(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, Eigen::Matrix<double, 12, 12>& thetahessian);

 } 
#endif