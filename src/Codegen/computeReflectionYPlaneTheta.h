#ifndef CODEGEN_COMPUTEREFLECTIONYPLANETHETA_H
#define CODEGEN_COMPUTEREFLECTIONYPLANETHETA_H

#include <Eigen/Core>

namespace Codegen { 
void computeReflectionYPlaneTheta(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double& theta);
void computeReflectionYPlaneThetaGradient(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 1>& thetagradient);
void computeReflectionYPlaneThetaHessian(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 9>& thetahessian);

 } 
#endif