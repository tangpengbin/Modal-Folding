#ifndef CODEGEN_COMPUTEREFLECTIONXPLANETHETA_H
#define CODEGEN_COMPUTEREFLECTIONXPLANETHETA_H

#include <Eigen/Core>

namespace Codegen { 
void computeReflectionXPlaneTheta(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, double& theta);
void computeReflectionXPlaneThetaGradient(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 1>& thetagradient);
void computeReflectionXPlaneThetaHessian(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 9>& thetahessian);

 } 
#endif