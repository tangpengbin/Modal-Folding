#ifndef CODEGEN_COMPUTESPRINGENERGY_H
#define CODEGEN_COMPUTESPRINGENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeSpringEnergy(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double& energy);
void computeSpringEnergyGradient(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, Eigen::Matrix<double, 6, 1>& energygradient);
void computeSpringEnergyHessian(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, Eigen::Matrix<double, 6, 6>& energyhessian);

 } 
#endif