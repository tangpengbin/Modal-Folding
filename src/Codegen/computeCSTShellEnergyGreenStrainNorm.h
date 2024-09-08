#ifndef CODEGEN_COMPUTECSTSHELLENERGYGREENSTRAINNORM_H
#define CODEGEN_COMPUTECSTSHELLENERGYGREENSTRAINNORM_H

#include <Eigen/Core>

namespace Codegen { 
void computeCSTShellEnergyGreenStrainNorm(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, double& greenStrainNorm);
void computeCSTShellEnergyGreenStrainNormGradient(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 1>& greenStrainNormgradient);
void computeCSTShellEnergyGreenStrainNormHessian(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& greenStrainNormhessian);

 } 
#endif