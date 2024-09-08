#ifndef CODEGEN_COMPUTEDSBENDINGFIXEDRESTPOSECORRECTENERGY_H
#define CODEGEN_COMPUTEDSBENDINGFIXEDRESTPOSECORRECTENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeDSBendingFixedRestPoseCorrectEnergy(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double restTheta, 
	double& energy);
void computeDSBendingFixedRestPoseCorrectEnergyGradient(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double restTheta, 
	Eigen::Matrix<double, 12, 1>& gradient);
void computeDSBendingFixedRestPoseCorrectEnergyHessian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double restTheta, 
	Eigen::Matrix<double, 12, 12>& hessian);
void computeDSBendingFixedRestPoseCorrectEnergyGradientParameterJacobian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & X3, double restTheta, 
	Eigen::Matrix<double, 12, 1>& gradientParameterJacobian);

 } 
#endif