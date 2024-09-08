#ifndef CODEGEN_COMPUTEREFLECTIONYPLANEDSBENDINGFIXEDRESTPOSEENERGY_H
#define CODEGEN_COMPUTEREFLECTIONYPLANEDSBENDINGFIXEDRESTPOSEENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeReflectionYPlaneDSBendingFixedRestPoseEnergy(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, double& energy);
void computeReflectionYPlaneDSBendingFixedRestPoseEnergyGradient(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, Eigen::Matrix<double, 9, 1>& energygradient);
void computeReflectionYPlaneDSBendingFixedRestPoseEnergyHessian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, Eigen::Matrix<double, 9, 9>& energyhessian);

 } 
#endif