#ifndef CODEGEN_COMPUTEREFLECTIONXPLANEDSBENDINGFIXEDRESTPOSEENERGY_H
#define CODEGEN_COMPUTEREFLECTIONXPLANEDSBENDINGFIXEDRESTPOSEENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeReflectionXPlaneDSBendingFixedRestPoseEnergy(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, double& energy);
void computeReflectionXPlaneDSBendingFixedRestPoseEnergyGradient(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, Eigen::Matrix<double, 9, 1>& energygradient);
void computeReflectionXPlaneDSBendingFixedRestPoseEnergyHessian(double bendingStiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & X0, 
	const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, double restTheta, Eigen::Matrix<double, 9, 9>& energyhessian);

 } 
#endif