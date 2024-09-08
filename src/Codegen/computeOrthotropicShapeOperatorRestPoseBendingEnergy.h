#ifndef CODEGEN_COMPUTEORTHOTROPICSHAPEOPERATORRESTPOSEBENDINGENERGY_H
#define CODEGEN_COMPUTEORTHOTROPICSHAPEOPERATORRESTPOSEBENDINGENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeOrthotropicShapeOperatorRestPoseBendingEnergy(double Y0, double Y1, double Y01, double G01, double thickness, 
	const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,3,1> & v2, const Eigen::Matrix<double,3,1> & v3, const Eigen::Matrix<double,3,1> & v4, const Eigen::Matrix<double,3,1> & v5, 
	const Eigen::Matrix<double,3,1> & v6, const Eigen::Matrix<double,3,1> & undeformed_v1, const Eigen::Matrix<double,3,1> & undeformed_v2, const Eigen::Matrix<double,3,1> & undeformed_v3, const Eigen::Matrix<double,2,2> & target_shapeOperator, 
	const Eigen::Matrix<double,3,1> & binary_multiplier, double& energy);
void computeOrthotropicShapeOperatorRestPoseBendingEnergyGradient(double Y0, double Y1, double Y01, double G01, double thickness, 
	const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,3,1> & v2, const Eigen::Matrix<double,3,1> & v3, const Eigen::Matrix<double,3,1> & v4, const Eigen::Matrix<double,3,1> & v5, 
	const Eigen::Matrix<double,3,1> & v6, const Eigen::Matrix<double,3,1> & undeformed_v1, const Eigen::Matrix<double,3,1> & undeformed_v2, const Eigen::Matrix<double,3,1> & undeformed_v3, const Eigen::Matrix<double,2,2> & target_shapeOperator, 
	const Eigen::Matrix<double,3,1> & binary_multiplier, Eigen::Matrix<double, 18, 1>& energygradient);
void computeOrthotropicShapeOperatorRestPoseBendingEnergyHessian(double Y0, double Y1, double Y01, double G01, double thickness, 
	const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,3,1> & v2, const Eigen::Matrix<double,3,1> & v3, const Eigen::Matrix<double,3,1> & v4, const Eigen::Matrix<double,3,1> & v5, 
	const Eigen::Matrix<double,3,1> & v6, const Eigen::Matrix<double,3,1> & undeformed_v1, const Eigen::Matrix<double,3,1> & undeformed_v2, const Eigen::Matrix<double,3,1> & undeformed_v3, const Eigen::Matrix<double,2,2> & target_shapeOperator, 
	const Eigen::Matrix<double,3,1> & binary_multiplier, Eigen::Matrix<double, 18, 18>& energyhessian);

 } 
#endif