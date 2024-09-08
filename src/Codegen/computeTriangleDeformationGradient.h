#ifndef CODEGEN_COMPUTETRIANGLEDEFORMATIONGRADIENT_H
#define CODEGEN_COMPUTETRIANGLEDEFORMATIONGRADIENT_H

#include <Eigen/Core>

namespace Codegen { 
void computeTriangleDeformationGradient(const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 1>& F_rowwise);
void computeTriangleDeformationGradientJacobian(const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 9>& jacobian);
void computeTriangleDeformationGradientHessians(const Eigen::Matrix<double,3,1> & X0, const Eigen::Matrix<double,3,1> & X1, const Eigen::Matrix<double,3,1> & X2, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, Eigen::Matrix<double, 9, 9>& hessian0, Eigen::Matrix<double, 9, 9>& hessian1, Eigen::Matrix<double, 9, 9>& hessian2, Eigen::Matrix<double, 9, 9>& hessian3, 
	Eigen::Matrix<double, 9, 9>& hessian4, Eigen::Matrix<double, 9, 9>& hessian5, Eigen::Matrix<double, 9, 9>& hessian6, Eigen::Matrix<double, 9, 9>& hessian7, Eigen::Matrix<double, 9, 9>& hessian8);

 } 
#endif