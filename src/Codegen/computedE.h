#ifndef CODEGEN_COMPUTEDE_H
#define CODEGEN_COMPUTEDE_H

#include <Eigen/Core>

namespace Codegen { 
void computedE(const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & x4, const Eigen::Matrix<double,3,1> & x1Undef, 
	const Eigen::Matrix<double,3,1> & x2Undef, const Eigen::Matrix<double,3,1> & x3Undef, const Eigen::Matrix<double,3,1> & x4Undef, Eigen::Matrix<double, 3, 3>& E);
void computedEjacobians(const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & x4, const Eigen::Matrix<double,3,1> & x1Undef, 
	const Eigen::Matrix<double,3,1> & x2Undef, const Eigen::Matrix<double,3,1> & x3Undef, const Eigen::Matrix<double,3,1> & x4Undef, Eigen::Matrix<double, 3, 3>& jacobian0, Eigen::Matrix<double, 3, 3>& jacobian1, 
	Eigen::Matrix<double, 3, 3>& jacobian2, Eigen::Matrix<double, 3, 3>& jacobian3, Eigen::Matrix<double, 3, 3>& jacobian4, Eigen::Matrix<double, 3, 3>& jacobian5, Eigen::Matrix<double, 3, 3>& jacobian6, 
	Eigen::Matrix<double, 3, 3>& jacobian7, Eigen::Matrix<double, 3, 3>& jacobian8, Eigen::Matrix<double, 3, 3>& jacobian9, Eigen::Matrix<double, 3, 3>& jacobian10, Eigen::Matrix<double, 3, 3>& jacobian11);

 } 
#endif