#ifndef CODEGEN_COMPUTESHAPEOPERATOR_H
#define CODEGEN_COMPUTESHAPEOPERATOR_H

#include <Eigen/Core>

namespace Codegen { 
void computeShapeOperator(const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,3,1> & v2, const Eigen::Matrix<double,3,1> & v3, const Eigen::Matrix<double,3,1> & v4, const Eigen::Matrix<double,3,1> & v5, 
	const Eigen::Matrix<double,3,1> & v6, const Eigen::Matrix<double,3,1> & undeformed_v1, const Eigen::Matrix<double,3,1> & undeformed_v2, const Eigen::Matrix<double,3,1> & undeformed_v3, const Eigen::Matrix<double,3,1> & binary_multiplier, 
	Eigen::Matrix<double, 2, 2>& shapeOperator);
void computeShapeOperatorjacobians(const Eigen::Matrix<double,3,1> & v1, const Eigen::Matrix<double,3,1> & v2, const Eigen::Matrix<double,3,1> & v3, const Eigen::Matrix<double,3,1> & v4, const Eigen::Matrix<double,3,1> & v5, 
	const Eigen::Matrix<double,3,1> & v6, const Eigen::Matrix<double,3,1> & undeformed_v1, const Eigen::Matrix<double,3,1> & undeformed_v2, const Eigen::Matrix<double,3,1> & undeformed_v3, const Eigen::Matrix<double,3,1> & binary_multiplier, 
	Eigen::Matrix<double, 2, 2>& jacobian0, Eigen::Matrix<double, 2, 2>& jacobian1, Eigen::Matrix<double, 2, 2>& jacobian2, Eigen::Matrix<double, 2, 2>& jacobian3, Eigen::Matrix<double, 2, 2>& jacobian4, 
	Eigen::Matrix<double, 2, 2>& jacobian5, Eigen::Matrix<double, 2, 2>& jacobian6, Eigen::Matrix<double, 2, 2>& jacobian7, Eigen::Matrix<double, 2, 2>& jacobian8, Eigen::Matrix<double, 2, 2>& jacobian9, 
	Eigen::Matrix<double, 2, 2>& jacobian10, Eigen::Matrix<double, 2, 2>& jacobian11, Eigen::Matrix<double, 2, 2>& jacobian12, Eigen::Matrix<double, 2, 2>& jacobian13, Eigen::Matrix<double, 2, 2>& jacobian14, 
	Eigen::Matrix<double, 2, 2>& jacobian15, Eigen::Matrix<double, 2, 2>& jacobian16, Eigen::Matrix<double, 2, 2>& jacobian17);

 } 
#endif