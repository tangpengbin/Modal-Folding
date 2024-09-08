#ifndef CODEGEN_COMPUTEHINGEEDGEALIGHMENTENERGY_H
#define CODEGEN_COMPUTEHINGEEDGEALIGHMENTENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeHingeEdgeAlighmentEnergy(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & fit_X0, 
	const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, double& energy);
void computeHingeEdgeAlighmentEnergyGradient(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & fit_X0, 
	const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, Eigen::Matrix<double, 12, 1>& energygradient);
void computeHingeEdgeAlighmentEnergyHessian(const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, const Eigen::Matrix<double,3,1> & fit_X0, 
	const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, Eigen::Matrix<double, 12, 12>& energyhessian);

 } 
#endif