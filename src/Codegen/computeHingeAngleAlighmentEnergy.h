#ifndef CODEGEN_COMPUTEHINGEANGLEALIGHMENTENERGY_H
#define CODEGEN_COMPUTEHINGEANGLEALIGHMENTENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeHingeAngleAlighmentEnergy(double theta, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & fit_X0, const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, double& energy);
void computeHingeAngleAlighmentEnergyGradient(double theta, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & fit_X0, const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, Eigen::Matrix<double, 1, 1>& energygradient);
void computeHingeAngleAlighmentEnergyHessian(double theta, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x3, 
	const Eigen::Matrix<double,3,1> & fit_X0, const Eigen::Matrix<double,3,1> & fit_X1, const Eigen::Matrix<double,3,1> & fit_X2, const Eigen::Matrix<double,3,1> & fit_X3, Eigen::Matrix<double, 1, 1>& energyhessian);

 } 
#endif