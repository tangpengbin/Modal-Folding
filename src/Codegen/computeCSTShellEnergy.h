#ifndef CODEGEN_COMPUTECSTSHELLENERGY_H
#define CODEGEN_COMPUTECSTSHELLENERGY_H

#include <Eigen/Core>

namespace Codegen { 
void computeCSTShellEnergy(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, double& energy);
void computeCSTShellEnergyGradient(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 1>& gradient);
void computeCSTShellEnergyHessian(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& hessian);
void computeCSTShellEnergyHessianProductJacobian(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, const Eigen::Matrix<double,9,1> & v, 
	Eigen::Matrix<double, 9, 9>& hessianProductJacobian);
void computeCSTShellEnergyRestConfJac(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& rest_conf_jac);
void computeCSTShellEnergyParam_jac(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 2>& param_jac);
void computeCSTShellEnergyParam_thirdOrderDerivatives(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& d3_0, 
	Eigen::Matrix<double, 9, 9>& d3_1);
void computeCSTShellEnergyThirdOrderDerivatives(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& d3_0, 
	Eigen::Matrix<double, 9, 9>& d3_1, Eigen::Matrix<double, 9, 9>& d3_2, Eigen::Matrix<double, 9, 9>& d3_3, Eigen::Matrix<double, 9, 9>& d3_4, Eigen::Matrix<double, 9, 9>& d3_5, 
	Eigen::Matrix<double, 9, 9>& d3_6, Eigen::Matrix<double, 9, 9>& d3_7, Eigen::Matrix<double, 9, 9>& d3_8);
void computeCSTShellEnergyRestConfThirdOrderDerivatives(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 9, 9>& d3_0, 
	Eigen::Matrix<double, 9, 9>& d3_1, Eigen::Matrix<double, 9, 9>& d3_2, Eigen::Matrix<double, 9, 9>& d3_3, Eigen::Matrix<double, 9, 9>& d3_4, Eigen::Matrix<double, 9, 9>& d3_5, 
	Eigen::Matrix<double, 9, 9>& d3_6, Eigen::Matrix<double, 9, 9>& d3_7, Eigen::Matrix<double, 9, 9>& d3_8);
void computeCSTShellEnergyFourthOrderDerivatives(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, const Eigen::Matrix<double,9,1> & v1, 
	const Eigen::Matrix<double,9,1> & v2, Eigen::Matrix<double, 9, 9>& fourthOrderDerivative);
void computeCSTShellEnergyThirdOrderToParameterFourthOrderDerivatives(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, const Eigen::Matrix<double,9,1> & v1, 
	const Eigen::Matrix<double,9,1> & v2, Eigen::Matrix<double, 9, 9>& thirdOrderToParameterFourthOrderDerivative);

 } 
#endif