#ifndef CODEGEN_COMPUTECSTSHELLENERGYGREENSTRAIN_H
#define CODEGEN_COMPUTECSTSHELLENERGYGREENSTRAIN_H

#include <Eigen/Core>

namespace Codegen { 
void computeCSTShellEnergyGreenStrain(double youngsModulus, double poissonsRatio, double thickness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, 
	const Eigen::Matrix<double,3,1> & x2, const Eigen::Matrix<double,3,1> & x0Undef, const Eigen::Matrix<double,3,1> & x1Undef, const Eigen::Matrix<double,3,1> & x2Undef, Eigen::Matrix<double, 2, 2>& greenStrain);

 } 
#endif