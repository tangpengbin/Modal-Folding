#include "computeSpringEnergy.h"

namespace Codegen { 
void computeSpringEnergy(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, double& energy){
	double _i_var[11];
	_i_var[0] = (x0(1,0))-(x1(1,0));
	_i_var[1] = (x0(0,0))-(x1(0,0));
	_i_var[2] = (x0(2,0))-(x1(2,0));
	_i_var[3] = (_i_var[0])*(_i_var[0]);
	_i_var[4] = (_i_var[1])*(_i_var[1]);
	_i_var[5] = (_i_var[2])*(_i_var[2]);
	_i_var[6] = (_i_var[4])+(_i_var[3]);
	_i_var[7] = 0.5;
	_i_var[8] = (_i_var[6])+(_i_var[5]);
	_i_var[9] = (_i_var[7])*(stiffness);
	_i_var[10] = (_i_var[9])*(_i_var[8]);
	energy = _i_var[10];
}
void computeSpringEnergyGradient(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, Eigen::Matrix<double, 6, 1>& energygradient){
	double _i_var[16];
	_i_var[0] = 0.5;
	_i_var[1] = (x0(0,0))-(x1(0,0));
	_i_var[2] = (_i_var[0])*(stiffness);
	_i_var[3] = (x0(1,0))-(x1(1,0));
	_i_var[4] = (x0(2,0))-(x1(2,0));
	_i_var[5] = (_i_var[2])*(_i_var[1]);
	_i_var[6] = 2;
	_i_var[7] = (_i_var[2])*(_i_var[3]);
	_i_var[8] = (_i_var[2])*(_i_var[4]);
	_i_var[9] = -1;
	_i_var[10] = (_i_var[6])*(_i_var[5]);
	_i_var[11] = (_i_var[6])*(_i_var[7]);
	_i_var[12] = (_i_var[6])*(_i_var[8]);
	_i_var[13] = (_i_var[10])*(_i_var[9]);
	_i_var[14] = (_i_var[11])*(_i_var[9]);
	_i_var[15] = (_i_var[12])*(_i_var[9]);
	energygradient(0,0) = _i_var[10];
	energygradient(1,0) = _i_var[11];
	energygradient(2,0) = _i_var[12];
	energygradient(3,0) = _i_var[13];
	energygradient(4,0) = _i_var[14];
	energygradient(5,0) = _i_var[15];
}
void computeSpringEnergyHessian(double stiffness, const Eigen::Matrix<double,3,1> & x0, const Eigen::Matrix<double,3,1> & x1, Eigen::Matrix<double, 6, 6>& energyhessian){
	double _i_var[7];
	_i_var[0] = 0.5;
	_i_var[1] = 2;
	_i_var[2] = (_i_var[0])*(stiffness);
	_i_var[3] = (_i_var[2])*(_i_var[1]);
	_i_var[4] = -1;
	_i_var[5] = 0;
	_i_var[6] = (_i_var[4])*(_i_var[3]);
	energyhessian(0,0) = _i_var[3];
	energyhessian(1,0) = _i_var[5];
	energyhessian(2,0) = _i_var[5];
	energyhessian(3,0) = _i_var[6];
	energyhessian(4,0) = _i_var[5];
	energyhessian(5,0) = _i_var[5];
	energyhessian(0,1) = _i_var[5];
	energyhessian(1,1) = _i_var[3];
	energyhessian(2,1) = _i_var[5];
	energyhessian(3,1) = _i_var[5];
	energyhessian(4,1) = _i_var[6];
	energyhessian(5,1) = _i_var[5];
	energyhessian(0,2) = _i_var[5];
	energyhessian(1,2) = _i_var[5];
	energyhessian(2,2) = _i_var[3];
	energyhessian(3,2) = _i_var[5];
	energyhessian(4,2) = _i_var[5];
	energyhessian(5,2) = _i_var[6];
	energyhessian(0,3) = _i_var[6];
	energyhessian(1,3) = _i_var[5];
	energyhessian(2,3) = _i_var[5];
	energyhessian(3,3) = _i_var[3];
	energyhessian(4,3) = _i_var[5];
	energyhessian(5,3) = _i_var[5];
	energyhessian(0,4) = _i_var[5];
	energyhessian(1,4) = _i_var[6];
	energyhessian(2,4) = _i_var[5];
	energyhessian(3,4) = _i_var[5];
	energyhessian(4,4) = _i_var[3];
	energyhessian(5,4) = _i_var[5];
	energyhessian(0,5) = _i_var[5];
	energyhessian(1,5) = _i_var[5];
	energyhessian(2,5) = _i_var[6];
	energyhessian(3,5) = _i_var[5];
	energyhessian(4,5) = _i_var[5];
	energyhessian(5,5) = _i_var[3];
}

 } 

