#ifndef CODEGEN_COMPUTECYLINDERICALTRANSFORMVERTEX_H
#define CODEGEN_COMPUTECYLINDERICALTRANSFORMVERTEX_H

#include <Eigen/Core>

namespace Codegen { 
void computeCylindericalTransformVertex(const Eigen::Matrix<double,3,1> & rotatePoint, const Eigen::Matrix<double,2,1> & translation0, const Eigen::Matrix<double,2,1> & translation1, double bendingCurvature, const Eigen::Matrix<double,2,1> & cylinderDirection, 
	Eigen::Matrix<double, 3, 1>& transformedVertex);
void computeCylindericalTransformVertexJacobian(const Eigen::Matrix<double,3,1> & rotatePoint, const Eigen::Matrix<double,2,1> & translation0, const Eigen::Matrix<double,2,1> & translation1, double bendingCurvature, const Eigen::Matrix<double,2,1> & cylinderDirection, 
	Eigen::Matrix<double, 3, 7>& jacobian);
void computeCylindericalTransformVertexHessians(const Eigen::Matrix<double,3,1> & rotatePoint, const Eigen::Matrix<double,2,1> & translation0, const Eigen::Matrix<double,2,1> & translation1, double bendingCurvature, const Eigen::Matrix<double,2,1> & cylinderDirection, 
	Eigen::Matrix<double, 3, 7>& hessian0, Eigen::Matrix<double, 3, 7>& hessian1, Eigen::Matrix<double, 3, 7>& hessian2, Eigen::Matrix<double, 3, 7>& hessian3, Eigen::Matrix<double, 3, 7>& hessian4, 
	Eigen::Matrix<double, 3, 7>& hessian5, Eigen::Matrix<double, 3, 7>& hessian6);

 } 
#endif