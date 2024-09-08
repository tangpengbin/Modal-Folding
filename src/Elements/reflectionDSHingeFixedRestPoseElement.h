/*=====================================================================================*/
/*! \file		DSHingeElement.h
\author		bthomasz
\brief		Declaration of DSHingeElement class
*/
/*=====================================================================================*/

#ifndef REFLECTION_DS_HINGE_FIXED_REST_POST_ELEMENT_H
#define REFLECTION_DS_HINGE_FIXED_REST_POST_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif


#include <SimLib/Elements/DofElementBase.h>

class reflectionDSHingeFixedRestPoseElement: public DofElementBase
{
public:
	reflectionDSHingeFixedRestPoseElement(void);
	~reflectionDSHingeFixedRestPoseElement(void);

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness, bool useReflectionXPlane);

	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	//rest state assumed constant (use precomputed variables)
	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf); 
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) override;

	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A);

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) override;
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

	double getCurrentAngle(const Eigen::VectorXd& vx);
	void computedThetadx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& dThetadx);
	void computed2Thetadx2(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& d2Thetadx2);

	int nverts() const { return 4; }
	void setRestAngle(double angle) { m_theta0 = angle; }
	double getRestAngle() { return m_theta0; }
	double getInitialRestAngle() { return m_theta_initial; }
protected:
	SOMaterial const * m_mat;
	//the rest angle
	double m_theta0;
	double m_theta_initial;
	double m_thickness;
	double m_geomFac;

	bool m_useReflectionXPlane;
};

#endif