#ifndef DS_HINGE_FIXED_REST_POST_ELEMENT_H
#define DS_HINGE_FIXED_REST_POST_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif
#include <SimLib/Elements/DofElementBase.h>

class DSHingeFixedRestPoseElement: public DofElementBase
{
public:
	DSHingeFixedRestPoseElement(void);
	~DSHingeFixedRestPoseElement(void);

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness, bool hardEdge = false);

	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	//rest state assumed constant (use precomputed variables)
	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf); 
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) override;

	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A);

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) override;
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

	//compute dfdtheta_rest
	void computeGradientParameterJacobian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& gradientParameterJacobian);

	//this use the x to compute theta
	double getCurrentAngle(const Eigen::VectorXd& vx);
	void computedThetadx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& dThetadx);
	void computed2Thetadx2(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& d2Thetadx2);

	int nverts() const { return 4; }
	void setRestAngle(double angle) { m_theta0 = angle; }
	double getRestAngle() { return m_theta0; }
	double getInitialRestAngle() { return m_theta_initial; }
	bool getHardEdge() { return m_hardEdge; }
protected:
	SOMaterial const * m_mat;
	//the rest angle
	double m_theta0;
	double m_theta_initial;
	double m_thickness;
	double m_geomFac;
	bool m_hardEdge;
};

#endif