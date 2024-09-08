#ifndef PERIODIC_DS_HINGE_FIXED_REST_POST_ELEMENT_H
#define PERIODIC_DS_HINGE_FIXED_REST_POST_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif


#include <SimLib/Elements/DofElementBase.h>

class periodicDSHingeFixedRestPoseElement: public DofElementBase
{
public:
	periodicDSHingeFixedRestPoseElement(void);
	~periodicDSHingeFixedRestPoseElement(void);

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const std::vector<Eigen::Vector2d>& verticesTranslationMultipliers, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness);
	std::vector<V3D> computePositions(const Eigen::VectorXd& vx);

	void transformVertices(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);
	Eigen::Matrix<double, 12, 16> transformVerticesJacobian(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);
	std::vector<Eigen::Matrix<double, 12, 16>> transformVerticesHessians(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);

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

	std::vector<Eigen::Vector2d> m_verticesTranslationMultipliers;
	Eigen::Matrix<double, 4 * 3, 4> m_dxdt0t1;


	double m_bendingCurvature;
	Eigen::Vector2d m_bendingDirection;
};

#endif