#ifndef PERIODIC_CST_SHELL_ELEMENT_H
#define PERIODIC_CST_SHELL_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif
#include "DofElementBase.h"
class periodicCSTShellElement :
	public DofElementBase
{
public:
	periodicCSTShellElement();
	~periodicCSTShellElement();

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const std::vector<Eigen::Vector2d>& verticesTranslationMultipliers, const Eigen::VectorXd& vx, SOMaterial* pMat);
	void convertEandNuToLame(double E, double nu, double& lambda, double& mu);
	double computeArea(const Eigen::VectorXd& vx);

	void transformVertices(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);

	Eigen::Matrix<double, 9, 13> transformVerticesJacobian(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);

	std::vector<Eigen::Matrix<double, 9, 13>> transformVerticesHessians(std::vector<V3D>& x, std::vector<V3D>& X,
		const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
		const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest);

	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	void computeMaxStrain(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, double& maxStrain);

	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf);
	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A);
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) override;

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) override;
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

	double getThickness() { return m_h; }
	void setThickness(double thinckness) { m_h = thinckness; }

protected:
	double m_h;

	std::vector<Eigen::Vector2d> m_verticesTranslationMultipliers;
	Eigen::Matrix<double, 3 * 3, 4> m_dxdt0t1;

	double m_bendingCurvature;
	Eigen::Vector2d m_bendingDirection;
};

#endif