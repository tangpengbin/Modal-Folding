#ifndef CST_SHELL_ELEMENT_H
#define CST_SHELL_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif
#include "DofElementBase.h"
class CSTShellElement :
	public DofElementBase
{
public:
	CSTShellElement();
	~CSTShellElement();

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat);
	void convertEandNuToLame(double E, double nu, double& lambda, double& mu);
	double computeArea(const Eigen::VectorXd& vx);

	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	void computeMaxStrain(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, double& maxStrain);

	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf);
	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A);
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) override;

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) override;
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

	double computeGreenStrainNorm(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	void computeGreenStrainNormGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad);
	void computeGreenStrainNormHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes);

	double getThickness() { return m_h; }
	void setThickness(double thinckness) { m_h = thinckness; }

protected:
	double m_h;
};

#endif