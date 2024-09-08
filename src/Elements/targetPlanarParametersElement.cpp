#include "targetPlanarParametersElement.h"
#include <SimLib/Core/SOUtils.h>

void targetPlanarParametersElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs,
	double penaltyStiffness,
	Eigen::Vector2d targetTranslation,
	bool applyElement)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);
	m_applyElement = applyElement;
	m_penaltyStiffness = penaltyStiffness;

	//default
	m_targetTranslation = targetTranslation;
}


double targetPlanarParametersElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	if (!m_applyElement)
		return 0.0;

	Eigen::VectorXd t = gatherVector(vx);
	Eigen::Vector2d currentTranslation = t.head(2);


	double energy = 0.0;
	energy += 0.5 * m_penaltyStiffness * (currentTranslation - m_targetTranslation).squaredNorm();
	return energy;
}

void targetPlanarParametersElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf) {
	if (!m_applyElement)
		return;

	Eigen::VectorXd t = gatherVector(vx);
	Eigen::Vector2d currentTranslation = t.head(2);


	Eigen::Vector2d gradient; gradient.setZero();
	gradient = m_penaltyStiffness * (currentTranslation - m_targetTranslation);
	

	gradient *= -1.0;

	scatherAdd(gradient, vf);
}

void targetPlanarParametersElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A) {
	throw std::logic_error("targetPlanarParametersElement addHessian Element not implemented");
}

void targetPlanarParametersElement::computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) {
	if (!m_applyElement)
	{
		grad = Eigen::VectorXd(m_dofIdcs.size());
		grad.setZero();
		return;
	}

	Eigen::VectorXd t = gatherVector(vx);
	Eigen::Vector2d currentTranslation = t.head(2);


	Eigen::Vector2d gradient; gradient.setZero();
	gradient = m_penaltyStiffness * (currentTranslation - m_targetTranslation);


	grad = gradient;
}

void targetPlanarParametersElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) {
	if (!m_applyElement)
	{
		hes = Eigen::MatrixXd(m_dofIdcs.size(), m_dofIdcs.size());
		hes.setZero();
		return;
	}
	Eigen::VectorXd t = gatherVector(vx);
	Eigen::Vector2d currentTranslation = t.head(2);


	Eigen::Matrix2d hessian; hessian.setZero();
	hessian(0, 0) = m_penaltyStiffness;
	hessian(1, 1) = m_penaltyStiffness;
	

	
	hes = hessian;
}

void targetPlanarParametersElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets) {
	if (!m_applyElement)
		return;

	Eigen::VectorXd t = gatherVector(vx);
	Eigen::Vector2d currentTranslation = t.head(2);


	Eigen::Matrix2d hessian; hessian.setZero();
	hessian(0, 0) = m_penaltyStiffness;
	hessian(1, 1) = m_penaltyStiffness;


	hessian *= -1.0;

	scatherAdd(hessian, triplets);
}
