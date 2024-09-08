#pragma once
#ifndef TARGET_PLANAR_PARAMETERS_ELEMENT_H
#define TARGET_PLANAR_PARAMETERS_ELEMENT_H
#include <SimLib/Elements/DofElementBase.h>
class targetPlanarParametersElement :
	public DofElementBase
{
public:
	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs,
		double penaltyStiffness,
		Eigen::Vector2d targetTranslation,
		bool applyElement = true);

	void setTargetTranslation(Eigen::Vector2d targetTranslation) 
	{
		m_targetTranslation = targetTranslation;
	}
	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);

	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf);
	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A);
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad);
	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes);
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

private:
	double m_penaltyStiffness;

	Eigen::Vector2d m_targetTranslation;
};

#endif