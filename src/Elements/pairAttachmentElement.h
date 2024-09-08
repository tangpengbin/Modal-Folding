#pragma once
#ifndef PAIR_ATTACHMENT_ELEMENT_H
#define PAIR_ATTACHMENT_ELEMENT_H

#include "DofElementBase.h"
class pairAttachmentElement :
    public DofElementBase
{
public:
	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, double stiffness, bool applyElement = false);
	double getStiffness()
	{
		return m_stiffness;
	}

	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf);
	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A);
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad);
	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes);\
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);
private:
	double m_stiffness;
};

#endif