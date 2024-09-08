#include "pairAttachmentElement.h"
#include "../Codegen/computeSpringEnergy.h"
#include <iostream>

void pairAttachmentElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, double stiffness, bool applyElement)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);
	m_stiffness = stiffness;
	m_applyElement = applyElement;
}


double pairAttachmentElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	if(!m_applyElement)
		return 0.0;

	Eigen::VectorXd x = gatherVector(vx);
	Eigen::VectorXd X = gatherVector(vX);

    double energy = 0.0;// = 0.5 * m_stiffness * (x.segment<3>(0) - x.segment<3>(3)).squaredNorm();
    Codegen::computeSpringEnergy(m_stiffness, x.segment<3>(0), x.segment<3>(3), energy);

	return energy;
}

void pairAttachmentElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
{
	if (!m_applyElement)
		return;

	Eigen::VectorXd x = gatherVector(vx);
	Eigen::VectorXd X = gatherVector(vX);

    Eigen::Matrix<double, 6, 1> grad; 
    Codegen::computeSpringEnergyGradient(m_stiffness, x.segment<3>(0), x.segment<3>(3), grad);

    grad *= -1.0;

    scatherAdd(grad, vf);
}

void pairAttachmentElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat& A)
{
	throw std::logic_error("addHessian not implemented for pairAttachmentElement");
	if (!m_applyElement)
		return;
}

void pairAttachmentElement::computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad)
{
	if (!m_applyElement)
	{
        grad = Eigen::Matrix<double, 6, 1>();
		grad.setZero();
		return;
	}

    Eigen::VectorXd x = gatherVector(vx);
    Eigen::VectorXd X = gatherVector(vX);

    Eigen::Matrix<double, 6, 1> gradient;
    Codegen::computeSpringEnergyGradient(m_stiffness, x.segment<3>(0), x.segment<3>(3), gradient);
    grad = gradient;

}

void pairAttachmentElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
{
	if (!m_applyElement)
	{
        hes = Eigen::Matrix<double, 6, 6>();
		hes.setZero();
		return;
	}

    Eigen::VectorXd x = gatherVector(vx);
    Eigen::VectorXd X = gatherVector(vX);


    Eigen::Matrix<double, 6, 6> hessian;
    Codegen::computeSpringEnergyHessian(m_stiffness, x.segment<3>(0), x.segment<3>(3), hessian);
    hes = hessian;

}

void pairAttachmentElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
{
	if (!m_applyElement)
		return;

	Eigen::VectorXd x = gatherVector(vx);
	Eigen::VectorXd X = gatherVector(vX);

	Eigen::Matrix<double, 6, 6> hessian;
	Codegen::computeSpringEnergyHessian(m_stiffness, x.segment<3>(0), x.segment<3>(3), hessian);
	hessian *= -1.0;
	
	scatherAdd(hessian, triplets);
}
