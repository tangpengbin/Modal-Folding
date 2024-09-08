#include "DSHingeFixedRestPoseElement.h"
#include "../Codegen/computeDSBendingFixedRestPoseEnergy.h"
#include "../Codegen/computeDSBendingFixedRestPoseCorrectEnergy.h"
#include "../Codegen/computedTheta.h"
#include <iostream>

const int NV = 4;
bool useCorrectEnergy = true;
DSHingeFixedRestPoseElement::DSHingeFixedRestPoseElement(void)
{

}

DSHingeFixedRestPoseElement::~DSHingeFixedRestPoseElement(void)
{
}

void DSHingeFixedRestPoseElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness, bool hardEdge)
{
	assert(dofIdcs.size() == NV * 3);
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);

	m_applyElement = true;

	m_mat = pMat;
	m_thickness = thickness;
	
	m_theta0 = getCurrentAngle(vx);
	m_theta_initial = m_theta0;
	m_hardEdge = hardEdge;
	//printf("the rest angle is %f\n", m_theta0);
}

double DSHingeFixedRestPoseElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);

	double E = 0;
	if (useCorrectEnergy)
	{
		Codegen::computeDSBendingFixedRestPoseCorrectEnergy(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, E);
		/*double Eold = 0;
		Codegen::computeDSBendingFixedRestPoseEnergy(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, Eold);
		double times = E/ Eold;
		printf("Times %f", times);
		exit(0);*/
	}
	else
	{
		Codegen::computeDSBendingFixedRestPoseEnergy(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, E);

		
	}
	//E /= 1800.0;
	return E;
}




void DSHingeFixedRestPoseElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);

	
	Eigen::Matrix<double, 12, 1> fx;
	if (useCorrectEnergy)
		Codegen::computeDSBendingFixedRestPoseCorrectEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);
	else
		Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);
	fx *= -1.0;
	

	//fx /= 1800.0;
	scatherAdd(fx, vf);
}

void DSHingeFixedRestPoseElement::computeGradient(const Eigen::VectorXd & vx, const Eigen::VectorXd & vX, Eigen::VectorXd & grad)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);

	computeEnergy(vx, vX);

	Eigen::Matrix<double, 12, 1> fx;
	if (useCorrectEnergy)
		Codegen::computeDSBendingFixedRestPoseCorrectEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);
	else
		Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);

	//fx /= 1800.0;
	grad = fx;
}

void DSHingeFixedRestPoseElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{


}

void DSHingeFixedRestPoseElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);

	Eigen::Matrix<double, 12, 12> H;
	if (useCorrectEnergy)
		Codegen::computeDSBendingFixedRestPoseCorrectEnergyHessian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, H);
	else
		Codegen::computeDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, H);

	//H /= 1800.0;
	hes = H;
}

void DSHingeFixedRestPoseElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);


	
	Eigen::Matrix<double, 12, 12> H;
	if (useCorrectEnergy)
		Codegen::computeDSBendingFixedRestPoseCorrectEnergyHessian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, H);
	else
		Codegen::computeDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, H);

	//H /= 1800.0;
	scatherAdd(-H, triplets);
}

void DSHingeFixedRestPoseElement::computeGradientParameterJacobian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& gradientParameterJacobian)
{
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	double k_bend = m_mat->kB * (m_hardEdge ? 100.0 : 1.0);

	Eigen::Matrix<double, 12, 1> Jx;
	if (useCorrectEnergy)
		Codegen::computeDSBendingFixedRestPoseCorrectEnergyGradientParameterJacobian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, Jx);
	else
		Codegen::computeDSBendingFixedRestPoseEnergyGradientParameterJacobian(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, Jx);

	//Jx /= 1800.0;
	gradientParameterJacobian = Jx;
}


double DSHingeFixedRestPoseElement::getCurrentAngle(const Eigen::VectorXd& vx)
{
	//compute rest theta
	std::vector<V3D> x(NV);
	Eigen::VectorXd x_gather = gatherVector(vx);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
	}
	double theta = 0.0;
	Codegen::computedTheta(x[0], x[1], x[2], x[3], theta);
	
	return theta;
}



void DSHingeFixedRestPoseElement::computedThetadx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& dThetadx)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	Eigen::Matrix<double, 12, 1> element_dThetadx;
	Codegen::computedThetaGradient(x[0], x[1], x[2], x[3], element_dThetadx);
	
	dThetadx = element_dThetadx;
}

void DSHingeFixedRestPoseElement::computed2Thetadx2(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& d2Thetadx2)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	Eigen::Matrix<double, 12, 12> element_d2Thetadx2;
	Codegen::computedThetaHessian(x[0], x[1], x[2], x[3], element_d2Thetadx2);

	d2Thetadx2 = element_d2Thetadx2;
}