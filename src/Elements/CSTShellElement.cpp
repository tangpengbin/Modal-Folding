#include "CSTShellElement.h"
#include <SimLib/Core/SOUtils.h>
#include "../Codegen/computeCSTShellEnergy.h"
#include "../Codegen/computeCSTShellEnergyGreenStrain.h"

#include "../Codegen/computeCSTShellEnergyGreenStrainNorm.h"
using namespace soutil;
bool USE_THICKNESS = true;
const int NV = 3;

CSTShellElement::CSTShellElement()
{
	m_h = 0.001;
}


CSTShellElement::~CSTShellElement()
{
}


void CSTShellElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat)
{
	assert(dofIdcs.size() == NV * 3);
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);
	m_applyElement = true;
	m_mat = pMat;
	
}

double CSTShellElement::computeArea(const Eigen::VectorXd& vx)
{
	std::vector<V3D> X(NV);

	Eigen::VectorXd X_gather = gatherVector(vx);
	V3D X1 = X_gather.segment<3>(0);
	V3D X2 = X_gather.segment<3>(3);
	V3D X3 = X_gather.segment<3>(6);

	V3D e1 = X2 - X1;
	V3D e2 = X3 - X1;
	V3D N = e1.cross(e2);
	double A = N.norm()*0.5;
	return A;
}

void CSTShellElement::convertEandNuToLame(double E, double nu, double& lambda, double& mu)
{
	lambda = E * nu / ((1 + nu)*(1 - 2 * nu));
	mu = E / (2 * (1 + nu));
}

double CSTShellElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
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
	double W = 0.0;

	Codegen::computeCSTShellEnergy(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], W);
	
	
	return W;
}

void CSTShellElement::computeMaxStrain(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, double& maxStrain)
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
	Eigen::Matrix2d greenStrain;

	Codegen::computeCSTShellEnergyGreenStrain(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], greenStrain);
	
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(greenStrain);
	if (eigensolver.info() != Eigen::Success) abort();
	Eigen::VectorXd eigenvalues = eigensolver.eigenvalues();
	eigenvalues.cwiseAbs();
	maxStrain = eigenvalues.maxCoeff();
}

void CSTShellElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
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

	Eigen::Matrix<double, 9, 1> f;
	
	Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], f);
	
	f *= -1.0;//force

	scatherAdd(f, vf);
}

void CSTShellElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{}

void CSTShellElement::computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad)
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

	Eigen::Matrix<double, 9, 1> f;
	
	Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], f);
	grad = f;
}

void CSTShellElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
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


	Eigen::Matrix<double, 9, 9> J;
	
	Codegen::computeCSTShellEnergyHessian(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], J);
	hes = J;
}

void CSTShellElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
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

	
	Eigen::Matrix<double, 9, 9> hessian;
	
	Codegen::computeCSTShellEnergyHessian(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], hessian);
	
	hessian *= -1.0;//dforce / dx

	scatherAdd(hessian, triplets);
}


double CSTShellElement::computeGreenStrainNorm(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
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

	double greenStrainSquaredNorm = 0.0;
	Codegen::computeCSTShellEnergyGreenStrainNorm(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], greenStrainSquaredNorm);

	return greenStrainSquaredNorm;
}
void CSTShellElement::computeGreenStrainNormGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad)
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

	Eigen::Matrix<double, 9, 1> g;
	Codegen::computeCSTShellEnergyGreenStrainNormGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], g);

	grad = g;
}
void CSTShellElement::computeGreenStrainNormHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
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

	Eigen::Matrix<double, 9, 9> h;
	Codegen::computeCSTShellEnergyGreenStrainNormHessian(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], h);

	hes = h;
}
