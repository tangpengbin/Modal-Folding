#include "periodicCSTShellElement.h"
#include <SimLib/Core/SOUtils.h>
#include "../Codegen/computeCSTShellEnergy.h"
#include "../Codegen/computeCSTShellEnergyGreenStrain.h"

#include"../Codegen/computeCylindericalTransformVertex.h"

using namespace soutil;

periodicCSTShellElement::periodicCSTShellElement()
{
	m_h = 0.001;
}


periodicCSTShellElement::~periodicCSTShellElement()
{
}


void periodicCSTShellElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const std::vector<Eigen::Vector2d>& verticesTranslationMultipliers, const Eigen::VectorXd& vx, SOMaterial* pMat)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);
	m_applyElement = true;
	m_mat = pMat;
	
	m_verticesTranslationMultipliers = verticesTranslationMultipliers;

	m_dxdt0t1.setZero();
	for (int i = 0; i < 3; i++)
	{
		m_dxdt0t1(3 * i + 0, 0) = m_verticesTranslationMultipliers[i][0];
		m_dxdt0t1(3 * i + 1, 1) = m_verticesTranslationMultipliers[i][0];

		m_dxdt0t1(3 * i + 0, 2) = m_verticesTranslationMultipliers[i][1];
		m_dxdt0t1(3 * i + 1, 3) = m_verticesTranslationMultipliers[i][1];
	}

	m_bendingCurvature = 0.0;
	m_bendingDirection = Eigen::Vector2d(1.0, 0.0);
}

double periodicCSTShellElement::computeArea(const Eigen::VectorXd& vx)
{
	std::vector<V3D> X(3);

	Eigen::VectorXd X_gather = gatherVector(vx);
	X[0] = X_gather.segment<3>(0);
	X[1] = X_gather.segment<3>(3);
	X[2] = X_gather.segment<3>(6);

	Eigen::Vector2d translation0 = X_gather.segment<2>(9);
	Eigen::Vector2d translation1 = X_gather.segment<2>(11);

	for (int i = 0; i < 3; i++)
	{
		X[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0;
		X[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1;
	}


	V3D e1 = X[1] - X[0];
	V3D e2 = X[2] - X[0];
	V3D N = e1.cross(e2);
	double A = N.norm()*0.5;
	return A;
}

void periodicCSTShellElement::transformVertices(std::vector<V3D>& x, std::vector<V3D>& X,
	const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
	const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{
	if (m_bendingCurvature == 0.0)
	{
		for (int i = 0; i < 3; i++)
		{
			x[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0;
			x[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1;

			X[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0_rest;
			X[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1_rest;
		}
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector2d input_translation0 = m_verticesTranslationMultipliers[i][0] * translation0;
			Eigen::Vector2d input_translation1 = m_verticesTranslationMultipliers[i][1] * translation1;
			Eigen::Vector3d rotate_x = x[i];
			Codegen::computeCylindericalTransformVertex(rotate_x, input_translation0, input_translation1, m_bendingCurvature, m_bendingDirection, x[i]);


			X[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0_rest;
			X[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1_rest;
		}
	}

}

Eigen::Matrix<double, 9, 13> periodicCSTShellElement::transformVerticesJacobian(std::vector<V3D>& x, std::vector<V3D>& X,
	const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
	const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{
	Eigen::Matrix<double, 3 * 3, 3 * 3 + 4> jacobian;
	jacobian.setZero();
	if (m_bendingCurvature == 0.0)
	{
		jacobian.leftCols(9).setIdentity();
		jacobian.rightCols(4) = m_dxdt0t1;
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector2d input_translation0 = m_verticesTranslationMultipliers[i][0] * translation0;
			Eigen::Vector2d input_translation1 = m_verticesTranslationMultipliers[i][1] * translation1;
			Eigen::Vector3d rotate_x = x[i];
			Eigen::Matrix<double, 3, 7> jacobian_i;
			Codegen::computeCylindericalTransformVertexJacobian(rotate_x, input_translation0, input_translation1, m_bendingCurvature, m_bendingDirection, jacobian_i);

			jacobian.block(3 * i, 3 * i, 3, 3) = jacobian_i.leftCols(3);
			jacobian.block(3 * i, 9, 3, 4) = jacobian_i.rightCols(4);
		}
	}
	return jacobian;
}

std::vector<Eigen::Matrix<double, 9, 13>> periodicCSTShellElement::transformVerticesHessians(std::vector<V3D>& x, std::vector<V3D>& X, const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1, const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{
	
	std::vector<Eigen::Matrix<double, 9, 13>> hessians(13);
	for (int i = 0; i < 13; i++)
		hessians[i].setZero();
	
	if (m_bendingCurvature == 0.0)
	{
	}
	else
	{
		for (int i = 0; i < 3; i++)
		{
			Eigen::Vector2d input_translation0 = m_verticesTranslationMultipliers[i][0] * translation0;
			Eigen::Vector2d input_translation1 = m_verticesTranslationMultipliers[i][1] * translation1;
			Eigen::Vector3d rotate_x = x[i];
			std::vector<Eigen::Matrix<double, 3, 7>> hessian_i(7);
			Codegen::computeCylindericalTransformVertexHessians(rotate_x, input_translation0, input_translation1, m_bendingCurvature, m_bendingDirection, 
				hessian_i[0], hessian_i[1], hessian_i[2], hessian_i[3], hessian_i[4], hessian_i[5], hessian_i[6]);

			//dx^2
			soutil::setTensorBlock<3, 7, 9, 13>(hessian_i, 0, 0, 0, 3, 3, 3, hessians, 3 * i, 3 * i, 3 * i);

			//dx dt0t1
			soutil::setTensorBlock<3, 7, 9, 13>(hessian_i, 3, 0, 0, 4, 3, 3, hessians, 3 * 3, 3 * i, 3 * i);

			//dt0t1 dx
			soutil::setTensorBlock<3, 7, 9, 13>(hessian_i, 0, 0, 3, 3, 3, 4, hessians, 3 * i, 3 * i, 3 * 3);

			//dt0t1 dt0t1
			soutil::setTensorBlock<3, 7, 9, 13>(hessian_i, 3, 0, 3, 4, 3, 4, hessians, 3 * 3, 3 * i, 3 * 3);
		}
	}

	return hessians;
}

void periodicCSTShellElement::convertEandNuToLame(double E, double nu, double& lambda, double& mu)
{
	lambda = E * nu / ((1 + nu)*(1 - 2 * nu));
	mu = E / (2 * (1 + nu));
}

// computes jacobian of barycentric wrt p, its constant wrt to p so argument is dropped
inline Eigen::Matrix<double, 2, 3> barycentricJacobian(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c)
{
	Eigen::Vector3d v0 = b - a, v1 = c - a;
	//Vector3e v2 = p - a;
	double d00 = v0.dot(v0); // (b - a).dot(b - a) => d d00 / d p = 0
	double d01 = v0.dot(v1); // (b - a).dot(c - a) => d d01 / d p = 0
	double d11 = v1.dot(v1); // (c - a).dot(c - a) => d d11 / dp = 0
	//Expr d20 = v2.dot(v0); // (p - a).dot(b - a) => d d20 / dp = (b - a)^T
	//Expr d21 = v2.dot(v1); //  ( p - a).dot(c - a) = > d21 / dp = (c - a)^T
	double denom = d00 * d11 - d01 * d01; // => d00 * d11 constant in => drops out, d01 constant in p => derivative is 0
	//v = (d11 * d20 - d01 * d21) / denom;
	//w = (d00 * d21 - d01 * d20) / denom;
	//u = 1.0f - v - w;
	//Vector3e dvdp = (d11 * dd20 / dp - d01 * d d21 / dp) / denom;
	Eigen::Vector3d dvdp = (d11 * (b - a) - d01 * (c - a)) / denom;
	Eigen::Vector3d dwdp = (d00 * (c - a) - d01 * (b - a)) / denom;
	Eigen::Matrix<double, 2, 3> result;
	result.row(0) = dvdp.transpose();
	result.row(1) = dwdp.transpose();
	return result;
}

inline Eigen::Matrix<double, 3, 2> compute3DCSTDeformationGradient(
	const Eigen::Vector3d& x1Undef, const Eigen::Vector3d& x2Undef, const Eigen::Vector3d& x3Undef,
	const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const Eigen::Vector3d& x3)
{
	// defGrad = d x / d X
	// X(b) = X * N(b) = X * [ 1 - b1 - b2; b1; b2];
	// b(X) = [X2 - X1, X3 - X1]^-1 [X - X1]
	// x(X) = x * N(b(X)) then take the jacobian of this to get defgrad
	// d x / d X = x * d N / dX = x * dN/db * [X2 - X1, X3 - X1]^-1;
	// however here X means a 2 dimensional vector! so we get a 3x2 defgrad
	// x(X) = x * N(Barycentric(X, X1, X2, X3))
	// defGrad = dx / d X = x * dN/dB * dB /dX
	// that would work except that the defGradient is gonna be 3x3 and in the undef config. something like [1, 0, 0;0,0,0;0,0,1];
	// and then E = 0.5 * (F^T F - I)  is gonna give a non zero energy in the undef config.
	//instead we choose a 2D coordinate system X* in the undef configuration for which we compute the def grad
	// defGrad = d x / d X*
	// x(X*) = x * N(Barycentric(X(X*)));
	// X(X*) = X1 + t * X*[0] + q *X*[1]

	Eigen::Vector3d tUndef = (x2Undef - x1Undef).normalized();
	Eigen::Vector3d e2Undef = (x3Undef - x1Undef);
	Eigen::Vector3d qUndef = (e2Undef - tUndef * e2Undef.dot(tUndef)).normalized();

	Eigen::Matrix3d x;
	x << x1, x2, x3;

	//N(b) = [1 - b1 - b2, b1, b2]
	Eigen::Matrix<double, 3, 2> dNdb;
	dNdb << -1.0, -1.0,
		1.0, 0.0,
		0.0, 1.0;

	Eigen::Matrix<double, 2, 3> dBdX = barycentricJacobian(x1Undef, x2Undef, x3Undef);
	Eigen::Matrix<double, 3, 2> dXdXStar;
	dXdXStar << tUndef, qUndef;
	Eigen::Matrix<double, 3, 2> defGrad = x * dNdb * dBdX * dXdXStar; //note that this F is not very intuitive it can contain -1 for undef configuration, but its not a problem as long as only F^T*F is used
	return defGrad;
}

inline Eigen::Matrix2d computeGreenStrain(const Eigen::Matrix<double, 3, 2>& defGrad)
{
	Eigen::Matrix2d greenStrain = 0.5 * (defGrad.transpose() * defGrad - Eigen::Matrix2d::Identity());
	return greenStrain;
}


auto computeTestEnergy(
	double youngsModulus,
	double poissonsRatio, double thickness,
	const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const Eigen::Vector3d& x3,
	const Eigen::Vector3d& x1Undef3D, const Eigen::Vector3d& x2Undef3D, const Eigen::Vector3d& x3Undef3D)
{

	double stretchStiffness = 0.5 * youngsModulus * thickness / (1.0 - poissonsRatio * poissonsRatio);
	Eigen::Vector2d quadPoint;
	quadPoint << 1.0 / 3.0, 1.0 / 3.0;

	double undefArea = 0.5 * (x2Undef3D - x1Undef3D).cross(x3Undef3D - x1Undef3D).norm();

	Eigen::Matrix<double, 3, 2> defGrad = compute3DCSTDeformationGradient(
		x1Undef3D, x2Undef3D, x3Undef3D,
		x1, x2, x3
	);
	Eigen::Matrix2d greenStrain = computeGreenStrain(defGrad);
	//MatrixXe xUndef(3, 3);
	//MatrixXe x(3, 3);
	//xUndef << x1Undef3D, x2Undef3D, x3Undef3D;
	//x << x1, x2, x3;
	//Matrix3e greenStrain = computeShellMembraneStrainTensor(xUndef, x);

	double trE = greenStrain.trace();
	double trESq = greenStrain.squaredNorm();

	double energy = stretchStiffness * undefArea * ((1.0 - poissonsRatio) * trESq + poissonsRatio * (trE * trE));

	//Expr lambda = youngsModulus * poissonsRatio / (1.0 + poissonsRatio) * (1.0 - 2.0 * poissonsRatio);
	//Expr mu = youngsModulus / (2.0 * (1.0 + poissonsRatio));
	//Expr strainEnergyDensity = lambda * 0.5 * (trE * trE) + mu * trESq;
	//Expr energy = strainEnergyDensity * thickness * undefArea;
	return energy;
};



double periodicCSTShellElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	double W = 0.0;

	Codegen::computeCSTShellEnergy(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], W);
	
	return W;
}

void periodicCSTShellElement::computeMaxStrain(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, double& maxStrain)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix2d greenStrain;

	Codegen::computeCSTShellEnergyGreenStrain(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], greenStrain);
	
	Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(greenStrain);
	if (eigensolver.info() != Eigen::Success) abort();
	Eigen::VectorXd eige3alues = eigensolver.eigenvalues();
	eige3alues.cwiseAbs();
	maxStrain = eige3alues.maxCoeff();
}

void periodicCSTShellElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix<double, 9, 1> fx;
	
	Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], fx);

	Eigen::Matrix<double, 13, 1> gradient;
	Eigen::Matrix<double, 9,13> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	gradient = x_jacobian.transpose() * fx;
	//gradient.head(9) = fx;
	//gradient.tail(4) = m_dxdt0t1.transpose() * fx;


	gradient *= -1.0;//force

	scatherAdd(gradient, vf);
}

void periodicCSTShellElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{}

void periodicCSTShellElement::computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix<double, 9, 1> fx;
	Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], fx);
	

	Eigen::Matrix<double, 13, 1> gradient;
	Eigen::Matrix<double, 9, 13> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	gradient = x_jacobian.transpose() * fx;
	//gradient.head(9) = fx;
	//gradient.tail(4) = m_dxdt0t1.transpose() * fx;


	grad = gradient;
}

void periodicCSTShellElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix<double, 9, 9> Jx;
	Codegen::computeCSTShellEnergyHessian(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], Jx);
	
	Eigen::Matrix<double, 13, 13> hessian;
	Eigen::Matrix<double, 9, 13> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	
	if (m_bendingCurvature == 0.0)
	{
		hessian = x_jacobian.transpose() * Jx * x_jacobian;
	}
	else
	{
		Eigen::Matrix<double, 9, 1> fx;
		Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], fx);
		std::vector<Eigen::Matrix<double, 9, 13>> x_hessians = transformVerticesHessians(x, X, translation0, translation1, translation0_rest, translation1_rest);



		Eigen::Matrix<double, 13, 13> hessian2;
		soutil::vectorTransposeMulTensor<9, 13, 9, 13>(fx, x_hessians, hessian2);
		hessian = x_jacobian.transpose() * Jx * x_jacobian + hessian2;
	}
	/*hessian.topLeftCorner(9, 9) = Jx;
	Eigen::Matrix<double, 9, 4> dxdt0t1 = Jx * m_dxdt0t1;
	hessian.topRightCorner(9, 4) = dxdt0t1;
	hessian.bottomLeftCorner(4, 9) = dxdt0t1.transpose();
	hessian.bottomRightCorner(4, 4) = m_dxdt0t1.transpose() * Jx * m_dxdt0t1;*/

	hes = hessian;
}

void periodicCSTShellElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
{
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(9);
	Eigen::Vector2d translation1 = x_gather.segment<2>(11);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(9);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(11);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	
	Eigen::Matrix<double, 9, 9> Jx;
	Codegen::computeCSTShellEnergyHessian(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], Jx);

	Eigen::Matrix<double, 13, 13> hessian;
	Eigen::Matrix<double, 9, 13> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);

	if (m_bendingCurvature == 0.0)
	{
		hessian = x_jacobian.transpose() * Jx * x_jacobian;
	}
	else
	{
		Eigen::Matrix<double, 9, 1> fx;
		Codegen::computeCSTShellEnergyGradient(m_mat->k1, m_mat->k2, m_h, x[0], x[1], x[2], X[0], X[1], X[2], fx);
		std::vector<Eigen::Matrix<double, 9, 13>> x_hessians = transformVerticesHessians(x, X, translation0, translation1, translation0_rest, translation1_rest);



		Eigen::Matrix<double, 13, 13> hessian2;
		soutil::vectorTransposeMulTensor<9, 13, 9, 13>(fx, x_hessians, hessian2);
		hessian = x_jacobian.transpose() * Jx * x_jacobian + hessian2;
	}
	/*hessian.topLeftCorner(9, 9) = Jx;
	Eigen::Matrix<double, 9, 4> dxdt0t1 = Jx * m_dxdt0t1;
	hessian.topRightCorner(9, 4) = dxdt0t1;
	hessian.bottomLeftCorner(4, 9) = dxdt0t1.transpose();
	hessian.bottomRightCorner(4, 4) = m_dxdt0t1.transpose() * Jx * m_dxdt0t1;*/

	hessian *= -1.0;//dforce / dx

	scatherAdd(hessian, triplets);
}