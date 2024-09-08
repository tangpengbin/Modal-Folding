#include "periodicDSHingeFixedRestPoseElement.h"
#include <SimLib/Core/SOUtils.h>
#include "../Codegen/computeDSBendingFixedRestPoseEnergy.h"
#include "../Codegen/computedTheta.h"
#include"../Codegen/computeCylindericalTransformVertex.h"


periodicDSHingeFixedRestPoseElement::periodicDSHingeFixedRestPoseElement(void)
{

}

periodicDSHingeFixedRestPoseElement::~periodicDSHingeFixedRestPoseElement(void)
{
}

void periodicDSHingeFixedRestPoseElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const std::vector<Eigen::Vector2d>& verticesTranslationMultipliers, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);

	m_applyElement = true;

	m_mat = pMat;
	m_thickness = thickness;

	m_verticesTranslationMultipliers = verticesTranslationMultipliers;

	m_theta0 = getCurrentAngle(vx);
	m_theta_initial = m_theta0;
	//printf("the rest angle is %f\n", m_theta0);

	m_dxdt0t1.setZero();
	for (int i = 0; i < 4; i++)
	{
		m_dxdt0t1(3 * i + 0, 0) = m_verticesTranslationMultipliers[i][0];
		m_dxdt0t1(3 * i + 1, 1) = m_verticesTranslationMultipliers[i][0];

		m_dxdt0t1(3 * i + 0, 2) = m_verticesTranslationMultipliers[i][1];
		m_dxdt0t1(3 * i + 1, 3) = m_verticesTranslationMultipliers[i][1];
	}

	m_bendingCurvature = 0.0;
	m_bendingDirection = Eigen::Vector2d(1.0, 0.0);
}

std::vector<V3D> periodicDSHingeFixedRestPoseElement::computePositions(const Eigen::VectorXd& vx)
{
	std::vector<V3D> x(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
	}

	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);

	for (int i = 0; i < 4; i++)
	{
		x[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0;
		x[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1;
	}
	return x;
}



void periodicDSHingeFixedRestPoseElement::transformVertices(std::vector<V3D>& x, std::vector<V3D>& X,
	const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
	const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{
	if (m_bendingCurvature == 0.0)
	{
		for (int i = 0; i < 4; i++)
		{
			x[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0;
			x[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1;

			X[i].head(2) += m_verticesTranslationMultipliers[i][0] * translation0_rest;
			X[i].head(2) += m_verticesTranslationMultipliers[i][1] * translation1_rest;
		}
	}
	else
	{
		for (int i = 0; i < 4; i++)
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

Eigen::Matrix<double, 12, 16> periodicDSHingeFixedRestPoseElement::transformVerticesJacobian(std::vector<V3D>& x, std::vector<V3D>& X,
	const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1,
	const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{
	Eigen::Matrix<double, 3 * 4, 3 * 4 + 4> jacobian;
	jacobian.setZero();
	if (m_bendingCurvature == 0.0)
	{
		jacobian.leftCols(12).setIdentity();
		jacobian.rightCols(4) = m_dxdt0t1;
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			Eigen::Vector2d input_translation0 = m_verticesTranslationMultipliers[i][0] * translation0;
			Eigen::Vector2d input_translation1 = m_verticesTranslationMultipliers[i][1] * translation1;
			Eigen::Vector3d rotate_x = x[i];
			Eigen::Matrix<double, 3, 7> jacobian_i;
			Codegen::computeCylindericalTransformVertexJacobian(rotate_x, input_translation0, input_translation1, m_bendingCurvature, m_bendingDirection, jacobian_i);

			jacobian.block(3 * i, 3 * i, 3, 3) = jacobian_i.leftCols(3);
			jacobian.block(3 * i, 12, 3, 4) = jacobian_i.rightCols(4);
		}
	}
	return jacobian;
}

std::vector<Eigen::Matrix<double, 12, 16>> periodicDSHingeFixedRestPoseElement::transformVerticesHessians(std::vector<V3D>& x, std::vector<V3D>& X, const Eigen::Vector2d& translation0, const Eigen::Vector2d& translation1, const Eigen::Vector2d& translation0_rest, const Eigen::Vector2d& translation1_rest)
{

	std::vector<Eigen::Matrix<double, 12, 16>> hessians(16);
	for (int i = 0; i < 16; i++)
		hessians[i].setZero();

	if (m_bendingCurvature == 0.0)
	{
	}
	else
	{
		for (int i = 0; i < 4; i++)
		{
			Eigen::Vector2d input_translation0 = m_verticesTranslationMultipliers[i][0] * translation0;
			Eigen::Vector2d input_translation1 = m_verticesTranslationMultipliers[i][1] * translation1;
			Eigen::Vector3d rotate_x = x[i];
			std::vector<Eigen::Matrix<double, 3, 7>> hessian_i(7);
			Codegen::computeCylindericalTransformVertexHessians(rotate_x, input_translation0, input_translation1, m_bendingCurvature, m_bendingDirection,
				hessian_i[0], hessian_i[1], hessian_i[2], hessian_i[3], hessian_i[4], hessian_i[5], hessian_i[6]);

			//dx^2
			soutil::setTensorBlock<3, 7, 12, 16>(hessian_i, 0, 0, 0, 3, 3, 3, hessians, 3 * i, 3 * i, 3 * i);

			//dx dt0t1
			soutil::setTensorBlock<3, 7, 12, 16>(hessian_i, 3, 0, 0, 4, 3, 3, hessians, 3 * 4, 3 * i, 3 * i);

			//dt0t1 dx
			soutil::setTensorBlock<3, 7, 12, 16>(hessian_i, 0, 0, 3, 3, 3, 4, hessians, 3 * i, 3 * i, 3 * 4);

			//dt0t1 dt0t1
			soutil::setTensorBlock<3, 7, 12, 16>(hessian_i, 3, 0, 3, 4, 3, 4, hessians, 3 * 4, 3 * i, 3 * 4);
		}
	}

	return hessians;
}


double periodicDSHingeFixedRestPoseElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);


	double k_bend = m_mat->kB;

	//double testAngle = computePeriodicHingeAngle(x[0], x[1], x[2], x[3], x[4], x[5]);

	double E = 0;
	Codegen::computeDSBendingFixedRestPoseEnergy(k_bend, x[0], x[1], x[2], x[3],
		X[0], X[1], X[2], X[3], m_theta0, E);

	return E;
}




void periodicDSHingeFixedRestPoseElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
{
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);


	double k_bend = m_mat->kB;

	
	Eigen::Matrix<double, 12, 1> fx;
	Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
		X[0], X[1], X[2], X[3], m_theta0, fx);

	Eigen::Matrix<double, 16, 1> gradient;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	gradient = x_jacobian.transpose() * fx;
	//gradient.head(12) = fx;
	//gradient.tail(4) = m_dxdt0t1.transpose() * fx;

	gradient *= -1.0;

	scatherAdd(gradient, vf);
}

void periodicDSHingeFixedRestPoseElement::computeGradient(const Eigen::VectorXd & vx, const Eigen::VectorXd & vX, Eigen::VectorXd & grad)
{
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	double k_bend = m_mat->kB;

	Eigen::Matrix<double, 12, 1> fx;
	Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
		X[0], X[1], X[2], X[3], m_theta0, fx);


	Eigen::Matrix<double, 16, 1> gradient;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	gradient = x_jacobian.transpose() * fx;
	//gradient.head(12) = fx;
	//gradient.tail(4) = m_dxdt0t1.transpose() * fx;


	grad = gradient;
}

void periodicDSHingeFixedRestPoseElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{


}

void periodicDSHingeFixedRestPoseElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
{
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	double k_bend = m_mat->kB;

	Eigen::Matrix<double, 12, 12> Jx;
	Codegen::computeDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2], x[3], 
		X[0], X[1], X[2], X[3], m_theta0, Jx);

	Eigen::Matrix<double, 16, 16> hessian;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);

	if (m_bendingCurvature == 0.0)
	{
		hessian = x_jacobian.transpose() * Jx * x_jacobian;
	}
	else
	{
		Eigen::Matrix<double, 12, 1> fx;
		Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);
		std::vector<Eigen::Matrix<double, 12, 16>> x_hessians = transformVerticesHessians(x, X, translation0, translation1, translation0_rest, translation1_rest);

		Eigen::Matrix<double, 16, 16> hessian2;
		soutil::vectorTransposeMulTensor<12, 16, 12, 16>(fx, x_hessians, hessian2);
		hessian = x_jacobian.transpose() * Jx * x_jacobian + hessian2;
	}
	/*hessian.topLeftCorner(12, 12) = Jx;
	Eigen::Matrix<double, 12, 4> dxdt0t1 = Jx * m_dxdt0t1;
	hessian.topRightCorner(12, 4) = dxdt0t1;
	hessian.bottomLeftCorner(4, 12) = dxdt0t1.transpose();
	hessian.bottomRightCorner(4, 4) = m_dxdt0t1.transpose() * Jx * m_dxdt0t1;*/

	hes = hessian;
}

void periodicDSHingeFixedRestPoseElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
{
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	double k_bend = m_mat->kB;


	
	Eigen::Matrix<double, 12, 12> Jx;
	Codegen::computeDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2], x[3], 
		X[0], X[1], X[2], X[3], m_theta0, Jx);


	Eigen::Matrix<double, 16, 16> hessian;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);

	if (m_bendingCurvature == 0.0)
	{
		hessian = x_jacobian.transpose() * Jx * x_jacobian;
	}
	else
	{
		Eigen::Matrix<double, 12, 1> fx;
		Codegen::computeDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2], x[3],
			X[0], X[1], X[2], X[3], m_theta0, fx);
		std::vector<Eigen::Matrix<double, 12, 16>> x_hessians = transformVerticesHessians(x, X, translation0, translation1, translation0_rest, translation1_rest);

		Eigen::Matrix<double, 16, 16> hessian2;
		soutil::vectorTransposeMulTensor<12, 16, 12, 16>(fx, x_hessians, hessian2);
		hessian = x_jacobian.transpose() * Jx * x_jacobian + hessian2;
	}
	/*hessian.topLeftCorner(12, 12) = Jx;
	Eigen::Matrix<double, 12, 4> dxdt0t1 = Jx * m_dxdt0t1;
	hessian.topRightCorner(12, 4) = dxdt0t1;
	hessian.bottomLeftCorner(4, 12) = dxdt0t1.transpose();
	hessian.bottomRightCorner(4, 4) = m_dxdt0t1.transpose() * Jx * m_dxdt0t1;*/


	if (isnan(hessian.norm()))
	{
		printf("hessian norm is nan");
	}

	scatherAdd(-hessian, triplets);
}


double periodicDSHingeFixedRestPoseElement::getCurrentAngle(const Eigen::VectorXd& vx)
{
	//compute rest theta

	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vx);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);



	double theta = 0.0;
	Codegen::computedTheta(x[0], x[1], x[2], x[3], theta);
	return theta;
}



void periodicDSHingeFixedRestPoseElement::computedThetadx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& dThetadx)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix<double, 12, 1> element_dThetadx;
	Codegen::computedThetaGradient(x[0], x[1], x[2], x[3],element_dThetadx);


	Eigen::Matrix<double, 16, 1> gradient;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);
	gradient = x_jacobian.transpose() * element_dThetadx;
	/*gradient.head(12) = element_dThetadx;
	gradient.tail(4) = m_dxdt0t1.transpose() * element_dThetadx;*/

	dThetadx = gradient;
}

void periodicDSHingeFixedRestPoseElement::computed2Thetadx2(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& d2Thetadx2)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(4);
	std::vector<V3D> X(4);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 4; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	Eigen::Vector2d translation0 = x_gather.segment<2>(12);
	Eigen::Vector2d translation1 = x_gather.segment<2>(14);
	Eigen::Vector2d translation0_rest = X_gather.segment<2>(12);
	Eigen::Vector2d translation1_rest = X_gather.segment<2>(14);

	transformVertices(x, X, translation0, translation1, translation0_rest, translation1_rest);

	Eigen::Matrix<double, 12, 12> element_d2Thetadx2;
	Codegen::computedThetaHessian(x[0], x[1], x[2], x[3], element_d2Thetadx2);


	Eigen::Matrix<double, 16, 16> hessian;
	Eigen::Matrix<double, 12, 16> x_jacobian = transformVerticesJacobian(x, X, translation0, translation1, translation0_rest, translation1_rest);

	if (m_bendingCurvature == 0.0)
	{
		hessian = x_jacobian.transpose() * element_d2Thetadx2 * x_jacobian;
	}
	else
	{
		Eigen::Matrix<double, 12, 1> element_dThetadx;
		Codegen::computedThetaGradient(x[0], x[1], x[2], x[3], element_dThetadx);
		std::vector<Eigen::Matrix<double, 12, 16>> x_hessians = transformVerticesHessians(x, X, translation0, translation1, translation0_rest, translation1_rest);

		Eigen::Matrix<double, 16, 16> hessian2;
		soutil::vectorTransposeMulTensor<12, 16, 12, 16>(element_dThetadx, x_hessians, hessian2);
		hessian = x_jacobian.transpose() * element_d2Thetadx2 * x_jacobian + hessian2;
	}
	/*hessian.topLeftCorner(12, 12) = element_d2Thetadx2;
	Eigen::Matrix<double, 12, 4> dxdt0t1 = element_d2Thetadx2 * m_dxdt0t1;
	hessian.topRightCorner(12, 4) = dxdt0t1;
	hessian.bottomLeftCorner(4, 12) = dxdt0t1.transpose();
	hessian.bottomRightCorner(4, 4) = m_dxdt0t1.transpose() * element_d2Thetadx2 * m_dxdt0t1;*/


	d2Thetadx2 = hessian;
}