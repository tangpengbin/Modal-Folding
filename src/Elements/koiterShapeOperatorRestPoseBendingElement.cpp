#include "koiterShapeOperatorRestPoseBendingElement.h"
#include "../Codegen/computeKoiterShapeOperatorRestPoseBendingEnergy.h"
#include "../Codegen/computeShapeOperator.h"

#include "../Codegen/computeOrthotropicShapeOperatorRestPoseBendingEnergy.h"
#include <iostream>

koiterShapeOperatorRestPoseBendingElement::koiterShapeOperatorRestPoseBendingElement(void)
{
	orthotropicBending = false;
}

koiterShapeOperatorRestPoseBendingElement::~koiterShapeOperatorRestPoseBendingElement(void)
{
}

Eigen::Vector3d reflectPoint(const Eigen::Vector3d& A, const Eigen::Vector3d& B, const Eigen::Vector3d& C)
{
	Eigen::Vector3d AB = B - A;
	Eigen::Vector3d AC = C - A;

	//D the projection of A to line AB
	Eigen::Vector3d AD = AB.dot(AC) / AB.norm() * AB;

	Eigen::Vector3d CD = -AC + AD;

	return C + 2 * CD;
}

void koiterShapeOperatorRestPoseBendingElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);

	m_applyElement = true;

	m_mat = pMat;
	m_thickness = thickness;
	
	int NV = dofIdcs.size() / 3;
	binary_multiplier.setZero();
	if (NV == 4)
		binary_multiplier[0] = 1.0;
	else if (NV == 5)
	{
		binary_multiplier[0] = 1.0;
		binary_multiplier[1] = 1.0;
	}
	else
		binary_multiplier.setOnes();

	m_shapeOperator0 = getCurrentShapeOperator(vx, vx);
	m_shapeOperator_initial = m_shapeOperator0;
	//printf("the rest angle is %f\n", m_theta0);



	double yongsModulus = m_mat->k1;
	double poissonsRatio = m_mat->k2;
	//Y01 = v01 * Y1 = v10 * Y0
	//Y0 = 2e6;   v01 = 0.4
	//Y1 = 1e6;   v10 = 0.2
	//Y01 = 4e5
	//G01 = 5e5

	//version1
	Y1 = yongsModulus;
	v01 = poissonsRatio;
	Y0 = Y1 / 10.0;
	v10 = v01 * Y1 / Y0;
	Y01 = v01 * Y1;
	G01 = Y1 / 8.0;

	//version2
	/*Y1 = yongsModulus;
	v01 = poissonsRatio;
	Y0 = Y1 / 10.0;
	v10 = v01 * Y1 / Y0;
	Y01 = v01 * Y1;
	G01 = Y1 / 5.0;*/

	//version3
	/*Y1 = yongsModulus;
	v01 = poissonsRatio;
	Y0 = Y1 / 10.0;
	v10 = v01 * Y1 / Y0;
	Y01 = v01 * Y1;
	G01 = Y1 / 2.5;*/
}


auto computeTriangleNormal(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3)
{
	Eigen::Vector3d v12 = v2 - v1;
	Eigen::Vector3d v13 = v3 - v1;

	Eigen::Vector3d normal = v12.cross(v13);
	return normal.normalized();
}

auto AngleSigned(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& n)
{
	//Vector3e m = a.cross(n);

	//Expr y = -m.dot(b);
	//Expr x = a.dot(b);
	//return atan2(y, x);

	return double(2.0) * atan(a.cross(b).dot(n) / (double(1.0) + (a.dot(b))));
}



Eigen::Matrix2d computeTriangleAveragedShapeOperator(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2, const Eigen::Vector3d& v3,
	const Eigen::Vector3d& v4, const Eigen::Vector3d& v5, const Eigen::Vector3d& v6,
	const Eigen::Vector3d& undeformed_v1, const Eigen::Vector3d& undeformed_v2, const Eigen::Vector3d& undeformed_v3,
	double& theta1, double& theta2, double& theta3)
{
	Eigen::Vector3d n = computeTriangleNormal(v1, v2, v3);

	Eigen::Vector3d ne1 = computeTriangleNormal(v2, v4, v3);
	Eigen::Vector3d ne2 = computeTriangleNormal(v1, v3, v5);
	Eigen::Vector3d ne3 = computeTriangleNormal(v1, v6, v2);

	Eigen::Vector3d V32 = undeformed_v3 - undeformed_v2;
	Eigen::Vector3d V13 = undeformed_v1 - undeformed_v3;
	Eigen::Vector3d V21 = undeformed_v2 - undeformed_v1;

	Eigen::Vector2d ta1(V32[1], -V32[0]);
	Eigen::Vector2d ta2(V13[1], -V13[0]);
	Eigen::Vector2d ta3(V21[1], -V21[0]);

	theta1 = AngleSigned(n, ne1, (v3 - v2).normalized());
	theta2 = AngleSigned(n, ne2, (v1 - v3).normalized());
	theta3 = AngleSigned(n, ne3, (v2 - v1).normalized());

	double area = 0.5 * (V21.cross(V13)).norm();
	Eigen::Matrix2d shape1 = (theta1 / 2.0 / (area * V32.norm())) * ta1 * ta1.transpose();
	Eigen::Matrix2d shape2 = (theta2 / 2.0 / (area * V13.norm())) * ta2 * ta2.transpose();
	Eigen::Matrix2d shape3 = (theta3 / 2.0 / (area * V21.norm())) * ta3 * ta3.transpose();

	Eigen::Matrix2d symmetricShapeOperator = shape1 + shape2 + shape3;
	return symmetricShapeOperator;
}


Eigen::Matrix2d koiterShapeOperatorRestPoseBendingElement::computeAngleAndShapeOperator(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX,
	double& theta1, double& theta2, double& theta3)
{
	int NV = m_dofIdcs.size() / 3;

	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	return computeTriangleAveragedShapeOperator(x[0], x[1], x[2], x[3], x[4], x[5],
		X[0], X[1], X[2],
		theta1, theta2, theta3);
}

//the order is following, we will duplicate vertices
//		 x3
//     /   \
//    x2---x1
//   /  \ /  \
//  x4---x0---x5

double koiterShapeOperatorRestPoseBendingElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	int NV = m_dofIdcs.size() / 3;

	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}


	double E = 0;
	if (orthotropicBending)
		Codegen::computeOrthotropicShapeOperatorRestPoseBendingEnergy(Y0, Y1, Y01, G01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, E);
	else
		Codegen::computeKoiterShapeOperatorRestPoseBendingEnergy(Y1, v01, m_thickness,
			x[0], x[1], x[2], x[3], x[4],x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, E);

	return E;
}

void koiterShapeOperatorRestPoseBendingElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
{
	int NV = m_dofIdcs.size() / 3;
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}
	
	double yongsModulus = m_mat->k1;
	double poissonsRatio = m_mat->k2;

	Eigen::Matrix<double, 18, 1> fx;
	if (orthotropicBending)
		Codegen::computeOrthotropicShapeOperatorRestPoseBendingEnergyGradient(Y0, Y1, Y01, G01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, fx);
	else
		Codegen::computeKoiterShapeOperatorRestPoseBendingEnergyGradient(Y1, v01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, fx);
	
	fx *= -1.0;
	
	scatherAdd(fx.head(NV*3), vf);
}

void koiterShapeOperatorRestPoseBendingElement::computeGradient(const Eigen::VectorXd & vx, const Eigen::VectorXd & vX, Eigen::VectorXd & grad)
{
	int NV = m_dofIdcs.size() / 3;
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}


	double yongsModulus = m_mat->k1;
	double poissonsRatio = m_mat->k2;

	Eigen::Matrix<double, 18, 1> fx;
	if (orthotropicBending)
		Codegen::computeOrthotropicShapeOperatorRestPoseBendingEnergyGradient(Y0, Y1, Y01, G01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, fx);
	else
		Codegen::computeKoiterShapeOperatorRestPoseBendingEnergyGradient(Y1, v01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, fx);

	grad = fx.head(NV * 3);
	
}

void koiterShapeOperatorRestPoseBendingElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{


}

void koiterShapeOperatorRestPoseBendingElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
{
	int NV = m_dofIdcs.size() / 3;
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}

	Eigen::Matrix<double, 18, 18> H;
	if (orthotropicBending)
		Codegen::computeOrthotropicShapeOperatorRestPoseBendingEnergyHessian(Y0, Y1, Y01, G01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, H);
	else
		Codegen::computeKoiterShapeOperatorRestPoseBendingEnergyHessian(Y1, v01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, H);

	hes = H.topLeftCorner(3 * NV, 3 * NV);
	
}

void koiterShapeOperatorRestPoseBendingElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
{
	int NV = m_dofIdcs.size() / 3;
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}

	Eigen::Matrix<double, 18, 18> H;
	if (orthotropicBending)
		Codegen::computeOrthotropicShapeOperatorRestPoseBendingEnergyHessian(Y0, Y1, Y01, G01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, H);
	else
		Codegen::computeKoiterShapeOperatorRestPoseBendingEnergyHessian(Y1, v01, m_thickness,
			x[0], x[1], x[2], x[3], x[4], x[5],
			X[0], X[1], X[2], m_shapeOperator0, binary_multiplier, H);

	//H /= 1800.0;
	scatherAdd(-H.topLeftCorner(NV * 3, NV * 3), triplets);
}


Eigen::Matrix2d koiterShapeOperatorRestPoseBendingElement::getCurrentShapeOperator(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
{
	int NV = m_dofIdcs.size() / 3;
	//compute rest theta
	std::vector<V3D> x(NV);
	std::vector<V3D> X(NV);
	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vx);
	for (int i = 0; i < NV; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}

	Eigen::Matrix2d shapeOperator;
	Codegen::computeShapeOperator(x[0], x[1], x[2], x[3], x[4], x[5],
		X[0], X[1], X[2], binary_multiplier, shapeOperator);
	
	return shapeOperator;
}



void koiterShapeOperatorRestPoseBendingElement::computedShapeOperatordx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, std::vector<Eigen::Matrix2d>& dShapeOperatordx)
{
	int NV = m_dofIdcs.size() / 3;
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
	if (NV == 4)
	{//we will add two duplicate vertices
		Eigen::Vector3d x4 = reflectPoint(x[0], x[2], x[1]); //reflect x1 to obtain x4
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x4);
		x.push_back(x5);
	}
	else if (NV == 5)
	{//we will add one duplicate vertex
		Eigen::Vector3d x5 = reflectPoint(x[0], x[1], x[2]); //reflect x2 to obtain x5
		x.push_back(x5);
	}

	std::vector<Eigen::Matrix2d> Ejacobian(6 * 3);
	Codegen::computeShapeOperatorjacobians(x[0], x[1], x[2], x[3], x[4], x[5],
		X[0], X[1], X[2], binary_multiplier,
		Ejacobian[0], Ejacobian[1], Ejacobian[2], Ejacobian[3], Ejacobian[4], Ejacobian[5], 
		Ejacobian[6], Ejacobian[7], Ejacobian[8], Ejacobian[9], Ejacobian[10], Ejacobian[11], 
		Ejacobian[12], Ejacobian[13], Ejacobian[14], Ejacobian[15], Ejacobian[16], Ejacobian[17]);
	
	if (NV == 4)
		Ejacobian.erase(Ejacobian.end() - 6, Ejacobian.end());
	else if (NV == 5)
		Ejacobian.erase(Ejacobian.end() - 3, Ejacobian.end());

	dShapeOperatordx = Ejacobian;
}
