#include "reflectionDSHingeFixedRestPoseElement.h"
#include "../Codegen/computeReflectionXPlaneDSBendingFixedRestPoseEnergy.h"
#include "../Codegen/computeReflectionYPlaneDSBendingFixedRestPoseEnergy.h"
#include "../Codegen/computeReflectionXPlaneTheta.h"
#include "../Codegen/computeReflectionYPlaneTheta.h"

reflectionDSHingeFixedRestPoseElement::reflectionDSHingeFixedRestPoseElement(void)
{

}

reflectionDSHingeFixedRestPoseElement::~reflectionDSHingeFixedRestPoseElement(void)
{
}

void reflectionDSHingeFixedRestPoseElement::init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness, bool useReflectionXPlane)
{
	setDofIndices(dofIdcs);
	setParameterIndices(pIdcs);

	m_applyElement = true;

	m_mat = pMat;
	m_thickness = thickness;
	
	m_useReflectionXPlane = useReflectionXPlane;
	m_theta0 = getCurrentAngle(vx);
	m_theta_initial = m_theta0;
	//printf("the rest angle is %f\n", m_theta0);
}


auto computeHingeAngle(const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2, const Eigen::Vector3d& x3)
{
	//vertex numbering : x0, x1 are part of the common edge, x2, x3 are on separate triangles
	//idea : compute theta and its derivative through 1. tan a = sin a / cos a, 2. a = arctan(tan a)
	Eigen::Vector3d x20 = x0 - x2;
	Eigen::Vector3d x21 = x1 - x2;
	Eigen::Vector3d x30 = x0 - x3;
	Eigen::Vector3d x31 = x1 - x3;

	Eigen::Vector3d E = x1 - x0;

	//face normals for the two triangles in deformed and undeformed configuration
	Eigen::Vector3d N2 = x30.cross(x31);
	Eigen::Vector3d N1 = x21.cross(x20);

	double nrm_E = sqrt(E.dot(E));

	double nrm_N1_2 = N1.dot(N1);
	double nrm_N2_2 = N2.dot(N2);

	Eigen::Vector3d En = E / nrm_E;
	Eigen::Vector3d N1n = N1 / sqrt(nrm_N1_2);
	Eigen::Vector3d N2n = N2 / sqrt(nrm_N2_2);

	//compute sin a = through cross(N1, N2) * E_normalized = | N1 | *| N2 | *sin a
	double sin_theta = (N1n.cross(N2n)).dot(En);
	double cos_theta = N1n.dot(N2n);
	//use tangent rule to compute tan(theta / 2) = sin(theta) / (cos(theta) + 1)
	double theta = 2.0 * atan(sin_theta / (1.0 + cos_theta));
	return theta;
}

auto computeReflectionXPlaneTheta(const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2)
{
	double reflectionX = x0[0];
	double x3_x = x2[0] + 2.0 * (reflectionX - x2[0]);
	Eigen::Vector3d x3(x3_x, x2[1], x2[2]);

	return computeHingeAngle(x0, x1, x2, x3);
}

auto computeReflectionXPlaneEnergy (double bendingStiffness, const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2,
	const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2,
	const double& restTheta)
{
	//x3 is the reflection of x2 w.r.t.  plane
	//the plane has two lines x0x1 and z axis

	Eigen::Vector3d E = X1 - X0;
	double E_norm = E.norm();

	double angleDiff = computeReflectionXPlaneTheta(x0, x1, x2) - restTheta;

	double energy = 0.5 * bendingStiffness * angleDiff * angleDiff * E_norm;
	return energy;
}

auto computeReflectionYPlaneTheta(const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2)
{
	double reflectionY = x0[1];
	double x3_y = x2[1] + 2.0 * (reflectionY - x2[1]);
	Eigen::Vector3d x3(x2[0], x3_y, x2[2]);

	return computeHingeAngle(x0, x1, x2, x3);
};
auto computeReflectionYPlaneEnergy(double bendingStiffness, const Eigen::Vector3d& x0, const Eigen::Vector3d& x1, const Eigen::Vector3d& x2,
	const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, const Eigen::Vector3d& X2,
	const double& restTheta)
{
	//x3 is the reflection of x2 w.r.t.  plane
	//the plane has two lines x0x1 and z axis

	Eigen::Vector3d E = X1 - X0;
	double E_norm = E.norm();

	double angleDiff = computeReflectionYPlaneTheta(x0, x1, x2) - restTheta;

	double energy = 0.5 * bendingStiffness * angleDiff * angleDiff * E_norm;
	return energy;
};

double reflectionDSHingeFixedRestPoseElement::computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX)
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

	double k_bend = m_mat->kB;

	double E = 0;

	

	if (m_useReflectionXPlane)
	{
		double testE = computeReflectionXPlaneEnergy(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0);
		Codegen::computeReflectionXPlaneDSBendingFixedRestPoseEnergy(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, E);
	}
	else
	{
		double testE = computeReflectionYPlaneEnergy(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0);
		Codegen::computeReflectionYPlaneDSBendingFixedRestPoseEnergy(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, E);
	}

	return E;
}




void reflectionDSHingeFixedRestPoseElement::addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf)
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

	double k_bend = m_mat->kB;

	
	Eigen::Matrix<double, 9, 1> fx;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, fx);
	else
		Codegen::computeReflectionYPlaneDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, fx);
	fx *= -1.0;

	scatherAdd(fx, vf);
}

void reflectionDSHingeFixedRestPoseElement::computeGradient(const Eigen::VectorXd & vx, const Eigen::VectorXd & vX, Eigen::VectorXd & grad)
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

	double k_bend = m_mat->kB;

	Eigen::Matrix<double, 9, 1> fx;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, fx);
	else
		Codegen::computeReflectionYPlaneDSBendingFixedRestPoseEnergyGradient(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, fx);

	grad = fx;
}

void reflectionDSHingeFixedRestPoseElement::addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A)
{


}

void reflectionDSHingeFixedRestPoseElement::computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes)
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

	double k_bend = m_mat->kB;

	Eigen::Matrix<double, 9, 9> Jx;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, Jx);
	else
		Codegen::computeReflectionYPlaneDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, Jx);

	hes = Jx;
}

void reflectionDSHingeFixedRestPoseElement::addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets)
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

	double k_bend = m_mat->kB;

	double energy = computeEnergy(vx, vX);

	Eigen::Matrix<double, 9, 9> Jx;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, Jx);
	else
		Codegen::computeReflectionYPlaneDSBendingFixedRestPoseEnergyHessian(k_bend, x[0], x[1], x[2],
			X[0], X[1], X[2], m_theta0, Jx);

	scatherAdd(-Jx, triplets);
}



double reflectionDSHingeFixedRestPoseElement::getCurrentAngle(const Eigen::VectorXd& vx)
{
	//compute rest theta
	std::vector<V3D> x(3);
	Eigen::VectorXd x_gather = gatherVector(vx);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
	}
	double theta = 0.0;
	if (m_useReflectionXPlane)
	{
		//double testangle = computeReflectionXPlaneTheta(x[0], x[1], x[2]);
		Codegen::computeReflectionXPlaneTheta(x[0], x[1], x[2], theta);
	}
	else
		Codegen::computeReflectionYPlaneTheta(x[0], x[1], x[2], theta);

	return theta;
}



void reflectionDSHingeFixedRestPoseElement::computedThetadx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& dThetadx)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	Eigen::Matrix<double, 9, 1> element_dThetadx;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneThetaGradient(x[0], x[1], x[2], element_dThetadx);
	else
		Codegen::computeReflectionYPlaneThetaGradient(x[0], x[1], x[2], element_dThetadx);

	dThetadx = element_dThetadx;
}

void reflectionDSHingeFixedRestPoseElement::computed2Thetadx2(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& d2Thetadx2)
{
	//compute theta of this hinge element w.r.t x
	std::vector<V3D> x(3);
	std::vector<V3D> X(3);

	Eigen::VectorXd x_gather = gatherVector(vx);
	Eigen::VectorXd X_gather = gatherVector(vX);
	for (int i = 0; i < 3; i++)
	{
		x[i] = x_gather.segment<3>(i * 3);
		X[i] = X_gather.segment<3>(i * 3);
	}

	Eigen::Matrix<double, 9, 9> element_d2Thetadx2;
	if (m_useReflectionXPlane)
		Codegen::computeReflectionXPlaneThetaHessian(x[0], x[1], x[2], element_d2Thetadx2);
	else
		Codegen::computeReflectionYPlaneThetaHessian(x[0], x[1], x[2], element_d2Thetadx2);

	d2Thetadx2 = element_d2Thetadx2;
}