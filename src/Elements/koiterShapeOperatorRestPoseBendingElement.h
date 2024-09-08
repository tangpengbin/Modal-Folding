#ifndef KOITER_SHAPE_OPERATOR_REST_POST_BENDING_ELEMENT_H
#define KOITER_SHAPE_OPERATOR_REST_POST_BENDING_ELEMENT_H

#if _MSC_VER > 1000
#pragma once
#endif


#include "DofElementBase.h"

class koiterShapeOperatorRestPoseBendingElement : public DofElementBase
{
public:
	koiterShapeOperatorRestPoseBendingElement(void);
	~koiterShapeOperatorRestPoseBendingElement(void);

	void init(const std::vector<int>& dofIdcs, const std::vector<int>& pIdcs, const Eigen::VectorXd& vx, SOMaterial* pMat, double thickness);

	Eigen::Matrix2d computeAngleAndShapeOperator(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX,
		double& theta1, double& theta2, double& theta3);
	double computeEnergy(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	//rest state assumed constant (use precomputed variables)
	void addGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& vf); 
	void computeGradient(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::VectorXd& grad) override;

	void addHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, SpMat & A);

	void computeHessian(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, Eigen::MatrixXd& hes) override;
	void addTriplets(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, TripVec& triplets);

	//this use the x to compute theta
	Eigen::Matrix2d getCurrentShapeOperator(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX);
	void computedShapeOperatordx(const Eigen::VectorXd& vx, const Eigen::VectorXd& vX, std::vector<Eigen::Matrix2d>& dShapeOperatordx);

	int nverts() const { return 4; }
	void setRestShapeOperator(const Eigen::Matrix2d& shapeOperator) { m_shapeOperator0 = shapeOperator; }
	Eigen::Matrix2d getRestShapeOperator() { return m_shapeOperator0; }
	Eigen::Matrix2d getInitialRestShapeOperator() { return m_shapeOperator_initial; }
protected:
	SOMaterial const * m_mat;
	Eigen::Vector3d binary_multiplier;
	//the rest angle
	Eigen::Matrix2d m_shapeOperator0;
	Eigen::Matrix2d m_shapeOperator_initial;
	double m_thickness;
	double m_geomFac;

	bool orthotropicBending;
	//orthotropic parameters
	double Y1;
	double v01;
	double Y0;
	double v10;
	double Y01;
	double G01;
};

#endif