#ifndef  SQUASH_TRIANGLE
#define SQUASH_TRIANGLE

#include <src/PolyMesh/include/PolyMesh/PolyMesh.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stack>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>

class SquashTriangle
{
public:
	acamcad::polymesh::PolyMesh* polyMesh;
	acamcad::polymesh::PolyMesh* storedMesh;
	SquashTriangle(acamcad::polymesh::MPolyFace* polyFace);
	~SquashTriangle();
	int x0Id, x1Id, x2Id;
	Eigen::Matrix2d restShape;//D_m 
	Eigen::Matrix2d restShapeInv;
	Eigen::Matrix2d poseShape;// D_s
	Eigen::Matrix2d F;// Deformation Gradients
	Eigen::Matrix2d U,V;// SVD
	Eigen::Vector2d Sigma;
	Eigen::Matrix2d S, R;// Polar decomposition
	double I1, I2, I3;// Smith invarients
	double vol;

	Eigen::Matrix<double, 4, 6> PFPx;// vec(\patial F / \patial x),F∈R^{2*2},x∈R^{2*3}

	void UpdateDeformationGradient(acamcad::polymesh::MPolyFace* polyFace);
	void UpdateDeformationGradient(Eigen::VectorXd& verticesX);

	// ARAP
	double GetARAP();
	Eigen::Vector4d GetARAPGradient();
	Eigen::Matrix4d GetARAPEigenSystem();

	//// SymmetricDirichlet
	double GetSD();
	Eigen::Vector4d GetSDGradient();
	Eigen::Matrix4d GetSDEigenSystem();


	double GetMinStep();

private:
	void TransferToLocalPlane(acamcad::polymesh::MPolyFace* polyFace, Eigen::Matrix2d& x);
	void UpdateInvarients();
	void UpdatePFPx();

};


#endif