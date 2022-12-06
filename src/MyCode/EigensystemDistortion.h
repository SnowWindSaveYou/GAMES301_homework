
#ifndef  EIGENSYSTEM_DISTORTION_ALG
#define EIGENSYSTEM_DISTORTION_ALG
#include <src/PolyMesh/include/PolyMesh/PolyMesh.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stack>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>
#include "SquashTriangle.h"
class EigensystemDistortion
{
public:
	enum EnergyType{
		ARAP = 0,
		SD = 1
	};
	acamcad::polymesh::PolyMesh* polyMesh;
	acamcad::polymesh::PolyMesh* storedMesh;
	Eigen::VectorXd vecX;
	int maxStepCount=64;
	int qCount;
	int vertexCount;
	EnergyType energyType = SD;

	std::vector<std::shared_ptr< SquashTriangle>> vecSquashTriengles;

	EigensystemDistortion(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh);
	~EigensystemDistortion();

	 void DoForward();
	 void Init();
	 void UpdateMesh();

	 double GetEnergy();
	 double GetEnergy(Eigen::VectorXd& vertexVec);
	 

private:
	void CalcDeformationGradient();
	void GetEnergyGradient(Eigen::VectorXd& g);
	double LinearSearch(
		const Eigen::VectorXd& d,
		const Eigen::VectorXd& g,
		double step_size,
		double decay
	);
	void ProjectHessian(Eigen::SparseMatrix<double>& H);
	void ProjectNewton();
	double GetMinStep(
		const Eigen::VectorXd& x,
		const Eigen::VectorXd& d
	);

};




#endif