
#ifndef  LSCM_H
#define LSCM_H
#include <src/PolyMesh/include/PolyMesh/PolyMesh.h>

//#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
//#include <stack>
//#include <cmath>
//#include <vector>
//#include <iostream>
//#include <random>

class LSCM
{
public:
	LSCM();
	~LSCM();


	static void LSCMParametrization(acamcad::polymesh::PolyMesh* polyMesh);

private:

};

#endif