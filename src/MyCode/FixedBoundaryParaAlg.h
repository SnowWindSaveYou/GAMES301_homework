
#ifndef  FIXED_BOUNDARY_PARA_ALG
#define FIXED_BOUNDARY_PARA_ALG

#include <src/PolyMesh/include/PolyMesh/PolyMesh.h>
enum BoundaryType
{
	Circle = 0,
	Square = 1,
	Current= 2
};

/**
* Homework 01
* implement of <<Parametrization and smooth approximation of surface triangulations>>
* Tutle's method and Floater's method
*/ 
class FixedBoundaryParaAlg
{
public:
	FixedBoundaryParaAlg();
	~FixedBoundaryParaAlg();


	static void TuttleParametrization(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh, BoundaryType boundaryType);
	static void FloaterParametrization(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh, BoundaryType boundaryType);

private:

};




#endif // ! FIXED_BOUNDARY_PARA_ALG
