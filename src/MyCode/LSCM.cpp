#include "LSCM.h"



std::pair<Eigen::Vector2d, Eigen::Vector2d >RotateToPlane(acamcad::polymesh::MVert* x_0, acamcad::polymesh::MVert* x_1, acamcad::polymesh::MVert* x_2) {
	// 摆放到以q为原点的坐标系中以去除平移
	auto vec_x_1 = x_1->position() - x_0->position();
	auto vec_x_2 = x_2->position() - x_0->position();

	Eigen::Vector2d v1, v2;
	v1 << vec_x_1.norm(), 0;
	v2<< vec_x_1.dot(vec_x_2) / vec_x_1.norm(), vec_x_1.cross(vec_x_2).norm() / vec_x_1.norm();
	return std::make_pair(v1, v2);
}

void LSCM::LSCMParametrization(acamcad::polymesh::PolyMesh* polyMesh) {
	using namespace acamcad::polymesh;
	std::vector< MPolyFace*> polyfaces = polyMesh->polyfaces();
	std::vector< MVert*> boundVertices = polyMesh->boundaryVertices();
	int qCount = polyMesh->numPolygons();
	int vCount = polyMesh->numVertices();

	// 随便固定两个顶点
	int fv1Id = 0;
	int fv2Id = boundVertices.size() / 2;
	MVert *fv1 = boundVertices[fv1Id];
	MVert *fv2 = boundVertices[fv2Id];


	// 虚数矩阵M转为实数形势处理 
	// M = A+iB, U = u+iv
	// C(T) = U(MM)U = [u^t,v^t] M2^tM2 [u,v] 
	// dC =  M2^tM2 [u,v]  = 0
	// M2 =| A,B,| -B,A|

	Eigen::SparseMatrix<double> M2Sparse(2 * qCount,  2 * qCount);
	Eigen::SparseMatrix<double> BSparse(2 * qCount, 4);
	std::vector<Eigen::Triplet<double>> M2Triplet;
	std::vector<Eigen::Triplet<double>> BTriplet;

	for (int q = 0;q < polyfaces.size();q++) {
		auto face = polyfaces[q];
		//assert(face->index() == q);
		auto v0 = face->halfEdge()->fromVertex();
		auto v1 = face->halfEdge()->next()->fromVertex();
		auto v2 = face->halfEdge()->next()->next()->fromVertex();

		// 转为平面
		auto planeV = RotateToPlane(v0, v1, v2);
		auto v1p = planeV.first;
		auto v2p = planeV.second;
		Eigen::Matrix2d localX;
		localX << v1p, v2p;
		double vol = localX.determinant();
		double sqrtvol = sqrt(vol);
		assert(vol > 0);

		//对边向量w
		Eigen::Vector2d w0 = v2p - v1p;//v2-v1
		Eigen::Vector2d w1 = - v2p;//v0-v2
		Eigen::Vector2d w2 = v1p;//v1-v0

		// Mij = Wij/sqrt(A)
		// M2 =| A,B|
		//     |-B,A|
		if (v0 != fv1 && v0 != fv2) 
		{
			Eigen::Vector2d m = w0 / sqrtvol;
			int vId = v0->index();
			M2Triplet.emplace_back(q, vId, m.x());//A
			M2Triplet.emplace_back(q + qCount, vId, -m.y());//-B
			M2Triplet.emplace_back(q, vId + vCount, m.y());//B
			M2Triplet.emplace_back(q + qCount, vId + vCount, m.x());//A
		}
		if (v1 != fv1 && v1 != fv2) 
		{
			Eigen::Vector2d m = w1 / sqrtvol;
			int vId = v1->index();
			M2Triplet.emplace_back(q, vId, m.x());//A
			M2Triplet.emplace_back(q + qCount, vId, -m.y());//-B
			M2Triplet.emplace_back(q, vId + vCount, m.y());//B
			M2Triplet.emplace_back(q + qCount, vId + vCount, m.x());//A
		}		
		if (v2 != fv1 && v2 != fv2) 
		{
			Eigen::Vector2d m = w2 / sqrtvol;
			int vId = v2->index();
			M2Triplet.emplace_back(q, vId, m.x());//A
			M2Triplet.emplace_back(q + qCount, vId, -m.y());//-B
			M2Triplet.emplace_back(q, vId + vCount, m.y());//B
			M2Triplet.emplace_back(q + qCount, vId + vCount, m.x());//A
		}

		// 固定点
		if (v0 == fv1||v1==fv1||v2 == fv1) {
			Eigen::Vector2d m;
			if(v0==fv1)			m = w0 / sqrtvol;
			else if (v1 == fv1)	m = w1 / sqrtvol;
			else				m = w2 / sqrtvol;

			BTriplet.emplace_back(q, 0, -m.x());
			BTriplet.emplace_back(q + qCount, 2, m.x());
			BTriplet.emplace_back(q, 2, m.y());
			BTriplet.emplace_back(q + qCount, 0, -m.y());

		}
		if (v0 == fv2 || v1 == fv2 || v2 == fv2) {
			Eigen::Vector2d m;
			if (v0 == fv2)		m = w0 / sqrtvol;
			else if (v1 == fv2)	m = w1 / sqrtvol;
			else				m = w2 / sqrtvol;

			BTriplet.emplace_back(q, 1, -m.x());
			BTriplet.emplace_back(q + qCount, 3, m.x());
			BTriplet.emplace_back(q, 3, m.y());
			BTriplet.emplace_back(q + qCount, 1, -m.y());
		}
	}
	M2Sparse.setFromTriplets(M2Triplet.begin(), M2Triplet.end());
	BSparse.setFromTriplets(BTriplet.begin(), BTriplet.end());
	M2Sparse.makeCompressed();
	BSparse.makeCompressed();

	Eigen::Vector4d Fixed;
	Fixed <<	0, 0, 
				1, 0;
	Eigen::MatrixXd B = BSparse*Fixed;

	// solve
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<double> > Solver_sparse;
	Solver_sparse.setTolerance(0.000001f);
	Solver_sparse.compute(M2Sparse);
	Eigen::MatrixXd newVtx = Solver_sparse.solve(B.sparseView());

	//Eigen::SparseLU<Eigen::SparseMatrix<double>>sparseSolver;
	//sparseSolver.analyzePattern(M2Sparse);
	//sparseSolver.factorize(M2Sparse);
	//Eigen::MatrixXd newVtx = sparseSolver.solve(B.sparseView());

	//std::cout << result << std::endl;
	// set vtx
	for (int q = 0;q < polyfaces.size();q++) {
		auto face = polyfaces[q];
		auto v0 = face->halfEdge()->fromVertex();
		auto v1 = face->halfEdge()->next()->fromVertex();
		auto v2 = face->halfEdge()->next()->next()->fromVertex();

		v0->setPosition(newVtx(v0->index()), newVtx(v0->index() + vCount), 0);
		v1->setPosition(newVtx(v1->index()), newVtx(v1->index() + vCount), 0);
		v2->setPosition(newVtx(v2->index()), newVtx(v2->index() + vCount), 0);
	}
	fv1->setPosition(0,-1, 0);
	fv2->setPosition(0,0, 0);
}