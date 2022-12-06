#include "FixedBoundaryParaAlg.h"

//#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <stack>
#include <cmath>
#include <vector>
#include <iostream>
#include <random>

float det(Eigen::Vector2f a, Eigen::Vector2f b) {
	return a.x() * b.y() - b.x() * a.y();
}

Eigen::Vector3f getBaryCenter(Eigen::Vector2f a, Eigen::Vector2f b, Eigen::Vector2f c) {

	//float totArea = det(a - b, a - c);
	float pbc = abs(det(b, c)) / 2;
	float pca = abs(det(c, a)) / 2;
	float pab = abs(det(a, b)) / 2;
	return Eigen::Vector3f(pbc, pca, pab) / (pbc + pca + pab);
}

void setBoundaryOnCircle(
	std::vector< acamcad::polymesh::MVert*> &boundVertices, 
	std::vector<Eigen::Triplet<float>> &boundarylist,
	std::unordered_set<int> &visitedVertices) {
	double circleAngleStep = 2 * std::_Pi / (boundVertices.size());
	for (int i = 0; i < boundVertices.size(); i++) {
		int idx = boundVertices[i]->index();
		boundarylist.push_back(Eigen::Triplet<float>(idx, 0, std::cos(i * circleAngleStep)));
		boundarylist.push_back(Eigen::Triplet<float>(idx, 1,  std::sin(i * circleAngleStep)));
		visitedVertices.insert(idx);
	}
}

void setBoundaryOnSquare(
	std::vector< acamcad::polymesh::MVert*>& boundVertices,
	std::vector<Eigen::Triplet<float>>& boundarylist,
	std::unordered_set<int>& visitedVertices) {
	float edgeSize =(float) (boundVertices.size() / 4);
	const float edgeLength = 1;
	for (int i = 0; i < boundVertices.size(); i++) {
		int idx = boundVertices[i]->index();
		float px, py;
		if (i < edgeSize) {
			 px = edgeLength * (i / edgeSize);
			 py = 0;
		}
		else if (i <= edgeSize*2) {
			 px = edgeLength;
			 py = edgeLength * ((i-edgeSize) / edgeSize);
		}
		else if (i <= edgeSize*3) {
			px = edgeLength - edgeLength * ((i - edgeSize * 2) / edgeSize);
			py = edgeLength;
		}
		else {
			px = 0;
			py = edgeLength - edgeLength * ((i - edgeSize * 3) / edgeSize);
		}
		px -= edgeLength / 2;
		py -= edgeLength / 2;
		boundarylist.push_back(Eigen::Triplet<float>(idx, 0,px));
		boundarylist.push_back(Eigen::Triplet<float>(idx, 1,py));
		visitedVertices.insert(idx);
	}
}

void setBoundaryCurrent(std::vector< acamcad::polymesh::MVert*>& boundVertices,
	std::vector<Eigen::Triplet<float>>& boundarylist,
	std::unordered_set<int>& visitedVertices) {
	float edgeSize = (float)(boundVertices.size() / 4);
	const float edgeLength = 1;
	for (int i = 0; i < boundVertices.size(); i++) {
		int idx = boundVertices[i]->index();
		auto pos = boundVertices[i]->position();
		boundarylist.push_back(Eigen::Triplet<float>(idx, 0, pos.x()));
		boundarylist.push_back(Eigen::Triplet<float>(idx, 1, pos.y()));
		visitedVertices.insert(idx);
	}
}

void FixedBoundaryParaAlg::TuttleParametrization(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh, BoundaryType boundaryType) {
	std::cout << "tuttle start------" << std::endl;
	std::cout << "setting bundary" << std::endl;
	using namespace acamcad::polymesh;
	//find first boundary
	MHalfedge* boundHE = nullptr;
	for (const auto& eh : storedMesh->edges()) {
		if (storedMesh->isBoundary(eh)) {
			boundHE = eh->halfEdge();
			if (!boundHE->isBoundary()) {
				boundHE = boundHE->pair();
			}
			break;
		}
	}
	if (boundHE == nullptr) {
		std::cout << "tuttle flatten err: input not a disk homeomorphism mesh" << std::endl;
		return;
	}

	// set mat

	// right mat, 除了确定的边界都是0
	std::vector< MVert*> boundVertices = polyMesh->boundaryVertices();// 使用形变过的边界作输入

	Eigen::SparseMatrix<float> boundarySparse(polyMesh->numVertices(), 2);
	std::vector<Eigen::Triplet<float>> boundarylist;//构建稀松矩阵
   // 边界放在圆上
	std::unordered_set<int> visitedVertices = std::unordered_set<int>(boundVertices.size());

	if (boundaryType == BoundaryType::Circle) {
		setBoundaryOnCircle(boundVertices, boundarylist, visitedVertices);
	}
	else if(boundaryType == BoundaryType::Square){
		setBoundaryOnSquare(boundVertices, boundarylist, visitedVertices);
	}
	else {
		setBoundaryCurrent(boundVertices, boundarylist, visitedVertices);
	}
	
	boundarySparse.setFromTriplets(boundarylist.begin(), boundarylist.end());
	boundarySparse.makeCompressed();

	//std::cout << boundryMatrix << std::endl;
	//// λ矩阵，链接边的权重
	//Eigen::MatrixXd lambdaMatrix= Eigen::MatrixXd::Zero(polyMesh->numVertices(), polyMesh->numVertices());
	std::vector<Eigen::Triplet<float>> lambdalist;//构建稀松矩阵
	Eigen::SparseMatrix<float> lambdaSparse(storedMesh->numVertices(), storedMesh->numVertices());
	for (auto vIter = storedMesh->vertices_begin(); vIter != storedMesh->vertices_end(); vIter++) {
		auto currVtx = (*vIter);
		int idx = currVtx->index();
		//lambdaMatrix(idx, idx)=1;//a_ij = 1, 权和为1
		//非边界点
		auto bbb = visitedVertices.end();
		if (visitedVertices.find(idx) == visitedVertices.end()) {
			int degree = 0;
			auto he = currVtx->halfEdge();
			do {//计算度数d_i
				degree++;
				he = he->rotateNext();
			} while (he != currVtx->halfEdge());

			float weight = -1.0f / float(degree);
			do {//每个链接边权重为1/d_i
				//lambdaMatrix(idx, he->toVertex()->index()) = -1/degree;//-λ_ij
				int toIdx = he->toVertex()->index();
				lambdalist.push_back(Eigen::Triplet<float>(idx, toIdx,weight));
				he = he->rotateNext();
			} while (he != currVtx->halfEdge());
		}
		lambdalist.push_back(Eigen::Triplet<float>(idx, idx, 1));
	}
	lambdaSparse.setFromTriplets(lambdalist.begin(), lambdalist.end());
	lambdaSparse.makeCompressed();

	//求解
	Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > sparseSolver;
	sparseSolver.setTolerance(0.000001f);
	sparseSolver.compute(lambdaSparse);
	Eigen::MatrixXf result = sparseSolver.solve(boundarySparse);

	//std::cout << result << std::endl;
	auto vertices = polyMesh->vertices();
	for (int i = 0; i < polyMesh->numVertices(); i++) {
		vertices[i]->setPosition(result(i, 0), result(i, 1), 0);
	}

	//
	
}



void FixedBoundaryParaAlg::FloaterParametrization(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh, BoundaryType boundaryType) {
	std::cout << "floater start------" << std::endl;
	std::cout << "setting bundary" << std::endl;
	using namespace acamcad::polymesh;
	//find first boundary
	MHalfedge* boundHE = nullptr;
	for (const auto& eh : storedMesh->edges()) {
		if (storedMesh->isBoundary(eh)) {
			boundHE = eh->halfEdge();
			if (!boundHE->isBoundary()) {
				boundHE = boundHE->pair();
			}
			break;
		}
	}
	if (boundHE == nullptr) {
		std::cout << "floater flatten err: input not a disk homeomorphism mesh" << std::endl;
		return;
	}

	// find ordered bound vertices
	std::vector< MVert*> boundVertices = polyMesh->boundaryVertices();// 使用形变过的边界作输入
	// set mat
	// right mat, 除了确定的边界都是0
	Eigen::SparseMatrix<float> boundarySparse(storedMesh->numVertices(), 2);
	std::vector<Eigen::Triplet<float>> boundarylist;//构建稀松矩阵

	std::unordered_set<int> visitedVertices = std::unordered_set<int>(boundVertices.size());
	if (boundaryType == BoundaryType::Circle) {
		setBoundaryOnCircle(boundVertices, boundarylist, visitedVertices);
	}
	else if (boundaryType == BoundaryType::Square) {
		setBoundaryOnSquare(boundVertices, boundarylist, visitedVertices);
	}
	else {
		setBoundaryCurrent(boundVertices, boundarylist, visitedVertices);
	}
	boundarySparse.setFromTriplets(boundarylist.begin(), boundarylist.end());
	boundarySparse.makeCompressed();

	//std::cout << boundryMatrix << std::endl;
	//// λ矩阵，链接边的权重
	std::vector<Eigen::Triplet<float>> lambdalist;//构建稀松矩阵
	Eigen::SparseMatrix<float> lambdaSparse(storedMesh->numVertices(), storedMesh->numVertices());
	std::vector <std::pair< MVert*, Eigen::Vector2f>> flattenRing = std::vector <std::pair< MVert*, Eigen::Vector2f>>(10);
	std::vector <float> angleList = std::vector <float>(10);

	for (auto vIter = storedMesh->vertices_begin(); vIter != storedMesh->vertices_end(); vIter++) {
		auto currVtx = (*vIter);
		auto currVtxPos = currVtx->position();
		int idx = currVtx->index();

		float sumWeight = 1;
		//非边界点
		if (visitedVertices.find(idx) == visitedVertices.end()) {
			int degree = 0;
			float angleSum = 0;
			sumWeight = 0;
			flattenRing.clear();
			angleList.clear();
			auto he = currVtx->halfEdge();

			//计算度数d_i & 顶点周围角的和
			angleList.push_back(0);
			do {
				auto toVtx = he->toVertex();
				auto nextToVtx = he->rotateNext()->toVertex();

				auto a = (toVtx->position() - currVtxPos).normalized();
				auto b = (nextToVtx->position() - currVtxPos).normalized();
				auto angle = acos(a.dot(b));
				assert(angle > 0);
				angleList.push_back(angle);
				degree++;
				angleSum += angle;

				he = he->rotateNext();
			} while (he != currVtx->halfEdge());

			//将局部扁平化
			float currAngleSum = 0;
			for (int j = 0; j < degree; j++) {
				float angle = angleList[j];
				angle = 2 * std::_Pi * angle / angleSum;//角度展到平面
				currAngleSum += angle;
				// 极坐标位置
				auto sideLength = (he->toVertex()->position() - currVtxPos).norm();
				Eigen::Vector2f vec(sideLength * cos(currAngleSum), sideLength * sin(currAngleSum));
				flattenRing.push_back(std::make_pair(he->toVertex(), vec));
				he = he->rotateNext();
			}
			he = currVtx->halfEdge();
			assert(flattenRing.size() == degree);

			// 重心坐标作为权重
			if (degree == 3) {
				auto baryCenter = getBaryCenter(flattenRing[0].second, flattenRing[1].second, flattenRing[2].second);
				lambdalist.push_back(Eigen::Triplet<float>(idx, flattenRing[0].first->index(), -baryCenter.x()));
				lambdalist.push_back(Eigen::Triplet<float>(idx, flattenRing[1].first->index(), -baryCenter.y()));
				lambdalist.push_back(Eigen::Triplet<float>(idx, flattenRing[2].first->index(), -baryCenter.z()));
				sumWeight += baryCenter.sum();
			}
			else {
				//当超过3度时取每个p_l与p直线相交边线段的两顶点构成的重心坐标均值和。
				std::vector<float> barySum = std::vector<float>(flattenRing.size());
				for (int j = 0; j < flattenRing.size(); j++) {

					auto p = flattenRing[j].second;

					for (int k = 0; k < flattenRing.size(); k++) {
						if (k == j || k + 1 == j)continue;
						auto p1 = flattenRing[k].second;
						auto p2 = flattenRing[(k + 1) % degree].second;
						// 通过判两点是否在p与原点构成的直线两侧来判断直线与线段相交
						float c1 = det(p, p1 - p);
						float c2 = det(p, p2 - p);
						if (c1 * c2 <= 0) {// 非同侧
							auto baryCenter = getBaryCenter(p, p1, p2);
							barySum[j] += baryCenter.x();
							barySum[k] += baryCenter.y();
							barySum[(k + 1) % degree] += baryCenter.z();
							break;
						}
					}
				}

				for (int j = 0; j < flattenRing.size(); j++) {
					float m = (1.0f / degree) * barySum[j];
					assert(m > 0);
					sumWeight += m;
					lambdalist.push_back(Eigen::Triplet<float>(idx, flattenRing[j].first->index(), -m));
				}
			}
			 
		}
		lambdalist.push_back(Eigen::Triplet<float>(idx, idx, sumWeight));
	}
	lambdaSparse.setFromTriplets(lambdalist.begin(), lambdalist.end());
	lambdaSparse.makeCompressed();

	//求解
	//Eigen::LeastSquaresConjugateGradient<Eigen::SparseMatrix<float> > sparseSolver;
	//sparseSolver.setTolerance(0.000001f);
	//sparseSolver.compute(lambdaSparse);
	//Eigen::MatrixXf result = sparseSolver.solve(boundarySparse);
	Eigen::SparseLU<Eigen::SparseMatrix<float>> sparseSolver;
	sparseSolver.analyzePattern(lambdaSparse);
	sparseSolver.factorize(lambdaSparse);
	Eigen::MatrixXf result = sparseSolver.solve(boundarySparse);
	//std::cout << result << std::endl;
	auto vertices = polyMesh->vertices();
	for (int i = 0; i < polyMesh->numVertices(); i++) {
		vertices[i]->setPosition(result(i, 0), result(i, 1), 0);
	}
	// calc mips
	//auto faces = polyMesh->polyfaces();
	//for (int i = 0; i < faces.size(); i++) {
	//	auto o = faces[i]->halfEdge()->fromVertex();
	//	auto a = faces[i]->halfEdge()->toVertex();
	//	auto b = faces[i]->halfEdge()->next()->toVertex();

	//	auto no = vertices[o->index()];
	//	auto na = vertices[a->index()];
	//	auto nb = vertices[b->index()];

	//	Eigen::Matrix2f A(1,1,1,1);
	//}
}