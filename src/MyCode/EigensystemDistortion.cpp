#include "EigensystemDistortion.h"



EigensystemDistortion::EigensystemDistortion(acamcad::polymesh::PolyMesh* polyMesh, acamcad::polymesh::PolyMesh* storedMesh)
{
		this->polyMesh = polyMesh;
		this->storedMesh = storedMesh;
		//this->maxStepCount = 10;
		this->vertexCount= storedMesh->vertices().size();
		this->qCount = storedMesh->polyfaces().size();
		//this->vecSquashTriengles = std::vector<std::shared_ptr< SquashTriangle>>(qCount);
		this->vecSquashTriengles.resize(qCount);
		// 初始化所有顶点的F
		for (int q = 0; q < qCount; q++) {
			auto face = storedMesh->polyfaces()[q];
			this->vecSquashTriengles[q].reset(new SquashTriangle(face));
			//this->vecSquashTriengles.push_back(std::shared_ptr< SquashTriangle>(new SquashTriangle(face)));
		}
}

EigensystemDistortion::~EigensystemDistortion() {

}

void EigensystemDistortion::Init() {
	// 记录旧的顶点
	using namespace acamcad::polymesh;
	std::vector< MVert*> XOverVertices = polyMesh->vertices();
	
	// vec(X) = x =  xy.x1y1...xnyn
	Eigen::VectorXd x(polyMesh->vertices().size() * 2);
	for (int i = 0; i < XOverVertices.size(); i++) {
		int idx = XOverVertices[i]->index();
		auto pos = XOverVertices[i]->position();
		x(idx * 2) = pos.x();
		x(idx * 2 +1) = pos.y();
	}
	this->vecX = x;

}
void EigensystemDistortion::DoForward() {
	std::cout << "EEA start------" << std::endl;
	using namespace acamcad::polymesh;
	this->ProjectNewton();

}
void EigensystemDistortion::CalcDeformationGradient() {
	for (int q = 0; q < qCount; q++) {
		auto tri = this->vecSquashTriengles[q];
		tri->UpdateDeformationGradient(vecX);
	}
}


double EigensystemDistortion::GetEnergy() {
	double energy = 0;

	if (this->energyType == ARAP) {
		for (int q = 0; q < qCount; q++) {
			energy += vecSquashTriengles[q]->vol * vecSquashTriengles[q]->GetARAP();
		}
	}
	else {
		for (int q = 0; q < qCount; q++) {
			energy += vecSquashTriengles[q]->vol * vecSquashTriengles[q]->GetSD();
		}
	}

	return energy;
}
double EigensystemDistortion::GetEnergy(Eigen::VectorXd& vertexVec) {
	double energy = 0;

	if (this->energyType == ARAP) {
		for (int q = 0; q < qCount; q++) {
			vecSquashTriengles[q]->UpdateDeformationGradient(vertexVec);
			energy += vecSquashTriengles[q]->vol* vecSquashTriengles[q]->GetARAP();
		}
	}
	else {
		for (int q = 0; q < qCount; q++) {
			vecSquashTriengles[q]->UpdateDeformationGradient(vertexVec);
			energy += vecSquashTriengles[q]->vol * vecSquashTriengles[q]->GetSD();
		}
	}

	return energy;
}
void EigensystemDistortion::GetEnergyGradient(Eigen::VectorXd& g) {
	g.setZero();
	for (int q = 0; q < qCount; q++) {
		auto tri = this->vecSquashTriengles[q];

		Eigen::Vector4d energyGrad;
		if (this->energyType == ARAP) {
			energyGrad = tri->GetARAPGradient();
		}
		else {
			energyGrad = tri->GetSDGradient();
		}
		

		Eigen::Vector<double, 6> gLocal;
		gLocal = tri->vol* tri->PFPx.transpose()* energyGrad;
		assert(!std::isnan(gLocal.norm()));
		// fill the local to global
		int indexMap[6] = {
			2 * tri->x0Id ,	2 * tri->x0Id + 1,
			2 * tri->x1Id , 2 * tri->x1Id + 1,
			2 * tri->x2Id , 2 * tri->x2Id + 1
		};
		for (int i = 0;i < 6;i++) {
			g[indexMap[i]] += gLocal[i];
		}
	}

}

void EigensystemDistortion::ProjectHessian(Eigen::SparseMatrix<double>& H) {
	Eigen::SparseMatrix<double> PFPx(qCount, 6);
	H.setZero();
	Eigen::Matrix4d Hq;
	std::vector<Eigen::Triplet<double>> HTriplets(6*qCount);

	for (int q = 0; q < qCount; q++) {
		auto tri = this->vecSquashTriengles[q];

		if (energyType == ARAP) {
			Hq = tri->GetARAPEigenSystem();
		}
		else {
			Hq = tri->GetSDEigenSystem();
		}
		
		auto PFPxLocal = tri->PFPx;
		Eigen::Matrix<double, 6, 6> Hi = tri->vol * PFPxLocal.transpose() * Hq * PFPxLocal;
		// fill the local to global
		int indexMap[6] = {
			2 * tri->x0Id ,	2 * tri->x0Id + 1,
			2 * tri->x1Id , 2 * tri->x1Id + 1,
			2 * tri->x2Id , 2 * tri->x2Id + 1
		};
		for (int i = 0;i < 6;i++) {
			for (int j = 0;j < 6;j++) {
				HTriplets.emplace_back(indexMap[i], indexMap[j], Hi(i, j));
			}
		}
	}
	H.setFromTriplets(HTriplets.begin(), HTriplets.end());
	H.makeCompressed();
};


// 来自hw2最佳作业
double EigensystemDistortion::GetMinStep(
	const Eigen::VectorXd& x,
	const Eigen::VectorXd& d
) {
	double minStep = 1.0;
	Eigen::Matrix2d Ds;
	Eigen::Matrix2d D;
	for (int q = 0; q < qCount; q++) {
		auto tri = this->vecSquashTriengles[q];
		 
		Ds << x.segment(tri->x1Id, 2) - x.segment(tri->x0Id, 2),
				x.segment(tri->x2Id, 2) - x.segment(tri->x0Id, 2);

		//Ds = tri->restShape;
		D << d.segment(tri->x1Id, 2) - d.segment(tri->x0Id, 2),
			d.segment(tri->x2Id, 2) - d.segment(tri->x0Id, 2);
		// 发生退化的条件为行列式为0 -》 |D_{s+1}| = |Ds+ aD| = 0
		// 通过分量展开转为一元两次方程，代入求根公式
		double A = D.determinant();
		double B = D.determinant() + Ds.determinant() - (Ds - D).determinant();
		double C = Ds.determinant();

		double t1 = 0.0;
		double t2 = 0.0;
		if (std::abs(A) > 1.0e-10)
		{
			double Delta = B * B - 4 * A * C;
			if (Delta <= 0)	continue;

			double delta = std::sqrt(Delta); // delta >= 0
			if (B >= 0)
			{
				double bd = -B - delta;
				t1 = 2 * C / bd;
				t2 = bd / (2 * A);
			}
			else
			{
				double bd = -B + delta;
				t1 = bd / (2 * A);
				t2 = (2 * C) / bd;
			}

			assert(std::isfinite(t1));
			assert(std::isfinite(t2));

			if (A < 0) std::swap(t1, t2); // make t1 > t2
			minStep = (t1 > 0) ? (std::min(minStep, t2 > 0 ? t2 : t1)) : minStep;
		}
		else
		{
			t1 = -C / B;
			//    avoid divide-by-zero
			minStep = (B == 0) ? minStep : ((t1 > 0) ? t1 : minStep);
		}
	}
	return minStep;
}

// from libigl
double EigensystemDistortion::LinearSearch(
	const Eigen::VectorXd& d,
	const Eigen::VectorXd& g,
	double stepSize,
	double decay
) {
	int MAX_STEP_SIZE_ITER = 12;

	double old_energy = GetEnergy(vecX);
	double new_energy = old_energy;
	int cur_iter = 0; 
	float ddotg = d.dot(g);

	while (new_energy >= old_energy && cur_iter < MAX_STEP_SIZE_ITER)
	{
		Eigen::VectorXd new_x = vecX + stepSize * d;

		double cur_e = GetEnergy(new_x);
		if (cur_e >= old_energy + 1.0e-4 * stepSize * ddotg)//使充分下降 Armijo 
		{
			stepSize *= decay;
		}
		else
		{
			vecX = new_x;
			new_energy = cur_e;
		}
		cur_iter++;
	}
	return new_energy;

};


void EigensystemDistortion::ProjectNewton() {
	Eigen::SparseMatrix<double> H(vertexCount*2, vertexCount*2);// xy per v
	Eigen::VectorXd bi(vertexCount * 2);// G
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> LDLTSolver;
	LDLTSolver.setShift(std::numeric_limits<float>::epsilon());

	double stepSize = 1;
	//bool hasAnalyzePattern = false;
	double energy = -1;
	const double decayRate = 0.8;

	for (int i = 0;i < maxStepCount;i++) {
		
		this->CalcDeformationGradient();
		this->GetEnergyGradient(bi);

		double maxGi = bi.cwiseAbs().maxCoeff();
		std::cout << "max |Gi| = " << maxGi << "\t";
		if (maxGi < 1.0e-4)
		{
			std::cout << "================CONVERGE!!!================\n";
			break;
		}
		this->ProjectHessian(H);
		LDLTSolver.compute(H);
		Eigen::VectorXd di = LDLTSolver.solve(-bi);
		//auto di = -H.cwiseInverse() * bi;

		// LineSearch
		stepSize =std::min( GetMinStep(vecX, di)*100, 1.0);//TODO 不知为何太小
		energy = this->LinearSearch( di, bi, stepSize, 0.8);
		std::cout << "iter energy: " << energy << std::endl;
	}
	
	UpdateMesh();
};

void EigensystemDistortion::UpdateMesh() {
	auto vertices = this->polyMesh->vertices();
	for (int i = 0;i < vertexCount;i++) {
		acamcad::MPoint3 a(vecX(i * 2), vecX(i * 2 + 1), 0);

		vertices[i]->setPosition(a);
	}
}