#include "SquashTriangle.h"


Eigen::Vector4d Vec(Eigen::Matrix2d &m) {
	return Eigen::Vector4d(m(0,0), m(1,0), m(0,1), m(1,1));
}

SquashTriangle::SquashTriangle(acamcad::polymesh::MPolyFace* polyFace) {
	auto x_0 = polyFace->halfEdge()->fromVertex();// q
	auto x_1 = polyFace->halfEdge()->next()->fromVertex();
	auto x_2 = polyFace->halfEdge()->next()->next()->fromVertex();
	this->x0Id = x_0->index();
	this->x1Id = x_1->index();
	this->x2Id = x_2->index();

	//this->TransferToLocalPlane(polyFace, restShape);
	// 摆放到以q为原点的坐标系中以去除平移
	auto vec_x_1 = x_1->position() - x_0->position();
	auto vec_x_2 = x_2->position() - x_0->position();

	// 把3d向量转到平面坐标系中
	restShape << vec_x_1.norm(), vec_x_1.dot(vec_x_2) / vec_x_1.norm(),
		0.0, vec_x_1.cross(vec_x_2).norm() / vec_x_1.norm();
	vol = restShape.determinant() / 2.0;
	restShapeInv = restShape.inverse();
	UpdatePFPx();

	
}
SquashTriangle::~SquashTriangle() {

}

void SquashTriangle::UpdateDeformationGradient(acamcad::polymesh::MPolyFace* polyFace) {
	this->TransferToLocalPlane(polyFace, poseShape);
	F = poseShape * restShapeInv;// F = D_s/ D_m^{-1}
	UpdateInvarients();
}

void SquashTriangle::UpdateDeformationGradient(Eigen::VectorXd &verticesX) {

	Eigen::Vector2d x_0(verticesX(this->x0Id * 2), verticesX(this->x0Id * 2 + 1));
	Eigen::Vector2d x_1(verticesX(this->x1Id * 2), verticesX(this->x1Id * 2 + 1));
	Eigen::Vector2d x_2(verticesX(this->x2Id * 2), verticesX(this->x2Id * 2 + 1));
	x_1 -= x_0;
	x_2 -= x_0;

	poseShape << x_1, x_2;

	F = poseShape * restShapeInv;// F = D_s/ D_m^{-1}
	UpdateInvarients();
}


void SquashTriangle::TransferToLocalPlane(acamcad::polymesh::MPolyFace* polyFace, Eigen::Matrix2d& x) {
	auto x_0 = polyFace->halfEdge()->fromVertex();// q
	auto x_1 = polyFace->halfEdge()->next()->fromVertex();
	auto x_2 = polyFace->halfEdge()->next()->next()->fromVertex();

	// 摆放到以q为原点的坐标系中以去除平移
	auto vec_x_1 = x_1->position() - x_0->position();
	auto vec_x_2 = x_2->position() - x_0->position();

	// 把3d向量转到平面坐标系中
	//auto normal = vec_x_1.cross(vec_x_2);
	//double area = normal.norm();
	//normal /= area;
	//auto x_basis = vec_x_1 / vec_x_1.norm();
	//auto y_basis = normal.cross(x_basis);

	//x<< vec_x_1.norm(), vec_x_2.dot(x_basis), 
	//	0.0,			vec_x_2.dot(y_basis);
	x << vec_x_1.norm(), vec_x_1.dot(vec_x_2) / vec_x_1.norm(),
		0.0, vec_x_1.cross(vec_x_2).norm() / vec_x_1.norm();
}
void SquashTriangle::UpdateInvarients() {
	Eigen::JacobiSVD<Eigen::Matrix2d> svd(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
	svd.computeU();
	svd.computeV();
	U = svd.matrixU();
	V = svd.matrixV();
	Sigma = svd.singularValues();


	 //fix reflection in rotation
	Eigen::Matrix2d adjustMat;
	adjustMat.setIdentity();
	adjustMat(1,1) = (U * V.transpose()).determinant()>0?1:-1; // det(R),不判断可能会丢精度的样子
	double detU = U.determinant();
	double detV = V.determinant();

	if(	detU < 0 && detV>0)U *= adjustMat;
	else if (detU > 0 && detV < 0)V *= adjustMat;

	Eigen::Matrix2d sigmaMat = Sigma.asDiagonal() * adjustMat;
	Sigma = sigmaMat.diagonal();

	// polar to extract rotation
	R = U * V.transpose();
	S = V * sigmaMat * V.transpose();

	//I1 \sum{sig} stretches
	I1 = S.trace();
	//I2 \sum{sig^2}distance
	I2 = S.squaredNorm();
	//I3 det(sig) volume change
	I3 = S.determinant();
}




void SquashTriangle::UpdatePFPx() {
	// PFPx = vec(patial F / patial x), F∈R^{2*2},x∈R^{2*3}
	// (patial F / patial x) = (patial D_s/ patial x) * D_m^{-1}
	// 由于(patial D_s/ patial x)为常数，所以PFPx只和rest shape有关

	// sum of cols
	double sum1 = restShapeInv(0, 0) + restShapeInv(1, 0) ;
	double sum2 = restShapeInv(0, 1) + restShapeInv(1, 1) ;
	PFPx.setZero();
	auto x = restShapeInv;

	PFPx << -sum1, 0, x(0, 0), 0, x(1, 0), 0,
			0, -sum1, 0, x(0, 0), 0, x(1, 0),
			-sum2, 0, x(0, 1), 0, x(1, 1), 0,
			0, -sum2, 0, x(0, 1), 0, x(1, 1);
}
/**
* ARAP = ||F-2||^2
* = ||F||^2 - 2tr(S)+||R||^2 
* = I2 - 2I1 + d^2
*/
double SquashTriangle::GetARAP() {// patial phi/ patial f
	return  I2 - 2 * I1 + 3;
}

// PI2/PF - 2I1/PF
Eigen::Vector4d SquashTriangle::GetARAPGradient() {
	Eigen::Vector4d PI1Pf = Vec(R);
	Eigen::Vector4d PI2Pf = 2 * Vec(F);
	//Eigen::Vector4d PI3Pf(F(1,1),-F(1,0),- F(0, 1), F(0,0));//U奇异值反过来V^-1
	Eigen::Vector4d PphiPf = PI2Pf -  2*PI1Pf;
	
	return PphiPf;
}

Eigen::Matrix4d SquashTriangle::GetARAPEigenSystem() {//PatialR/PatialF
	Eigen::Matrix2d twistMat;
	twistMat << 0, -1, 1, 0;
	twistMat *= 1 / sqrt(2.0);
	Eigen::Matrix2d T = U* twistMat* V.transpose();
	const Eigen::Vector4d e = Vec(T);
	double filtered = I1 >= 2.0 ? 2.0 / I1 : 1;

	Eigen::Matrix4d H;
	H.setIdentity();
	H -= filtered * (e * e.transpose());
	H *= 2.0;
	return H;

	// 
	//double lambda1 = 2 - 4 / pow(Sigma(0), Sigma(1));
	//double lambda2 = 2;
	//double lambda3 = 2;
	//double lambda4 = 2;

	//Eigen::Matrix2d flipMat;
	//flipMat << 0,1,1,0;

	//Eigen::Matrix2d twistMat;
	//twistMat << 0, -1, 1, 0;


	//Eigen::Matrix2d D1 = U * Eigen::Vector2d(1,0).asDiagonal() * V.transpose();
	//Eigen::Matrix2d D2 = U * Eigen::Vector2d(0,1).asDiagonal() * V.transpose();
	//Eigen::Matrix2d L =  1 / sqrt(2) * U * flipMat * V.transpose();
	//Eigen::Matrix2d T = 1 / sqrt(2) * U * twistMat * V.transpose();

	//auto e1 = Vec(T);
	//auto e2= Vec(L);
	//auto e3 = Vec(D1);
	//auto e4 = Vec(D2);
	//
	//
	//// patial R/ patial F = patial^2 I1 / patial F^2 = \sum{lambda*vec(1)*vect(q)^T}
	//// Hq
	//Eigen::Matrix4d vecPRPF;
	//vecPRPF << std::max(lambda1, 0.0) * e1 * e1.transpose()
	//	+ std::max(lambda2, 0.0) * e2 * e2.transpose()
	//	+ std::max(lambda3, 0.0) * e3 * e3.transpose()
	//	+ std::max(lambda4, 0.0) * e4 * e4.transpose();
	//return vecPRPF;
}
/**
* SD = (||F||^2+||F^{-1}||^2)/2
* SD_2D = (I2+I2/I3^2)/2
*/
double SquashTriangle::GetSD() {
	double sd = (I2+I2/pow(I3,2))/2.0;
	return sd;
}
Eigen::Vector4d SquashTriangle::GetSDGradient() {
	Eigen::Matrix2d G = U * Sigma.reverse().asDiagonal() * V.transpose();
	Eigen::Vector4d g = (1 + 1.0 / pow(I3 ,2.0))* Vec(F) //PpsiPI2 * 2f 
						- I2 / pow(I3,3.0)		* Vec(G);		// Ppsi/PI3*g
	return g;
}
Eigen::Matrix4d SquashTriangle::GetSDEigenSystem() {
	double lambda1 = 1 + 3 / pow(Sigma(0), 4.0);
	double lambda2 = 1 + 3 / pow(Sigma(1), 4.0);
	double lambda3 = 1 + 1 / pow(I3, 2.0) + I2 / pow(I3, 3);
	double lambda4 = 1 + 1 / pow(I3, 2.0) - I2 / pow(I3, 3);

	Eigen::Matrix2d flipMat;
	flipMat << 0, 1, 1, 0;

	Eigen::Matrix2d twistMat;
	twistMat << 0, -1, 1, 0;


	Eigen::Matrix2d D1 = U * Eigen::Vector2d(1, 0).asDiagonal() * V.transpose();
	Eigen::Matrix2d D2 = U * Eigen::Vector2d(0, 1).asDiagonal() * V.transpose();
	Eigen::Matrix2d L = U * flipMat * V.transpose() / sqrt(2.0);
	Eigen::Matrix2d T = U * twistMat * V.transpose()/ sqrt(2.0);

	auto e1 = Vec(D1);
	auto e2 = Vec(D2);
	auto e3 = Vec(L);
	auto e4 = Vec(T);
	// patial R/ patial F = patial^2 I1 / patial F^2 = \sum{lambda*vec(1)*vect(q)^T}
	// Hq
	Eigen::Matrix4d vecPRPF = std::max(lambda1, 0.0) * e1 * e1.transpose()
		+ std::max(lambda2, 0.0) * e2 * e2.transpose()
		+ std::max(lambda3, 0.0) * e3 * e3.transpose()
		+ std::max(lambda4, 0.0) * e4 * e4.transpose();
	return vecPRPF;
}

double SquashTriangle::GetMinStep() {
	// 计算全局和计算局部是否一致？

	return 1;

}