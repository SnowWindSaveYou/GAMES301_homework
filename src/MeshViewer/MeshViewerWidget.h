#pragma once
#include <QString>
#include "QGLViewerWidget.h"
#include "../PolyMesh/include/PolyMesh/IOManager.h"
#include "src/MyCode/FixedBoundaryParaAlg.h"
#include "src/MyCode/EigensystemDistortion.h"
#include "src/MyCode/SquashTriangle.h"

class MeshViewerWidget : public QGLViewerWidget
{
	Q_OBJECT
	//Q_ENUM(BoundaryType)
public:
	MeshViewerWidget(QWidget* parent = 0);
	virtual ~MeshViewerWidget(void);
	bool LoadMesh(const std::string & filename);
	void Clear(void);
	void UpdateMesh(void);
	bool SaveMesh(const std::string & filename);
	bool ScreenShot(void);
	void SetDrawBoundingBox(bool b);
	void SetDrawBoundary(bool b);
	void EnableLighting(bool b);
	void EnableDoubleSide(bool b);
	void ResetView(void);
	void ViewCenter(void);
	void CopyRotation(void);
	void LoadRotation(void);
	
signals:
	void LoadMeshOKSignal(bool, QString);
public slots:
	void PrintMeshInfo(void);
	void DoTuttleFlatten(void);
	void DoFloaterFlatten(void);
	void ChangeBoundaryType(int);
	void DoARAPOpt(void);
	void DoSDOpt(void);
	void DoLSCM(void);
	void DoDOF(void);
protected:
	virtual void DrawScene(void) override; 
	void DrawSceneMesh(void);

private:
	void DrawPoints(void) const;
	void DrawWireframe(void) const;
	void DrawHiddenLines(void) const;
	void DrawFlatLines(void) const;
	void DrawFlat(void) const;
	void DrawBoundingBox(void) const;
	void DrawBoundary(void) const;
protected:
	acamcad::polymesh::PolyMesh* polyMesh = new acamcad::polymesh::PolyMesh();

	// 保留个备份用于和初始mesh做计算或显示
	acamcad::polymesh::PolyMesh* storedMesh = new acamcad::polymesh::PolyMesh();

	QString strMeshFileName;
	QString strMeshBaseName;
	QString strMeshPath;
	acamcad::MPoint3 ptMin;
	acamcad::MPoint3 ptMax;
	
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool isDrawBoundingBox;
	bool isDrawBoundary;
	bool isParamerized = false;

	BoundaryType boundaryType = BoundaryType::Circle;

	std::shared_ptr< EigensystemDistortion> eigenAnalysis;
};
