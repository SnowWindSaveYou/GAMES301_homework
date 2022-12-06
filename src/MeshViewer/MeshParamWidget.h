#pragma once

#include <QWidget>
#include <QtGui>
#include <QtWidgets>
#include "src/MyCode/FixedBoundaryParaAlg.h"
class MeshParamWidget : public QWidget
{
	Q_OBJECT

public:
	MeshParamWidget(QWidget *parent = 0);
	~MeshParamWidget(void);
private:
	void CreateTabWidget(void);
	void CreateLayout(void);
signals:
	void PrintInfoSignal();
	void DoTuttleSignal();
	void DoFloaterSignal();
	void ChangeBoundaryTypeSignal(int boundaryType);
	void DoARAPOptSignal();
	void DoSDOptSignal();
	void DoLSCMSignal();
	void DoDOFSignal();
private:
	QTabWidget *twParam;
	QWidget *wParam;
	QScrollArea *saParam;
	QPushButton *pbPrintInfo;
	QPushButton *tuttleParamBTN;
	QPushButton* floaterParamBTN;
	QComboBox* boundryTypeComboBox;
	QPushButton* arapOptBTN;
	QPushButton* sdOptBTN;
	QPushButton* lscmBTN;
	QPushButton* dofBTN;

	
};
