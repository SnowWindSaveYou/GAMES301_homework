#include "MeshParamWidget.h"

MeshParamWidget::MeshParamWidget(QWidget *parent)
	: QWidget(parent)
{
	CreateTabWidget();
	CreateLayout();
}

MeshParamWidget::~MeshParamWidget()
{
}

void MeshParamWidget::CreateTabWidget(void)
{
	QVBoxLayout* layout = new QVBoxLayout();

	pbPrintInfo = new QPushButton(tr("Print Mesh Info"));
	connect(pbPrintInfo, SIGNAL(clicked()), SIGNAL(PrintInfoSignal()));
	layout->addWidget(pbPrintInfo);

	//TAG 在这添加组件
	tuttleParamBTN = new QPushButton(tr("Tuttle"));
	connect(tuttleParamBTN, SIGNAL(clicked()), SIGNAL(DoTuttleSignal()));
	layout->addWidget(tuttleParamBTN);

	floaterParamBTN = new QPushButton(tr("Floater"));
	connect(floaterParamBTN, SIGNAL(clicked()), SIGNAL(DoFloaterSignal()));
	layout->addWidget(floaterParamBTN);
	

	boundryTypeComboBox = new QComboBox();
	boundryTypeComboBox->addItem("Circle",(int) BoundaryType::Circle);
	boundryTypeComboBox->addItem("Square", (int)BoundaryType::Square);
	boundryTypeComboBox->addItem("Current", (int)BoundaryType::Current);
	boundryTypeComboBox->setCurrentIndex(0);
	connect(boundryTypeComboBox, SIGNAL(currentIndexChanged(int)), SIGNAL(ChangeBoundaryTypeSignal(int)));
	layout->addWidget(boundryTypeComboBox);

	arapOptBTN = new QPushButton(tr("Do ARAP"));
	connect(arapOptBTN, SIGNAL(clicked()), SIGNAL(DoARAPOptSignal()));
	layout->addWidget(arapOptBTN);

	sdOptBTN = new QPushButton(tr("Do SD"));
	connect(sdOptBTN, SIGNAL(clicked()), SIGNAL(DoSDOptSignal()));
	layout->addWidget(sdOptBTN);


	lscmBTN = new QPushButton(tr("Do LSCM"));
	connect(lscmBTN, SIGNAL(clicked()), SIGNAL(DoLSCMSignal()));
	layout->addWidget(lscmBTN);

	dofBTN = new QPushButton(tr("Do One Form"));
	connect(dofBTN, SIGNAL(clicked()), SIGNAL(DoDOFSignal()));
	layout->addWidget(dofBTN);
	
	layout->addStretch();
	wParam = new QWidget();
	wParam->setLayout(layout);
	saParam = new QScrollArea();
	saParam->setFocusPolicy(Qt::NoFocus);
	saParam->setFrameStyle(QFrame::NoFrame);
	saParam->setWidget(wParam);
	saParam->setWidgetResizable(true);
}

void MeshParamWidget::CreateLayout(void)
{
	twParam = new QTabWidget();
	twParam->addTab(saParam, "Tab");
	QGridLayout *layout = new QGridLayout();
	layout->addWidget(twParam, 0, 0, 1, 1);
	this->setLayout(layout);
}
