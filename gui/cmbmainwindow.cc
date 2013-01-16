#include "cmbmainwindow.h"

#include "plotcontrol.h"
#include "parameterdock.h"
#include "modelinspectordock.h"
#include "automaticrunner.h"

#include "plotcontrol.h"
#include "parameterdock.h"
#include "cmbeasywidget.h"
#include "controlpanel.h"
#include "plotwidget.h"
#include "likeliworker.h"
#include "helpwindow.h"
#include "tipdialog.h"
#include "configurationdialog.h"
#include "q4buttongroup.h"

#include "cmbeasy.h"
#include "lensing.h"
#include "gauge.h"
#include "postscriptplot.h"
#include "chainshop.h"
#include "cmbcalc.h"
#include "model.h"
#include "useful.h"
#include "data.h"
#include "miscmath.h"
#include "quintcosmos.h"

#include <qtabwidget.h>
#include <qprinter.h>
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <qspinbox.h>
#include <qcombobox.h>
#include <qradiobutton.h>



#include <QPixmap>
#include <QStatusBar>
#include <QPrinter>
#include <QMessageBox>
#include <QMenu>
#include <QMenuBar>
#include <QDesktopWidget>
#include <QDockWidget>
#include <QToolTip>
#include <QToolBar>
#include <QTimer>
#include <QSettings>
#include <QAction>
#include <QFileDialog>
#include <QWhatsThis>
#include <QProgressBar>
#include <QTextDocument>
#include <QTextCursor>
#include <QTextTable>
#include <QTextCharFormat>
#include <QInputDialog>
#include <QTableWidget>
#include <QTextStream>

#include <QDebug>

#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <algorithm>

#define INITMINMAX(a,b) a=1e20;b=-1e20;
#define SETMINMAX(a,b,c) if (a < b) b=a; if (a > c) c=a; 


CmbMainWindow::CmbMainWindow()  : QMainWindow(), Mathobject(), analyze(new AnalyzeThis()) , TipD(0), interrupt(false), fullPostscript(false),  postscriptFile("") ,  infoModel(-1),  InitializationOk(true) {
  mIgnoreBadErrors = false;
  mGenerating2dPlotTable = false;
  m2dPlotTableX = 0;

  try {

    qDebug() << endl << endl << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl << endl;
    PixDir = cmbeasyDir("/pix");
    noPingPong=false;

    IsLockedIn1d = false;
    IsKeep2dActive = false;

    initCosmos();
    statusBar();

    // get a cmbcalc Object and connect it's integrate() signals with our slots
    cmbcalc = new CmbCalc();  
    if (cmbcalc == 0) throw Bad_Error("No memory left. Error code 4");
    cmbcalc->setGauge(Gauge::speedyDEInvariant);
    //    cmbcalc->setGauge(Gauge::synchronous);

    cw = new CmbEasyWidget(*cosmos,this,"cmbeasy");
    if (cw == 0) throw Bad_Error("No memory left. Error code 5");
    setCentralWidget(cw);

    mParameterDock = new ParameterDockWidget(*cosmos, this);
    mPlotControlDock = new PlotControlDockWidget(this);
    mModelInspectorDock = new ModelInspectorWidget(this);


    mParamDockWidget = new QDockWidget(this);
    mPlotControlDockWidget = new QDockWidget(this);
    mModelInspectorDockWidget = new QDockWidget(this);
    mParamDockWidget->setWindowTitle("Cosmological Parameters");
    mParamDockWidget->setWidget(mParameterDock);
    mPlotControlDockWidget->setWindowTitle("Plotting Functions");
    mPlotControlDockWidget->setWidget(mPlotControlDock);
    mModelInspectorDockWidget->setWindowTitle("Inspector");
    mModelInspectorDockWidget->setWidget(mModelInspectorDock);

    addDockWidget(Qt::LeftDockWidgetArea, mParamDockWidget);
    addDockWidget(Qt::RightDockWidgetArea, mModelInspectorDockWidget);
    addDockWidget(Qt::BottomDockWidgetArea, mPlotControlDockWidget);

    //mModelInspectorDockWidget->resize(mModelInspectorDock->width(), 250);
    //mModelInspectorDock->resize(mModelInspectorDock->width(), 250);
    resize(750, 550);

    // Bundle some plots together
    map<int,bool>* ValidIds_Common = new map<int,bool>();
    cw->cmbplot->validIds = ValidIds_Common;
    cw->structureplot->validIds = ValidIds_Common;
    cw->TensorPlot->validIds = ValidIds_Common;
    cw->PolarEPlot->validIds = ValidIds_Common;
    cw->PolarCPlot->validIds = ValidIds_Common;
    // the 1-Dim likeliplot get's its own id's
    cw->likeliplot1d->validIds  = new map<int,bool>();

    QuintTypeGroup = new Q4ButtonGroup(cw);
    QuintTypeGroup->addButton(mParameterDock->LambdaType, 0);
    QuintTypeGroup->addButton(mParameterDock->LKTType, 1);
    QuintTypeGroup->addButton(mParameterDock->IPLType, 2);

    mPlotControlDock->AutoScaleButton->setIcon(QIcon(QString::fromStdString(
                 ControlPanel::cmbeasyDir("/pix/autoscale.png"))));


    mPlotControlDock->LikeliBox->setMaxCount( 20 );
    mPlotControlDock->LikeliBox_1d->setMaxCount( 20 );

//X     qDebug() << "sizehint:" << mPlotControlDock->sizeHint();
//X     qDebug() << "sizehint:" << mParameterDock->sizeHint();
//X     mPlotControlDockWidget->resize(mPlotControlDock->sizeHint());
//X     mParamDockWidget->resize(mParameterDock->sizeHint());

    // finally show the widget
    cw->show();


    cfg = new ConfigurationDialog(PrtProfiles,this);
    connect(cfg,SIGNAL(saveProfile()),this,SLOT(savePrtProfile()));
    connect(this,SIGNAL(prtProfilesChanged()),cfg,SLOT(updateProfiles()));

    tool = new QToolBar("Cmbeasy Toolbar",this);
    const int metric = style()->pixelMetric( QStyle::PM_SmallIconSize );
    QSize sz( metric, metric );
    tool->setIconSize( sz );

    addToolBar(tool);

    whatsThisAction = QWhatsThis::createAction(this);
    tool->addAction(whatsThisAction);

    initAction();
    initControl();
    initMenu();

    prt = new QPrinter();

    // get a LikeliWorker for Likelihood

    //Likeli = new LikeliWorker(this, &CmbMainWindow::oneLikeli);
    KWorker = new LikeliWorker(this, &CmbMainWindow::spool);
    RasterizeWorker = new LikeliWorker(this,&CmbMainWindow::drawLikeli);
    ReadWMAPWorker = new LikeliWorker(this,&CmbMainWindow::readWMAP);
    DistillWorker = new LikeliWorker(this,&CmbMainWindow::distill);

    startReadWMAP();

    QLayout* dataBoxLayout = mPlotControlDock->DataBox->layout();

    QCheckBox* wmap5Data = new QCheckBox(mPlotControlDock);
    wmap5Data->setObjectName("cmbeasy_data_wmap05");
    wmap5Data->setProperty("dataFileName", cmbeasyDir("/resources/gui_data/wmap_binned_tt_spectrum_5yr_v3.txt"));
    wmap5Data->setProperty("dataType", "xminmax_table" /*or newdat*/);
    wmap5Data->setText("Wmap 5 year data");
    dataBoxLayout->addWidget(wmap5Data);
    connect(wmap5Data, SIGNAL(toggled(bool)), this, SLOT(readData()));

    QCheckBox* acbar08Data = new QCheckBox(mPlotControlDock);
    acbar08Data->setObjectName("cmbeasy_data_acbar08");
    acbar08Data->setProperty("dataFileName", cmbeasyDir("/resources/acbar08/data/acbar2007.newdat"));
    acbar08Data->setProperty("dataType", "newdat");
    acbar08Data->setText("Acbar 08 data");
    dataBoxLayout->addWidget(acbar08Data);
    connect(acbar08Data, SIGNAL(toggled(bool)), this, SLOT(readData()));

    // Get models and Data
    connect(mPlotControlDock->WMAPBox,SIGNAL(toggled(bool)),this,SLOT(readData()));
    connect(mPlotControlDock->BoomBox,SIGNAL(toggled(bool)),this,SLOT(readData()));
    connect(mPlotControlDock->TegBox,SIGNAL(toggled(bool)),this,SLOT(readData()));

    connect(this,SIGNAL(repaintPlot()),cw->cmbplot,SLOT(repaint()));
    connect(this,SIGNAL(repaintPlot()),cw->TensorPlot,SLOT(repaint()));
    connect(this,SIGNAL(repaintPlot()),cw->structureplot,SLOT(repaint()));
    connect(this,SIGNAL(repaintPlot()),cw->likeliplot,SLOT(repaint()));
    connect(this,SIGNAL(repaintPlot()),cw->likeliplot1d,SLOT(repaint()));


    mPlotControlDock->DrawGroup->setEnabled(!mPlotControlDock->Analytic->isChecked());

    cw->structureplot->setLog(true,true);

    readData();

    qDebug() << "finished reading data" << endl;


    cosmos->setOptDistanceLss(0.12);

    // Parameters p(cosmos);
    /*
    qDebug() << "=====" << endl;
    p.printStatus("meins");
    mMin.printStatus("min");
    mMax.printStatus("max");
    // Model*m = findCloseModel(p);
    qDebug() << "+++++" << endl;
    
    
      Sn1aCorrList = new Data("correlation.dat","correlation");
      Sn1aCorr = Sn1aCorrList->fill2d();
      if (Sn1aCorr == 0) {
      QMessageBox::critical(this,"Sn1a correlation data corrupt","The correlation data is not ok", QMessageBox::Ok,QMessageBox::NoButton);
      }
    */
    /*
    SplineWeb *s = new SplineWeb("likelihood",0, 10,10);
  
    for (double x = -3; x <= 3; x+=0.25) {
      for (double y = -3; y <= 3 ;y+=0.25) {
	double sigma = 0.1; 
	double f = exp(- (x*x+y*y)/(2*sigma));
	s->set(x,y, f);
      }
    }
  
    cw->likeliplot->setWeb(s);
    */

 
  // connect Plotwidget
    connect(cw->cmbplot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->structureplot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->likeliplot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->likeliplot,SIGNAL(leftButtonDoubleClicked(QPointF)),this,SLOT(updateModelInspector(QPointF)));
    connect(cw->likeliplot1d,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->TensorPlot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->PolarEPlot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));
    connect(cw->PolarCPlot,SIGNAL(pressedRMB(int, const QPoint&)),this,SLOT(plotPressedRMB(int, const QPoint&)));

    connect(cw->cmbplot,SIGNAL(highlighted(int)),this,SLOT(plotCursorOverCurve(int)));
    connect(cw->structureplot,SIGNAL(highlighted(int)),this,SLOT(plotCursorOverCurve(int)));
    connect(cw->TensorPlot,SIGNAL(highlighted(int)),this,SLOT(plotCursorOverCurve(int)));

    connect(cw->cmbplot,SIGNAL(scaleChanged()),this,SLOT(syncPlot()));
    connect(cw->structureplot,SIGNAL(scaleChanged()),this,SLOT(syncPlot()));
    connect(cw->likeliplot,SIGNAL(scaleChanged()),this,SLOT(syncPlot()));
    connect(cw->TensorPlot,SIGNAL(scaleChanged()),this,SLOT(syncPlot()));

    connect(mPlotControlDock->AntiAlias,SIGNAL(toggled(bool)),cw->cmbplot,SLOT(setAntiAlias(bool)));
    connect(mPlotControlDock->AntiAlias,SIGNAL(toggled(bool)),cw->structureplot,SLOT(setAntiAlias(bool)));
    connect(mPlotControlDock->AntiAlias,SIGNAL(toggled(bool)),cw->TensorPlot,SLOT(setAntiAlias(bool)));

    connect(mPlotControlDock->XLogBox,SIGNAL(toggled(bool)),this,SLOT(resyncPlot()));
    connect(mPlotControlDock->YLogBox,SIGNAL(toggled(bool)),this,SLOT(resyncPlot()));

    connect(cw->PlotTab,SIGNAL(currentChanged ( int )),this,SLOT(syncPlot()));
    connect(mPlotControlDock->PlotXMin,SIGNAL(textChanged(const QString&)), this, SLOT(resyncPlot()));
    connect(mPlotControlDock->PlotXMax,SIGNAL(textChanged(const QString&)), this, SLOT(resyncPlot()));


    connect(mPlotControlDock->AutoScaleButton,SIGNAL(clicked()),this, SLOT(autoScale()));
    connect(mPlotControlDock->LikeliBox_1d, SIGNAL(currentIndexChanged(const QString&)),
                                      this, SLOT(updateParameterNames(const QString&)));
    connect(mPlotControlDock->LikeliBox, SIGNAL(currentIndexChanged(const QString&)),
                                      this, SLOT(updateParameterNames(const QString&)));


    // connect to synchronize the papersize changes
    //   connect(cw->PaperSize_1d,SIGNAL(highlighted(int)),this,SLOT(syncPaperSizeTo1(int)));
    //    connect(cw->PaperSize,SIGNAL(highlighted(int)),this,SLOT(syncPaperSizeTo2(int)));


    // connect info
    connect(mPlotControlDock->AskCosmos,SIGNAL(clicked()),this,SLOT(askCosmos()));

    // connect likelidraw
    connect(mPlotControlDock->LikeliDrawButton,SIGNAL(clicked()),this,SLOT(startDrawLikeli_2d()));
    connect(mPlotControlDock->LikeliAutoButton,SIGNAL(clicked()),this,SLOT(autoAdapt()));
    //connect(cw->Likeli2Dim,SIGNAL(toggled(bool)),this,SLOT(toggledDims(bool)));
    connect(mPlotControlDock->InvokeConfigButton_2d,SIGNAL(clicked()),this,SLOT(advancedPlotSettings()));

    connect(mPlotControlDock->LockIn1d,SIGNAL(clicked()),this,SLOT(lockIn1d()));
    connect(mPlotControlDock->keep2d(), SIGNAL( stateChanged(int) ), SLOT( keep2dStateChanged(int) ) );

    connect(mPlotControlDock->LikeliDrawButton_1d,SIGNAL(clicked()),this,SLOT(startDrawLikeli_1d()));
    connect(mPlotControlDock->LikeliAutoButton_1d,SIGNAL(clicked()),this,SLOT(autoAdapt_1d()));
    connect(mPlotControlDock->InvokeConfigButton_1d,SIGNAL(clicked()),this,SLOT(advancedPlotSettings()));


    // connect quintessence type checkboxes
    connect(QuintTypeGroup,SIGNAL(buttonClicked(int)),this,SLOT(quintTypeChanged(int)));
    connect(mParameterDock->TuneQuint,SIGNAL(toggled(bool)),this,SLOT(tuneQuintChanged()));
    connect(mParameterDock->TuneQuint,SIGNAL(toggled(bool)),mParameterDock->QuintLastS,SLOT(setEnabled(bool)));

    // connect Tensor / Scalar Ratio 
//    connect(cw->Scalar,SIGNAL(toggled(bool)),this,SLOT(scalarOrTensorToggle()));
//    connect(cw->Tensor,SIGNAL(toggled(bool)),this,SLOT(scalarOrTensorToggle()));

//X     resize(800,600);
  
    cmbcalc->initjl(ControlPanel::cmbeasyDir("/resources/jlgen.dat") ,1500);  // max nr of l's
  
  }  catch (Bad_Error x) {
    cerr << "BAD: ERROR" << endl;
    cerr << x.s << endl;
    QMessageBox::critical(this,"Initialization Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
    InitializationOk = false;
  }
  qDebug() << "CMBMAINWINDOW is up and running";
}


/*!
  Creates ControlPanel for use with cmbcalc etc.
*/
void CmbMainWindow::initControl() {
  control = new ControlPanel();
  control->power_cdm = true;
  control->cmb = true;
  control->setInitialConditions(ControlPanel::adiabatic);
}

void CmbMainWindow::initMenu() {
  PlotRMB = new QMenu(this);
  PlotRMB->addAction(plotInfoAction);
  PlotRMB->addAction(plotCutAction);
  PlotRMB->addAction(cosmosAction);
  PlotRMB->addAction(saveCurveAction);

  Plot2ndRMB = new QMenu(this);
  //Plot2ndRMB->addAction(plotPrtGnuplotAction);
  Plot2ndRMB->addAction(printAction);
  Plot2ndRMB->addAction(saveAllCurvesAction);

  PlotLikeliRMB = new QMenu(this);
  //PlotLikeliRMB->addAction(plotPrtGnuplotAction);
  PlotLikeliRMB->addAction(drawLikeliAction);
  PlotLikeliRMB->addAction(flipAxisAction);
  PlotLikeliRMB->addAction(autoAdaptAction);
  PlotLikeliRMB->addAction(PlotLayoutAction);

  PlotLikeliHoverRMB = new QMenu(this);
  PlotLikeliHoverRMB->addAction(toFrontAction);
  PlotLikeliHoverRMB->addAction(toBackAction);
  PlotLikeliHoverRMB->addAction(deleteRegionAction);
  PlotLikeliHoverRMB->addAction(colorRegionAction);
  PlotLikeliHoverRMB->addAction(stickyRegionAction);

  FileMenu = new QMenu("&File");
  FileMenu->addAction(exitAction);

  menuBar()->addMenu(FileMenu);

  PlotMenu = new QMenu("&Plot");
  //  PlotMenu->addAction(plotPrtGnuplotAction);
  PlotMenu->addAction(printAction);
  PlotMenu->addAction(autoScaleAction);
  menuBar()->addMenu(PlotMenu);

  CalculationMenu = new QMenu("&Calculation");
  CalculationMenu->addAction(calcAction);
  CalculationMenu->addAction(stopAction);
  menuBar()->addMenu(CalculationMenu);

  LikeliMenu= new QMenu("&Likelihood");
  LikeliMenu->addAction(QIcon(QPixmap(PixDir+"/fileopen.png")),"load likelihood",this,SLOT(loadLikeli()));
  LikeliMenu->addAction(QIcon(QPixmap(PixDir+"/filesave.png")),"distill likelihood",this,SLOT(startDistill()));

  LikeliMenu->addAction(drawLikeliAction);
  LikeliMenu->addAction(flipAxisAction);

  menuBar()->addMenu(LikeliMenu);

  QMenu *ViewMenu = new QMenu("&View");
  ViewMenu->addAction(mParamDockWidget->toggleViewAction());
  ViewMenu->addAction(mPlotControlDockWidget->toggleViewAction());
  ViewMenu->addAction(mModelInspectorDockWidget->toggleViewAction());

  menuBar()->addMenu(ViewMenu);

  HelpMenu = new QMenu("&Help");
  HelpMenu->addAction("Manual",this,SLOT(help()), Qt::Key_F1);
  HelpMenu->addAction(whatsThisAction);
  HelpMenu->addAction( "About", this , SLOT(about()));

  menuBar()->addMenu(HelpMenu);
}

QString CmbMainWindow::cmbeasyDir(QString addon) {
  return QString(ControlPanel::cmbeasyDir(addon.toStdString()).c_str());
}

void CmbMainWindow::help() {
  QString home(cmbeasyDir("/resources/help/index.html")); 

  HelpWindow *help = new HelpWindow(home, ".", 0, "help viewer");
  help->setWindowTitle("CMBEASY Manual");
  if ( QApplication::desktop()->width() > 400
       && QApplication::desktop()->height() > 500 )
    help->show();
  else
    help->showMaximized();
}

void CmbMainWindow::about() {
  QString home(cmbeasyDir("/resources/help/about.html"));

  HelpWindow *help = new HelpWindow(home, ".", 0, "help viewer");
  help->setWindowTitle("About");
  if ( QApplication::desktop()->width() > 300
       && QApplication::desktop()->height() > 300 ) {
    help->resize(500,400);
    help->show();
  }
  else
    help->showMaximized();
}


/*
void CmbMainWindow::marginalizeHelp() {
   QString home(cmbeasyDir("/resources/help/likelihood.html"));
  
   if (dia->TextBrowser->isVisible()) {
     dia->TextBrowser->hide();
     QSize si = dia->minimumSize();
     dia->resize(si);
     dia->setGeometry(QApplication::desktop()->width()/2-si.width()/2,QApplication::desktop()->height()/2-si.height()/2,si.width(),si.height());
   }
   else { 
     dia->TextBrowser->show();
     //dia->resize(600,600);
     dia->setGeometry(QApplication::desktop()->width ()/2-300,QApplication::desktop()->height()/2-300,600,600);
     dia->TextBrowser->setSource(home);
   }
}
*/

void CmbMainWindow::initCosmos() {
  cosmos = new QuintCosmos();
  if (cosmos == 0) throw Bad_Error("No memory left. Error code 2");
  cosmos->setT_cmb(2.726);
  cosmos->setY_he(0.24);
  cosmos->setNuR(3.04); 
  cosmos->setNuNR(0);

  cosmos->InitialPower[0] = 1.0;
}


/*!
  inits actions for menus and toolbar
*/
void CmbMainWindow::initAction() { 
  calcAction = new QAction(QIcon(QPixmap(PixDir+"/int.png")),"&Start calculation",this);
  calcAction->setShortcut(Qt::AltModifier+Qt::Key_S);
  tool->addAction(calcAction);


  stopAction = new QAction(QIcon(QPixmap(PixDir+"/stop.png")), "Stop &Calculation",this);
  stopAction->setShortcut(Qt::AltModifier+Qt::Key_C);
  stopAction->setEnabled(false);
  tool->addAction(stopAction);

  plotInfoAction= new QAction(QIcon(QPixmap(PixDir+"/info.png")), "&Look up Info",this);
  plotInfoAction->setShortcut(Qt::AltModifier+Qt::Key_L);

  plotCutAction =  new QAction(QIcon(QPixmap(PixDir+"/editcut.png")), "&Cut",this);
  plotCutAction->setShortcut(Qt::ControlModifier+Qt::Key_X);

  printAction = new QAction(QIcon(QPixmap(PixDir+"/fileprint.png")),"Print",this);
  printAction->setShortcut(Qt::ControlModifier+Qt::Key_P);
  tool->addAction(printAction);

  //  plotPrtGnuplotAction =  new QAction(QIconSet(QPixmap(PixDir+"/fileprint.png")),"Print Plot using Gnuplot",Qt::ControlModifier + Qt::Key_T,this);


  cosmosAction = new QAction(QIcon(QPixmap(PixDir+"/wizard.png")),"Invoke Cosmos",this);
  cosmosAction->setShortcut(Qt::AltModifier+Qt::Key_I);
  //tool->addAction(cosmosAction);

  autoScaleAction = new QAction(QIcon(QPixmap(PixDir+"/autoscale.png")),"Auto Scale",this);
  autoScaleAction->setShortcut(Qt::ControlModifier + Qt::Key_A);
  tool->addAction(autoScaleAction);

  saveCurveAction = new QAction(QIcon(QPixmap(PixDir+"/filesave.png")),"Save curve",this);
  saveCurveAction->setShortcut(Qt::ControlModifier + Qt::Key_D);

  saveAllCurvesAction  = new QAction(QIcon(QPixmap(PixDir+"/filesave.png")),"Save all curves",this);
  saveAllCurvesAction->setShortcut(Qt::ControlModifier + Qt::Key_T);

  exitAction = new QAction(QIcon(QPixmap(PixDir+"/exit.png")),"Exit",this);
  exitAction->setShortcut(Qt::ControlModifier+Qt::Key_Q);

  drawLikeliAction = new QAction("draw likelihood",this);
  drawLikeliAction->setShortcut(Qt::ControlModifier+Qt::Key_M);
  drawLikeliAction->setEnabled(true);

  flipAxisAction = new QAction("flip axis",this);
  flipAxisAction->setShortcut(Qt::ControlModifier+Qt::Key_F);
  flipAxisAction->setEnabled(false);

  autoAdaptAction = new QAction("auto adapt ranges",this);
  autoAdaptAction->setShortcut(Qt::ControlModifier+Qt::Key_R);
  autoAdaptAction->setEnabled(true);


  revertZoomAction = new QAction(QIcon(QPixmap(PixDir+"/viewmag-.png")),"revert Zoom",this);
  revertZoomAction->setShortcut(Qt::ControlModifier+Qt::Key_Z);
  tool->addAction(revertZoomAction);

  PlotLayoutAction = new QAction("Plot Layout",this);
  PlotLayoutAction->setShortcut(Qt::ControlModifier+Qt::Key_L);
  PlotLayoutAction->setEnabled(true);

  generateStatisticsFileAction = new QAction( "Generate Statistics File", this );
  generateStatisticsFileAction->setShortcut( Qt::ControlModifier+Qt::Key_W );
  addAction( generateStatisticsFileAction );

  generateStatisticsFileAction->setShortcutContext( Qt::ApplicationShortcut );
  generateStatisticsFileAction->setEnabled( true );
  connect(generateStatisticsFileAction, SIGNAL(triggered()), this, SLOT(generateStatisticsFile()));

  generatePlotTableAction = new QAction( "Generate 1d Plot Table", this );
  generatePlotTableAction->setShortcut( Qt::ControlModifier+Qt::Key_Y );
  addAction( generatePlotTableAction );

  generatePlotTableAction->setShortcutContext( Qt::ApplicationShortcut );
  generatePlotTableAction->setEnabled( true );
  connect(generatePlotTableAction, SIGNAL(triggered()), this, SLOT(generatePlotTable()));

  connect(calcAction,SIGNAL(triggered()),this,SLOT(calc()));
  connect(stopAction,SIGNAL(triggered()),this,SLOT(stop()));
  connect(plotInfoAction,SIGNAL(triggered()),this,SLOT(plotInfo()));
  connect(plotCutAction,SIGNAL(triggered()),this,SLOT(plotCut()));
  connect(printAction,SIGNAL(triggered()),this,SLOT(print()));
  connect(cosmosAction,SIGNAL(triggered()),this,SLOT(askCosmos()));
  connect(autoScaleAction,SIGNAL(triggered()),this,SLOT(autoScale()));
  //  connect(plotPrtGnuplotAction,SIGNAL(triggered()),this,SLOT(printGnuplot()));
  connect(saveCurveAction,SIGNAL(triggered()),this,SLOT(saveCurve()));
  connect(saveAllCurvesAction,SIGNAL(triggered()),this,SLOT(saveAllCurves()));

  toFrontAction = new QAction(QIcon(QPixmap(PixDir+"/up.png")),"Bring to Front",this);
  toFrontAction->setShortcut( Qt::AltModifier+Qt::Key_F);
  toFrontAction->setEnabled(true);

  toBackAction = new QAction(QIcon(QPixmap(PixDir+"/down.png")),"Send to Back",this);
  toBackAction->setShortcut( Qt::AltModifier+Qt::Key_B);
  toBackAction->setEnabled(true);

  deleteRegionAction = new QAction(QIcon(QPixmap(PixDir+"/delete.png")),"Delete Region",this);
  deleteRegionAction->setShortcut( Qt::AltModifier+Qt::Key_D);
  deleteRegionAction->setEnabled(true);

  colorRegionAction = new QAction(QIcon(QPixmap(PixDir+"/color.png")),"Color Region",this);
  colorRegionAction->setShortcut( Qt::AltModifier+Qt::Key_C );
  colorRegionAction->setEnabled(true);

  stickyRegionAction = new QAction(QIcon(QPixmap(PixDir+"/pin4c.png")),"Sticky Region",this);
  stickyRegionAction->setShortcut( Qt::AltModifier+Qt::Key_S );
  stickyRegionAction->setEnabled(true);


  connect(flipAxisAction,SIGNAL(triggered()),this,SLOT(flipAxis()));
  connect(autoAdaptAction,SIGNAL(triggered()),this,SLOT(autoAdapt()));
  connect(revertZoomAction,SIGNAL(triggered()),this,SLOT(revertZoom()));
  
  connect(drawLikeliAction,SIGNAL(triggered()),this,SLOT(startDrawLikeli()));
  connect(PlotLayoutAction,SIGNAL(triggered()),this,SLOT(advancedPlotSettings()));
  
  connect(toFrontAction,SIGNAL(triggered()),this,SLOT(toFront()));
  connect(toBackAction,SIGNAL(triggered()),this,SLOT(toBack()));
  connect(deleteRegionAction,SIGNAL(triggered()),this,SLOT(deleteRegion()));
  connect(colorRegionAction,SIGNAL(triggered()),this,SLOT(colorRegion()));
  connect(stickyRegionAction,SIGNAL(triggered()),this,SLOT(stickyRegion()));	  

}
 

void CmbMainWindow::readData() {
  try {
    // get rid of the previous data, if exists
    for (list<Data*>::iterator i=measured.begin(); i!= measured.end();i++) delete *i;
    measured.clear(); 

    QList<const QCheckBox*> dataList = mPlotControlDock->findChildren<const QCheckBox*>(QRegExp("^cmbeasy_data"));
    foreach (const QCheckBox * dataSet, dataList) {
      if (dataSet->isChecked()) {
        //measured.push_back(new Data(cmbeasyDir("/resources/gui_data/wmap_tt.dat").toLatin1(),"WMAP",Data::Yerror));
        QVariant dataFileNameVariant = dataSet->property("dataFileName");
        if (dataFileNameVariant.isValid()) {
          QString dataFileName = dataFileNameVariant.toString();
          if (dataSet->property("dataType").toString() == "xminmax_table") {
            measured.push_back(new Data(dataFileName.toStdString().c_str(),
                  "WMAP-5yr", Data::AsymXYerror));
          } else if (dataSet->property("dataType").toString() == "newdat") {
            measured.push_back(new Data(dataFileName.toStdString().c_str(),
                  "Acbar 08", Data::Newdat));
          } else {
            qWarning() << "CmbMainWindow::readData() - wrong data type for: " << dataSet->objectName();
          }
        } else {
          qWarning() << "CmbMainWindow::readData() - no data for: " << dataSet->objectName();
        }
      }
    }

    // the old, manual way of doing things:
    if (mPlotControlDock->WMAPBox->isChecked())
      measured.push_back(new Data(cmbeasyDir("/resources/gui_data/wmap_tt.dat").toLatin1(),"WMAP",Data::Yerror));
    if (mPlotControlDock->BoomBox->isChecked())
      measured.push_back(new Data(cmbeasyDir("/resources/gui_data/boomerang.dat").toLatin1(),"boomerang '01"));
    if (mPlotControlDock->TegBox->isChecked())
      measured.push_back(new Data(cmbeasyDir("/resources/gui_data/tegmark.dat").toLatin1(),"Tegmark '01"));

    qDebug() << "read in" << endl;

    cw->cmbplot->setDataList(&measured);
    qDebug() << "set data list" << endl;
    // cw->cmbplot->drawData();
    cw->cmbplot->autoScale(false,false);
    qDebug() << "autoscaled" << endl;
    emit(repaintPlot());
    qDebug() << "emitted repaint (1)" << endl;
  } catch (Bad_Error e) {
    QMessageBox::warning(this, "Error Loading Data File", QString::fromStdString(e.s));
  }
}


double CmbMainWindow::sigma8(Model *m) {
  Anchor anchor;
  return sigma8(m->power_cdm(&anchor), m->h);
}

double CmbMainWindow::sigma8(Spline* s,double h) {
  tmpSig8 = s;
  tmpH = h;
  if (!s->isArmed()) s->arm();
  // s->dump("sigma8.dat");
  return sqrt(Miscmath::rombint((moSingle)&CmbMainWindow::sig8Integrator, this, s->start()*h, s->stop()*h,1e-6));
}

double CmbMainWindow::sig8Integrator(const double k) const  {  
  static double tpi =1.0/(2*M_PI*M_PI); //   4*3.14159    
  double r=8/tmpH;
  double x = k*r;
  double win = ( sin(x) - x*cos(x) ) / (x*x*x);
  
  double power = tmpSig8->fastY(k/tmpH) * pow(tmpH,-3);  // -4
  return tpi*power*pow(k*3*win,2);
}



void CmbMainWindow::calc()
{
  interrupt=false;
  calcAction->setEnabled(false);
  stopAction->setEnabled(true);

  mParameterDock->fillCosmos();

  control->power_cdm = mParameterDock->PowerBox->isChecked();

  if (control->power_cdm) { 
    control->highPrecisionTransfer = true;  // if at all cdm, high precision?
    control->transferMaxK=5*cosmos->h();  // maximal k for cdm  
    control->transferPerLog=5;  // k-values per log interval
  }


  cosmos->reset(); //! doesn't change physical things, just resets splines etc
  int n = cosmos->InitialPower.size();
  cl.clear();
  cl.resize(n);
  syncControl();

  mParameterDock->ParameterTab->setCurrentIndex(1);
  bool haveDE = (cosmos->omega_q()!=0);
  bool syncGauge = mParameterDock->SyncGauge->isChecked();
  bool invariantGauge = mParameterDock->InvariantGauge->isChecked();
  if (syncGauge && invariantGauge) {  // can't happen
    QMessageBox::warning(this, "Two Gauges Selected", "... will use invariant gauge.",
                                QMessageBox::Ok);
    syncGauge = false;
  }

  if (invariantGauge && !haveDE) {
    cmbcalc->setGauge(Gauge::speedyInvariant);
  } else if (invariantGauge && haveDE) {
    cmbcalc->setGauge(Gauge::speedyDEInvariant);
  } else if (syncGauge && !haveDE) {
    cmbcalc->setGauge(Gauge::synchronous);
  } else if (syncGauge && haveDE) {
    cmbcalc->setGauge(Gauge::quintSynchronous);
  }


  if (control->tensor) {
    if (cosmos->InitialTensorPower[0] == 0) {
      QMessageBox::warning(this,"Initial Tensor Spectrum 0",
          "You request a tensor spectrum, however,\n"
          "you've set the tensor spectral index  to 0."
          " This will lead to a vanishing tensor spectrum.",
          QMessageBox::Ok, QMessageBox::Abort);
      mParameterDock->ParameterTab->setCurrentIndex(1);
      calcAction->setEnabled(true);
      stopAction->setEnabled(false);
      return;
    }
    if (control->power_cdm && !control->scalar) { 
      QMessageBox::warning(this,"Power-spectrum without scalars but with tensors",
                                "At the moment, you can't have the combination: scalars no, power yes, tensor yes. Sorry. Disable power or enable scalar.",
                                QMessageBox::Ok, QMessageBox::Abort);
      calcAction->setEnabled(true);
      stopAction->setEnabled(false);
      return;
    }
  }

  statusBar()->showMessage("Preparing Cosmos");
  TimeStart = QTime::currentTime();
  try {
    if (mParameterDock->TuneQuint->isChecked())
      cosmos->tuneQuintessence(mParameterDock->QuintLastS->text().toDouble(),0.7);

    allocateWatchThread(KWorker);
  } catch (Bad_Error x) {
    cerr << "BAD: ERROR" << endl;
    cerr << x.s << endl;
    QMessageBox::critical(this,"Initialization Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
    calcAction->setEnabled(true);
    stopAction->setEnabled(false);
  }
  qDebug() << "~CALC" << endl;
}


void CmbMainWindow::allocateWatchThread(LikeliWorker* worker) {
  Active[worker] = false; // to prevent very unlikelik clashes while not being properly initialized
  if (Progress.find(worker)==Progress.end())  Progress[worker] = new QProgressBar(statusBar());
  
  Timer[worker] = new QTimer( this );
  
  connect(Timer[worker], SIGNAL(timeout()), this,SLOT(watchdogCentral()) );
  Step[worker]=0;
  TotalSteps[worker]=100;

  OldStep[worker]=100;
  OldTotalSteps[worker]=100;
  
  StatusMsg[worker] = "";
  OldMsg[worker]="";
  CriticalMsg[worker]="";
  CriticalHead[worker]=""; 
  


  Timer[worker]->start(300);
  statusBar()->addPermanentWidget(Progress[worker],0);
  
  Progress[worker]->reset();
  Progress[worker]->show();
  Progress[worker]->setMaximum(100);
  Progress[worker]->setValue(0);

  
  /*
  if (WatchThread==0) WatchThread = worker; 
  else WatchThread2 = worker;
  */
  Active[worker] = true; // switch on supervision
  worker->ende=true;  // single shot, always
  worker->start();
}



void CmbMainWindow::watchdogCentral() {
  for (map<LikeliWorker*,QProgressBar*>::iterator i = Progress.begin(); i != Progress.end(); i++) {
    LikeliWorker* worker = i->first;
    QProgressBar* progress = i->second;
    if (Active[worker]) {

      if (StatusMsg[worker] != OldMsg[worker]) statusBar()->showMessage(StatusMsg[worker]);
      if (OldStep[worker] != Step[worker]) progress->setValue(Step[worker]);
      if (OldTotalSteps[worker] != TotalSteps[worker]) progress->setMaximum(TotalSteps[worker]);

      OldMsg[worker] = "";
      OldStep[worker] = Step[worker];
      OldTotalSteps[worker] = TotalSteps[worker];

      if (CriticalMsg[worker] != "") {
        cerr << "watchdog, caught critical message: " << CriticalMsg[worker].toStdString() << endl;
        delete Timer[worker];
        progress->hide();
        statusBar()->removeWidget(progress);
        QMessageBox::critical(this,CriticalHead[worker],CriticalMsg[worker],  QMessageBox::Ok,QMessageBox::NoButton);
        worker->wait();
        CriticalMsg[worker]="";
        if (worker== RasterizeWorker) endDrawLikeli();
        Active[worker] = false;
      } else if (worker->isFinished()) {
        qDebug() << "worker finished" << endl;
        delete Timer[worker];
        progress->hide();
        statusBar()->removeWidget(progress);
        Active[worker] = false;
        if (! interrupt) {
          if (worker == KWorker) endSpool();
          if (worker == RasterizeWorker) endDrawLikeli();
          if (worker == ReadWMAPWorker) endReadWMAP();
        }
      }
    }
  }
}

void CmbMainWindow::spool() {
  try {
    int nk = cmbcalc->go(cosmos,"dudei",*control,true);
    qDebug() << "TOTAL SOURCES: " << nk << endl;
    StatusMsg[KWorker] = "Integrating CMB Sources in k-space";
    TotalSteps[KWorker] = nk;
    Step[KWorker]=0;
    for (int k = 1; k <= nk; ++k) {
      if (interrupt) return;
      cmbcalc->oneK(k,cosmos,*control);
      Step[KWorker] = k;
    }
    if (control->power_cdm) {
      StatusMsg[KWorker] = "Integrating remaining modes for Transfer functions";
      int nkt= cmbcalc->goTransfer(cosmos,"dudei",*control,true);
      TotalSteps[KWorker] = nkt-nk-1;
      Step[KWorker]=0;
      for (int k = nk + 1; k <= nkt; ++k) {
	if (interrupt) return;    
	cmbcalc->oneKTransfer(k, cosmos,*control);
	Step[KWorker] = k - nk -1;
      }  
    }
    if (interrupt) return;
    
    StatusMsg[KWorker] = "Interpolating sources";
    int maxJ = cmbcalc->prepareCl(cosmos,*control,cl,true);
    if (control->scalar) {
      StatusMsg[KWorker] = "Integrating Scalar Cl's";
      TotalSteps[KWorker] = maxJ;
      for (int j = 0; j < maxJ ; ++j) {
	if (interrupt) return;
	Step[KWorker] = j;
	cmbcalc->oneCl(cosmos,j,cl,*control);
      }
      cmbcalc->ridOfScalarIntegrator();
    }
    if (control->tensor) {
      StatusMsg[KWorker] = "Integrating Tensor Cl's";
      TotalSteps[KWorker] = maxJ;
      for (int j = 0; j < maxJ ; ++j) {
	if (interrupt) return;
	Step[KWorker] = j;
	cmbcalc->oneTensorCl(cosmos,j,cl,*control);
      }
      cmbcalc->ridOfTensorIntegrator();
    }
  } catch  (Bad_Error x) {
    cerr << "BAD: ERROR" << endl;
    cerr << x.s << endl;
    CriticalHead[KWorker] = "cmb calculation error";
    CriticalMsg[KWorker] = x.s.c_str();
    interrupt=true;
  }

}

void CmbMainWindow::endSpool() {
  try {
    cmbcalc->finish(cosmos,*control,cl);
    
    calcAction->setEnabled(true);
    stopAction->setEnabled(false);
 
    if (control->cmb) cl.ts[0]->setChildrensN();
    
    if (control->cmb) {
      // cl.tt[0]->dump("tensor");
      vector<WMAPNorm> norm;
      // for the very first run, the user may be faster to press "calc" 
      // than the conversion to binary of the wmap data can be accomplished
      // hence, we only make a quick WMAP normalize in that case  (i.e. without chi^2)
      analyze->scaleCls(cl,0,pow(2.725e6,2));  // bring the cl's to muK^2 
      // New in cmbeasy v4.0, you can set the initial amplitudes for A_s and A_t
      // or use the inflationary result 
      // and apply or not apply the wmap normalization
      vector<double> A_s, A_t;  // the Amplitudes. For each spectral index an Amplitude A_s and A_t
      analyze->fiducialAmplitudes(*cosmos,A_s,A_t); // Initialize the vectors A_s and A_t (convenience)
      // For our single spectral index, you can choose A_s & A_t (not needed if WMAP normalized later)
      A_s[0] = 20e-10; A_t[0] = 20e-10; 
      // Apply the A_t  = -8 n_t A_s inflationary formula (if you want this, if not, comment out)
      analyze->applyInflationaryTensorRatio(*cosmos,A_s,A_t); 
      // In all cases, you need to call rescaleSpectra() which applies the necessary factors of pi etc
      // to the spectra and also calculates sigma_8 [upon wmapnormalizing, sigma_8 will be recalcualted
      // automatically]
      analyze->rescaleSpectra(*cosmos,*control,cl,A_s,A_t);
      // finally WMAP normalize (if you like to, if not, comment out);
     
     
      if (analyze->WMAPNotYetInitialized) analyze->quickWMAPNormalize(*cosmos,*control,cl);
      else  norm = analyze->WMAPNormalize(*cosmos,*control,cl); // WMAP normalize s

      cl.ts[0]->arm(Spline::all);      // arm all output splines
      cl.tt[0]->dump("aftertensor");
      //Model neu(*cosmos,cl,sigma8[0],*control);
      if (control->tensor) { 
	cw->TensorPlot->setSpline(*(cl.tt[0])); 
	cw->TensorPlot->autoScale();
      }
      else cw->TensorPlot->incSplineId();

      if ( control->isLensing() )
      {
        Lensing lens( *cosmos,*control,cl,*cmbcalc );
        lens.doLensing( Lensing::AllSky );
      }

      cl.ts[0]->disarm(Spline::all);


      cl.ts[0]->arm(Spline::all);
      int id = cw->cmbplot->setSpline(*(cl.ts[0]));
      plotModel[id] =  new Model(*cosmos,cl,*control,0,1500); // Model //  neu;

      ExtendedInfo info;
      fillExtendedInfo(info); 
      
      info.wmap_tt = norm[0].chi2_tt;	
      info.wmap_te = norm[0].chi2_te;
      extendedInfo[id] = info;
      

      qDebug() << "analyze says: MAP_TT: " << norm[0].chi2_tt << endl;	

      
      updateInfo(id);
      cw->cmbplot->autoScale();
      cw->PolarEPlot->setSpline(*(cl.es[0]));
      cw->PolarCPlot->setSpline(*(cl.cs[0]));
      cw->PolarEPlot->autoScale(true,false);
      cw->PolarCPlot->autoScale(true,false);

      cw->cmbplot->printStatus();
   
      if (control->power_cdm) {
	cosmos->dumpPower(0,"cdm");
	Spline *s =  cosmos->createPower(0,"cmbmainwindow_cdm");
	s->arm();
	//s->dump("save",false);
	cw->structureplot->setSpline(*s); 
	delete s;
	cw->structureplot->autoScale(true,false);
	qDebug() << cw->structureplot->Co.x << endl;
	qDebug() << cw->structureplot->Co.X << endl;
	qDebug() << cw->structureplot->Co.y << endl;
	qDebug() << cw->structureplot->Co.Y << endl;
      } else cw->structureplot->incSplineId();  // skip this id in structure plot for keeping all id's in sync


      //cw->plot->drawSplines();
      emit(repaintPlot());
    } else {
      cw->cmbplot->incSplineId();
      cw->PolarEPlot->incSplineId();
      cw->PolarCPlot->incSplineId();
    }

    if (control->power_cdm) {
      string base = ControlPanel::cmbeasyDir("/output/");
      cosmos->dumpPower(0,base+"cdm");
      Spline *s =  cosmos->createPower(0,"cmbmainwindow_cdm");
      cw->structureplot->setSpline(*s); 
      delete s;
      cw->structureplot->autoScale(true,false);   
    } else cw->structureplot->incSplineId();  // skip this id in structure plot for keeping all id's in sync
   
       
    statusBar()->showMessage("finished in " +  toStr(TimeStart.secsTo(QTime::currentTime())) + " s",2000);   
    emit(repaintPlot());
    cosmos->printStatus();

  }  catch (Bad_Error x) {
      qWarning() << "endspool bad error" << endl;
      qWarning() << QString::fromStdString(x.s) << endl;
      QMessageBox::critical(this,"Bad Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
  }   
}


void CmbMainWindow::saveCurve() {
  if (RMBPlot_id) {
    PlotWidget *p = currentPlotWidget();
    if (p) { 
      QString saveName;
      saveName = QFileDialog::getSaveFileName(this, "Save Likelihood file as...", PlotDir, "Plot data (*.plt)");
      if (saveName == QString::null) return;
      if (saveName.right(4) == ".plt") saveName = saveName.left(saveName.length()-4);
      PlotDir = getDir(saveName);
      ofstream o((saveName+".plt").toLatin1());
      p->saveToFile(RMBPlot_id,o);
    }
  }
}



void CmbMainWindow::saveAllCurves() {
  PlotWidget *p = currentPlotWidget();
  if (p) { 
    QString saveName;
    saveName = QFileDialog::getSaveFileName(this, "Save Plot as...", PlotDir,"Plot data (*.plt)");
    if (saveName == QString::null) return;
    if (saveName.right(4) == ".plt") saveName = saveName.left(saveName.length()-4);
    PlotDir = getDir(saveName);
    ofstream o((saveName+".plt").toLatin1());
    p->saveToFile(o);
  }
}

/*!
  Slot connected to PlotWidget RMB press - detection, QPoint pos not needed right
  now 
*/
void CmbMainWindow::plotPressedRMB(int id, const QPoint& pos) {
  qDebug() << "CMBMAINWINDOW: plotPressed" << endl;
  RMBPlot_currentPage = cw->PlotTab->currentIndex();
   if (RMBPlot_currentPage== 3) { // i.e. likelihood plot
     if (id != -1) {
       RMBPlot_id = id;
       PlotLikeliHoverRMB->popup(QCursor::pos());
     }
     else PlotLikeliRMB->popup(QCursor::pos());
     
   } else {
     if (id) {  // so if it has been over a spline
       RMBPlot_id = id;
       cosmosAction->setEnabled(extendedInfo.find(id) == extendedInfo.end() && RMBPlot_currentPage != 4);
       PlotRMB->popup(QCursor::pos());
     } else {
       Plot2ndRMB->popup(QCursor::pos());
       // to do
     }
   }
}

void CmbMainWindow::plotCursorOverCurve(int id) {
  updateInfo(id);
  RMBPlot_id = id;  // this is sensible, cause we want this for Keyboard
  // shortcuts also (which invoke actions and need an id) 
}


void  CmbMainWindow::plotInfo() {
  if (RMBPlot_id) {
    updateInfo(RMBPlot_id);
    mPlotControlDock->BigTab->setCurrentIndex(0);
  }
}

void  CmbMainWindow::updateInfo(int id) {
  qDebug() << "updateinfo: " << id << endl;
  if (plotModel.find(id) == plotModel.end())  return;
  infoModel = id;
  Model* m = plotModel[id];
  qDebug() << "m: " << m << endl;
 

  QString t;
  bool quintessence = false;
  
  ExtendedInfo e;


  if (extendedInfo.find(id) != extendedInfo.end()) {
    e = extendedInfo[id];
    //qDebug() << "THERE IS AN EXTENDED INFO" << endl;
    //qDebug() << e.QuintType << "   " << e.o_q << endl;
    if (e.QuintType != "Base") quintessence = true;
  }

  //qDebug() << "#################################################" << endl;
  // qDebug() << "INFO" << endl;
/*  Qt 4.0.0 has problems with nested tables in HTML; 4.0.1 is supposedly handling them ok though :-/
  t+= "<html><body><table><tr><td>";

  t += "<table>";
  t += "<tr><td>hubble h:</td><td><font color=#aa0000><b> " + toStr(m->h,4) + "</b></font></td></tr>\n";
  t += "<tr><td>Omega_b:</td><td><font color=#aa0000><b> " + toStr(m->o_b,4) + "</b></font></td></tr>\n";
  t += "<tr><td>Omega_cdm:</td><td><font color=#aa0000><b> " + toStr(m->o_cdm,4) + "</b></font></td></tr>\n";
  qDebug() << "deb1" << endl;
  if (quintessence) 
    t += "<tr><td>Omega_q (" +e.QuintType +"):</td><td><font color=#00aa00><b> " + toStr(e.o_q ,4) + "</b></font></td></tr>\n";
  else 
    qDebug() << "deb2" << endl;
    t += "<tr><td>Omega_vac:</td><td><font color=#aa0000><b> " + toStr(m->o_lambda,4) + "</b></font></td></tr>\n";
  t += "<tr><td>spectral index:</td><td><font color=#aa0000><b> " + toStr(m->n,4) + "</b></font></td></tr>\n";
  t += "<tr><td>optdlss:</td><td><font color=#aa0000><b> " + toStr(m->optdlss,4) + "</b></font></td></tr>\n";
  t += "<tr><td>sigma8:</td><td><font color=#aa0000><b> " + toStr(m->sigma8,4) + "</b></font></td></tr>\n";
  // qDebug() << "taking on the sigma8" << endl;
  "<tr><td>sigma8 O_m^gamma (0.5 +- 0.1):</td><td><font color=#aa0000><b> " + toStr(cosmos->sigma8Omega(m->n,m->sigma8,QuintCosmos::logAWeff),4) + "</b></font></td></tr>\n";
  
  //qDebug() << "done" <<endl;

  //  t += "<tr><td>Omega_vac:</td><td><font color=#aa0000><b> " + toStr(m->o_lambda,4) + "</b></font></td></tr>\n";

  t += "</table>";


  if (extendedInfo.find(id) != extendedInfo.end()) {
    ExtendedInfo &e = extendedInfo[id];
    t+="</td><td>";
    t += "<table>";
    t+= "<tr><td>tau<sub>0</sub>:</td><td><font color=#0000aa><b> " + toStr(e.tau0,4) + "</b></font></td></tr>\n";
    t+= "<tr><td>tau<sub>ls</sub>:</td><td><font color=#0000aa><b> " + toStr(e.tauls,4) + "</b></font></td></tr>\n";
    t+= "<tr><td>tau<sub>equ</sub>:</td><td><font color=#0000aa><b> " + toStr(e.taueq,4) + "</b></font></td></tr>\n";
    t+= "<tr><td>z<sub>ls</sub>:</td><td><font color=#0000aa><b> " + toStr(e.zls,4) + "</b></font></td></tr>\n";
    t+= "<tr><td>z<sub>equ</sub>:</td><td><font color=#0000aa><b> " + toStr(e.zeq,4) + "</b></font></td></tr>\n";


    t += "</table>";

    t+="</td><td>";
    t += "<table>";
    t+="<tr><td>WMAP chi<sup>2</sup>: </td><td><font color=#008833><b> " + toStr(e.wmap_tt + e.wmap_te,4) + "</b></font></td></tr>\n";
    t+="<tr><td>WMAP (TT ) chi<sup>2</sup>: </td><td><font color=#008833><b> " + toStr(e.wmap_tt,4) + "</b></font></td></tr>\n";
    t+="<tr><td>WMAP (TE) chi<sup>2</sup>: </td><td><font color=#008833><b> " +  toStr(e.wmap_te,4) + "</b></font></td></tr>\n";
    t += "</table>";
  }

  t += "</td></tr></table></body></html>";
*/

  QColor c;
  c.setNamedColor( "#aa0000");
  QTextCharFormat redFormat;
  redFormat.setForeground( c );
  redFormat.setFontWeight( QFont::Bold );
  QTextCharFormat blueFormat;
  c.setNamedColor( "#0000aa");
  blueFormat.setForeground( c );
  blueFormat.setFontWeight( QFont::Bold );
  QTextCharFormat greenFormat;
  c.setNamedColor( "#008833");
  greenFormat.setForeground( c );
  greenFormat.setFontWeight( QFont::Bold );
  QTextCharFormat subScriptFormat;
  subScriptFormat.setVerticalAlignment( QTextCharFormat::AlignSubScript );
  QTextCharFormat superScriptFormat;
  superScriptFormat.setVerticalAlignment( QTextCharFormat::AlignSuperScript );

  QTextDocument *doc = mPlotControlDock->InfoView->document();
  doc->clear();
  QTextCursor cur(doc);
  QTextTable *table = cur.insertTable( 8, 8 );
  QTextTableFormat tableFormat;
  tableFormat.setBorder( 0 );
  table->setFormat( tableFormat );

  int i = 0;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "hubble h:" );
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->h, 4 ), redFormat );
  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "Omega_b:" );
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->o_b, 4 ), redFormat );
  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "Omega_cdm:" );
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->o_cdm, 4 ), redFormat );

  if ( quintessence )
  {
    ++i;
    table->cellAt( i, 0 ).firstCursorPosition().insertText( "Omega_q (" + e.QuintType + "):" + toStr( e.o_q, 4 ) );
    table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->o_cdm, 4 ), redFormat );
  }

  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "Omega_vac:");
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->o_lambda, 4 ), redFormat );
  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "spectral index:");
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->n, 4 ), redFormat );
  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "optdlss:");
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->optdlss, 4 ), redFormat );
  ++i;
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "sigma8:");
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( m->sigma8, 4 ), redFormat );
  // qDebug() << "taking on the sigma8" << endl;
  ++i;
/*
  table->cellAt( i, 0 ).firstCursorPosition().insertText( "sigma8 O_m^gamma (0.5 +- 0.1):" );
  table->cellAt( i, 1 ).firstCursorPosition().insertText( toStr( cosmos->sigma8Omega(m->n,m->sigma8,QuintCosmos::logAWeff), 4 ) );
*/

  //qDebug() << "done" <<endl;

  if (extendedInfo.find(id) != extendedInfo.end()) {
    ExtendedInfo &e = extendedInfo[id];
    i = 0;
    table->cellAt( i, 3 ).firstCursorPosition().insertText( "   tau");
    table->cellAt( i, 3 ).lastCursorPosition().insertText( "0", subScriptFormat );
    table->cellAt( i, 3 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 4 ).lastCursorPosition().insertText( toStr( e.tau0, 4 ), blueFormat );
    i++;
    table->cellAt( i, 3 ).firstCursorPosition().insertText( "   tau");
    table->cellAt( i, 3 ).lastCursorPosition().insertText( "ls", subScriptFormat );
    table->cellAt( i, 3 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 4 ).lastCursorPosition().insertText( toStr( e.tauls, 4 ), blueFormat );
    i++;
    table->cellAt( i, 3 ).firstCursorPosition().insertText( "   tau");
    table->cellAt( i, 3 ).lastCursorPosition().insertText( "equ", subScriptFormat );
    table->cellAt( i, 3 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 4 ).lastCursorPosition().insertText( toStr( e.taueq, 4 ), blueFormat );
    i++;
    table->cellAt( i, 3 ).firstCursorPosition().insertText( "   z");
    table->cellAt( i, 3 ).lastCursorPosition().insertText( "ls", subScriptFormat );
    table->cellAt( i, 3 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 4 ).lastCursorPosition().insertText( toStr( e.zls, 4 ), blueFormat );

    table->cellAt( i, 3 ).firstCursorPosition().insertText( "   z");
    table->cellAt( i, 3 ).lastCursorPosition().insertText( "equ", subScriptFormat );
    table->cellAt( i, 3 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 4 ).lastCursorPosition().insertText( toStr( e.zeq, 4 ), blueFormat );

    i = 0;
    table->cellAt( i, 5 ).firstCursorPosition().insertText( "   WMAP chi");
    table->cellAt( i, 5 ).lastCursorPosition().insertText( "2", superScriptFormat );
    table->cellAt( i, 5 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 6 ).lastCursorPosition().insertText( toStr( e.wmap_tt + e.wmap_te, 4 ), greenFormat );
    i++;
    table->cellAt( i, 5 ).firstCursorPosition().insertText( "   WMAP (TT ) chi");
    table->cellAt( i, 5 ).lastCursorPosition().insertText( "2", superScriptFormat );
    table->cellAt( i, 5 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 6 ).lastCursorPosition().insertText( toStr( e.wmap_tt, 4 ), greenFormat );

    i++;
    table->cellAt( i, 5 ).firstCursorPosition().insertText( "   WMAP (TE) chi");
    table->cellAt( i, 5 ).lastCursorPosition().insertText( "2", superScriptFormat );
    table->cellAt( i, 5 ).lastCursorPosition().insertText( ":", QTextCharFormat() );
    table->cellAt( i, 6 ).lastCursorPosition().insertText( toStr( e.wmap_te, 4 ), greenFormat );
  }

  //qDebug() << "END INFO" << endl;
  //qDebug() << "######################################################" << endl;
  //cw->InfoView->setHtml(t);
}

void  CmbMainWindow::plotCut() {
  qDebug() << "plotCut " << RMBPlot_id << endl;
  if (RMBPlot_currentPage != 4) { // not likelihood 1-D
    if (RMBPlot_id) {
      // tell plot to get rid of it 
      cw->cmbplot->ridOf(RMBPlot_id);
      cw->structureplot->ridOf(RMBPlot_id);
      cw->TensorPlot->ridOf(RMBPlot_id);
      cw->PolarEPlot->ridOf(RMBPlot_id);
      cw->PolarCPlot->ridOf(RMBPlot_id);


      // get rid of the model ourselves
      plotModel.erase(plotModel.find(RMBPlot_id));

      if (extendedInfo.find(RMBPlot_id) != extendedInfo.end())
	extendedInfo.erase(extendedInfo.find(RMBPlot_id));   qDebug() << "deb4" << endl;
    }
  } else {  // likelihood 1-D
    cw->likeliplot1d->ridOf(RMBPlot_id);
  }
}

static QMap<int, ColumnPlotInfo> readParameterNameFile(const QString& nameFileName)
{
  QMap<int, ColumnPlotInfo> nameMap;
  int colNr;
  QString colName;
  QString enabled;
  QFile nameFile(nameFileName);
  if (nameFile.open(QFile::ReadOnly)) {
    qDebug() << "opened parameter file " << nameFileName;
    QTextStream nameStream(&nameFile);
    bool firstColRead = false;
    while (!nameStream.atEnd())   {
      nameStream >> colNr;
      // if we're only reading in junk (like a single eol before eof),
      // operator >> converts it to 0
      if (firstColRead && colNr == 0)
        break;
      firstColRead = true;

      nameStream >> enabled;
      colName = nameStream.readLine();
      qDebug() << "colNr:" << colNr << "enabled" << enabled;
      ColumnPlotInfo info;
      info.enabled = false;
      if (enabled.trimmed() == QLatin1String("x")
          && !colName.trimmed().isEmpty())    //  we have a parameter named "x" w/o enabled/disabled info
        info.enabled=true;
      else if (enabled.trimmed() != QLatin1String("o"))
        colName.prepend(enabled);

      info.columnLabel = colName;
      info.columnNumber = colNr;
      nameMap.insert(colNr, info);
    }
  }
  return nameMap;
}

static void saveParameterNameFile(const QString& nameFileName, const QMap<int, ColumnPlotInfo>& nameMap)
{
  QFile nameFile(nameFileName);
  if (nameFile.open(QFile::WriteOnly)) {
    QTextStream nameStream(&nameFile);
    QMapIterator<int, ColumnPlotInfo> i(nameMap);
    while (i.hasNext())   {
      ColumnPlotInfo info = i.next().value();
      nameStream << info.columnNumber;
      nameStream << (info.enabled?"  x ":"  o ");
      nameStream << info.columnLabel << '\n';
    }
  }
}

bool CmbMainWindow::onlyZeroesInColumn(const int col)
{
  mPlotControlDock->LikeXBox_1d->setValue( col );
  return !autoAdapt_1d();
}

static bool columnIsIncluded(const int col, const std::vector<ColumnPlotInfo>& infos)
{
  std::vector<ColumnPlotInfo>::const_iterator it;
  it = std::find_if(infos.begin(), infos.end(), std::bind1st(std::equal_to<int>(), col));
  if (it != infos.end() && (*it).enabled)
    return true;
  return false;
}

bool isPlottingEnabledForColumn(const ColumnPlotInfo& info)
{
  return info.enabled;
}

static int nrOfEnabledColumns(const std::vector<ColumnPlotInfo>& infos)
{
  return std::count_if(infos.begin(), infos.end(), isPlottingEnabledForColumn);
}

static std::string labelForColumn(const int col, const std::vector<ColumnPlotInfo>& infos)
{
  std::vector<ColumnPlotInfo>::const_iterator it;
  it = std::find_if(infos.begin(), infos.end(), std::bind1st(std::equal_to<int>(), col));
  if (it != infos.end() && (*it).enabled)
    return (*it).columnLabel.toStdString();
  return std::string();
}

std::vector<ColumnPlotInfo> CmbMainWindow::columnPlotInfoFromParameterNamesFile(const QString& nameFileName)
{
  QMap<int, ColumnPlotInfo> nameMap = readParameterNameFile(nameFileName);

  QString colName;

  std::vector<ColumnPlotInfo> retVal;
  QMapIterator<int, ColumnPlotInfo> it(nameMap);
  while (it.hasNext()) {
    ColumnPlotInfo info = it.next().value();
    if (!onlyZeroesInColumn(info.columnNumber))
        retVal.push_back(info);
  }

  return retVal;
}

std::vector<ColumnPlotInfo> CmbMainWindow::askForColumnsToPlot()
{
  QDialog* dialog = new QDialog(this);

  QTableWidget* tw = new QTableWidget(mNumberOfColumns, 2, dialog);

  QStringList labels;
  labels << "Include?" << "Name";
  tw->setHorizontalHeaderLabels(labels);

  QString proposedName = QFileInfo(currentMccFile()).absoluteFilePath();
  proposedName.chop(4); //remove the ".mcc"
  proposedName += "parameterNames";
  QString nameFileName = QFileDialog::getOpenFileName(this, "Choose file with column names",
                                              proposedName, "Parameter Name Files (*.parameterNames)", 0,
                                              QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  QMap<int, ColumnPlotInfo> nameMap = readParameterNameFile(nameFileName);

  QString colName;

  for (int col = 0; col < mNumberOfColumns; ++col) {
    QCheckBox* cb = new QCheckBox(tw);
    tw->setCellWidget(col, 0, cb);

    if (nameMap.contains(col)) {
      colName = nameMap.value(col).columnLabel;
      cb->setChecked(nameMap.value(col).enabled);
    } else {
      colName = "Column " + QString::number(col);
    }

    QTableWidgetItem* item = new QTableWidgetItem(colName);
    tw->setItem(col, 1, item);
  }


  QPushButton* ok = new QPushButton(dialog);
  ok->setDefault(true);
  ok->setText("Generate Plots and Table");

  connect(ok, SIGNAL(clicked()), dialog, SLOT(accept()));

  QVBoxLayout* layout = new QVBoxLayout();
  layout->addWidget(tw);
  layout->addWidget(ok);
  dialog->setLayout(layout);

  std::vector<ColumnPlotInfo> retVal;
  if (dialog->exec()) {
    for (int col = 0; col < mNumberOfColumns; ++col) {
      QCheckBox* box = qobject_cast<QCheckBox*>(tw->cellWidget(col,0));
      if (box){
        ColumnPlotInfo info;
        info.columnNumber = col;
        info.columnLabel = tw->item(col, 1)->text();
        bool isAllZeroes = onlyZeroesInColumn(info.columnNumber);
        info.enabled = (box->isChecked() && !isAllZeroes);
        if ( info.enabled )
          retVal.push_back(info);
        nameMap[col] = info;
      }
    }
  }
  else
  {
  }

  saveParameterNameFile( nameFileName, nameMap);

  return retVal;
}

QString CmbMainWindow::currentMccFile()
{
  if (mGenerating2dPlotTable)
    return  mPlotControlDock->LikeliBox->currentText()+".mcc";
  return  mPlotControlDock->LikeliBox_1d->currentText()+".mcc";
}

int CmbMainWindow::currentLikelihoodColumn()
{
  if (mGenerating2dPlotTable)
    return  mPlotControlDock->LikeLikeBox->value();
  return  mPlotControlDock->LikeLikeBox_1d->value();
}

void CmbMainWindow::generatePlotTable(int col, QString dir, QString parameterNamesFile, int dimensions)
{
  mPlotControlDock->Analytic->setChecked(false);
//X   qDebug() << "starting to plot columns:" << m2dPlotTableX << col << endl;
//X   QString message = "starting to plot columns:" + QString::number(m2dPlotTableX) + '-' + QString::number(col);
//X   QMessageBox::information(this, "", message);

  if ( m2dPlotTableX == 0 && col == 0) {
    if (dir.isEmpty()) {
    mAutomaticPrintDirectory =
            QFileDialog::getExistingDirectory(this, "Choose a Directory for the Plots and tex-Table",
                                              QFileInfo(currentMccFile()).absolutePath(),
                                              QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    } else if (mAutomaticPrintDirectory.isEmpty()) {
      mAutomaticPrintDirectory = dir;
    }

    if (mAutomaticPrintDirectory.isEmpty())
      return;

    bool ok;
    int dims;
    if (dimensions == 1 || dimensions == 2) {
      dims = dimensions;
      ok = true;
    } else {
      dims = QInputDialog::getInteger(this, "", "Create a 1d or 2d plot?",  1, 1, 2, 1, &ok);
    }
    if (!ok)
      return;

    if (dims == 2) {
      mGenerating2dPlotTable = true;
    } else if (dims != 1) {
      QMessageBox::critical(this, "", "Cannot create a table of plots for the currently active plot widget.");
      return;
    }


    mIgnoreBadErrors = true;

    if (mGenerating2dPlotTable)
      cw->PlotTab->setCurrentIndex(3);
    else
      cw->PlotTab->setCurrentIndex(4);

    QString mccfile = currentMccFile();
    mNumberOfColumns = ChainShop::numberOfColumns(mccfile.toStdString());

    if (parameterNamesFile.isEmpty())
      mColumnsToPlot = askForColumnsToPlot();
    else
      mColumnsToPlot = columnPlotInfoFromParameterNamesFile(parameterNamesFile);

    //debug: create only the table, if the plots are already there
//X     outputPlotTable();
//X     return;

    if (mColumnsToPlot.empty()) {
      qDebug() << "no columns to plot";
      return;
    }
  }
  if ( col < mNumberOfColumns ) {
    std::vector<ColumnPlotInfo>::const_iterator it;
    bool columnIncluded = columnIsIncluded(col, mColumnsToPlot);

    if (columnIncluded && !mGenerating2dPlotTable) {
      mPlotControlDock->LikeXBox_1d->setValue( col );
      autoAdapt_1d();
      connect( this, SIGNAL( likelihoodDrawingFinished() ), this, SLOT( outputPlotForTable() ) );
      startDrawLikeli_1d();
      return;
    }
    else if (columnIncluded && mGenerating2dPlotTable && (m2dPlotTableX != col)) {
      mPlotControlDock->LikeXBox->setValue( m2dPlotTableX );
      mPlotControlDock->LikeYBox->setValue( col );
      autoAdapt();
      connect( this, SIGNAL( likelihoodDrawingFinished() ), this, SLOT( outputPlotForTable() ) );
      startDrawLikeli_2d();
      return;
    } else {
      qDebug() << "not plotting " << m2dPlotTableX << col;
      if (mGenerating2dPlotTable) {
        mPlotControlDock->LikeYBox->stepUp();
        generatePlotTable(mPlotControlDock->LikeYBox->value());
        return;
      } else {
        mPlotControlDock->LikeXBox_1d->stepUp();
        generatePlotTable(mPlotControlDock->LikeXBox_1d->value());
        return;
      }
    }
  }
  if ( col == mNumberOfColumns ) {
      if (mGenerating2dPlotTable && (m2dPlotTableX < mNumberOfColumns)) {
        ++m2dPlotTableX;
        while ( (m2dPlotTableX < mNumberOfColumns)
                && !columnIsIncluded(m2dPlotTableX, mColumnsToPlot))
          ++m2dPlotTableX;

        if (m2dPlotTableX == mNumberOfColumns)
          generatePlotTable(mPlotControlDock->LikeYBox->value()); // i.e. go do outputPlotTable()

        mPlotControlDock->LikeXBox->setValue(m2dPlotTableX);
        col = 0;
        mPlotControlDock->LikeYBox->setValue(col);
        generatePlotTable(mPlotControlDock->LikeYBox->value());
        return;
      } else {
        mIgnoreBadErrors = false;
        outputPlotTable();
        qDebug() << "done" << endl;
        mGenerating2dPlotTable = false;
      }
  }
}

void CmbMainWindow::outputPlotForTable()
{
  disconnect( this, SIGNAL( likelihoodDrawingFinished() ), this, SLOT( outputPlotForTable() ) );
  printAutomatic();
  if (mGenerating2dPlotTable) {
    mPlotControlDock->LikeYBox->stepUp();
    generatePlotTable(mPlotControlDock->LikeYBox->value());
  } else {
    mPlotControlDock->LikeXBox_1d->stepUp();
    generatePlotTable(mPlotControlDock->LikeXBox_1d->value());
  }
}

QString CmbMainWindow::automaticPrintFileName(int colX, int colY)
{
  QString fileName = QFileInfo(currentMccFile()).fileName();
  fileName.chop(4); // remove the ".mcc"
  QString thisParameter;
  if (mGenerating2dPlotTable) {
    if (colX == -1) {
      thisParameter = mPlotControlDock->LikeXBox->cleanText();
    } else {
      thisParameter = mPlotControlDock->LikeXBox->textFromValue(colX).trimmed();
    }
    fileName += "_2d_" + thisParameter + '-';
  } else {
    fileName += "_1d_";
  }
  if (colX == -1 && colY == -1) {
    if (mGenerating2dPlotTable) {
      fileName +=  mPlotControlDock->LikeYBox->cleanText();
    } else {
      fileName +=  mPlotControlDock->LikeXBox_1d->cleanText();
    }
  } else {
    if (mGenerating2dPlotTable) {
      fileName +=  mPlotControlDock->LikeYBox->textFromValue(colY).trimmed();
    } else {
      qDebug() << "1_d_textfrom value" <<  colX <<  mPlotControlDock->LikeXBox_1d->textFromValue(colX);
      fileName +=  mPlotControlDock->LikeXBox_1d->textFromValue(colX).trimmed();
    }
  }
  qDebug() << "automaticPrintFileName:" << colX << colY << "returning" << fileName;
  return mAutomaticPrintDirectory + fileName + ".eps";
}

void CmbMainWindow::printAutomatic()
{
  PlotWidget *p = currentPlotWidget();
  if (p)  {
    try {
      int idx =cw->PlotTab->currentIndex();
      PostscriptPlot post(PrtProfiles[idx]);
      post.setResourceFile(ControlPanel::cmbeasyDir("/resources/advancedpostscript.eps"));

      p->autoTick(post);

      ThingsToPlot Things;
      //check whether we want to plot the star marking the best fit model
      if ( mPlotControlDock->BestFitBox->isChecked() ) {
        CoordPoint pos( cw->likeliplot->BestFitModel );
        SpecialPoint p( SpecialPoint::Star, pos, 24e-3 );
        Things.SpecialPoints.push_back(p);
      }

      if ( mPlotControlDock->MarkReferenceModelBox->isChecked() ) {
        CoordPoint pos( cw->likeliplot->ModelToMark );
        SpecialPoint p( SpecialPoint::Star, pos, 24e-3 );
        Things.SpecialPoints.push_back(p);
      }

      std::string xAxisLabelText = PrtProfiles[idx].XLabelText;

      int col = mPlotControlDock->LikeXBox_1d->value();
      if (mGenerating2dPlotTable)
        col = mPlotControlDock->LikeXBox->value();

      bool columnIncluded = columnIsIncluded(col, mColumnsToPlot);

      if (columnIncluded)
        xAxisLabelText = labelForColumn(col, mColumnsToPlot);
      else {
        stringstream columnNumberStr;
        columnNumberStr << col;
        throw Bad_Error("CmbMainWindow::printAutomatic - no label text available for col no " + columnNumberStr.str());
      }

      AxisLabel xaxis(post.parseTexLabel(xAxisLabelText));  
      xaxis.offset = post.XLabelOffset;
      xaxis.size= 0.08*0.01*post.AxisLabelSize; 
      Things.XAxisLabel = xaxis;

      std::string yAxisLabelText = PrtProfiles[idx].YLabelText;
      if (mGenerating2dPlotTable) {
        int ycol = mPlotControlDock->LikeYBox->value();
        bool columnIncluded = columnIsIncluded(ycol, mColumnsToPlot);

        if (columnIncluded)
          yAxisLabelText = labelForColumn(ycol, mColumnsToPlot);
        else {
          stringstream columnNumberStr;
          columnNumberStr << ycol;
          throw Bad_Error("CmbMainWindow::printAutomatic - no label text available for y-col no " + columnNumberStr.str());
        }
      }

      AxisLabel yaxis(post.parseTexLabel(yAxisLabelText));  
      yaxis.offset = post.YLabelOffset;
      yaxis.size= 0.08*0.01*post.AxisLabelSize; 
      Things.YAxisLabel = yaxis;

      p->printPostscript(post,Things, automaticPrintFileName());
    } catch (Bad_Error x) {
      qWarning() << "bad error during the creation of a table of the 1d plots." << endl;
      qWarning() << QString::fromStdString(x.s)  << endl;
      QMessageBox::critical(this,"Bad Error", x.s.c_str(),  QMessageBox::Ok);
    }
  }
}


void CmbMainWindow::outputPlotTable()
{
  QString info;

  QString currentFile = currentMccFile();
  currentFile.chop(4); // remove the ".mcc"
  QString name = QFileInfo(currentFile).fileName();

  info += "\\documentclass[]{revtex4}\n";
  info += "\\usepackage[]{graphicx}\n";
  info += "\\begin{document}\n";

  if (mGenerating2dPlotTable)
    info += "Summary information (2d-plots) for the monte carlo file: ";
  else
    info += "Summary information (1d-plots) for the monte carlo file: ";

  info += QFileInfo(name).fileName().replace('_', QString("\\_"));
  info += '\n';
  info += '\n';

  int nrOfColumns = 3;
  if (mGenerating2dPlotTable)
    nrOfColumns = nrOfEnabledColumns(mColumnsToPlot)-1;
  info += "\\begin{tabular}{c";
  for (int c=0; c < nrOfColumns-2; ++c)
    info += "c";
  info += "c}\n";

  int colsPrinted = 1;
  for (int col = 0; col < mNumberOfColumns; ++col) {
    bool columnIncluded = columnIsIncluded(col, mColumnsToPlot);

    if (!columnIncluded)
      continue;

    m2dPlotTableX = 0;
    int nrOfXColumns = 0;
    if (mGenerating2dPlotTable)
      nrOfXColumns = mNumberOfColumns-1;
    int xcolsPrinted=0;
    for ( ; m2dPlotTableX <= nrOfXColumns; ++m2dPlotTableX) {
      qDebug() << "at" << col << m2dPlotTableX;
      bool columnIncluded2 = columnIsIncluded(m2dPlotTableX, mColumnsToPlot);

         if (mGenerating2dPlotTable && (!columnIncluded2 || (m2dPlotTableX==col))) {
           continue;
         }
      QString fileName = QFileInfo(automaticPrintFileName(col, m2dPlotTableX)).fileName();
      QString width;
      width.setNum(1./(nrOfColumns+1), 'f', 3);
      if (!mGenerating2dPlotTable || col < m2dPlotTableX)
        info += QLatin1String("    ") + "\\includegraphics[angle=270, width=" + width + "\\textwidth]{" + fileName + "}";
      ++colsPrinted; ++xcolsPrinted;
      if ( (!mGenerating2dPlotTable && ((colsPrinted%nrOfColumns)&1))
           || (mGenerating2dPlotTable && (xcolsPrinted>=nrOfColumns)))
        info += QLatin1String("\\\\") + "\n";
      else
        info += "& \n";
    }
  }
  info += '\n';

  info += "\\end{tabular}";
  info += '\n';
  info += "\\end{document}\n";

  QString outputFileName = mAutomaticPrintDirectory + QFileInfo(currentMccFile()).fileName();
  outputFileName.chop(4); //remove the .mcc
  if (mGenerating2dPlotTable)
    outputFileName += "-2dplots.tex";
  else
    outputFileName += "-1dplots.tex";

  QFile file( outputFileName );
  if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) ) {
    QMessageBox::critical(this, "Error", "Could not write to " +  outputFileName, QMessageBox::Ok );
    return;
  }

  QTextStream out( &file );
  out << info;

  mGenerating2dPlotTable = false;
  m2dPlotTableX = 0;

  if (!AutomaticRunner::runningNonInteractively) {
    QMessageBox::information(this,"Done", "Table written to: " + outputFileName , QMessageBox::Ok);
  } else {
    qDebug() << "Table written to: " + outputFileName;
    emit summaryDone();
  }
}


void CmbMainWindow::generateStatisticsFile(int col, QString dir, QString parameterNames)
{
  qDebug() << "CmbMainWindow::generateStatisticsFile " << dir << mAutomaticPrintDirectory;
  if (dir.isEmpty() && mAutomaticPrintDirectory.isEmpty()) {
    mAutomaticPrintDirectory =
      QFileDialog::getExistingDirectory(this, "Choose a Directory for the statistics file",
          QFileInfo(currentMccFile()).absolutePath(),
          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
  } else if (mAutomaticPrintDirectory.isEmpty()) {
    mAutomaticPrintDirectory = dir;
  }

  if (mAutomaticPrintDirectory.isEmpty())
    return;


  mIgnoreBadErrors = true;
  cw->PlotTab->setCurrentIndex(4);
  if ( col == 0 ) {
    mNumberOfColumns = ChainShop::numberOfColumns(currentMccFile().toStdString());
    if (parameterNames.isEmpty())
      m1dStatisticsColumns = askForColumnsToPlot();
    else
      m1dStatisticsColumns = columnPlotInfoFromParameterNamesFile(parameterNames);
    mLikelihoodInfos.clear();
    mPlotControlDock->LikeXBox_1d->setValue(0);
  }
  if ( col < mNumberOfColumns ) {
    mPlotControlDock->LikeXBox_1d->setValue( col );
    autoAdapt_1d();
    connect( this, SIGNAL( likelihoodDrawingFinished() ), this, SLOT( collectLikelihoodInfos() ) );
    mLastLikelihoodInfos.clear();
    QTimer::singleShot(0, this, SLOT(startDrawLikeli_1d()));
    return;
  }
  if ( col == mNumberOfColumns ) {
    mIgnoreBadErrors = false;
    askForStatisticsFileFormat();
  }
}

void CmbMainWindow::askForStatisticsFileFormat()
{
  if (AutomaticRunner::runningNonInteractively) {
    outputStatisticsTexFile();
    emit summaryDone();
    return;
  }
  QMessageBox msgBox;
  msgBox.setWindowTitle("Cmbeasy - Statistics Output Format");
  msgBox.setText("Which output format should be used for the summary? ");

  msgBox.setStandardButtons(QMessageBox::Cancel);
  QPushButton *plainTextButton = msgBox.addButton("Plain Text", QMessageBox::AcceptRole);
  QPushButton *texButton = msgBox.addButton("Tex Table", QMessageBox::AcceptRole);


  msgBox.exec();

  if (msgBox.clickedButton() == plainTextButton) {
    outputStatisticsFile();
  } else if (msgBox.clickedButton() == texButton) {
    outputStatisticsTexFile();
  }
}

void CmbMainWindow::collectLikelihoodInfos()
{
  disconnect( this, SIGNAL( likelihoodDrawingFinished() ), this, SLOT( collectLikelihoodInfos() ) );
  if ( !mLastLikelihoodInfos.empty() ) {
    mLikelihoodInfos.push_back( mLastLikelihoodInfos );
  } else {
    vector<LikelihoodInfo> v;
    v.push_back( LikelihoodInfo( LikelihoodInfo::Error ) );
    mLikelihoodInfos.push_back( v );
  }
  mLastLikelihoodInfos.clear();
  mPlotControlDock->LikeXBox_1d->stepUp();
  generateStatisticsFile( mPlotControlDock->LikeXBox_1d->value() );
}

//size of one field in the statistics output
#define TABSIZE 14

static QString makeString( double value, QString prefix = QString() )
{
  QString s = QString::number ( value ).prepend( prefix ).rightJustified( TABSIZE, ' ' );
  s.prepend( ' ' );
  return s;
}

static QString makeString( const QString& str )
{
  QString s = str.rightJustified( TABSIZE, ' ' );
  s.prepend( ' ' );
  return s;
}

void CmbMainWindow::outputStatisticsFile() 
{
  QString info;


  QString name = mPlotControlDock->LikeliBox_1d->currentText();
  int likepos = mPlotControlDock->LikeLikeBox_1d->value();
  vector<float> bestModel = bestFitModel(name+QLatin1String(".mcc"),likepos);
  typedef vector<LikelihoodInfo> LikeInfoVec;
  typedef vector<LikeInfoVec> LikeInfosCols;

  LikeInfosCols::iterator col = mLikelihoodInfos.begin();
  int colNr = 0;
  LikeInfosCols::iterator lastCol = mLikelihoodInfos.end();
  info += name;
  info += '\n';
  info += "-------------------------------------------------------------------------------";
  info += '\n';

  QString line;
  QString emptyTab;
  emptyTab.fill( ' ', TABSIZE );
  line += emptyTab;
  const QString colNames[4] = { "Chi2", "68% under curve: ", "95% under curve: ", "Best"};
  for ( int i = 0; i < 4; ++i ) {
    QString curColName = colNames[ i ].leftJustified( TABSIZE, ' ' );
    line += emptyTab + curColName  + emptyTab;
  }
  info +=  line + '\n';
  for ( ; col != lastCol; ++col, ++colNr ) {
    std::vector<ColumnPlotInfo>::const_iterator plotInfoIt;
    plotInfoIt = std::find(m1dStatisticsColumns.begin(), m1dStatisticsColumns.end(), colNr);
    bool shouldBeIncluded = (plotInfoIt != m1dStatisticsColumns.end());
    if (!shouldBeIncluded)
      continue;
    ColumnPlotInfo curPlotInfo = *plotInfoIt;
    if (!curPlotInfo.enabled)
      continue;

    QString chi2Str, underCurve95Str, underCurve68Str;
    LikeInfoVec::iterator it = col->begin();
    LikeInfoVec::iterator end = col->end();
    qDebug() << colNr << " has " << col->size() << " infos. " << endl;
    for ( ; it != end; ++it ) {
      QString currentInfoStr;
      LikelihoodInfo i = *it;
      currentInfoStr += ( i.hasMaximum )?  makeString( i.maximum() ) : makeString( "n.a." );
      currentInfoStr += ( i.hasLeft )?  makeString( i.left(), "(-)" ) : makeString( "n.a." );
      currentInfoStr +=  ( i.hasRight )?  makeString( i.right(), "(+)" ) : makeString( "n.a." );
      switch ( i.type ) {
        case LikelihoodInfo::Chi2: chi2Str = currentInfoStr; break;
        case LikelihoodInfo::UnderCurve95: underCurve95Str =  currentInfoStr; break;
        case LikelihoodInfo::UnderCurve68: underCurve68Str = currentInfoStr; break;
        case LikelihoodInfo::None: /*fall through*/;
        default: /* shouldn't happen */;
      }
    }
    qDebug() << "column nr: " << colNr << endl;

    QString rowDesc = curPlotInfo.columnLabel.leftJustified( TABSIZE, ' ' );
    QString valueStr = chi2Str + underCurve68Str + underCurve95Str;
    if ( valueStr.isEmpty() ) valueStr.fill( ' ', 10 * TABSIZE );
    valueStr += makeString(bestModel[colNr]);
    //if ( QString errorStr( "Error during assembling likelihood infos... " );
    info += rowDesc + valueStr+ '\n';
  }
  qDebug() << "info: " << info << endl;

  QFile file( name + ".statistics" );
  if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) )
    return;

  QTextStream out( &file );
  out << info;
}

void CmbMainWindow::outputStatisticsTexFile()
{
  QString info;

  QString name = mPlotControlDock->LikeliBox_1d->currentText();
  int likepos = mPlotControlDock->LikeLikeBox_1d->value();
  vector<float> bestModel = bestFitModel(name+QLatin1String(".mcc"),likepos);
  typedef vector<LikelihoodInfo> LikeInfoVec;
  typedef vector<LikeInfoVec> LikeInfosCols;

  info += "\\documentclass[draft]{revtex4}\n";
  info += "\\begin{document}\n";
  info += "Summary information for the monte carlo file: ";
  info += QFileInfo(name).fileName().replace('_', QString("\\_"));
  info += '\n';
  info += '\n';
  info += "\\begin{tabular}{c|c|c|c|c|c|}";
  info += '\n';

  QString line;
  QString emptyTab;
  emptyTab.fill( ' ', TABSIZE );
  line += " - ";
  const QString colNames[4] = { "Chi2", "68\\% under curve: ", "95\\% under curve: ", "Best"};
  for ( int i = 1; i < 4; ++i ) {
    QString curColName = colNames[ i ].leftJustified( TABSIZE, ' ' );
    line += '&'+ emptyTab + curColName  + emptyTab ;
  }
  info +=  line + "\\\\ \\hline \n";
  int colNr = 0;
  LikeInfosCols::iterator col = mLikelihoodInfos.begin();
  LikeInfosCols::iterator lastCol = mLikelihoodInfos.end();
  for ( ; col != lastCol; ++col, ++colNr ) {
    std::vector<ColumnPlotInfo>::const_iterator plotInfoIt;
    plotInfoIt = std::find(m1dStatisticsColumns.begin(), m1dStatisticsColumns.end(), colNr);
    bool shouldBeIncluded = (plotInfoIt != m1dStatisticsColumns.end());
    if (!shouldBeIncluded)
      continue;
    ColumnPlotInfo curPlotInfo = *plotInfoIt;
    if (!curPlotInfo.enabled)
      continue;


    QString chi2Str, underCurve95Str, underCurve68Str;
    LikeInfoVec::iterator it = col->begin();
    LikeInfoVec::iterator end = col->end();
    qDebug() << colNr << " has " << col->size() << " infos. " << endl;
    for ( ; it != end; ++it ) {
      QString currentInfoStr;
      LikelihoodInfo i = *it;

      if (!i.hasMaximum && !i.hasRight && i.hasLeft)
        currentInfoStr += " > " + makeString(i.left());
      else if (!i.hasMaximum && !i.hasLeft && i.hasRight)
        currentInfoStr += " < " + makeString(i.right());
      else
      {
        currentInfoStr += ( i.hasMaximum )?  makeString( i.maximum() ) : makeString( "n.a." );
        currentInfoStr += "^{" + (( i.hasLeft )?  makeString( i.left(), "(-)" ) : makeString( "n.a." )) + "}";
        currentInfoStr += "_{" +  (( i.hasRight )?  makeString( i.right(), "(+)" ) : makeString( "n.a." )) + "}";
      }
      switch ( i.type ) {
        case LikelihoodInfo::Chi2: chi2Str = currentInfoStr; break;
        case LikelihoodInfo::UnderCurve95: underCurve95Str =  currentInfoStr; break;
        case LikelihoodInfo::UnderCurve68: underCurve68Str = currentInfoStr; break;
        case LikelihoodInfo::None: /*fall through*/;
        default: /* shouldn't happen */;
      }
    }
    qDebug() << "column nr: " << colNr << endl;
    QString rowDesc = curPlotInfo.columnLabel.leftJustified( TABSIZE, ' ' );
    rowDesc.prepend('$');
    rowDesc.append('$');
    QString valueStr =/* " & $" +chi2Str + " $ " + */ " & $ " + underCurve68Str +  " $ & $ " + underCurve95Str + " $ & ";
    if ( valueStr.isEmpty() ) valueStr.fill( ' ', 10 * TABSIZE );
    valueStr += makeString(bestModel[colNr]);
    //if ( QString errorStr( "Error during assembling likelihood infos... " );
    if (!(/*chi2Str.isEmpty() &&*/ underCurve95Str.isEmpty() && underCurve68Str.isEmpty()))
      info += rowDesc + valueStr+ " \\\\ \n";
  }

  info += "\\end{tabular}";
  info += '\n';
  info += "\\end{document}\n";
  qDebug() << "info:" << info;

  if (!mAutomaticPrintDirectory.endsWith('/'))
      mAutomaticPrintDirectory.append('/');
  QFileInfo fileInf(name);
  QFile file( mAutomaticPrintDirectory + fileInf.fileName() + "-statistics.tex" );
  if ( !file.open( QIODevice::WriteOnly | QIODevice::Text ) ) {
    qDebug() << "Error writing to: " << file.fileName();
    return;
  }

  QTextStream out( &file );
  out << info;
  qDebug() << "wrote tex table to: " << file.fileName();
}

void CmbMainWindow::savePrtProfile()
{
	bool ok;
	QString text = QInputDialog::getText( this, "Please enter a name for this profile",
			"profile name:", QLineEdit::Normal,
			"unnamed profile", &ok );
	if ( ok && !text.isEmpty() ){
		cfg->syncProfile();
		LowLevelPlot pr( PrtProfiles[ cfg->CurrentProfile ] );
		pr.Name = text.toStdString();
		PrtProfiles.push_back( pr );
	}
	cfg->updateProfiles();
}

void CmbMainWindow::print() {
  PlotWidget *p = currentPlotWidget();
  //  if (p != cw->likeliplot) { // most of the plots are easy...
  if (true) {
    if (p)  {
      try { 
        //PostscriptPlot post(ControlPanel::cmbeasyDir("/resources/advancedpostscript.eps"));
        cfg->syncProfile();
        if (cfg->CurrentProfile == -1)
          cfg->setProfile(cw->PlotTab->currentIndex());
        int idx = cfg->CurrentProfile;
        PostscriptPlot post(PrtProfiles[idx]);
        qDebug() << "printing using profile: " << cfg->CurrentProfile << QString::fromStdString(post.Name);
        post.setResourceFile(ControlPanel::cmbeasyDir("/resources/advancedpostscript.eps"));

        if (post.AutomaticTick) {
          qDebug() << "start autoticking" << endl;
          p->autoTick(post);
        } else { 
          p->semiAutoTick(post);
        }


        ThingsToPlot Things;
        //check whether we want to plot the star marking the best fit model
        if ( mPlotControlDock->BestFitBox->isChecked() ) {
          CoordPoint pos( cw->likeliplot->BestFitModel );
          SpecialPoint p( SpecialPoint::Star, pos, 24e-3 );
          Things.SpecialPoints.push_back(p);
        }
        if ( mPlotControlDock->MarkReferenceModelBox->isChecked() ) {
          CoordPoint pos( cw->likeliplot->ModelToMark );
          SpecialPoint p( SpecialPoint::Star, pos, 24e-3 );
          Things.SpecialPoints.push_back(p);
        }

        AxisLabel xaxis(post.parseTexLabel(PrtProfiles[idx].XLabelText));  
        xaxis.offset = post.XLabelOffset;
        xaxis.size= 0.08*0.01*post.AxisLabelSize; 
        Things.XAxisLabel = xaxis;
        AxisLabel yaxis(post.parseTexLabel(PrtProfiles[idx].YLabelText));  
        yaxis.offset = post.YLabelOffset;
        yaxis.size= 0.08*0.01*post.AxisLabelSize; 
        Things.YAxisLabel = yaxis;

        p->printPostscript(post,Things);
      } catch (Bad_Error x) {
        qDebug() << "endspool bad error" << endl;
        qDebug() << QString::fromStdString(x.s)  << endl;
        QMessageBox::critical(this,"Bad Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
      }
    }
  } else { // if it is the likeli-plot
    throw Bad_Error("why here");
    /*
    QSettings s;
    prt->setOutputToFile(s.value("/cmbeasy/plot/printqt/tofile",false).toBool());
    prt->setOutputFileName(s.value("/cmbeasy/plot/printqt/name","print.ps").toString());
    prt->setOrientation((QPrinter::Orientation)s.readNumEntry("/cmbeasy/plot/printqt/orientation",QPrinter::Landscape));
    
    if (prt->setup()) {
      s.setValue("/cmbeasy/plot/printqt/tofile",prt->outputToFile());
      s.setValue("/cmbeasy/plot/printqt/name",prt->outputFileName());
      s.setValue("/cmbeasy/plot/printqt/orientation",prt->orientation());
      
      if (prt->outputToFile()) postscriptFile = prt->outputFileName().latin1();
      
      // a4 paper:   x = 595.3    y = 841.9 
      double scaleX=595.3,scaleY=841.9;
      double t;
      switch (prt->pageSize()) {
      case QPrinter::A5: t=scaleX; scaleX=scaleY /2; scaleY=t; break;
      case QPrinter::A3: t=scaleX; scaleX=scaleY; scaleY=2*t; break;
      default: {};
      }
    }     
    printLikeli();  
    */
  } 
}

PlotWidget* CmbMainWindow::idx2plot(int idx) {
  PlotWidget *p=0;
  if (idx == 0) p = cw->cmbplot;  //cmb
  if (idx == 1) p = cw->TensorPlot;
  if (idx == 2) p = cw->structureplot; // cdm_power
  if (idx == 3) p = cw->likeliplot; // likelihood 2-dim
  if (idx == 4) p = cw->likeliplot1d; // likelihood 1-dim
  if (idx == 5) p = cw->PolarEPlot;
  if (idx == 6) p = cw->PolarCPlot;
  return p;
}

int CmbMainWindow::plot2idx(PlotWidget* p) {
  int idx=0;
  if (p == cw->cmbplot) idx = 0;
  if (p == cw->TensorPlot) idx = 1;
  if (p == cw->structureplot) idx =2;
  if (p == cw->likeliplot) idx=3;
  if (p == cw->likeliplot1d) idx=4;
  if (p == cw->PolarEPlot) idx =5;
  if (p == cw->PolarCPlot) idx=6;
  return idx;
}

void CmbMainWindow::toFront() {
  PlotWidget *p = currentPlotWidget();
  if (!p) return;
  p->toFront(RMBPlot_id);
  p->update();
}

void CmbMainWindow::toBack() {
  PlotWidget *p = currentPlotWidget();
  if (!p) return;
  p->toBack(RMBPlot_id);
  p->update();
}

void CmbMainWindow::deleteRegion() {
  PlotWidget *p = currentPlotWidget();
  if (!p) return;
  p->deleteRegionGroup(RMBPlot_id);
  p->update();
}

void CmbMainWindow::colorRegion() {
  PlotWidget *p = currentPlotWidget();
  if (!p) return;
  p->colorRegionGroup(RMBPlot_id);
  p->update();
}

void CmbMainWindow::stickyRegion() {
  PlotWidget *p = currentPlotWidget();
  if (!p) return;
  p->setRegionGroupSticky(RMBPlot_id);
  p->update();
}


PlotWidget* CmbMainWindow::currentPlotWidget() {
  return  idx2plot(cw->PlotTab->currentIndex());
}

void CmbMainWindow::printGnuplot() {
  PlotWidget *p = currentPlotWidget();
  if (p)  p->printGnuplot();
}


void CmbMainWindow::autoScale() {
  currentPlotWidget()->autoScale();
  emit(repaintPlot());
}

void CmbMainWindow::syncPlot() {
  //  qDebug() << "sync plot" << endl;
  int idx = cw->PlotTab->currentIndex();
  PlotWidget *p=currentPlotWidget();
  //qDebug() << "index: "<<idx << endl;

  noPingPong=true;
  mPlotControlDock->PlotXMin->setText(toStr(p->Co.x,5));
  mPlotControlDock->PlotXMax->setText(toStr(p->Co.X,5));
  
  mPlotControlDock->XLogBox->setChecked(p->logX());
  mPlotControlDock->YLogBox->setChecked(p->logY());
  noPingPong=false;

  if (idx != 0) mPlotControlDock->DataBox->hide();
  if (idx == 0) mPlotControlDock->DataBox->show();
}

void CmbMainWindow::resyncPlot() {
  qDebug() << "Resyncplot" << endl;
  if (noPingPong) return;
  //return;
  int idx = cw->PlotTab->currentIndex();
  PlotWidget *p = idx2plot(idx);
 
  qDebug() << "idx: " << idx << endl;

  

  if (mPlotControlDock->PlotXMin->text().length() > 0) p->Co.x = mPlotControlDock->PlotXMin->text().toDouble();
  if (mPlotControlDock->PlotXMax->text().length() > 0) p->Co.X = mPlotControlDock->PlotXMax->text().toDouble();

  p->setLog(mPlotControlDock->XLogBox->isChecked(), mPlotControlDock->YLogBox->isChecked());
  
  qDebug() << "AUTOSCALING IN RESYNC()" << endl;
  // autoscale without notify
  p-> autoScale(false,false); // cause scaleChanged() is connected to us
  p->repaint();
}


void CmbMainWindow::askCosmos() {
  
  if (plotModel.find(infoModel) == plotModel.end())  return;
  Model* m = plotModel[infoModel];


  cosmos->reset();
  cosmos->seth(m->h);
  cosmos->setOmega_b(m->o_b);
  cosmos->setOmega_cdm(m->o_cdm);
  cosmos->setInitialPower(m->n);
  cosmos->setOptDistanceLss(m->optdlss);

  ControlPanel *tmp = new ControlPanel(*control);

  control->cmb = false;
  control->power_cdm = false;
  cmbcalc->go(cosmos,"dudei",*control);
  
  ExtendedInfo e;
  extendedInfo[infoModel] =   fillExtendedInfo(e); 

  updateInfo(infoModel);
  mPlotControlDock->BigTab->setCurrentIndex(0);
  delete control;
  control = tmp;
}

ExtendedInfo& CmbMainWindow::fillExtendedInfo(ExtendedInfo& e) {

 e.tau0 = cosmos->tau_0();
 e.tauls = cosmos->tau_ls();
 e.taueq = cosmos->tau_equ();
 
 e.zls = cosmos->z_ls();
 e.zeq = cosmos->z_equ(); 

 e.o_q = cosmos->omega_q();
 e.QuintType = cosmos->quintessence()->name().c_str();
 
 
 return e;

}


void CmbMainWindow::syncControl() {
  control->power_cdm = mParameterDock->PowerBox->isChecked();
  control->scalar = mParameterDock->Scalar->isChecked() || control->power_cdm;
  control->tensor  = mParameterDock->Tensor->isChecked();
  control->cmb = mParameterDock->Scalar->isChecked() || control->tensor;
  if (mParameterDock->Lensed->isChecked() && mParameterDock->Lensed->isEnabled()) 
    control->setLensing(ControlPanel::linear);
  else  control->setLensing(ControlPanel::none);
}


bool CmbMainWindow::autoAdapt() {
  qDebug() << "auto adapt splines" << endl;
  QString mccfile = mPlotControlDock->LikeliBox->currentText()+".mcc";
  int xpos = mPlotControlDock->LikeXBox->value();
  int ypos = mPlotControlDock->LikeYBox->value();
  int likepos = mPlotControlDock->LikeLikeBox->value();

  vector<float> bestModel;
  if ( mPlotControlDock->BestFitBox->isChecked() )
    bestModel = bestFitModel(mccfile,likepos);
  vector<float> referenceModel;
  if ( mPlotControlDock->MarkReferenceModelBox->isChecked() )
    referenceModel = modelToMark();

  vector<float> v;
  try {
    v= ChainShop::autoScale(mccfile.toStdString(),xpos,ypos);
    if ( mPlotControlDock->BestFitBox->isChecked() )
      adjustRangeForModel(xpos, ypos, bestModel, v);
    if ( mPlotControlDock->MarkReferenceModelBox->isChecked() )
      adjustRangeForModel(xpos, ypos, referenceModel, v);
    int precision = 3;
    while (QString::number(v[0], 'f' /*format as 99.9.., not 9.9e1*/, precision)
               == QString::number(v[1], 'f', precision++)) {}
    QString xminString = QString::number(v[0], 'f', precision);
    QString xmaxString = QString::number(v[1], 'f', precision);
    precision = 3;
    while (QString::number(v[2], 'f', precision)
               == QString::number(v[3], 'f', precision++)) {}
    QString yminString = QString::number(v[2], 'f', precision);
    QString ymaxString = QString::number(v[3], 'f', precision);
    mPlotControlDock->LikeMinX->setText(xminString);
    mPlotControlDock->LikeMaxX->setText(xmaxString);
    mPlotControlDock->LikeMinY->setText(yminString);
    mPlotControlDock->LikeMaxY->setText(ymaxString);
  }   catch(Bad_Error x) {
    qDebug() << "BAD ERROR IN DRAWLIKE: " << endl;
    qDebug() << QString::fromStdString(x.s) << endl;
    QMessageBox::critical(this,"AutoAdapt Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
    return false;
  }

  if (v[0]==0 && v[1]==0)
    return false;
  if (v[2]==0 && v[3]==0)
    return false;
  return true;
}

bool CmbMainWindow::autoAdapt_1d() {
  QString mccfile = mPlotControlDock->LikeliBox_1d->currentText()+".mcc";
  int xpos = mPlotControlDock->LikeXBox_1d->value();
  int likepos = mPlotControlDock->LikeLikeBox_1d->value();
  vector<float> bestModel;
  if ( mPlotControlDock->BestFitBox->isChecked() )
    bestModel = bestFitModel(mccfile,likepos);
  vector<float> referenceModel;
  if ( mPlotControlDock->MarkReferenceModelBox->isChecked() )
    referenceModel = modelToMark();
  vector<float> v;
  try {
    v= ChainShop::autoScale(mccfile.toStdString(),xpos,xpos);
    if ( mPlotControlDock->BestFitBox->isChecked() )
      adjustRangeForModel(xpos, xpos, bestModel, v);
    if ( mPlotControlDock->MarkReferenceModelBox->isChecked() )
      adjustRangeForModel(xpos, xpos, referenceModel, v);
    qDebug() << "autoadapt: " << v[0] << " to: " << v[1] << endl;
    int precision = 3;
    while (QString::number(v[0], 'f' /*format as 99.9.., not 9.9e1*/, precision)
               == QString::number(v[1], 'f', precision++)) {}
    mPlotControlDock->LikeMinX_1d->setText(QString::number(v[0], 'f', precision));
    mPlotControlDock->LikeMaxX_1d->setText(QString::number(v[1], 'f', precision));
  }   catch(Bad_Error x) {
    qDebug() << QString::fromStdString(x.s);
    QMessageBox::critical(this,"AutoAdapt Error", x.s.c_str(),  QMessageBox::Ok,QMessageBox::NoButton);
  }
  if (v[0]==0 && v[1]==0)
    return false;
  return true;
}

void CmbMainWindow::adjustRangeForModel(int xPos, int yPos, const std::vector<float>& model, std::vector<float>& range) const
{
  const float xRange = range[1]-range[0];
  const float yRange = range[3]-range[2];
  range[0] = qMin(range[0], model[xPos]-0.05f*xRange);
  range[1] = qMax(range[1], model[xPos]+.05f*xRange);
  range[2] = qMin(range[2], model[yPos]-.05f*yRange);
  range[3] = qMax(range[3], model[yPos]+0.05f*yRange);
}

/*
void CmbMainWindow::syncPaperSizeTo1(int id) {
  cw->PaperSize->setCurrentItem(id);
}
void CmbMainWindow::syncPaperSizeTo2(int id) {
  //  cw->PaperSize_1d->setCurrentItem(id);
}
*/

void CmbMainWindow::advancedPlotSettings() {
  int idx = cw->PlotTab->currentIndex();
  if (cfg->CurrentProfile == -1)
    cfg->setProfile(idx);
  cfg->show();
}

void CmbMainWindow::lockIn1d() {
  if (IsLockedIn1d) {
    InformLike1d = "";
  } else {
    InformLike1d = "<font color=#aa0000>I locked in on: " +  mPlotControlDock->LikeliBox_1d->currentText();
    InformLike1d += "</font><br>";
    LockedExpPoly = CurrentExpPoly;
  }
  IsLockedIn1d=!IsLockedIn1d;
  mPlotControlDock->LikeText_1d->setHtml(InformLike1d);
}

void CmbMainWindow::keep2dStateChanged( int state )
{
  IsKeep2dActive = ( state > 0 );
  qDebug() << "Will " << ( (IsKeep2dActive )?"": " not" ) << " keep the 2d Plot when drawing next time..." << endl;
}

/*!
  Exchange x and y axis in Likelihood
*/
void CmbMainWindow::flipAxis() {
  qDebug() << "FLIPAXIS" << endl;
  int i = mPlotControlDock->LikeXBox->value();
  mPlotControlDock->LikeXBox->setValue(mPlotControlDock->LikeYBox->value());
  mPlotControlDock->LikeYBox->setValue(i);
  QString f = mPlotControlDock->LikeMinX->text();
  mPlotControlDock->LikeMinX->setText(mPlotControlDock->LikeMinY->text());
  mPlotControlDock->LikeMinY->setText(f);
  f = mPlotControlDock->LikeMaxY->text();
  mPlotControlDock->LikeMaxY->setText(mPlotControlDock->LikeMaxX->text());
  mPlotControlDock->LikeMaxX->setText(f);
  startDrawLikeli_2d();
}


void CmbMainWindow::bigTabPageChanged(QWidget *w) {
  if (mPlotControlDock->BigTab->currentIndex() == 2)
    mParameterDock->ParameterTab->setCurrentIndex(2);
}

void CmbMainWindow::startReadWMAP() {
  allocateWatchThread(ReadWMAPWorker);
}

void CmbMainWindow::readWMAP() {
  try {
    TotalSteps[ReadWMAPWorker] = 100;
    StatusMsg[ReadWMAPWorker] = "Converting to binary:: Temperature covariance matrix";
    analyze-> initWMAPCommonTT(&Step[ReadWMAPWorker]);
    StatusMsg[ReadWMAPWorker] = "Converting to binary:: Polarization covariance matrix";
    analyze-> initWMAPCommonTE(&Step[ReadWMAPWorker]);
    analyze->WMAPNotYetInitialized=false;
  } catch (Bad_Error x) {
    CriticalHead[ReadWMAPWorker] = "Likelihood Error";
    CriticalMsg[ReadWMAPWorker] = x.s.c_str();
  }
}


void CmbMainWindow::endReadWMAP() { 
  qDebug() << "end READWMAP" << endl;
  //  qDebug() << "enddrawlikeli" << endl;
  calcAction->setEnabled(true);
  stopAction->setEnabled(false);
  drawLikeliAction->setEnabled(true);
  mPlotControlDock->LikeliDrawButton->setEnabled(true);
 
  statusBar()->showMessage("Finished reading of WMAP data");
 
}

QString CmbMainWindow::getDir(QString n) {
  QString dir;
  int pos = n.lastIndexOf('/');
  if (pos == -1) return QString::null;
  return n.left(pos+1);
}

QString CmbMainWindow::getMccFileName()
{
  QString fileName;
  ifstream in;
  do {
    fileName = QFileDialog::getOpenFileName(this, "Open Likelihood file", LikeDir,"Likelihood data (*.mcc)");
    if (fileName.isEmpty())
      break;
    qDebug() << "LOADING FILENAME: " << fileName << endl;
    in.open(fileName.toLatin1());
    if (fileName.endsWith(".mcc"))
      fileName.chop(4);
    LikeDir = getDir(fileName);
    if (!in.good()) {
      QMessageBox::warning(this,"Montecarlo  files should end in .mcc",
            "Please choose an existing Montecarlo file,\ndenoted by the suffix .mcc");
    }
  } while (!in.good());
  return fileName;
}

void CmbMainWindow::loadLikeli(const QString& mccFileName)
{
  QString fileName = mccFileName;
  if (mccFileName.isEmpty())
    fileName = getMccFileName();
  if (fileName.isEmpty())
    return;
  if (fileName.endsWith(".mcc"))
      fileName.chop(4);
  mPlotControlDock->LikeliBox->insertItem(0, fileName);
  mPlotControlDock->LikeliBox_1d->insertItem(0, fileName);
  mPlotControlDock->LikeliBox->setCurrentIndex(0);
  mPlotControlDock->LikeliBox_1d->setCurrentIndex(0);
  LikeliFiles.push_front(fileName);
  updateParameterNames(fileName);
}

QMap<int, ColumnPlotInfo> readParameterInfos(CmbMainWindow* w, const QString& paramNameFile)
{
  QString fileName=paramNameFile;
  if (fileName.endsWith(".mcc")) {
    fileName.chop(4);
  }
  QString nameFileName = QFileInfo(fileName).path()+"/automatic.parameterNames";
  qDebug() << "checking file:" << nameFileName << "for parameter names";
  QMap<int, ColumnPlotInfo> nameMap;
  if (QFile::exists(nameFileName)) {
    nameMap = readParameterNameFile(nameFileName);
  }
  return nameMap;
}


void CmbMainWindow::updateParameterNames(const QString& text)
{
  QMap<int, ColumnPlotInfo> nameMap = readParameterInfos(this, text);
  QList<ParameterSpinBox*> l = findChildren<ParameterSpinBox*>();
  foreach (ParameterSpinBox* s, l) {
    s->setParameterNames(nameMap);
  }
}


void CmbMainWindow::startDistill() {
  QStringList files = QFileDialog::getOpenFileNames(this, "Open several raw likelihood files", LikeDir, "Likelihood raw data (*.dat)");
  if (files.empty())
    return;
  QString proposedName = files.front();
  int suffixLength=QFileInfo(proposedName).suffix().length();
  proposedName.chop(suffixLength);
  proposedName+="mcc";
  QString name = QFileDialog::getSaveFileName(this, "Save Likelihood file as...", proposedName,"Likelihood data (*.mcc)");
  if (name.isEmpty())
    return;

  startDistill(files, name);
}

void CmbMainWindow::startDistill(const QStringList& files, const QString& distilledName)
{
  DistillSaveName = distilledName;

  if (DistillSaveName.endsWith(".mcc"))
    DistillSaveName.chop(4);

  LikeDir = getDir(DistillSaveName);

  QStringList::ConstIterator it = files.begin();
  DistillNames.resize(files.size());
  int count =0;
  while( it != files.end() ) {
    QString s = *it;
    DistillNames[count++] = s.toStdString();
    ++it;
  }
  allocateWatchThread(DistillWorker);
}

void CmbMainWindow::distill() {
  try {
    StatusMsg[DistillWorker] = "Distilling";
    // new version, drop distill obsolete
    double average = ChainShop::distillChain(DistillNames, DistillSaveName.toStdString(), &Step[DistillWorker],true,false);
    qDebug() << "Average stay at model: " << average << endl;
    mPlotControlDock->LikeliBox->insertItem(0,DistillSaveName);
    mPlotControlDock->LikeliBox->setCurrentIndex(0);
    mPlotControlDock->LikeliBox_1d->insertItem(0,DistillSaveName);
    mPlotControlDock->LikeliBox_1d->setCurrentIndex(0);
    LikeliFiles.push_front(DistillSaveName);
    emit summaryDone();
  } catch (Bad_Error x) {
    CriticalHead[DistillWorker] = "Distill Error";
    CriticalMsg[DistillWorker] = x.s.c_str();
  }
}

void CmbMainWindow::printLikeli() {
  postscriptFile = "";

  if (!prt->outputFileName().isEmpty()) {
    QString command = "cp  /tmp/likeli_nix.eps " + prt->outputFileName();
    system(command.toLatin1());
    return; // stop here, if only to file
  }

  QString PrinterName = prt->printerName();
  QString command = "eps2eps /tmp/likeli_nix.eps /tmp/likeli_nix_distilled.eps";
  system(command.toLatin1());

  command =  "lpr -P"+PrinterName + " /tmp/likeli_nix_distilled.eps";
  system(command.toLatin1());
  return;
}

void CmbMainWindow::startDrawLikeli() {
  int likepos;
  QString likeBoxText;
  if (LikelihoodDimension_is1) {
   likepos = mPlotControlDock->LikeLikeBox_1d->value();
   likeBoxText=mPlotControlDock->LikeLikeBox_1d->text();
  } else {
   likepos = mPlotControlDock->LikeLikeBox->value();
   likeBoxText=mPlotControlDock->LikeLikeBox->text();
  }
  if (!likeBoxText.isEmpty() && likeBoxText!=QLatin1String("TotalLogLike")) {
    int ret = QMessageBox::warning(this, "Likelihood column has suspicous name",
        "The column in the chainfile used to read in the likelihood of a point is named:\n"
        + likeBoxText +"\n" + "Do you want to continue anyways?",
        QMessageBox::Yes | QMessageBox::Abort);
    if (ret==QMessageBox::Abort) {
      return;
    }
  }
  allocateWatchThread(RasterizeWorker);
  drawLikeliAction->setEnabled(false);
  mPlotControlDock->LikeliDrawButton->setEnabled(false);
  statusBar()->showMessage("Rasterizing");
}

void CmbMainWindow::startDrawLikeli_1d() {
  LikelihoodDimension_is1 = true;
  startDrawLikeli();
}

#ifndef PRERELEASE
static ChainShop::DistributionType debugDistributionType;
static ChainShop::AveragingType debugAveragingType;
#endif

void CmbMainWindow::startDrawLikeli_2d() {
  if ( mPlotControlDock->MarkReferenceModelBox->isChecked() ) {
    QString dir = QFileInfo(currentMccFile()).absolutePath();
    QString fileName = "/modeltomark.txt";
    QFile nameFile(dir+fileName);
    if (!nameFile.open(QFile::ReadOnly)) {
      QMessageBox::warning(this, "Could not find reference model on disc.",
          (QLatin1String("could not open file: ")
           +dir+fileName+"\n Disabling option 'mark reference model'"));
      mPlotControlDock->MarkReferenceModelBox->setChecked(false);
    }
  }

#ifndef PRERELEASE
  if(!mPlotControlDock->Analytic->isChecked()) {
    QMessageBox msgBox, msgBox1;
    QPushButton *likeButton = msgBox.addButton("LnLike", QMessageBox::ActionRole);
    QPushButton *weightButton = msgBox.addButton("Weight", QMessageBox::ActionRole);
    msgBox.setText("Weight or Chi2?");
    msgBox.exec();
    if(msgBox.clickedButton()==likeButton)
      debugDistributionType=ChainShop::UseChi2;
    else if(msgBox.clickedButton()==weightButton)
      debugDistributionType=ChainShop::UseWeight;
    else
      return;
    QPushButton *maxButton, *totalButton, *avgButton;
    maxButton = msgBox1.addButton("max", QMessageBox::ActionRole);
    totalButton = msgBox1.addButton("total", QMessageBox::ActionRole);
    avgButton = msgBox1.addButton("average", QMessageBox::ActionRole);
    msgBox1.exec();
    if(msgBox1.clickedButton()==maxButton)
      debugAveragingType=ChainShop::UseMaximum;
    else if(msgBox1.clickedButton()==totalButton)
      debugAveragingType=ChainShop::UseTotal;
    else if(msgBox1.clickedButton()==avgButton)
      debugAveragingType=ChainShop::UseAverage;
    else
      return;
    } else {
      debugDistributionType=ChainShop::UseWeight;
      debugAveragingType=ChainShop::UseTotal;
    }
#endif

  LikelihoodDimension_is1 = false;
  startDrawLikeli();
}


void CmbMainWindow::endDrawLikeli() {
  calcAction->setEnabled(true);
  stopAction->setEnabled(false);
  drawLikeliAction->setEnabled(true);
  flipAxisAction->setEnabled(true);
  mPlotControlDock->LikeliDrawButton->setEnabled(true);
  if (LikelihoodDimension_is1)
    cw->PlotTab->setCurrentIndex(4); else
      cw->PlotTab->setCurrentIndex(3);
  mPlotControlDock->LikeText_1d->setHtml(Like1Text);

  statusBar()->showMessage("Finished Rasterizing");

  emit(repaintPlot());
  emit(likelihoodDrawingFinished());
}

vector<double> CmbMainWindow::adjustAxis(double start, double stop) {
  
  double dx = stop-start;
  dx /= 10.0;  // approximate order of ticks

  double order = rint( log(dx)/log(10.0) ); // order of magnitude of range
  
  double startorder = log(start)/log(10.0);
  double stoporder = log(stop)/log(10.0);

  //  pair<double,double> result = make_pair<double,double>(start,stop);
  vector<double> result(3);
  result[0] = start;
  result[1] = stop;
  //result[2] = pow(10, floor(stop-start)-1);

  order = rint(order);
  stoporder = floor(stoporder);
  startorder = ceil(startorder);
  if (start !=0 ) {
    double t = start * pow(10, -order);
    t = (double)(int)(t);
    t *= pow(10,order);
    result[0] = t;
  }
  if (stop !=0 ) {
    double t = stop * pow(10, -order);
    t = (double)(int)(t);
    t *= pow(10,order);
    result[1] = t;
  }
  result[2] = pow(10,order);
  
  while (fabs((stop-start)/result[2]) > 10) { 
    double r;
    r =  result[2] * 2.0;
    if (fabs((stop-start)/r) < 10) { result[2] = r; break;}
    r =  result[2] * 3.0;
    if (fabs((stop-start)/r) < 10) { result[2] = r; break;}
    result [2]*= 5.0;
  }
 return result;
}




void CmbMainWindow::drawLikeli() {
  qDebug() << "drawLikeli())" << endl;
  if (LikelihoodDimension_is1) {
    drawLikeli_1d();  // it's not the 2-D, but the little one :-) 
    return;
  }

  // we will change what the plotwidget should draw
  // but are called in a different thread than the main gui
  // thread where the plotwidget does its drawing
  // so this is necessary to avoid crashes
  cw->likeliplot->setUpdatesEnabled(false);
  // and make sure we re-enable them when we're done
  struct OnScopeExit
  {
    CmbEasyWidget* cW;
    OnScopeExit(CmbEasyWidget* w): cW(w) {}
    ~OnScopeExit() { cW->likeliplot->setUpdatesEnabled(true); }
  } ReEnableDrawing(cw);

  TotalSteps[RasterizeWorker] = 100;
  Step[RasterizeWorker] = 0;
  try {
    float xstart = mPlotControlDock->LikeMinX->text().toFloat();
    float xstop = mPlotControlDock->LikeMaxX->text().toFloat();
    int xpos = mPlotControlDock->LikeXBox->value();
    int likepos = mPlotControlDock->LikeLikeBox->value();
    QString mccfile = mPlotControlDock->LikeliBox->currentText()+".mcc";
    float smear = mPlotControlDock->LikeSmear->text().toFloat();

    ChainShop::Grid matrix;

    int ypos = mPlotControlDock->LikeYBox->value(); 
    float ystart = mPlotControlDock->LikeMinY->text().toFloat();
    float ystop = mPlotControlDock->LikeMaxY->text().toFloat();
    int steps = mPlotControlDock->LikeResolution->value();
    bool grid = mPlotControlDock->CreateGridFile->isChecked();
    bool Analytic = mPlotControlDock->Analytic->isChecked();
    if (Analytic) {
      StatusMsg[RasterizeWorker] = "Binning data points" ;
      smear = 0.2; // almost no smear for analytic determination, only binning
      steps = 30;
    }  else StatusMsg[RasterizeWorker] = "Smearing out the discrete data points" ;

    qDebug() << "BEFORE 2DIM " << endl;
    ChainShop shop;
    ChainShop::GridInfo gridInfo;
#ifndef PRERELEASE
    shop.setAveragingType(debugAveragingType);
    shop.setDistributionType(debugDistributionType);
#endif
    shop.get2DimMCDistribution(mccfile.toStdString(), matrix,
        xpos, ypos, likepos,
        gridInfo,
        steps, smear,
        &Step[RasterizeWorker], grid);
    mModelInspectorDock->allData()->setGrid(gridInfo);
    mModelInspectorDock->selectedData()->setGrid(gridInfo);
    vector<vector<float> > v;
    ChainShop::readMCC(mccfile.toStdString(), v);
    mModelInspectorDock->allData()->setPoints(v);
    mModelInspectorDock->selectedData()->setPoints(v);
    QMap<int, ColumnPlotInfo> ci = readParameterInfos(this, mccfile);
    mModelInspectorDock->allData()->setParameterInfos(ci);
    mModelInspectorDock->selectedData()->setParameterInfos(ci);
    mModelInspectorDock->showCols(xpos, ypos);
    mModelInspectorDock->showAllModels->setChecked(true);

    qDebug() << "DONE 2DIM" << endl;

    StatusMsg[RasterizeWorker] = "Calculating confidence regions"; 
    Step[RasterizeWorker] = 0;
    cw->likeliplot->Co = CoordRect(xstart,xstop,ystart,ystop);
  
    //  RasterizeReturn *r = analyze->rasterize(matrix, out, prt,true,&Step[RasterizeWorker]) ;
    // cw->likeliplot->setBlocks(r);  UNCOMMENT FOR RASTERIZE

    if (Analytic) {
      int persigma = mPlotControlDock->PerSigmaSpin->value();
      float maxsigma = mPlotControlDock->SigmaSpin->value();
      bool shade =true;
      if (persigma == 0) { persigma =1; shade = false;}  // If PerSigma =0 , then we don't want shading
      ChainShop::ConfidenceInference inference = ChainShop::Bayesian;
      if (mPlotControlDock->DeltaChi2->isChecked()) inference = ChainShop::DeltaChi2;
      int order = mPlotControlDock->Order2d->value();
      int quality = mPlotControlDock->Quality->value();
      quality = max(quality,1);
      quality = min(quality,9);
      list<ConfidenceRegion> *regions = ChainShop::getConfidenceRegions(matrix,maxsigma,persigma,shade,inference,order,quality, &Step[RasterizeWorker] );
      cw->likeliplot->addRegionGroup(regions);


      cw->likeliplot->setBlocks(0);
    } else {
      //      RasterizeReturn *r = analyze->rasterize(matrix, out, prt,true,&Step[RasterizeWorker]) ;
      list<Block>* b;
      bool old = false;
//X       qDebug() << "old or new? "; cin >> old;
#ifndef PRERELEASE
      if (debugDistributionType==ChainShop::UseChi2)
        old =true;
#endif
      if (old) {
        qDebug() << endl << "using old one." << endl;
        b = shop.rasterize(matrix,&Step[RasterizeWorker]) ;
      } else {
        qDebug() << endl << "using new one." << endl;
        b = shop.rasterize_Bayesian2(matrix, gridInfo, &Step[RasterizeWorker]) ;
      }
      cw->likeliplot->setBlocks(b);  
      cw->likeliplot->clearRegionGroups();
    }

    if (mPlotControlDock->ThirdBox->isChecked()) {
      ColorPoints cp;
      int zpos = mPlotControlDock->LikeZBox->value(); 
      ChainShop::getThirdDimension(mccfile.toStdString(),xpos,ypos,zpos,cp);
      cp.MinColor = Color(0,0,1);
      cp.MiddleColor = Color(0,1,0);
      cp.MaxColor = Color(1,0,0);
      cp.Radius = 6e-3;
      cp.Width = 0.03;
      cp.Height = 0.3;
      cp.LeftLower = CoordPoint(0.96,0.04);
      cw->likeliplot->attachColorPoints(cp);
    }

    if ( mPlotControlDock->BestFitBox->isChecked() ) {
      cw->likeliplot->setDrawBestFitModel( true );
      vector<float> bestModel = bestFitModel(mccfile,likepos);
      cw->likeliplot->BestFitModel = pair<float,float>( bestModel[ xpos ], bestModel[ ypos ] );
    } else {
      cw->likeliplot->setDrawBestFitModel( false );
    }
    if ( mPlotControlDock->MarkReferenceModelBox->isChecked() ) {
      cw->likeliplot->setDrawModelToMark( true );
      vector<float> model = modelToMark();
      cw->likeliplot->ModelToMark = pair<float,float>( model[ xpos ], model[ ypos ] );
    } else {
      cw->likeliplot->setDrawModelToMark( false );
    }
  } catch(Bad_Error x) {
    CriticalHead[RasterizeWorker] = "Likelihood Error";
    CriticalMsg[RasterizeWorker] = x.s.c_str();
  }
}

vector<float> CmbMainWindow::bestFitModel(QString mccfile, int likelipos) {
  vector<float> v = ChainShop::bestFitModel(mccfile.toLatin1().data(),likelipos); 
  qDebug() << endl<<endl;
  qDebug() << "================================" << endl;
  qDebug() << " Best Fit Model" << endl;
  qDebug() << "================================" << endl;
  for (unsigned int i = 0; i < v.size(); i ++) {
    qDebug() << "column[" << i << "]  :  " << v[i] << endl;
  }
  qDebug() << "================================" << endl;
  return v;
}

vector<float> CmbMainWindow::modelToMark()
{
  vector<float> v;

  QString dir = QFileInfo(currentMccFile()).absolutePath();
  QString fileName = "/modeltomark.txt";
  unsigned int lastCol = 0;
  QFile nameFile(dir+fileName);

  QTextStream nameStream(&nameFile);
  while (!nameStream.atEnd())   {
    unsigned int colNr;
    float parameterValue;
    nameStream >> colNr >> parameterValue;

    qDebug() << "read " << parameterValue << "at" << colNr;
    qDebug() << "for lastCol:" << lastCol;
    while(lastCol<colNr) {
      lastCol++;
      qDebug() << "inserting 0 for: " << lastCol;
      v.push_back(0.);
    }

    qDebug() << "inserting " << parameterValue << "at lastCol" <<lastCol;
    v.push_back(parameterValue);
    lastCol++;
  }

  unsigned int maxCol = ChainShop::numberOfColumns(currentMccFile().toStdString());

  while (lastCol < maxCol) {
    lastCol++;
    v.push_back(0.);
  }

  qDebug() << endl<<endl;
  qDebug() << "================================" << endl;
  qDebug() << " Model to mark in plot" << endl;
  qDebug() << "================================" << endl;
  for (unsigned int i = 0; i < v.size(); i ++) {
    qDebug() << "column[" << i << "]  :  " << v[i] << endl;
  }
  qDebug() << "================================" << endl;

  return v;
}


void CmbMainWindow::drawLikeli_1d() {
  ChainShop shop;
  TotalSteps[RasterizeWorker] = 100;
  Step[RasterizeWorker] = 0;
  try {
    if (mPlotControlDock->LikeResolution_1d->value() < 8) throw Bad_Error("CmbMainWindow::drawLikeli_1d() 8 bins is the minimum number");

    float xstart = mPlotControlDock->LikeMinX_1d->text().toFloat();
    float xstop = mPlotControlDock->LikeMaxX_1d->text().toFloat();
    int xpos = mPlotControlDock->LikeXBox_1d->value();
    int likepos = mPlotControlDock->LikeLikeBox_1d->value();
    QString mccfile = mPlotControlDock->LikeliBox_1d->currentText()+".mcc";

    map<double, vector<double> > distribution;

    shop.getMCDistribution(mccfile.toStdString(), distribution,xpos, xstart, xstop,
                             mPlotControlDock->LikeResolution_1d->value(),&Step[RasterizeWorker]);
    Spline s;
    mLastLikelihoodInfos.clear();

    //    s.set(analyze->fitExpPoly(distribution));    

    int order=mPlotControlDock->Order1d->value();
    CurrentExpPoly = ChainShop::fitExpPoly(distribution,order);
    map<double,double> current = ChainShop::evaluateExpPoly(CurrentExpPoly,xstart,xstop,20);
    if (IsLockedIn1d) {
      map<double,double> locked = ChainShop::evaluateExpPoly(LockedExpPoly,xstart,xstop,20);
      for (map<double,double>::iterator i=locked.begin(); i != locked.end(); i++) {
	qDebug() << "checking: " << i->first << "  sizes: " << current.size() << "  :: " << locked.size() << endl;
	if (current.find(i->first) == current.end()) throw Bad_Error("Strange, Current and Locked do not share same data");
	current[i->first] *= i->second;  
      }

      Spline tmp;
      tmp.set(current); // temporarily use one
      tmp.arm();

      double peakvalue = tmp(tmp.maximum());
      tmp *= 1.0/peakvalue;

      // replace the rectangles
      for (map<double,vector<double> >::iterator i =distribution.begin(); i != distribution.end(); i++) {
	double x = i->first;
	qDebug() << "left position is: " << x;
	x += 0.5*i->second[0];
	qDebug() << " plus the step width, we have: " << x << endl;
	i->second[1] = tmp(x);
      } 
      for (map<double,double>::iterator i=current.begin(); i != current.end(); i++) i->second /= peakvalue;
    }

    s.set(current);



    int id = cw->likeliplot1d->setSpline(s);
    cw->likeliplot1d->ridOf(id-1, false /*don't repaint*/);


    Color FillColorSigma1(157.0/255.0,179.0/255.0,  234.0/255.0);
    Color FillColorSigma2(221.0/255.0,244.0/255.0,  215.0/255.0);
    Color FillColorSigma3(249.0/255.0,246.0/255.0,  237.0/255.0);


    //Miscmath misc;

    double Threshold1 = 1.0;// misc.DeltaChi2CorrespondingToSigma(1.0,1.0);
    double Threshold2 = 4.0; //misc.DeltaChi2CorrespondingToSigma(1.0,2.0);

    qDebug() << "1-D Likelihood chi2 corresponding to 1 and 2 sigma: " << Threshold1 << " and " << Threshold2 << endl;

    Threshold1 = exp(-0.5*Threshold1);
    Threshold2 = exp(-0.5*Threshold2);
 
    list<Rectangle> *rects  = new list<Rectangle>; 
    for (map<double,vector<double> >::iterator i =distribution.begin(); i != distribution.end(); i++) {
      double like = i->second[1];
      Rectangle r(i->first,i->first + i->second[0],0,like);
      r.FillColor = FillColorSigma3;
      if (like > Threshold2) r.FillColor = FillColorSigma2;
      if (like > Threshold1) r.FillColor = FillColorSigma1;
      rects->push_back(r);

      qDebug()  << "Pushing rectangle: x= " << i->first << " X: " << i->first + i->second[0] << "  height: " << like << endl;

    }

    if ( mPlotControlDock->BestFitBox->isChecked() ) {
      QString name = mPlotControlDock->LikeliBox_1d->currentText();
      vector<float> model = bestFitModel(name+QLatin1String(".mcc"),likepos);
      int indexToMark = mPlotControlDock->LikeXBox_1d->value();
      double valueToMark = model[indexToMark];
      double stepWidth = distribution.begin()->second[0]/6.;
      qDebug() << "stepWidth: " << stepWidth << endl;
      Rectangle r(valueToMark-stepWidth,valueToMark+stepWidth,0., 1.);
      r.FillColor = Color(0.9,0.1,0.0);
      rects->push_back(r);
    }
    if ( mPlotControlDock->MarkReferenceModelBox->isChecked() ) {
      vector<float> model = modelToMark();
      if (!model.empty()) {
        int indexToMark = mPlotControlDock->LikeXBox_1d->value();
        double valueToMark = model[indexToMark];
        double stepWidth = distribution.begin()->second[0]/6.;
        qDebug() << "stepWidth: " << stepWidth << endl;
        Rectangle r(valueToMark-stepWidth,valueToMark+stepWidth,0., 1.);
        r.FillColor = Color(0.0,0.1,0.9);
        rects->push_back(r);
      }
    }


    // Re-scale the plotwidget
    cw->likeliplot1d->setRectangles(rects);
    cw->likeliplot1d->Co =  CoordRect(xstart,xstop,0,1.05);

    // Calculate maximum and sigmas in two ways:
    // The first one is a delta-chi2 test
    // The second one is meassuring the area under the curve

    // Delta Chi2 -test
    s.arm();
    //s.dump("sdump");
    Threshold1 *= s(s.maximum()); // re-normalize to fitted curve absolute normalisation
    vector<double> OneSigma = s.getZeroVector(Threshold1);
    vector<double> PlusMinus;

    if (OneSigma.size() == 2) {
      PlusMinus.resize(2);
      PlusMinus[0] = s.maximum() - OneSigma[0];
      PlusMinus[1] = OneSigma[1] - s.maximum();
    }

    // Area under curve
    Spline sint(1000,"sint");
    double sum=0;
    for (int i =0; i < s.size()-1; i++) {
      sint.set(s.x(i),sum); 
      sum += s.integrate(s.x(i),s.x(i+1));
    }
    sint.set(s.stop(),sum);
    sint.arm();
    sint.dump("sint");


    QString t(InformLike1d);
    LikelihoodInfo chi2Info( LikelihoodInfo::Chi2 );
    t+= "<table>";
    t += "<tr><td>Max.&nbsp;likl.&nbsp;&nbsp;at:&nbsp;</td><td><font color=#aa0000><b> " + toStr(s.maximum(),4) + "</b>"; 
    chi2Info.setMaximum( s.maximum() );
    if (PlusMinus.size()>1) {
      t += "&nbsp;+&nbsp;" + toStr(PlusMinus[1],4)+ "&nbsp;-&nbsp;" + toStr(PlusMinus[0],4);
      chi2Info.setLeft( PlusMinus[ 0 ] ); chi2Info.setRight( PlusMinus[ 1 ] );
    }
    mLastLikelihoodInfos.push_back( chi2Info );
    qDebug() << "deb1" << endl;
    t += "</font></td>";
    t += "<td><font color=#00aa00><font>delta&nbsp;chi<sup>2</sup></font></td></tr>\n";

    double Area = 68;
    for (;;) {
      double divide = 0.5*s(s.maximum());  // put it on half 
      double step= divide*0.5;
      double left=1e100;
      double right=-1e100;
      for (;;) {
        vector<double> cut = s.getZeroVector(divide);
        qDebug() << "CUTSIZE: "<< cut.size() << endl;
        //      if (cut.size() < 2) break;
        if (cut.size() < 1) break;
        if (cut.size() > 1) {
          left = cut[0];
          right = cut[cut.size()-1];
        }
        if (cut.size() == 1) {
          left = cut[0];
          right = left;
          if (left > s.maximum()) left = s.start();
          if (right < s.maximum()) right = s.stop();
        }
        double area = sint(right)-sint(left); 

        qDebug() << "left: "<< left << " right: " << right << "  area: " << area << " sum: " << sum << endl;
        if (area < 0.01*Area * sum) divide -= step; else divide += step;
        step *=0.5;
        if (step < 1e-3) break;
      }

      LikelihoodInfo likeInfo( LikelihoodInfo::None );
      if ( Area == 68 ) likeInfo.type = LikelihoodInfo::UnderCurve68;
      if ( Area == 95 ) likeInfo.type = LikelihoodInfo::UnderCurve95;
      if ( likeInfo.type == LikelihoodInfo::None ) throw Bad_Error( "drawlikeli_1d: couldn't determine likelihood type." );

      if (left != s.start() && right != s.stop()) {
        t += "<tr><td>Max.&nbsp;likl.&nbsp;&nbsp;at:&nbsp;</td><td><font color=#aa0000><b> " + toStr(s.maximum(),4) + "</b>"; 
        t += "&nbsp;+&nbsp;" + toStr(right-s.maximum(),4)+ "&nbsp;-&nbsp;" + toStr(s.maximum()-left,4);
        likeInfo.setMaximum( s.maximum() );
        likeInfo.setLeft( s.maximum()-left );
        likeInfo.setRight( right-s.maximum() );
      } else {
        if (left == s.start()) {
          t += "<tr><td>Constraint to: </td><td><font color=#aa0000><b>&#60;&nbsp;" +  toStr(right,4) + "</b>";
          likeInfo.setRight( right );
        } else 
          t += "<tr><td>Constraint to: </td><td><font color=#aa0000><b>&#62;&nbsp;" +  toStr(left,4) + "</b>";
        likeInfo.setLeft( left );
      }
      qDebug() << "deb2" << endl;
      t += "</font></td>\n";
      t += "<td><font color=#00aa00>"+toStr(Area,2) +"%&nbsp;area&nbsp;under&nbsp;curve</font></td></tr>\n";
      mLastLikelihoodInfos.push_back( likeInfo );
      if (Area == 95) break;
      if (Area == 68) Area=95;
    }
    t += "</table>";
    Like1Text =t; // will be displayed shortly. For some reasons this crashed in newer versions of qt - sure, because it might cause a repaint from a non gui-thread :)
    //    cw->LikeText_1d->setText("hallo");
    qDebug() << "Text is: \n" << t << endl ;
    qDebug() << "deb3" << endl;
  } catch(Bad_Error x) {
    if ( mIgnoreBadErrors ) return;
    qDebug() << "caught " << QString::fromStdString(x.s) << endl;
//X     QMessageBox::critical(this, "", QString::fromStdString("caught " + x.s), QMessageBox::Ok);
    CriticalHead[RasterizeWorker] = "Likelihood Error";
    CriticalMsg[RasterizeWorker] = x.s.c_str();
  }
  qDebug() << "~drawlikeli_1d" << endl;
}

void CmbMainWindow::toggledDims(bool state) {
  if (state) { 
    cw->PlotTab->setCurrentIndex(3); 
  } else {
    cw->PlotTab->setCurrentIndex(4);
  }
}


void CmbMainWindow::saveSettings() {
  qDebug() << "************** SAVING SETTINGS ************" << endl;
  QSettings s;
  s.setValue("/cmbeasy/background/omega_b",mParameterDock->ob->text());
  s.setValue("/cmbeasy/background/omega_bh2",mParameterDock->obh2->text());
  s.setValue("/cmbeasy/background/omega_cdm",mParameterDock->ocdm->text());
  s.setValue("/cmbeasy/background/omega_cdmh2",mParameterDock->ocdmh2->text());
  s.setValue("/cmbeasy/background/hubble",mParameterDock->hubble->text());
  s.setValue("/cmbeasy/background/olambda",mParameterDock->olambda->text());
  s.setValue("/cmbeasy/background/onu",mParameterDock->onu->text());
  s.setValue("/cmbeasy/background/nnu",mParameterDock->nnu->text());
  s.setValue("/cmbeasy/background/reionize",mParameterDock->reionize->text());
 
  s.setValue("/cmbeasy/spectral/spectral",mParameterDock->spectral->text());
  s.setValue("/cmbeasy/spectral/spectralTensor",mParameterDock->spectralTensor->text());
  s.setValue("/cmbeasy/spectral/stRatio",mParameterDock->stRatio->text());

  s.setValue("/cmbeasy/basic/adjust",mParameterDock->AdjustLambda->isChecked());
  s.setValue("/cmbeasy/basic/scalar",mParameterDock->Scalar->isChecked());
  s.setValue("/cmbeasy/basic/power",mParameterDock->PowerBox->isChecked());
  //  s.setValue("/cmbeasy/basic/morphing",mParameterDock->Morphing->isChecked()); MORPHHERE
  s.setValue("/cmbeasy/basic/tensor",mParameterDock->Tensor->isChecked());
  s.setValue("/cmbeasy/basic/lensed",mParameterDock->Lensed->isChecked());

  s.setValue("/cmbeasy/dark/tune",mParameterDock->TuneQuint->isChecked());
  s.setValue("/cmbeasy/dark/type",QuintTypeGroup->checkedId());

  s.setValue("/cmbeasy/likeli/third",mPlotControlDock->ThirdBox->isChecked());
  s.setValue("/cmbeasy/likeli/xcol",mPlotControlDock->LikeXBox->value());
  s.setValue("/cmbeasy/likeli/ycol",mPlotControlDock->LikeYBox->value());
  s.setValue("/cmbeasy/likeli/zcol",mPlotControlDock->LikeZBox->value());
  s.setValue("/cmbeasy/likeli/likecol",mPlotControlDock->LikeLikeBox->value());
  s.setValue("/cmbeasy/likeli/minx",mPlotControlDock->LikeMinX->text());
  s.setValue("/cmbeasy/likeli/maxx",mPlotControlDock->LikeMaxX->text());
  s.setValue("/cmbeasy/likeli/miny",mPlotControlDock->LikeMinY->text());
  s.setValue("/cmbeasy/likeli/maxy",mPlotControlDock->LikeMaxY->text());
  s.setValue("/cmbeasy/likeli/file",mPlotControlDock->LikeliBox->currentText());
  s.setValue("/cmbeasy/likeli/resolution",mPlotControlDock->LikeResolution->value());
  s.setValue("/cmbeasy/likeli/smear",mPlotControlDock->LikeSmear->text());
  s.setValue("/cmbeasy/likeli/likedir",LikeDir);
  //s.setValue("/cmbeasy/likeli/labelx",mPlotControlDock->LikeLabelX->text());
  //s.setValue("/cmbeasy/likeli/labely",mPlotControlDock->LikeLabelY->text());
  //  s.setValue("/cmbeasy/likeli/papersize",mPlotControlDock->PaperSize->currentItem());

  s.setValue("/cmbeasy/likeli/persigma",mPlotControlDock->PerSigmaSpin->value());
  s.setValue("/cmbeasy/likeli/maxsigma",mPlotControlDock->SigmaSpin->value());
  s.setValue("/cmbeasy/likeli/order2d",mPlotControlDock->Order2d->value());
  s.setValue("/cmbeasy/likeli/markbestfit", mPlotControlDock->BestFitBox->isChecked());
  s.setValue("/cmbeasy/likeli/markreference", mPlotControlDock->MarkReferenceModelBox->isChecked());

  s.setValue("/cmbeasy/likeli1d/bins",mPlotControlDock->LikeResolution_1d->value());
  s.setValue("/cmbeasy/likeli1d/xcol",mPlotControlDock->LikeXBox_1d->value());
  s.setValue("/cmbeasy/likeli1d/likecol",mPlotControlDock->LikeLikeBox_1d->value());
  s.setValue("/cmbeasy/likeli1d/minx",mPlotControlDock->LikeMinX_1d->text());
  s.setValue("/cmbeasy/likeli1d/maxx",mPlotControlDock->LikeMaxX_1d->text());
  s.setValue("/cmbeasy/likeli1d/derived",mPlotControlDock->Derived_1d->isChecked());
  s.setValue("/cmbeasy/likeli1d/order1d",mPlotControlDock->Order1d->value());

  s.setValue("/cmbeasy/plot/plotdir",PlotDir);
  s.setValue("/cmbeasy/printing/profiles/number", (int) PrtProfiles.size());


  int files = mPlotControlDock->LikeliBox->count();
  s.setValue("/cmbeasy/likeli/files",files);
 
  for (int i=0; i < files; i++) {
    QString num;
    s.setValue("/cmbeasy/likeli/filename" + num.setNum(i) , mPlotControlDock->LikeliBox->itemText(i) );  
  }


  for (unsigned int i = 0; i < PrtProfiles.size(); i++) {
    QString num, loc="/cmbeasy/printing/profiles/profile_";
    num.setNum(i);
    loc += num;
    //    s.setValue(loc + "/name",k->first); // map<name,printing>
    //    s.setValue(loc + "/scaleX",PrtProfiles[i].scaleX);
    //s.setValue(loc + "/scaleY",PrtProfiles[i].scaleY);
    s.setValue(loc + "/name",QString::fromStdString(PrtProfiles[i].Name)); // map<name,printing>
    s.setValue(loc +  "/xoffset",PrtProfiles[i].XLabelOffset);
    s.setValue(loc +  "/yoffset",PrtProfiles[i].YLabelOffset);
    s.setValue(loc + "/axislabelsize",PrtProfiles[i].AxisLabelSize);
    s.setValue(loc + "/ticklabelsize",PrtProfiles[i].TickLabelSize);
    s.setValue(loc + "/AutoTick",PrtProfiles[i].AutomaticTick);
    //    s.setValue(loc + "/ProtectedSize",PrtProfiles[i].ProtectedSize);


    s.setValue(loc + "/xlabeltext"   ,PrtProfiles[i].XLabelText.c_str());
    s.setValue(loc + "/ylabeltext"   ,PrtProfiles[i].YLabelText.c_str());
    s.setValue(loc + "/significantx"   ,PrtProfiles[i].Significant_x);
    s.setValue(loc + "/significanty"   ,PrtProfiles[i].Significant_y);
    s.setValue(loc + "/xstyle"  ,ConfigurationDialog::labelStyle2int(PrtProfiles[i].XTickLabelStyle));
    s.setValue(loc + "/ystyle"  ,ConfigurationDialog::labelStyle2int(PrtProfiles[i].YTickLabelStyle));
    s.setValue(loc +   "/leftmargin"  ,PrtProfiles[i].LeftMargin);
    s.setValue(loc +   "/rightmargin"  ,PrtProfiles[i].RightMargin);
    s.setValue(loc +   "/bottommargin"  ,PrtProfiles[i].BottomMargin       );
    s.setValue(loc +   "/topmargin"  ,PrtProfiles[i].TopMargin       );
    s.setValue(loc +   "/starttickx"  ,PrtProfiles[i].StartTick_x       );
    s.setValue(loc +   "/startticky"  ,PrtProfiles[i].StartTick_y       );
    s.setValue(loc +   "/steptickx"  ,PrtProfiles[i].StepTick_x       );
    s.setValue(loc +   "/stepticky"  ,PrtProfiles[i].StepTick_y       );
    s.setValue(loc + "/framelinewidth",PrtProfiles[i].FrameLineWidth);
    s.setValue(loc + "/curvelinewidth",PrtProfiles[i].CurveLineWidth);
    s.setValue(loc + "/tickslinewidth",PrtProfiles[i].TicksLineWidth);

    s.setValue(loc + "/custompaperx",PrtProfiles[i].Width);
    s.setValue(loc + "/custompapery",PrtProfiles[i].Height);

  }

  s.setValue("/cmbeasy/plot/data/wmap",mPlotControlDock->WMAPBox->isChecked());
  s.setValue("/cmbeasy/plot/data/boomerang",mPlotControlDock->BoomBox->isChecked());
  s.setValue("/cmbeasy/plot/data/tegmark",mPlotControlDock->TegBox->isChecked());

  s.setValue("/cmbeasy/tip/nr",(int)TipNr);
  s.setValue("/cmbeasy/tip/again",TipAgain);
}

void CmbMainWindow::loadSettings() {
  QSettings s;
  mParameterDock->ob->setText(s.value("/cmbeasy/background/omega_b","0.05").toString());
  mParameterDock->obh2->setText(s.value("/cmbeasy/background/omega_bh2","0.021125").toString());
  mParameterDock->ocdm->setText(s.value("/cmbeasy/background/omega_cdm","0.35").toString());
  mParameterDock->ocdmh2->setText(s.value("/cmbeasy/background/omega_cdmh2","0.147875").toString());
  mParameterDock->hubble->setText(s.value("/cmbeasy/background/hubble","0.65").toString());
  mParameterDock->olambda->setText(s.value("/cmbeasy/background/olambda","0.6").toString());
  mParameterDock->onu->setText(s.value("/cmbeasy/background/onu","0").toString());
  mParameterDock->nnu->setText(s.value("/cmbeasy/background/nnu","0").toString());
  mParameterDock->reionize->setText(s.value("/cmbeasy/background/reionize","0").toString());

  mParameterDock->spectral->setText(s.value("/cmbeasy/spectral/spectral","0.97").toString());
  mParameterDock->spectralTensor->setText(s.value("/cmbeasy/spectral/spectralTensor","-0.03").toString());
  mParameterDock->stRatio->setText(s.value("/cmbeasy/spectral/stRatio","0").toString());

  mParameterDock->AdjustLambda->setChecked(s.value("/cmbeasy/basic/adjust",true).toBool());
  mParameterDock->Scalar->setChecked(s.value("/cmbeasy/basic/scalar",true).toBool());
  mParameterDock->Tensor->setChecked(s.value("/cmbeasy/basic/tensor",false).toBool());
  mParameterDock->PowerBox->setChecked(s.value("/cmbeasy/basic/power",false).toBool());
  //  mParameterDock->Morphing->setChecked(s.value("/cmbeasy/basic/morphing",false)); // MORPHHERE
  mParameterDock->Lensed->setChecked(s.value("/cmbeasy/basic/lensed",false).toBool());

  mParameterDock->TuneQuint->setChecked(s.value("/cmbeasy/dark/tune",true).toBool());
  QuintTypeGroup->button(s.value("/cmbeasy/dark/type",0).toInt())->setChecked(true);


  mPlotControlDock->ThirdBox->setChecked(s.value("/cmbeasy/likeli/third",false).toBool());
  mPlotControlDock->LikeXBox->setValue(s.value("/cmbeasy/likeli/xcol",0).toInt());
  mPlotControlDock->LikeYBox->setValue(s.value("/cmbeasy/likeli/ycol",0).toInt());
  mPlotControlDock->LikeZBox->setValue(s.value("/cmbeasy/likeli/zcol",0).toInt());
  mPlotControlDock->LikeLikeBox->setValue(s.value("/cmbeasy/likeli/likecol",6).toInt());

  mPlotControlDock->LikeMinX->setText(s.value("/cmbeasy/likeli/minx","0.0").toString());
  mPlotControlDock->LikeMaxX->setText(s.value("/cmbeasy/likeli/maxx","0.0").toString());
  mPlotControlDock->LikeMinY->setText(s.value("/cmbeasy/likeli/miny","0.0").toString());
  mPlotControlDock->LikeMaxY->setText(s.value("/cmbeasy/likeli/maxy","0.0").toString());
  mPlotControlDock->BestFitBox->setChecked(s.value("/cmbeasy/likeli/markbestfit", true).toBool());
  mPlotControlDock->MarkReferenceModelBox->setChecked(s.value("/cmbeasy/likeli/markreference", true).toBool());

  int files = s.value( "/cmbeasy/likeli/files",0).toInt();
  files = min(files,5);

  for (int i=0; i < files; i++) {
    QString key = QLatin1String("/cmbeasy/likeli/filename") + QString::number(i);
    QString temp = s.value(key,"").toString();
    mPlotControlDock->LikeliBox_1d->insertItem(i,temp);
    mPlotControlDock->LikeliBox->insertItem(i,temp);
  }

 
  mPlotControlDock->LikeResolution->setValue(s.value("/cmbeasy/likeli/resolution",100).toInt());
  mPlotControlDock->LikeSmear->setText(s.value("/cmbeasy/likeli/smear","2").toString());

  mPlotControlDock->PerSigmaSpin->setValue(s.value("/cmbeasy/likeli/persigma",3).toInt());
  mPlotControlDock->SigmaSpin->setValue(s.value("/cmbeasy/likeli/maxsigma",2).toInt());

  mPlotControlDock->Order2d->setValue(s.value("/cmbeasy/likeli/order2d",4).toInt());

  LikeDir = s.value("/cmbeasy/likeli/likedir",QString()).toString();
  PlotDir = s.value("/cmbeasy/plot/plotdir",QString()).toString();
  //  mPlotControlDock->PaperSize->setCurrentItem(s.readNumEntry("/cmbeasy/likeli/papersize",0));

  mPlotControlDock->LikeResolution_1d->setValue(s.value("/cmbeasy/likeli1d/bins",12).toInt());
  mPlotControlDock->LikeXBox_1d->setValue(s.value("/cmbeasy/likeli1d/xcol",0).toInt());
  mPlotControlDock->LikeMinX_1d->setText(s.value("/cmbeasy/likeli1d/minx","0").toString());
  mPlotControlDock->LikeMaxX_1d->setText(s.value("/cmbeasy/likeli1d/maxx","0").toString());
  mPlotControlDock->Derived_1d->setChecked(s.value("/cmbeasy/likeli1d/derived",false).toBool());
  mPlotControlDock->LikeLikeBox_1d->setValue(s.value("/cmbeasy/likeli1d/likecol",6).toInt());
  mPlotControlDock->Order1d->setValue(s.value("/cmbeasy/likeli1d/order1d",4).toInt());



  //  mPlotControlDock->LikeLabelX->setText(s.value("/cmbeasy/likeli/labelx",QString()).toString());
  // mPlotControlDock->LikeLabelY->setText(s.value("/cmbeasy/likeli/labely",QString()).toString());

  mPlotControlDock->WMAPBox->setChecked(s.value("/cmbeasy/plot/data/wmap", true).toBool());
  mPlotControlDock->BoomBox->setChecked(s.value("/cmbeasy/plot/data/boomerang", false).toBool());
  mPlotControlDock->TegBox->setChecked(s.value("/cmbeasy/plot/data/tegmark", false).toBool());

  TipNr = s.value("/cmbeasy/tip/nr",0).toUInt();
  TipAgain = s.value("/cmbeasy/tip/again",true).toBool();


  int psize = s.value("/cmbeasy/printing/profiles/number",0).toInt();
  if (psize >= 7) {
    qDebug() << "NUMBER OF PROFILES: " << psize << endl;
    PrtProfiles.resize(psize);
    for (int i = 0; i < psize; i++,k++) {
      QString num, loc="/cmbeasy/printing/profiles/profile_";
      num.setNum(i);
      loc += num;
      LowLevelPlot prt;
      //      prt.scaleX = s.value(loc + "/scaleX",0).toDouble();
      // prt.scaleY = s.value(loc + "/scaleY",0).toDouble();
      prt.Name = s.value(loc + "/name",QLatin1String("unnamed")).toString().toStdString();
      prt.XLabelOffset =  s.value(loc + "/xoffset",0.06).toDouble();
      prt.YLabelOffset =  s.value(loc + "/yoffset",0.11).toDouble();
      prt.AxisLabelSize = s.value(loc + "/axislabelsize",100).toInt();
      prt.TickLabelSize  = s.value(loc + "/ticklabelsize",100).toInt();
      
      prt.AutomaticTick = s.value(loc + "/AutoTick",true).toBool(); 

      prt.XLabelText = s.value(loc + "/xlabeltext","").toString().toStdString();
      prt.YLabelText = s.value(loc + "/ylabeltext","").toString().toStdString();

      prt.Significant_x = s.value(loc + "/significantx",2).toInt();
      prt.Significant_y = s.value(loc + "/significanty",2).toInt();

      prt.XTickLabelStyle = ConfigurationDialog::int2labelStyle(s.value( loc + "/xstyle",1).toInt());
      prt.YTickLabelStyle = ConfigurationDialog::int2labelStyle(s.value( loc + "/ystyle",1).toInt());
      
      prt.LeftMargin = s.value(loc + "/leftmargin",0.15).toDouble();
      prt.RightMargin = s.value(loc + "/rightmargin",0.5).toDouble();
      prt.BottomMargin = s.value(loc + "/bottommargin",0.15).toDouble();
      prt.TopMargin = s.value(loc + "/topmargin",0.05).toDouble();


      prt.StartTick_x = s.value(loc +   "/starttickx",0).toDouble();
      prt.StartTick_y = s.value(loc +   "/startticky",0).toDouble();
      
      prt.StepTick_x = s.value(loc +   "/steptickx",0).toDouble();
      prt.StepTick_y = s.value(loc +   "/stepticky",0).toDouble();

      prt.FrameLineWidth = s.value(loc + "/framelinewidth",100).toInt();
      prt.CurveLineWidth = s.value(loc + "/curvelinewidth",100).toInt();
      prt.TicksLineWidth = s.value(loc + "/tickslinewidth",100).toInt();

      prt.Width = s.value(loc + "/custompaperx",90).toDouble();
      prt.Height = s.value(loc + "/custompapery",100).toDouble();
  
      //      prt.ProtectedSize =  s.value(loc + "/ProtectedSize",true).toBool(); 
      PrtProfiles[i] = prt;
      //qDebug() << "I DID READ: " << name << endl;
    }
  } else if ( psize < 7 ){ // updating from something very old...
    PrtProfiles.resize(7);
    for (int i = 0; i < 7; i++) {
      LowLevelPlot prt;
      PrtProfiles[i]=prt;
    }
  }

  for (int i = 0; i < 7; i++) { //make sure the standard profiles have the correct names
    string n;
    switch (i) {
    case 0: n = "Scalar <TT> plot"; break;
    case 1: n = "Tensor <TT> plot"; break;
    case 2: n = "Cold Dark Matter powerspectrum plot "; break;
    case 3: n = "2-D likelihood contour plot"; break;
    case 4: n = "1-D likelihood contour plot"; break;
    case 5: n = "Scalar <EE> plot"; break;
    case 6: n = "Scalar <TE> plot"; break;
    default: n = "Scalar <TT> plot";
    }
    PrtProfiles[i].Name = n;
  }
  emit prtProfilesChanged();
  qDebug() << "PRTPROFILESIZE: " << PrtProfiles.size() << endl;
}

void CmbMainWindow::quintTypeChanged(int id) {
  vector<QLabel*> label(3);
  label[0] = mParameterDock->QParamText1;
  label[1] = mParameterDock->QParamText2;
  label[2] = mParameterDock->QParamText3;
  vector<QLineEdit*> edit(3);
  edit[0] = mParameterDock->QuintParam1;
  edit[1] = mParameterDock->QuintParam2;
  edit[2] = mParameterDock->QuintParam3;

  if (id > -1) {
    switch (id) {
    case 1:
      cosmos->setQuintessence(Quintessence::leaping);
      break;
    case 2:
      cosmos->setQuintessence(Quintessence::ipl);
      break;
    default:
      cosmos->setQuintessence(Quintessence::none);
    }
  }
 
  // next, we  force cosmos to create the quintessence
  // the value of omega_q is not important, it will never be
  // used in a calculation cause this is always set before the actual calcul.
  cosmos->setOmega_quintessence(1e-3);

  vector<QPName> names = cosmos->quintessence()->parameterNames();
  vector<double> param =  cosmos->quintessence()->parameters();
  int given = names.size();  // store it

  if (names.size() < 3) names.resize(3);
  

  for (int i=0; i < 3; i++) {
    QPName n = names[i];
    if (n.name == "") n.name = "void";
    label[i]->setText((n.name + ":").c_str());
    if (n.tooltip != "") {
      edit[i]->setToolTip(QString::fromStdString(n.tooltip));
    }
    if (n.whatsthis != "") {
      edit[i]->setWhatsThis(QString::fromStdString(n.whatsthis));
    }
    // if it is to be tuned and it is determined and it is really a parameter (i < given), 
    // we disable the parameter
    if (i < given) 
      qDebug() << "i < given " << n.determined << "   " <<  mParameterDock->TuneQuint->isChecked() << endl;


    edit[i]->setText(toStr(param[i],3));

    if (i < given) {
      edit[i]->setHidden(n.determined && mParameterDock->TuneQuint->isChecked()); 
      label[i]->setHidden(n.determined && mParameterDock->TuneQuint->isChecked()); 
    }
    else {
      edit[i]->setHidden(true);
      label[i]->setHidden(true);
    }
  }
}

void CmbMainWindow::tuneQuintChanged() { quintTypeChanged(-1); }

QString CmbMainWindow::toStr(double x, int post) {
  LowLevelPlot low;
  string s = low.double2string(x,post);
  return s.c_str();
}
/*
  QString tmp;
  tmp.setNum(x);
  
  int k = tmp.indexOf('.');
  int e = tmp.indexOf('e');

  if (k == -1 && e == -1) return tmp;   // not . not e
  if (e == -1) return tmp.left(k+1+post);  //  .  but no e
 
  // both
  
  return tmp.left(k+1+post) + tmp.right(tmp.length()-e);
}
*/

void CmbMainWindow::showTip() {
  if (TipAgain) {
    if (TipD) delete TipD;
    TipD = new TipDialog(&TipNr,&TipAgain,this);
    TipD->show();
  }
}

void CmbMainWindow::showTipNow() {
  TipAgain = true;
  showTip();
}

/*! 
  After loading tjhe settings, this method updates whatever
  necessary
*/
void CmbMainWindow::syncWithSettings() {
  quintTypeChanged(QuintTypeGroup->checkedId());
  /*
for (map<QString,Printing>::iterator i = PrtProfiles.begin(); i != PrtProfiles.end(); i++) {
    cw->PaperSize->insertItem(i->first);
    cw->PaperSize_1d->insertItem(i->first);
  }
  */
}

void CmbMainWindow::revertZoom() {
  PlotWidget *p = currentPlotWidget();
  if (p) p->revertZoom();
  p->repaint();
}

/*!
  Connected to the scalar and tensor checkbox in basic
  control. If both are checked then enable Tensor / scalar etc
*/
/*
void CmbMainWindow::scalarOrTensorToggle() {
  qDebug() << "SCALAR TENSOR TOGGLE" << endl;
  if (cw->Scalar->isChecked() && cw->Tensor->isChecked()) qDebug() << "should be enabled" << endl;
  //cw->TensorScalarFrame->setEnabled(cw->Scalar->isChecked() && cw->Tensor->isChecked());
  cw->stgiven->setEnabled(cw->Scalar->isChecked() && cw->Tensor->isChecked());
  cw->stRatio->setEnabled(cw->Scalar->isChecked() && cw->Tensor->isChecked());
}
*/

void CmbMainWindow::stop() { 
  interrupt=true; 
  calcAction->setEnabled(true);
  stopAction->setEnabled(false);
  KWorker->wait();
  RasterizeWorker->wait();
}

void CmbMainWindow::updateModelInspector(const QPointF& p)
{
  mModelInspectorDock->selectedData()->updateDataSelection(mPlotControlDock->LikeXBox->value(),
                                                           mPlotControlDock->LikeYBox->value(),
                                                           p.x(), p.y());
  mModelInspectorDock->reSortSelectedData();
}

