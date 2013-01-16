#ifndef CMBMAINWINDOW_H
#define CMBMAINWINDOW_H

#include "spline.h"
#include "safevector.h"
#include "cl.h"
#include "model.h"
#include "mathobject.h"
#include "data.h"
#include "analyzethis.h"
#include "lowlevelplot.h"
#include "chainshop.h"
#include "plotcontrol.h"
#include "parameterspinbox.h"

#include <QMainWindow>
#include <QString>
#include <QTime>
#include <QDockWidget>

class CmbEasyWidget;
class PlotWidget;
class ConfigurationDialog;
class ParameterDockWidget;
class ModelInspectorWidget;
class TipDialog;
class LikeliWorker;

class Data;
class ControlPanel;
class QuintCosmos;
class CmbCalc;
class AnalyzeThis;

class Q4ButtonGroup;

class QToolBar;
class QAction;
class QPopupMenu;
class QPoint;
class QProgressBar;
class QMutex;
class QPrinter;
class QProcess;
class QFileDialog;

#define index_type long

struct ExtendedInfo {
  double tau0, taueq, tauls;
  double zls,zeq;
  QString QuintType;
  double o_q;
  double wmap_tt;
  double wmap_te;
};

struct LikelihoodInfo {
  bool hasMaximum, hasLeft, hasRight;
  enum Type {None, Error, Chi2, UnderCurve68, UnderCurve95};
  Type type;
  LikelihoodInfo( Type t ): hasMaximum( false ), hasLeft( false ), hasRight( false ), type( t ) {}
  void setLeft( double l ) { mLeft = l; hasLeft = true; }
  void setRight( double r ) { mRight = r; hasRight = true; }
  void setMaximum( double m ) { mMaximum = m; hasMaximum = true; }
  double left() { if ( !hasLeft ) throw Bad_Error( "LikelihoodInfo: invalid access" ); return mLeft; }
  double right() { if ( !hasRight ) throw Bad_Error( "LikelihoodInfo: invalid access" ); return mRight; }
  double maximum() { if ( !hasMaximum ) throw Bad_Error( "LikelihoodInfo: invalid access" ); return mMaximum; }
  bool isEmpty() { return (!hasMaximum && !hasLeft && !hasRight); }
private:
  LikelihoodInfo() {}
  double mMaximum, mLeft, mRight;
};


/*
struct Parameters {
  Parameters(Model*);
  Parameters() {}
  QString toStr(double,int);
  float o_b, lambda, h, n,obh2,optdlss,o_qls;
  Parameters(QuintCosmos*);
  void printStatus(string="");
  QString statusString();
  QString printToLine();
  void readLine(ifstream &in); // set parameters from input stream
  double operator() (const int);
};
*/
class CmbMainWindow : public QMainWindow, Mathobject {
  Q_OBJECT
    private:
 
  enum CMBParam { baryons, cold, hubble, spectral,optical};
 
  CmbEasyWidget * cw;
  AnalyzeThis *analyze;
  bool LikelihoodDimension_is1; // true if 1-D, false if 2-D 

 public:
  QAction* calcAction;
  QAction* stopAction;
  QAction* plotInfoAction; //!< RMB plot entry
  QAction* plotCutAction; //!< RMB plot entry
  QAction* printAction; //!< print
  QAction* exitAction;
  QAction *drawLikeliAction; //!< draw likelihood
  QAction* cosmosAction; //!< like klicking on the cosmos button
  QAction* autoScaleAction; //!< autoscaling of plots
  QAction* plotPrtGnuplotAction; //!< for printing the plot via paint devices
  QAction* saveCurveAction; //!< save the curve to disk
  QAction* saveAllCurvesAction; //!< save all curves of one plot to disk
  QAction* flipAxisAction; //!< flip axis for likelihood
  QAction *autoAdaptAction; //!< given the columns for x and y find out about the ranges
  QAction* revertZoomAction; //!< revert Zoom
  QAction* PlotLayoutAction; //!< Detailed Plotting options for printing
  QAction* whatsThisAction; //!< enter "What's This?" mode
  QAction* toFrontAction; //!< bring Region to front
  QAction* toBackAction; //!< send Region to background
  QAction* deleteRegionAction; //!< delete Region
  QAction* colorRegionAction; //!<  assign different color to Region
  QAction* stickyRegionAction; //!<  assign different color to Region
  QAction* generateStatisticsFileAction; //!< generate a file with max and 1sigma for each column
  QAction* generatePlotTableAction; //!< automatically generate 1d plots for all columns and a tex file that includes them


 private:


  ParameterDockWidget*   mParameterDock;
  PlotControlDockWidget* mPlotControlDock;
  ModelInspectorWidget*  mModelInspectorDock;

  QDockWidget *mParamDockWidget;
  QDockWidget *mPlotControlDockWidget;
  QDockWidget *mModelInspectorDockWidget;

  QToolBar *tool;
  QMenu *PlotRMB; //!< for the plots if mouse is over a spline
  QMenu *Plot2ndRMB; //!< also for the plots but when no spline is selected
  QMenu *PlotLikeliRMB; //!< also for the plots but only when likelihood plot is the current page
  QMenu *PlotLikeliHoverRMB; //!< also for the plots but only when 2d likelihood plot is the current page and we hover over a region
  QMenu *FileMenu;
  QMenu *CurveaMenu;
  QMenu *PlotMenu;
  QMenu *CalculationMenu;
  QMenu *HelpMenu;
  QMenu *LikeliMenu;

  Q4ButtonGroup *QuintTypeGroup;

  TipDialog *TipD;
  unsigned int TipNr;
  bool TipAgain;
  ConfigurationDialog *cfg;
  vector<LowLevelPlot> PrtProfiles;

       
   map<LikeliWorker*,QProgressBar*> Progress;
   map<LikeliWorker*,QString> CriticalHead, CriticalMsg; //!< if this is none-zero, watchdog will display a critical msg
   map<LikeliWorker*,QString> StatusMsg,OldMsg;    //!< watchdog publishes this, if changed
   map<LikeliWorker*,int> Step, TotalSteps, OldStep,OldTotalSteps;
   map<LikeliWorker*,QTimer*> Timer;
   map<LikeliWorker*,bool> Active;
   
  QTime TimeStart; //!< for time information
  //  QMutex *mutex; //!<

  // RMB Plot Menu and stuff
  int RMBPlot_id; //<! id of the line 
  int RMBPlot_currentPage; //<! page on which RMB has been pressed

  QString PixDir; // Path to icon files
  //  QString ModelFile;
  // QString ModelInfFile;
  ofstream *debugstream;

  // CMB Stuff

  QuintCosmos *cosmos;
  ControlPanel *control;
  SafeVector<double> InitialPower; // replaces intialps 
  CmbCalc* cmbcalc;
  CL cl;
  vector<int> lval; // The l - values


  int QuintTypeId; //! For keeping track of the selected type of quintessence in the button group

  //vector<Model*> models;  //!< all models models
  map<index_type,Model*> cacheModel; //!< models already read in
  map<index_type,long> hashModel; //!< position of model with hash_key in file

  map<int, Model*> plotModel; //!< the models displayed as plots
  map<int,ExtendedInfo> extendedInfo;



  bool interrupt; // if true, then stop() has been pressed
 
  //LikeliWorker *Likeli; //!< A thread that calls our routines :-)
  LikeliWorker *KWorker;
  LikeliWorker *RasterizeWorker;
  LikeliWorker *ReadWMAPWorker;
  LikeliWorker *DistillWorker;
  bool fullPostscript; //!< evaluated by drawLikeli. If true it generates the plot as postscript to file
  string postscriptFile; //!< use this file to write the data to...
  QPrinter *prt;
  QString LikeDir,PlotDir;

  list<QString> LikeliFiles;
  vector<string> DistillNames; //!< a list of filenames to distill
  QString DistillSaveName; //!< The file name to save the distilled data to

  QString InformLike1d; //!< A string that is to be displayed always in the 1-dim panel
  QString Like1Text; //!< The total string to be displayed in the 1-dim panel
  ExpPoly LockedExpPoly; //!< If locked in, this is the ExpPoly
  ExpPoly CurrentExpPoly; //!< the Current ExpPoly
  bool IsLockedIn1d;
  bool IsKeep2dActive; //!< experimental/internal

  int k, nk;
  int nkTransfer; //!< only needed for progressbar. this is the k value from which on the transfers are calculated 
  int j,maxJ;

  vector<LikelihoodInfo> mLastLikelihoodInfos; //!< contains info about the calculated 1d likelihood (only the one that was calculated last, that is)
  vector< vector<LikelihoodInfo> > mLikelihoodInfos; //!< contains info about the calculated 1d likelihoods (for all columns that have been done so far)
  std::vector<ColumnPlotInfo> m1dStatisticsColumns; //!< keeps the info about the 1d statistics columns that should be plotted
  int mNumberOfColumns;
  bool mIgnoreBadErrors;


  // morphing variables
  /*
    Parameters mMin, mMax,mDelta;
    bool DatabaseExists;

  map<double,int> obmap, lambdamap,hmap,nmap,optdlssmap;
  
  Data* Sn1a, *Sn1aCorrList, *Calan;
  vector< vector< DataEntry> >* Sn1aCorr;

  bool addDeltaSuccess; //!< set to false if the closest or sandwich model are not in the database
  */
  list<Data*> measured; 
  bool noPingPong; //!< if true, resyncPlot() will immediatly return. This is used by syncPlot() which can safely setText() to PlotXMin etc
  
  int infoModel; //!< store the id of the model that is currently displayed in info (or -1) for nothing

  /*
  void readInf(bool rebuilt=false); //!< read model space and the position of each model in models.dat from file models.inf. 
  void readMap(ifstream&,map<double, int> &); //!< read in the maps of parameters and their ordering number 1,2,3,4...
  void writeMap(ofstream&,map<double, int> &); //!< write map of parameters to file (models.inf usually)
 
  Model* findCloseModel(Parameters&);
  Model* findSandwichModel(CMBParam);
  index_type modelIndex(Parameters&);
  index_type modelIndex(Model* m) { 
    Parameters p(m);
    return modelIndex(p);
  }
  bool inDataBaseRange(Parameters&); //!< if the parameters are within, i.e [..], not (..) the parameters of the model data base, return true.

  bool valid(Model*);
 

  Model* displayMorph(bool gfx=true);
  void addDelta(Spline*,CMBParam,bool cmb);  
  */

  double sigma8(Model*); //!< calculate sigma 8 from the spectrum of this model (for fun)
  double sigma8(Spline*,double); //!< as above, for fun
  double sig8Integrator(const double) const;
  Spline* tmpSig8; //!< needed for integration temporarily
  double tmpH; //!< needed for integration, this hubble h


  PlotWidget *idx2plot(int); //!< get PlotWidget for tab number int
  PlotWidget *currentPlotWidget(); //!< get PlotWidget of currently displayed tab
  int plot2idx(PlotWidget*p); //!< return index of plotwidget

  void initAction();
  void initCosmos(); 
  void initControl();
  void initMenu();

  QString mAutomaticPrintDirectory; //!< holds the directory where automatically generated plots should be placed
  std::vector<ColumnPlotInfo> mColumnsToPlot;
  bool mGenerating2dPlotTable;
  int  m2dPlotTableX;

 public:
  CmbMainWindow();
  bool InitializationOk; //!< if something in creator goes wrong, this will tell you
  void updateInfo(int id); //!< updates info for model with id 
  ExtendedInfo& fillExtendedInfo(ExtendedInfo&);

  void syncWithSettings();
  static QString cmbeasyDir(QString addon=""); //!< converts ControlPanel::cmbeasyDir() 's string to QString. This is for convenience.
  void allocateWatchThread(LikeliWorker*);

  PlotControlDockWidget* plotControlDock() const { return mPlotControlDock; }

  public slots:
    // Slots for calculation
    void calc();  //!< start cmb calculation
  //void intK(int);  //!< start integration of k-values
  void spool(); 
  //void intKTransfer(int,int); //!< start integration of higher k-values
  // void intCl(int); //!< start integration of Cl - spectrum
  //void callOneK(); //!< ask for one k- value to be calculated
  //void callOneKTransfer();  //!< same as callOneK() for higher k
  //void callOneCl(); //!< ask for one l- value to be calculated
  void endSpool();
  void watchdogCentral();
  void stop(); //! interrupt calculation
  void readData(); //!< read measured Data points to Data structure
  void plotPressedRMB(int, const QPoint&); //!< connected to the plotwidgets.
  void plotInfo(); //!< slot for RMB plot info request
  void plotCursorOverCurve(int); //!< called by plot if curve is touched by cursor
  void plotCut(); //!< slot for RMB plot cut request

  void print(); //!< called by printAction
  void printGnuplot(); //!< called by plotPrtQTAction
  void syncPlot(); //!< take current plot and synchronise the min_x etc in the BigTab
  void resyncPlot(); //!< the other direction, take values of BigTab and ask PlotWidget to redraw
  
  void showTip(); //!< Create and show tip dialog (if TipAgain is true)
  void showTipNow(); //!< Show the tip dialog, invoked from the help menu
  //  void scalarOrTensorToggle();

  void askCosmos(); //!< called by Cosmos button in info tab. Get Background evolution.

  void syncControl(); //!< set control.power_cdm to this value
  
  void setRMBPlotID(int id) { RMBPlot_id = id;}; //!< We (mis)-use the RMBPlot_id variable for ocasions like detected collision in a plot (to enable sensible shortcuts)

  static QString toStr(double,int post=0); //!< convert to string
 
  void saveCurve(); //!< save curve to disk
  void saveAllCurves(); //!< save all curves of one plot

  void autoScale(); //!< connected to the AutoScale button, calls autoScale() for the corresponding plotwidget
  
  void toFront(); //!< only in 2-d likelihood: bring regiongroup to front
  void toBack(); //!< only in 2-d likelihood: send regiongroup to back
  void deleteRegion(); //!< only in 2-d likelihood: delete RegionGroup
  void colorRegion(); //!< only in 2-d likelihood: assing new color to  RegionGroup
  void stickyRegion(); //!< only in 2-d likelihood: make region sticky
  

  void saveSettings();
  void loadSettings();
  void quintTypeChanged(int); //!< if user chooses other quinttype, this slot is called
  void tuneQuintChanged();
  void help();  //!< oline manual
  void about(); //!< about this
  void bigTabPageChanged(QWidget*); //!< used to check whenever the likelihood page comes up

  
  void startReadWMAP(); //!< setting up the thread for WMAP initialization
  void readWMAP(); //!< call AnalyzeThis::initMAPCommon() as a thread during startup
  void endReadWMAP(); //!< cleaning up after WMAP intitialization

  
  QString getDir(QString);
  vector <double> adjustAxis(double start, double stop);
  QString getMccFileName();
  void loadLikeli(const QString& mccFileName = QString());
  void updateParameterNames(const QString& fileName);
  void startDrawLikeli_2d();
  void startDrawLikeli_1d();
  void startDrawLikeli();
  void drawLikeli();
  void drawLikeli_1d();  
  void endDrawLikeli();
  void printLikeli();
  void distill(); //!< distill
  void startDistill(); //!< get information for distilling from user
  void startDistill(const QStringList& files, const QString& distilledName); //!< set up distilling of montecarlo files
  vector<float> bestFitModel(QString mccfile, int likelipos);
  vector<float> modelToMark(); //!< reads in parameters from a file that can be marked in a 2d-plot

  //! write out a file that has max and 1sigma for each column (1st part)
  void generateStatisticsFile(int col=0, QString dir=QString(), QString parameterNames=QString());
  void collectLikelihoodInfos(); //!< helper slot for generateStatisticsFile();
  void askForStatisticsFileFormat(); //!< ask the user in which format to use for the statistics file
  void outputStatisticsFile(); //!< write out a file that has max and 1sigma for each column (2nd part)
  void outputStatisticsTexFile(); //!< same as the outputStatisticsFile(), but formatted as a tex table

  void generatePlotTable(int col=0, QString dir=QString(), QString parameterNamesFile=QString(), int dimensions=0); //!< generate a tex file with a table that includes 1d plots for all parameters of a mc run
  void outputPlotForTable(); //!< write the graphics for a table of plots to disc
  void printAutomatic(); //!< print current plot without user intervention
  void outputPlotTable(); //!< output the tex table with all plots
  QString automaticPrintFileName(int colX=-1, int colY=-1); //!< automatically create a file name for plot output
  std::vector<ColumnPlotInfo> askForColumnsToPlot(); //!< show Dialog with checkboxes for which columns to include in summary
  //! read a file to decide which columns to include in summary infos
  std::vector<ColumnPlotInfo> columnPlotInfoFromParameterNamesFile(const QString& nameFileName);

  QString currentMccFile();
  int currentLikelihoodColumn();

  void flipAxis(); //!< flip axis of likelihood plot + recalc the likelihoods
  bool autoAdapt(); //!< given the columns for x and y of likelihood in .mcc file, auto find the ranges of the parameters
  bool autoAdapt_1d(); //!< version for 1d
  //! make sure that the values at xPos and yPos of model are within range
  void adjustRangeForModel(int xPos, int yPos, const std::vector<float>& model, std::vector<float>& range) const;
  bool onlyZeroesInColumn(const int col);  // checks wheter column col is all zeros by calling autoAdapt_1d()

  void advancedPlotSettings();

  void lockIn1d();
  void keep2dStateChanged(int);

 private slots:
   //  void readModels();  //!< read models, get model index and write them to file models.inf, then call readInf() 
  void revertZoom(); //!< if reverting of zoom for a plot is asked for
  void savePrtProfile();
  void toggledDims(bool);

  //!< update the model inspector widget, connected to left mouse button click on the plotwidget
  void updateModelInspector(const QPointF&);

 // void syncPaperSizeTo1(int);
 // void syncPaperSizeTo2(int);


 signals:
 void likelihoodDrawingFinished();
 void repaintPlot();
 void prtProfilesChanged();
 void summaryDone();

  // Slots for cw->handling
};

#endif
