#include "automaticrunner.h"

#include <QMessageBox>
#include <QTimer>
#include <QApplication>
#include <QDir>
#include <QDebug>
#include <QDockWidget>

bool AutomaticRunner::runningNonInteractively = false;

AutomaticRunner::AutomaticRunner(CmbMainWindow *w, const QStringList& args)
            : QObject(), mMainWindow(w)
{
  connect(mMainWindow, SIGNAL(summaryDone()), qApp, SLOT(closeAllWindows()));
  if (args.count() > 2 && args.at(1)==QLatin1String("--distillTo")) {
    mTask = TaskDistill;
    mTargetFile = args.at(2);
    mMccDataFiles = args.mid(3);
    QTimer::singleShot(0, this, SLOT(distill()));
    runningNonInteractively = true;
    return;
  }
  if (args.count() != 5) {
    QMessageBox::warning(w, "Wrong number of arguments",
             "The commandline for generating a summary should be of the form:\n cmbeasy <action> file targetDir parameternamesFiles");
    QApplication::exit(0);
  }

  if (args.at(1)==QLatin1String("--generateStatistics"))
    mTask = TaskStatistics;
  if ( args.at(1)==QLatin1String("--generate1dPlots"))
    mTask = TaskPlots1d;
  if (args.at(1)==QLatin1String("--generate2dPlots"))
    mTask = TaskPlots2d;

  mLikeliFile = args.at(2);

  mTargetDir = args.at(3);
  QDir target;
  if (target.exists(mTargetDir)) {
    QMessageBox::warning(w, "Target Directory exists",
             "The target directory for automatically genereated plots must not exist, it will be created automatically");
    QTimer::singleShot(0, qApp, SLOT(quit()));
    return;
  }
  if (!target.mkpath(mTargetDir)) {
    QMessageBox::warning(w, "Error creating directory:", "Could not create directory:\n"+mTargetDir);
    QTimer::singleShot(0, qApp, SLOT(quit()));
    return;
  }
  mParameterNamesFile = args.at(4);

  runningNonInteractively = true;
  mMainWindow->plotControlDock()->DeltaChi2->setChecked(true);
  QTimer::singleShot(0, this, SLOT(run()));
}

void AutomaticRunner::run()
{
  mMainWindow->loadLikeli(mLikeliFile);
  if (!mTargetDir.endsWith('/'))
    mTargetDir += '/';
  switch (mTask) {
   case TaskStatistics: mMainWindow->generateStatisticsFile(0, mTargetDir, mParameterNamesFile); break;
   case TaskPlots1d:    mMainWindow->generatePlotTable(0, mTargetDir, mParameterNamesFile, 1); break;
   case TaskPlots2d:    mMainWindow->generatePlotTable(0, mTargetDir, mParameterNamesFile, 2); break;
   default: qDebug() << "unknown task in AutomaticRunner::run()";
  }
}

void AutomaticRunner::distill()
{
  qDebug() << "distilling: " << mMccDataFiles << "target: " << mTargetFile;
  mMainWindow->startDistill(mMccDataFiles, mTargetFile);
}
