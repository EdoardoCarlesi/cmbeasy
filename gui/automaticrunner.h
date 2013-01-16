#ifndef AUTOMATICRUNNER_H
#define AUTOMATICRUNNER_H

#include "cmbmainwindow.h"

#include <QObject>

class AutomaticRunner: public QObject
{
  Q_OBJECT
  public:
    AutomaticRunner(CmbMainWindow *w, const QStringList& args);

    static bool runningNonInteractively;

  public slots:
    void run();
    void distill();

  private:
    CmbMainWindow *mMainWindow;
    enum { TaskDistill, TaskStatistics, TaskPlots1d, TaskPlots2d } mTask;
    QString mTargetDir;
    QString mTargetFile;
    QString mParameterNamesFile;
    QStringList mMccDataFiles;
    QString mLikeliFile;
};

#endif // AUTOMATICRUNNER_H
