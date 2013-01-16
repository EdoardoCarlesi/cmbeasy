#include "cmbeasy.h"
#include "controlpanel.h"
#include "cmbmainwindow.h"
#include "automaticrunner.h"

#include "gsl/gsl_errno.h"

#include <QApplication>
#include <QFont>
#include <QObject>
#include <QPixmap>
#include <QAction>
#include <QPaintEvent>
#include <QSplashScreen>
#include <QDesktopWidget>
#include <QPainter>
#include <QFont>
#include <QtDebug>

#include <unistd.h>

void cmbeasy_gui_gsl_error_handler( const char * reason, const char * file, int line, int gsl_errno )
{
  std::string errorString;
  errorString += "An error occured in the GSL library: ";
  if (file) errorString += file;
  errorString += ", line ";
  errorString += line;
  errorString += ".\nThe reason was: ";
  if (reason) errorString += reason;
  errorString += "\nGSL error no: ";
  errorString += QString::number(gsl_errno).toStdString();
  errorString += '\n';
  throw Bad_Error( errorString );
}

int main(int argc, char *argv[]) {

  gsl_set_error_handler ( cmbeasy_gui_gsl_error_handler );

  try {
  QApplication app(argc,argv);

  QStringList args = app.arguments();
  bool automaticRun = false;
  AutomaticRunner*  automaticRunner = 0;
  if (args.size() >= 2
      && (args.at(1)==QLatin1String("--generateStatistics")
          || args.at(1)==QLatin1String("--generate1dPlots")
          || args.at(1)==QLatin1String("--generate2dPlots")
          || args.at(1)==QLatin1String("--distillTo"))) {
    automaticRun = true;
  }

  QPixmap pm( CmbMainWindow::cmbeasyDir( "/resources/help/splash.png" ) );
  QSplashScreen *splash = 0;
  QString m("by Michael Doran,\n  Georg Robbers\n and Christian Mueller");
  if (!automaticRun) {
    splash = new QSplashScreen( pm );
    splash->show();
    splash->showMessage(m, Qt::AlignBottom | Qt::AlignRight, Qt::white);
    qApp->processEvents();
  }

  CmbMainWindow *w = new CmbMainWindow();
  if (automaticRun) {
    automaticRunner = new AutomaticRunner(w, args);
  }


  if (w->InitializationOk) {
    QObject::connect(&app,SIGNAL(aboutToQuit()),w,SLOT(saveSettings()));
    QObject::connect(&app,SIGNAL(lastWindowClosed()),&app,SLOT(quit()));
    QObject::connect(w->exitAction,SIGNAL(triggered()),&app,SLOT(quit()));
    w->loadSettings();
    w->syncWithSettings();
    w->show();


    if (!automaticRun) {
      sleep(1);
      w->showTip();
      splash->finish(w);
      delete splash;
    }
    return app.exec();
  }
  } catch (Bad_Error e) {
    qCritical("Bad_Error occured");
    qCritical(e.s.c_str());
  }
}
