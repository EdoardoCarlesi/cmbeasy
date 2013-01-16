#ifndef CONFIGURATIONDIALOG_H
#define CONFIGURATIONDIALOG_H

class QPushButton;

using namespace std;

#include <vector>
#include "ui_design_configurationdialog.h"
#include "analyzethis.h"
#include "lowlevelplot.h"
class ConfigurationDialog : public QDialog, private Ui::DesignConfigurationDialog {
  Q_OBJECT
  vector<QString> ConfigurationString;

  LowLevelPlot safety;
  vector<LowLevelPlot> &PrtProfiles;
  double Point;

 public:
  int CurrentProfile;
  ConfigurationDialog(vector<LowLevelPlot>& ,QWidget* parent = 0, const char* name = 0, bool modal = FALSE);
  void setProfile(int id);
  void syncDialog(LowLevelPlot&);
  void syncProfile();
  virtual void reject();
  virtual void accept();
  static LowLevelPlot::LabelStyle int2labelStyle(int);
  static int labelStyle2int(LowLevelPlot::LabelStyle);

 signals:
  void saveProfile();

 public slots:
  void updateProfiles();
  void deleteProfile();
  void activateProfile( QListWidgetItem* item );
};
#endif
