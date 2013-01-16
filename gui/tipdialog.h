#ifndef TIPDIALOG_H
#define TIPDIALOG_H

using namespace std;

#include <vector>
#include "ui_design_tipdialog.h"

class TipDialog : public QDialog, private Ui::DesignTipDialog {
  Q_OBJECT
  vector<QString> TipString;
  
 public:
  TipDialog(unsigned int *tip, bool *again, QWidget * parent = 0, const char * name = 0, bool modal = FALSE);
  unsigned int *Tip;
  bool *ShowAgain;

  bool close();
  
  public slots:
  
    void  closeIt();
  void showTip();

  //signals:
  //  void closed(int);
};

#endif
