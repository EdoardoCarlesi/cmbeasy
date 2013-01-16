#ifndef CMBEASYWIDGET_H
#define CMBEASYWIDGET_H

#include "ui_design_cmbeasywidget.h"

#include <QObject>

class QuintCosmos;

class CmbEasyWidget : public QWidget, public Ui::DesignCmbEasyWidget {
  Q_OBJECT

  public:
    CmbEasyWidget(QuintCosmos&,QWidget* parent = 0, const char* name = 0 );
    virtual ~CmbEasyWidget();
};
#endif // CMBEASYWIDGET_H
