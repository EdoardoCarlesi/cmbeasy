#ifndef GUI_PLOTCONTROL_H
#define GUI_PLOTCONTROL_H

#include "ui_design_plotcontrol.h"

class PlotControlDockWidget: public QWidget, public Ui::PlotControlForm
{
  Q_OBJECT

  public:
    PlotControlDockWidget(QWidget *parent=0);

    QSize minimumSizeHint() const { return QSize(0,0); }

    QCheckBox* keep2d() { return keep2dButton; }

  private:
    QCheckBox *keep2dButton;

};
#endif // GUI_PLOTCONTROL_H
