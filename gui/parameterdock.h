#ifndef GUI_PARAMETERDOCK_H
#define GUI_PARAMETERDOCK_H


#include "ui_design_parameterdock.h"

class QuintCosmos;


class ParameterDockWidget: public QWidget, public Ui::ParameterDockForm
{
  Q_OBJECT

  public:
    ParameterDockWidget(QuintCosmos& c, QWidget* parent=0);

    QSize minimumSizeHint() const { return QSize(0,0); }
    QSize sizeHint() const { return QSize(400,400); }

  public slots:
    void obh2Return();
    void ocdmh2Return();
    void obReturn();
    void ocdmReturn();
    void hubbleReturn();
    void lambdaReturn();

    void adjustLambda();
    void fillCosmos();
    QString toStr(double, int=0);
    //void olambdaReturn();

  private:
    QuintCosmos& cosmos;
};

#endif // GUI_PARAMETERDOCK_H
