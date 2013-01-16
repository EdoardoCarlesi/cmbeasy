#include "plotcontrol.h"

PlotControlDockWidget::PlotControlDockWidget(QWidget *parent)
  :QWidget(parent)
{
  setupUi(this);

  bool hidelock = true;
  Confidence_Inference->setHidden(true);  // "Chi2" is only for testing

#ifndef PRERELEASE
  hidelock=false;
  Confidence_Inference->setHidden(false);
  Confidence_Inference->setEnabled(true);
  DeltaChi2->setEnabled(true);
  Bayesian->setEnabled(true);
#endif

  Derived_1d->setHidden(true); // not needed

  LockIn1d->setHidden(hidelock);

  keep2dButton = new QCheckBox( buttonFrame2d );
  buttonFrame2d->layout()->addWidget( keep2dButton );
  keep2dButton->setText( "Lock 2d" );

#ifndef EXPERIMENTAL
  keep2dButton->setHidden( true );
#endif
}
