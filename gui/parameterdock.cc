#include "parameterdock.h"

#include "cmbmainwindow.h"

#include "quintcosmos.h"

#include <QDoubleValidator>

ParameterDockWidget::ParameterDockWidget(QuintCosmos& c, QWidget *parent)
               : QWidget(parent), cosmos(c)
{
  setupUi(this);


  QString pixDir = CmbMainWindow::cmbeasyDir("/pix/");
  Obh2Label->setPixmap( QPixmap( pixDir + QString::fromUtf8( "../pix/omega_bh2.png" ) ) );
  ObLabel->setPixmap( QPixmap( pixDir + QString::fromUtf8( "../pix/omega_b.png" ) ) );
  OmLabel->setPixmap( QPixmap( pixDir + QString::fromUtf8( "../pix/omega_cdm.png" ) ) );
  OlambdaLabel->setPixmap( QPixmap( pixDir + QString::fromUtf8( "../pix/omega_lambda.png" ) ) );
  OnuLabel->setPixmap( QPixmap( pixDir + QString::fromUtf8( "../pix/omega_nu.png" ) ) );

  connect(obh2,SIGNAL(returnPressed()),this,SLOT(obh2Return()));
  connect(ocdmh2,SIGNAL(returnPressed()),this,SLOT(ocdmh2Return()));
  connect(ob,SIGNAL(returnPressed()),this,SLOT(obReturn()));
  connect(ocdm,SIGNAL(returnPressed()),this,SLOT(ocdmReturn()));
  connect(hubble,SIGNAL(textChanged(const QString&)),this,SLOT(hubbleReturn()));

  connect(olambda,SIGNAL(textChanged(const QString&)),this,SLOT(lambdaReturn()));

  obh2->setValidator(new QDoubleValidator(0.0,1.0,12,obh2));
  ob->setValidator(new QDoubleValidator(0.0,1.0,12,ob));
  ocdmh2->setValidator(new QDoubleValidator(0,1,12,ocdmh2));
  ocdm->setValidator(new QDoubleValidator(0,1,12,ocdm));
  hubble->setValidator(new QDoubleValidator(0,1,12,hubble));
  spectral->setValidator(new QDoubleValidator(0.6,1.5,12,spectral));
  reionize->setValidator(new QDoubleValidator(0,1,12,reionize));

  ocdm->setText(toStr(cosmos.omega_cdm(),3));  
  ob->setText(toStr(cosmos.omega_b(),3));  
  olambda->setText(toStr(cosmos.omega_v(),3)) ;
  hubble->setText(toStr(cosmos.h(),3)) ;
  reionize->setText(toStr(cosmos.optDistanceLss(),3));
  spectral->setText(toStr(cosmos.InitialPower[0]));


  cosmos.setReionizationFraction(1);


  ocdmReturn();
  obReturn();
}


void ParameterDockWidget::obh2Return() { 
  double h = hubble->text().toDouble();
  double o = obh2->text().toDouble() / (h*h);
  //cosmos.setOmegaH2_b(obh2->text().toDouble());
  ob->setText(toStr(o ,3));  
  adjustLambda();
}

void ParameterDockWidget::ocdmh2Return() {
  double h = hubble->text().toDouble();
  double o = ocdmh2->text().toDouble() / (h*h);
  //cosmos.setOmegaH2_cdm(ocdmh2->text().toDouble());
  ocdm->setText(toStr(o,3));  
  adjustLambda();
}

void ParameterDockWidget::obReturn() { 
  double h = hubble->text().toDouble();
  double o = ob->text().toDouble() * (h*h);
  //cosmos.setOmega_b(ob->text().toDouble());
  obh2->setText(toStr(o,3));  
  adjustLambda();
}

void ParameterDockWidget::ocdmReturn() { 
  double h = hubble->text().toDouble();
  double o = ocdm->text().toDouble() * (h*h);
  //cosmos.setOmega_cdm(ocdm->text().toDouble());
  ocdmh2->setText(toStr(o,3));  
  adjustLambda();
}

void ParameterDockWidget::hubbleReturn() { 
  cosmos.seth(hubble->text().toDouble());
  ocdmReturn();
  obReturn(); // misuse these two such that obH2 etc are recalculated ... 
}


void ParameterDockWidget::lambdaReturn() {
  if (AdjustLambda->isChecked()) {
    double o = olambda->text().toDouble();
    cosmos.setOmega_vacuum(o);
    o += ob->text().toDouble();
    cosmos.setOmega_cdm(1.0 - o );
    ocdm->setText(toStr(1.0 - o,3));  
    ocdmh2->setText(toStr((1-o)*pow(hubble->text().toDouble(),2),3));     
  }
}

void ParameterDockWidget::adjustLambda() {
  if (AdjustLambda->isChecked()) {
    double o = ob->text().toDouble() + ocdm->text().toDouble();
    olambda->setText(toStr(1.0 - o,3)); 
  }
}

void ParameterDockWidget::fillCosmos() {
  //cout << "fill Cosmos" << endl;
  double lam = olambda->text().toDouble();
  double qunt = 0; 
  if (cosmos.quintessenceType() != Quintessence::none) {
    lam =0;
    qunt = olambda->text().toDouble();
  }
  cosmos.setOmega_quintessence(qunt);
  cosmos.setOmega_cdm(ocdm->text().toDouble());
  cosmos.setOmega_b(ob->text().toDouble());
  cosmos.setOmega_vacuum(lam);
  cosmos.seth(hubble->text().toDouble());
  cosmos.setOptDistanceLss( reionize->text().toDouble() );
  cosmos.InitialPower[0]=spectral->text().toDouble();
  cosmos.InitialTensorPower[0]=spectralTensor->text().toDouble();
  cosmos.setOmega_nuNR(onu->text().toDouble());
  cosmos.setNuNR(nnu->text().toDouble());


  // now fill quintessence parameters
  vector<QLineEdit*> edit(3);
  edit[0] = QuintParam1;
  edit[1] = QuintParam2;
  edit[2] = QuintParam3;

  vector<QPName> names = cosmos.quintessence()->parameterNames();
  vector<double> p = cosmos.quintessence()->parameters();
  for (unsigned int i = 0; i < names.size(); i++) {
    if (!names[i].determined ||   ! TuneQuint->isChecked()  ) {   // if tuning does not determine it or tuning is not on
      p[i] = edit[i]->text().toDouble();
    }
  }
  cosmos.setQParameters(p);
}


QString ParameterDockWidget::toStr(double x, int post) {
  //cout << "cmbeasywidget::tostr " << x << " post: " << post << endl;
  QString tmp;
  tmp.setNum(x);
  
  //cout << "Ceas1" << endl;
  int k = tmp.indexOf('.');
  int e = tmp.indexOf('e');
  //cout << "Ceas2 " << e << "  " << k <<  endl;
  if (k == -1 && e == -1) return tmp;   // not . not e
  if (e == -1) return tmp.left(k+1+post);  //  .  but no e
  //cout << "Ceas3" << endl;
  // both
  
  return tmp.left(k+1+post) + tmp.right(tmp.length()-e);
}


