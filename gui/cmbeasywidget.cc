#include "cmbeasywidget.h"
#include "plotwidget.h"
#include "cmbmainwindow.h"
#include "quintcosmos.h"

CmbEasyWidget::CmbEasyWidget( QuintCosmos& c,QWidget* parent, const char* name )
  : QWidget( parent ), DesignCmbEasyWidget()
{
  setupUi(this);
}

CmbEasyWidget::~CmbEasyWidget()
{
}

