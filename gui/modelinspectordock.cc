#include "modelinspectordock.h"

ModelInspectorWidget::ModelInspectorWidget(QWidget *parent)
    : QWidget(parent), mSelectedDataCount(0), mAllDataCount(0)
{
  setupUi(this);

  connect(showAllModels, SIGNAL(toggled(bool)), this, SLOT(showAllData(bool)));
  connect(selectedModelsOnly, SIGNAL(toggled(bool)), this, SLOT(showCurrentData(bool)));
  connect(&mAllDataModel, SIGNAL(numberShownChanged(unsigned int)), this, SLOT(setAllDataCount(unsigned int)));
  connect(&mSelectedDataModel, SIGNAL(numberShownChanged(unsigned int)), this, SLOT(setSelectedDataCount(unsigned int)));
  connect(&mAllDataModel, SIGNAL(layoutChanged()), this, SLOT(updateHeaders()));
  connect(&mSelectedDataModel, SIGNAL(layoutChanged()), this, SLOT(updateHeaders()));

  setModelsShown();

  modelTable->setSortingEnabled(true);
  connect(modelTable->horizontalHeader(), SIGNAL(sortIndicatorChanged(int,Qt::SortOrder)),
          this, SLOT(sortByColumn(int, Qt::SortOrder)));
  modelTable->horizontalHeader()->setClickable(true);
  modelTable->horizontalHeader()->setContextMenuPolicy(Qt::ActionsContextMenu);

  mSignalMapper = new QSignalMapper(this);
}

void ModelInspectorWidget::updateHeaders()
{
  const QMap<int, ColumnPlotInfo> *pI=0;
  if(showAllModels->isChecked()) {
    pI = &mAllDataModel.paramInfos();
  } else if(selectedModelsOnly->isChecked()) {
    pI = &mSelectedDataModel.paramInfos();
  }
  if(!pI)
    return;

  delete mSignalMapper;
  mSignalMapper = new QSignalMapper(this);
  foreach(QAction *a, modelTable->horizontalHeader()->actions()) {
    delete a;
  }
  QMapIterator<int, ColumnPlotInfo> it(*pI);
  while(it.hasNext()) {
    it.next();
    QAction *action = new QAction(it.value().columnLabel, this);
    action->setText(it.value().columnLabel);
    action->setCheckable(true);
    modelTable->horizontalHeader()->addAction(action);
    mSignalMapper->setMapping(action, it.key());
    bool shown = false;
    if(mHeaderShown.contains(it.key())) {
      shown = mHeaderShown.value(it.key());
    }
    //modelTable->horizontalHeader()->setSectionHidden(it.key(), !shown);
    modelTable->setColumnHidden(it.key(), !shown);
    action->setChecked(shown);
    connect(action, SIGNAL(toggled(bool)), mSignalMapper, SLOT(map()));
  }
  connect(mSignalMapper, SIGNAL(mapped(int)),
            this, SLOT(setHeaderShown(int)));
  qDebug() << "updated headers: "<< mHeaderShown;

  modelTable->reset();
}
