#ifndef GUI_MODELINSPECTORDOCK_H
#define GUI_MODELINSPECTORDOCK_H


#include "ui_design_modelinspector.h"

#include "chainshop.h"
#include "parameterspinbox.h" // for ColumnPlotInfo

#include <QAbstractTableModel>
#include <QMutex>
#include <QSignalMapper>


#include <vector>
#include <algorithm>
#include <stdexcept>

typedef std::vector<std::vector<float> > ChainData;

struct VecVecCompareByCol
{
  int col;
  VecVecCompareByCol(int c): col(c) {}
  bool operator()(const std::vector<float>& v1, const std::vector<float>& v2) const {
                 return v1[col]<v2[col];
  }
};

class AllChainDataModel: public QAbstractTableModel
{
  Q_OBJECT

  public:
    virtual int rowCount(const QModelIndex & parent = QModelIndex() ) const {
      return mDataPoints.size();
    }

    virtual int columnCount ( const QModelIndex & parent = QModelIndex() ) const {
      return mDataPoints.begin()->size();
    }

    virtual QVariant data ( const QModelIndex & index, int role = Qt::DisplayRole ) const {
      QMutex *mMutableMutex = const_cast<QMutex*>(&mMutex);
      QMutexLocker locker(mMutableMutex);
      if (!index.isValid()
          || role != Qt::DisplayRole
          || index.row()>=mDataPoints.size()
          || index.column()>=mDataPoints.begin()->size())
        return QVariant();

      return mDataPoints[index.row()][index.column()];
    }

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const {
      if (role != Qt::DisplayRole
          || orientation != Qt::Horizontal
          || section>mParamInfos.size()) {
        return QAbstractTableModel::headerData(section, orientation, role);
      }
      return mParamInfos[section].columnLabel;
    }

    void setPoints(const ChainData& dataPoints) {
      mMutex.lock();
      mDataPoints=dataPoints;
      mMutex.unlock();
      reset();
      emit numberShownChanged(mDataPoints.size());
    }

    const QMap<int, ColumnPlotInfo>& paramInfos() const { return mParamInfos; }
    void setParameterInfos(const QMap<int, ColumnPlotInfo>& i) {
      mMutex.lock();
      mParamInfos=i;
      mMutex.unlock();
      reset();
    }

    void setGrid(const ChainShop::GridInfo& gi) { mGridInfo = gi; reset(); }

    virtual void sort(int column, Qt::SortOrder order=Qt::AscendingOrder) {
      switch(order) {
        case Qt::AscendingOrder:
               std::sort(mDataPoints.begin(), mDataPoints.end(),
                          VecVecCompareByCol(column));
               break;
        case Qt::DescendingOrder:
               std::sort(mDataPoints.rbegin(), mDataPoints.rend(),
                          VecVecCompareByCol(column));
               break;
        default:
               ; // doesn't happen
      }
      reset();
    }
  signals:
     void numberShownChanged(unsigned int);


  private:
    ChainData mDataPoints;
    ChainShop::GridInfo mGridInfo;
    QMap<int, ColumnPlotInfo> mParamInfos;
    QMutex mMutex;
};

class SelectedChainDataModel: public QAbstractTableModel
{
  Q_OBJECT

  public:
    virtual int rowCount(const QModelIndex & parent = QModelIndex() ) const {
      return mSelectedData.size();
    }

    virtual int columnCount ( const QModelIndex & parent = QModelIndex() ) const {
      if(mSelectedData.size()==0)
        return 0;
      return mSelectedData.begin()->size();
    }

    virtual QVariant data ( const QModelIndex & index, int role = Qt::DisplayRole ) const {
      QMutex *mMutableMutex = const_cast<QMutex*>(&mMutex);
      QMutexLocker locker(mMutableMutex);
      if (!index.isValid()
          || role != Qt::DisplayRole
          || index.row()>=mSelectedData.size()
          || index.column()>=mSelectedData.begin()->size())
        return QVariant();
      return mSelectedData[index.row()][index.column()];
    }

    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const {
      if (role != Qt::DisplayRole
          || orientation != Qt::Horizontal
          || section>mParamInfos.size()) {
        return QAbstractTableModel::headerData(section, orientation, role);
      }
      return mParamInfos[section].columnLabel;
    }

    void setPoints(const ChainData& dataPoints) {
      mMutex.lock();
      mDataPoints=dataPoints;
      mMutex.unlock();
      reset();
    }

    void setParameterInfos(const QMap<int, ColumnPlotInfo>& i) {
      mMutex.lock();
      mParamInfos=i;
      mMutex.unlock();
      reset();
    }

    void setGrid(const ChainShop::GridInfo& gi) { mGridInfo = gi; reset(); }

    struct IsOutsideGridBox
    {
      unsigned int _xBin, _yBin, _xCol, _yCol;
      const ChainShop::GridInfo& _g;

      IsOutsideGridBox(float x, float y, unsigned int xCol, unsigned int yCol,
                   const ChainShop::GridInfo& g)
        : _xCol(xCol), _yCol(yCol), _g(g) {
        _xBin = _g.x2bin(x);
        _yBin = _g.y2bin(y);
        //qDebug() << "bins are:" << _xBin << _yBin;
      }

      bool operator()(const vector<float>& p) {
         return (_xBin!=_g.x2bin(p[_xCol]) || _yBin!=_g.y2bin(p[_yCol]));
      }
    };

    void updateDataSelection(unsigned int xCol, unsigned int yCol, float x, float y) {
      emit layoutAboutToBeChanged();
      //QModelIndex topLeft = index(0,0);
      //QModelIndex bottomRight = index(mSelectedData.size(),mParamInfos.size());
      //int oldColumnCount=mSelectedData.size();
      mSelectedData.clear();
      remove_copy_if(mDataPoints.begin(), mDataPoints.end(),
                     back_inserter(mSelectedData),
                     IsOutsideGridBox(x, y, xCol, yCol, mGridInfo));
      emit numberShownChanged(mSelectedData.size());
      //emit dataChanged(topLeft, bottomRight);
      emit layoutChanged();
    }

    const QMap<int, ColumnPlotInfo>& paramInfos() const { return mParamInfos; }

    virtual void sort(int column, Qt::SortOrder order=Qt::AscendingOrder) {
      emit layoutAboutToBeChanged();
      qDebug() << "SelectedChainDataModel::sort()";
      switch(order) {
        case Qt::AscendingOrder:
               std::sort(mSelectedData.begin(), mSelectedData.end(),
                          VecVecCompareByCol(column));
               break;
        case Qt::DescendingOrder:
               std::sort(mSelectedData.rbegin(), mSelectedData.rend(),
                          VecVecCompareByCol(column));
               break;
        default:
               qWarning() << "SelectedChainDataModel::sort() - default case shouldn't happen";
               ; // doesn't happen
      }
      emit layoutChanged();
    }

  signals:
     void numberShownChanged(unsigned int);

  private:
    ChainData mDataPoints, mSelectedData;
    ChainShop::GridInfo mGridInfo;
    QMap<int, ColumnPlotInfo> mParamInfos;
    QMutex mMutex;
};


class ModelInspectorWidget: public QWidget, public Ui::ModelInspectorWidget
{
  Q_OBJECT

  public:
    ModelInspectorWidget(QWidget* parent=0);

    AllChainDataModel* allData() { return &mAllDataModel; }
    SelectedChainDataModel* selectedData() { return &mSelectedDataModel; }
    QAbstractTableModel* currentModel() {
      if(showAllModels->isChecked())
        return allData();
      else if(selectedModelsOnly->isChecked())
        return selectedData();
      return 0;
    }

    void reSortSelectedData() {
       if(modelTable->horizontalHeader()->isSortIndicatorShown()) {
         selectedData()->sort(modelTable->horizontalHeader()->sortIndicatorSection(),
                              modelTable->horizontalHeader()->sortIndicatorOrder());
       }
     }

   void showCols(int x, int y) { mHeaderShown[x] = mHeaderShown[y] = true; }

  public slots:
    void showAllData(bool enabled=true) {
       if(!enabled)
         return;
       if(modelTable->horizontalHeader()->isSortIndicatorShown()) {
         allData()->sort(modelTable->horizontalHeader()->sortIndicatorSection(),
                         modelTable->horizontalHeader()->sortIndicatorOrder());
       }
       QItemSelectionModel *m = modelTable->selectionModel();
       modelTable->setModel(allData());
       delete m;
       setModelsShown();
       updateHeaders();
     }

    void showCurrentData(bool enabled=true) {
       if(!enabled)
         return;
       reSortSelectedData();
       QItemSelectionModel *m = modelTable->selectionModel();
       modelTable->setModel(selectedData());
       delete m;
       setModelsShown();
       updateHeaders();
    }

    void setAllDataCount(unsigned int i) { mAllDataCount=i; setModelsShown(); }
    void setSelectedDataCount(unsigned int i) { mSelectedDataCount=i; setModelsShown(); }

    void updateHeaders();
    void setHeaderShown(int i) {
      mHeaderShown[i]=qobject_cast<QAction*>(mSignalMapper->mapping(i))->isChecked();
      modelTable->setColumnHidden(i, !mHeaderShown[i]);
      qDebug() << "setHeaderShown()" << mHeaderShown;
    }

  void sortByColumn (int column, Qt::SortOrder order) {
    qDebug() << "ModelInspectorWidget::sortByColumn()";
    currentModel()->sort(column, order);
  }

  protected:
    void setModelsShown() {
        unsigned int count=0;
      if(showAllModels->isChecked())
        count = mAllDataCount;
      else if(selectedModelsOnly->isChecked())
        count = mSelectedDataCount;

      statusLine->setText(QString("No of Models Listed: %1").arg(count));
    }

  private:
    AllChainDataModel mAllDataModel;
    SelectedChainDataModel mSelectedDataModel;
    unsigned int mSelectedDataCount, mAllDataCount;
    QSignalMapper* mSignalMapper;
    QMap<int, bool> mHeaderShown;
};



#endif // GUI_MODELINSPECTORDOCK_H
