#ifndef PARAMETERSPINBOX_H
#define PARAMETERSPINBOX_H

#include <QSpinBox>
#include <QMap>
#include <QDebug>

struct ColumnPlotInfo {
  int columnNumber;
  bool enabled;
  QString columnLabel;
  operator int () { return columnNumber; }
  operator int () const { return columnNumber; }
};

class ParameterSpinBox: public QSpinBox
{
  Q_OBJECT
  public:
    ParameterSpinBox(QWidget* parent=0): QSpinBox(parent) {}
    void setParameterNames(const QMap<int, ColumnPlotInfo>& n) {
      mNameMap = n;
      // workaround to force updating
      int oldVal = value();
      setValue(minimum());
      setValue(oldVal);
    }

    virtual QString textFromValue(int value) const {
      if (mNameMap.contains(value)) {
        return mNameMap[value].columnLabel;
      } else {
        return QSpinBox::textFromValue(value);
      }
    }

  private:
    QMap<int, ColumnPlotInfo> mNameMap;
};
#endif // PARAMETERSPINBOX_H
