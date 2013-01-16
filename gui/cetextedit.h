#ifndef CETEXTEDIT_H
#define CETEXTEDIT_H

#include <QTextEdit>
//! A QTextEdit with a smaller standard sizeHint(), so that the
//! PlotControllDockWidget can have a reasonable standard size
class CETextEdit: public QTextEdit
{
  public:
    CETextEdit(QWidget* parent): QTextEdit(parent) {}
    QSize sizeHint() const
    {
      int f = 2 * frameWidth();
      QSize sz(f, f);
      int h = fontMetrics().height();
      sz += QSize(4*h, 4*h);
      return sz.boundedTo(QSize(5*h, 5*h));

    }
};

#endif // CETEXTEDIT_H
