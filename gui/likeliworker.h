#ifndef LIKELIWORKER_H
#define LIKELIWORKER_H

#include <QtCore/QThread>

#include "cmbmainwindow.h"
typedef void (CmbMainWindow::*Worker) ();

/*!
  LikeliWorker is a Thread which calls one function of CmbMainWindow.
  The name is legacy and in cmbeasy, there are currently three objects
  of type LikeliWorker: One for calculating anisotropies, one for likelihood
  plotting and one for initializing the WMAP data.

  \param c The CmbMainWindow
  \param w The function to call
  \param e If true, call function only once. Otherwise, call forever.

*/
class LikeliWorker : public QThread { 
  CmbMainWindow *cmb;
  Worker worker;

 public:
  bool ende;
  LikeliWorker(CmbMainWindow *c, Worker w, bool e=false) : QThread() ,cmb(c) , worker(w),  ende(e)  {};

 protected:
  void run() {
    for (;;) {
    (cmb->*worker)();
    if (ende) break;
    }
  }
};

#endif
