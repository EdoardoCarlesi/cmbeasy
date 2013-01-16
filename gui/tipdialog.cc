#include "tipdialog.h"
#include "cmbmainwindow.h"
#include <iostream>

TipDialog::TipDialog(unsigned int *tip, bool *again,QWidget *parent, const char * name, bool modal)  : QDialog(parent), Ui::DesignTipDialog() , Tip(tip) , ShowAgain(again) {

  setupUi(this);

  AgainBox->setChecked(*ShowAgain);
  ifstream in(CmbMainWindow::cmbeasyDir("/resources/gui_data/tips.txt").toLatin1().data());
  
  unsigned int Size;
  char buffer[1000];
  unsigned int count = 0;
  in >> Size;
  in.getline(buffer,999);
  TipString.resize(Size);
  while (in.getline(buffer,999) && count < Size) {
    if (buffer[0] == '#') count++;
    else {
      TipString[count] += buffer;
      TipString[count] += " ";
    }
  }
  connect(NextButton,SIGNAL(clicked()),this,SLOT(showTip()));
  connect(DoneButton,SIGNAL(clicked()),this,SLOT(closeIt()));
  
  cout << "-========---------=================="<<endl;
  cout << Size << endl;
  for (unsigned int i = 0; i < Size; i++) cout << TipString[i].toStdString() << "\n" << endl;

  showTip();
}

void TipDialog::showTip() {
  cout << "SHOWTIP: " << *Tip <<endl;
  if (*Tip >= TipString.size()) cout << "TipDialog:: Tip Nr ( " << *Tip << " ) too large" << endl;
  else {
    TipEdit->setPlainText(TipString[(*Tip)++]);
    if (*Tip >= TipString.size()) *Tip = 0;
  }
}

void TipDialog::closeIt() { close(); }

bool TipDialog::close() {  
  *ShowAgain = AgainBox->isChecked();
  QDialog::close();
  return true;
}

