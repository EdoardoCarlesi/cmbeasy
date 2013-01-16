/****************************************************************************
** $Id:  qt/helpwindow.h   3.0.1   edited Oct 12 12:18 $
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#ifndef HELPWINDOW_H
#define HELPWINDOW_H

#include <QMainWindow>
#include <QTextBrowser>
#include <QString>
#include <QMap>

class QComboBox;
class QPopupMenu;
class QAction;

class HelpWindow : public QMainWindow
{
    Q_OBJECT
public:
    HelpWindow( const QString& home_,  const QString& path, QWidget* parent = 0, const char *name=0 );
    ~HelpWindow();

private slots:
    void setBackwardAvailable( bool );
    void setForwardAvailable( bool );

    void textChanged();
    void about();
    void aboutQt();
    void openFile();
    void newWindow();
    void print();

    void pathSelected( const QString & );
    void histChosen( QAction* );
    void bookmChosen( QAction* );
    void addBookmark();

private:
    void readHistory();
    void readBookmarks();

    QTextBrowser* browser;
    QComboBox *pathCombo;
    QAction* backwardId, *forwardId;
    QString selectedURL;
    QStringList history, bookmarks;
    QMap<QAction*, QString> mHistory, mBookmarks; // ###FIXME: use QAction->data() instead
    QMenu *hist, *bookm;
};





#endif

