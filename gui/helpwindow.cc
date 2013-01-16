/****************************************************************************
** $Id:  qt/helpwindow.cpp   3.0.1   edited Oct 31 20:35 $
**
** Copyright (C) 1992-2000 Trolltech AS.  All rights reserved.
**
** This file is part of an example program for Qt.  This example
** program may be used, distributed and modified without limitation.
**
*****************************************************************************/

#include "helpwindow.h"
#include <qstatusbar.h>
#include <QPixmap>
#include <QIcon>
#include <QToolButton>
#include <QToolBar>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <QPrinter>
#include <QPrintDialog>
#include <QFileDialog>
#include <QApplication>
#include <QComboBox>
#include <QEvent>
#include <QPainter>
#include <QMenuBar>
#include <QFrame>
#include <QMenu>
#include "cmbmainwindow.h"

#include <ctype.h>

HelpWindow::HelpWindow( const QString& home_, const QString& _path,
			QWidget* parent, const char *name )
    : QMainWindow( parent ),
      pathCombo( 0 ), selectedURL()
{
    setObjectName( name );
    setAttribute( Qt::WA_DeleteOnClose );

    readHistory();
    readBookmarks();

    browser = new QTextBrowser( this );

    browser->setSearchPaths( QStringList() << _path );
    browser->setFrameStyle( QFrame::Panel | QFrame::Sunken );
    connect( browser, SIGNAL( textChanged() ),
	     this, SLOT( textChanged() ) );

    setCentralWidget( browser );

    if ( !home_.isEmpty() )
	browser->setSource( home_ );

    connect( browser, SIGNAL( highlighted( const QString&) ),
	     statusBar(), SLOT( showMessage( const QString&)) );
//    connect( browser, SIGNAL( anchorClicked( const QUrl& ) ),
//	     browser, SLOT( setSource( const QUrl& ) ) );


    resize( 640,700 );

    QMenu* file = new QMenu(tr("&File"), this);
    file->addAction( tr("&New Window"), this, SLOT( newWindow() ), Qt::ControlModifier+Qt::Key_N );
    file->addAction( tr("&Open File"), this, SLOT( openFile() ), Qt::ControlModifier+Qt::Key_O );
    file->addAction( tr("&Print"), this, SLOT( print() ), Qt::ControlModifier+Qt::Key_P );
    file->addSeparator();
    file->addAction( tr("&Close"), this, SLOT( close() ), Qt::ControlModifier+Qt::Key_Q );
    file->addAction( tr("E&xit"), qApp, SLOT( closeAllWindows() ), Qt::ControlModifier+Qt::Key_X );

    // The same three icons are used twice each.
    QPixmap back(CmbMainWindow::cmbeasyDir("/pix/back.png"));
    QPixmap forward(CmbMainWindow::cmbeasyDir("/pix/forward.png"));
    QPixmap home(CmbMainWindow::cmbeasyDir("/pix/gohome.png"));
    QIcon icon_back(back);
    QIcon icon_forward(forward);
    QIcon icon_home(home); 

    QMenu* go = new QMenu( tr("&Go"), this );
    backwardId = go->addAction( icon_back,
				 tr("&Backward"), browser, SLOT( backward() ),
				 Qt::ControlModifier+Qt::Key_Left );
    forwardId = go->addAction( icon_forward,
				tr("&Forward"), browser, SLOT( forward() ),
				Qt::ControlModifier+Qt::Key_Right );
    go->addAction( icon_home, tr("&Home"), browser, SLOT( home() ) );

    QMenu* help = new QMenu( tr("&Help"), this );
    help->addAction( tr("&About ..."), this, SLOT( about() ) );
    help->addAction( tr("About &Qt ..."), this, SLOT( aboutQt() ) );

    hist = new QMenu( tr( "History" ), this );
    QStringList::Iterator it = history.begin();
    for ( ; it != history.end(); ++it )
	mHistory[ hist->addAction( *it ) ] = *it;
    connect( hist, SIGNAL( triggered( QAction* ) ),
	     this, SLOT( histChosen( QAction* ) ) );

    bookm = new QMenu( tr( "Bookmarks" ), this );
    bookm->addAction( tr( "Add Bookmark" ), this, SLOT( addBookmark() ) );
    bookm->addSeparator();

    QStringList::Iterator it2 = bookmarks.begin();
    for ( ; it2 != bookmarks.end(); ++it2 )
	mBookmarks[ bookm->addAction( *it2 ) ] = *it2;
    connect( bookm, SIGNAL( triggered( QAction* ) ),
	     this, SLOT( bookmChosen( QAction* ) ) );

    menuBar()->addMenu( file );
    menuBar()->addMenu( go );
    menuBar()->addMenu( hist );
    menuBar()->addMenu( bookm );
    menuBar()->addSeparator();
    menuBar()->addMenu( help );

    forwardId->setEnabled(FALSE);
    backwardId->setEnabled(FALSE);
    connect( browser, SIGNAL( backwardAvailable( bool ) ),
	     this, SLOT( setBackwardAvailable( bool ) ) );
    connect( browser, SIGNAL( forwardAvailable( bool ) ),
	     this, SLOT( setForwardAvailable( bool ) ) );


    QToolBar* toolbar = new QToolBar( this );
    addToolBar( toolbar );
    QToolButton* button;

    button = new QToolButton( toolbar );
    button->setDefaultAction( new QAction( icon_back, tr("Backward"), this ) );
    connect( button->defaultAction(), SIGNAL( triggered() ), browser, SLOT(backward()) );
    toolbar->addWidget( button );
    connect( browser, SIGNAL( backwardAvailable(bool) ), button, SLOT( setEnabled(bool) ) );
    button->setEnabled( FALSE );

    button = new QToolButton( toolbar );
    button->setDefaultAction( new QAction( icon_forward, tr("Forward"), this ) );
    connect( button->defaultAction(), SIGNAL( triggered() ), browser, SLOT(forward() ) ),
    toolbar->addWidget( button );
    connect( browser, SIGNAL( forwardAvailable(bool) ), button, SLOT( setEnabled(bool) ) );
    button->setEnabled( FALSE );

    button = new QToolButton( toolbar );
    button->setDefaultAction( new QAction( icon_home, tr("Home"), this ) );
    connect( button->defaultAction(), SIGNAL( triggered() ), browser, SLOT(home() ) );
    toolbar->addWidget( button );

    toolbar->addSeparator();

    pathCombo = new QComboBox( toolbar );
    pathCombo->setEditable( true );
    toolbar->addWidget( pathCombo );
    connect( pathCombo, SIGNAL( activated( const QString & ) ),
	     this, SLOT( pathSelected( const QString & ) ) );

    pathCombo->addItem( home_ );
}


void HelpWindow::setBackwardAvailable( bool b)
{
    backwardId->setEnabled(b);
}

void HelpWindow::setForwardAvailable( bool b)
{
    forwardId->setEnabled(b);
}


void HelpWindow::textChanged()
{
    if ( browser->documentTitle().isNull() ) {
	setWindowTitle( "Helpviewer - " + browser->source().path() );
	selectedURL = browser->source().path();
    }
    else {
	setWindowTitle( "Helpviewer - " + browser->documentTitle() ) ;
	selectedURL = browser->documentTitle();
    }

    if ( !selectedURL.isEmpty() && pathCombo ) {
	bool exists = FALSE;
	int i;
	for ( i = 0; i < pathCombo->count(); ++i ) {
	    if ( pathCombo->itemText( i ) == selectedURL ) {
		exists = TRUE;
		break;
	    }
	}
	if ( !exists ) {
	    pathCombo->insertItem( 0, selectedURL );
	    pathCombo->setCurrentIndex( 0 );
	    mHistory[ hist->addAction( selectedURL ) ] = selectedURL;
	} else
	    pathCombo->setCurrentIndex( i );
	selectedURL = QString();
    }
}

HelpWindow::~HelpWindow()
{
    history.clear();
    QMap<QAction*, QString>::Iterator it = mHistory.begin();
    for ( ; it != mHistory.end(); ++it )
	history.append( *it );

    QFile f( QDir::currentPath() + "/.history" );
    f.open( QIODevice::WriteOnly );
    QDataStream s( &f );
    s << history;
    f.close();

    bookmarks.clear();
    QMap<QAction*, QString>::Iterator it2 = mBookmarks.begin();
    for ( ; it2 != mBookmarks.end(); ++it2 )
	bookmarks.append( *it2 );

    QFile f2( QDir::currentPath() + "/.bookmarks" );
    f2.open( QIODevice::WriteOnly );
    QDataStream s2( &f2 );
    s2 << bookmarks;
    f2.close();
}

void HelpWindow::about()
{
    QMessageBox::about( this, "HelpViewer for Cmbeasy",
			"<p>This is a simple browser "
			"for the documentation of Cmbeasy's graphical user interface</p>"
			);
}


void HelpWindow::aboutQt()
{
    QMessageBox::aboutQt( this, "QBrowser" );
}

void HelpWindow::openFile()
{
#ifndef QT_NO_FILEDIALOG
    QString fn = QFileDialog::getOpenFileName( this );
    if ( !fn.isEmpty() )
	browser->setSource( fn );
#endif
}

void HelpWindow::newWindow()
{
    ( new HelpWindow(browser->source().path(), "qbrowser") )->show();
}

void HelpWindow::print()
{
#ifndef QT_NO_PRINTER
    QPrinter printer;//(QPrinter::HighResolution );
    printer.setFullPage(TRUE);
    QTextDocument *document = browser->document();

    QPrintDialog *dlg = new QPrintDialog( &printer, this );
    if ( dlg->exec() != QDialog::Accepted )
        return;
    document->print( &printer );
#endif
}

void HelpWindow::pathSelected( const QString &_path )
{
    browser->setSource( _path );
    QMap<QAction*, QString>::Iterator it = mHistory.begin();
    bool exists = FALSE;
    for ( ; it != mHistory.end(); ++it ) {
	if ( *it == _path ) {
	    exists = TRUE;
	    break;
	}
    }
    if ( !exists )
	mHistory[ hist->addAction( _path ) ] = _path;
}

void HelpWindow::readHistory()
{
    if ( QFile::exists( QDir::currentPath() + "/.history" ) ) {
	QFile f( QDir::currentPath() + "/.history" );
	f.open( QIODevice::ReadOnly );
	QDataStream s( &f );
	s >> history;
	f.close();
	while ( history.count() > 20 )
	    history.erase( history.begin() );
    }
}

void HelpWindow::readBookmarks()
{
    if ( QFile::exists( QDir::currentPath() + "/.bookmarks" ) ) {
	QFile f( QDir::currentPath() + "/.bookmarks" );
	f.open( QIODevice::ReadOnly );
	QDataStream s( &f );
	s >> bookmarks;
	f.close();
    }
}

void HelpWindow::histChosen( QAction* i )
{
    if ( mHistory.contains( i ) )
	browser->setSource( mHistory[ i ] );
}

void HelpWindow::bookmChosen( QAction* i )
{
    if ( mBookmarks.contains( i ) )
	browser->setSource( mBookmarks[ i ] );
}

void HelpWindow::addBookmark()
{
   mBookmarks[ bookm->addAction( windowTitle() ) ] = browser->source().path();
}
