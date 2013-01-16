#ifndef Q4BUTTONGROUP_H
#define Q4BUTTONGROUP_H
// Simple class to add the Qt3 mapping functionality
// back to Qt4 QButtonGroup
// written by: Matthias Ettrich <ettrich@trolltech.com
// email to: qt4-preview-feedback@trolltech.com
// dated: 2005-06-03 14:11
// "Here's a little class that you can use: [..] Hope this helps until 4.1"

#include <QButtonGroup>
#include <QMap>

class Q4ButtonGroup : public QButtonGroup
{
    Q_OBJECT
public:
    Q4ButtonGroup( QObject *parent = 0 ):QButtonGroup( parent ){
        connect( this, SIGNAL( buttonClicked( QAbstractButton* ) ),
                this, SLOT( mapButton( QAbstractButton* ) ) );
    }
    void addButton( QAbstractButton *button, int id = -1 ) {
        QButtonGroup::addButton( button );
        if ( id != -1 )
            mapping[ button ] = id;
    }

    QAbstractButton *button( int id ) const { return mapping.key( id ); }
    int id( QAbstractButton *button ) const { return mapping.value( button ); }

    int checkedId() const { return mapping.value( checkedButton(), -1 ); }

signals:
        void buttonClicked( int );

private slots:
        void mapButton( QAbstractButton *button ) {
            emit buttonClicked( id( button ) );
        }
private:
    QMap<QAbstractButton*, int> mapping;
};
#endif //Q4BUTTONGROUP_H
