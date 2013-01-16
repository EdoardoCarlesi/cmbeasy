#ifndef CL_H
#define CL_H

#include "cleanvector.h"

class Spline;
class Anchor;

/*!
  The multipole spectrum is stored as
  Splines. In order to pass the various
  spectra (possibly for several spectral
  indices) around, CL has been
  invented.

  Accessing for instance the scalar 
  temperature anisotropy is quite simple:

  \code
  CL mycl;
  mycl.ts[0]->dump("scalar anisotropy for first spectral index");
  \endcode
*/
struct CL
{
  typedef CleanVector<Spline*> SplineCleanVector;
  SplineCleanVector ts, tt, es, et, bt,
                    cs, ct, kk, tk, bs;

  CL* clone();
  void resize(unsigned int n); //!< resize all Cl splines to n
  void clear(); //!< clear all splines
  /*! Multiply all splines for spectral index number i with
   *  factor x */
  void mul(const double x, const int i);
  Spline* createTotalTTSpline(const int i, const std::string name=std::string(), Anchor* a=0) const;
  Spline* createTotalEESpline(const int i, const std::string name=std::string(), Anchor* a=0) const;
  Spline* createTotalTESpline(const int i, const std::string name=std::string(), Anchor* a=0) const;
  Spline* createTotalBBSpline(const int i, const std::string name=std::string(), Anchor* a=0) const;

private:
  Spline* createTotalXXSpline(const int i, const SplineCleanVector& scalar, const SplineCleanVector& tensor,
                              const std::string name=std::string(), Anchor* a=0) const;
};


#endif

//==========================================
//
// Some online documentation following.
// I attach it here, as CL is just a small class and
// here it should not disturb too much
//
//==========================================

/*! 
\page classesobj Classes, Objects and Sub-Classing

Sometimes, it makes sense to group some data and functions
acting on it in one class. There are three obvious advantages of 
doing this: 

1. The code is more readable, as a logic block of the code is 
contained within one class and often within one file.

2. As the class defines an interface through which other
parts of the code can interact with the class, one can change
the implementation of certain functions, without changing 
other code, for as long as the interface (i.e. the arguments
of the functions) stay the same.

3. One can use sub-classing to inherit from an existing class.
The sub-class has the same functionality of the partent class,
yet it can extend this functionality if needed.  In addition,
the sub-class can re-implement functions of the base-class.
Which will be executed instead of the base-class'es versions.

Usually, the definition of the interface (and layout) of the
class resides in a so called header file. A small example class,
could look like this:

\code
class Square { 
  private:
  double x;
  public:
  Square(double argument);
  virtual double getSquare();
};
\endcode

This code would go into the file "square.h". The "private:" and
"public:" statements indicate regions within the class layout
that are for internal class use only or for everyone to access.
Please note the "virtual" in front of double getSquare(). This
tells the compiler that this function can (but doesn't have to) be re-implemented by
a sub-class. 
The implementation of this class would go into "square.cc" (or ".cpp") 
and could look like this:

\code
#include "square.h"

Square::Square(double argument) : x(argument) {}

double Square::getSquare() {
  return x*x;
}
\endcode

In the first line, we include the header file (otherwise the
compiler would not know about a class called Square). The next
line is a so called constructor. Whenever you create an
object of this class, a constructor will be called. There may
me many constructors (often provided for convenience), see e.g. Spline.
In our case, you can only construct a Square object, if you
provide an argument. This argument will be stored in the variable x (the
statement after the : in the constructor initializes x with the argument). 

The only function of our class is called getSquare. The fact that
it belongs to the class Square is indicated by the Square:: in front of getSquare()
(there can be many classes with the same function-name and the Square:: in
front of the cuntion distinguishes between those).
We return a double x*x, just as promised in the definition of "square.h"

Say, we wanted to extend this class. Lets call it SuperSquare. The header file
could look like this:
\code
#include "square.h"

class SuperSquare : public Square {

public:
  SuperSquare(double argument); 
  double getSquare();
};
\endcode
The function getSquare() is not virtual this time. Hence, a sub-class of SuperSquare
cannot re-implement it. We could however, easily declare it virtual, the any
sub-class could re-implement it.

The implementation would be in "supersquare.cc" and could look like this
\code
#include "supersquare.h"

SuperSquare::SuperSquare(double a) : Square(a) {}

double getSquare() {
  return pow(x,2);
}
\endcode

Here, the constructor calls its base class with the necessary argument. The
  base-class in turn, will be constructed and initialize its variable "x" with
this argument. 

  Our re-implementation of getSquare() uses the math-library pow() function.
I have to admit that the SuperClass wouldn't be as fast. 

However this
is certainly a valid implementation of the task: "take the variable x stored in
Square and return the square of it!"

In some other part of the code, we can now use the Square class (we have
to #include "supersquare.h" of course!)
\code
Square old(3);
SuperSquare better(3), hundred(10);

double x = old.getSquare();
double y = better.getSquare();
double z = hundred.getSquare();

\endcode

Here, we have three objects, one of type Square and two of type SuperSquare.
Both old and better will give the same result (namely 9). This is important
when using pointers, see \ref pointers.
*/

/*!
\page virt Virtual functions
Within a class definition, you can specify functions as being "virtual". Virtual functions
can be  re-implemented by sub-classes. Say you declare the function within the
class Synchronous:

\code
virtual void getReady();
\endcode

Whatever your class needs to do to get ready will be done in getReady(). 
Now, a sub-class (say QuintSynchronous) can also delcare this function:

\code
virtual void getReady();
\endcode

However, it may do something else to get ready. It may, however still wish to call
the getReady() function of its parent class, this can always be done by specifying
the name of the parent class to distinguish between the own implementation of 
getReady() and the one of the parent class:
\code
Invariant::getReady();
\endcode

\subsection pure Purely Virtual functions.
Sometimes, a base class doesn't implement any functionality, but defines a function
that the sub-classes necessarily need to provide. The rest of the code can then
  rely on using this class, no matter what kind of sub-class is used. Such functions
  are called purely virtual and defined using "=0":
\code
  void getReady()=0;  // I don't implement anything, but the sub-classes need to !
\endcode

This is extremely useful when using pointers to an object, as the rest of the
  code can stay umodified when substituting one sub-class for another, if the
rest of the code only depends on the functions declared in the parent-class
    (for instance Quintessence or Perturbation).
*/

/*!
\page objp Objects, Pointers and References
    The easiest way to create an object is to declare it:
\code
Cosmos cosmos;
Cosmos alternative;
\endcode

Here, two objects of type Cosmos are created. However,
when calling subroutines, it is usually not a clever
idea to hand over something as sophisticated as
a cosmos:
\code
void DoSomething(Cosmos c) {
...
}
\endcode

A call of DoSomething(alternative) would try to <b>copy</b>
the entire cosmos object. This is very expensive. Hence, it 
is much cheaper to provide a <b>pointer</b>, or a 
<b>reference</b> to an object. A pointer holds the
address at which a certain object resides in memory. The same
is true for a reference. The distinction between these two
is that a pointer can change its value, while a reference can't.
So this would work:
\code 
Cosmos *cosmos = new Cosmos();  // construct a cosmos, return a pointer 
cosmos = new Cosmos();  // new cosmos, re-use the pointer
\endcode

The "*" before cosmos (or after Cosmos) indicates that cosmos is
in this case not of type Cosmos, but a pointer to an object with that
type.
When calling functions acting on the object, we have to either
use the "->" or "*" operator. These two calls are equivalent:
\code
cosmos->printStatus();
*cosmos.printStatus();
\endcode

The "*" in front of the pointer cosmos de-references it.

As mentioned, a reference can only be initialized once (until it is destroyed and created again):
\code
Cosmos alternative;
Cosmos &cosmos = alternative;
\endcode
Whereas alternative is a honest Cosmos object, the cosmos reference below (indicated by "&") is technically speaking only a constant pointer. Whatever you
do to the cosmos, will act on alternative. They are two names for one and
only one Cosmos object, namely alternative in our case.

\section poly Polymorphism
One big advantage of using references and pointers is polymorphism: You can for instance define a function that has a pointer to Perturbation as argument:
\code
void mylittlefunction(Perturbation *p) {
  p->getReady();
}
\endcode

Which function getReady() will be executed ? The answer is: it depends!
If you call 
\code
QuintSynchronous perturb;
mylittlefunction(&perturb);
\endcode
where the "&" returns a pointer to the perturb object, the code
of QuintSynchronous::getReady() will be executed. If you plug
in Invariant for QuintSynchronous in the above, it will be
Invariant::getReady that is called. As long as your function
"mylittlefunction()" only cares about the fact that the perturbation
object gets ready, it doesn't mind which one it asks to
get ready. Hence, when you come up with a new cosmology
and a new Perturbation class, most parts of CMBEASY will
just not care and you don't have to change the code there. All
you need to do is provide a sub-class of Perturbation and 
that's it !

\section summa Syntax Summary
Here, we summarize the usage of objects, pointers and references:
\code
Cosmos object;
Cosmos *pointer = new Cosmos();
Cosmos &reference1 = object;
Cosmos &reference2 = *pointer;

object.reset();

pointer->reset();
*pointer.reset();  // equivalent

reference1.reset();
reference2.reset();
\endcode
*/

/*!
\page templ Templates
Sometimes, a function (or class) can act on many different
data types in exactly the same manner. Instead of writing 
several functions for each data - type, one uses templates.
The main usage of templates within CMBEASY occurs for
C++ standard library classes, such as vector, list or map.
However, to have a realy world CMBEASY example, let us
look at the Model::write() function:
\code
template<class T> static void Model::write(ostream& o,const T& s) {
    o.write((const char*) &s,sizeof(T));
}
\endcode

Its purpose is writing the binary content of a variable to a file.
The template<class T> indicates that the compiler should
create a function Model::write() for each data type specified
by the user. Hence, code like this:
\code
ofstream file("myfile.dat");
double x = 1.234;
int y = 5;
Model::write<double>(file, x);
Model::write<int>(file,y);
\endcode

would write a double and an int in binary format (i.e. 8 bytes for the double and 4 bytes for the int) to the file. 

\section \tempusage Main usage: standard library
As mentioned, the main usage for templates within CMBEASY arises
within the standard library. CMBEASY frequently uses data types
like vector, list and map. All of these are template classes. Hence, one
needs to specify the type of objects that will be stored in a certain
vector, list or map:
\code
vector<double> one_dim_array;
vector< vector< double> > two_dim_array;
map<string, double> look_up;

look_up["hello world"] = 1.032; // a mapping of string -> double

\endcode

Please note that the data types used as templates can be user-defined
and in principle very complex (and templates themselves). 
*/

/*!
\page itrat Iterators
Quite often, it is useful to perform a calculation for all variables stored in a standard container like vector, list or map.
For all of these containers, the syntax of looping through the content is exactly the same:
\code
vector<double> my_vector;
int c = 0;
for (vector<double>::iterator i = my_vector.begin(); i != my_vector.end(); i++) {
  double x = *i;
  double y = my_vector[c]; // equivalent
  c++;
} 
\endcode

Instead of a vector<double> we could have used any other standard container, e.g. a 
\code
map<string, Spline*> alternative;
\endcode

The iterator here is initialized with the begin() of the container. The loop stops,
when the iterator equals the .end() [the end() is not part of the container, i.e. the last element is begore end()]. Increasing the iterator by i++ goes to the next element. In our example, x and y will have exactly the same
value.  
Whereas for a vector this is not a very sensible thing to do, it is very useful for map and list. Elements of a map are stored in sorted manner, i.e.
\code
map<double, Spline*> splines;
for (double k = 3; k > 0; k -= 0.03) {
  splines[k] = new Spline();
}
\endcode

would store some pointers to splines each indexed by a k-value in the
map splines. The order in which this examples initializes the splines 
is from large k to low k. Yet, within map, it will be stored from low k to large k.
Hence, when traversing a map using an iterator, one will always start
at the lowest k-value.
*/
/*! \page subcls Modifying CMBEASY through sub-classing 


If you would like to make modifications, there are two things, that will probably interest you most:

\section mintro1 1. How are my changes surviving updates ? 

It is certainly quite cumbersome to keep track of changes in between releases. However, there
is a clever way to avoid a lot of work: <b>use sub-classing!</b>

The nicest solution is to inherit from an existing class and re-implement a virtual function. 
For instance, if you are unhappy with some handling in Synchronous, say getReady(), inherit 
from Synchronous and re-implement it

The header would look like this:
\code
\  #include "synchronous.cc"
    MoreHappySync : public Synchronous {
  public:
    MoreHappySync(Cosmos *c);
    void getReady(const ControlPanel&);
};
\endcode

The implementation, in  morehappysync.cc would look like:
\code
MoreHappySync::MoreHappySync(Cosmos *c) : Synchronous(c) {}  // constructor

void MoreHappySync::getReady(const ControlPanel& control) {
  double x = 3.141;
  double y = 2*x;  // make me more happy
  Synchronous::getReady(control);  // call Base-Class function (at least if you like)
}
\endcode

In our case, it's been easy to re-implement, cause getReady() is declared
virtual (i.e. you may substitute your solution for the existing one) in synchronous.h.

<b>If you find that a function should be virtual and is not (I didn't think all cases through),
please let me know. I'll make it virtual then and you can sub-class!</b>

\section mintro2 2. What, if the change is something that no base class would need, can my changes still survive ?

This is not so nice from the aesthetic point of view. However, no problem in principle.

Let us say, you were to add a flag to ControlPanel. Let us call this bool MyFlag. Say you would
add this to ControlPanel. If someone (like 
me) changes ControlPanel and you would like to update, you either would have to compare
your version to the update version,  or install the update version and add your flag (ok, this example
is rather trivial, but it is very important for more complicated stuff).

The way around this is to define your own ControlPanel, let us call it MyPanel. The header file could
look like this:

\code
#include "controlpanel.h"
class MyPanel : public ControlPanel {
public:
bool MyFlag;
};
\endcode

In the driver, you would instead of ControlPanel use a MyControl:

\code
MyControl control;
\endcode

Now, as MyControl is certainly a control, the whole package just won't care about
your added Flag. And if ControlPanel gets some more features, MyControl will also know
about it immediately (cause it inherits). 

Say you would like to access your new control in a function f:
\code  
f(ControlPanel& control) { ...... }
\endcode
which takes ControlPanel as argument. Now you can just
substitute ControlPanel for MyControl. Change the code to
\code
f(MyControl& control)
\endcode

Another possibility would be to make an explicit cast from ControlPanel to MyControl.
In this case, f keeps its old argument (ControlPanel& control). 
Within the function body, you can then upgrad control to MyControl:
\code
MyControl& myversion = (MyControl&) control;
bool this_is_my_flag = control.MyFlag;
\endcode

The (MyControl&) tells the compiler that the reference to control which he thought
of must be of type ControlPanel& is indeed not only a ControlPanel, but a MyControl.
You can, of course only do this, if you are sure that the ControlPanel given
to f is indeed of type MyControl and not just ControlPanel. 
What you are more or less doing here is explicitly telling the compiler that
the object will indeed be a MyControl. 
As mentioned above, using virtual functions that are added to the base classes are
by far the nicer option. 
*/

/*!
  \page uss Using Splines
  
  Within CMBEASY, Splines are very useful. I suggest using 
  Splines or a SplineWeb, whenever you need to store a quantity
  that you will need to process later.  The basic usage of 
  Splines is explained in the documentation of the Spline class. 
  However, to summarize some frequently used functions, I repeat
  them here (assuming there is a Spline s with data and arm()ed):
\code

// Interpolating
Spline *pointer = &s;  // for instructions
double y = s(0.3);  interpolate at 0.3
double y2 = pointer->fastY(0.3); // the same
double y3 = (*pointer)(0.3); // the same

// Deriving and integrating

Spline drv;
s.derive(&drv);  // derive s and store result in drv

double y = s.integrate(); // integrate spline s 
double y2  = s.integrate(-1,1);  // integrate from x = -1... 1 


// Finding extrema

double x = s.maximum();   // x-value of spline's maximum
double y = s(x);              // interpolated y-value of maximum


// Dump to file
s.dump("debug_it");  // dump x and y data to file "debug_it.dat"
s.dump("debug_it",false); // dump data to "debug_it.dat" *and* interpolated values (10 by default per data-interval) to "debug_it.plt". Spline needs to be armed for this


// bounds and positions

int n = s.last(); // position of the last point in the data arrays
int size = s.size(); // size of the array (equal n+1)

double x_min = s.start();  
double x_min = s.stop();

double y_first = s.front();  // y-data of first point
double y_last = s.back();  // y-data of last point

double x = s.stop(3); // x-data of the fourth point from the back: n-3

// convenient sampling of spline in linear fashion, in our 
// case making steps of size (stop()-start())/100:
for (double x = s.start(); x <= s.stop(); x += s.range(100) ) .....

// accessing and manipulating data
double y = s[5]; // / ydata[5], i.e. point number 6 (we count from 0 on)
double y = s.y(5);  // ydata[5], i.e. point number 6 (we count from 0 on), same as operator []
double x = s.x(5);  // x-data at point number 6 

s += 1.0; // add 1.0 to the spline y-data
s *= 3.0;  // multiply the y-data by 3.0.
s.mul(7,0.3); // multiply y-data point number 7 by 0.3
\endcode
*/

/*!
\page modifying How to modify the code

The purpose of this documentation is to provide you with a working
knowledge of ways to modify the code. As you may not be familiar with
some features of C++, I will give a short introduction to frequently used
C++ elements  in CMBEASY. 

\subsection gcpps General C++ stuff 

1. \ref classesobj

2. \ref virt

3. \ref objp

4. \ref templ 

5. \ref itrat

\subsection cmbeasystuff cmbeasy stuff

6. \ref subcls

7. \ref uss
*/



/*!
\mainpage CMBEASY

This is the documentation of CMBEASY. Go to the <a href="www.cmbeasy.org">homepage</a> for more information. 

\section intro Introduction
 
CMBEASY is a software package for calculating
	the evolution of density fluctuations in the universe. Most notably, the Cosmic Microwave background temperature anisotropies. 
       The code is based on the CMBFAST package by Uros Seljak and Mattias Zaldarriaga. 
Without their donation of the CMBFAST package to the public domain, CMBEASY would not exist. 
      
           Even though its ancestor is a Fortran program, CMBEASY is fully object oriented C++. This considerably simplifies manipulations and extensions of the code. In addition, a powerful Spline class can be used to easily store and visualize data. 
    Many features of the CMBEASY package are also accessible via a graphical user interface. This may be helpful for gaining intuition, as well as educational purposes.


\section installation Installation

Download the appropriate source code package, uncompress it using 'gzip', then untar it 
using 'tar -xvf package_file', where package_file should be replaced by the name of the
package you've downloaded. After un-taring, there will be a directory called
'cmbeasy'. Read the file 'INSTALLATION' on how to proceed.
 
\section howto How to find out about the classes

In order to get a quick overview of the classes, click on the 'Class Hierarchy' link
on top of this page. Don't be afraid of the numerous classes displayed ! Most of them
are quite small, having only a very limited purpose. 

Click on the class you would like to learn more about. Then scroll down to the
section 'detailed description'. Quite often, you will even find examples on how to 
use this class. 
In addition, CMBEASY itself is an example of how the classes work.

\section first run

For some first runs, I suggest calling "cmb" or "xcmb" with a control-file as argument.

> cmb resources/parameters.cmb

or

> xcmb resources/configuration.cfg

The xcmb driver is more versatile than the cmb driver. Choose
whichever you like more.  

Yet, to gain full control, call "cmb" without an argument and change the settings in detailed() in
driver.cc

\section modify Modifying the code

As a first step, I would suggest to change settings in 'driver.cc'. 
For modifications of the background cosmology, QuintCosmos may serve as 
an example.  The same is true for the perturbation classes. It is also very convenient,
to use Spline  and SplineWeb, especially, if you would like to keep track of 
the evolution of some quantity. 

After modifying, just type 'make' in cmbeasy's directory.  This will rebuild the 
package, including your changes.

Please see also \ref modifying on how to modify the code.



*/
