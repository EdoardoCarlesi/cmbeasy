#ifndef PORTING_HELPER_H
#define PORTING_HELPER_H
/****************************************************************
 *   This file is part of the C++ port of the Fortran likelihood
 *   code for the WMAP 7yr release provided by the WMAP team
 *   at http://lambda.gsfc.nasa.gov/ .
 *
 *   The code was ported by Georg Robbers for easier
 *   interfacing with cmbeasy (http://www.cmbeasy.org).
 *   Bugs in this port should be reported to the
 *   cmbeasy authors (bugs@cmbeasy.org).
 ****************************************************************/

//X map the fortran data types
#define REAL_4 float
#define REAL float
#define REAL_8 double
#define COMPLEX std::complex<float>
#define COMPLEX_4 std::complex<float>

#include <complex>
#include <algorithm>
#include <numeric>
#include <functional>
#include <string>

#define GSL_RANGE_CHECK_OFF

#ifdef WMAP_STANDALONE
/* these should be caught when using the standalone version*/
struct Bad_Error {
  std::string s;
  Bad_Error(std::string q="") { s=q; }
};
#endif

template<typename T> class Array1D;
template<typename T> class Matrix;
template<typename T> class Array3D;

typedef Array1D<int> int_1d;

typedef Array1D<REAL_4> real4_1d;
typedef Matrix<REAL_4> real4_2d;
typedef Array3D<REAL_4> real4_3d;

typedef Array1D<REAL> real_1d;
typedef Matrix<REAL> real_2d;
typedef Array3D<REAL> real_3d;

typedef Array1D<REAL_8> real8_1d;
typedef Matrix<REAL_8> real8_2d;
typedef Array3D<REAL_8> real8_3d;

typedef Array1D<double> double_1d;
typedef Matrix<double> double_2d;
typedef Array3D<double> double_3d;

typedef Array1D<COMPLEX> complex_1d;
typedef Matrix<COMPLEX> complex_2d;

/*! \class Range
 *  \brief little helper class holding a range of positive integers
 *
 *  This class is used to indicate the first and last index of each
 *  dimension of the Matrix, Array1D and Array3D classes.
 *  Convenient for simulating fortran-like arrays with first indices
 *  different from 0.
 */
class Range
{
public:
  unsigned int firstIndex; //!< the starting integer of the range
  unsigned int lastIndex;  //!< the last integer included in the range

  //! Construct a range from first to last
  Range( unsigned int first, unsigned int last )
      : firstIndex( first ), lastIndex( last ) {}
};

/*! \class Array1D
 *  \brief yet another 1D Array class
 *
 *   Rather simplistic one-dimensional Array class to make the C++ port of the
 *   WMAP-likelihood look a bit more like the Fortran90 original.
 *
 *   Array1D can encapsulate an Healpix-Array arr, in order to avoid converting
 *   back and forth when interfacing with Healpix.
 *
 *   \warning There is only _very_ limited error checking, and this class is supposed
 *   to only be used for the WMAP 3-year likelihood port.
 */
template<typename T>
class Array1D
{
private:
  unsigned int mFirstIndex; //!< starting index of the Array, 0 by default
  unsigned int mLength;     //!< Array length
  T* mData;                 //!< pointer the the underlying data array
public:
  /*! \brief  Construct an Array1D of length len
   *
   * @param len the length of the Array, 0 by default
   */
  Array1D( unsigned int len = 0 )
      : mFirstIndex(0), mLength(len),
      mData(NULL)
  {
    if ( len > 0 )
      mData = new T[len];
  }

  /*! \brief  Copy constructor
   *
   * @param m the Array to be copied
   */
  Array1D(const Array1D<T>& m)
  {
    *this = m;
  }

  /*! \brief Construct an array with specific Range
   *
   *  Construct an array with base index r.firstIndex,
   *  and size r.lastIndex-r.firstIndex+1.
   *
   *  @param r the range of indices for this Array
   */
  Array1D( Range r )
      : mFirstIndex( r.firstIndex ),
      mLength(r.lastIndex-r.firstIndex+1),
      mData(NULL)
  {
    if ( mLength > 0 )
      mData = new T[mLength];
  }

  ~Array1D()
  {
    delete[] mData;
  }

  //! free the memory held by this Array.
  //! equivalent to resize(0)
  void free()
  {
    resize( 0 );
  }

  /*! \class unsupported_action
   *  \brief exception class
   *
   *  This exception is thrown when trying to
   *  resize an array that operates on top of a
   *  Healpix Array arr.
   */
  class unsupported_action {};

  //! Resize the array to length len.
  //! Old data is deleted.
  void resize( unsigned int len )
  {
    mLength = len;
    mFirstIndex = 0;
    delete[] mData;
    mData = new T[len];
  }

  //! Resize the array to have the range r.
  //! Old data is deleted.
  void resize( Range r )
  {
    mFirstIndex = r.firstIndex;
    mLength = r.lastIndex-r.firstIndex+1;
    delete[] mData;
    mData = new T[mLength];
  }

  /*! \class BoundsViolation
   *  \brief exception thrown when accessing an array out of bounds
   */
  class BoundsViolation {};

  //! return a reference to the element at index idx
  T&       operator() (unsigned int idx)
  {
    //if (idx >= mLength+mFirstIndex) throw BoundsViolation();
    return mData[-mFirstIndex+idx];
  }

  //! return a reference to the element at index idx
  const T& operator() (unsigned int  idx) const
  {
    //if (idx >= mLength+mFirstIndex) throw BoundsViolation();
    return mData[-mFirstIndex+idx];
  }

  //! Returns an Array1D containing the elements in the range r.
  //! The elements in the range are copied to a new array with
  //! range r.
  Array1D<T> operator() (Range r)
  {
    //if (r.lastIndex >= mLength+mFirstIndex || r.firstIndex > mFirstIndex) throw BoundsViolation();
    Array1D out(Range( r.firstIndex, r.lastIndex ) );
    std::copy( mData-mFirstIndex+r.firstIndex, mData-mFirstIndex+r.lastIndex+1, out.data() );
    return out;
  }

  T* data() const { return mData; }                       //!< returns a pointer to the underlying data
  unsigned int firstIndex() const { return mFirstIndex; } //!< returns the base index
  void reindexSelf( unsigned int i) { mFirstIndex = i; }  //!< set the base index to i
  unsigned int length() const { return mLength; }         //!< returns the length of this array
  unsigned int size() const { return mLength; }           //!< returns the size of this Array

  //! assignment operator; the data is copied
  const T& operator= (const T& t)
  {
    std::fill( mData, mData+mLength, t );
    return t;
  }

  /* Make this Array a copy of m. The data is copied.
   *
   */
  Array1D<T>& operator= (const Array1D<T>& m)
  {
    if ( this == &m )
      return *this;
    mLength = m.length();
    mFirstIndex = m.firstIndex();
    delete[] mData;
    mData = new T[mLength];

    std::copy( m.data(), m.data()+m.length(), mData );
    return *this;
  }

  /* Make this Array a copy of the data pointed to by m. The data is copied.
   *
   * @warning: m is assumed to have the same size as this Array1D
   */
  Array1D<T>& operator= (const T* const m)
  {
    std::copy( m, m+mLength, mData);
    return *this;
  }

  //! divide all elements in the array by t
  Array1D<T>& operator/ (const T& t)
  {
    using namespace std;
    transform( mData, mData+mLength,  mData, bind2nd( divides<T>(), t ));
    return *this;
  }

  //! add the array m to this array (elementwise)
  Array1D<T>& operator+= (const Array1D& m)
  {
    using namespace std;
    transform( m.data(), m.data()+m.length(),  mData, mData, plus<T>() );
    return *this;
  }

  //! subtract the array m from this array (elementwise)
  Array1D<T> operator- (const Array1D& m) const
  {
    using namespace std;
    Array1D out(Range( m.firstIndex(), m.length()+m.firstIndex()-1 ) );
    transform( mData, mData+mLength,  m.data(), out.data(), minus<T>() );
    return out;
  }

  //! multiply all elements in the array by t
  Array1D<T>& operator*= (const T& t)
  {
    using namespace std;
    transform( mData, mData+mLength,  mData, bind2nd( multiplies<T>(), t ));
    return *this;
  }

  //! multiply this array by the array m (elementwise)
  Array1D<T>& operator*= (const Array1D& m)
  {
    using namespace std;
    transform( mData, mData+mLength,  m.data(), mData, multiplies<T>() );
    return *this;
  }

  //! set the range r to t
  void set (const Range& r, const T& t)
  {
    std::fill( mData-mFirstIndex+r.firstIndex, mData-mFirstIndex+r.lastIndex+1, t );
  }

  //! return the sum of all the elements in the array
  T sum() const
  {
    return std::accumulate( mData, mData+mLength, 0. );
  }

  //! return the sum of the elements in the range r
  T sum(const Range& r) const
  {
    return std::accumulate( mData-mFirstIndex+r.firstIndex, mData-mFirstIndex+r.lastIndex+1, 0. );
  }
};


/*! \class Matrix
 *  \brief simple Matrix class
 *
 *   A tiny Matrix class to make the C++ port of the
 *   WMAP-likelihood look a bit more like the Fortran90 original.
 *
 *   A Matrix can have differenct base indices for each dimension.
 *
 *   \warning There is no error checking, and this class is supposed
 *   to only be used for the WMAP 3-year likelihood port.
 */

template<typename T>
class Matrix
{

private:
  unsigned int  mNRows, mNCols, mRowsFirstIndex, mColsFirstIndex;
  T* mData;
public:
  //! Construct a Matrix with nrows rows and ncols cols.
  Matrix( unsigned int nrows = 0, unsigned int ncols = 0 )
      : mNRows (nrows), mNCols(ncols),
      mRowsFirstIndex(0), mColsFirstIndex(0),
      mData(NULL)
  {
    if ( nrows > 0 && ncols > 0 )
      mData = new T[nrows * ncols];
  }

  //! Construct a Matrix with indices in the range r1 (first dimension)
  //! and r2 (second dimension)
  Matrix( Range r1, Range r2 )
      : mNRows (r1.lastIndex-r1.firstIndex+1),
      mNCols(r2.lastIndex-r2.firstIndex+1),
      mRowsFirstIndex(r1.firstIndex),
      mColsFirstIndex(r2.firstIndex),
      mData(NULL)
  {
    if ( mNRows > 0 && mNCols > 0 )
      mData = new T[mNRows * mNCols];
  }

  //! Copy constructor; the data is copied - there is no refcounting or somesuch.
  Matrix(const Matrix<T>& m)
  {
    *this = m;
  }

  ~Matrix()
  {
    delete[] mData;
  }

  //! free the memory held by this matrix. Simply calls resize(0,0)
  void free()
  {
    resize( 0,0 );
  }

  //! resize this matrix to be nrows x ncols; the old data is deleted.
  void resize( unsigned int nrows, unsigned int ncols )
  {
    mNRows = nrows;
    mNCols = ncols;
    mRowsFirstIndex = 0;
    mColsFirstIndex = 0;
    delete[] mData;
    if ( nrows == 0 && ncols == 0 )
      mData = NULL;
    else
      mData = new T[nrows * ncols];
  }

  //! see above
  void resize( Range r1, Range r2 )
  {
    mNRows=r1.lastIndex-r1.firstIndex+1;
    mNCols=r2.lastIndex-r2.firstIndex+1;
    mRowsFirstIndex=r1.firstIndex;
    mColsFirstIndex=r2.firstIndex;
    delete[] mData;
    mData = new T[mNCols*mNRows];
  }

  //! change the base index to i and j, respectively.
  void reindexSelf( unsigned int i, unsigned int j )
  {
    mRowsFirstIndex = i;
    mColsFirstIndex = j;
  }

  /*! \class BoundsViolation
   *  \brief exception thrown when accessing an index out of bounds
   */
  class BoundsViolation {};

  //! return a reference to the element at (row, col)
  T&       operator() (unsigned int row, unsigned int col)
  {
    //if (row >= mNRows+mRowsFirstIndex || col >= mNCols+mColsFirstIndex)
    //  throw BoundsViolation();
    return mData[(row-mRowsFirstIndex)*mNCols + col-mColsFirstIndex];
  }

  //! return a reference to the element at (row, col)
  const T& operator() (unsigned int  row, unsigned int  col) const
  {
    //if (row >= mNRows+mRowsFirstIndex || col >= mNCols+mColsFirstIndex)
    //  throw BoundsViolation();
    return mData[(row-mRowsFirstIndex)*mNCols + col-mColsFirstIndex];
  }

  //! set all elements in the matrix to t
  const T& operator= (const T& t)
  {
    std::fill( mData, mData+mNRows*mNCols, t );
    return t;
  }

  //! delete the own data and copy the data from m
  //! resizes to the shape of m
  Matrix<T>& operator= (const Matrix<T>& m)
  {
    if ( this == &m )
      return *this;
    delete[] mData;
    mRowsFirstIndex = m.rowsFirstIndex();
    mColsFirstIndex = m.colsFirstIndex();
    mNRows = m.rows();
    mNCols = m.cols();
    mData = new T[m.size()];
    std::copy( m.data(), m.data()+m.size(), mData );
    return *this;
  }

  //! elementwise addition of matrices
  Matrix<T> operator+ (const Matrix<T>& m)
  {
    Matrix<T> out( mNRows, mNCols );
    std::transform( mData, mData+mNRows*mNCols, m.data(), out.data(), std::plus<T>() );
    return out;
  }

  T* data() const { return mData; }                                  //!< return a pointer to the beginning of the data of the matrix
  unsigned int size() const { return mNCols*mNRows; }                //!< return the size of the matrix
  unsigned int rows() const { return mNRows; }                       //!< return number of rows
  unsigned int cols() const { return mNCols; }                       //!< return number of columns
  unsigned int rowsFirstIndex() const { return mRowsFirstIndex; }    //!< return the starting index of the rows
  unsigned int colsFirstIndex() const { return mColsFirstIndex; }    //!< return the starting index of the columns
};



/*! \class Array3D
 *  \brief a three dimensional array
 *
 *   A three dimensional array class to make the C++ port of the
 *   WMAP-likelihood look a bit more like the Fortran90 original.
 *
 *   An array can have differenct base indices for each dimension.
 *
 *   \warning There is no error checking, and this class is supposed
 *   to only be used for the WMAP 3-year likelihood port.
 */

/*! \enum Dim
 *  \brief the three possible dimensions of a three dimensional array
 */
enum Dim {firstDim, secondDim, thirdDim};
template<typename T>
class Array3D
{

protected:
  unsigned int  mNDim1, mNDim2, mNDim3,
  mDim1FirstIndex, mDim2FirstIndex,
  mDim3FirstIndex;
  T* mData;
public:
  //! construc an 3d array with given size
  Array3D( unsigned int nDim1 = 0, unsigned int nDim2 = 0, unsigned int nDim3 = 0 )
      : mNDim1 (nDim1), mNDim2(nDim2), mNDim3(nDim3),
      mDim1FirstIndex(0), mDim2FirstIndex(0),
      mDim3FirstIndex(0), mData(NULL)
  {
    if ( nDim1 > 0 && nDim2 > 0 && nDim3 > 0 )
      mData = new T[nDim1 * nDim2 * nDim3];
  }

  Array3D( Range r1, Range r2, Range r3 )
      : mNDim1 (r1.lastIndex-r1.firstIndex+1), mNDim2(r2.lastIndex-r2.firstIndex+1),
        mNDim3(r3.lastIndex-r3.firstIndex+1),
        mDim1FirstIndex(r1.firstIndex), mDim2FirstIndex(r2.firstIndex),
        mDim3FirstIndex(r3.firstIndex),
        mData(NULL)
  {
    if ( mNDim1 > 0 && mNDim2 > 0 && mNDim3 > 0 )
      mData = new T[mNDim1 * mNDim2 * mNDim3];
  }

  //! Copy constructor; the data is copied - there is no refcounting or somesuch.
  Array3D(const Array3D<T>& m) { *this = m; }

  //! delete the own data and copy the data from m
  //! resizes to the shape of m
  Array3D<T>& operator= (const Array3D<T>& m)
  {
    if ( this == &m )
      return *this;
    delete[] mData;

    mNDim1 = m.dim1Size();
    mNDim2 = m.dim2Size();
    mNDim3 = m.dim3Size();

    mDim1FirstIndex = m.dim1FirstIndex();
    mDim2FirstIndex = m.dim2FirstIndex();
    mDim3FirstIndex = m.dim3FirstIndex();

    mData = new T[m.size()];
    std::copy( m.data(), m.data()+m.size(), mData );
    return *this;
  }

  ~Array3D()
  {
    delete[] mData;
  }

  //! change the base indices of each dimension to a, b and c
  void reindexSelf( unsigned int a, unsigned int b, unsigned int c )
  {
    mDim1FirstIndex = a;
    mDim2FirstIndex = b;
    mDim3FirstIndex = c;
  }

  unsigned int dim1FirstIndex() const { return mDim1FirstIndex; }  //!< return base index for first dimension
  unsigned int dim2FirstIndex() const { return mDim2FirstIndex; }  //!< return base index for second dimension
  unsigned int dim3FirstIndex() const { return mDim3FirstIndex; }  //!< return base index for third dimension

  unsigned int dim1Size() const { return mNDim1; }  //!< return length in first dimension
  unsigned int dim2Size() const { return mNDim2; }  //!< return length in second dimension
  unsigned int dim3Size() const { return mNDim3; }  //!< return length in third dimension

  //! simply calls resize(0,0,0);
  void free() { resize(0,0,0); }

  //! resize the array. The old data is deleted.
  void resize( unsigned int nDim1, unsigned int nDim2, unsigned int nDim3 )
  {
    mNDim1 = nDim1;
    mNDim2 = nDim2;
    mNDim3 = nDim3;
    mDim1FirstIndex = 0;
    mDim2FirstIndex = 0;
    mDim3FirstIndex = 0;
    delete[] mData;
    if ( nDim1 == 0 && nDim2 == 0 && nDim3 == 0 )
      mData = NULL;
    else
      mData = new T[nDim1 * nDim2 * nDim3];
  }

  //! resize the array. The old data is deleted.
  void resize( Range r1, Range r2, Range r3 )
  {
   mNDim1 = r1.lastIndex-r1.firstIndex+1;
   mNDim2 = r2.lastIndex-r2.firstIndex+1;
   mNDim3 = r3.lastIndex-r3.firstIndex+1;
   mDim1FirstIndex = r1.firstIndex;
   mDim2FirstIndex = r2.firstIndex;
   mDim3FirstIndex = r3.firstIndex;

   if ( mNDim1 == 0 && mNDim2 == 0 && mNDim3 == 0 )
     mData = NULL;
   else
     mData = new T[mNDim1 * mNDim2 * mNDim3];
  }

  /*! \class BoundsViolation
   *  \brief exception thrown when accessing an index out of bounds
   */
  class BoundsViolation {};

  //! return a reference to the element at (i,l,k)
  T&       operator() (unsigned int i, unsigned int j, unsigned int k)
  {
    //X                     if (i >= mNDim1+mDim1FirstIndex || j >= mNDim2+mDim2FirstIndex
    //X                         || k >= mNDim3+mDim3FirstIndex )
    //X                       throw BoundsViolation();
    return mData[i-mDim1FirstIndex + ( j-mDim2FirstIndex )*mNDim1 + (k-mDim3FirstIndex )*mNDim1*mNDim2 ];
  }

  //! return a reference to the element at (i,l,k)
  const T&       operator() (unsigned int i, unsigned int j, unsigned int k) const
  {
    //X                     if (i >= mNDim1+mDim1FirstIndex || j >= mNDim2+mDim2FirstIndex
    //X                         || k >= mNDim3+mDim3FirstIndex )
    //X                       throw BoundsViolation();
    return mData[i-mDim1FirstIndex + ( j-mDim2FirstIndex )*mNDim1 + (k-mDim3FirstIndex )*mNDim1*mNDim2 ];
  }

  //!  assign t to all elements in the array
  const T& operator= (const T& t)
  {
    std::fill( mData, mData+mNDim1*mNDim2*mNDim3, t );
    return t;
  }

  //! return a pointer to the beginning of the data held by the array
  T* data() const { return mData; }
  unsigned int size() const { return mNDim1*mNDim2*mNDim3; }  //!< return the size of the array

  //! returns the stride of dimension d
  unsigned int stride( Dim d ) const
  {
    switch ( d )
    {
      case firstDim: return ( &((*this)(1,0,0 ))-&((*this )( 0,0,0 )) );
      case secondDim: return ( &((*this)(0,1,0 ))-&((*this )( 0,0,0 )) );
      case thirdDim: return ( &((*this)(0,0,1 ))-&((*this )( 0,0,0 )) );
    }
    return 0; //never happens
  }

};


//! helper class to convert float arrays to double arrays
//! (could be handled a lot more elegantly, I know)
class FloatToDoubleArray3D: public Array3D<double>
{
public:
  FloatToDoubleArray3D& fromFloat (const Array3D<float>& m)
  {
    delete[] mData;

    mNDim1 = m.dim1Size();
    mNDim2 = m.dim2Size();
    mNDim3 = m.dim3Size();

    mDim1FirstIndex = m.dim1FirstIndex();
    mDim2FirstIndex = m.dim2FirstIndex();
    mDim3FirstIndex = m.dim3FirstIndex();

    mData = new double[m.size()];
    std::copy( m.data(), m.data()+m.size(), mData );
    return *this;
  }
};

#endif // PORTING_HELPER_H
