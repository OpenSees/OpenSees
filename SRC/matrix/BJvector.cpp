                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/BJvector.cpp,v $
                                                                        
                                                                        
//############################################################################
//#                                                                          #
//#             /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/~~\              #
//#            |                                          |____|             #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                 B A S E                  |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |              C L A S S E S               |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |          C + +     S O U R C E           |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#            |                                          |                  #
//#         /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/   |                  #
//#        |                                         |    |                  #
//#         \_________________________________________\__/                   #
//#                                                                          #
//#                                                                          #
//############################################################################
//
//   "C makes it easy to shoot yourself in the foot, C++ makes it harder,
//   but when you do, it blows away your whole leg" -- Bjarne Stroustrup
//
///*
//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:             Vector class                                              #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic ( with help from ARKoenig in JOOP )         #
//# PROGRAMMER(S):     Boris Jeremic ( with help from ARKoenig in JOOP )         #
//#                                                                              #
//#                                                                              #
//# DATE:              Nov. 14. 1992.                                            #
//# UPDATE HISTORY:    05 - __ avgust '93.  derived class from BJmatrix class      #
//#                                         which is derived from                #
//#                                         nDarray class                        #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

// All of this inheritance idioms are after
// Jim Coplien : "Advanced C++ programing styles and idioms".

#ifndef VECTOR_CC
#define VECTOR_CC

//  #include "basics.hh"
//  #include "ndarray.hh"
#include "BJvector.h"

//##############################################################################
BJvector::BJvector(int order_n, double initvalue):
  BJmatrix( 2, order_n, 1, initvalue)  {  }  // default constructor
// rank 2 ^ just to be consistent with rank of BJmatrix
//##############################################################################
BJvector::BJvector(int order_n, double *initval):
  BJmatrix(2, order_n, 1, initval)  {  }
//rank 2 ^ just to be consistent with rank of BJmatrix

//##############################################################################
BJvector::BJvector( const nDarray & x):
  BJmatrix( x )   {  } // copy-initializer



//.... // IT IS NOT INHERITED so must be defined in all derived classes
//.... // See ARM page 277.
//.... //##############################################################################
//.... BJvector::~BJvector()
//.... {
//....   if (reference_count(-1) == 0)  // if reference count  goes to 0
//....     {
//.... // DEallocate memory of the actual nDarray
//.... //    delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
//.... // nema potrebe za brojem clanova koji se brisu## see ELLIS & STROUSTRUP $18.3
//.... //                                                and note on the p.65($5.3.4)
//.... //  and the page 276 ($12.4)
//....     delete [] data();
//....     delete [] dim();
//....     delete pc_nDarray_rep;
//....   }
//.... }




//#############################################################################
BJvector& BJvector::operator=( const BJvector & rval)
{
    rval.pc_nDarray_rep->n++; // we're adding another reference.
//    rval.reference_count(+1);  // tell the rval it has another reference
//   /*  It is important to increment the reference_counter in the new
//       BJtensor before decrementing the reference_counter in the
//       old BJtensor_rep to ensure proper operation when assigning a
//       BJtensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */
// clean up current value;
    if( reference_count(-1) == 0)  // if nobody else is referencing us.
      {
        delete [] data();
        delete [] dim();
        delete pc_nDarray_rep;
      }
 // connect to new value
    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval nDarray_rep
    return *this;
}

//..//#############################################################################
//..BJvector& BJvector::operator=( const BJmatrix & rval)
//..{
//..    rval.pc_nDarray_rep->n++; // we're adding another reference.
//..//    rval.reference_count(+1);  // tell the rval it has another reference
//..//   /*  It is important to increment the reference_counter in the new
//..//       BJtensor before decrementing the reference_counter in the
//..//       old BJtensor_rep to ensure proper operation when assigning a
//..//       BJtensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */
//..
//.. // clean up current value;
//..    if( reference_count(-1) == 0)  // if nobody else is referencing us.
//..      {
//..        delete [] data();
//..        delete [] dim();
//..        delete pc_nDarray_rep;
//..      }
//..
//..// set back rank to 1 for BJvector instead of 2 as in BJmatrix case
//..      rval.pc_nDarray_rep->nDarray_rank = 1;
//..//    rval.rank(1);
//..// connect to new value
//..    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval nDarray_rep
//..    return *this;
//..}
//..
//..//#############################################################################
//..BJvector& BJvector::operator=( const nDarray & rval)
//..{
//..    rval.pc_nDarray_rep->n++; // we're adding another reference.
//..//    rval.reference_count(+1);  // tell the rval it has another reference
//..//   /*  It is important to increment the reference_counter in the new
//..//       BJtensor before decrementing the reference_counter in the
//..//       old BJtensor_rep to ensure proper operation when assigning a
//..//       BJtensor_rep to itself ( after ARKoenig JOOP May/June '90 )  */
//..
//.. // clean up current value;
//..    if( reference_count(-1) == 0)  // if nobody else is referencing us.
//..      {
//..        delete [] data();
//..        delete [] dim();
//..        delete pc_nDarray_rep;
//..      }
//..
//.. // connect to new value
//..    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval nDarray_rep
//..    return *this;
//..}
//..

//#######  //#############################################################################
//#######  // I had to overload the operator* from BJmatrix class
//#######  // because in the case of BJvectors s*sT and so on
//#######  BJmatrix BJvector::operator*( BJvector & arg)
//#######    {
//#######  //    if( cols() != arg.rows())
//#######  //      error("# rows of second mat must equal "
//#######  //               "# cols of first for multiply#");
//#######      BJmatrix result(rows(),arg.cols());
//#######      for( int row=0 ; row<rows() ; row++ )
//#######        for( int col=0 ; col<arg.cols() ; col++ )
//#######          {
//#######            double sum = 0;
//#######            for( int i=0 ; i<cols() ; i++ )
//#######              sum += mval(row,i)*arg.mval(i,col);
//#######            result.mval(row,col) = sum;
//#######          }
//#######      return result; // Returning a local variable?
//#######      // copy-initializer happens before the destructor,
//#######      // so reference count is 2 when destructor is called,
//#######      // thus destructor doesn't free the memory.
//#######    }
//#######  
//....//#############################################################################
//....BJvector BJvector::operator*( double arg)
//....  {
//....    BJvector result(rows());
//....    for ( int i=0 ; i<rows() ; i++ )
//....      result.val(i) = cval(i) * arg;
//....    return result;
//....  }



//##############################################################################
//
// CODE BOUND checking routine ( slower but safer )
//
double BJvector::cval(int subscript, ... )  const
  {
    return (this->BJmatrix::cval(subscript,1));
  }
double & BJvector::val(int subscript, ... )
  {
    return (this->BJmatrix::val(subscript,1));
  }

#endif
