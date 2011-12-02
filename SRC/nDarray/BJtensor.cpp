                                                                        
// $Revision: 1.3 $                                                              
// $Date: 2004-07-26 14:26:10 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/nDarray/BJtensor.cpp,v $                                                                
                                                                        
                                                                        
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
//# CLASS:             BJtensor                                                    #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              May 28. - July 21 '93                                     #
//# UPDATE HISTORY:    july 8. '93. BJtensor02 - BJtensor multiplication             #
//#                                 inner and outer products                     #
//#                    august 17-19 '93 fixed default constructor that wasted    #
//#                                     memory ###                               #
//#                    october 11 '93 added transpose0110, transpose0101,        #
//#                                   transpose0111 so the creation of           #
//#                                   isotropic BJtensor is much easer and         #
//#                                   understandable !                           #
//#                    januar 06 '93  added BJtensor2BJmatrix_1, BJtensor2BJmatrix_2     #
//#                                   inverse_1, inverse_2, inverse_3            #
//#                    januar 20 '93  added inverse  TRUE ONE                    #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                    28June2004     added val4 for efficiency still            #
//#                                   to be worked on                            #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
//
#ifndef TENSOR_CC
#define TENSOR_CC

//  #include "basics.h"
//  #include "nDarray.h"
#include "BJtensor.h"

// just send appropriate arguments to the base constructor
//##############################################################################

BJtensor::BJtensor(int rank_of_BJtensor, double initval):
  nDarray(rank_of_BJtensor, initval)  // default constructor
    {
      indices1 = (char *) NULL; //since NULL is (void *)  02june98
      indices2 = (char *) NULL; //since NULL is (void *)  02june98 
    }


//##############################################################################
BJtensor::BJtensor(int rank_of_BJtensor, const int *pdim, double *values):
  nDarray(rank_of_BJtensor, pdim, values)
    {
      indices1 = (char *) NULL; //since NULL is (void *)  02june98 
      indices2 = (char *) NULL; //since NULL is (void *)  02june98 
    }


//##############################################################################
BJtensor::BJtensor(int rank_of_BJtensor, const int *pdim, double initvalue):
  nDarray(rank_of_BJtensor, pdim, initvalue)
    {
      indices1 = (char *) NULL; //since NULL is (void *)  02june98 
      indices2 = (char *) NULL; //since NULL is (void *)  02june98 
    }


//##############################################################################
BJtensor::BJtensor(char *flag, int rank_of_BJtensor, const int *pdim):
  nDarray( flag, rank_of_BJtensor, pdim) // create a unit nDarray
    {
      indices1 = (char *) NULL; //since NULL is (void *)  02june98 
      indices2 = (char *) NULL; //since NULL is (void *)  02june98 
    }


//##############################################################################
BJtensor::BJtensor(char *flag):
  nDarray(flag) { } //this one used to send "NO" message



//##############################################################################
BJtensor::BJtensor(const BJtensor & x):   // copy initializer
  nDarray("NO")           // with base class constructor cancelation
  {
    x.pc_nDarray_rep->n++;  // we're adding another reference.
//    x.reference_count(+1); // we're adding another reference.
    pc_nDarray_rep = x.pc_nDarray_rep;  // point to the new BJtensor_rep.
// add the indices
    indices1 = x.indices1;
    indices2 = x.indices2;
 }



//##############################################################################
BJtensor::BJtensor(const nDarray & x):
  nDarray( x ) {  } // copy-initializer
                     // DO NOT delete ( NULL ) indices because return
                     // from operator*( BJtensor & arg ) this one
                     // is invoked and so I need this indices.




 // IT IS NOT INHERITED so must be defined in all derived classes
 // See ARM page 277.
 //##############################################################################
// BJtensor::~BJtensor()
// {
//   if (reference_count(-1) == 0)  // if reference count  goes to 0
//     {
// // DEallocate memory of the actual nDarray
// // delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
// // nema potrebe za brojem clanova koji se brisu!! see ELLIS & STROUSTRUP $18.3
// //                                                and note on the p.65($5.3.4)
// //  and the page 276 ($12.4)
//     delete data();
//     delete dim();
//     delete pc_nDarray_rep;
//   }
// }




// IT IS NOT INHERITED so must be defined in all derived classes
// See ARM page 306.
//##############################################################################
BJtensor& BJtensor::operator=( const BJtensor & rval)
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
    pc_nDarray_rep = rval.pc_nDarray_rep;  // point at the rval BJtensor_rep
// null indices in the rval AND in the this
// because rval is temporary anyway and I need *this
//    rval.null_indices();
    this->null_indices();
    return *this;
  }




//##############################################################################
// this is supposed to fill in the char * indices1 or the char * indices2
// array in BJtensor object
// so that they can be checked for matching later on when operations like
// single contraction (.), double contraction (:), dyadic product (otimes)
// are performed the object can choose the right operator.
// Also when multiplying BJtensor by itself ( like for example:
//               D("ij")*D("kl")
// the other char * indices2 is defined to take the other array of
// indices
// so that mulitplication will be O.K.
// Since only two BJtensors can be multiplied at the time ( binary operation )
// only indices1 and indices2 are needed #
// WATCH OUT THIS IS NOT STANDRAD AS YOU CANNOT QUARANTY THE ORDER OF EXECUTION!!!!
BJtensor & BJtensor::operator()(char *indices_from_user)
  {
    if ( this->indices1 == NULL )
      {
        this->indices1 = indices_from_user;
        return *this;
      }
    else
      {
        this->indices2 = indices_from_user;
        return *this;
      }
  }



//##############################################################################
void BJtensor::null_indices()
  {
    this->indices1 = (char *) NULL; //since NULL is (void *)  02june98 
    this->indices2 = (char *) NULL; //since NULL is (void *)  02june98 
  }




//++//##############################################################################
//++// BJtensor addition
//++BJtensor BJtensor::operator+( BJtensor & rval)
//++  {
//++    int this_rank_of_BJtensor = this->rank();  //pc_nDarray_rep->nDarray_rank;
//++    int rval_rank_of_BJtensor =  rval.rank();  //pc_nDarray_rep->nDarray_rank;
//++
//++    if(this_rank_of_BJtensor != rval_rank_of_BJtensor)
//++      {
//++        ::printf("\a\nTensors of different ranks: addition not possible\n");
//++        ::exit ( 1 );
//++      }
//++
//++    for ( int i=0 ; i<this_rank_of_BJtensor ; i++ )
//++      if (this->dim()[i] != rval.dim()[i] )
//++//      if (this->dim(i) != rval.dim(i) )
//++        {
//++          ::fprintf(stderr,"\a\nDimension discrepancy in operator+\n",
//++                   "this->dim()[%d]=%d\n",
//++                   "arg.dim()[%d]=%d\n",
//++                    i,this->dim()[i],
//++                    i,rval.dim()[i] );
//++          ::exit(1);
//++        }
//++// construct BJtensor using the same control numbers as for the
//++// original one .
//++//      BJtensor add(pc_nDarray_rep->nDarray_rank, pc_nDarray_rep->dim, 0.0);
//++      BJtensor add( this->rank(), dim(), 0.0);
//++
//++      add.indices1 = this->indices1;
//++
//++//      switch(pc_nDarray_rep->nDarray_rank)
//++      switch(this->rank())
//++        {
//++          case 0:
//++            {
//++              add.val(1) = val(1) + rval.val(1);
//++              break;
//++            }
//++
//++          case 1:
//++            {
//++              for ( int i1=1 ; i1<=3 ; i1++ )
//++                add.val(i1) = val(i1) + rval.val(i1);
//++              break;
//++            }
//++
//++          case 2:
//++            {
//++              for ( int i2=1 ; i2<=3 ; i2++ )
//++                for ( int j2=1 ; j2<=3 ; j2++ )
//++                  add.val(i2, j2) = val(i2, j2) + rval.val(i2, j2);
//++              break;
//++            }
//++
//++          case 3:
//++            {
//++              for ( int i3=1 ; i3<=3 ; i3++ )
//++                for ( int j3=1 ; j3<=3 ; j3++ )
//++                  for ( int k3=1 ; k3<=3 ; k3++ )
//++                    add.val(i3, j3, k3) = val(i3, j3, k3) +
//++                                                 rval.val(i3, j3, k3);
//++              break;
//++            }
//++
//++          case 4:
//++            {
//++              for ( int i4=1 ; i4<=3 ; i4++ )
//++                for ( int j4=1 ; j4<=3 ; j4++ )
//++                  for ( int k4=1 ; k4<=3 ; k4++ )
//++                    for ( int l4=1 ; l4<=3 ; l4++ )
//++                      add.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)+
//++                                                  rval.val(i4,j4,k4,l4);
//++              break;
//++            }
//++        }
//++
//++    null_indices();
//++    rval.null_indices();
//++
//++    return add;
//++  }
//++
//````````//##############################################################################
//````````// BJtensor addition
//````````BJtensor& BJtensor::operator+=(const BJtensor & rval)
//````````  {
//````````    int this_rank_of_BJtensor = this->rank();  //pc_nDarray_rep->nDarray_rank;
//````````    int rval_rank_of_BJtensor =  rval.rank();  //pc_nDarray_rep->nDarray_rank;
//````````
//````````    if(this_rank_of_BJtensor != rval_rank_of_BJtensor)
//````````      {
//````````        ::printf("\a\nTensors of different ranks: += not possible\n");
//````````        ::exit ( 1 );
//````````      }
//````````
//````````    for ( int i=0 ; i<this_rank_of_BJtensor ; i++ )
//````````      if (this->dim()[i] != rval.dim()[i] )
//````````//      if (this->dim(i) != rval.dim(i) )
//````````        {
//````````::fprintf(stderr,"\a\nDimension discrepancy in operator+=\n this->dim()[%d]=%d\n arg.dim()[%d]=%d\n",
//````````i,this->dim()[i],
//````````i,rval.dim()[i] );
//````````::exit(1);
//````````        }
//````````// Copy *this if necessary
//````````    if ( this->pc_nDarray_rep->n > 1 )// see ARK in JOOP may/june '90
//````````      {                               // "Letter From a Newcomer"
//````````//..............................................................................
//````````      // create the structure:
//````````        nDarray_rep * New_pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
//````````        New_pc_nDarray_rep->nDarray_rank = this->pc_nDarray_rep->nDarray_rank;
//````````// in the case of nDarray_rank=0 add one to get right thing from the
//````````// operator new
//````````        int one_or0 = 0;
//````````        if(!New_pc_nDarray_rep->nDarray_rank) one_or0 = 1;
//````````        New_pc_nDarray_rep->dim = new int[New_pc_nDarray_rep->nDarray_rank+one_or0];
//````````                                  // array for dimensions
//````````        New_pc_nDarray_rep->total_numb = 1;
//````````        for( int idim = 0 ; idim < New_pc_nDarray_rep->nDarray_rank ; idim++ )
//````````          {
//````````            New_pc_nDarray_rep->dim[idim] = this->pc_nDarray_rep->dim[idim];
//````````            New_pc_nDarray_rep->total_numb *= New_pc_nDarray_rep->dim[idim];
//````````          }
//````````// allocate memory for the actual nDarray as nDarray
//````````        New_pc_nDarray_rep->pd_nDdata = new double [(size_t)New_pc_nDarray_rep->total_numb];
//````````          if (!New_pc_nDarray_rep->pd_nDdata)
//````````            {
//````````              ::fprintf(stderr,"\a\nInsufficient memory for array\n");
//````````              ::exit(1);
//````````            }
//````````         New_pc_nDarray_rep->n = 1;  // so far, there's one reference
//````````         for ( i=0 ; i<New_pc_nDarray_rep->total_numb ; i++ )
//````````           New_pc_nDarray_rep->pd_nDdata[i] = this->pc_nDarray_rep->pd_nDdata[i];
//````````//.........
//````````         this->pc_nDarray_rep->total_numb--;
//````````         this->pc_nDarray_rep = New_pc_nDarray_rep;
//````````//..............................................................................
//````````      }
//````````// it appears that I can add this two BJtensors just as a simple BJvectors:
//````````    for (int j=0 ; j<this->pc_nDarray_rep->total_numb ; j++)
//````````      this->pc_nDarray_rep->pd_nDdata[j] += rval.pc_nDarray_rep->pd_nDdata[j];
//````````
//````````    this->indices1 = rval.indices1;
//````````    return *this;
//````````  }
//````````

//##############################################################################
// BJtensor addition
BJtensor operator+(const BJtensor & lval, const BJtensor & rval)
  {
    BJtensor result(lval);
    result += rval;
    return result;
  }




////##############################################################################
//// scalar addition
//BJtensor  BJtensor::operator+(double rval)
// {
//// construct BJtensor using the same control numbers as for the
//// original one.
//      BJtensor add(this->rank(), dim(), 0.0);
//
//      add.indices1 = this->indices1;
//
//      switch(this->rank())
//        {
//          case 0:
//            {
//              add.val(1) = val(1) + rval;
//              break;
//            }
//
//          case 1:
//            {
//              for ( int i1=1 ; i1<=this->dim()[0] ; i1++ )
//                add.val(i1) = val(i1) + rval;
//              break;
//            }
//
//          case 2:
//            {
//              for ( int i2=1 ; i2<=this->dim()[0] ; i2++ )
//                for ( int j2=1 ; j2<=this->dim()[1] ; j2++ )
//                  add.val(i2, j2) = val(i2, j2) + rval;
//              break;
//            }
//
//          case 3:
//            {
//              for ( int i3=1 ; i3<=this->dim()[0] ; i3++ )
//                for ( int j3=1 ; j3<=this->dim()[1] ; j3++ )
//                  for ( int k3=1 ; k3<=this->dim()[2] ; k3++ )
//                    add.val(i3, j3, k3) = val(i3, j3, k3) + rval;
//              break;
//            }
//
//          case 4:
//            {
//              for ( int i4=1 ; i4<=this->dim()[0] ; i4++ )
//                for ( int j4=1 ; j4<=this->dim()[1] ; j4++ )
//                  for ( int k4=1 ; k4<=this->dim()[2] ; k4++ )
//                    for ( int l4=1 ; l4<=this->dim()[3] ; l4++ )
//                      add.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)+rval;
//              break;
//            }
//        }
//
//    null_indices();
//
//    return add;
// }



//..//++//##############################################################################
//..//++// BJtensor substraction
//..//++BJtensor BJtensor::operator-(BJtensor & rval)
//++ {
//++    int this_rank_of_BJtensor = this->rank();  //pc_nDarray_rep->nDarray_rank;
//++    int rval_rank_of_BJtensor = rval.rank();   //pc_nDarray_rep->nDarray_rank;
//++
//++    if(this_rank_of_BJtensor != rval_rank_of_BJtensor)
//++      {
//++        ::printf("\a\nTensors of different ranks:",
//++                 " substraction not possible\n");
//++        ::exit ( 1 );
//++      }
//++
//++    for ( int i=0 ; i<this_rank_of_BJtensor ; i++ )
//++      if (this->dim()[i] != rval.dim()[i] )
//++        {
//++          ::fprintf(stderr,"\a\nDimension discrepancy in operator  -  \n",
//++                   "this->dim()[%d]=%d\n",
//++                   "arg.dim()[%d]=%d\n",
//++                    i,this->dim()[i],
//++                    i,rval.dim()[i]);
//++          ::exit(1);
//++        }
//++// construct BJtensor using the same control numbers as for the
//++// original one.
//++      BJtensor sub(this->rank(), dim(), 0.0);
//++
//++      sub.indices1 = this->indices1;
//++
//++      switch(this->rank())
//++        {
//++          case 0:
//++            {
//++              sub.val(1) = val(1) - rval.val(1);
//++              break;
//++            }
//++
//++          case 1:
//++            {
//++              for ( int i1=1 ; i1<=3 ; i1++ )
//++                sub.val(i1) = val(i1) - rval.val(i1);
//++              break;
//++            }
//++
//++          case 2:
//++            {
//++              for ( int i2=1 ; i2<=3 ; i2++ )
//++                for ( int j2=1 ; j2<=3 ; j2++ )
//++                  sub.val(i2, j2) = val(i2, j2) -
//++                                           rval.val(i2, j2);
//++              break;
//++            }
//++
//++          case 3:
//++            {
//++              for ( int i3=1 ; i3<=3 ; i3++ )
//++                for ( int j3=1 ; j3<=3 ; j3++ )
//++                  for ( int k3=1 ; k3<=3 ; k3++ )
//++                    sub.val(i3, j3, k3) = val(i3, j3, k3) -
//++                                                 rval.val(i3, j3, k3);
//++              break;
//++            }
//++
//++          case 4:
//++            {
//++              for ( int i4=1 ; i4<=3 ; i4++ )
//++                for ( int j4=1 ; j4<=3 ; j4++ )
//++                  for ( int k4=1 ; k4<=3 ; k4++ )
//++                    for ( int l4=1 ; l4<=3 ; l4++ )
//++                      sub.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)-
//++                                                  rval.val(i4,j4,k4,l4);
//++              break;
//++            }
//++        }
//++
//++    null_indices();
//++    rval.null_indices();
//++
//++    return sub;
//++ }

//````````//##############################################################################
//````````// BJtensor substraction
//````````BJtensor& BJtensor::operator-=(const BJtensor & rval)
//````````  {
//````````    int this_rank_of_BJtensor = this->rank();  //pc_nDarray_rep->nDarray_rank;
//````````    int rval_rank_of_BJtensor =  rval.rank();  //pc_nDarray_rep->nDarray_rank;
//````````
//````````    if(this_rank_of_BJtensor != rval_rank_of_BJtensor)
//````````      {
//````````        ::printf("\a\nTensors of different ranks: -= not possible\n");
//````````        ::exit ( 1 );
//````````      }
//````````
//````````    for ( int i=0 ; i<this_rank_of_BJtensor ; i++ )
//````````      if (this->dim()[i] != rval.dim()[i] )
//````````//      if (this->dim(i) != rval.dim(i) )
//````````        {
//````````::fprintf(stderr,"\a\nDimension discrepancy in operator-=\n this->dim()[%d]=%d\n arg.dim()[%d]=%d\n",
//````````i,this->dim()[i],
//````````i,rval.dim()[i] );
//````````::exit(1);
//````````        }
//````````// Copy *this if necessary
//````````    if ( this->pc_nDarray_rep->n > 1 )// see ARK in JOOP may/june '90
//````````      {                               // "Letter From a Newcomer"
//````````//..............................................................................
//````````      // create the structure:
//````````        nDarray_rep * New_pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
//````````        New_pc_nDarray_rep->nDarray_rank = this->pc_nDarray_rep->nDarray_rank;
//````````// in the case of nDarray_rank=0 add one to get right thing from the
//````````// operator new
//````````        int one_or0 = 0;
//````````        if(!New_pc_nDarray_rep->nDarray_rank) one_or0 = 1;
//````````        New_pc_nDarray_rep->dim = new int[New_pc_nDarray_rep->nDarray_rank+one_or0];
//````````                                  // array for dimensions
//````````        New_pc_nDarray_rep->total_numb = 1;
//````````        for( int idim = 0 ; idim < New_pc_nDarray_rep->nDarray_rank ; idim++ )
//````````          {
//````````            New_pc_nDarray_rep->dim[idim] = this->pc_nDarray_rep->dim[idim];
//````````            New_pc_nDarray_rep->total_numb *= New_pc_nDarray_rep->dim[idim];
//````````          }
//````````// allocate memory for the actual nDarray as nDarray
//````````        New_pc_nDarray_rep->pd_nDdata = new double [(size_t)New_pc_nDarray_rep->total_numb];
//````````          if (!New_pc_nDarray_rep->pd_nDdata)
//````````            {
//````````              ::fprintf(stderr,"\a\nInsufficient memory for array\n");
//````````              ::exit(1);
//````````            }
//````````         New_pc_nDarray_rep->n = 1;  // so far, there's one reference
//````````         for ( i=0 ; i<New_pc_nDarray_rep->total_numb ; i++ )
//````````           New_pc_nDarray_rep->pd_nDdata[i] = this->pc_nDarray_rep->pd_nDdata[i];
//````````//.........
//````````         this->pc_nDarray_rep->total_numb--;
//````````         this->pc_nDarray_rep = New_pc_nDarray_rep;
//````````//..............................................................................
//````````      }
//````````// it appears that I can add this two BJtensors just as a simple BJvectors:
//````````    for (int j=0 ; j<this->pc_nDarray_rep->total_numb ; j++)
//````````      this->pc_nDarray_rep->pd_nDdata[j] -= rval.pc_nDarray_rep->pd_nDdata[j];
//````````
//````````    this->indices1 = rval.indices1;
//````````    return *this;
//````````  }
//````````

//##############################################################################
// BJtensor substraction
BJtensor operator-(const BJtensor & lval, const BJtensor & rval)
  {
    BJtensor result(lval);
    result -= rval;
    return result;
  }




////##############################################################################
//// scalar substraction
//BJtensor  BJtensor::operator-( double rval)
//  {
//// construct BJtensor using the same control numbers as for the
//// original one.
//    BJtensor sub(this->rank(), dim(), 0.0);
//
//    sub.indices1 = this->indices1;
//
//    switch(this->rank())
//      {
//        case 0:
//          {
//            sub.val(1) = val(1) - rval;
//            break;
//          }
//
//        case 1:
//          {
//            for ( int i1=1 ; i1<=this->dim()[0] ; i1++ )
//              sub.val(i1) = val(i1) - rval;
//            break;
//          }
//
//        case 2:
//          {
//            for ( int i2=1 ; i2<=this->dim()[0] ; i2++ )
//              for ( int j2=1 ; j2<=this->dim()[1] ; j2++ )
//                sub.val(i2, j2) = val(i2, j2) - rval;
//            break;
//          }
//
//        case 3:
//          {
//            for ( int i3=1 ; i3<=this->dim()[0] ; i3++ )
//              for ( int j3=1 ; j3<=this->dim()[1] ; j3++ )
//                for ( int k3=1 ; k3<=this->dim()[2] ; k3++ )
//                  sub.val(i3, j3, k3) = val(i3, j3, k3) - rval;
//            break;
//          }
//
//        case 4:
//          {
//            for ( int i4=1 ; i4<=this->dim()[0] ; i4++ )
//              for ( int j4=1 ; j4<=this->dim()[1] ; j4++ )
//                for ( int k4=1 ; k4<=this->dim()[2] ; k4++ )
//                  for ( int l4=1 ; l4<=this->dim()[3] ; l4++ )
//                    sub.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)-rval;
//            break;
//          }
//      }
//
//    null_indices();
//
//    return sub;
//  }
//

//##############################################################################
//// scalar multiplication (BJtensor * scalar)
//check what happens if you are have tenosr A=B and then you do B*5.13
// does it change A as well !!!!!!!!!!!!!!!!!!!!!!!!!!
// TS KR!!!!
//
BJtensor BJtensor::operator*( const double rval)
 {
// construct BJtensor using the same control numbers as for the
// original one.
    BJtensor mult(this->rank(), dim(), 0.0);

    mult.indices1 = this->indices1;

    switch(this->rank())
      {
        case 0:
          {
            mult.val(1) = val(1) * rval;
            break;
          }

        case 1:
          {
            for ( int i1=1 ; i1<=this->dim()[0] ; i1++ )
              mult.val(i1) = val(i1) * rval;
            break;
          }

        case 2:
          {
            for ( int i2=1 ; i2<=this->dim()[0] ; i2++ )
              for ( int j2=1 ; j2<=this->dim()[1] ; j2++ )
                mult.val(i2, j2) = val(i2, j2) * rval;
            break;
          }

        case 3:
          {
            for ( int i3=1 ; i3<=this->dim()[0] ; i3++ )
              for ( int j3=1 ; j3<=this->dim()[1] ; j3++ )
                for ( int k3=1 ; k3<=this->dim()[2] ; k3++ )
                  mult.val(i3, j3, k3) = val(i3, j3, k3) * rval;
            break;
          }

        case 4:
          {
            for ( int i4=1 ; i4<=this->dim()[0] ; i4++ )
              for ( int j4=1 ; j4<=this->dim()[1] ; j4++ )
                for ( int k4=1 ; k4<=this->dim()[2] ; k4++ )
                  for ( int l4=1 ; l4<=this->dim()[3] ; l4++ )
                    mult.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)*rval;
            break;
          }
      }

    null_indices();

    return mult;
 }

//    BJtensor operator*( double lval, nDarray & rval);  // REVIEWER global
//##############################################################################
// scalar multiplication
BJtensor  operator*( const double lval, BJtensor & rval)
  {
    return rval*lval;
  }

//##############################################################################
//##############################################################################
//##############################################################################
//##############################################################################
//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
//
//                                _____..---========+^+=========---.._____
//   ______________________ __,_='=====____  ================== _____=====`=
//  (._____________________I__) _ __=_/    `--------=+=--------'
//      /      /--...---===='---+----'
//     '------'---.--- _  _ =   _._'    "Make it so..."
//                    `--------'                  Captain Jean-Luc Picard
//                                                USS Enterprise, NCC-1701D
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
BJtensor BJtensor::operator*( BJtensor & arg)
  {
//DEBUGprint::printf("before block in operator* coreleft is: %lu bytes\n", (unsigned long) ::coreleft());
//   const int MAX_TENS_ORD = 8;
   const int MAX_TENS_ORD = 4;
//....   int results_rank = 0;

//DEBUGprint::printf("this->rank()=%d\n",
//DEBUGprint          this->rank());
//DEBUGprint::printf("arg.rank()=%d\n",
//DEBUGprint          arg.rank());

// space for CONTRACTED indices
   int *this_contr = new int[MAX_TENS_ORD];
   int *arg_contr  = new int[MAX_TENS_ORD];

// space for UN-CONTRACTED indices
   int *this_uncontr = new int[MAX_TENS_ORD];
   int *arg_uncontr  = new int[MAX_TENS_ORD];

   for(int this_ic=0 ; this_ic<MAX_TENS_ORD ; this_ic++ )
     {
       this_contr[this_ic] = 0;
       this_uncontr[this_ic] = 0;
     }
   for( int arg_ic=0 ; arg_ic<MAX_TENS_ORD ; arg_ic++ )
     {
       arg_contr[arg_ic]  = 0;
       arg_uncontr[arg_ic]  = 0;
     }

// now make the right choice for indices for this and arg BJtensor !!!
//
//                  |||~          |||~          |||~
//                 (0 0)         (0 0)         (0 0)
//--------------ooO-(_)-Ooo---ooO-(_)-Ooo---ooO-(_)-Ooo----------------------
//
// WATCH OUT THIS IS NOT STANDRAD AS YOU CANNOT QUARANTY THE ORDER OF EXECUTION!!!!

   char * this_indices = this->indices1;
   char * arg_indices  = arg.indices1;

//fprintf(stdout,"\n\n\n\n  this->indices1 = %s   ; this->indices2 = %s\n", this->indices1,this->indices2);
// if the BJtensors are same then split indices :-)
   if ( this->pc_nDarray_rep == arg.pc_nDarray_rep )
     {
// WATCH OUT THIS IS NOT STANDRAD AS YOU CANNOT QUARANTY THE ORDER OF EXECUTION!!!!
       this_indices = this->indices2; // this is changed because
       arg_indices  = arg.indices1;  // it reads indices from
                                    // right so "in" indices1 are
                                   // the right ( arg ) and "in"
                                  // indices2 are the
                                 // left indices
// WATCH OUT THIS IS NOT STANDRAD AS YOU CANNOT QUARANTY THE ORDER OF EXECUTION!!!!
//       this_indices = this->indices1; // Now it turns out that it is 
//       arg_indices  = arg.indices2;  // reading indices from left!!!!!!
//                                    // Not any more as per g++ (12 may 2000
//                                   // 
//                                  // 
//                                 //  
// WATCH OUT THIS IS NOT STANDRAD AS YOU CANNOT QUARANTY THE ORDER OF EXECUTION!!!!
     }
//DEBUGprint::printf("this_indices=%s\n", this_indices);
//DEBUGprint::printf("arg_indices=%s\n",  arg_indices);


//   int this_indices_number = ::strlen(this->pc_nDarray_rep->indices1);
//   int  arg_indices_number = ::strlen(arg.pc_nDarray_rep->indices1);

   int this_indices_number = ::strlen(this_indices);
   int  arg_indices_number = ::strlen(arg_indices);

//DEBUGprint::printf("this_indices_number=%d\n",
//DEBUGprint          this_indices_number);
//DEBUGprint::printf("arg_indices_number=%d\n",
//DEBUGprint          arg_indices_number);

// check for indices number: should be EQUAL than rank of BJtensor
   if ( this_indices_number != this->rank() )
     {
::fprintf(stderr,"\a\n 'this' has more/less indices than expected:\n \
this_indices_number = %d\n \
this_indices = %s\n \
this->rank() = %d\n",
this_indices_number,
this_indices,
this->rank());
::exit(1);
     }

   if ( arg_indices_number != arg.rank() )
     {
::fprintf(stderr,"\a\n 'arg' has more/less indices than expected:\n \
arg_indices_number = %d\n \
arg_indices = %s\n \
arg.rank() = %d\n",
arg_indices_number,
arg_indices,
arg.rank());
::exit(1);
     }

// counter for contracted indices
   int contr_counter = contracted_ind(this_indices,
                                      arg_indices,
                                      this_contr,
                                      arg_contr,
                                      this_indices_number,
                                      arg_indices_number  );

   int this_uni_count = uncontracted_ind( this_uncontr,
                                          this_contr,
                                          this_indices_number);

   int arg_uni_count = uncontracted_ind( arg_uncontr,
                                         arg_contr,
                                         arg_indices_number);

//number of UNcontractions
   int uncontr_counter = this_uni_count + arg_uni_count;

//TEMP  this is just the test because right now only up to order = 4
   if ( uncontr_counter > MAX_TENS_ORD )
     {
       ::fprintf(stderr,"\a\a\n\n  OOOPS product of multiplication has order %d\n",
                           uncontr_counter);
       ::exit(1);
     }

//DEBUGprint ::printf("..... this_uni_count = %d, arg_uni_count = %d\n",
//DEBUGprint                 this_uni_count,      arg_uni_count);
//DEBUGprint ::printf(" uncontr_counter %d\n",uncontr_counter);
//DEBUGprint for( int pthis_ic=0 ;
//DEBUGprint      pthis_ic<(this->rank()+1) ;
//DEBUGprint      pthis_ic++ )
//DEBUGprint   ::printf("this_uncontr[%d]=%d\n",pthis_ic,this_uncontr[pthis_ic]);
//DEBUGprint for( int parg_ic=0 ;
//DEBUGprint      parg_ic<(arg.rank()+1) ;
//DEBUGprint      parg_ic++ )
//DEBUGprint   ::printf("             arg_uncontr[%d]=%d\n",parg_ic,arg_uncontr[parg_ic]);

/////////////////////////////////////////////////////
// let's make result BJtensor

   int t = 0;
   int a = 0;
   static int results_dims[MAX_TENS_ORD];
   for( t=0 ; t < this_uni_count ; t++ )
     results_dims[t]=this->dim()[this_uncontr[t]-1];
   for( a=0 ; a < arg_uni_count ; a++ )
     results_dims[a+t]=arg.dim()[arg_uncontr[a]-1];
// do kraja (till the end) . . .
   for( ; a < MAX_TENS_ORD-t ; a++ )
     results_dims[a+t]=1;


// Indices  pa sada indeksi ( index )
   static char results_indices[MAX_TENS_ORD+1];
   for( t=0 ; t < this_uni_count ; t++ )
     results_indices[t]=this_indices[this_uncontr[t]-1];
   for( a=0 ; a < arg_uni_count ; a++ )
     results_indices[a+t]=arg_indices[arg_uncontr[a]-1];
// the last one . . .
   results_indices[uncontr_counter] = '\0';


//DEBUGprint  ::printf("\n\n......  uncontr_counter+1=%d\n", uncontr_counter+1);
//DEBUGprintfor( int pa=0 ; pa < MAX_TENS_ORD ; pa++ )
//DEBUGprint  ::printf("   results_dims[%d]=%d\n",pa,results_dims[pa]);

   BJtensor result(uncontr_counter, results_dims, 0.0);
   result(results_indices);  // initialize indices in result
//DEBUGprint  ::printf("\n\n..................   results_indices=%s\n",results_indices);

//////////////////////////////////////////////////
// check dimensions in each BJtensor order
   for( t=0 ; t<contr_counter ; t++ )
     if ( this->dim()[this_contr[t]-1] !=
          arg.dim()[arg_contr[t]-1] )
       {
         ::fprintf(stderr,"\a\a\n t = %d \n"
                  "this_contr[t]-1  = %d\n"
                  "this->dim()[%d] = %d\n"
                  "arg_contr[t]-1  = %d\n"
                  "arg.dim()[%d] = %d\n",
                   t,
                   this_contr[t]-1,
                   this_contr[t]-1,this->dim()[this_contr[t]-1],
                   arg_contr[t]-1,
                   arg_contr[t]-1,arg.dim()[arg_contr[t]-1]);
         ::exit(1);
       }

/////////////////////////////////////////////////////#*%
//////////////////////////////////////////////////#*%
///////////////////////////////////////////////#*%
////////////////////////////////////////////#*%
/////////////////////////////////////////#*%
//////////////////////////////////////#*%
///////////////////////////////////#*%
////////////////////////////////#*%
/////////////////////////////#*%
//////////////////////////#*%
///////////////////////#*%
////////////////////#*%
/////////////////#*%
//////////////#*%
///////////#*%
////////#*%
/////#*%
//#*%

   int *inerr_dims = new int[MAX_TENS_ORD];
   for( t=0 ; t<contr_counter ; t++ )
     {
       inerr_dims[t] = this->dim()[this_contr[t]-1];
//DEBUGprint       ::printf("    inerr_dims[%d] = %d\n",t,inerr_dims[t]);
     }
   for( ; t<MAX_TENS_ORD ; t++ )
     {
       inerr_dims[t] = 1;
//DEBUGprint       ::printf("          the rest of inerr_dims[%d] = %d\n",t,inerr_dims[t]);
     }


   int *lid = new int[MAX_TENS_ORD];
   for( t=0 ; t<MAX_TENS_ORD ; t++ )
     {
       lid[t] = 1;
//DEBUGprint       ::printf("    lid[%d] = %d\n",t,lid[t]);
     }

   int *rid = new int[MAX_TENS_ORD];
   for( t=0 ;  t<MAX_TENS_ORD ; t++ )
     {
       rid[t] = 1;
//DEBUGprint       ::printf("    rid[%d] = %d\n",t,rid[t]);
     }

   int *cd = new int[MAX_TENS_ORD];
   int *rd = new int[MAX_TENS_ORD];

   int *resd = new int[MAX_TENS_ORD];
   for( t=0 ;  t<MAX_TENS_ORD ; t++ )
     {
       resd[t] = 1;
//DEBUGprint       ::printf("    resd[%d] = %d\n",t,resd[t]);
     }

  double inerr=0.0;

// obrnuti red zbog brzeg izvrsenja ( rd[0] ide do najveceg broja . . . )
//  for ( rd[7]=1 ; rd[7]<=results_dims[7] ; rd[7]++ )
//   for ( rd[6]=1 ; rd[6]<=results_dims[6] ; rd[6]++ )
//    for ( rd[5]=1 ; rd[5]<=results_dims[5] ; rd[5]++ )
//     for ( rd[4]=1 ; rd[4]<=results_dims[4] ; rd[4]++ )
      for ( rd[3]=1 ; rd[3]<=results_dims[3] ; rd[3]++ )
       for ( rd[2]=1 ; rd[2]<=results_dims[2] ; rd[2]++ )
        for ( rd[1]=1 ; rd[1]<=results_dims[1] ; rd[1]++ )
         for ( rd[0]=1 ; rd[0]<=results_dims[0] ; rd[0]++ )
          {
            inerr = 0.0;

//            for ( cd[7]=1 ; cd[7]<=inerr_dims[7] ; cd[7]++ )
//             for ( cd[6]=1 ; cd[6]<=inerr_dims[6] ; cd[6]++ )
//              for ( cd[5]=1 ; cd[5]<=inerr_dims[5] ; cd[5]++ )
//               for ( cd[4]=1 ; cd[4]<=inerr_dims[4] ; cd[4]++ )
                for ( cd[3]=1 ; cd[3]<=inerr_dims[3] ; cd[3]++ )
                 for ( cd[2]=1 ; cd[2]<=inerr_dims[2] ; cd[2]++ )
                  for ( cd[1]=1 ; cd[1]<=inerr_dims[1] ; cd[1]++ )
                   for ( cd[0]=1 ; cd[0]<=inerr_dims[0] ; cd[0]++ )
                    {

// contracted indices
                      for ( t=0 ; t<contr_counter ; t++ )
                        {
                          lid[ this_contr[t]-1 ]  = cd[t];
//DEBUGprint                          ::printf("this contracted  t=%d  lid[%d]=%d     ",
//DEBUGprint                                    t,this_contr[t]-1,lid[this_contr[t]-1]);
                        }
//DEBUGprint                      ::printf("\n");
                      for ( t=0 ; t<contr_counter ; t++ )
                        {
                          rid[ arg_contr[t]-1 ]  = cd[t];
//DEBUGprint                          ::printf("arg contracted  t=%d  rid[%d]=%d    ",
//DEBUGprint                                    t,arg_contr[t]-1,rid[arg_contr[t]-1]);
                        }
//DEBUGprint                      ::printf("\n");


// uncontracted indices
                      for ( t=0 ; t<this_uni_count ; t++ )
                        {
                          lid[ this_uncontr[t]-1 ]  = rd[t];
//DEBUGprint                          ::printf("this uncontracted ----- t=%d lid[%d]=%d  ",
//DEBUGprint                                    t,this_uncontr[t]-1,lid[this_uncontr[t]-1]);
                        }
//DEBUGprint                      ::printf("\n");

//                      for ( t=0 ; t<arg_uni_count ; t++ )
                      for ( ; t<(this_uni_count+arg_uni_count) ; t++ )
                        {
                          rid[ arg_uncontr[t-this_uni_count]-1 ] = rd[t];
//DEBUGprint                          ::printf("arg uncontracted ----- t=%d  rid[%d]=%d  ",
//DEBUGprint                                    t,
//DEBUGprint                                    arg_uncontr[t-this_uni_count]-1,
//DEBUGprint                                    rid[arg_uncontr[t-this_uni_count]-1]);
                        }
//DEBUGprint                      ::printf("\n");



//                      inerr =
//                      this->val(lid[0],lid[1],lid[2],lid[3],
//                                       lid[4],lid[5],lid[6],lid[7]) *
//                      arg.val(rid[0],rid[1],rid[2],rid[3],
//                                     rid[4],rid[5],rid[6],rid[7]);

                      inerr = inerr +
                      this->val(lid[0],lid[1],lid[2],lid[3]) *
                      arg.val(rid[0],rid[1],rid[2],rid[3]);
//::printf(" inerr = %12.4lf \n   this[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %12.4lf arg[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %12.4lf\n",
//      inerr,lid[0],lid[1],lid[2],lid[3],lid[4],lid[5],lid[6],lid[7],
//      this->val(lid[0],lid[1],lid[2],lid[3],lid[4],lid[5],lid[6],lid[7]),
//      rid[0],rid[1],rid[2],rid[3],rid[4],rid[5],rid[6],rid[7],
//      arg.val(rid[0],rid[1],rid[2],rid[3],rid[4],rid[5],rid[6],rid[7]));

//DEBUGprint::printf("inerr=%6.2lf this[%1d,%1d,%1d,%1d]=%6.2lf arg[%1d,%1d,%1d,%1d]=%6.2lf\n",
//DEBUGprint      inerr,lid[0],lid[1],lid[2],lid[3],
//DEBUGprint      this->val(lid[0],lid[1],lid[2],lid[3]),
//DEBUGprint      rid[0],rid[1],rid[2],rid[3],
//DEBUGprint      arg.val(rid[0],rid[1],rid[2],rid[3]));

                    }

                                                           // might be optimized 
            for ( t=0 ; t<this_uni_count ; t++ )           // put first on in second
              {                                            // loop TS
//..                resd[ this_uncontr[t]-1 ] = rd[t];
                resd[t] = rd[t];
//DEBUGprint                ::printf("##### t=%d  resd[%d]=%d\n",
//DEBUGprint                          t,this_uncontr[t]-1,resd[t]);
//DEBUGprint                          t,t,resd[t]);
              }
            for ( ; t<(this_uni_count+arg_uni_count) ; t++ )
              {
//..                resd[ arg_uncontr[t]-1 ] = rd[t];
                resd[t] = rd[t];
//DEBUGprint                ::printf("##### t=%d  resd[%d]=%d\n",
//DEBUGprint                          t,arg_uncontr[t]-1,resd[t]);
//DEBUGprint                          t,t,resd[t]);
              }

//          result(resd[0],resd[1],resd[2],resd[3],
//                 resd[4],resd[5],resd[6],resd[7]) = inerr ;
//::printf(" -----------result[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %lf\n",



//      resd[0],resd[1],resd[2],resd[3],resd[4],resd[5],resd[6],resd[7],
//      result(resd[0],resd[1],resd[2],resd[3],resd[4],resd[5],resd[6],resd[7]));

          result.val(resd[0],resd[1],resd[2],resd[3]) = inerr ;
//DEBUGprint::printf("##################################### result[%1d,%1d,%1d,%1d]=%8.2lf\n",
//DEBUGprint             resd[0],resd[1],resd[2],resd[3],
//DEBUGprint      result(resd[0],resd[1],resd[2],resd[3]));


          }

// deleting dinamically allocated arrays
                             
    delete  [] this_contr;   
    delete  [] arg_contr;    
                             
    delete  [] this_uncontr; 
    delete  [] arg_uncontr ; 
                             
    delete  [] inerr_dims;   
                             
    delete  [] lid;          
    delete  [] rid;          
                             
    delete  [] cd;           
    delete  [] rd;           
                             
    delete  [] resd;         

// NULLification of indices

    null_indices();
    arg.null_indices();

//DEBUGprint::printf("     after delete in operator* coreleft is: %lu bytes\n", (unsigned long) ::coreleft());

    return result; // Returning a local variable ??
                   // copy-initializer happens before the destructor,
                   // so reference count is 2 when destructor is called,
                   // thus destructor doesn't free the memory.
  }



//......//##############################################################################
//......//##############################################################################
//......//##############################################################################
//......//##############################################################################
//......//-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//......//
//......//
//......//                                _____..---========+^+=========---.._____
//......//   ______________________ __,_='=====____  ================== _____=====`=
//......//  (._____________________I__) _ __=_/    `--------=+=--------'
//......//      /      /--...---===='---+----'
//......//     '------'---.--- _  _ =   _._'    "Make it so..."
//......//                    `--------'                  Captain Jean-Luc Picard
//......//                                                USS Enterprise, NCC-1701D
//......//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//......BJtensor operator*(BJtensor& lval, BJtensor& rval)
//......  {
//......//DEBUGprint::printf("before block in operator* coreleft is: %lu bytes\n", (unsigned long) ::coreleft());
//......//   const int MAX_TENS_ORD = 8;
//......   const int MAX_TENS_ORD = 4;
//......//....   int results_rank = 0;
//......
//......//DEBUGprint::printf("lval.rank()=%d\n",
//......//DEBUGprint          lval.rank());
//......//DEBUGprint::printf("rval.rank()=%d\n",
//......//DEBUGprint          rval.rank());
//......
//......// space for CONTRACTED indices
//......   int *lval_contr = new int[MAX_TENS_ORD];
//......   int *rval_contr  = new int[MAX_TENS_ORD];
//......
//......// space for UN-CONTRACTED indices
//......   int *lval_uncontr = new int[MAX_TENS_ORD];
//......   int *rval_uncontr = new int[MAX_TENS_ORD];
//......
//......   for(int lval_ic=0 ; lval_ic<MAX_TENS_ORD ; lval_ic++ )
//......     {
//......       lval_contr[lval_ic] = 0;
//......       lval_uncontr[lval_ic] = 0;
//......     }
//......   for( int rval_ic=0 ; rval_ic<MAX_TENS_ORD ; rval_ic++ )
//......     {
//......       rval_contr[rval_ic]  = 0;
//......       rval_uncontr[rval_ic]  = 0;
//......     }
//......
//......
//......
//......// now make the right choice for indices for this and rval BJtensor !!!
//......//
//......//                  |||~          |||~          |||~
//......//                 (0 0)         (0 0)         (0 0)
//......//--------------ooO-(_)-Ooo---ooO-(_)-Ooo---ooO-(_)-Ooo----------------------
//......//
//......
//......   char * lval_indices = lval.indices1;
//......   char * rval_indices = rval.indices1;
//......
//......// if the BJtensors are same then split indices :-)
//......   if ( lval.pc_nDarray_rep == rval.pc_nDarray_rep )
//......     {
//......       lval_indices = lval.indices2; // this is changed because
//......       rval_indices  = rval.indices1;  // it reads indices from
//......                                    // right so "in" indices1 are
//......                                   // the right ( rval ) and "in"
//......                                  // indices2 are the
//......                                 // left indices
//......     }
//......//DEBUGprint::printf("lval_indices=%s\n", lval_indices);
//......//DEBUGprint::printf("rval_indices=%s\n",  rval_indices);
//......
//......
//......
//......//   int lval_indices_number = ::strlen(lval.pc_nDarray_rep->indices1);
//......//   int  rval_indices_number = ::strlen(rval.pc_nDarray_rep->indices1);
//......
//......   int lval_indices_number = ::strlen(lval_indices);
//......   int rval_indices_number = ::strlen(rval_indices);
//......
//......//DEBUGprint::printf("lval_indices_number=%d\n",
//......//DEBUGprint          lval_indices_number);
//......//DEBUGprint::printf("rval_indices_number=%d\n",
//......//DEBUGprint          rval_indices_number);
//......
//......// check for indices number: should be EQUAL than rank of BJtensor
//......   if ( lval_indices_number != lval.rank() )
//......     {
//......       ::fprintf(stderr,"\a\n 'this' has more/less indices than expected:\n",
//......                " lval_indices_number = %d\n",
//......                " lval.rank() = %d\n",
//......                  lval_indices_number,
//......                  lval.rank());
//......       ::exit(1);
//......     }
//......
//......   if ( rval_indices_number != rval.rank() )
//......     {
//......       ::fprintf(stderr,"\a\n 'rval' has more/less indices than expected:\n",
//......                " rval_indices_number = %d\n",
//......                " rval.rank() = %d\n",
//......                  rval_indices_number,
//......                  rval.rank());
//......       ::exit(1);
//......     }
//......
//......// counter for contracted indices
//......   int contr_counter = lval.contracted_ind(lval_indices,
//......                                           rval_indices,
//......                                           lval_contr,
//......                                           rval_contr,
//......                                           lval_indices_number,
//......                                           rval_indices_number  );
//......
//......   int lval_uni_count = lval.uncontracted_ind( lval_uncontr,
//......                                               lval_contr,
//......                                               lval_indices_number);
//......
//......   int rval_uni_count = rval.uncontracted_ind( rval_uncontr,
//......                                               rval_contr,
//......                                               rval_indices_number);
//......
//......//number of UNcontractions
//......   int uncontr_counter = lval_uni_count + rval_uni_count;
//......
//......//TEMP  this is just the test because right now only up to order = 4
//......   if ( uncontr_counter > 4 )
//......     {
//......       ::fprintf(stderr,"\a\a\n\n  OOOPS product of multiplication has order %d\n",
//......                           uncontr_counter);
//......       ::exit(1);
//......     }
//......
//......
//......//DEBUGprint ::printf("..... lval_uni_count = %d, rval_uni_count = %d\n",
//......//DEBUGprint                 lval_uni_count,      rval_uni_count);
//......//DEBUGprint ::printf(" uncontr_counter %d\n",uncontr_counter);
//......//DEBUGprint for( int plval_ic=0 ;
//......//DEBUGprint      plval_ic<(lval.rank()+1) ;
//......//DEBUGprint      plval_ic++ )
//......//DEBUGprint   ::printf("lval_uncontr[%d]=%d\n",plval_ic,lval_uncontr[plval_ic]);
//......//DEBUGprint for( int prval_ic=0 ;
//......//DEBUGprint      prval_ic<(rval.rank()+1) ;
//......//DEBUGprint      prval_ic++ )
//......//DEBUGprint   ::printf("             rval_uncontr[%d]=%d\n",prval_ic,rval_uncontr[prval_ic]);
//......
//....../////////////////////////////////////////////////////
//......// let's make result BJtensor
//......
//......   static int results_dims[MAX_TENS_ORD];
//......   for( int t=0 ; t < lval_uni_count ; t++ )
//......     results_dims[t]=lval.dim()[lval_uncontr[t]-1];
//......   for( int a=0 ; a < rval_uni_count ; a++ )
//......     results_dims[a+t]=rval.dim()[rval_uncontr[a]-1];
//......// do kraja . . .
//......   for( ; a < MAX_TENS_ORD-t ; a++ )
//......     results_dims[a+t]=1;
//......
//......
//......
//......
//......
//......
//......// pa sada indeksi
//......   static char results_indices[MAX_TENS_ORD];
//......   for( t=0 ; t < lval_uni_count ; t++ )
//......     results_indices[t]=lval_indices[lval_uncontr[t]-1];
//......   for( a=0 ; a < rval_uni_count ; a++ )
//......     results_indices[a+t]=rval_indices[rval_uncontr[a]-1];
//......// the last one . . .
//......   results_indices[uncontr_counter] = '\0';
//......
//......
//......//DEBUGprint  ::printf("\n\n......  uncontr_counter+1=%d\n", uncontr_counter+1);
//......//DEBUGprintfor( int pa=0 ; pa < MAX_TENS_ORD ; pa++ )
//......//DEBUGprint  ::printf("   results_dims[%d]=%d\n",pa,results_dims[pa]);
//......
//......   BJtensor result(uncontr_counter, results_dims, 0.0);
//......   result(results_indices);  // initialize indices in result
//......//DEBUGprint  ::printf("\n\n..................   results_indices=%s\n",results_indices);
//......
//......//////////////////////////////////////////////////
//......// check dimensions in each tesnor order
//......   for( t=0 ; t<contr_counter ; t++ )
//......     if ( lval.dim()[lval_contr[t]-1] !=
//......          rval.dim()[rval_contr[t]-1] )
//......       {
//......         ::fprintf(stderr,"\a\a\n t = %d \n"
//......                  "lval_contr[t]-1  = %d\n"
//......                  "lval.dim()[%d] = %d\n"
//......                  "rval_contr[t]-1  = %d\n"
//......                  "rval.dim()[%d] = %d\n",
//......                   t,
//......                   lval_contr[t]-1,
//......                   lval_contr[t]-1,lval.dim()[lval_contr[t]-1],
//......                   rval_contr[t]-1,
//......                   rval_contr[t]-1,rval.dim()[rval_contr[t]-1]);
//......         ::exit(1);
//......       }
//......
//....../////////////////////////////////////////////////////#*%
//......//////////////////////////////////////////////////#*%
//......///////////////////////////////////////////////#*%
//......////////////////////////////////////////////#*%
//....../////////////////////////////////////////#*%
//......//////////////////////////////////////#*%
//......///////////////////////////////////#*%
//......////////////////////////////////#*%
//....../////////////////////////////#*%
//......//////////////////////////#*%
//......///////////////////////#*%
//......////////////////////#*%
//....../////////////////#*%
//......//////////////#*%
//......///////////#*%
//......////////#*%
//....../////#*%
//......//#*%
//......
//......   int *inerr_dims = new int[MAX_TENS_ORD];
//......   for( t=0 ; t<contr_counter ; t++ )
//......     {
//......       inerr_dims[t] = lval.dim()[lval_contr[t]-1];
//......//DEBUGprint       ::printf("    inerr_dims[%d] = %d\n",t,inerr_dims[t]);
//......     }
//......   for( ; t<MAX_TENS_ORD ; t++ )
//......     {
//......       inerr_dims[t] = 1;
//......//DEBUGprint       ::printf("          the rest of inerr_dims[%d] = %d\n",t,inerr_dims[t]);
//......     }
//......
//......
//......   int *lid = new int[MAX_TENS_ORD];
//......   for( t=0 ; t<MAX_TENS_ORD ; t++ )
//......     {
//......       lid[t] = 1;
//......//DEBUGprint       ::printf("    lid[%d] = %d\n",t,lid[t]);
//......     }
//......
//......   int *rid = new int[MAX_TENS_ORD];
//......   for( t=0 ;  t<MAX_TENS_ORD ; t++ )
//......     {
//......       rid[t] = 1;
//......//DEBUGprint       ::printf("    rid[%d] = %d\n",t,rid[t]);
//......     }
//......
//......   int *cd = new int[MAX_TENS_ORD];
//......   int *rd = new int[MAX_TENS_ORD];
//......
//......   int *resd = new int[MAX_TENS_ORD];
//......   for( t=0 ;  t<MAX_TENS_ORD ; t++ )
//......     {
//......       resd[t] = 1;
//......//DEBUGprint       ::printf("    resd[%d] = %d\n",t,resd[t]);
//......     }
//......
//......  double inerr=0.0;
//......
//......// obrnuti red zbog brzeg izvrsenja ( rd[0] ide do najveceg broja . . . )
//......//  for ( rd[7]=1 ; rd[7]<=results_dims[7] ; rd[7]++ )
//......//   for ( rd[6]=1 ; rd[6]<=results_dims[6] ; rd[6]++ )
//......//    for ( rd[5]=1 ; rd[5]<=results_dims[5] ; rd[5]++ )
//......//     for ( rd[4]=1 ; rd[4]<=results_dims[4] ; rd[4]++ )
//......      for ( rd[3]=1 ; rd[3]<=results_dims[3] ; rd[3]++ )
//......       for ( rd[2]=1 ; rd[2]<=results_dims[2] ; rd[2]++ )
//......        for ( rd[1]=1 ; rd[1]<=results_dims[1] ; rd[1]++ )
//......         for ( rd[0]=1 ; rd[0]<=results_dims[0] ; rd[0]++ )
//......          {
//......            inerr = 0.0;
//......
//......//            for ( cd[7]=1 ; cd[7]<=inerr_dims[7] ; cd[7]++ )
//......//             for ( cd[6]=1 ; cd[6]<=inerr_dims[6] ; cd[6]++ )
//......//              for ( cd[5]=1 ; cd[5]<=inerr_dims[5] ; cd[5]++ )
//......//               for ( cd[4]=1 ; cd[4]<=inerr_dims[4] ; cd[4]++ )
//......                for ( cd[3]=1 ; cd[3]<=inerr_dims[3] ; cd[3]++ )
//......                 for ( cd[2]=1 ; cd[2]<=inerr_dims[2] ; cd[2]++ )
//......                  for ( cd[1]=1 ; cd[1]<=inerr_dims[1] ; cd[1]++ )
//......                   for ( cd[0]=1 ; cd[0]<=inerr_dims[0] ; cd[0]++ )
//......
//......
//......
//......
//......
//......                    {
//......
//......// contracted indices
//......                      for ( t=0 ; t<contr_counter ; t++ )
//......                        {
//......                          lid[ lval_contr[t]-1 ]  = cd[t];
//......//DEBUGprint                          ::printf("this contracted  t=%d  lid[%d]=%d     ",
//......//DEBUGprint                                    t,lval_contr[t]-1,lid[lval_contr[t]-1]);
//......                        }
//......//DEBUGprint                      ::printf("\n");
//......                      for ( t=0 ; t<contr_counter ; t++ )
//......                        {
//......                          rid[ rval_contr[t]-1 ]  = cd[t];
//......//DEBUGprint                          ::printf("rval contracted  t=%d  rid[%d]=%d    ",
//......//DEBUGprint                                    t,rval_contr[t]-1,rid[rval_contr[t]-1]);
//......                        }
//......//DEBUGprint                      ::printf("\n");
//......
//......
//......// uncontracted indices
//......                      for ( t=0 ; t<lval_uni_count ; t++ )
//......                        {
//......                          lid[ lval_uncontr[t]-1 ]  = rd[t];
//......//DEBUGprint                          ::printf("this uncontracted ----- t=%d lid[%d]=%d  ",
//......//DEBUGprint                                    t,lval_uncontr[t]-1,lid[lval_uncontr[t]-1]);
//......                        }
//......//DEBUGprint                      ::printf("\n");
//......
//......//                      for ( t=0 ; t<rval_uni_count ; t++ )
//......                      for ( ; t<(lval_uni_count+rval_uni_count) ; t++ )
//......                        {
//......                          rid[ rval_uncontr[t-lval_uni_count]-1 ] = rd[t];
//......//DEBUGprint                          ::printf("rval uncontracted ----- t=%d  rid[%d]=%d  ",
//......//DEBUGprint                                    t,
//......//DEBUGprint                                    rval_uncontr[t-lval_uni_count]-1,
//......//DEBUGprint                                    rid[rval_uncontr[t-lval_uni_count]-1]);
//......                        }
//......//DEBUGprint                      ::printf("\n");
//......
//......
//......
//......//                      inerr =
//......//                      lval.val(lid[0],lid[1],lid[2],lid[3],
//......//                                       lid[4],lid[5],lid[6],lid[7]) *
//......//                      rval.val(rid[0],rid[1],rid[2],rid[3],
//......//                                     rid[4],rid[5],rid[6],rid[7]);
//......
//......                      inerr = inerr +
//......                      lval.val(lid[0],lid[1],lid[2],lid[3]) *
//......                      rval.val(rid[0],rid[1],rid[2],rid[3]);
//......//::printf(" inerr = %12.4lf \n   this[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %12.4lf rval[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %12.4lf\n",
//......//      inerr,lid[0],lid[1],lid[2],lid[3],lid[4],lid[5],lid[6],lid[7],
//......//      lval.val(lid[0],lid[1],lid[2],lid[3],lid[4],lid[5],lid[6],lid[7]),
//......//      rid[0],rid[1],rid[2],rid[3],rid[4],rid[5],rid[6],rid[7],
//......//      rval.val(rid[0],rid[1],rid[2],rid[3],rid[4],rid[5],rid[6],rid[7]));
//......
//......//DEBUGprint::printf("inerr=%6.2lf this[%1d,%1d,%1d,%1d]=%6.2lf rval[%1d,%1d,%1d,%1d]=%6.2lf\n",
//......//DEBUGprint      inerr,lid[0],lid[1],lid[2],lid[3],
//......//DEBUGprint      lval.val(lid[0],lid[1],lid[2],lid[3]),
//......//DEBUGprint      rid[0],rid[1],rid[2],rid[3],
//......//DEBUGprint      rval.val(rid[0],rid[1],rid[2],rid[3]));
//......
//......                    }
//......
//......
//......            for ( t=0 ; t<lval_uni_count ; t++ )
//......              {
//......//..                resd[ lval_uncontr[t]-1 ] = rd[t];
//......                resd[t] = rd[t];
//......//DEBUGprint                ::printf("##### t=%d  resd[%d]=%d\n",
//......//DEBUGprint                          t,lval_uncontr[t]-1,resd[t]);
//......//DEBUGprint                          t,t,resd[t]);
//......              }
//......            for ( ; t<(lval_uni_count+rval_uni_count) ; t++ )
//......              {
//......//..                resd[ rval_uncontr[t]-1 ] = rd[t];
//......                resd[t] = rd[t];
//......//DEBUGprint                ::printf("##### t=%d  resd[%d]=%d\n",
//......//DEBUGprint                          t,rval_uncontr[t]-1,resd[t]);
//......//DEBUGprint                          t,t,resd[t]);
//......              }
//......
//......//          result(resd[0],resd[1],resd[2],resd[3],
//......//                 resd[4],resd[5],resd[6],resd[7]) = inerr ;
//......//::printf(" -----------result[%1d,%1d,%1d,%1d,%1d,%1d,%1d,%1d] = %lf\n",
//......
//......
//......
//......//      resd[0],resd[1],resd[2],resd[3],resd[4],resd[5],resd[6],resd[7],
//......//      result(resd[0],resd[1],resd[2],resd[3],resd[4],resd[5],resd[6],resd[7]));
//......
//......          result.val(resd[0],resd[1],resd[2],resd[3]) = inerr ;
//......//DEBUGprint::printf("##################################### result[%1d,%1d,%1d,%1d]=%8.2lf\n",
//......//DEBUGprint             resd[0],resd[1],resd[2],resd[3],
//......//DEBUGprint      result(resd[0],resd[1],resd[2],resd[3]));
//......
//......
//......          }
//......
//......// deleting dinamically allocated arrays
//......
//......    delete [] lval_contr;
//......    delete [] rval_contr;
//......
//......    delete [] lval_uncontr;
//......    delete [] rval_uncontr ;
//......
//......    delete [] inerr_dims;
//......
//......    delete [] lid;
//......    delete [] rid;
//......
//......    delete [] cd;
//......    delete [] rd;
//......
//......    delete [] resd;
//......
//......// NULLification of indices
//......
//......    lval.null_indices();
//......    rval.null_indices();
//......
//......//DEBUGprint::printf("     after delete in operator* coreleft is: %lu bytes\n", (unsigned long) ::coreleft());
//......
//......    return BJtensor(result); // Returning a local variable ??
//......                   // copy-initializer happens before the destructor,
//......                   // so reference count is 2 when destructor is called,
//......                   // thus destructor doesn't free the memory.
//......  }
//......



//##############################################################################
// counter for contracted indices
int BJtensor::contracted_ind(char  * argleft_indices,
                           char  * argright_indices,
                           int   * argleft_contr,
                           int   * argright_contr,
                           int     argleft_ind_numb,
                           int     argright_ind_numb)
  {
    int contr_counter = 0;
    int argleft_i_count = 0;
    int argright_i_count  = 0;

    for ( int argleft_indices_counter = 1 ;
          argleft_indices_counter<=argleft_ind_numb ;
          argleft_indices_counter++ )
      {
        int t_i_c = argleft_indices_counter - 1;
        for ( int argright_indices_counter = 1 ;
              argright_indices_counter<= argright_ind_numb ;
              argright_indices_counter++ )
          {
            int a_i_c = argright_indices_counter - 1;
            if ( argleft_indices[t_i_c] ==
                 argright_indices[a_i_c] )
              {
                argleft_contr[argleft_i_count] = argleft_indices_counter;
                argright_contr[argright_i_count] = argright_indices_counter;

                argleft_i_count++;
                argright_i_count++;

                contr_counter++;
              }
          }
      }


//DEBUGprint ::printf("..... argleft_i_count = %d, argright_i_count = %d\n",
//DEBUGprint                 argleft_i_count,      argright_i_count);
//DEBUGprint ::printf(" contr_counter %d\n",contr_counter);
//DEBUGprint for(int pthis_ic=0 ;
//DEBUGprint     pthis_ic<(argleft->rank()+1) ;
//DEBUGprint     pthis_ic++ )
//DEBUGprint   ::printf("argleft_contr[%d]=%d\n",pthis_ic,argleft_contr[pthis_ic]);
//DEBUGprint for( int parg_ic=0 ;
//DEBUGprint      parg_ic<(argright->rank()+1) ;
//DEBUGprint      parg_ic++ )
//DEBUGprint   ::printf("                 argright_contr[%d]=%d\n",
//DEBUGprint                              parg_ic,argright_contr[parg_ic]);


    return contr_counter;
  }

//##############################################################################
// counter for UNcontracted indices
int BJtensor::uncontracted_ind(int *tens_uncontr,
                             int *tens_contr,
                             int tens_ind_numb)
  {
// counter for UNcontracted indices
    int tens_uni_count = 0;

    for ( int tens_unindices_counter = 1 ;
          tens_unindices_counter<=tens_ind_numb ;
          tens_unindices_counter++ )
      {
        tens_uncontr[tens_uni_count] = tens_unindices_counter;
        for ( int tens_indices_counter = 1 ;
              tens_indices_counter<=tens_ind_numb ;
              tens_indices_counter++ )
          {
            int t_i_c = tens_indices_counter - 1;
//DEBUGprint::printf("CHECK tens_uncontr[%d]=%d",tens_uni_count, tens_uncontr[tens_uni_count]);
//DEBUGprint::printf("  WITH tens_contr[%d]=%d\n", t_i_c,tens_contr[t_i_c]);
            if ( tens_uncontr[tens_uni_count] == tens_contr[t_i_c] )
              {
                tens_uncontr[tens_uni_count] = 0;
                tens_uni_count--;
                tens_indices_counter=tens_ind_numb;//out of inner for loop
              }
          }
        tens_uni_count++;
      }
    return tens_uni_count;
  }

//##############################################################################
// BJtensorial division THE rval MUST BE 0 ORDER BJtensor
BJtensor BJtensor::operator/( BJtensor & rval)
 {
// construct BJtensor using the same control numbers as for the
// original one.
    BJtensor div(this->rank(), dim(), 0.0);
    double rvalDouble = rval.trace();
    if ( rval.rank() != 1 ) ::printf("rval.rank() != 1 for BJtensor BJtensor::operator/( BJtensor & rval)\n");

    div.indices1 = this->indices1;

    switch(this->rank())
      {
        case 0:
          {
            div.val(1) = val(1)/rvalDouble;
            break;
          }

        case 1:
          {
            for ( int i1=1 ; i1<=3 ; i1++ )
              div.val(i1) = val(i1)/rvalDouble;
            break;
          }

        case 2:
          {
            for ( int i2=1 ; i2<=3 ; i2++ )
              for ( int j2=1 ; j2<=3 ; j2++ )
                div.val(i2, j2) = val(i2, j2)/rvalDouble;
            break;
          }

        case 3:
          {
            for ( int i3=1 ; i3<=3 ; i3++ )
              for ( int j3=1 ; j3<=3 ; j3++ )
                for ( int k3=1 ; k3<=3 ; k3++ )
                  div.val(i3, j3, k3) = val(i3, j3, k3)/rvalDouble;
            break;
          }

        case 4:
          {
            for ( int i4=1 ; i4<=3 ; i4++ )
              for ( int j4=1 ; j4<=3 ; j4++ )
                for ( int k4=1 ; k4<=3 ; k4++ )
                  for ( int l4=1 ; l4<=3 ; l4++ )
                    div.val(i4,j4,k4,l4)=val(i4,j4,k4,l4)/rvalDouble;
            break;
          }
      }

    null_indices();

    return div;
 }



//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> ikjl    */
BJtensor BJtensor::transpose0110() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, k, j, l);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> ikjl    */
BJtensor BJtensor::transposeoverbar() const // same as transpose0110
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, k, j, l);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }
//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> ilkj    */
BJtensor BJtensor::transpose0101() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, l, k, j);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> iljk    */
BJtensor BJtensor::transpose0111() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, l, j, k);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }
//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> iljk    */
BJtensor BJtensor::transposeunderbar() const   // transpose0111() 
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, l, j, k);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> jikl
BJtensor BJtensor::transpose1100() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(j, i, k, l);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> ijlk
BJtensor BJtensor::transpose0011() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(i, j, l, k);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* transpose function for 4th rank BJtensors:
//     ijkl  -->> iilk
BJtensor BJtensor::transpose1001() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       for ( int k=1 ; k<=dim()[2] ; k++ )
         for ( int l=1 ; l<=dim()[3] ; l++ )
           trans1.val(i, j, k, l) = cval(l, j, k, i);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }




//##############################################################################
///* transpose function for 2th rank BJtensors:
//     ij  -->> ji    */
BJtensor BJtensor::transpose11() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than transpose it.
   BJtensor trans1(this->rank(), dim(), 0.0);

   for ( int i=1 ; i<=dim()[0] ; i++ )
     for ( int j=1 ; j<=dim()[1] ; j++ )
       trans1.val(i, j) = cval(j, i);

    trans1.null_indices();
//    null_indices();
    return trans1;
 }

//##############################################################################
///* symmterize function for 2th rank BJtensors:
//     ij  
BJtensor BJtensor::symmetrize11() const
  {
// construct BJtensor using the same control numbers as for the
// original one and than symmetrize it.
//   BJtensor sym(this->rank(), dim(), 0.0);


   BJtensor temp(*this);

   BJtensor sym = (temp + temp.transpose11())*0.5;

   sym.null_indices();
   return sym;
 }

//..//##############################################################################
//..///* symmterize function for 2th rank BJtensors:
//..//     ij  
//..BJtensor BJtensor::symmetrize11() const
//..  {
//..// construct BJtensor using the same control numbers as for the
//..// original one and than symmetrize it.
//..//   BJtensor sym(this->rank(), dim(), 0.0);
//..
//..
//..   BJtensor temp(*this);
//..
//..   BJtensor sym = (temp + temp.transpose11())*0.5;
//..
//..   sym.null_indices();
//..   return sym;
//.. }
//..


// moved to ndarray april 3 1995
//--//##############################################################################
//--// trace function: trace of a second rank BJtensor
//--// what about fourth ( 4th ) rank BJtensor trace or any other rank ??
//--double BJtensor::trace() const 
//--  {
//--    double tr = 0.0;
//--// ovaj case ne treba vec moze sve do 4-tog reda ( ili kasnije osmog # )
//--    switch(this->rank())
//--      {
//--        case 0:
//--          {
//--//            ::printf(" trace=%.4e  ", val(1));
//--//            ::printf("\n");
//--            tr = cval(1);
//--            break;
//--          }
//--
//--        case 1:
//--          {
//--            if(dim()[0] != 1)
//--              {
//--::printf("\a\nERROR in trace function : not a squared 1-st rank BJtensor\n");
//--::exit( 1 );
//--              }
//--            tr = cval(1);
//--            break;
//--          }
//--
//--        case 2:
//--          {
//--            if(dim()[0] != dim()[1])
//--              {
//--::printf("\a\nERROR in trace function : not a sqared 2nd-rank BJtensor\n");
//--::exit( 1 );
//--              }
//--            for ( int i2=1 ; i2<=dim()[0] ; i2++ )
//--              tr += cval(i2, i2);
//--            break;
//--          }
//--
//--        case 3:
//--          {
//--            if( dim()[0] != dim()[1] ||
//--                dim()[1] != dim()[2] ||
//--                dim()[2] != dim()[0]    )
//--              {
//--::printf("\a\nERROR in trace function : not a sqared 3nd-rank BJtensor\n");
//--::exit( 1 );
//--              }
//--            for ( int i3=1 ; i3<=dim()[0] ; i3++ )
//--              tr += cval(i3, i3, i3);
//--            break;
//--          }
//--
//--        case 4:
//--          {
//--            if( dim()[0] != dim()[1] ||
//--                dim()[1] != dim()[2] ||
//--                dim()[2] != dim()[3] ||
//--                dim()[3] != dim()[0]    )
//--              {
//--::printf("\a\nERROR in trace function : not a sqared 4nd-rank BJtensor\n");
//--::exit( 1 );
//--              }
//--            for ( int i4=1 ; i4<=dim()[0] ; i4++ )
//--              tr += cval(i4, i4, i4, i4);
//--            break;
//--          }
//--      }
//--
//..........    null_indices();
//--
//--    return tr;
//--  }


//##############################################################################
// determinant function
double BJtensor::determinant() const
  {
    double det = 0.0;
//    switch(pc_nDarray_rep->nDarray_rank)
    switch(this->rank())
      {
        case 0:
          {
            det = cval(1);
            break;
          }

        case 1:
          {
            if(dim()[0] != 1)
              {
                ::fprintf(stderr,"\a\n ERROR: only for square BJvector (1) \n");
                ::exit(1);
              }
            det = cval(1);
            break;
          }

        case 2:
          {
            BJmatrix converted = this->BJtensor2BJmatrix_2();
            det = converted.determinant();
            break;
          }

        case 3:
          {
            ::fprintf(stderr,"\a\n SORRY: not implemented \n");
            ::exit(1);
            break;
          }

        case 4:
          {
            ::fprintf(stderr,"\a\n SORRY: not implemented \n");
            ::exit(1);
            break;
          }
      }

//..........    null_indices();

    return det;
  }


//##############################################################################
BJmatrix BJtensor::BJtensor2BJmatrix_1() const // convert BJtensor of even order to BJmatrix
                                // to be used in inversion process
                                // in the begining only BJtensor of
                                // 2 and 4 order
                                // I_ijkl scheme
  {
    int BJtensor_order = this->rank();
    static int BJmatrix_dims[] = {0,0};

    if ( BJtensor_order == 2 )
      {
        BJmatrix_dims[0] = this->dim(1);
        BJmatrix_dims[1] = this->dim(2);
      }
    else if ( BJtensor_order == 4 )
      {
        BJmatrix_dims[0] = (this->dim(1))*(this->dim(3));
        BJmatrix_dims[1] = (this->dim(2))*(this->dim(4));
      }
    else
      {
        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
        ::exit(1);
      }

    BJmatrix converted(BJmatrix_dims[0], BJmatrix_dims[1], 0.0);


   int m41 = 0;
   int m42 = 0;
// filling out the converted BJmatrix
    if ( BJtensor_order == 2 )
      {
        for ( int c22=1 ; c22<=this->dim(2) ; c22++ )
          for ( int c21=1 ; c21<=this->dim(1) ; c21++ )
            converted.val(c21,c22)=this->cval(c21,c22);
      }
    else if ( BJtensor_order == 4 )
      {
        for ( int c44=1 ; c44<=this->dim(4) ; c44++ )
          for ( int c43=1 ; c43<=this->dim(3) ; c43++ )
            for ( int c42=1 ; c42<=this->dim(2) ; c42++ )
              for ( int c41=1 ; c41<=this->dim(1) ; c41++ )
                {
                  m41 = this->dim(1)*(c41-1)+c43;
                  m42 = this->dim(2)*(c42-1)+c44;

                  converted.val( m41, m42 ) = this->cval( c41, c42, c43, c44 );

                }
      }
    else
      {
        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
        ::exit(1);
      }
//    converted.print("t","\n\n converted \n");
    return converted;
  }


//##############################################################################
BJmatrix BJtensor::BJtensor2BJmatrix_2() const // convert BJtensor of even order to BJmatrix
                                // to be used in inversion process
                                // in the begining only BJtensor of
                                // 2 and 4 order
                                // I_ikjl scheme
  {
    int BJtensor_order = this->rank();
    static int BJmatrix_dims[] = {0,0};

    if ( BJtensor_order == 2 )
      {
        BJmatrix_dims[0] = this->dim(1);
        BJmatrix_dims[1] = this->dim(2);
      }
    else if ( BJtensor_order == 4 )
      {
        BJmatrix_dims[0] = (this->dim(1))*(this->dim(2));
        BJmatrix_dims[1] = (this->dim(3))*(this->dim(4));
      }
    else
      {
        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
        ::exit(1);
      }

    BJmatrix converted(BJmatrix_dims[0], BJmatrix_dims[1], 0.0);


   int m41 = 0;
   int m42 = 0;
// filling out the converted BJmatrix
    if ( BJtensor_order == 2 )
      {
        for ( int c22=1 ; c22<=this->dim(2) ; c22++ )
          for ( int c21=1 ; c21<=this->dim(1) ; c21++ )
            converted.val(c21,c22)=this->cval(c21,c22);
      }
    else if ( BJtensor_order == 4 )
      {
        for ( int c44=1 ; c44<=this->dim(4) ; c44++ )
          for ( int c43=1 ; c43<=this->dim(3) ; c43++ )
            for ( int c42=1 ; c42<=this->dim(2) ; c42++ )
              for ( int c41=1 ; c41<=this->dim(1) ; c41++ )
                {
                  m41 = this->dim(1)*(c41-1)+c42;
                  m42 = this->dim(3)*(c43-1)+c44;

                  converted.val( m41, m42 ) = this->cval( c41, c42, c43, c44 );

                }
      }
    else
      {
        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
        ::exit(1);
      }

//    converted.print("t","\n\n converted \n");

    return converted;
  }


//~~~~//##############################################################################
//~~~~BJmatrix BJtensor::BJtensor2BJmatrix_3()  // convert BJtensor of even order to BJmatrix
//~~~~                                // to be used in inversion process
//~~~~                                // in the begining only BJtensor of
//~~~~                                // 2 and 4 order
//~~~~                                // I_iljk scheme
//~~~~  {
//~~~~    int BJtensor_order = this->rank();
//~~~~    static int BJmatrix_dims[] = {0,0};
//~~~~
//~~~~    if ( BJtensor_order == 2 )
//~~~~      {
//~~~~        BJmatrix_dims[0] = this->dim(1);
//~~~~        BJmatrix_dims[1] = this->dim(2);
//~~~~      }
//~~~~    else if ( BJtensor_order == 4 )
//~~~~      {
//~~~~        BJmatrix_dims[0] = (this->dim(1))*(this->dim(4));
//~~~~        BJmatrix_dims[1] = (this->dim(2))*(this->dim(3));
//~~~~      }
//~~~~    else
//~~~~      {
//~~~~        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
//~~~~        ::exit(1);
//~~~~      }
//~~~~
//~~~~    BJmatrix converted(BJmatrix_dims[0], BJmatrix_dims[1], 0.0);
//~~~~
//~~~~
//~~~~   int m41 = 0;
//~~~~   int m42 = 0;
//~~~~// filling out the converted BJmatrix
//~~~~    if ( BJtensor_order == 2 )
//~~~~      {
//~~~~        for ( int c22=1 ; c22<=this->dim(2) ; c22++ )
//~~~~          for ( int c21=1 ; c21<=this->dim(1) ; c21++ )
//~~~~            converted.val(c21,c22)=this->val(c21,c22);
//~~~~      }
//~~~~    else if ( BJtensor_order == 4 )
//~~~~      {
//~~~~        for ( int c44=1 ; c44<=this->dim(4) ; c44++ )
//~~~~          for ( int c43=1 ; c43<=this->dim(3) ; c43++ )
//~~~~            for ( int c42=1 ; c42<=this->dim(2) ; c42++ )
//~~~~              for ( int c41=1 ; c41<=this->dim(1) ; c41++ )
//~~~~                {
//~~~~                  m41 = this->dim(1)*(c41-1)+c44;
//~~~~                  m42 = this->dim(2)*(c42-1)+c43;
//~~~~
//~~~~                  converted.val( m41, m42 ) = this->val( c41, c42, c43, c44 );
//~~~~
//~~~~                }
//~~~~      }
//~~~~    else
//~~~~      {
//~~~~        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
//~~~~        ::exit(1);
//~~~~      }
//~~~~
//~~~~//    converted.print("t","\n\n converted \n");
//~~~~
//~~~~    return converted;
//~~~~  }
//~~~~
//~~~~
//~~~~

//..//##############################################################################
//..BJtensor BJtensor::inverse()  // invert BJtensor of even rank by
//..                          // converting it to BJmatrix
//..                          // and than inverting the BJmatrix
//..                          // and at the end puting it back
//..                          // again ( back conversion )
//..  {
//..    static int BJmatrix_dims[] = {0,0};
//..    if ( this->rank() == 2 )
//..      {
//..        BJmatrix_dims[0] = this->dim(1);
//..        BJmatrix_dims[1] = this->dim(2);
//..      }
//..    else if ( this->rank() == 4 )
//..      {
//..        BJmatrix_dims[0] = (this->dim(1))*(this->dim(2));
//..        BJmatrix_dims[1] = (this->dim(3))*(this->dim(4));
//..      }
//..    else
//..      {
//..        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
//..        ::exit(1);
//..      }
//..    BJmatrix converted(BJmatrix_dims[0], BJmatrix_dims[1], 0.0);
//..
//..   int m41 = 0;
//..   int m42 = 0;
//..
//..// filling out the converted BJmatrix
//..    if ( this->rank() == 2 )
//..      {
//..        for ( int c22=1 ; c22<=this->dim(2) ; c22++ )
//..          for ( int c21=1 ; c21<=this->dim(1) ; c21++ )
//..            converted.val(c21,c22)=this->val(c21,c22);
//..      }
//..    else if ( this->rank() == 4 )
//..      {
//..        for ( int c44=1 ; c44<=this->dim(4) ; c44++ )
//..          for ( int c43=1 ; c43<=this->dim(3) ; c43++ )
//..            for ( int c42=1 ; c42<=this->dim(2) ; c42++ )
//..              for ( int c41=1 ; c41<=this->dim(1) ; c41++ )
//..                {
//..                  m41 = this->dim(1)*(c41-1)+c42;
//..                  m42 = this->dim(3)*(c43-1)+c44;
//..
//..                  converted.val( m41, m42 ) = this->val( c41, c42, c43, c44 );
//..                }
//..      }
//..    else
//..      {
//..        ::fprintf(stderr,"\a\n only BJtensor rank 2 and 4 #\n");
//..        ::exit(1);
//..      }
//..
//..//  converted.print("c","\n Converted:\n");
//..  BJmatrix converted_inverse = converted.inverse();
//..
//..// building up a BJtensor result
//..  BJtensor result(this->rank(), this->dim(), 0.0);
//..
//..// filling back the inverted BJmatrix to BJtensor
//..    if ( this->rank() == 2 )
//..      {
//..        for ( int c22=1 ; c22<=this->dim(2) ; c22++ )
//..          for ( int c21=1 ; c21<=this->dim(1) ; c21++ )
//..            result.val(c21,c22)=converted_inverse.val(c21,c22);
//..      }
//..    else if ( this->rank() == 4 )
//..      {
//..        for ( int c44=1 ; c44<=this->dim(4) ; c44++ )
//..          for ( int c43=1 ; c43<=this->dim(3) ; c43++ )
//..            for ( int c42=1 ; c42<=this->dim(2) ; c42++ )
//..              for ( int c41=1 ; c41<=this->dim(1) ; c41++ )
//..                {
//..                  m41 = this->dim(1)*(c41-1)+c42;
//..                  m42 = this->dim(3)*(c43-1)+c44;
//..
//..                  result.val(c41,c42,c43,c44) = converted_inverse.val(m41,m42);
//..                }
//..      }
//..//    result.print("t","back to BJtensor result");
//..    return result;
//..  }
//..

//~~~~//##############################################################################
//~~~~BJtensor BJtensor::inverse_1()  // invert BJtensor of even rank by
//~~~~                            // converting it to BJmatrix
//~~~~                            // and than inverting the BJmatrix
//~~~~                            // and at the end puting it back
//~~~~                            // again ( back conversion )
//~~~~                            // I_ijkl scheme
//~~~~
//~~~~  {
//~~~~
//~~~~  BJmatrix converted = this->BJtensor2BJmatrix_1();
//~~~~
//~~~~//  converted.print("c","\n Converted:\n");
//~~~~
//~~~~  BJmatrix converted_inverse = converted.inverse();
//~~~~
//~~~~//  converted_inverse.print("t","\n\n converted \n");
//~~~~
//~~~~  BJtensor result = converted_inverse.BJmatrix2BJtensor_1();
//~~~~
//~~~~//    result.print("t","back to BJtensor result");
//~~~~
//~~~~    return result;
//~~~~  }


//##############################################################################
BJtensor BJtensor::inverse_2()  const // invert BJtensor of even rank by
                            // converting it to BJmatrix
                            // and than inverting the BJmatrix
                            // and at the end puting it back
                            // again ( back conversion )
                            // I_ikjl scheme

{

  BJmatrix converted = this->BJtensor2BJmatrix_2();

//  converted.print("c","\n Converted:\n");

  BJmatrix converted_inverse = converted.inverse();

//  converted_inverse.print("t","\n\n converted_inverse \n");

  BJtensor result = converted_inverse.BJmatrix2BJtensor_2();

//  result.print("t","back to BJtensor result");

  return result;
}


//##############################################################################
BJtensor BJtensor::inverse()  const // invert BJtensor of even rank by
                          // converting it to BJmatrix
                          // and than inverting the BJmatrix
                          // and at the end puting it back
                          // again ( back conversion )
                          // I_ikjl scheme   TRUE ONE

  {

  int BJmatrix_or_BJtensor4 = this->rank();

  BJmatrix converted = this->BJtensor2BJmatrix_2();

  BJtensor result;

//converted.print("c","\n Converted:\n");

  BJmatrix converted_inverse = converted.inverse();

//converted_inverse.print("t","\n\n Converted_inverse \n");

  if ( BJmatrix_or_BJtensor4 == 4 ) 
    {
      result = converted_inverse.BJmatrix2BJtensor_2();
    }
  else if ( BJmatrix_or_BJtensor4 == 2 ) 
    {
      result = converted_inverse.BJmatrix2BJtensor_22();
//result.print("r","\n\n result  = converted_inverse.BJmatrix2BJtensor_2();\n");
    }
  else 
    {
      ::fprintf(stderr,"only 2 and 4 BJtensors \n");
      ::exit(1);
    }

//    result.print("t","back to BJtensor result");
#ifdef SASA
//***** Dodo Sasa
//  result.pc_nDarray_rep->dim=this->dim();
//************** Dosta prljavo ( ako se ubije treba videti da se
// iskopira niz a ne pointer jer ovako moze da bude problema ako 
// se unisti jedan od nizova !!!
// Takodje ako su dimenzije tenzora razlicite bi ce problema
 if(BJmatrix_or_BJtensor4==4){
  result.pc_nDarray_rep->dim[0]=this->dim()[2];
  result.pc_nDarray_rep->dim[1]=this->dim()[3];
  result.pc_nDarray_rep->dim[2]=this->dim()[0];
  result.pc_nDarray_rep->dim[3]=this->dim()[1];
  }
 else {
  result.pc_nDarray_rep->dim[0]=this->dim()[1];
  result.pc_nDarray_rep->dim[1]=this->dim()[0];
  }
#endif
// ***************** Kraj izmena koje je dodao Sasa
    return result;
  }


//~~~~//##############################################################################
//~~~~BJtensor BJtensor::inverse_3()  // invert BJtensor of even rank by
//~~~~                            // converting it to BJmatrix
//~~~~                            // and than inverting the BJmatrix
//~~~~                            // and at the end puting it back
//~~~~                            // again ( back conversion )
//~~~~                            // I_ilkj scheme
//~~~~
//~~~~  {
//~~~~
//~~~~  BJmatrix converted = this->BJtensor2BJmatrix_3();
//~~~~
//~~~~//  converted.print("c","\n Converted:\n");
//~~~~
//~~~~  BJmatrix converted_inverse = converted.inverse();
//~~~~
//~~~~//  converted_inverse.print("t","\n\n converted \n");
//~~~~
//~~~~  BJtensor result = converted_inverse.BJmatrix2BJtensor_3();
//~~~~
//~~~~//    result.print("t","back to BJtensor result");
//~~~~
//~~~~    return result;
//~~~~  }
//~~~~


//.~~~~      //##############################################################################
//.~~~~      // use "from" and initialize already allocated BJtensor from "from" values
//.~~~~      void BJtensor::Initialize( const BJtensor & from )
//.~~~~        {
//.~~~~      // copy only data because everything else has already been defined
//.~~~~      // WATCH OUT IT HAS TO BE DEFINED BEFORE THIS FUNCTIONS IS CALLED
//.~~~~          for ( int i=0 ; i<pc_nDarray_rep->total_numb ; i++ )
//.~~~~            this->pc_nDarray_rep->pd_nDdata[i] = from.pc_nDarray_rep->pd_nDdata[i] ;
//.~~~~        }


//##############################################################################
char * BJtensor::f_indices1( void ) const
  {
    return( this->indices1 );
  }
//##############################################################################
char * BJtensor::f_indices2( void ) const
  {
    return( this->indices2 );
  }


#endif

