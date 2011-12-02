                                                                        
// $Revision: 1.2 $                                                              
// $Date: 2003-02-14 23:01:48 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/nDarray/BJmatrix.cpp,v $                                                                

                                                                        
                                                                        
///*
//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:             BJmatrix                                                  #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                    Bruce Eckel { changed/improved a lot by BJ }              #
//#                                                                              #
//# DATE:              November '92                                              #
//# UPDATE HISTORY:    05 - __ avgust '93.  redefined as derived class from      #
//#                                 nDarray class                                #
//#                    january 06 '93  added BJmatrix2BJtensor_1, BJmatrix2BJtensor_2    #
//#                                   BJmatrix2BJtensor_3                        #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                    May 10 2000 Xiaoyan found a bug with 			 #
//#        BJmatrix::BJmatrix(int rows, int columns, double *initvalues):        #
//#           nDarray( 2, rows, columns, initvalues){ } 			 #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
#ifndef MATRIX_CC
#define MATRIX_CC
//  #include "basics.h"
//  #include "nDarray.h"
#include "BJmatrix.h"
//#include "BJvector.h"
#include <OPS_Stream.h>

#define TINY  1e-20
void BJmatrix::error(char * msg1, char * msg2)
  {
    ::fprintf(stderr,"BJmatrix error: %s %s\n", msg1, msg2);
    exit( 1 );
  }

//##############################################################################
BJmatrix::BJmatrix(int rows, int columns, double initval):
  nDarray( 2, rows, columns, initval){ } // calling the appropriate
                                          // base constructor
//##############################################################################
BJmatrix::BJmatrix(int rows, int columns, double *initvalues):
  nDarray( 2, rows, columns, initvalues){ } // calling the appropriate
                                             // base constructor
//Xiaoyan found a bug it should be initvalue instead of *initvalue
//##############################################################################
//  special for vector
BJmatrix::BJmatrix(int rank, int rows, int columns, double *initvalues):
  nDarray( rank, rows, columns, initvalues){ } // calling the appropriate
                                             // base constructor
//##############################################################################
//  special for vector
BJmatrix::BJmatrix(int rank, int rows, int columns, double initvalues):
  nDarray( rank, rows, columns, initvalues){ } // calling the appropriate
                                             // base constructor
//##############################################################################
BJmatrix::BJmatrix(char *flag, int dimension ): // create an ident BJmatrix
  nDarray("NO")           // with base class constructor cancelation
  {
    if ( flag[0] != 'I' && flag[0] != 'e' )
     {
      ::printf("\nTo create a 2nd rank Kronecker delta type: nDarray (\"I\",2);\n");
      ::exit( 1 );
     }
// create the structure:
     pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
     rank(2);  //rank_of_nDarray here BJmatrix = 2

// in the case of nDarray_rank=0 add one to get right thing from the
// operator new
     set_dim_pointer( new int[rank()] );  // array for dimensions
     long int number = dimension*dimension;
     total_number(number);
     for( int idim = 1 ; idim <= rank() ; idim++ )
       {
         dim()[idim-1] = dimension;
       }

// allocate memory for the actual nDarray as nDarray
     set_data_pointer( new double [ (size_t) total_number() ]);
       if (!data())
         {
           ::printf("\a\nInsufficient memory for array\n");
           ::exit(1);
         }

     set_reference_count(+1);  // so far, there's one reference

     for ( int i2=1 ; i2<=dim(1) ; i2++ )
       for ( int j2=1 ; j2<=dim(2) ; j2++ )
         val(i2,j2) = (i2 == j2 ? 1  : 0);
}

//#############################################################################
// error message when trying to read a 'standard'
// BJmatrix file:
static char nonstandard[] =
"is a 'non-standard'file. A 'standard' BJmatrix file must\n"
"start with the dimensions of the BJmatrix, i.e.:\n"
"\t rows= 12 columns= 14\n or abbreviated: r= 12 c= 14\n"
"Notice rows appear before columns, and chars are lowercase\n"
"comments follow '!' signs to end of line, data follows :::\n";

//#############################################################################
// read from "standard" BJmatrix file:
BJmatrix::BJmatrix(char *initfile):
  nDarray("NO")           // with base class constructor cancelation
{
#define BSIZE 120
  FILE *from;

  int rows    = 0;
  int columns = 0;

  if ((from = fopen(initfile,"r")) == 0)
    error("cannot open BJmatrix initializer file", initfile);
  char buf[BSIZE], *cp, *cp2;
  int rfound = 0, cfound = 0, colonsfound = 0;
//        /*Parse file initialization header  */
  while( fgets(buf, BSIZE, from) )         // for each header line
    {
// Remove comments with ANSI C library function 'strp
      if( ( cp = strpbrk(buf,"!")) != NULL ) // look for comments
        *cp = '\0';  // terminate string at comment
      if(( cp = strpbrk(buf,"r")) != NULL )
        if (( cp2 = strpbrk(cp, "=")) != NULL )
          if (( cp = strpbrk( cp2, "0123456789")) != NULL )
            {
              rows = atoi(cp);
              rfound++;  // flag to say rows were found
            }
      if( ( cp = strpbrk(buf,"c") ) != NULL )
         if ( ( cp2 = strpbrk(cp, "=")) != NULL )
           if ( ( cp = strpbrk(cp2, "O123456789")) != NULL )
             {
               columns = atoi(cp);
               cfound++;  // flag to say cols were found
             }
      if ( strstr(buf,":::") != NULL )
        {
          colonsfound++;
          break; //... out of "while" loop
        }
    }
  if ( !rfound || !cfound || !colonsfound )
    {
      fprintf(stderr, " %s%s", initfile, nonstandard);
      exit(1);
    }

// let's construct a BJmatrix

// create the structure:
  pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
  rank(2);  // rank_of_nDarray here BJmatrix = 2

  set_dim_pointer( (new int[ rank() ]) );// array for dimensions
  long int number = rows*columns;
  total_number(number);

  dim()[0] = rows;
  dim()[1] = columns;

// allocate memory for the actual nDarray as nDarray
  set_data_pointer( new double [ (size_t) total_number() ] );
    if (!data())
      {
        ::printf("\a\nInsufficient memory for array\n");
        ::exit(1);
      }

  set_reference_count(+1);  // so far, there's one reference

  for ( int row=1 ; row <= rows ; row++ )
    for ( int col=1 ; col <= columns ; col++ )
      {
        char nb[20];
        fscanf(from, "%s", nb); // scan for space-delimited string
        val(row,col) = (double) atof(nb); // convert it to a double
        if(ferror(from))
          error("problem with BJmatrix initializer file ", initfile);
      }
  fclose(from);
}

//#############################################################################
// read from flat BJmatrix file: and write to test output
BJmatrix::BJmatrix(char *initfile, char *outfile):
  nDarray("NO")           // with base class constructor cancelation
{
//#define BSIZE 120
  FILE *from;
//Kai  FILE *to;

  int rows    = 0;
  int columns = 0;

  if ((from = fopen(initfile,"r")) == 0)
    error("cannot open BJmatrix initializer file", initfile);
//  char buf[BSIZE];
  char  cp[32];
  fscanf(from, "%s", cp); // scan for space-delimited string
//  cp = strpbrk( buf, "0123456789")
  rows = columns = atoi(cp);
// create the structure:
  pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
  rank(2);  // rank_of_nDarray here BJmatrix = 2

  set_dim_pointer( (new int[ rank() ]) );// array for dimensions
  long int number = rows*columns;
  total_number(number);

  dim()[0] = rows;
  dim()[1] = columns;

// allocate memory for the actual nDarray as nDarray
  set_data_pointer( new double [ (size_t) total_number() ] );
    if (!data())
      {
        ::printf("\a\nInsufficient memory for array\n");
        ::exit(1);
      }

  set_reference_count(+1);  // so far, there's one reference

  for ( int row=1 ; row <= rows ; row++ )
    for ( int col=1 ; col <= columns ; col++ )
      {
        char nb[20];
        fscanf(from, "%s", nb); // scan for space-delimited string
        val(row,col) = (double) atof(nb); // convert it to a double
        if(ferror(from))
          error("problem with BJmatrix initializer file ", initfile);
      }
  fclose(from);

//Kai  if ((to = fopen(outfile,"w")) == 0)
//Kai    error("cannot open BJmatrix output file", outfile);
  write_standard(outfile, "this is a test BJmatrix output"); 
}

//##############################################################################
BJmatrix::BJmatrix(const BJmatrix & x): // copy-initializer
  nDarray("NO")     // with base class constructor cancelation
    {
      x.pc_nDarray_rep->n++;  // we're adding another reference.
//      x.reference_count(+1); // we're adding another reference.
      pc_nDarray_rep = x.pc_nDarray_rep;  // point to the new BJtensor_rep.
    }



//##############################################################################
BJmatrix::BJmatrix(const nDarray & x):
  nDarray( x ) { }



//--// IT IS NOT INHERITED so must be defined in all derived classes
//--// See ARM page 277.
//--//##############################################################################
//--BJmatrix::~BJmatrix()
//--{
//--  if (reference_count(-1) == 0)  // if reference count  goes to 0
//--    {
//--// DEallocate memory of the actual nDarray
//--//    delete [pc_nDarray_rep->pc_nDarray_rep->total_numb] pc_nDarray_rep->pd_nDdata;
//--// nema potrebe za brojem clanova koji se brisu## see ELLIS & STROUSTRUP $18.3
//--//                                                and note on the p.65($5.3.4)
//--//  and the page 276 ($12.4)
//--    delete [] data();
//--    delete [] dim();
//--    delete pc_nDarray_rep;
//--  }
//--}


//#############################################################################
int BJmatrix::rows( ) const  // rows in BJmatrix
  {
    return this->dim(1);
  }



//#############################################################################
int BJmatrix::cols( ) const       // cols in BJmatrix
  {
    return this->dim(2);
  }



//#############################################################################
BJmatrix& BJmatrix::operator=( const BJmatrix & rval)
{
      rval.pc_nDarray_rep->n++;  // we're adding another reference.
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
    return *this;
}



//#############################################################################
void BJmatrix::write_standard(char *filename, char *msg)
  {
    FILE *to;
    if ((to = fopen(filename,"w")) == NULL)
      error("cannot open or create BJmatrix output file",filename);
    fprintf(to, "# %s: BJmatrix file written in \"standard\" format\n", filename);
    time_t clock;
    time(&clock);
    fprintf(to, "# %s", asctime(localtime(&clock)));
    fprintf(to, "# %s", msg);
    fprintf(to, "rows= %d columns= %d", rows(), cols());
    fprintf(to, ":::\n");
    for ( int row=0; row<rows() ; row++ )
      {
        for( int col=0; col<cols() ; col++ )
          {
            fprintf(to, "%6.6g  ",mval(row,col));
            if(ferror(to))
              error("problem with BJmatrix output file ", filename);
          }
        fprintf(to, "\n");
      }
}

//....//#############################################################################
//....// This is from JOOP May/June 1990 after ARKoenig ( I did the rest )
//....BJmatrix & BJmatrix::operator*=( const BJmatrix & rval )
//....  {
//....    if( this->cols() != rval.rows())
//....      error("# rows of second mat must equal "
//....               "# cols of first for multiply#");
//....
//....   int this_r = this->pc_nDarray_rep->dim[0];
//....   int this_c = this->pc_nDarray_rep->dim[1];
//....   int rval_r = rval.pc_nDarray_rep->dim[0];
//....   int rval_c = rval.pc_nDarray_rep->dim[1];
//....// temporary BJmatrix that will keep result
//....   BJmatrix result( this_r, rval_c, 0.0 );
//....
//....// copy *this if necessery
//....// i.e. if somebody else is pointing to the same nDarray_rep class
//....// then allocate new memory for the this and disconect this from
//....// the old one!
//....    if ( this->pc_nDarray_rep->n > 1 )// see ARK in JOOP may/june '90
//....      {                               // "Letter From a Newcomer"
//....//..//..............................................................................
//....//..// create the structure:
//....//..        nDarray_rep * New_pc_nDarray_rep = new nDarray_rep; // this 'new' is overloaded
//....//..        New_pc_nDarray_rep->nDarray_rank = this->pc_nDarray_rep->nDarray_rank;
//....//..// in the case of nDarray_rank=0 add one to get right thing from the
//....//..// operator new
//....//..        int one_or0 = 0;
//....//..        if(!New_pc_nDarray_rep->nDarray_rank) one_or0 = 1;
//....//..        New_pc_nDarray_rep->dim = new int[New_pc_nDarray_rep->nDarray_rank+one_or0];
//....//..                                  // array for dimensions
//....//..        New_pc_nDarray_rep->total_numb = 1;
//....//..        for( int idim = 0 ; idim < New_pc_nDarray_rep->nDarray_rank ; idim++ )
//....//..          {
//....//..            New_pc_nDarray_rep->dim[idim] = this->pc_nDarray_rep->dim[idim];
//....//..            New_pc_nDarray_rep->total_numb *= New_pc_nDarray_rep->dim[idim];
//....//..          }
//....//..// allocate memory for the actual nDarray as nDarray
//....//..        New_pc_nDarray_rep->pd_nDdata = new double [(size_t)New_pc_nDarray_rep->total_numb];
//....//..          if (!New_pc_nDarray_rep->pd_nDdata)
//....//..            {
//....//..              ::printf("\a\nInsufficient memory for array\n");
//....//..              ::exit(1);
//....//..            }
//....//..         New_pc_nDarray_rep->n = 1;  // so far, there's one reference
//....//..         for ( int i=0 ; i<New_pc_nDarray_rep->total_numb ; i++ )
//....//..           New_pc_nDarray_rep->pd_nDdata[i] = this->pc_nDarray_rep->pd_nDdata[i];
//....//.........
//....         this->pc_nDarray_rep->total_numb--;
//....         this->pc_nDarray_rep = New_pc_nDarray_rep;
//....//..............................................................................
//....      }
//....
//....// multiply rval with *this
//....   for ( int row=0 ; row<this_r ; row++ )
//....      for ( int col=0 ; col<rval_c ; col++ )
//....        {
//....          double sum = 0.0;
//....          for( int i=0 ; i<this_c ; i++ )
//....            sum += this->pc_nDarray_rep->m[row][i]*rval.pc_nDarray_rep->m[i][col];
//....          result.val(row,col) = sum;
//....        }
//....//..// copy result to *this
//....//..   for ( int i=0 ; i<this_r ; i++ )
//....//..      for ( int j=0 ; j<rval_c ; j++ )
//....//..         this->pc_nDarray_rep->m[i][j] = result.val(i,j);
//....// copy result to *this actually take the whole class pointer and just
//....// assign it to the this
//....
//....   this->pc_nDarray_rep = result.pc_nDarray_rep;
//....
//....   return *this;
//....  }
//....
//....//#############################################################################
//....// This is from JOOP May/June 1990 after ARKoenig
//....inline BJmatrix operator*(const BJmatrix & left_val, const BJmatrix & right_val)
//....  {
//....    BJmatrix result(left_value);
//....    result *= right_val;
//....    return result;
//....  }
//....
//....




//#############################################################################
// OLD good version
BJmatrix BJmatrix::operator*( BJmatrix & arg)
  {
    if( cols() != arg.rows())
      error("# rows of second mat must equal "
               "# cols of first for multiply#");
    BJmatrix result(rows(),arg.cols());
    for( int row=0 ; row<rows() ; row++ )
      for( int col=0 ; col<arg.cols() ; col++ )
        {
          double sum = 0;
          for( int i=0 ; i<cols() ; i++ )
            sum += mval(row,i)*arg.mval(i,col);
          result.mval(row,col) = sum;
        }
    return result; // Returning a local variable?
    // copy-initializer happens before the destructor,
    // so reference count is 2 when destructor is called,
    // thus destructor doesn't free the memory.
  }

//**sky**// this is COUPLING between ordinary BJmatrix and SKYMATRIX
//**sky**// it will be usefull in multiplying vector with skyBJmatrix
//**sky**BJmatrix BJmatrix::operator*(const skyBJmatrix & arg)
//**sky**  {
//**sky**    if( cols() != arg.dimension_of_sky_M())
//**sky**      error("# dimension_of_sky_M must equal "
//**sky**               "# cols of first for multiply#");
//**sky**    BJmatrix result(rows(),arg.dimension_of_sky_M());
//**sky**    for( int row=0 ; row<rows() ; row++ )
//**sky**      for( int col=0 ; col<arg.dimension_of_sky_M() ; col++ )
//**sky**        {
//**sky**          double sum = 0;
//**sky**          for( int i=0 ; i<cols() ; i++ )
//**sky**            {
//**sky**//            double prvi = val(row,i);
//**sky**//            double drugi = arg.full_val(i+1,col+1);
//**sky**//            double temp = prvi*drugi; // in skyBJmatrix
//**sky**//            sum += temp;
//**sky**              sum += val(row,i)*arg.full_val(i+1,col+1);
//**sky**            }
//**sky**          result.val(row,col) = sum;                   // start from 1
//**sky**        }
//**sky**    return result; // Returning a local variable?
//**sky**    // copy-initializer happens before the destructor,
//**sky**    // so reference count is 2 when destructor is called,
//**sky**    // thus destructor doesn't free the memory.
//**sky**  }
//**sky**



//#############################################################################
BJmatrix BJmatrix::operator*( double arg)
  {
    BJmatrix result(rows(),cols());
    for ( int i=0 ; i<rows() ; i++ )
     for ( int j=0 ; j<cols() ; j++ )
       result.mval(i,j) = mval(i,j) * arg;
    return result;
  }

//....//#############################################################################
//....// I had to overload the operator* from BJmatrix class
//....// because in the case of vectors s*sT and so on
//....vector BJmatrix::operator*( vector & arg)
//....  {
//....    if( cols() != arg.rows())
//....      error("# rows of second mat must equal "
//....               "# cols of first for multiply in BJmatrix::operator*( vector & arg) ");
//....    vector result(arg.rows());
//....      for( int row=0 ; row<arg.rows() ; row++ )
//....        {
//....          double sum = 0;
//....          for( int i=0 ; i<cols() ; i++ )
//....            sum += mval(row,i)*arg.cval(i);
//....          result.val(row) = sum;
//....        }
//....    return result; // Returning a local variable?
//....    // copy-initializer happens before the destructor,
//....    // so reference count is 2 when destructor is called,
//....    // thus destructor doesn't free the memory.
//....  }


//#############################################################################
BJmatrix BJmatrix::transpose()
  {
// This request is relaxed for vectors
//    if(rows() != cols())
//    error("BJmatrix must be square to transpose#\n");
    BJmatrix trans(cols(), rows());
    for ( int row=0 ; row<rows() ; row++ )
      for ( int col=0 ; col<cols() ; col++ )
        trans.mval(col,row) = mval(row,col);
    return trans;
  }

//#############################################################################
double BJmatrix::mmin()
  {
    double temp;
    if ( rows()<=0 || cols()<=0)
      error("bad BJmatrix size for min ()");
    double minimum = mval(0,0);
    for ( int row=0 ; row<rows() ; row++ )
      for ( int col=0 ; col<cols() ; col++ )
        if ( (temp=mval(row,col)) < minimum )
          minimum = temp;
    return minimum;
  }

//#############################################################################
double BJmatrix::mmax()
  {
    double temp;
    if( rows()<=0 || cols()<=0)
      error("bad BJmatrix size for max()");
    double maximum = mval(0,0);
    for ( int row=0 ; row<rows() ; row++ )
      for ( int col=0 ; col<cols() ; col++ )
        {
          if ( (temp=mval(row,col)) > maximum )
          maximum = temp;
        }
    return maximum;
  }


//#############################################################################
double BJmatrix::mean()
  {
    int col = 0;
    int row = 0;
    double sum = 0;
    for ( row=0 ; row<rows() ; row++ )
      for ( col=0 ; col<cols() ; col++ )
        sum += fabs(mval(row,col));
    return sum/(row*col);
  }


//#############################################################################
double BJmatrix::sum()
  {
    int col;
    double sum = 0;
    for ( int row=0 ; row<rows() ; row++ )
      for ( col=0 ; col<cols() ; col++ )
        sum += mval(row,col);
    return sum;
  }



//#############################################################################
double BJmatrix::variance()
  {
    int col = 0;
    int row = 0;
    double s_squared = 0;
    double mn = mean();
    for ( row=0 ; row<rows() ; row++ )
      for ( col=0 ; col<cols() ; col++ )
        {
          double temp = mval(row,col) - mn;
          temp *= temp;
          s_squared += temp;
        }
    s_squared /= row*col-1; // number of elements minus one
    return s_squared;
  }




//#############################################################################
double BJmatrix::determinant()
  {
    if( rows() != cols() )
      error("BJmatrix must be square for determinant()");
    BJmatrix indx(cols()); // create the "index vector"
//    BJmatrix B(cols()); // see pp 38. in Numerical Recipes
    int d;
    // perform the decomposition once:
    BJmatrix decomp = lu_decompose(indx,d);
//    decomp.print("decomposed BJmatrix");
    double determinant = d;
      for( int i=0 ; i < cols() ; i++ )
        determinant *= decomp.mval(i,i);
    return determinant;
  }


//#############################################################################
BJmatrix BJmatrix::inverse()
  {
    if( rows()!=cols() )
      error("BJmatrix must be square for inverse()");
    BJmatrix Y("I", rows()); // create an identity BJmatrix
    BJmatrix indx(cols());   // create the "index vector"
    BJmatrix B(cols());      // see Press & Flannery
    int d;
    // perform the decomposition once:
    BJmatrix decomp = lu_decompose(indx,d);
    for(int col=0 ; col<cols() ; col++ )
      {
        B.copy_column(Y,col,0);
        decomp.lu_back_subst(indx, B);
        Y.copy_column(B, 0, col);
      }
    return Y.transpose();
  }




//~~~~//##############################################################################
//~~~~BJtensor BJmatrix::BJmatrix2BJtensor_1()  // convert BJmatrix of to 4th order BJtensor
//~~~~  {                               // I_ijkl scheme
//~~~~
//~~~~// building up a BJtensor back_conv
//~~~~// this is special case just for 4th order with def_dim_4 !!
//~~~~  BJtensor back_conv( 4, def_dim_4, 0.0);
//~~~~
//~~~~  int m41 = 0;
//~~~~  int m42 = 0;
//~~~~
//~~~~// filling back the BJmatrix to BJtensor
//~~~~  for ( int c44=1 ; c44<=back_conv.dim(4) ; c44++ )
//~~~~    for ( int c43=1 ; c43<=back_conv.dim(3) ; c43++ )
//~~~~      for ( int c42=1 ; c42<=back_conv.dim(2) ; c42++ )
//~~~~        for ( int c41=1 ; c41<=back_conv.dim(1) ; c41++ )
//~~~~          {
//~~~~            m41 = back_conv.dim(1)*(c41-1)+c43;
//~~~~            m42 = back_conv.dim(2)*(c42-1)+c44;
//~~~~
//~~~~            back_conv.val(c41,c42,c43,c44) = this->val(m41,m42);
//~~~~          }
//~~~~
//~~~~//    back_conv.print("t","back to BJtensor back_conv");
//~~~~    return back_conv;
//~~~~  }

//##############################################################################
BJtensor BJmatrix::BJmatrix2BJtensor_2()  // convert BJmatrix of to 4th order BJtensor
  {                               // I_ikjl scheme

// building up a BJtensor back_conv
// this is special case just for 4th order with def_dim_4 !!
  BJtensor back_conv( 4, def_dim_4, 0.0);
#ifdef SASA
  BJtensor Back_conv(4,def_dim_4_2,0.0);
 if(cols()==4) back_conv=Back_conv; // dimenzije su po 2
#endif
  int m41 = 0;
  int m42 = 0;

// filling back the BJmatrix to BJtensor
  for ( int c44=1 ; c44<=back_conv.dim(4) ; c44++ )
    for ( int c43=1 ; c43<=back_conv.dim(3) ; c43++ )
      for ( int c42=1 ; c42<=back_conv.dim(2) ; c42++ )
        for ( int c41=1 ; c41<=back_conv.dim(1) ; c41++ )
          {
            m41 = back_conv.dim(1)*(c41-1)+c42;
            m42 = back_conv.dim(3)*(c43-1)+c44;

            back_conv.val(c41,c42,c43,c44) = this->val(m41,m42);
          }

//    back_conv.print("t","back to BJtensor back_conv");
    return back_conv;
  }

//##############################################################################
BJtensor BJmatrix::BJmatrix2BJtensor_22()  // convert BJmatrix of to 2th order BJtensor
  {                               // I_ikjl scheme

// building up a BJtensor back_conv
// this is special case just for 4th order with def_dim_4 !!

//  BJtensor back_conv( 2, def_dim_2, 0.0);
  BJtensor back_conv( 2, def_dim_2, this->data());

//..  //  int m41 = 0;
//..  //  int m42 = 0;
//..  
//..  // filling back the BJmatrix to BJtensor
//..     for ( int c22=1 ; c22<=back_conv.dim(2) ; c22++ )
//..       for ( int c21=1 ; c21<=back_conv.dim(1) ; c21++ )
//..         {
//..  //         m41 = back_conv.dim(1)*(c41-1)+c42;
//..  //         m42 = back_conv.dim(3)*(c43-1)+c44;
//..       
//..           back_conv.val(c21,c22) = this->val(m41,m42);
//..         }
//..  
//..  //    back_conv.print("t","back to BJtensor back_conv");

    return back_conv;
  }
//~~~~//##############################################################################
//~~~~BJtensor BJmatrix::BJmatrix2BJtensor_3()  // convert BJmatrix of to 4th order BJtensor
//~~~~  {                               // I_ilkj scheme
//~~~~
//~~~~// building up a BJtensor back_conv
//~~~~// this is special case just for 4th order with def_dim_4 !!
//~~~~  BJtensor back_conv( 4, def_dim_4, 0.0);
//~~~~
//~~~~  int m41 = 0;
//~~~~  int m42 = 0;
//~~~~
//~~~~// filling back the BJmatrix to BJtensor
//~~~~  for ( int c44=1 ; c44<=back_conv.dim(4) ; c44++ )
//~~~~    for ( int c43=1 ; c43<=back_conv.dim(3) ; c43++ )
//~~~~      for ( int c42=1 ; c42<=back_conv.dim(2) ; c42++ )
//~~~~        for ( int c41=1 ; c41<=back_conv.dim(1) ; c41++ )
//~~~~          {
//~~~~            m41 = back_conv.dim(1)*(c41-1)+c44;
//~~~~            m42 = back_conv.dim(2)*(c42-1)+c43;
//~~~~
//~~~~            back_conv.val(c41,c42,c43,c44) = this->val(m41,m42);
//~~~~          }
//~~~~
//~~~~//    back_conv.print("t","back to BJtensor back_conv");
//~~~~    return back_conv;
//~~~~  }
//~~~~
//~~~~
//#############################################################################
//##BJmatrix BJmatrix::compress_col(int col1, int col2, int to_col)
//##  {
//##    BJmatrix mat_temp( rows(), cols()-1);
//##//    mat_temp.print("mat_temp first appearence\n");
//##    for( int row=1 ; row<=rows() ; row++ )
//##      {
//##        // make temp vector to be put in the to_col :
//##// val STARTS from 0
//##        mat_temp(row,to_col) = val(row-1,col1-1)+val(row-1,col2-1);
//##      }
//##//    mat_temp.print("mat_temp second appearence\n");
//##
//##   int temp_row=1;
//##   int temp_col=1;
//##
//##    for( int col=1 ; col<=cols() ; col++ )
//##      {
//##        if( col!=col1 && col!=col2 )
//##          {
//##            temp_row=1;
//##            if( temp_col==to_col )
//##              {
//##                temp_col++;
//##              }
//##
//##          for( row=1 ; row<=rows() ; row++ )
//##            {
//##              mat_temp(temp_row,temp_col) = val(row-1,col-1);
//##              temp_row++;
//##            }
//##            temp_col++;
//##          }
//##      }
//##//    mat_temp.print("mat_temp third appearence\n");
//##
//##    return mat_temp;
//##  }
//##


//#############################################################################
//##BJmatrix BJmatrix::compress_row(int row1, int row2, int to_row)
//##  {
//##    BJmatrix mat_temp( rows()-1, cols());
//##//    mat_temp.print("mat_temp first appearence\n");
//##    for( int col=1 ; col<=cols() ; col++ )
//##      {
//##        // make temp vector to be put in the to_row :
//##// val STARTS from 0
//##        mat_temp(to_row,col) = val(row1-1,col-1)+val(row2-1,col-1);
//##      }
//##//    mat_temp.print("mat_temp second appearence\n");
//##
//##   int temp_col=1;
//##   int temp_row=1;
//##
//##    for( int row=1 ; row<=rows() ; row++ )
//##      {
//##        if( row!=row1 && row!=row2 )
//##          {
//##            temp_col=1;
//##            if( temp_row==to_row )
//##              {
//##                temp_row++;
//##              }
//##
//##          for( col=1 ; col<=cols() ; col++ )
//##            {
//##              mat_temp(temp_row,temp_col) = val(row-1,col-1);
//##              temp_col++;
//##            }
//##            temp_row++;
//##          }
//##      }
//##//    mat_temp.print("mat_temp third appearence\n");
//##
//##    return mat_temp;
//##  }
//##
//##


///*************************************************************
//The private support functions for determinant & inverse.
//**************************************************************/
//

//#############################################################################
// copy the from_col of mm to the to  col of "this"
void BJmatrix::copy_column( BJmatrix & mm, int from_col, int to_col)
  {
    if( rows()!=mm.rows() )
      error("number of rows must be equal for copy_column()");
    for( int row=0 ; row<rows() ; row++ )
      mval( row, to_col) = mm.mval(row,from_col);
  }


//#############################################################################
void BJmatrix::switch_columns(int col1, int col2 )
  {
    int row = 0;
    BJmatrix temp( rows() );
    for( row=0 ; row<rows() ; row++ )
      // temporarily store col 1:
      temp.mval(row,0) = mval(row,col1);
    for( row=0 ; row<rows() ; row++ )
      mval(row,col1 ) = mval(row,col2); // move col2 to col1
    for( row=0 ; row<rows() ; row++ )
      mval(row,col2) = temp.mval(row,0); // move temp to col2
  }


//#############################################################################
// make an image of a BJmatrix (used in L-U decomposition)
void BJmatrix::deepcopy(BJmatrix & from, BJmatrix & to)
  {
    if( from.rows()!=to.rows() || from.cols()!=to.cols() )
      error("matrices must be equal dimensions for deepcopy()");
    for( int row=0 ; row<from.rows() ; row++ )
      for( int col=0 ; col<from.cols() ; col++ )
        to.mval(row,col) = from.mval(row,col);
}



//#############################################################################
// scale a BJmatrix (used in L-U decomposition)
BJmatrix BJmatrix::scale( )
  {
//    ::fprintf(stderr,"BJmatrix BJmatrix::scale( )\n"); 
    double double_eps = d_macheps();
    double temp;
    if( rows()<=0 || cols()<=0 )
      error("bad BJmatrix size for scale()");
    if( rows()!=cols() )
      error("BJmatrix must be square for scale()");
    BJmatrix scale_vector(rows());
    for ( int col=0 ; col<cols() ; col++ )
      {
        double maximum = 0.0;
        for( int row=0 ; row<rows() ; row++ )
          if ( (temp=(double)fabs(mval(row,col))) > maximum )
            maximum = temp; // find max column magnitude in this row
          if( fabs(maximum) <= double_eps )
            error("singular BJmatrix in scale()");
          scale_vector.mval(col,0) = 1/maximum; // save the scaling
      }
    return scale_vector;
  }



//#############################################################################
BJmatrix BJmatrix::lu_decompose(BJmatrix & indx, int & d )
  {
///*
// Returns the L-U decomposition of a BJmatrix. indx is an output
// vector that records the row permutation effected by the
// partial pivoting, d is output as +-1 depending on whether the
// number of row interchanges was even or odd, respectively.
// This routine is used in combination with lu_back_subst to
// solve linear equations or invert a BJmatrix.
//*/
    if( rows()!=cols() )
      error("Matrix must be square to L-U decompose#\n");
    d = 1; // parity check
    int row=0, col=0, k=0, col_max=0; // counters
    double dum=0.0; // from the book - I don't know significance (Authors remark##)
    double sum=0.0;
    double maximum=0.0;
    BJmatrix lu_decomp(rows(),cols());
 // make a direct copy of the original BJmatrix.
    deepcopy( *this, lu_decomp );
    BJmatrix scale_vector = lu_decomp.scale(); // scale the BJmatrix

 // The loop over columns of Crout's method:
    for( row=0 ; row<rows() ; row++ )
      {
        if ( row>0 )        // eqn ±.3. 1 ± except for row=col:
          {
            for ( col=0 ; col <= row-1 ; col++ )
              {
                sum = lu_decomp.mval(row,col);
                if( col > 0 )
                  {
                    for( k = 0 ; k <= col-1 ; k++ )
                      sum -= lu_decomp.mval(row,k)*lu_decomp.mval(k, col);
                    lu_decomp.mval(row,col) = sum;
                  }
              }
          }

  // Initialize for the search for the largest pivot element
        maximum = 0;
  // i=j of eq 2.3.12 & i=j+ 1.. N of 2.3.13.
        for( col=row ; col <= cols()-1 ; col++ )
          {
            sum = lu_decomp.mval(row,col);
            if( row > 0 )
              {
                for( k=0 ; k <= row-1 ; k++ )
                  sum -=  lu_decomp.mval(k,col)*lu_decomp.mval(row,k);
                lu_decomp.mval(row,col) = sum;
              }
   // figure of merit for pivot:
            dum = scale_vector.mval(col,0)*fabs(sum);
            if ( dum >= maximum )   // is it better than the best so far?
              {
                col_max = col;
                maximum = dum;
              }
          }

  // Do we need to interchange rows?
        if( row != col_max)
          {
            lu_decomp.switch_columns(col_max,row); // Yes, do so ...
            d *= -1; //... and change the parity of d
  // also interchange the scale factor:
            dum = scale_vector.mval(col_max,0);
	    // scale_vector.mval(col_max,0) = scale_vector.mval(col,0); is wrong
	    // Bug! busted by Frank McKenna (fmckenna@ce.berkeley.edu)  

	    scale_vector.mval(col_max,0) = scale_vector.mval(row,0); // this is OK
            scale_vector.mval(row,0) = dum;
          }

        indx.mval(row,0) = col_max;
  // Now, finally, divide by the pivot element:
        if( row != rows() - 1 )
          {
            if( lu_decomp.mval(row,row) == 0 )
              lu_decomp.mval(row,row) = TINY;
       // If the pivot element is zero the BJmatrix is
       // singular (at least to the precision of the
       // algorithm).  For some applications on singular
       // matrices, it is desirable to substitute TINY for zero

            dum = 1/lu_decomp.mval(row,row);

            for( col=row+1 ; col <= cols()-1 ; col++ )
              lu_decomp.mval(row,col) *= dum;

          }

      }
    
    if( lu_decomp.mval(rows()-1, cols()-1 ) == 0 )
      lu_decomp.mval(rows()-1, cols()-1 ) = TINY;

    return lu_decomp;
  }




//#############################################################################
void BJmatrix::lu_back_subst(BJmatrix & indx, BJmatrix & b)
  {
///*
// Solves the set of N linear equations A*X = B.  Here "this"
// is the LU-decomposition of the BJmatrix A, determined by the
// routine lu_decompose(). Indx is input as the permutation
// vector returned  by lu_decompose().  B is input as the
// right-hand side vector B,  and returns with the solution
// vector X.  This routine takes into  account the possibility
// that B will begin with many zero elements,  so it is efficient
// for use in BJmatrix inversion.  See pp 36-37 in
// Press & Flannery.
//*/
    if( rows() != cols() )
      error ("non-square lu_decomp BJmatrix in lu_back_subst()");
    if( rows() != b.rows() )
      error("wrong size B vector passed to lu_back_subst()");
    if( rows() != indx.rows() )
      error("wrong size indx vector passed to lu_back_subst()");
    int row, col, ll;
    int ii = 0;
    double sum;
    for( col=0 ; col < cols() ; col++ )
      {
        ll = (int)indx.mval(col,0);
        sum = b.mval(ll,0);
        b.mval(ll,0) = b.mval(col,0);
        if ( ii >= 0 )
          for( row = ii ; row <= col-1 ; row++ )
            sum -= mval(row,col)*b.mval(row,0);
        else if( sum != 0 )
        ii = col;
        b.mval(col,0) = sum;
      }
    for( col = cols() -1 ; col >= 0 ; col-- )
      {
        sum = b.mval(col,0);
        if ( col < cols() - 1 )
          for ( row = col + 1 ; row <= rows()-1 ; row++ )
        sum -= mval(row,col)*b.mval(row,0);
  // store a component of the soln vector X:
        b.mval(col,0) = sum/mval(col,col);
      }
  }



//#############################################################################
double & BJmatrix::mval (int row, int col) // I am still keeping mval
  {                                      // operator for compatibility
     return( this->val(row+1,col+1) );   // with old BJmatrix class members
  }                                      // and they start from 0 !
    // used by BJmatrix functions which KNOW they aren't
    // exceeding the boundaries

// // Moved to  nDarray 14 Oct. 1996
// //#############################################################################
// vector BJmatrix::eigenvalues(void)
//   {
//     int rows = this->rows();
//     int cols = this->cols();
//     if ( rows != cols )
//       {
//         ::printf("rows!=cols in eigenvalues\n");
//         ::exit(1);
//       }
// 
//     vector EV((rows), 0.0);
// 
// // najbezbolinije da stvarno napravim dvodimenzioni niz pa da  ga kopiram u 'a'
// // PAZI oni u NRC rade kao u FORTRANU dakle nizovi od 1 do n
//     double ** a = new double *[rows+1];
//     if ( !a ) {::printf("memory exhausted for **a \n"); ::exit(1);}
//     for ( int i=0 ; i<(rows+1) ; i++ )
//       {
//         a[i] = new double[rows+1];
//         if ( !a[i] ) {::printf("memory exhausted for *a \n"); ::exit(1);}
//       }
//     double * d = new double [rows+1];
//     double * e = new double [rows+1];
// 
//     for ( int j=0 ; j<rows ; j++) 
//       for ( int k=0 ; k<rows ; k++) 
//         {
//           a[j+1][k+1] = this->cval(j+1,k+1);
//         }
// 
// 
// // Householder reduction of a real symmetric BJmatrix see NRC page 373
//     tred2( a, rows, d, e);
// // QL algorithm with implicit shift see NRC page 380
//     tqli( d, e, rows, a);
// // sort eigenvalues see NRC page 366
//     eigsrt( d, a, rows);
// 
//     for ( int l=0 ; l<rows ; l++ )
//       {
//         EV.val(l+1) = d[l+1];
//       }
// 
// // delete new - ed arrays
//     for ( int i1=0 ; i1<rows ; i1++ )
//       {
//         delete a[i1];
//       }
//     delete a;
//     delete d;
//     delete e;
// 
//     return EV;
//   }
// 
// //#############################################################################
// BJmatrix BJmatrix::eigenvectors(void)
//   {
//     int rows = this->rows();
//     int cols = this->cols();
//     if ( rows != cols )
//       {
//         ::printf("rows!=cols in eigenvectors\n");
//         ::exit(1);
//       }
// 
//     BJmatrix EV(rows, rows, 0.0);
// //    BJmatrix temp( rows, rows, rows, this->data() );
// 
// // najbezbolinije da stvarno napravim dvodimenzioni niz pa da  ga kopiram u 'a'
// // PAZI oni u NRC rade kao u FORTRANU dakle nizovi od 1 - n
//     double ** a = new double *[rows+1];
//     if ( !a ) {::printf("memory exhausted for **a \n"); ::exit(1);}
//     for ( int i=0 ; i<(rows+1) ; i++ )
//       {
//         a[i] = new double[rows+1];
//         if ( !a[i] ) {::printf("memory exhausted for *a \n"); ::exit(1);}
//       }
//     double * d = new double [rows+1];
//     double * e = new double [rows+1];
// 
//     for ( int j=0 ; j<rows ; j++) 
//       for ( int k=0 ; k<rows ; k++) 
//         {
//           a[j+1][k+1] = this->cval(j+1,k+1);
//         }
// // Householder reduction of a real symmetric BJmatrix see NRC page 373
//     tred2( a, rows, d, e);
// // QL algorithm with implicit shift see NRC page 380
//     tqli( d, e, rows, a);
// // sort eigenvalues see NRC page 366
//     eigsrt( d, a, rows);
// 
//     for ( int l=0 ; l<rows ; l++ )
//       for ( int l1=0; l1<rows ; l1++ )
//         {
//           EV.val(l+1,l1+1) = a[l+1][l1+1];
//         }
// 
// // delete new - ed arrays
//     for ( int i1=0 ; i1<(rows+1) ; i1++ )
//       {
//         delete a[i1];
//       }
//     delete a;
//     delete d;
//     delete e;
// 
//     return EV;
//   }
//


// // prebacen u nDarray 14 oktobra 1996
// //#############################################################################
// #define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))
// 
// void BJmatrix::tqli(double * d, double * e, int n, double ** z)
// {
//         int m,l,iter,i,k;
//         double s,r,p,g,f,dd,c,b;
// //	void nrerror();
// 
//         for (i=2;i<=n;i++) e[i-1]=e[i];
//         e[n]=0.0;
//         for (l=1;l<=n;l++) {
//                 iter=0;
//                 do {
//                         for (m=l;m<=n-1;m++) {
//                                 dd=fabs(d[m])+fabs(d[m+1]);
//                                 if (fabs(e[m])+dd == dd) break;
//                         }
//                         if (m != l) {
//                                 if (iter++ == 30) { ::printf("Too many iterations in TQLI\n"); ::exit(1); }
//                                 g=(d[l+1]-d[l])/(2.0*e[l]);
//                                 r=sqrt((g*g)+1.0);
//                                 g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
//                                 s=c=1.0;
//                                 p=0.0;
//                                 for (i=m-1;i>=l;i--) {
//                                         f=s*e[i];
//                                         b=c*e[i];
//                                         if (fabs(f) >= fabs(g)) {
//                                                 c=g/f;
//                                                 r=sqrt((c*c)+1.0);
//                                                 e[i+1]=f*r;
//                                                 c *= (s=1.0/r);
//                                         } else {
//                                                 s=f/g;
//                                                 r=sqrt((s*s)+1.0);
//                                                 e[i+1]=g*r;
//                                                 s *= (c=1.0/r);
//                                         }
//                                         g=d[i+1]-p;
//                                         r=(d[i]-g)*s+2.0*c*b;
//                                         p=s*r;
//                                         d[i+1]=g+p;
//                                         g=c*r-b;
//                                         /* Next loop can be omitted if eigenvectors not wanted */
//                                         for (k=1;k<=n;k++) {
//                                                 f=z[k][i+1];
//                                                 z[k][i+1]=s*z[k][i]+c*f;
//                                                 z[k][i]=c*z[k][i]-s*f;
//                                         }
//                                 }
//                                 d[l]=d[l]-p;
//                                 e[l]=g;
//                                 e[m]=0.0;
//                         }
//                 } while (m != l);
//         }
// }
// 
// 
// //#############################################################################
// void BJmatrix::tred2(double ** a, int n, double * d, double * e)
// {
//         int l,k,j,i;
//         double scale,hh,h,g,f;
// 
//         for (i=n;i>=2;i--) {
//                 l=i-1;
//                 h=scale=0.0;
//                 if (l > 1) {
//                         for (k=1;k<=l;k++)
//                                 scale += fabs(a[i][k]);
//                         if (scale == 0.0)
//                                 e[i]=a[i][l];
//                         else {
//                                 for (k=1;k<=l;k++) {
//                                         a[i][k] /= scale;
//                                         h += a[i][k]*a[i][k];
//                                 }
//                                 f=a[i][l];
//                                 g = f>0 ? -sqrt(h) : sqrt(h);
//                                 e[i]=scale*g;
//                                 h -= f*g;
//                                 a[i][l]=f-g;
//                                 f=0.0;
//                                 for (j=1;j<=l;j++) {
//                                 /* Next statement can be omitted if eigenvectors not wanted */
//                                         a[j][i]=a[i][j]/h;
//                                         g=0.0;
//                                         for (k=1;k<=j;k++)
//                                                 g += a[j][k]*a[i][k];
//                                         for (k=j+1;k<=l;k++)
//                                                 g += a[k][j]*a[i][k];
//                                         e[j]=g/h;
//                                         f += e[j]*a[i][j];
//                                 }
//                                 hh=f/(h+h);
//                                 for (j=1;j<=l;j++) {
//                                         f=a[i][j];
//                                         e[j]=g=e[j]-hh*f;
//                                         for (k=1;k<=j;k++)
//                                                 a[j][k] -= (f*e[k]+g*a[i][k]);
//                                 }
//                         }
//                 } else
//                         e[i]=a[i][l];
//                 d[i]=h;
//         }
//         /* Next statement can be omitted if eigenvectors not wanted */
//         d[1]=0.0;
//         e[1]=0.0;
//         /* Contents of this loop can be omitted if eigenvectors not
//                         wanted except for statement d[i]=a[i][i]; */
//         for (i=1;i<=n;i++) {
//                 l=i-1;
//                 if (d[i]) {
//                         for (j=1;j<=l;j++) {
//                                 g=0.0;
//                                 for (k=1;k<=l;k++)
//                                         g += a[i][k]*a[k][j];
//                                 for (k=1;k<=l;k++)
//                                         a[k][j] -= g*a[k][i];
//                         }
//                 }
//                 d[i]=a[i][i];
//                 a[i][i]=1.0;
//                 for (j=1;j<=l;j++) a[j][i]=a[i][j]=0.0;
//         }
// }
// 
// 
// 
// //#############################################################################
// void BJmatrix::eigsrt(double * d, double ** v, int n)
// {
//         int k,j,i;
//         double p;
// 
//         for (i=1;i<n;i++) {
//                 p=d[k=i];
//                 for (j=i+1;j<=n;j++)
//                         if (d[j] >= p) p=d[k=j];
//                 if (k != i) {
//                         d[k]=d[i];
//                         d[i]=p;
//                         for (j=1;j<=n;j++) {
//                                 p=v[j][i];
//                                 v[j][i]=v[j][k];
//                                 v[j][k]=p;
//                         }
//                 }
//         }
// }
// 

 
// //#############################################################################
double*  BJmatrix::BJmatrixtoarray(int& num)
  {
        num=pc_nDarray_rep->total_numb;
/*
        double *data = new double [num];

        for(int i=0 ; i< num ; i++)
        {
           data[i]=pc_nDarray_rep->pd_nDdata[i];
        }

        return data;
*/
        return pc_nDarray_rep->pd_nDdata;
  }



#endif 







