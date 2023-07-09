/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2008-01-16 00:19:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/TestDataOutputFileHandler.cpp,v $
                                                                        
                                                                        
// Written: fmk 11/04
// Revised:
//
// Purpose: This file is a driver to test the DataOutputHandler concrete subclasses.

#include <stdlib.h>

#include <OPS_Globals.h>
#include <Domain.h>
#include <Vector.h>
#include <FEM_ObjectBroker.h>
#include <StandardStream.h>

#include <DataOutputFileHandler.h>
#include <DataOutputFileHandler.h>

// global variables

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;

double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;


int main(int argc, char **argv)
{
  //
  // local variables used in test
  //
  int numColumns = 5;
  char *columns[numColumns];
  char column1[] = "col1";
  char column2[] = "col2";
  char column3[] = "col3";
  char column4[] = "col4";
  char column5[] = "col5";
  columns[0] = column1;
  columns[1] = column2;
  columns[2] = column3;
  columns[3] = column4;
  columns[4] = column5;
  int result, result1, result2, result3, result4;

  Vector data1(numColumns);
  Vector data2(numColumns-1);
  for (int i=0; i<numColumns; i++)
    data1(i) = i+1;
  for (int i=0; i<numColumns-1; i++)
    data2(i) = i+1;



  opserr << " *******************************************************************\n";
  opserr << "                        DataOutputFileHandler unit test\n";
  opserr << " *******************************************************************\n\n";
  opserr << " The DataOutputFile class provides implements 2 methods \n";
  opserr << "    open() - to open a connection to the o/p stream";
  opserr << "    write() - to send the contents of a vector to this stream\n\n";
  opserr << " *******************************************************************\n";

  //
  //  first test we can create DataOutputFile object to test.
  //

  opserr << "TEST: Creation of DataOutputFile object - will create 3 objects 1,2, & 3\n";
  
  DataOutputFileHandler *theHandler1  = new DataOutputFileHandler("UnitTest1");
  DataOutputFileHandler *theHandler2  = new DataOutputFileHandler("UnitTest2", DATA_FILE);
  DataOutputFileHandler *theHandler3  = new DataOutputFileHandler("UnitTest3", XML_FILE, APPEND);
  DataOutputFileHandler *theHandler4  = new DataOutputFileHandler(0);

  if (theHandler1 != 0 && theHandler2 != 0 && theHandler3 != 0 && theHandler4 != 0) 
    opserr << "PASS: Creation of DataOutputFileHandler object\n\n";    
  else {
    opserr << "FAIL: Creation of DataOutputFileHandler object\n\n";
    exit(-1);
  }

  //
  // test open() method
  //

  opserr << "Now We Test open()\n\t 1. invoke open() on our 3 objects (only 2nd will print the headings)\n";

  opserr << "Invoking on handler 1\n";
  result1 = theHandler1->open(columns, numColumns);
  opserr << "Invoking on handler 2\n";
  result2 = theHandler2->open(columns, numColumns);
  opserr << "Invoking on handler 3\n";
  result3 = theHandler3->open(columns, numColumns-1);
  
  if (result1 == 0 && result2 == 0 && result3 == 0) 
    opserr << "PASS: open() method  with VALID args\n\n";    
  else {
    opserr << "FAIL: open() method  with VALID args\n\n";    
  }

  //
  //  a number of tests to test open() method will fail if invalid args passed
  //     expected result: should fail, print out an error message and return != 0 for all following tests.
  //  

  opserr << "TEST: open() method with bad column descriptions\n";
  result = theHandler1->open(0, 2);
  if (result != 0) 
    opserr << "PASS: open() method with bad column descriptors\n\n"; 
  else {
    opserr << "FAIL: open() method with bad column descriptors\n\n"; 
  }

  opserr << "TEST: open() method with negative number of columns\n";
  result = theHandler1->open(columns, -1);
  if (result != 0) 
    opserr << "PASS: open() method with negative columns\n\n"; 
  else {
    opserr << "FAIL: open() method with negative columns\n\n"; 
  }


  //
  // test write() method
  //

  opserr << "Now We Test write()\n\t 1. invoke write() 3 times with valid args on both handler1, handler2 and handler3\n";
  opserr << "\t 2. invoke write() with invalid args a number of times\n\n";

  //
  //  test we can insert data into the file
  //     expected result: should return 0 and append line to <tableName>.out file
  //


  for (int i=0; i<3; i++) {
    result1 += theHandler1->write(data1);
    result1 += theHandler2->write(data1);
    result3 += theHandler3->write(data2);
    data1 += numColumns;
    data2 += numColumns-1;
  }
  
  // check the files exist and of correct size, UnitTest1.out should have 3 lines.
  if (result1 == 0 && result2 == 0 && result3 == 0) 
    opserr << "PASS: write() method  with VALID args\n\n";    
  else {
    opserr << "FAIL: write() method  with VALID args\n\n";    
  }

  //
  //  a number of tests to test the insert data will fail if data invalied
  //     expected result: all should return negative values & possible warning messages
  //

  opserr << "TEST: call write() method with incorrect vector size\n";
  result = theHandler3->write(data1);

  if (result != 0) 
    opserr << "PASS: write() method  with incorrect vector size\n\n";
  else {
    opserr << "FAIL: write() method  with incorrect vector size\n\n";
  }


  opserr << "TEST: call write() method on handler which has not been set\n";
  result = theHandler4->write(data1);

  if (result != 0) 
    opserr << "PASS: write() method to handler where open has not been called\n\n";
  else {
    opserr << "FAIL: write() method to handler where open has not been called\n\n";
  }

  //
  //  finally test we can destroy a DataOutputFileHandler object.
  //

  opserr << "TEST: Destruction of DataOutputFileHandler object\n";
  delete theHandler1;
  delete theHandler2;
  delete theHandler3;
  delete theHandler4;

  // if get to here without a segmentation fault destructor passes
  opserr << "PASS: Destruction of DataOutputFileHandler object\n\n";    
  exit(0);
}	
	
