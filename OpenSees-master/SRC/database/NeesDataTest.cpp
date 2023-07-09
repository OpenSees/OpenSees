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
                                                                        
// $Revision: 1.1 $
// $Date: 2004-07-13 23:13:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/NeesDataTest.cpp,v $
                                                                        
                                                                        
// Written: fmk 12/95
// Revised:
//
// Purpose: This file is a driver to unitetest the NEESData class.

#include <stdlib.h>

#include <OPS_Globals.h>
#include <NEESData.h>
#include <Domain.h>
#include <Vector.h>
#include <FEM_ObjectBroker.h>
#include <StandardStream.h>

// global variables
StandardStream sserr;
OPS_Stream &opserr = sserr;
double   ops_Dt =0;                
Domain  *ops_TheActiveDomain  =0;   
Element *ops_TheActiveElement =0;  

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
  int result, result1, result2;

  Vector data1(numColumns);
  Vector data2(numColumns-1);
  for (int i=0; i<numColumns; i++)
    data1(i) = i+1;
  for (int i=0; i<numColumns-1; i++)
    data2(i) = i+1;



  opserr << " *******************************************************************\n";
  opserr << "                        NEESData unit test\n";
  opserr << " *******************************************************************\n\n";
  opserr << " The NEESData class provides implements 2 methods \n";
  opserr << "    createTable() - to create a table, for NEESData this is 2 files,\n";
  opserr << "                    a <tableName>.out and a <tableName>.xml\n";
  opserr << "    insertData()  - to place data into the <tableName>.out file\n\n";
  opserr << " *******************************************************************\n";

  //
  //  first test we can create NEESData object to test.
  //

  opserr << "TEST: Creation of NEESData object\n";
  FEM_ObjectBroker theBroker;
  Domain *theDomain = new Domain();
  
  NEESData *theDatabase  = new NEESData("UnitTest",*theDomain,theBroker);
  if (theDatabase != 0) 
    opserr << "PASS: Creation of NEESData object\n\n";    
  else {
    opserr << "FAIL: Creation of NEESData object\n\n";
    exit(-1);
  }


  //
  // test createTable() method
  //

  opserr << "Now We Test createTable()\n\t 1. invoke createTable() twice to create tables UnitTest1 and UnitTest2\n";
  opserr << "\t 2. invoke createTable() with invalid args a number of times\n\n";

  //
  //  test we can create table; try it twice with different file names
  //     expected result: should return 0 and create 2 files for each createTable
  //

  opserr << "TEST: createTable() method  with VALID args (2 tables should be created, UnitTest1 and UnitTest2)\n";
  result1 = theDatabase->createTable("UnitTest1", numColumns, columns);
  result2 = theDatabase->createTable("UnitTest2", numColumns-1, columns);
  
  if (result1 == 0 && result2 == 0) 
    opserr << "PASS: createTable() method  with VALID args\n\n";    
  else {
    opserr << "FAIL: createTable() method  with VALID args\n\n";    
  }

  //
  //  a number of tests to test createTable() method will fail if invalid args passed
  //     expected result: should fail, print out an error message and return != 0 for all following tests.
  //  

  opserr << "TEST: createTable() method with bad fileName\n";
  result = theDatabase->createTable("/UnitTest", 3, columns);
  if (result != 0) 
    opserr << "PASS: createTable() method with bad fileName\n\n"; 
  else {
    opserr << "FAIL: createTable() method with bad fileName\n\n"; 
    opserr << " *** ARE YOU RUNNING THIS AS ROOOT -- YOU SHOULDN'T BE, su AND TRY AGAIN *** \n\n";    
  }

  opserr << "TEST: createTable() method with duplicate fileName\n";
  result = theDatabase->createTable("UnitTest1", 3, columns);
  if (result != 0) 
    opserr << "PASS: createTable() method with duplicate tableName\n\n"; 
  else {
    opserr << "FAIL: createTable() method with duplicate tableName\n\n"; 
  }

  opserr << "TEST: createTable() method with 0 columns\n";
  result = theDatabase->createTable("UnitTestFail", 0, columns);
  if (result != 0) 
    opserr << "PASS: createTable() method with 0 columns\n\n"; 
  else {
    opserr << "FAIL: createTable() method with 0 columns\n\n"; 
  }


  //
  // test insertData() method
  //

  opserr << "Now We Test insertData()\n\t 1. invoke insertData() 3 times with valid args on both UnitTest1 and UnitTest2\n";
  opserr << "\t 2. invoke insertData() with invalid args a number of times\n\n";

  //
  //  test we can insert data into the file
  //     expected result: should return 0 and append line to <tableName>.out file
  //

  opserr << "TEST: call insertData() method  with VALID args 3 times on both UNitTest1 and UnitTest2\n";
  for (int i=0; i<3; i++) {
    result1 = theDatabase->insertData("UnitTest1", columns, 1, data1);
    result2 = theDatabase->insertData("UnitTest2", columns, 1, data2);
    data1 += numColumns;
    data2 += numColumns-1;
  }
  
  // check the files exist and of correct size, UnitTest1.out should have 3 lines.
  if (result1 == 0 && result2 == 0) {
    
  }
  if (result1 == 0 && result2 == 0) 
    opserr << "PASS: insertData() method  with VALID args\n\n";    
  else {
    opserr << "FAIL: insertData() method  with VALID args\n\n";    
  }

  //
  //  a number of tests to test the insert data will fail if data invalied
  //     expected result: all should retuen != 0 and <tableName>.out should be unmodified
  //

  opserr << "TEST: call insertData() method with incorrect vector size\n";
  result = theDatabase->insertData("UnitTest2", columns, 1, data1);

  if (result != 0) 
    opserr << "PASS: insertData() method  with incorrect vector size\n\n";
  else {
    opserr << "FAIL: insertData() method  with incorrect vector size\n\n";
  }

  opserr << "TEST: call insertData() method to a table not yet created\n";
  result = theDatabase->insertData("UnitTest", columns, 1, data1);

  if (result != 0) 
    opserr << "PASS: insertData() method to a table not yet created\n\n";
  else {
    opserr << "FAIL: insertData() method to a table not yet created\n\n";
  }

  //
  //  finally test we can destroy a NEESData object.
  //

  opserr << "TEST: Destruction of NEESData object\n";
  delete theDatabase;

  // if get to here without a segmentation fault destructor passes
  opserr << "PASS: Destruction of NEESData object\n\n";    
  exit(0);
}	
	
