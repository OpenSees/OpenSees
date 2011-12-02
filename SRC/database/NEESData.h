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
                                                                        
// $Revision: 1.3 $
// $Date: 2006-01-10 18:16:20 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/NEESData.h,v $
                                                                        
                                                                        
#ifndef NeesData_h
#define NeesData_h

// Written: fmk 
// Created: 03/04
//
// Description: This file contains the class definition for NeesData.
// NeesData is a concretes subclass of FE_Datastore. NEESData datastore
// is used to store data from the recorders in a file format that can
// be stored in the NEES data repository.
//
// What: "@(#) NeesData.h, revA"

#include <FE_Datastore.h>

#include <fstream>
using std::fstream;

#include <bool.h>

class NEES_table
{
 public:
  char     *name;
  char    **columns;
  int       numColumns;
  NEES_table    *next;
  bool      hasOutExtension;
};

class FEM_ObjectBroker;

class NEESData: public FE_Datastore
{
  public:
  /** @brief
   *
   * @param database name
   * @param domain
   * @param object broker
   *
   * @author fmk
   */
    NEESData(const char *dataBase,
	     Domain &theDomain, 
	     FEM_ObjectBroker &theBroker);    
    
    ~NEESData();

    // method to get a database tag
    int getDbTag(void);

    // methods for sending and receiving the data
    int sendMsg(int dbTag, int commitTag, 
		const Message &, 
		ChannelAddress *theAddress =0);    
    int recvMsg(int dbTag, int commitTag, 
		Message &, 
		ChannelAddress *theAddress =0);        

    int sendMatrix(int dbTag, int commitTag, 
		   const Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    int recvMatrix(int dbTag, int commitTag, 
		   Matrix &theMatrix, 
		   ChannelAddress *theAddress =0);
    
    int sendVector(int dbTag, int commitTag, 
		   const Vector &theVector, 
		   ChannelAddress *theAddress =0);
    int recvVector(int dbTag, int commitTag, 
		   Vector &theVector, 
		   ChannelAddress *theAddress =0);
    
    int sendID(int dbTag, int commitTag,
	       const ID &theID,
	       ChannelAddress *theAddress =0);
    int recvID(int dbTag, int commitTag,
	       ID &theID,
	       ChannelAddress *theAddress =0);

    int createTable(const char *tableName, int numColumns, char *columns[]);
    int insertData(const char *tableName, char *columns[], 
		   int commitTag, const Vector &data);

    int getData(const char *tableName, char *columns[], int commitTag, Vector &data);

    // the commitState method
    int commitState(int commitTag);        
    
  protected:

  private:
    char         *dataBase;
    NEES_table    *tables;
    int          numTables;
};



#endif

