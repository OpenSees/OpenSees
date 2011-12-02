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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/tagged/storage/main.cpp,v $
                                                                        
                                                                        
#include <Node.h>
#include <ArrayOfTaggedObjects.h>
#include <MapOfTaggedObjects.h>
#include <MapOfTaggedObjectsIter.h>

#include <stdlib.h>
#include <iostream.h>

int main(int argc, char **argv)
{
    int startID = 1;
    int i;

    MapOfTaggedObjects *theHolder = new MapOfTaggedObjects();
    cerr << "adding 5 components starting with id 1 to an object of initial size 5\n";

    for (i=0; i<5; i++) {
	Node *theNode = new Node(i+startID,1,i+startID);
	theHolder->addComponent(theNode);
    }
    
    cerr << "the contents\n";    
    theHolder->Print(cerr);

    TaggedObject *theNode;
    
    cerr << "removing node with id 4\n";
    theNode = theHolder->removeComponent(4); 
    if (theNode != 0) cerr << *theNode; else cerr << "Not There\n";        

    cerr << "the contents\n";    
    theHolder->Print(cerr);

    cerr << "adding 3 more components 10,11,12\n";        
    for (i=0; i<3; i++) {
	Node *theNode = new Node(i+10,1,i+10);
	theHolder->addComponent(theNode);
    }

    cerr << "the contents\n";    
    theHolder->Print(cerr);
    
    // destroying the holder
    cerr << "invoking the object destructor\n";        
    delete theHolder;
    
    cerr << "end of test1\n";
    double a; cin >> a;
    
    startID = 0;
    theHolder = new MapOfTaggedObjects();
    cerr << "adding 5 components starting with id 0 to an object of initial size 5\n";

    for (i=0; i<5; i++) {
	Node *theNode = new Node(i+startID,1,i+startID);
	theHolder->addComponent(theNode);
    }
    
    cerr << "the contents\n";    
    theHolder->Print(cerr);

    cerr << "removing node with id 4\n";
    theNode = theHolder->removeComponent(4); 
    if (theNode != 0) cerr << *theNode; else cerr << "Not There\n";        

    cerr << "the contents\n";    
    theHolder->Print(cerr);    
    
    cerr << "adding 3 more components 10,11,12\n";        
    for (i=0; i<3; i++) {
	Node *theNode = new Node(i+startID,1,i+10);
	theHolder->addComponent(theNode);
    }

    cerr << "the contents using print\n";    
    theHolder->Print(cerr);    

    cerr << "now check that the iterator works - should see contents\n";
    TaggedObject *theItem;
    TaggedObjectIter &theItems = theHolder->getComponents();
    while ((theItem = theItems()) != 0)
	cerr << *theItem;
    
    theHolder->clearAll();
    cerr << "the contetnts after clearAll\n";            
    theHolder->Print(cerr);
}
