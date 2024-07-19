//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
// Author: fmk, cmp
//
#include <Logging.h>
#include <Parsing.h>
#include <Domain.h>
#include <Node.h>
#include <MP_Constraint.h>
#include <ID.h>
#include <Matrix.h>
#include <Vector.h>
//
#include <assert.h>
#include <tcl.h>
#include <BasicModelBuilder.h>

#define CONSTRAINT_ERROR TCL_ERROR
#define CONSTRAINT_OK    TCL_OK

static int createLinearRigidBeam(Domain &theDomain, int ret_tag, int con_tag);
static int createLinearRigidRod (Domain &theDomain, int ret_tag, int con_tag);
static int createLinearRigidDiaphragm(Domain &theDomain, int ret_tag, ID &nC, int perpPlaneConstrained);


//
// Written: fmk 12/99
//
static int
createLinearRigidBeam(Domain &theDomain, int ret_tag, int con_tag)
{
 
    // get a pointer to the retained and constrained nodes - make sure they exist
    Node *nodeR = theDomain.getNode(ret_tag);
    if (nodeR == nullptr) {
      opserr << G3_ERROR_PROMPT 
             << "retained Node" <<  ret_tag <<  " not in domain\n";
      return CONSTRAINT_ERROR;
    }

    Node *nodeC = theDomain.getNode(con_tag);
    if (nodeR == nullptr) {
      opserr << G3_ERROR_PROMPT 
             << "constrained node " <<  con_tag <<  " not in domain\n";
      return CONSTRAINT_ERROR;
    }

    // get the coordinates of the two nodes - check dimensions are the same
    const Vector &crdR = nodeR->getCrds();
    const Vector &crdC = nodeC->getCrds();
    int dimR = crdR.Size();
    int dimC = crdC.Size();
    if (dimR != dimC) {
      opserr << G3_ERROR_PROMPT 
             << "mismatch in dimension between constrained node " 
             <<  con_tag <<  " and retained node " << ret_tag << "\n";
      return CONSTRAINT_ERROR;
    }
    
    // check the number of dof at each node is the same
    int numDOF = nodeR->getNumberDOF();
    if (numDOF != nodeC->getNumberDOF()){ 
      opserr << G3_ERROR_PROMPT 
             << "mismatch in numDOF between constrained node " 
             <<  con_tag <<  " and retained node " << ret_tag << "\n";
      return CONSTRAINT_ERROR;
    }

    // check the number of dof at the nodes >= dimension of problem
    if(numDOF < dimR){    
      opserr << G3_ERROR_PROMPT 
             << "numDOF at nodes " << ret_tag << " and " 
             <<  con_tag <<  " must be >= dimension of problem\n";
      return CONSTRAINT_ERROR;
    }
    
    // create the ID to identify the constrained dof 
    ID id(numDOF);

    // construct the transformation matrix Ccr, where  Uc = Ccr Ur & set the diag, Ccr = I
    Matrix mat(numDOF,numDOF);
    mat.Zero();

    // set the values
    for (int i=0; i<numDOF; ++i) {
      mat(i,i) = 1.0;
      id(i) = i;
    }

    // if there are rotational dof - we must modify Ccr DONE ASSUMING SMALL ROTATIONS
    if (dimR != numDOF) {
      if (dimR == 2 && numDOF == 3) {
	double deltaX = crdC(0) - crdR(0);
	double deltaY = crdC(1) - crdR(1);	    
	mat(0,2) = -deltaY;
	mat(1,2) = deltaX;

      } else if (dimR == 3 && numDOF == 6) {
	double deltaX = crdC(0) - crdR(0);
	double deltaY = crdC(1) - crdR(1);	    
	double deltaZ = crdC(2) - crdR(2);

	// rotation about z/3 axis
	mat(0,5) = -deltaY;
	mat(1,5) =  deltaX;	

	// rotation about y/2 axis
	mat(0,4) =  deltaZ;
	mat(2,4) = -deltaX;

	// rotation about x/1 axis
	mat(1,3) = -deltaZ;
	mat(2,3) =  deltaY;

      } else { // not valid
	opserr << G3_ERROR_PROMPT 
               << " for nodes " << ret_tag << " and " << con_tag 
               << " nodes do not have valid numDOF for their dimension\n";
	return CONSTRAINT_ERROR;
      }
    }	
    
    // create the MP_Constraint
    MP_Constraint *newC = new MP_Constraint(ret_tag, con_tag, mat, id, id);

    // add the constraint to the domain
    if (theDomain.addMP_Constraint(newC) == false) {
      delete newC;
      opserr << G3_ERROR_PROMPT 
             << "nodes " << con_tag << " and " << ret_tag << ", could not add to domain\n";
      return CONSTRAINT_ERROR;
    }

    return CONSTRAINT_OK;
}

// File: ~/model/constraints/RigidRod.C
//
// Written: fmk 12/99
//
static int
createLinearRigidRod(Domain &theDomain, int ret_tag, int con_tag)
{

    // get a pointer to the retained node and constrained nodes - ensure these exist
    Node *nodeR = theDomain.getNode(ret_tag);
    if (nodeR == 0) {
      opserr << G3_ERROR_PROMPT
             << "retained node " <<  ret_tag <<  " not in domain\n";
      return CONSTRAINT_ERROR;
    }
    Node *nodeC = theDomain.getNode(con_tag);
    if (nodeR == 0) {
      opserr << G3_ERROR_PROMPT 
             << "constrained node " <<  con_tag <<  " not in domain\n";
      return CONSTRAINT_ERROR;
    }

    // get the coordinates of the two nodes - check dimensions are the same
    const Vector &crdR = nodeR->getCrds();
    const Vector &crdC = nodeC->getCrds();
    int dimR = crdR.Size();
    int dimC = crdC.Size();
    if (dimR != dimC) {
      opserr << G3_ERROR_PROMPT 
             << "mismatch in dimension between constrained node " 
             <<  con_tag <<  " and retained node " << ret_tag << "\n";
      return CONSTRAINT_ERROR;
    }
    
    // check the number of dof at each node is the same 
    int numDOF = nodeR->getNumberDOF();
    if (numDOF != nodeC->getNumberDOF()){ 
      opserr << G3_ERROR_PROMPT 
             << "mismatch in numDOF " << "between constrained node " 
             <<  con_tag <<  " and retained node " << ret_tag << "\n";
      return CONSTRAINT_ERROR;
    }

    // check the number of dof at the nodes >= dimension of problem
    if(numDOF < dimR){    
      opserr << G3_ERROR_PROMPT 
             << "numDOF at nodes " << ret_tag << " and " << con_tag 
             << " must be >= dimension of problem\n";
      return CONSTRAINT_ERROR;
    }
 
    // create the ID to identify the constrained dof 
    ID id(dimR);

    // construct the transformation matrix Ccr, where  Uc = Ccr Ur & set the diag
    Matrix mat(dimR,dimR);
    mat.Zero();

    // set the values
    for (int i=0; i<dimR; ++i) {
      mat(i,i) = 1.0;
      id(i) = i;
    }

    // create the MP_Constraint
    MP_Constraint *newC = new MP_Constraint(ret_tag, con_tag, mat, id, id);

    // add the constraint to the domain
    if (theDomain.addMP_Constraint(newC) == false) {
      delete newC;
      opserr << G3_ERROR_PROMPT << "for nodes " << con_tag << " and " << ret_tag << " could not add to domain\n";
      return CONSTRAINT_ERROR;
    }
    return CONSTRAINT_OK;
}


//
// Written: fmk 1/99
//
static int
createLinearRigidDiaphragm(Domain &theDomain, int ret_tag, ID &nC, 
			       int perpPlaneConstrained)
{
    // check plane is valid, i.e. perpPlaneConstrained must be 0, 1 or 2
    if (perpPlaneConstrained < 0 || perpPlaneConstrained > 2) {
      opserr << G3_ERROR_PROMPT << 
	"the dirn of perpendicular to constrained plane " << perpPlaneConstrained <<  " not valid\n";
      return CONSTRAINT_ERROR;
    }

    // check constrainedNodes ID does not contain the retained node
    if (nC.getLocation(ret_tag) >= 0) {
      opserr << G3_ERROR_PROMPT 
             << "retained node " << ret_tag << " is in constrained node list\n";
      return CONSTRAINT_ERROR;
    }	
    
    // get a pointer to the retained node and check node in 3d with 6 DOFs
    Node *nodeR = theDomain.getNode(ret_tag);
    if (nodeR == nullptr) {
      opserr << G3_ERROR_PROMPT 
             << "retained Node " <<  ret_tag <<  " not in domain\n";
      return CONSTRAINT_ERROR;
    }

    const Vector &crdR = nodeR->getCrds();
    if ((nodeR->getNumberDOF() != 6) || (crdR.Size() != 3)){
      opserr << G3_ERROR_PROMPT 
             << "retained Node " << ret_tag << " not in 3d space with 6 DOFs\n";
      return CONSTRAINT_ERROR;
    }	

    // 
    // create some objects which will be passed to the MP_Constraint 
    // constructor, elements of objects are filled in later
    // 
    // create the ID to identify the constrained DOFs 
    ID id(3);

    // construct the transformation matrix Ccr, where  Uc = Ccr Ur & set the diag
    Matrix mat(3,3);
    mat.Zero();
    mat(0,0) = 1.0; 
    mat(1,1) = 1.0; 
    mat(2,2) = 1.0;


    // now for each of the specified constrained DOFs we:
    // 1. check it's in the plane of the retained node, 
    // 2. set the ID and transformation matrix,
    // 3. create the MP_Constrainet and add it to the domain

    for (int i=0; i<nC.Size(); ++i) {

      // get the constrained node
      int ndC = nC(i);
      Node *nodeC = theDomain.getNode(ndC);

      // ensure node exists
      if (nodeC == nullptr) {
        // TODO: Document change; this used to be accepted without error.
        opserr << G3_ERROR_PROMPT 
               << "cannot constrain node " << ndC << " as no node in domain\n";
        return CONSTRAINT_ERROR;
      }

      // get node coordinates
      const Vector &crdC = nodeC->getCrds();

      // check constrained node has correct dim and number of DOFs
      if ((nodeR->getNumberDOF() == 6) && (crdR.Size() == 3)){

        // determine delta Coordintaes
        double deltaX = crdC(0) - crdR(0);
        double deltaY = crdC(1) - crdR(1);	    
        double deltaZ = crdC(2) - crdR(2);
        
        // rigid diaphragm in xy plane
        if (perpPlaneConstrained == 2) { 

          // check constrained node in xy plane with retained node
          if (deltaZ == 0.0) {

            // DOFs corresponding to dX, dY and theta Z (0,1,5)
            id(0) = 0; id(1) = 1; id(2) = 5;

            // set up transformation matrix
            mat(0,2) = - deltaY;
            mat(1,2) = deltaX;

          } else {
            opserr << G3_ERROR_PROMPT 
                   << "ignoring constrained node " << ndC << ", not in xy plane\n";
            return CONSTRAINT_ERROR;
          }

        // rigid diaphragm in xz plane
        } else if (perpPlaneConstrained == 1) { 

          // check constrained node in xy plane with retained node
          if (deltaY == 0.0) {

            // DOFs corresponding to dX, dZ and theta Y (0,2,4)
            id(0) = 0; id(1) = 2; id(2) = 4;

            // set up transformation matrix
            mat(0,2) = deltaZ;
            mat(1,2) = -deltaX;

          } else {
            opserr << G3_ERROR_PROMPT 
                   << "ignoring constrained node " << ndC << ", not in xz plane\n";
            return CONSTRAINT_ERROR;
          }

        // rigid diaphragm in yz plane
        } else {	  

          // check constrained node in xy plane with retained node
          if (deltaX == 0.0) {

            // DOFs corresponding to dY, dZ and theta X (1,2,3)
            id(0) = 1; id(1) = 2; id(2) = 3;

            // set up transformation matrix
            mat(0,2) = -deltaZ;
            mat(1,2) = deltaY;

          } else {
            opserr << G3_ERROR_PROMPT 
                   << "ignoring constrained node " << ndC << ", not in xz plane\n";
            return CONSTRAINT_ERROR;
          }
        }
            
        // create the MP_Constraint
        MP_Constraint *newC = new MP_Constraint(ret_tag, ndC, mat, id, id);

        // add the constraint to the domain
        if (theDomain.addMP_Constraint(newC) == false) {
          opserr << G3_ERROR_PROMPT 
                 << "ignoring constrained node " << ndC << ", failed to add\n";
          delete newC;
          return CONSTRAINT_ERROR;
        }

      } else  {
        opserr << G3_WARN_PROMPT 
               << "ignoring constrained node  " << ndC << ", not 3D node\n";
        return CONSTRAINT_OK;
      }
      
    } // for each node in constrained nodes
    return CONSTRAINT_OK;
}

int
TclCommand_RigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  // TODO: Change RigidDiaphragm to take Domain as clientData
  Domain *theTclDomain = ((BasicModelBuilder*)clientData)->getDomain();

  if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "rigidLink perpDirn? rNode? <cNodes?>\n";
      return TCL_ERROR;
  }

  int rNode, perpDirn;
  if (Tcl_GetInt(interp, argv[1], &perpDirn) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "rigidLink perpDirn rNode cNodes - could not read perpDirn? \n";
      return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "rigidLink perpDirn rNode cNodes - could not read rNode \n";
      return TCL_ERROR;
  }

  // read in the constrained nodes
  int numConstrainedNodes = argc - 3;
  ID constrainedNodes(numConstrainedNodes);
  for (int i=0; i<numConstrainedNodes; ++i) {
      int cNode;
      if (Tcl_GetInt(interp, argv[3+i], &cNode) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "rigidLink perpDirn rNode cNodes - could not read a cNode\n";
          return TCL_ERROR;
      }
      constrainedNodes(i) = cNode;
  }

  //RigidDiaphragm theLink(*theTclDomain, rNode, constrainedNodes, perpDirn-1);
  return createLinearRigidDiaphragm(*theTclDomain, rNode, constrainedNodes, perpDirn-1);
}



int
TclCommand_RigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theTclDomain = ((BasicModelBuilder*)clientData)->getDomain();

  if (argc < 4) {
      opserr << G3_ERROR_PROMPT << "rigidLink linkType? rNode? cNode?\n";
      return TCL_ERROR;
  }

  int rNode, cNode;
  if (Tcl_GetInt(interp, argv[2], &rNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "rigidLink linkType? rNode? cNode? - could not read rNode \n"; 
      return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &cNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "rigidLink linkType? rNode? cNode? - could not read CNode \n"; 
      return TCL_ERROR;
  }

  // construct a rigid rod or beam depending on 1st arg
  if ((strcmp(argv[1],"-bar") == 0) || (strcmp(argv[1],"bar") == 0)) {
    // RigidRod theLink(*theTclDomain, rNode, cNode);
    createLinearRigidRod(*theTclDomain, rNode, cNode);

  } else if ((strcmp(argv[1],"-beam") == 0) || (strcmp(argv[1],"beam") == 0)) {
    createLinearRigidBeam(*theTclDomain, rNode, cNode);
    //RigidBeam theLink(*theTclDomain, rNode, cNode);

  } else {
      opserr << G3_ERROR_PROMPT << "rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) \n"; 
      return TCL_ERROR;
  }
  return TCL_OK;
}

