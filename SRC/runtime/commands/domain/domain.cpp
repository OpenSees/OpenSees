//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Description: Domain manipulation commands which do not require
// additional namespacing.
//
#include <string.h>
#ifdef WIN32
  #define strdup _strdup
#endif
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>     // std::sort

#include <tcl.h>
#include <FileStream.h>
#include <Logging.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>

#include <Pressure_Constraint.h>
#include <Element.h>
#include <ElementIter.h>

#include <Node.h>
#include <Vector.h>

#define MAX_NDF 6

namespace OpenSees {
namespace DomainCommands {
int
domainChange(ClientData clientData, Tcl_Interp *interp, int argc,
             Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  ((Domain*)clientData)->domainChange();
  return TCL_OK;
}


int
removeObject(ClientData clientData, Tcl_Interp *interp, int argc,
             Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;

  // make sure at least one other argument to contain type of object
  if (argc < 2) {
    opserr << "WARNING want - remove objectType?\n";
    return TCL_ERROR;
  }
  const char* remove_type = Tcl_GetString(objv[1]);

  int tag;
  if ((strcmp(remove_type, "element") == 0) || 
      (strcmp(remove_type, "ele") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove element eleTag?\n";
      return TCL_ERROR;
    }

    if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove element tag? failed to read tag: " 
             << Tcl_GetString(objv[2]) << "\n";
      return TCL_ERROR;
    }
    Element *theEle = the_domain->removeElement(tag);
    if (theEle != nullptr) {

      // finally invoke the destructor on the element
      delete theEle;
    }
  }

  else if ((strcmp(remove_type, "loadPattern") == 0) ||
           (strcmp(remove_type, "pattern") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove loadPattern patternTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove loadPattern tag? failed to read tag: "
             << Tcl_GetString(objv[2]) << "\n";
      return TCL_ERROR;
    }
    LoadPattern *thePattern = the_domain->removeLoadPattern(tag);
    if (thePattern != nullptr) {
      thePattern->clearAll();
      delete thePattern;
    }
  }

  else if (strcmp(remove_type, "parameter") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove parameter paramTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove parameter tag? failed to read tag: " << Tcl_GetString(objv[2])
             << "\n";
      return TCL_ERROR;
    }
    Parameter *theParameter = the_domain->removeParameter(tag);
    if (theParameter != nullptr) {
      delete theParameter;
    }
  }

  else if (strcmp(remove_type, "node") == 0) {
    if (argc < 3) {
      opserr << "WARNING want - remove node nodeTag?\n";
      return TCL_ERROR;
    }
    if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
      opserr << "WARNING remove node tag? failed to read tag: " 
             << Tcl_GetString(objv[2]) << "\n";
      return TCL_ERROR;
    }
    Node *theNode = the_domain->removeNode(tag);
    if (theNode != nullptr) {
      delete theNode;
    }
    Pressure_Constraint *thePC = the_domain->removePressure_Constraint(tag);
    if (thePC != nullptr) {
      delete thePC;
    }
  }

  else if (strcmp(remove_type, "recorders") == 0) {
    the_domain->removeRecorders();
  }

  else if (strcmp(remove_type, "recorder") == 0) {
    if (argc < 3) {
      opserr << OpenSees::PromptValueError << "want - remove recorder recorderTag?\n";
      return TCL_ERROR;
    }
    int tag;
    if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "remove recorder tag? failed to read tag: " << Tcl_GetString(objv[2])
               << "\n";
      return TCL_ERROR;
    }
    if (the_domain->removeRecorder(tag) != 0) {
      opserr << OpenSees::PromptValueError << "No recorder found with tag " << tag << "\n";
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  else if ((strcmp(remove_type, "SPconstraint") == 0) ||
           (strcmp(remove_type, "sp") == 0)) {

    if (argc < 3) {
      opserr << "WARNING want - remove SPconstraint spTag? -or- remove SPconstraint nodeTag? dofTag? <patternTag?>\n";
      return TCL_ERROR;
    }    
    if (argc == 3) {
      if (Tcl_GetIntFromObj(interp, objv[2], &tag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read tag: " << objv[2] << "\n";
        return TCL_ERROR;
      }      

      SP_Constraint *theSPconstraint = the_domain->removeSP_Constraint(tag);
      if (theSPconstraint != nullptr)
        delete theSPconstraint;

    } else {
      int nodeTag, dofTag;
      int patternTag = -1;
      
      if (Tcl_GetIntFromObj(interp, objv[2], &nodeTag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read node tag: " << objv[2] << "\n";
        return TCL_ERROR;
      }      
      if (Tcl_GetIntFromObj(interp, objv[3], &dofTag) != TCL_OK) {
        opserr << "WARNING remove sp tag? failed to read dof tag: " << objv[3] << "\n";
        return TCL_ERROR;
      }      
      
      if (argc == 5) {
        if (Tcl_GetIntFromObj(interp, objv[4], &patternTag) != TCL_OK) {
          opserr << "WARNING remove sp tag? failed to read pattern tag: " << objv[4] << "\n";
          return TCL_ERROR;
        }
      }
      dofTag--;  // one for C++ indexing of dof
      
      the_domain->removeSP_Constraint(nodeTag, dofTag, patternTag);

      return TCL_OK;
    }
    //

//  const char** const args = new const char*[argc+1];
//  args[0] = objv[1];
//  args[1] = objv[0];
//  args[0] = strdup(objv[1]);
//  args[1] = strdup(objv[0]);
//  opserr << args[0] << " " << args[1] << " ";
//  for (int i=2; i<argc; ++i) {
//    args[i] = strdup(objv[i]);
//    opserr << args[i] << " ";
//  }
//  args[argc] = nullptr;
 // opserr << "\n";

//  Tcl_CmdInfo info;
//  assert(Tcl_GetCommandInfo(interp, args[0], &info) == 1);
//  int status = info.proc(info.clientData, interp, argc, args);
//  for (int i = 0; i < argc; ++i)
//    free((void*)args[i]);
//  delete[] args;
//  return status;
  }

  else if ((strcmp(Tcl_GetString(objv[1]), "MPconstraint") == 0) ||
           (strcmp(Tcl_GetString(objv[1]), "mp") == 0)) {
    if (argc < 3) {
      opserr << "WARNING want - remove MPconstraint nNodeTag? -or- remove "
                "MPconstraint -tag mpTag\n";
      return TCL_ERROR;
    }
    int nodTag = 0;
    if (argc == 3) {
      if (Tcl_GetIntFromObj(interp, objv[2], &nodTag) != TCL_OK) {
        opserr << "WARNING remove mp nodeTag? failed to read nodeTag: "
               << Tcl_GetString(objv[2]) << "\n";
        return TCL_ERROR;
      }

      the_domain->removeMP_Constraints(nodTag);
      return TCL_OK;
    }
    if (strcmp(Tcl_GetString(objv[2]), "-tag") == 0 && argc > 3) {
      if (Tcl_GetIntFromObj(interp, objv[3], &nodTag) != TCL_OK) {
        opserr << "WARNING remove mp -tag mpTag? failed to read mpTag: "
               << Tcl_GetString(objv[3]) << "\n";
        return TCL_ERROR;
      }

      the_domain->removeMP_Constraint(nodTag);
      return TCL_OK;
    }
  }

  else
    opserr << "WARNING remove " 
           << Tcl_GetString(objv[1]) << " not supported" << "\n";

  return TCL_OK;
}


int
fixedNodes(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  Domain * domain = (Domain*)clientData;

  SP_Constraint *theSP;
  SP_ConstraintIter &spIter = domain->getDomainAndLoadPatternSPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theSP = spIter()) != 0) {
    tag = theSP->getNodeTag();
    tags.insert(tag);
  }
  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());
  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
fixedDOFs(ClientData clientData, Tcl_Interp *interp, int argc, Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  Domain * theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - fixedDOFs fNode?\n";
    return TCL_ERROR;
  }

  int fNode;
  if (Tcl_GetIntFromObj(interp, objv[1], &fNode) != TCL_OK) {
    opserr << "WARNING fixedDOFs fNode? - could not read fNode? \n";
    return TCL_ERROR;
  }

  SP_Constraint *theSP;
  SP_ConstraintIter &spIter = theDomain->getDomainAndLoadPatternSPs();

  Node *node = theDomain->getNode(fNode);
  if (node == nullptr) {
    opserr << OpenSees::PromptValueError << " fixedDOFs fNode? - could not find node with tag " << fNode << "\n";
    return TCL_ERROR;
  }

  Vector fixed(node->getNumberDOF());
  while ((theSP = spIter()) != nullptr) {
    int tag = theSP->getNodeTag();
    if (tag == fNode) {
      fixed(theSP->getDOF_Number()) = 1;
    }
  }

  char buffer[20];
  for (int i = 0; i < 6; ++i) {
    if (fixed(i) == 1) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
constrainedNodes(ClientData clientData, Tcl_Interp *interp, int argc,
                 Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  Domain * theDomain = (Domain*)clientData;

  bool all = true;
  int rNode;
  if (argc > 1) {
    if (Tcl_GetIntFromObj(interp, objv[1], &rNode) != TCL_OK) {
      opserr << "WARNING constrainedNodes <rNode?> - could not read rNode? \n";
      return TCL_ERROR;
    }
    all = false;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = theDomain->getMPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  while ((theMP = mpIter()) != nullptr) {
    int tag = theMP->getNodeConstrained();
    if (all || rNode == theMP->getNodeRetained()) {
      tags.insert(tag);
    }
  }
  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());
  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
constrainedDOFs(ClientData clientData, Tcl_Interp *interp, int argc,
                Tcl_Obj *const *objv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - constrainedDOFs cNode? <rNode?> <rDOF?>\n";
    return TCL_ERROR;
  }

  int cNode;
  if (Tcl_GetIntFromObj(interp, objv[1], &cNode) != TCL_OK) {
    opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
              "read cNode?\n";
    return TCL_ERROR;
  }

  int rNode;
  bool allNodes = true;
  if (argc > 2) {
    if (Tcl_GetIntFromObj(interp, objv[2], &rNode) != TCL_OK) {
      opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
                "read rNode? \n";
      return TCL_ERROR;
    }
    allNodes = false;
  }

  int rDOF;
  bool allDOFs = true;
  if (argc > 3) {
    if (Tcl_GetIntFromObj(interp, objv[3], &rDOF) != TCL_OK) {
      opserr << "WARNING constrainedDOFs cNode? <rNode?> <rDOF?> - could not "
                "read rDOF? \n";
      return TCL_ERROR;
    }
    rDOF--;
    allDOFs = false;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = theDomain->getMPs();

  bool constrained[MAX_NDF];
  while ((theMP = mpIter()) != nullptr) {
    int tag = theMP->getNodeConstrained();
    if (tag == cNode) {
      if (allNodes || rNode == theMP->getNodeRetained()) {
        const ID &cDOFs = theMP->getConstrainedDOFs();
        int n = cDOFs.Size();
        if (allDOFs) {
          for (int i = 0; i < n; ++i)
            constrained[cDOFs(i)] = true;

        } else {
          const ID &rDOFs = theMP->getRetainedDOFs();
          for (int i = 0; i < n; ++i)
            if (rDOF == rDOFs(i))
              constrained[cDOFs(i)] = true;

        }
      }
    }
  }
  char buffer[20];
  for (int i = 0; i < MAX_NDF; ++i) {
    if (constrained[i]) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }
  return TCL_OK;
}

int
retainedDOFs(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  //
  // retainedDOFs rNode? <cNode?> <cDOF?>
  //
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "want - retainedDOFs rNode? <cNode?> <cDOF?>\n";
    return TCL_ERROR;
  }

  int rNode;
  if (Tcl_GetInt(interp, argv[1], &rNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "could not read rNode? \n";
    return TCL_ERROR;
  }

  int cNode;
  bool allNodes = 1;
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &cNode) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "could not read cNode? \n";
      return TCL_ERROR;
    }
    allNodes = 0;
  }

  int cDOF;
  bool allDOFs = 1;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &cDOF) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "could not read cDOF? \n";
      return TCL_ERROR;
    }
    cDOF--;
    allDOFs = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = domain->getMPs();

  int tag;
  int i;
  int n;
  Vector retained(6);
  while ((theMP = mpIter()) != nullptr) {
    tag = theMP->getNodeRetained();
    if (tag == rNode) {
      if (allNodes || cNode == theMP->getNodeConstrained()) {
        const ID &rDOFs = theMP->getRetainedDOFs();
        n = rDOFs.Size();
        if (allDOFs) {
          for (i = 0; i < n; ++i) {
            retained(rDOFs(i)) = 1;
          }
        } else {
          const ID &cDOFs = theMP->getConstrainedDOFs();
          for (int i = 0; i < n; ++i) {
            if (cDOF == cDOFs(i))
              retained(rDOFs(i)) = 1;
          }
        }
      }
    }
  }
  char buffer[20];
  for (int i = 0; i < 6; ++i) {
    if (retained(i) == 1) {
      sprintf(buffer, "%d ", i + 1);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
updateElementDomain(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv)
{
  // Need to "setDomain" to make the change take effect.
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  ElementIter &theElements = domain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != nullptr)
    theElement->setDomain(domain);

  return 0;
}
}
}
