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

// Written: Minjie

// Description: misc commands

#include <elementAPI.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <ElementalLoad.h>
#include <ElementalLoadIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <Parameter.h>
#include <RandomVariable.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <TimeSeries.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>
#include <NodalLoad.h>
#include <NodalLoadIter.h>
#include <Matrix.h>
#include <MeshRegion.h>
#include <StringContainer.h>
#include <fstream>
#include <string>
#include <InitialStateParameter.h>
#include <RigidRod.h>
#include <RigidBeam.h>
#include <RigidDiaphragm.h>
#include <vector>
#include <TriMesh.h>
#include <TetMesh.h>
#include <BackgroundMesh.h>
#include <Damping.h>

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#include <metis.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

int OPS_loadConst()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    theDomain->setLoadConstant();

    if (OPS_GetNumRemainingInputArgs() == 2) {
	const char* opt = OPS_GetString();
	if (strcmp(opt, "-time") == 0) {
	    double newTime;
	    int numdata = 1;
	    if (OPS_GetDoubleInput(&numdata, &newTime) < 0) {
		opserr<<"WARNING readingvalue - loadConst -time value\n";
		return -1;
	    }
	    theDomain->setCurrentTime(newTime);
	    theDomain->setCommittedTime(newTime);
	}
    }

    return 0;
}

int OPS_calculateNodalReactions()
{
    // make sure at least one other argument to contain type of system
    int incInertia = 0;

    if (OPS_GetNumRemainingInputArgs() == 1)  {
	const char* type = OPS_GetString();

	if ((strcmp(type,"-incInertia") == 0)
	    || (strcmp(type,"-dynamical") == 0)
	    || (strcmp(type,"-Dynamic") == 0)
	    || (strcmp(type,"-dynamic") == 0))

	    incInertia = 1;

	else if ((strcmp(type,"-rayleigh") == 0))

	    incInertia = 2;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    theDomain->calculateNodalReactions(incInertia);

    return 0;
}

int OPS_rayleighDamping() {
    if (OPS_GetNumRemainingInputArgs() < 4) {
        opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - "
                  "not enough arguments to command\n";
        return -1;
    }

    // double alphaM, betaK, betaK0, betaKc;
    double data[4];
    int numdata = 4;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
        opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - "
                  "could not read ? \n";
        return -1;
    }

    Domain *theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
        const char *option = OPS_GetString();
        bool isEle = true;
        if (strcmp(option, "-ele") == 0) {
            isEle = true;
        } else if (strcmp(option, "-node") == 0) {
            isEle = false;
        } else {
            opserr << "WARNING: valid options are -ele or -node\n";
            return -1;
        }
        numdata = OPS_GetNumRemainingInputArgs();
        std::vector<int> tags(numdata);
        if (OPS_GetIntInput(&numdata, &tags[0]) < 0) {
            opserr << "WARNING: failed to get element tags\n";
            return -1;
        }

        for (int i = 0; i < (int)tags.size(); ++i) {
            if (isEle) {
                Element *ele = theDomain->getElement(tags[i]);
                if (ele == 0) {
                    opserr << "WARNING: element " << tags[i]
                           << " does not exist\n";
                    return -1;
                }
                ele->setRayleighDampingFactors(data[0], data[1],
                                               data[2], data[3]);
            } else {
                Node *node = theDomain->getNode(tags[i]);
                if (node == 0) {
                    opserr << "WARNING: node " << tags[i]
                           << " does not exist\n";
                    return -1;
                }
                node->setRayleighDampingFactor(data[0]);
            }
        }

    } else {
        theDomain->setRayleighDampingFactors(data[0], data[1],
                                             data[2], data[3]);
    }

    return 0;
}

int OPS_setCreep()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING illegal command - setCreep value? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int newFlag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &newFlag) < 0) {
	opserr << "WARNING reading creep value - setCreep value? \n";
	return -1;
    } else {
	theDomain->setCreep(newFlag);
    }
    return 0;
}

int OPS_setTime()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING illegal command - time pseudoTime? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    double newTime;
    int numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &newTime) < 0) {
	opserr << "WARNING reading time value - time pseudoTime? \n";
	return -1;
    } else {
	theDomain->setCurrentTime(newTime);
	theDomain->setCommittedTime(newTime);
    }
    return 0;
}

int OPS_removeObject()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - remove objectType?\n";
	return -1;
    }
    const char* type = OPS_GetString();

    int tag;
    if ((strcmp(type,"element") == 0) || (strcmp(type,"ele") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove element eleTag?\n";
	    return -1;
	}

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove element tag? failed to read tag\n ";
	    return -1;
	}
	Element *theEle = theDomain->removeElement(tag);
	if (theEle != 0) {
	    // we also have to remove any elemental loads from the domain
	    LoadPatternIter &theLoadPatterns = theDomain->getLoadPatterns();
	    LoadPattern *thePattern;

	    // go through all load patterns
	    while ((thePattern = theLoadPatterns()) != 0) {
		ElementalLoadIter theEleLoads = thePattern->getElementalLoads();
		ElementalLoad *theLoad;

		// go through all elemental loads in the pattern
		while ((theLoad = theEleLoads()) != 0) {

		    // remove & destroy elemental from elemental load if there
		    // note - if last element in load, remove the load and delete it

		    /* *****************
		       int numLoadsLeft = theLoad->removeElement(tag);
		       if (numLoadsLeft == 0) {
		       thePattern->removeElementalLoad(theLoad->getTag());
		       delete theLoad;
		       }
		    *********************/
		}
	    }

	    // finally invoke the destructor on the element
	    delete theEle;
	}
    }

    else if (strcmp(type,"loadPattern") == 0 || strcmp(type,"pattern") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove loadPattern patternTag?\n";
	    return -1;
	}

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove loadPattern tag? failed to read tag\n ";
	    return -1;
	}
	LoadPattern *thePattern = theDomain->removeLoadPattern(tag);
	if (thePattern != 0) {
	    thePattern->clearAll();
	    delete thePattern;
	}
    }

    else if (strcmp(type,"parameter") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove parameter paramTag?\n";
	    return -1;
	}

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove parameter tag? failed to read tag\n ";
	    return -1;
	}
	Parameter *theParameter = theDomain->removeParameter(tag);
	if (theParameter != 0) {
	    delete theParameter;
	}
    }

    else if (strcmp(type,"node") == 0) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove node nodeTag?\n";
	    return -1;
	}
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove node tag? failed to read tag\n";
	    return -1;
	}
	Node *theNode = theDomain->removeNode(tag);
	if (theNode != 0) {
	    delete theNode;
	}
	Pressure_Constraint* thePC = theDomain->removePressure_Constraint(tag);
	if(thePC != 0) {
	    delete thePC;
	}
    }


    else if (strcmp(type,"recorders") == 0) {
	theDomain->removeRecorders();
    }

    else if ((strcmp(type,"recorder") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove recorder recorderTag?\n";
	    return -1;
	}
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove recorder tag? failed to read tag\n";
	    return -1;
	}
	return theDomain->removeRecorder(tag);
    }

    else if ((strcmp(type,"timeSeries") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove timeSeries $tag\n";
	    return -1;
	}
	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
	    opserr << "WARNING remove timeSeries tag? failed to read tag\n";
	    return -1;
	}
	return OPS_removeTimeSeries(tag);
    }


    else if ((strcmp(type,"SPconstraint") == 0) || (strcmp(type,"sp") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove SPconstraint spTag? -or- remove SPconstraint nodeTag? dofTag? <patternTag?>\n";
	    return -1;
	}
	if (OPS_GetNumRemainingInputArgs() == 1) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING remove sp tag? failed to read tag\n";
		return -1;
	    }

	    SP_Constraint *theSPconstraint = theDomain->removeSP_Constraint(tag);
	    if (theSPconstraint != 0) {
		delete theSPconstraint;
	    }
	} else {

	    // nodeTag, dofTag patternTag
	    int tags[3] = {0,0,-1};
	    int numdata = OPS_GetNumRemainingInputArgs();
	    if (numdata > 3) numdata = 3;
	    if (OPS_GetIntInput(&numdata, tags) < 0) {
		opserr << "WARNING remove sp tag? failed to read tags\n";
		return -1;
	    }

	    tags[1]--;  // one for C++ indexing of dof

	    theDomain->removeSP_Constraint(tags[0],tags[1],tags[2]);

	    return 0;
	}
    }

    else if ((strcmp(type,"MPconstraint") == 0) || (strcmp(type,"mp") == 0)) {
	if (OPS_GetNumRemainingInputArgs() < 1) {
	    opserr << "WARNING want - remove MPconstraint nNodeTag? -or- remove MPconstraint -tag mpTag\n";
	    return -1;
	}
	int nodTag = 0;
	if (OPS_GetNumRemainingInputArgs() == 1) {
	    int numdata = 1;
	    if (OPS_GetIntInput(&numdata, &nodTag) < 0) {
		opserr << "WARNING remove mp tag? failed to read tag\n";
		return -1;
	    }

	    theDomain->removeMP_Constraints(nodTag);
	    return 0;
	}
	if (OPS_GetNumRemainingInputArgs() > 1) {
	    const char* type = OPS_GetString();
	    if (strcmp(type,"-tag") == 0) {
		int numdata = 1;
		if (OPS_GetIntInput(&numdata, &nodTag) < 0) {
		    opserr << "WARNING remove mp -tag mpTag? failed to read mpTag\n";
		    return -1;
		}
	    }

	    theDomain->removeMP_Constraint(nodTag);
	    return 0;
	}
    }


    else
	opserr << "WARNING remove " << type << " not supported\n";

    return 0;
}

int OPS_addNodalMass()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - mass nodeTag? <mass values>?\n";
	return -1;
    }

    int nodeTag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &nodeTag) < 0) {
	opserr << "WARNING invalid nodeTag\n";
	return -1;
    }

    int ndf = OPS_GetNDF();
    Matrix mass(ndf, ndf);
    double theMass;
    for (int i=0; i<ndf; i++) {

	if (OPS_GetNumRemainingInputArgs() < 1) {
	    break;
	}

	if (OPS_GetDoubleInput(&numdata, &theMass) < 0) {
	    opserr << "WARNING invalid mass value\n";
	    return -1;
	}
	mass(i,i) = theMass;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (theDomain->setMass(mass, nodeTag) != 0) {
	opserr << "WARNING failed to set mass at node " << nodeTag << "\n";
	return -1;
    }

    return 0;

}

int OPS_buildModel()
{
    return 0;
}

int OPS_setNodeDisp()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - setNodeDisp nodeTag? dof? value? <-commit>\n";
        return -1;
    }

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING setNodeDisp nodeTag? dof? - could not read nodeTag? \n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
        opserr << "WARNING setNodeDisp -- node with tag " << tag << " not found" << endln;
        return -1;
    }

    if (OPS_GetIntInput(&numdata, &dof) < 0) {
        opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
        opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
        const char* optArg = OPS_GetString();
        if (strcmp(optArg, "-commit") == 0)
            commit = true;
    }

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
        Vector disp(numDOF);
        disp = theNode->getDisp();
        disp(dof) = value;
        theNode->setTrialDisp(disp);
    }
    if (commit)
        theNode->commitState();

    return 0;
}

int OPS_setNodeVel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - setNodeVel nodeTag? dof? value? <-commit>\n";
        return -1;
    }

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING setNodeVel nodeTag? dof? - could not read nodeTag? \n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
        opserr << "WARNING setNodeVel -- node with tag " << tag << " not found" << endln;
        return -1;
    }

    if (OPS_GetIntInput(&numdata, &dof) < 0) {
        opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
        opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
        const char* optArg = OPS_GetString();
        if (strcmp(optArg, "-commit") == 0)
            commit = true;
    }

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
        Vector vel(numDOF);
        vel = theNode->getVel();
        vel(dof) = value;
        theNode->setTrialVel(vel);
    }
    if (commit)
        theNode->commitState();

    return 0;
}

int OPS_setNodeAccel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
        opserr << "WARNING want - setNodeAccel nodeTag? dof? value? <-commit>\n";
        return -1;
    }

    int tag;
    int dof = -1;
    double value = 0.0;
    bool commit = false;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING setNodeAccel nodeTag? dof? - could not read nodeTag? \n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
        opserr << "WARNING setNodeAccel -- node with tag " << tag << " not found" << endln;
        return -1;
    }

    if (OPS_GetIntInput(&numdata, &dof) < 0) {
        opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
        opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read dof? \n";
        return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
        const char* optArg = OPS_GetString();
        if (strcmp(optArg, "-commit") == 0)
            commit = true;
    }

    dof--;

    int numDOF = theNode->getNumberDOF();

    if (dof >= 0 && dof < numDOF) {
        Vector accel(numDOF);
        accel = theNode->getAccel();
        accel(dof) = value;
        theNode->setTrialAccel(accel);
    }
    if (commit)
        theNode->commitState();

    return 0;
}

int OPS_setElementRayleighDampingFactors()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
	opserr << "WARNING setElementRayleighDampingFactors eleTag? alphaM? betaK? betaK0? betaKc? - not enough arguments to command\n";
	return -1;
    }
    int eleTag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
	opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read eleTag? \n";
	return -1;
    }

    numdata = 4;
    double data[4];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING rayleigh alphaM? betaK? betaK0? betaKc? - could not read double inputs? \n";
	return -1;
    }

    double alphaM = data[0];
    double betaK = data[1];
    double betaK0 = data[2];
    double betaKc = data[3];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theEle = theDomain->getElement(eleTag);
    theEle->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);

    return 0;
}

int OPS_MeshRegion()
{
    int tag;
    double alphaM = 0.0;
    double betaK  = 0.0;
    double betaK0 = 0.0;
    double betaKc = 0.0;
    Damping *theDamping = 0;

    ID *theNodes = 0;
    ID *theElements = 0;
    int numNodes = 0;
    int numElements = 0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // first get tag for region
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING region tag? - no tag specified\n";
	return -1;
    }

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING region tag? .. - invalid tag " << endln;
	return -1;
    }

    // now continue until end of command
    bool only = false;
    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* flag = OPS_GetString();

	if (strcmp(flag,"-ele") == 0 || strcmp(flag,"-eleOnly") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING region tag? .. -ele tag1? .. - no ele tags specified\n";
		return -1;
	    }

	    //
	    // read in a list of ele until end of command or other flag
	    //
	    if (theElements == 0) {
		theElements = new ID(0, 64);
	    }
	    int eleTag;
	    while (OPS_GetNumRemainingInputArgs() > 0) {
		if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
		    // back one arg
		    OPS_ResetCurrentInputArg(-1);
		    break;
		}

		(*theElements)[numElements++] = eleTag;
	    }

	    if (strcmp(flag,"-eleOnly") == 0) {
		only = true;
	    }

	} else if (strcmp(flag,"-eleRange") == 0 || strcmp(flag,"-eleOnlyRange") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "WARNING region tag? .. -eleRange start? end?  .. - no ele tags specified\n";
		return -1;
	    }

	    //
	    // read in start and end tags of two elements & add set [start,end]
	    //

	    int start, end;
	    if (OPS_GetIntInput(&numdata, &start) < 0) {
		opserr << "WARNING region tag? -eleRange start? end? - invalid start " << endln;
		return -1;
	    }
	    if (OPS_GetIntInput(&numdata, &end) < 0) {
		opserr << "WARNING region tag? -eleRange start? end? - invalid end " << endln;
		return -1;
	    }
	    if (start > end) {
		int swap = end;
		end = start;
		start = swap;
	    }
	    int numEle = end-start+1;

	    if (theElements == 0) {
		theElements = new ID(0, numEle);
	    }
	    for (int i=start; i<=end; i++) {
		(*theElements)[numElements++] = i;
	    }

	    if (strcmp(flag,"-eleOnlyRange") == 0) {
		only = true;
	    }

	} else if (strcmp(flag,"-node") == 0 || strcmp(flag,"-nodeOnly") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING region tag? .. -node tag1? .. - no node tags specified\n";
		return -1;
	    }

	    // read in list of nodes
	    if (theNodes == 0) {
		theNodes = new ID(0, 64);
	    }

	    int nodTag;
	    while (OPS_GetNumRemainingInputArgs() > 0) {
		if (OPS_GetIntInput(&numdata, &nodTag) < 0) {
		    // back one arg
		    OPS_ResetCurrentInputArg(-1);
		    break;
		}

		(*theNodes)[numNodes++] = nodTag;
	    }

	    if (strcmp(flag,"-nodeOnly") == 0) {
		only = true;
	    }

	} else if (strcmp(flag,"-nodeRange") == 0 || strcmp(flag,"-nodeOnlyRange") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 2) {
		opserr << "WARNING region tag? .. -nodeRange start? end?  .. - no node tags specified\n";
		return -1;
	    }

	    // read in start and end ele tags
	    int start, end;
	    if (OPS_GetIntInput(&numdata, &start) < 0) {
		opserr << "WARNING region tag? -nodeRange start? end? - invalid start " << endln;
		return -1;
	    }
	    if (OPS_GetIntInput(&numdata, &end) < 0) {
		opserr << "WARNING region tag? -nodeRange start? end? - invalid end " << endln;
		return -1;
	    }
	    if (start > end) {
		int swap = end;
		end = start;
		start = swap;
	    }
	    int numNode = end-start+1;

	    if (theNodes == 0) {
		theNodes = new ID(0, numNode);
	    }
	    for (int i=start; i<=end; i++) {
		(*theNodes)[numNodes++] = i;
	    }

	    if (strcmp(flag,"-nodeOnlyRange") == 0) {
		only = true;
	    }

	} else if (strcmp(flag,"-rayleigh") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 4) {
		opserr << "WARNING region tag? .. -rayleigh aM? bK? bK0?  .. - not enough factors\n";
		return -1;
	    }

	    // read in rayleigh damping factors
	    if (OPS_GetDoubleInput(&numdata, &alphaM) < 0) {
		opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid aM " << endln;
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &betaK) < 0) {
		opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bK " << endln;
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &betaK0) < 0) {
		opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bK0 " << endln;
		return -1;
	    }
	    if (OPS_GetDoubleInput(&numdata, &betaKc) < 0) {
		opserr << "WARNING region tag? .. -rayleigh aM bK bK0 - invalid bKc " << endln;
		return -1;
	    }
	}  else if (strcmp(flag,"-damp") == 0) {

	    // ensure no segmentation fault if user messes up
	    if (OPS_GetNumRemainingInputArgs() < 1) {
		opserr << "WARNING region tag? .. -damp dampingTag?\n";
		return -1;
	    }
		
		int dampingTag = 0;

	    // read in damping factors
	    if (OPS_GetIntInput(&numdata, &dampingTag) < 0) {
		opserr << "WARNING region tag? .. -damp dampingTag?" << endln;
		return -1;
	    }
	    theDamping = OPS_getDamping(dampingTag);
	    if(theDamping == 0) {
		opserr << "damping not found\n";
	    return -1;
		}
	}
    }

    MeshRegion *theRegion = theDomain->getRegion(tag);
    if (theRegion == 0) {
	theRegion = new MeshRegion(tag);
	if (theRegion == 0) {
	    opserr << "could not create region\n";
	    return -1;
	}
	if (theDomain->addRegion(*theRegion) < 0) {
	    opserr << "WARNING could not add to domain - region " << tag << endln;
	    delete theRegion;
	    return -1;
	}
    }

    // if elements or nodes have been set, set them in the Region
    if (theElements != 0) {
	if (only) {
	    theRegion->setElementsOnly(*theElements);
	} else {
	    theRegion->setElements(*theElements);
	}
    }

    if (theNodes != 0) {
	if (theElements == 0) {
	    if (only) {
		theRegion->setNodesOnly(*theNodes);
	    } else {
		theRegion->setNodes(*theNodes);
	    }

	} else {
	    opserr << "WARNING region - both elements & nodes set, ONLY set using elements\n";
	}
    }

    // if damping has been specified set the damping factors
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0) {
	theRegion->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    }
	
	if(theDamping) theRegion->setDamping(theDamping);

    if (theElements != 0) {
	delete theElements;
    }

    if (theNodes != 0) {
	delete theNodes;
    }

    return 0;
}

extern
int peerSearchNGA(const char *eq,
		  const char *soilType,
		  const char *fault,
		  const char *magLo,
		  const char *magHi,
		  const char *distLo,
		  const char *distHi,
		  const char *vsLo,
		  const char *vsHi,
		  const char *pgaLo,
		  const char *pgaHi,
		  const char *latSW,
		  const char *latNE,
		  const char *lngSW,
		  const char *lngNW,
		  StringContainer &recordNames);

int OPS_peerNGA()
{
    StringContainer ngaRecordNames;
    const char *eq =0;
    const char *soilType = 0;
    const char *fault =0;
    const char *magLo =0;
    const char *magHi =0;
    const char *distLo =0;
    const char *distHi =0;
    const char *vsLo =0;
    const char *vsHi =0;
    const char *pgaLo =0;
    const char *pgaHi =0;
    const char *latSW =0;
    const char *latNE =0;
    const char *lngSW =0;
    const char *lngNW =0;

    while (OPS_GetNumRemainingInputArgs() > 1) {

	const char* flag = OPS_GetString();

	if (strcmp(flag,"-eq") == 0) {
	    eq = OPS_GetString();
	} else if (strcmp(flag,"-fault") == 0) {
	    fault = OPS_GetString();
	} else if (strcmp(flag,"-soil") == 0) {
	    soilType = OPS_GetString();
	} else if (strcmp(flag,"-magLo") == 0) {
	    magLo = OPS_GetString();
	} else if (strcmp(flag,"-magHi") == 0) {
	    magHi = OPS_GetString();
	} else if (strcmp(flag,"-distLo") == 0) {
	    distLo  = OPS_GetString();
	} else if (strcmp(flag,"-distHi") == 0) {
	    distHi = OPS_GetString();
	} else if (strcmp(flag,"-vsLo") == 0) {
	    vsLo = OPS_GetString();
	} else if (strcmp(flag,"-vsHi") == 0) {
	    vsHi = OPS_GetString();
	} else if (strcmp(flag,"-pgaLo") == 0) {
	    pgaLo = OPS_GetString();
	} else if (strcmp(flag,"-pgaHi") == 0) {
	    pgaHi = OPS_GetString();
	} else if (strcmp(flag,"-latSW") == 0) {
	    latSW = OPS_GetString();
	} else if (strcmp(flag,"-latNE") == 0) {
	    latNE = OPS_GetString();
	} else if (strcmp(flag,"-lngSW") == 0) {
	    lngSW = OPS_GetString();
	} else if (strcmp(flag,"-lngNW") == 0) {
	    lngNW = OPS_GetString();
	}
    }

    peerSearchNGA(eq,
		  soilType,
		  fault,
		  magLo,
		  magHi,
		  distLo,
		  distHi,
		  vsLo,
		  vsHi,
		  pgaLo,
		  pgaHi,
		  latSW,
		  latNE,
		  lngSW,
		  lngNW,
		  ngaRecordNames);

    int numStrings = ngaRecordNames.getNumStrings();
    if (numStrings == 0) return 0;

    int len = 3;
    for (int i=0; i<numStrings; i++) {
	len += (int)strlen(ngaRecordNames.getString(i))+3;
    }

    char* result = new char[len];
    strcpy(result, " ");
    for (int i=0; i<numStrings; i++) {
	strcat(result, ngaRecordNames.getString(i));
	strcat(result, " ");
    }

    if (OPS_SetString(result) < 0) {
	opserr << "WARNING failed to set result string\n";
	delete [] result;
	return -1;
    }

    delete [] result;
    return 0;
}

int OPS_domainChange()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    theDomain->domainChange();

    return 0;
}

int OPS_record()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    theDomain->record(false);

    return 0;
}

int OPS_stripOpenSeesXML()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "ERROR incorrect # args - stripXML input.xml output.dat <output.xml>\n";
	return -1;
    }

    const char *inputFile = OPS_GetString();
    const char *outputDataFile = OPS_GetString();
    const char *outputDescriptiveFile = 0;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	outputDescriptiveFile = OPS_GetString();
    }

    // open files
    std::ifstream theInputFile;
    theInputFile.open(inputFile, std::ios::in);
    if (theInputFile.bad()) {
	opserr << "stripXML - error opening input file: " << inputFile << endln;
	return -1;
    }

    std::ofstream theOutputDataFile;
    theOutputDataFile.open(outputDataFile, std::ios::out);
    if (theOutputDataFile.bad()) {
	opserr << "stripXML - error opening input file: " << outputDataFile << endln;
	return -1;
    }

    std::ofstream theOutputDescriptiveFile;
    if (outputDescriptiveFile != 0) {
	theOutputDescriptiveFile.open(outputDescriptiveFile, std::ios::out);
	if (theOutputDescriptiveFile.bad()) {
	    opserr << "stripXML - error opening input file: " << outputDescriptiveFile << endln;
	    return -1;
	}
    }

    std::string line;
    bool spitData = false;
    while (! theInputFile.eof() ) {
	getline(theInputFile, line);
	const char *inputLine = line.c_str();

	if (spitData == true) {
	    if (strstr(inputLine,"</Data>") != 0)
		spitData = false;
        else {
            //theOutputDataFile << line << endln;
        }
	} else {
	    const char *inputLine = line.c_str();
	    if (strstr(inputLine,"<Data>") != 0)
		spitData = true;
        else if (outputDescriptiveFile != 0) {
            //theOutputDescriptiveFile << line << endln;
        }
	}
    }

    theInputFile.close();
    theOutputDataFile.close();

    if (outputDescriptiveFile != 0)
	theOutputDescriptiveFile.close();

    return 0;
}

extern int binaryToText(const char *inputFilename, const char *outputFilename);
extern int textToBinary(const char *inputFilename, const char *outputFilename);

int OPS_convertBinaryToText()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "ERROR incorrect # args - convertBinaryToText inputFile outputFile\n";
	return -1;
    }

    const char *inputFile = OPS_GetString();
    const char *outputFile = OPS_GetString();

    return binaryToText(inputFile, outputFile);
}

int OPS_convertTextToBinary()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "ERROR incorrect # args - convertTextToBinary inputFile outputFile\n";
	return -1;
    }

    const char *inputFile = OPS_GetString();
    const char *outputFile = OPS_GetString();

    return textToBinary(inputFile, outputFile);
}

int OPS_InitialStateAnalysis()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING: Incorrect number of arguments for InitialStateAnalysis command" << endln;
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const char* flag = OPS_GetString();
    if (strcmp(flag,"on") == 0) {
	opserr << "InitialStateAnalysis ON" << endln;

	// set global variable to true
	// FMK changes for parallel:
	// ops_InitialStateAnalysis = true;

	Parameter *theP = new InitialStateParameter(true);
	theDomain->addParameter(theP);
	delete theP;

	return 0;

    } else if (strcmp(flag,"off") == 0) {
	opserr << "InitialStateAnalysis OFF" <<endln;

	// call revert to start to zero the displacements
	theDomain->revertToStart();

	// set global variable to false
	// FMK changes for parallel
	// ops_InitialStateAnalysis = false;
	Parameter *theP = new InitialStateParameter(false);
	theDomain->addParameter(theP);
	delete theP;

	return 0;

    } else {
	opserr << "WARNING: Incorrect arguments - want InitialStateAnalysis on, or InitialStateAnalysis off" << endln;

	return -1;
    }

    return 0;
}

int OPS_maxOpenFiles()
{
    int maxOpenFiles;

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want maxNumFiles\n";
	return -1;
    }
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &maxOpenFiles) < 0) {
	return -1;
    }

#ifdef _WIN32
    int newMax = _setmaxstdio(maxOpenFiles);
    if (maxOpenFiles > 2048) {
	opserr << "setMaxOpenFiles: too many files specified (2048 max)\n";
    } else {
	if (newMax != maxOpenFiles) {
	    opserr << "setMaxOpenFiles FAILED: max allowed files: " << newMax;
	    return -1;
	}
    }
    return 0;
#endif

    opserr << "setMaxOpenFiles FAILED: - command not available on this machine\n";
    return 0;
}

int OPS_RigidLink()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING rigidLink linkType? rNode? cNode?\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const char* type = OPS_GetString();

    int rNode, cNode;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &rNode) < 0) {
	opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read rNode \n";
	return -1;
    }
    if (OPS_GetIntInput(&numdata, &cNode) < 0) {
	opserr << "WARNING rigidLink linkType? rNode? cNode? - could not read CNode \n";
	return -1;
    }

    // construct a rigid rod or beam depending on 1st arg
    if ((strcmp(type,"-bar") == 0) || (strcmp(type,"bar") == 0)) {
	RigidRod theLink(*theDomain, rNode, cNode);
    } else if ((strcmp(type,"-beam") == 0) || (strcmp(type,"beam") == 0)) {
	RigidBeam theLink(*theDomain, rNode, cNode);
    } else {
	opserr << "WARNING rigidLink linkType? rNode? cNode? - unrecognised link type (-bar, -beam) \n";
	return -1;
    }

    return 0;
}

int OPS_RigidDiaphragm()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING rigidDiaphragm perpDirn? rNode? <cNodes?>\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int rNode, perpDirn;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &perpDirn) < 0) {
	opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read perpDirn? \n";
	return -1;
    }

    if (OPS_GetIntInput(&numdata, &rNode) < 0) {
	opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read rNode \n";
	return -1;
    }

    // read in the constrained Nodes
    int numConstrainedNodes = OPS_GetNumRemainingInputArgs();
    ID constrainedNodes(numConstrainedNodes);
    for (int i=0; i<numConstrainedNodes; i++) {
	int cNode;
	if (OPS_GetIntInput(&numdata, &cNode) < 0) {
	    opserr << "WARNING rigidLink perpDirn rNode cNodes - could not read a cNode\n";
	    return -1;
	}
	constrainedNodes(i) = cNode;
    }

    RigidDiaphragm theLink(*theDomain, rNode, constrainedNodes, perpDirn-1);


    return 0;
}

int OPS_addElementRayleigh()
{
    // make sure correct number of arguments on command line
    if (OPS_GetNumRemainingInputArgs() < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: setElementRayleighFactors elementTag?  alphaM? $betaK? $betaKinit? $betaKcomm? \n";
	return -1;
    }

    int eleTag =0;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
	opserr << "WARNING: setElementRayleighFactors invalid eleTag: ";
	opserr << " \n";
	return -1;
    }

    double alphaM,betaK,betaKinit,betaKcomm;

    if (OPS_GetDoubleInput(&numdata, &alphaM) < 0) {
	opserr << "WARNING : setElementRayleighFactors invalid ";
	opserr << "alphaM: " << endln;
	return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &betaK) < 0) {
	opserr << "WARNING : setElementRayleighFactors invalid ";
	opserr << "betaK: " << endln;
	return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &betaKinit) < 0) {
	opserr << "WARNING : setElementRayleighFactors invalid ";
	opserr << "betaKinit: " << endln;
	return -1;
    }

    if (OPS_GetDoubleInput(&numdata, &betaKcomm) < 0) {
	opserr << "WARNING : setElementRayleighFactors invalid ";
	opserr << "betaKcomm: " << endln;
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Element* elePtr = theDomain->getElement(eleTag);

    if (elePtr == 0) {
	opserr << "WARNING : setElementRayleighFactors invalid eleTag: " << eleTag << " the element does not exist in the domain \n";
	return -1;
    }


    if ( elePtr->setRayleighDampingFactors(alphaM, betaK, betaKinit, betaKcomm) != 0 ) {
	opserr << "ERROR : setElementRayleighFactors: FAILED to add damping factors for element " << eleTag << "\n";
	return -1;
    }

    return 0;
}

extern int OPS_LineMesh();
extern int OPS_TriMesh();
extern int OPS_TetMesh();
extern int OPS_QuadMesh();
extern int OPS_BgMesh();
extern int OPS_ParticleGroup();
extern int OPS_Flume();
extern int OPS_BeamBrick();
extern int OPS_BeamDisk();
extern BackgroundMesh& OPS_getBgMesh();


int OPS_mesh()
{
    // make sure correct number of arguments on command line
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: mesh type? ...>\n";
        return -1;
    }

    // mesh type
    const char* type = OPS_GetString();
    int res = 0;
    if (strcmp(type, "line") == 0) {
        res = OPS_LineMesh();
    } else if (strcmp(type, "tri") == 0) {
        res = OPS_TriMesh();
    } else if (strcmp(type, "part") == 0) {
	res = OPS_ParticleGroup();
    } else if (strcmp(type, "bg") == 0) {
        res = OPS_BgMesh();
    } else if (strcmp(type, "tet") == 0) {
	res = OPS_TetMesh();
    } else if (strcmp(type, "quad") == 0) {
	res = OPS_QuadMesh();
    } else if (strcmp(type, "flume") == 0) {
	res = OPS_Flume();
    } else if (strcmp(type, "beamBrick") == 0) {
	res = OPS_BeamBrick();
    } else if (strcmp(type, "beamDisk") == 0) {
	res = OPS_BeamDisk();
    } else {
        opserr<<"WARNING: mesh type "<<type<<" is unknown\n";
        return -1;
    }

    if (res < 0) {
        return -1;
    }

    return 0;
}

int OPS_remesh()
{
    if (OPS_GetNumRemainingInputArgs() > 0) {
	double alpha = -1.0;
	int numdata = 1;
	if (OPS_GetDoubleInput(&numdata,&alpha) < 0) {
	    opserr << "WARNING: invalid alpha\n";
	    return -1;
	}

	int ndm = OPS_GetNDM();

	if (ndm == 2) {
	    if (TriMesh::remesh(alpha) < 0) {
		opserr << "WARNING: failed to remesh\n";
		return -1;
	    }
	} else if (ndm == 3) {
	    if (TetMesh::remesh(alpha) < 0) {
		opserr << "WARNING: failed to remesh\n";
		return -1;
	    }
	}

    } else {
	BackgroundMesh& bgmesh = OPS_getBgMesh();
	if (bgmesh.remesh() < 0) {
	    opserr << "WARNING: failed to remesh background\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_getPID()
{
    int pid = 0;

#ifdef _PARALLEL_INTERPRETERS
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
#endif

    int size = 1;
    if (OPS_SetIntOutput(&size, &pid, true) < 0) {
	opserr << "WARNING: failed to set pid\n";
	return -1;
    }

    return 0;
}

int OPS_getNP()
{
    int nump = 1;

#ifdef _PARALLEL_INTERPRETERS
    MPI_Comm_size(MPI_COMM_WORLD, &nump);
#endif

    int size = 1;
    if (OPS_SetIntOutput(&size, &nump, true) < 0) {
	opserr << "WARNING: failed to set np\n";
	return -1;
    }

    return 0;
}

int OPS_barrier()
{
#ifdef _PARALLEL_INTERPRETERS
    return MPI_Barrier(MPI_COMM_WORLD);
#endif

    return 0;
}

int OPS_send()
{
#ifdef _PARALLEL_INTERPRETERS

    // get ids
    int otherPID = -1;
    int myPID = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPID);

    int np = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (OPS_GetNumRemainingInputArgs()  < 1) {
    	opserr << "WARNING: need send <'-pid' pid> data\n";
    	return -1;
    }

    // get other pid
    const char *sdata = OPS_GetString();
    if (strcmp(sdata, "-pid") == 0) {
        if (OPS_GetNumRemainingInputArgs() > 0) {
            int num = 1;
            if (OPS_GetIntInput(&num, &otherPID) < 0) {
                opserr << "WARNING: failed to get pid\n";
                return -1;
            }
            if (otherPID<0 || otherPID>=np || otherPID==myPID) {
                opserr << "WARNING: invalid pid "<<otherPID<<"\n";
                return -1;
            }
        } else {
            opserr << "WARNING: need pid\n";
            return -1;
        }
    } else {
        opserr << "WARNING: -pid is required\n";
        return -1;
    }

    // get data type
    int num = OPS_GetNumRemainingInputArgs();
    if (num < 1) {
        opserr << "WARNING: need data\n";
        return -1;
    }

    std::vector<double> ddata(num);
    MPI_Datatype datatype = MPI_DOUBLE;
    if (OPS_GetDoubleInput(&num, &ddata[0]) < 0) {
        if (num > 1) {
            opserr << "WARNING: data is string and size must be 1\n";
            return -1;
        }
        datatype = MPI_CHAR;
        OPS_ResetCurrentInputArg(-1);
        sdata = OPS_GetString();
    }

    // data length, 1-double, 2-string
    int msgLength[2] = {0,0};
    void* buffer = 0;
    if (datatype == MPI_DOUBLE) {
        msgLength[0] = num;
        msgLength[1] = 1;
        buffer = (void *) &ddata[0];
    } else {
        msgLength[0] = strlen(sdata)+1;
        msgLength[1] = 2;
        buffer = (void *) sdata;
    }

    // send data to PID
    MPI_Send((void *)(&msgLength[0]), 2, MPI_INT, otherPID,
             0, MPI_COMM_WORLD);
    MPI_Send(buffer, msgLength[0], datatype, otherPID,
             1, MPI_COMM_WORLD);

#endif
    return 0;
}

int OPS_recv()
{
#ifdef _PARALLEL_INTERPRETERS

    // get ids
    int otherPID = -1;
    int myPID = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPID);

    int np = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (OPS_GetNumRemainingInputArgs()  < 1) {
        opserr << "WARNING: need recv '-pid' pid\n";
        return -1;
    }

    // get other pid
    const char *sdata = OPS_GetString();
    if (strcmp(sdata, "-pid") == 0) {
        if (OPS_GetNumRemainingInputArgs() < 1) {
            opserr << "WANRING: need pid\n";
            return -1;
        }
        int num = 1;
        if (OPS_GetIntInput(&num, &otherPID) == 0) {
            if (otherPID<0 || otherPID>=np || otherPID==myPID) {
                opserr << "WARNING: invalid pid "<<otherPID<<"\n";
                return -1;
            }
        } else {
            OPS_ResetCurrentInputArg(-1);

            sdata = OPS_GetString();
            if (strcmp(sdata, "ANY") == 0) {
                opserr << "WARNING: ANY pid is no longer available\n";
                return -1;
            } else {
                opserr << "WARNING: Invalid pid\n";
                return -1;
            }
        }
    } else {
        opserr << "WARNING: -pid is required\n";
        return -1;
    }

    // receive data
    MPI_Status status;

    // receive length and type
    int msgLength[2] = {0, 0};
    MPI_Recv((void *)(&msgLength[0]), 2, MPI_INT, otherPID, 0,
             MPI_COMM_WORLD, &status);

    // get data
    if (msgLength[0] > 0) {

        // get type
        MPI_Datatype datatype = MPI_DOUBLE;
        char* gMsg = new char[msgLength[0]];
        std::vector<double> ddata(msgLength[0]);
        void* buffer = 0;
        if (msgLength[1] == 1) {
            datatype = MPI_DOUBLE;
            buffer = (void *) &ddata[0];
        } else {
            datatype = MPI_CHAR;
            buffer = (void *) gMsg;
        }

        // receive data
        MPI_Recv(buffer, msgLength[0], datatype, otherPID, 1,
                 MPI_COMM_WORLD, &status);

        // set outputs
        int res = 0;
        if (datatype == MPI_DOUBLE) {

            res = OPS_SetDoubleOutput(&msgLength[0], &ddata[0], false);

	    } else {

            res = OPS_SetString(gMsg);

	    }
	    if (res < 0) {
            opserr << "WARNING: failed to set results\n";
            return -1;
	    }
        delete [] gMsg;
	}


#endif

    return 0;
}

int OPS_Bcast() {

#ifdef _PARALLEL_INTERPRETERS
    // get ids
    int myPID = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &myPID);

    int np = 0;
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    int msgLength[2] = {0,0};
    void* buffer = 0;
    MPI_Datatype datatype;
    if (myPID == 0) {

        // send
        int num = OPS_GetNumRemainingInputArgs();
        if (num < 1) {
            opserr << "WARNING: need data\n";
            return -1;
        }

        std::vector<double> ddata(num);
        datatype = MPI_DOUBLE;
        const char* sdata = 0;
        if (OPS_GetDoubleInput(&num, &ddata[0]) < 0) {
            if (num > 1) {
                opserr
                    << "WARNING: data is string and size must be 1\n";
                return -1;
            }
            datatype = MPI_CHAR;
            OPS_ResetCurrentInputArg(-1);
            sdata = OPS_GetString();
        }

        // data length, 1-double, 2-string
        if (datatype == MPI_DOUBLE) {

            msgLength[0] = num;
            msgLength[1] = 1;
            buffer = (void *) &ddata[0];

        } else {

            msgLength[0] = strlen(sdata)+1;
            msgLength[1] = 2;
            buffer = (void *) sdata;
        }

        MPI_Bcast((void *)(&msgLength[0]), 2, MPI_INT,  0, MPI_COMM_WORLD);
        MPI_Bcast(buffer, msgLength[0], datatype, 0, MPI_COMM_WORLD);
    } else {

        // receive
        MPI_Bcast((void *)(&msgLength[0]), 2, MPI_INT, 0, MPI_COMM_WORLD);
        if (msgLength[0] > 0) {

            // get type
            datatype = MPI_DOUBLE;
            char* gMsg = new char[msgLength[0]];
            std::vector<double> ddata(msgLength[0]);
            if (msgLength[1] == 1) {

                datatype = MPI_DOUBLE;
                buffer = (void *) &ddata[0];

            } else {

                datatype = MPI_CHAR;
                buffer = (void *) gMsg;
            }
            MPI_Bcast(buffer, msgLength[0], datatype, 0, MPI_COMM_WORLD);

            // set outputs
            int res = 0;
            if (datatype == MPI_DOUBLE) {

                res = OPS_SetDoubleOutput(&msgLength[0], &ddata[0], false);

            } else {

                res = OPS_SetString(gMsg);

            }
            if (res < 0) {
                opserr << "WARNING: failed to set results\n";
                return -1;
            }
            delete [] gMsg;
        }
    }
#endif
    return 0;
}

int OPS_getNumThreads()
{
#ifdef _OPENMP
    int num = omp_get_max_threads();
    int numdata = 1;

    if (OPS_SetIntOutput(&numdata,&num,true) < 0) {
	opserr << "WARNING: failed to set output -- getNumThreads\n";
	return -1;
    }
#endif

    return 0;
}

int OPS_setNumThreads()
{
#ifdef _OPENMP
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING: need setNumThreads num\n";
	return -1;
    }
    
    int num;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata,&num) < 0) {
	opserr << "WARNING: failed to set output -- getNumThreads\n";
	return -1;
    }

    omp_set_num_threads(num);
#endif

    return 0;
}

int OPS_setStartNodeTag() {
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: needs tag\n";
        return -1;
    }

    int tag;
    int num = 1;
    if (OPS_GetIntInput(&num, &tag) < 0) {
        opserr << "WARNING: failed to get tag\n";
        return -1;
    }

    Mesh::setStartNodeTag(tag);

    return 0;
}

int OPS_partition() {
#ifdef _PARALLEL_INTERPRETERS
    // domain
    Domain* domain = OPS_GetDomain();
    if (domain ==0) {
        return 0;
    }

    int pid = 0;
    int np = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &np);

    if (np == 1) return 0;

    // parameters
    int ncuts = 1;
    int niter = 10;
    int ufactor = 30;
    int info = 0;
    while (OPS_GetNumRemainingInputArgs() > 0) {
        int num = 1;
        auto opt = OPS_GetString();
        if (strcmp(opt, "-ncuts") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0 &&
                OPS_GetIntInput(&num, &ncuts) < 0) {
                opserr << "WARNING: failed to get ncuts\n";
                return -1;
            }
            if (ncuts < 1) {
                ncuts = 1;
            }
        } else if (strcmp(opt, "-niter") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0 &&
                OPS_GetIntInput(&num, &niter) < 0) {
                opserr << "WARNING: failed to get niter\n";
                return -1;
            }
            if (niter < 1) {
                niter = 10;
            }
        } else if (strcmp(opt, "-ufactor") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0 &&
                OPS_GetIntInput(&num, &ncuts) < 0) {
                opserr << "WARNING: failed to get ufactor\n";
                return -1;
            }
            if (ufactor < 1) {
                ufactor = 30;
            }
        } else if (strcmp(opt, "-info") == 0) {
            info = METIS_DBG_INFO;
        }
    }


    // map of node tag to index
    std::map<int, idx_t> nind;
    Node *nd = 0;
    auto &nodes = domain->getNodes();
    int index = 0;
    while ((nd = nodes()) != 0) {
        nind[nd->getTag()] = index++;
    }

    // check size
    if (nind.empty()) {
        opserr << "WARNING: no node is found on processorr "<<pid<<"\n";
        return -1;
    }

    // number of nodes
    idx_t nn = (idx_t)nind.size();

    // check if all processors have same number
    idx_t nn_max = 0;
    if (MPI_Reduce(&nn, &nn_max, 1, MPI_INT,
                   MPI_MAX, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        opserr << "WARNING: failed to get max nn\n";
        return -1;
    }
    if (pid == 0) {
        if (nn_max != nn) {
            opserr << "WARNING: all processors should "
                      "have the same mesh before partitioning";
            return -1;
        }
    }

    // pair of arrays storing the mesh
    std::map<int, idx_t> eindex;
    std::vector<idx_t> eptr;
    std::vector<idx_t> eind;
    std::vector<int> etag;

    Element *ele = 0;
    auto &eles = domain->getElements();
    eptr.push_back(0);
    index = 0;
    while ((ele = eles()) != 0) {
        const auto &elenodes = ele->getExternalNodes();
        for (int i = 0; i < elenodes.Size(); ++i) {
            if (nind.find(elenodes(i)) == nind.end()) {
                opserr << "Process " << pid << ": element " << ele->getTag()
                << " has a node " << elenodes(i) 
                << " that does not exist in domain \n";
                return -1;
            }
            eind.push_back(nind[elenodes(i)]);
        }
        eptr.push_back((idx_t)eind.size());
        eindex[ele->getTag()] = index++;
        etag.push_back(ele->getTag());
    }

    // check size
    if (eptr.size() == 1) {
        opserr << "WARNING: no element is found on processorr 0\n";
        return -1;
    }

    // number of elements
    idx_t ne = (idx_t)eptr.size() - 1;

    // check if all processors have same number
    idx_t ne_max = 0;
    if (MPI_Reduce(&ne, &ne_max, 1, MPI_INT,
                   MPI_MAX, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        opserr << "WARNING: failed to get max ne\n";
        return -1;
    }
    if (pid == 0) {
        if (ne_max != ne) {
            opserr << "WARNING: all processors should "
                      "have the same mesh before partitioning";
            return -1;
        }
    }

    // partition for elements
    std::vector<idx_t> epart(ne);

    // partition for nodes
    std::vector<idx_t> npart(nn);

    // do partition on P0
    if (pid == 0) {

        // options
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;
        options[METIS_OPTION_NCUTS] = ncuts;
        options[METIS_OPTION_NITER] = niter;
        options[METIS_OPTION_NUMBERING] = 0;
        options[METIS_OPTION_UFACTOR] = ufactor;
        options[METIS_OPTION_DBGLVL] = info;

        // number of parts to partition
        idx_t nparts = np;

        // total communications volume
        idx_t objval;

        // call metis
        auto res = METIS_PartMeshNodal(
            &ne, &nn, &eptr[0], &eind[0], NULL, NULL, &nparts,
            NULL, options, &objval, &epart[0], &npart[0]);


        // check errors 
        if (res == METIS_ERROR_INPUT) {
            opserr << "WARNING: metis input error\n";
            return -1;
        } else if (res == METIS_ERROR_MEMORY) {
            opserr << "WARNING: metis memory error\n";
            return -1;
        } else if (res == METIS_ERROR) {
            opserr << "WARNING: metis error\n";
            return -1;
        }
    }

    // broadcast element partitions
    if (MPI_Bcast(&epart[0], ne, MPI_INT, 0, MPI_COMM_WORLD) !=
        MPI_SUCCESS) {
        opserr << "WARNING: failed to broadcast epart\n";
        return -1;
    }

    // broadcast node partitions
    if (MPI_Bcast(&npart[0], nn, MPI_INT, 0, MPI_COMM_WORLD) !=
        MPI_SUCCESS) {
        opserr << "WARNING: failed to broadcast npart\n";
        return -1;
    }

    // get nodes in current processor and remove elements
    for (int i = 0; i < ne; ++i) {

        // remove element
        if (epart[i] != pid) {
            Element *ele = domain->removeElement(etag[i]);
            if (ele != 0) {
                delete ele;
            }
            continue;
        }

        // get elenodes
        Element *ele = domain->getElement(etag[i]);
        const auto &elenodes = ele->getExternalNodes();
        for (int j = 0; j < elenodes.Size(); ++j) {
            auto &id = nind[elenodes(j)];
            if (id >= 0) {
                id = -id - 1;
            }
        }
    }

    // get nodes in current processor and remove mp
    auto &mps = domain->getMPs();
    MP_Constraint *mp = 0;
    while ((mp = mps()) != 0) {
        int rtag = mp->getNodeRetained();
        int ctag = mp->getNodeConstrained();
        auto &rid = nind[rtag];
        auto &cid = nind[ctag];
        if (rid < 0 || cid < 0) {
            if (rid >= 0) {
                rid = -rid - 1;
            }
            if (cid >= 0) {
                cid = -cid - 1;
            }
        } else {
            domain->removeMP_Constraint(mp->getTag());
            delete mp;
        }
    }

    // remove nodes
    for (const auto& item: nind) {
        int ndtag = item.first;
        auto id = item.second;
        if (id >= 0) {
            auto node = domain->removeNode(ndtag);
            if (node != 0) {
                delete node;
            }
            auto pc = domain->removePressure_Constraint(ndtag);
            if (pc != 0) {
                delete pc;
            }
        }
    }

    // remove sps
    auto& sps = domain->getSPs();
    SP_Constraint* sp = 0;
    while ((sp = sps()) != 0) {
        int ndtag = sp->getNodeTag();
        auto id = nind[ndtag];
        if (id >= 0) {
            domain->removeSP_Constraint(sp->getTag());
            delete sp;
        }
    }

    // go through load patterns
    auto& lps = domain->getLoadPatterns();
    LoadPattern* lp = 0;
    while ((lp = lps()) != 0) {
        // remove nodal loads
        auto& nloads = lp->getNodalLoads();
        NodalLoad* nload = 0;
        while ((nload = nloads()) != 0) {
            int ndtag = nload->getNodeTag();
            auto id = nind[ndtag];
            if (id >= 0) {
                lp->removeNodalLoad(nload->getTag());
                delete nload;
            } else {
                // nodal load can only appear in one place
                id = -id - 1;
                if (npart[id] != pid) {
                    lp->removeNodalLoad(nload->getTag());
                    delete nload;
                }
            }
        }

        // remove elemental loads
        auto& eloads = lp->getElementalLoads();
        ElementalLoad* eload = 0;
        while ((eload = eloads()) != 0) {
            int e = eload->getElementTag();
            auto id = epart[eindex[e]];
            if (id >= 0) {
                lp->removeElementalLoad(eload->getTag());
                delete eload;
            }
        }

        // remove sps
        auto& sps2 = lp->getSPs();
        while ((sp = sps2()) != 0) {
            int ndtag = sp->getNodeTag();
            auto id = nind[ndtag];
            if (id >= 0) {
                lp->removeSP_Constraint(sp->getTag());
                delete sp;
            }
        }
    }

#endif
    return 0;
}
