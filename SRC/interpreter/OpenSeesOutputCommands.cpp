/* *****************************************************************************
Copyright (c) 2015-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

// Written: Minjie

// Description: commands to output

#include <elementAPI.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <DOF_Group.h>
#include <Matrix.h>
#include <LoadPattern.h>
#include <FileStream.h>
#include <ID.h>
#include <Element.h>
#include <ElementIter.h>
#include <map>
#include <Recorder.h>
#include <Pressure_Constraint.h>
#include <vector>
#include <Parameter.h>
#include <ParameterIter.h>
#include <DummyStream.h>
#include <Response.h>
#include <Mesh.h>
#include <BackgroundMesh.h>
#include <Parameter.h>
#include <ParameterIter.h>

void* OPS_NodeRecorder();
void* OPS_EnvelopeNodeRecorder();
void* OPS_ElementRecorder();
void* OPS_EnvelopeElementRecorder();
void* OPS_PVDRecorder();
BackgroundMesh& OPS_getBgMesh();

//void* OPS_DriftRecorder();
//void* OPS_PatternRecorder();

namespace {

    struct char_cmp {
        bool operator () (const char *a, const char *b) const
        {
            return strcmp(a, b)<0;
        }
    };

    typedef std::map<const char *, void *(*)(void), char_cmp> OPS_ParsingFunctionMap;

    static OPS_ParsingFunctionMap recordersMap;

    static int setUpRecorders(void) {
        recordersMap.insert(std::make_pair("Node", &OPS_NodeRecorder));
        recordersMap.insert(std::make_pair("EnvelopeNode", &OPS_EnvelopeNodeRecorder));
        recordersMap.insert(std::make_pair("Element", &OPS_ElementRecorder));
        recordersMap.insert(std::make_pair("EnvelopeElement", &OPS_EnvelopeElementRecorder));
	recordersMap.insert(std::make_pair("PVD", &OPS_PVDRecorder));
	recordersMap.insert(std::make_pair("BgPVD", &OPS_PVDRecorder));
        //recordersMap.insert(std::make_pair("Drift", &OPS_DriftRecorder));
        //recordersMap.insert(std::make_pair("Pattern", &OPS_PatternRecorder));

        return 0;
    }
}

int OPS_Recorder()
{
    static bool initDone = false;
    if (initDone == false) {
        setUpRecorders();
        initDone = true;
    }

    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING too few arguments: recorder type? tag? ...\n";
        return -1;
    }

    const char* type = OPS_GetString();
    OPS_ParsingFunctionMap::const_iterator iter = recordersMap.find(type);
    if (iter == recordersMap.end()) {
        opserr << "WARNING recorder type " << type << " is unknown\n";
        return -1;
    }

    Recorder* theRecorder = (Recorder*)(*iter->second)();
    if (theRecorder == 0) {
        opserr << "WARNING failed to create recorder\n";
        return -1;
    }

    if (strcmp(type,"BgPVD") == 0) {
	BackgroundMesh& bg = OPS_getBgMesh();
	bg.addRecorder(theRecorder);
    } else {

	// now add the recorder to the domain
	Domain* theDomain = OPS_GetDomain();
	if (theDomain == 0) return -1;

	if (theDomain->addRecorder(*theRecorder) < 0) {
	    opserr << "ERROR could not add to domain - recorder.\n";
	    delete theRecorder;
	    return -1;
	}
    }

    return 0;
}

int OPS_nodeDisp()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING insufficient args: nodeDisp nodeTag <dof ...>\n";
	return -1;
    }

    // tag and dof
    int data[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) {
	numdata = 2;
    }

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING nodeDisp - failed to read int inputs\n";
	return -1;
    }
    data[1]--;

    // get response
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(data[0], Disp);

    if (nodalResponse == 0) {
	opserr << "WARNING no response is found\n";
	return -1;
    }

    // set outputs
    int size = nodalResponse->Size();
    if (data[1] >= 0) {
	if (data[1] >= size) {
	    opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
	    return -1;
	}

	double value = (*nodalResponse)(data[1]);
	numdata = 1;

	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr<<"WARNING nodeDisp - failed to read double inputs\n";
	    return -1;
	}


    } else {

	std::vector<double> values(size);
	for (int i=0; i<size; i++) {
	    values[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0]) < 0) {
	    opserr<<"WARNING nodeDisp - failed to read double inputs\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_nodeReaction()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeReaction nodeTag? <dof?>\n";
	return -1;
   }

    int data[2] = {0, -1};
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 2) numdata = 2;

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING nodeReaction - failed to read int inputs\n";
	return -1;
    }

    data[1]--;

    // get response
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(data[0], Reaction);

    if (nodalResponse == 0) {
	return -1;
    }

    int size = nodalResponse->Size();

    if (data[1] >= 0) {

      if (data[1] >= size) {
	  opserr << "WARNING nodeReaction nodeTag? dof? - dofTag? too large\n";
	  return -1;
      }

      double value = (*nodalResponse)(data[1]);
      numdata = 1;

      // now we copy the value to the tcl string that is returned
      if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	  opserr<<"WARNING nodeReaction - failed to set double output\n";
	  return -1;
      }

    } else {

	std::vector<double> values(size);
	for (int i=0; i<size; i++) {
	    values[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, &values[0]) < 0) {
	    opserr<<"WARNING nodeReaction - failed to set double output\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_nodeEigenvector()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
	return -1;
    }

    // tag, eigenvector, dof
    if (numdata > 3) numdata = 3;
    int data[3] = {0, 0, -1};

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr<<"WARNING invalid int inputs\n";
	return -1;
    }
    int eigenvector = data[1];
    int dof = data[2];

    dof--;
    eigenvector--;

    // get eigen vectors
    Node* theNode = theDomain->getNode(data[0]);
    const Matrix &theEigenvectors = theNode->getEigenvectors();

    int size = theEigenvectors.noRows();
    int numEigen = theEigenvectors.noCols();

    if (eigenvector < 0 || eigenvector >= numEigen) {
	opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
	return -1;
    }

    if (dof >= 0) {
	if (dof >= size) {
	    opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
	    return -1;
	}

	double value = theEigenvectors(dof, eigenvector);
	size = 1;

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&size, &value) < 0) {
	    opserr<<"WARNING nodeEigenvector - failed to set double output\n";
	    return -1;
	}

    } else {

	Vector values(size);
	for (int i=0; i<size; i++) {
	    values(i) = theEigenvectors(i, eigenvector);
	}

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&size, &values(0)) < 0) {
	    opserr<<"WARNING nodeEigenvector - failed to set double output\n";
	    return -1;
	}
    }

    return 0;
}

int OPS_getTime()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    double time = theDomain->getCurrentTime();
    int numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &time) < 0) {
	opserr << "WARNING failed to get current time\n";
	return -1;
    }

    return 0;
}

int OPS_eleResponse()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 2) {
	opserr << "WARNING want - eleResponse eleTag? eleArgs...\n";
	return -1;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "could not read eleTag\n";
	return -1;
    }

    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 0) {
	const char** argv = new const char*[numdata];
	for (int i=0; i<numdata; i++) {
	    argv[i] = OPS_GetString();
	}
	const Vector* data = theDomain->getElementResponse(tag, argv, numdata);
	delete [] argv;

	if (data != 0) {
	    int size = data->Size();
	    double* newdata = new double[size];
	    for (int i=0; i<size; i++) {
		newdata[i] = (*data)(i);
	    }
	    if (OPS_SetDoubleOutput(&size, newdata) < 0) {
		opserr << "WARNING failed to et response\n";
		delete [] newdata;
		return -1;
	    }
	    delete [] newdata;

	}
    }
    return 0;

}

int OPS_getLoadFactor()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
	return -1;
    }

    int pattern;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &pattern) < 0) {
	opserr << "ERROR reading load pattern tag -- getLoadFactor\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    LoadPattern *thePattern = theDomain->getLoadPattern(pattern);
    if (thePattern == 0) {
	opserr << "ERROR load pattern with tag " << pattern << " not found in domain -- getLoadFactor\n";
	return -1;
    }

    double factor = thePattern->getLoadFactor();
    if (OPS_SetDoubleOutput(&numdata, &factor) < 0) {
	opserr << "WARNING failed to set load factor\n";
	return -1;
    }

    return 0;
}

// printNode():
// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// input: nodeArg: integer equal to arg count to node plus 1
//        output: output stream to which the results are sent
//
int printNode(OPS_Stream &output)
{
    int flag = 0; // default flag sent to a nodes Print() method
    int nodeArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // if just 'print <filename> node' print all the nodes - no flag
    if (argc == 0) {
        NodeIter &theNodes = theDomain->getNodes();
        Node *theNode;
        while ((theNode = theNodes()) != 0) {
            theNode->Print(output);
        }
        return 0;
    }

    // if 'print <filename> node flag int <int int ..>' get the flag
    const char* flagArg = OPS_GetString();
    if ((strcmp(flagArg, "flag") == 0) || (strcmp(flagArg, "-flag") == 0)) {
        // get the specified flag
        if (argc < 2) {
            opserr << "WARNING print <filename> node <flag int> no int specified \n";
            return -1;
        }
        int numdata = 1;
        if (OPS_GetIntInput(&numdata, &flag) < 0) {
            opserr << "WARNING print node failed to get integer flag: \n";
            return -1;
        }
        nodeArg += 2;
    }
    else {
        OPS_ResetCurrentInputArg(2);
    }

    // now print the nodes with the specified flag, 0 by default

    // if 'print <filename> node flag'
    //     print out all the nodes in the domain with flag
    if (argc == nodeArg) {
        NodeIter &theNodes = theDomain->getNodes();
        Node *theNode;
        while ((theNode = theNodes()) != 0) {
            theNode->Print(output, flag);
        }
        return 0;
    }
    else {
        // otherwise print out the specified nodes i j k .. with flag
        int numNodes = argc - nodeArg;
        ID *theNodes = new ID(numNodes);
        for (int i = 0; i < numNodes; i++) {
            int nodeTag;
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &nodeTag) < 0) {
                opserr << "WARNING print node failed to get integer: " << endln;
                delete theNodes;
                return -1;
            }
            (*theNodes)(i) = nodeTag;
            nodeArg++;
        }
        theDomain->Print(output, theNodes, 0, flag);
        delete theNodes;
    }

    return 0;
}

int printElement(OPS_Stream &output)
{
    int flag = 0; // default flag sent to a nodes Print() method
    int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // if just 'print <filename> node' print all the nodes - no flag
    if (argc == 0) {
        ElementIter &theElements = theDomain->getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output);
        }
        return 0;
    }

    // if 'print <filename> Element flag int <int int ..>' get the flag
    const char* eleflag = OPS_GetString();
    if ((strcmp(eleflag, "flag") == 0) || (strcmp(eleflag, "-flag")) == 0) {
        // get the specified flag
        if (argc < 2) {
            opserr << "WARNING print <filename> ele <flag int> no int specified \n";
            return -1;
        }
        int numdata = 1;
        if (OPS_GetIntInput(&numdata, &flag) < 0) {
            opserr << "WARNING print ele failed to get integer flag: \n";
            return -1;
        }
        eleArg += 2;
    }
    else {
        OPS_ResetCurrentInputArg(2);
    }

    // now print the Elements with the specified flag, 0 by default
    if (argc == eleArg) {
        ElementIter &theElements = theDomain->getElements();
        Element *theElement;
        while ((theElement = theElements()) != 0) {
            theElement->Print(output, flag);
        }
        return 0;
    }
    else {
        // otherwise print out the specified nodes i j k .. with flag
        int numEle = argc - eleArg;
        ID *theEle = new ID(numEle);
        for (int i = 0; i < numEle; i++) {
            int eleTag;
            int numdata = 1;
            if (OPS_GetIntInput(&numdata, &eleTag) < 0) {
                opserr << "WARNING print ele failed to get integer: " << endln;
                delete theEle;
                return -1;
            }
            (*theEle)(i) = eleTag;
        }
        theDomain->Print(output, 0, theEle, flag);
        delete theEle;
    }

    return 0;
}

int printAlgorithm(OPS_Stream &output)
{
    /*int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    EquiSolnAlgo** theAlgorithm = OPS_GetAlgorithm();
    if (theAlgorithm == 0) return -1;

    // if just 'print <filename> algorithm'- no flag
    if (argc == 0) {
        theAlgorithm->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &flag) < 0) {
        opserr << "WARNING print algorithm failed to get integer flag: \n";
        return -1;
    }
    theAlgorithm->Print(output, flag);*/

    return 0;
}

int printIntegrator(OPS_Stream &output)
{
    /*int eleArg = 0;
    int argc = OPS_GetNumRemainingInputArgs();

    StaticIntegrator** theStaticIntegrator = OPS_GetStaticIntegrator();
    TransientIntegrator** theTransientIntegrator = OPS_GetTransientIntegrator();

    if (theStaticIntegrator == 0 && theTransientIntegrator == 0)
        return 0;

    IncrementalIntegrator *theIntegrator;
    if (theStaticIntegrator != 0)
        theIntegrator = theStaticIntegrator;
    else
        theIntegrator = *theTransientIntegrator;

    // if just 'print <filename> algorithm'- no flag
    if (argc == 0) {
        theIntegrator->Print(output);
        return 0;
    }

    // if 'print <filename> Algorithm flag' get the flag
    int flag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &flag) < 0) {
        opserr << "WARNING print integrator failed to get integer flag: \n";
        return -1;
    }
    theIntegrator->Print(output, flag);*/

    return 0;
}

int OPS_printModelGID()
{
    // This function print's a file with node and elements in a format useful for GID
    int res = 0;
    bool hasLinear = 0;
    bool hasTri3  = 0;
    bool hasQuad4 = 0;
    bool hasQuad8 = 0;
    bool hasQuad9 = 0;
    bool hasBrick = 0;
    int startEle = 1;
    int endEle = 1;
    int eleRange = 0;
    int numdata = 1;

    FileStream outputFile;
    // OPS_Stream *output = &opserr;
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING printGID fileName? - no filename supplied\n";
	return -1;
    }
    const char* filename = OPS_GetString();
    openMode mode = OVERWRITE;
    while (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if (strcmp(flag,"-append") == 0) {
	    mode = APPEND;
	}
	if (strcmp(flag,"-eleRange") == 0 && OPS_GetNumRemainingInputArgs()>1) {
	    //opserr<<"WARNING:commands: eleRange defined"<<endln;
	    eleRange = 1;
	    if (OPS_GetIntInput(&numdata, &startEle) < 0) {
		opserr << "WARNING print node failed to get integer: "  << endln;
		return -1;
	    }
	    if (OPS_GetIntInput(&numdata, &endEle) < 0) {
		opserr << "WARNING print node failed to get integer: " << endln;
		return -1;
	    }
	    //opserr<<"startEle = "<<startEle<<" endEle = "<<endEle<<endln;
	}
    }

    if (outputFile.setFile(filename, mode) < 0) {
        opserr << "WARNING printGID " << filename << " failed to set the file\n";
	return -1;
    }

    // Cycle over Elements to understand what type of elements are there
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
	// int tag = theElement->getTag();

	// Check type of Element with Number of Nodes
	// if 2 Nodes print the Element
	int nNode = theElement->getNumExternalNodes();
	if (nNode == 2) {
	    hasLinear = 1;
	} else if (nNode == 4) {
	    hasQuad4 = 1;
	} else if (nNode == 3) {
	    hasTri3 = 1;
	} else if (nNode == 9) {
	    hasQuad9 = 1;
	} else if (nNode == 8) {
	    const char *name = theElement->getClassType();
	    if (strcmp(name,"Brick") == 0) {
		hasBrick = 1;
	    } else {
		hasQuad8 = 1;
	    }
	}
    }
    // **** Linear Elements - 2 Nodes
    if (hasLinear == 1) {
	// Print HEADER
	outputFile << "MESH \"2NMESH\" dimension 3 ElemType Linear Nnode 2" << endln;
	outputFile << "#color 0 0 255" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {
		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 2) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    Vector tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Quadrilateral Elements - 4 Nodes
    if (hasQuad4 == 1) {
	// Print HEADER
	outputFile << "MESH \"4NMESH\" dimension 3 ElemType Quadrilateral Nnode 4" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 4) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    Vector tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Triangular Elements - 3 Nodes
    if (hasTri3 == 1) {
	// Print HEADER
	outputFile << "MESH \"3NMESH\" dimension 3 ElemType Triangle Nnode 3" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {

		// Check type of Element with Number of Nodes
		// if 3 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 3) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    Vector tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Quadrilateral Elements - 9 Nodes
    if (hasQuad9 == 1) {
	// Print HEADER
	outputFile << "MESH \"9NMESH\" dimension 3 ElemType Linear Nnode 9" << endln;
	outputFile << "#color 0 255 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange ==0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 9) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    Vector tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }
    // **** Hexahedra Elements - 8 Nodes
    if (hasBrick == 1) {
	// Print HEADER
	outputFile << "MESH \"8NMESH\" dimension 3 ElemType Hexahedra Nnode 8" << endln;
	outputFile << "#color 255 0 0" << endln << endln;

	// Print node coordinates
	outputFile << "Coordinates" << endln;
	NodeIter &theNodes = theDomain->getNodes();
	// MeshRegion *myRegion = theDomain->getRegion(0);
	Node *theNode;
	while ((theNode = theNodes()) != 0) {
	    int tag = theNode->getTag();
	    const Vector &crds = theNode->getCrds();
	    //outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" << crds(2) << endln;
	    int l_tmp = crds.Size();
	    outputFile << tag << "\t\t";
	    for (int ii = 0; ii<l_tmp; ii++) {
		outputFile << crds(ii) << "\t";
	    };
	    for (int ii = l_tmp; ii<3; ii++) {
		outputFile << 0.0 << "\t";
	    };
	    outputFile << endln;
	}
	outputFile << "End coordinates" << endln << endln;

	// Print elements connectivity
	outputFile << "Elements" << endln;
	ElementIter &theElements = theDomain->getElements();
	Element *theElement;
	while ((theElement = theElements()) != 0) {
	    int tag = theElement->getTag();
	    // Check if element tag is inside theRange
	    if (((tag<=endEle) & (tag>=startEle)) || (eleRange == 0)) {

		// Check type of Element with Number of Nodes
		// if 2 Nodes print the Element
		int nNode = theElement->getNumExternalNodes();
		if (nNode == 8) {
		    Node **NodePtrs;
		    NodePtrs = theElement->getNodePtrs();
		    Vector tagNodes(nNode);
		    for (int i = 0; i < nNode; i++) {
			tagNodes(i)=NodePtrs[i]->getTag();
		    }
		    outputFile << tag << "\t\t";
		    for (int i = 0; i < nNode; i++) {
			outputFile << tagNodes(i) << "\t";
		    }
		    outputFile << endln;
		}
	    }
	}
	outputFile << "End elements" << endln;
    }

    outputFile.close();
    return res;
}

int OPS_eleForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    /*
      Element *theEle = theDomain.getElement(tag);
      if (theEle == 0)
      return TCL_ERROR;

      const Vector &force = theEle->getResistingForce();
    */

    const char *myArgv[1];
    char myArgv0[8];
    strcpy(myArgv0,"forces");
    myArgv[0] = myArgv0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *force = theDomain->getElementResponse(tag, &myArgv[0], 1);
    if (force != 0) {
	int size = force->Size();

	if (dof >= 0) {

	    if (size < dof) {
		opserr << "WARNING eleForce dof > size\n";
		return -1;
	    }

	    double value = (*force)(dof);

	    // now we copy the value to the tcl string that is returned
	    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
		opserr << "WARNING eleForce failed to set output\n";
		return -1;
	    }

	} else {

	    double* data = new double[size];
	    for (int i=0; i<size; i++) {
		data[i] = (*force)(i);
	    }
	    if (OPS_SetDoubleOutput(&size, data) < 0) {
		opserr << "WARNING eleForce failed to set outputs\n";
		delete [] data;
		return -1;
	    }

	    delete [] data;
	    return 0;
	}
    }

    return 0;
}

int OPS_eleDynamicalForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleForce eleTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theEle = theDomain->getElement(tag);
    if (theEle == 0) {
	opserr << "WARNING element "<<tag<<" does not exist\n";
	return -1;
    }

    const Vector &force = theEle->getResistingForceIncInertia();
    int size = force.Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING eleDyanmicalForce size < dof\n";
	    return -1;
	}

	double value = force(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = force(i);
	}
	if (OPS_SetDoubleOutput(&size, data) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;

}

int OPS_nodeUnbalance()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeUnbalance nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING eleForce eleTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Unbalance);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeUnbalance size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING nodeUnbalance failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data) < 0) {
	    opserr << "WARNING eleDyanmicalForce failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeVel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeVel nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeVel nodeTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING nodeVel nodeTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Vel);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeVel size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING nodeVel failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data) < 0) {
	    opserr << "WARNING nodeVel failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeAccel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeAccel nodeTag? <dof?>\n";
	return -1;
    }

    int tag;
    int dof = -1;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeAccel nodeTag? dof? - could not read nodeTag? \n";
	return -1;
    }

    if (OPS_GetNumRemainingInputArgs() > 0) {
	if (OPS_GetIntInput(&numdata, &dof) < 0) {
	    opserr << "WARNING nodeAccel nodeTag? dof? - could not read dof? \n";
	    return -1;
	}
    }

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, Accel);
    if (nodalResponse == 0) {
	opserr << "WARNING failed to get nodal response\n";
	return -1;
    }

    int size = nodalResponse->Size();

    if (dof >= 0) {

	if (size < dof) {
	    opserr << "WARNING nodeAccel size < dof\n";
	    return -1;
	}

	double value = (*nodalResponse)(dof);

	// now we copy the value to the tcl string that is returned
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING nodeAccel failed to set output\n";
	    return -1;
	}

    } else {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = (*nodalResponse)(i);
	}
	if (OPS_SetDoubleOutput(&size, data) < 0) {
	    opserr << "WARNING nodeAccel failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
	return 0;
    }

    return 0;
}

int OPS_nodeResponse()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
	return -1;
    }

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING nodeResponse - could not read int inputs \n";
	return -1;
    }

    int tag=data[0], dof=data[1], responseID=data[2];

    dof--;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *nodalResponse = theDomain->getNodeResponse(tag, (NodeResponseType)responseID);
    if (nodalResponse == 0 || nodalResponse->Size() < dof || dof < 0) {
	opserr << "WARNING errors in read response\n";
	return -1;
    }

    double value = (*nodalResponse)(dof);
    numdata = 1;

    // now we copy the value to the tcl string that is returned
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_nodeCoord()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeCoord nodeTag? <dim?>\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeCoord nodeTag? dim? - could not read nodeTag? \n";
	return -1;
    }

    int dim = -1;

    if (OPS_GetNumRemainingInputArgs() > 0) {
	const char* flag = OPS_GetString();
	if (strcmp(flag,"X") == 0 || strcmp(flag,"x") == 0) {
	    dim = 0;
	} else if (strcmp(flag,"Y") == 0 || strcmp(flag,"y") == 0) {
	    dim = 1;
	} else if (strcmp(flag,"Z") == 0 || strcmp(flag,"z") == 0) {
	    dim = 2;
	} else {
	    OPS_ResetCurrentInputArg(-1);
	    if (OPS_GetIntInput(&numdata, &dim) < 0) {
		opserr << "WARNING nodeCoord nodeTag? dim? - could not read dim? \n";
		return -1;
	    }
	    dim--;
	}
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Node *theNode = theDomain->getNode(tag);

    if (theNode == 0) {
	opserr << "WARNING node "<<tag<<" does not exist\n";
	return -1;
    }

    const Vector &coords = theNode->getCrds();

    int size = coords.Size();
    if (dim == -1) {

	double* data = new double[size];
	for (int i=0; i<size; i++) {
	    data[i] = coords(i);
	}
	if (OPS_SetDoubleOutput(&size, data) < 0) {
	    opserr << "WARNING failed to set output\n";
	    delete [] data;
	    return -1;
	}
	delete [] data;

    } else if (dim < size) {
	double value = coords(dim); // -1 for OpenSees vs C indexing
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING failed to set output\n";
	    return -1;
	}

    } else {
	opserr << "WARNING invalid dim\n";
	return -1;
    }

    return 0;
}

int OPS_setNodeCoord()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - setNodeCoord nodeTag? dim? value?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read nodeTag? \n";
	return -1;
    }

    int dim;
    double value;

    if (OPS_GetIntInput(&numdata, &dim) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read dim? \n";
	return -1;
    }
    if (OPS_GetDoubleInput(&numdata, &value) < 0) {
	opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read value? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    Node *theNode = theDomain->getNode(tag);

    if (theNode == 0) {
	opserr << "WARNING node "<< tag<<" does not exist\n";
	return -1;
    }

    Vector coords(theNode->getCrds());
    coords(dim-1) = value;
    theNode->setCrds(coords);

    return 0;
}

int OPS_updateElementDomain()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    // Need to "setDomain" to make the change take effect.
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
	theElement->setDomain(theDomain);
    }

    return 0;
}

int OPS_eleNodes()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - eleNodes eleTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING eleNodes eleTag? \n";
	return -1;
    }

    const char *myArgv[1];
    char myArgv0[80];
    strcpy(myArgv0,"nodeTags");
    myArgv[0] = myArgv0;

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    const Vector *tags = theDomain->getElementResponse(tag, &myArgv[0], 1);
    //  Element *theElement = theDomain.getElement(tag);
    if (tags != 0) {
	int numTags = tags->Size();
	int* data = new int[numTags];
	for (int i = 0; i < numTags; i++) {
	    data[i] = (int)(*tags)(i);
	}

	if (OPS_SetIntOutput(&numTags, data) < 0) {
	    opserr << "WARNING failed to set outputs\n";
	    delete [] data;
	    return -1;
	}

	delete [] data;
    }

    return 0;
}

int OPS_nodeDOFs()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - nodeDOFs nodeTag?\n";
	return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeDOFs nodeTag?\n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
	opserr << "WARNING nodeDOFs node " << tag << " not found" << endln;
	return -1;
    }
    int numDOF = theNode->getNumberDOF();

    DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
    if (theDOFgroup == 0) {
      opserr << "WARNING nodeDOFs DOF group null" << endln;
      return -1;
    }
    const ID &eqnNumbers = theDOFgroup->getID();
    int *data = new int[numDOF];
    for (int i = 0; i < numDOF; i++) {
      data[i] = eqnNumbers(i);
    }
    if (OPS_SetIntOutput(&numDOF, data) < 0) {
      opserr << "WARNING nodeDOFs failed to set outputs\n";
      delete [] data;
      return -1;
    }
    delete [] data;

    return 0;
}

int OPS_nodeMass()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - nodeMass nodeTag? nodeDOF?\n";
	return -1;
    }

    int tag, dof;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
	return -1;
    }
    if (OPS_GetIntInput(&numdata, &dof) < 0) {
	opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Node *theNode = theDomain->getNode(tag);
    if (theNode == 0) {
	opserr << "WARNING nodeMass node " << tag << " not found" << endln;
	return -1;
    }
    int numDOF = theNode->getNumberDOF();
    if (dof < 1 || dof > numDOF) {
	opserr << "WARNING nodeMass dof " << dof << " not in range" << endln;
	return -1;
    }
    else {
	const Matrix &mass = theNode->getMass();
	double value = mass(dof-1,dof-1);
	if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	    opserr << "WARNING nodeMass failed to set mass\n";
	}
    }

    return 0;
}

int OPS_nodePressure()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING: want - nodePressure nodeTag?\n";
        return -1;
    }

    int tag;
    int numdata = 1;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING: nodePressure invalid tag\n";
        return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    double pressure = 0.0;
    Pressure_Constraint* thePC = theDomain->getPressure_Constraint(tag);
    if(thePC != 0) {
        pressure = thePC->getPressure();
    }
    if (OPS_SetDoubleOutput(&numdata, &pressure) < 0) {
	opserr << "WARNING failed to get presure\n";
	return -1;
    }

    return 0;
}

int OPS_nodeBounds()
{

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;
    const Vector &bounds = theDomain->getPhysicalBounds();
    int size = bounds.Size();
    double* data = new double[size];
    for (int i = 0; i < size; i++)
      data[i] = bounds(i);

    if (OPS_SetDoubleOutput(&size, data) < 0) {
	opserr << "WARNING failed to get node bounds\n";
	delete [] data;
	return -1;
    }

    delete [] data;

    return 0;
}

int OPS_setPrecision()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING setPrecision precision? - no precision value supplied\n";
	return -1;
    }
    int precision;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &precision) < 0) {
	opserr << "WARNING setPrecision precision? - error reading precision value supplied\n";
	return -1;
    }
    opserr.setPrecision(precision);

    return 0;
}

int OPS_getEleTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    std::vector<int> eletags;
    if (OPS_GetNumRemainingInputArgs() < 1) {
	// return all elements

	Element *theEle;
	ElementIter &eleIter = theDomain->getElements();

	while ((theEle = eleIter()) != 0) {
	    eletags.push_back(theEle->getTag());
	}
    } else if (OPS_GetNumRemainingInputArgs() == 2) {

	// return nodes in mesh
	const char* type = OPS_GetString();
	if (strcmp(type,"-mesh") == 0) {
	    int tag;
	    int num = 1;
	    if (OPS_GetIntInput(&num, &tag) < 0) {
		opserr << "WARNING: failed to get mesh tag\n";
		return -1;
	    }
	    Mesh* msh = OPS_getMesh(tag);
	    if (msh == 0) {
		opserr << "WARNING: mesh "<<tag<<" does not exist\n";
		return -1;
	    }
	    const ID& tags = msh->getEleTags();
	    for (int i=0; i<tags.Size(); ++i) {
		eletags.push_back(tags(i));
	    }
	}
    }

    if (eletags.empty()) return 0;

    int size = (int)eletags.size();
    int* data = &eletags[0];

    if (OPS_SetIntOutput(&size, data) < 0) {
	opserr << "WARNING failed to set outputs\n";
	return -1;
    }

    return 0;
}

int OPS_getNodeTags()
{
    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    std::vector<int> nodetags;
    if (OPS_GetNumRemainingInputArgs() < 1) {

	// return all nodes
	Node *theNode;
	NodeIter &nodeIter = theDomain->getNodes();

	while ((theNode = nodeIter()) != 0) {
	    nodetags.push_back(theNode->getTag());
	}
    } else if (OPS_GetNumRemainingInputArgs() == 2) {

	// return nodes in mesh
	const char* type = OPS_GetString();
	if (strcmp(type,"-mesh") == 0) {
	    int tag;
	    int num = 1;
	    if (OPS_GetIntInput(&num, &tag) < 0) {
		opserr << "WARNING: failed to get mesh tag\n";
		return -1;
	    }
	    Mesh* msh = OPS_getMesh(tag);
	    if (msh == 0) {
		opserr << "WARNING: mesh "<<tag<<" does not exist\n";
		return -1;
	    }
	    const ID& tags = msh->getNodeTags();
	    for (int i=0; i<tags.Size(); ++i) {
		nodetags.push_back(tags(i));
	    }
	    const ID& newtags = msh->getNewNodeTags();
	    for (int i=0; i<newtags.Size(); ++i) {
		nodetags.push_back(newtags(i));
	    }
	}
    }

    if (nodetags.empty()) return 0;

    int size = (int)nodetags.size();
    int* data = &nodetags[0];

    if (OPS_SetIntOutput(&size, data) < 0) {
	opserr << "WARNING failed to set outputs\n";
	return -1;
    }

    return 0;
}

int OPS_sectionForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sectionForce eleTag? secNum? dof? \n";
	return -1;
    }

    //opserr << "sectionForce: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 3;
    int data[3];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionForce eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];
    int dof = data[2];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "force";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    if (dof <= 0 || dof > theVec.Size()) {
	opserr << "WARNING invalid dof "<<dof<<"\n";
	delete theResponse;
	return -1;
    }

    double value = theVec(dof-1);
    numdata = 1;

    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionDeformation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sectionDeformation eleTag? secNum? dof? \n";
	return -1;
    }

    //opserr << "sectionDeformation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 3;
    int data[3];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionDeformation eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];
    int dof = data[2];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionDeformation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "deformation";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    if (dof <= 0 || dof > theVec.Size()) {
	opserr << "WARNING invalid dof "<<dof<<"\n";
	delete theResponse;
	return -1;
    }

    double value = theVec(dof-1);
    numdata = 1;

    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionStiffness()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionStiffness eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionStiffness: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionStiffness eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionStiffness element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "stiffness";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMat = *(info.theMatrix);
    int nsdof = theMat.noCols();
    int size = nsdof*nsdof;
    if (size == 0) {
	delete theResponse;
	return 0;
    }

    std::vector<double> values;
    values.reserve(size);

    for (int i = 0; i < nsdof; i++) {
	for (int j = 0; j < nsdof; j++) {
	    values.push_back(theMat(i,j));
	}
    }

    if (OPS_SetDoubleOutput(&size, &values[0]) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionFlexibility()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionFlexibility eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionFlexibility: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionFlexibility eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionFlexibility element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 3;
    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "flexibility";
    const char *argvv[3];
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMat = *(info.theMatrix);
    int nsdof = theMat.noCols();
    int size = nsdof*nsdof;
    if (size == 0) {
	delete theResponse;
	return 0;
    }

    std::vector<double> values;
    values.reserve(size);

    for (int i = 0; i < nsdof; i++) {
	for (int j = 0; j < nsdof; j++) {
	    values.push_back(theMat(i,j));
	}
    }

    if (OPS_SetDoubleOutput(&size, &values[0]) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionLocation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionLocation eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionLocation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionLocation eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionLocation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "integrationPoints";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    if (secNum <= 0 || secNum > theVec.Size()) {
	opserr << "WARNING invalid secNum\n";
	delete theResponse;
	return -1;
    }

    double value = theVec(secNum-1);
    numdata = 1;

    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_sectionWeight()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sectionWeight eleTag? secNum? \n";
	return -1;
    }

    //opserr << "sectionWeight: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 2;
    int data[2];

    if (OPS_GetIntInput(&numdata, data) < 0) {
	opserr << "WARNING sectionWeight eleTag? secNum? dof? - could not read int input? \n";
	return -1;
    }

    int tag = data[0];
    int secNum = data[1];

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sectionWeight element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "integrationWeights";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    if (secNum <= 0 || secNum > theVec.Size()) {
	opserr << "WARNING invalid secNum\n";
	delete theResponse;
	return -1;
    }

    double value = theVec(secNum-1);
    numdata = 1;

    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_basicDeformation()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicDeformation eleTag? \n";
	return -1;
    }

    //opserr << "basicDeformation: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicDeformation eleTag? dofNum? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicDeformation element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicDeformation";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    int nbf = theVec.Size();

    std::vector<double> data(nbf);
    for (int i=0; i<nbf; i++) {
	data[i] = theVec(i);
    }

    if (OPS_SetDoubleOutput(&nbf, &data[0]) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_basicForce()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicForce eleTag? \n";
	return -1;
    }

    //opserr << "basicForce: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicForce eleTag? dofNum? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicForce";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	double res = 0.0;
	if (OPS_SetDoubleOutput(&numdata, &res) < 0) {
	    opserr << "WARNING: failed to set output\n";
	    return -1;
	}
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Vector &theVec = *(info.theVector);
    int nbf = theVec.Size();

    std::vector<double> data(nbf);
    for (int i=0; i<nbf; i++) {
	data[i] = theVec(i);
    }

    if (OPS_SetDoubleOutput(&nbf, &data[0]) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_basicStiffness()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 1) {
	opserr << "WARNING want - basicStiffness eleTag? \n";
	return -1;
    }

    //opserr << "basicStiffness: ";
    //for (int i = 0; i < argc; i++)
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = 1;
    int tag;

    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING basicStiffness eleTag? dofNum? - could not read eleTag? \n";
	return -1;
    }

    Domain* theDomain = OPS_GetDomain();
    if (theDomain == 0) return -1;

    Element *theElement = theDomain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING basicStiffness element with tag " << tag << " not found in domain \n";
	return -1;
    }

    int argcc = 1;
    char a[80] = "basicStiffness";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	return 0;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMatrix = *(info.theMatrix);
    int nbf = theMatrix.noCols();

    std::vector<double> values;
    int size = nbf*nbf;
    if (size == 0) return 0;
    values.reserve(size);


    for (int i = 0; i < nbf; i++) {
	for (int j = 0; j < nbf; j++) {
	    values.push_back(theMatrix(i,j));
	}
    }

    if (OPS_SetDoubleOutput(&size, &values[0]) < 0) {
	opserr << "WARNING failed to set output\n";
	delete theResponse;
	return -1;
    }

    delete theResponse;

    return 0;
}

int OPS_version()
{
    if (OPS_SetString(OPS_VERSION) < 0) {
	opserr << "WARNING failed to set version string\n";
	return -1;
    }

    return 0;
}

int OPS_logFile()
{
    if (OPS_GetNumRemainingInputArgs() < 1) { 
	opserr << "WARNING logFile fileName? - no filename supplied\n";
	return -1;
    }
    openMode mode = OVERWRITE;
    bool echo = true;

    const char* filename = OPS_GetString();
    if (strcmp(filename, "Invalid String Input!") == 0) {
	opserr << "WARNING: invalid string input\n";
	return -1;
    }

    while (OPS_GetNumRemainingInputArgs() > 0) {

	const char* opt = OPS_GetString();
	
	if (strcmp(opt,"-append") == 0) {
	    mode = APPEND;
	} else if (strcmp(opt,"-noEcho") == 0) {
	    echo = false;
	}
    }

    if (opserr.setFile(filename, mode, echo) < 0) {
	opserr << "WARNING logFile " << filename << " failed to set the file\n";
	return -1;
    }

    // const char *pwd = getInterpPWD(interp);  
    // simulationInfo.addOutputFile(argv[1], pwd);

    return 0;
}

// Sensitivity:BEGIN /////////////////////////////////////////////
int OPS_sensNodeDisp()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeDisp nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeDisp: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeDisp: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getDispSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensNodeVel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeVel nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeVel: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeVel: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getVelSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensNodeAccel()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensNodeAccel nodeTag? dof? paramTag?\n";
	return -1;
    }    

    int data[3];
    int numdata = 3;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag, dof or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;


    Node *theNode = domain->getNode(data[0]);
    if (theNode == 0) {
	opserr << "sensNodeAccel: node " << data[0] << " not found" << "\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[2]);
    if (theParam == 0) {
	opserr << "sensNodeAccel: parameter " << data[2] << " not found" << "\n";
	return -1;
    }

    int gradIndex = theParam->getGradIndex();

    double value = theNode->getAccSensitivity(data[1],gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &value) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensLambda()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING no load pattern supplied -- getLoadFactor\n";
	return -1;
    }

    int data[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to read patternTag or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    LoadPattern *thePattern = domain->getLoadPattern(data[0]);
    if (thePattern == 0) {
	opserr << "ERROR load pattern with tag " << data[0] << " not found in domain\n";
	return -1;
    }

    Parameter *theParam = domain->getParameter(data[1]);
    if (theParam == 0) {
	opserr << "sensLambda: parameter " << data[1] << " not found" << "\n";
	return -1;
    }
  
    int gradIndex = theParam->getGradIndex();
    double factor = thePattern->getLoadFactorSensitivity(gradIndex);

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &factor) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

int OPS_sensSectionForce()
{
    if (OPS_GetNumRemainingInputArgs() < 3) {
	opserr << "WARNING want - sensSectionForce eleTag? <secNum?> dof? paramTag?\n";
	return -1;
    }    
  
    //opserr << "sensSectionForce: ";
    //for (int i = 0; i < argc; i++) 
    //  opserr << argv[i] << ' ' ;
    //opserr << endln;

    int numdata = OPS_GetNumRemainingInputArgs();
    std::vector<int> data(numdata);
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to read input data\n";
	return -1;
    }

    int tag, dof, paramTag;
    int secNum = -1;
    if (numdata == 3) {
	tag = data[0];
	dof = data[1];
	paramTag = data[2]; 
    } else {
	tag = data[0];
	secNum = data[1];
	dof = data[2];
	paramTag = data[3]; 
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    ParameterIter &pIter = domain->getParameters();
    Parameter *theParam;
    while ((theParam = pIter()) != 0) {
	theParam->activate(false);
    }

    theParam = domain->getParameter(paramTag);
    int gradIndex = theParam->getGradIndex();
    theParam->activate(true);

    Element *theElement = domain->getElement(tag);
    if (theElement == 0) {
	opserr << "WARNING sensSectionForce element with tag " << tag << " not found in domain \n";
	return -1;
    }

    char a[80] = "section";
    char b[80];
    sprintf(b, "%d", secNum);
    char c[80] = "dsdh";
    const char *argvv[3];
    int argcc = 3;
    argvv[0] = a;
    argvv[1] = b;
    argvv[2] = c;
    if (secNum < 0) { // For zeroLengthSection
	argcc = 2;
	argvv[1] = c;
    }

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == 0) {
	numdata = 1;
	double res = 0.0;
	if (OPS_SetDoubleOutput(&numdata, &res) < 0) {
	    opserr<<"WARNING failed to set output\n";
	    return -1;
	}
	return 0;
    }

    theResponse->getResponseSensitivity(gradIndex);
    Information &info = theResponse->getInformation();

    Vector theVec = *(info.theVector);

    numdata = theVec.Size();
    if (OPS_SetDoubleOutput(&numdata, &theVec(dof-1)) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    theParam->activate(false);

    delete theResponse;

    return 0;
}

int OPS_sensNodePressure()
{
    // make sure at least one other argument to contain type of system
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr << "WARNING want - sensNodePressure nodeTag? paramTag?\n";
	return -1;
    }    

    int data[2];
    int numdata = 2;
    if (OPS_GetIntInput(&numdata, &data[0]) < 0) {
	opserr << "WARNING: failed to get tag or paramTag\n";
	return -1;
    }

    Domain* domain = OPS_GetDomain();
    if (domain == 0) return 0;

    double dp = 0.0;
    Pressure_Constraint* thePC = domain->getPressure_Constraint(data[0]);
    if(thePC != 0) {
        // int ptag = thePC->getPressureNode();
        // Node* pNode = theDomain.getNode(ptag);
        Node* pNode = thePC->getPressureNode();
        if(pNode != 0) {

            Parameter *theParam = domain->getParameter(data[1]);
            if (theParam == 0) {
                opserr << "sensNodePressure: parameter " << data[1] << " not found" << endln;
                return -1;
            }

            int gradIndex = theParam->getGradIndex();
            dp = pNode->getVelSensitivity(1,gradIndex);
        }
    }

    numdata = 1;
    if (OPS_SetDoubleOutput(&numdata, &dp) < 0) {
	opserr<<"WARNING failed to set output\n";
	return -1;
    }

    return 0;
}

// Sensitivity:END /////////////////////////////////////////////
