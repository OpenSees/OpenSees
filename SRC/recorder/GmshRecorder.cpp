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

// Written: Jose Abell
//
// Description: This file contains the class definition for 
// GmshRecorder. A GmshRecorder is used to store mesh and all responses in gmsh format.

#include "GmshRecorder.h"
#include <sstream>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <Matrix.h>
#include <classTags.h>
#include <NodeIter.h>
#include <Message.h>
#include <Channel.h>


#ifdef _PARALLEL_PROCESSING
#include <MachineBroker.h>
#include <Graph.h>
#include <Vertex.h>
#include <VertexIter.h>
#include <mpi.h>
#endif

#ifdef _PARALLEL_PROCESSING
extern MachineBroker *theMachineBroker;
#endif



// #ifdef GMSHRECORDER_DEBUG
// #define DEBUGSTREAM std::cerr
// #else 
// #define DEBUGSTREAM 
// #endif

#define GMSHRECORDER_DEBUG_DISABLED 1
#define DEBUGSTREAM \
    if (GMSHRECORDER_DEBUG_DISABLED) {} \
    else opserr

std::map<int, GmshRecorder::GmshType> GmshRecorder::gmshtypes;

void* OPS_GmshRecorder()
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 1) {
        opserr << "WARNING: insufficient number of arguments\n";
        return 0;
    }

    GmshRecorder::setGMSHType();

    // filename
    const char* name = OPS_GetString();

    // plotting options
    numdata = OPS_GetNumRemainingInputArgs();
    int indent = 2;
    int precision = 10;
    int write_graph_mesh = 0;
    int write_update_time = 0;
    int write_ele_updatetime = 0;
    GmshRecorder::NodeData nodedata;
    std::vector<GmshRecorder::EleData> eledata;
    while (numdata > 0) {
        std::string type = OPS_GetString();
        if (type == "disp") {
            nodedata.disp = true;
        } else if (type == "vel") {
            nodedata.vel = true;
        } else if (type == "accel") {
            nodedata.accel = true;
        } else if (type == "incrDisp") {
            nodedata.incrdisp = true;
        } else if (type == "reaction") {
            nodedata.reaction = true;
        } else if (type == "pressure") {
            nodedata.pressure = true;
        } else if (type == "graph") {
            write_graph_mesh = true;
        } else if (type == "updatetime") {
            write_update_time = true;
        } else if (type == "unbalancedLoad") {
            nodedata.unbalanced = true;
        } else if (type == "mass") {
            nodedata.mass = true;
        } else if (type == "eleupdatetime") {
            write_ele_updatetime = true;
        // } else if (type == "partition") {
            // write_partition = true;
        }  else if (type == "eigen") {
            numdata = OPS_GetNumRemainingInputArgs();
            if (numdata < 1) {
                opserr << "WARNING: eigen needs 'numEigenvector'\n";
                return 0;
            }
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &nodedata.numeigen) < 0) return 0;
        } else if (type == "-precision") {
            numdata = OPS_GetNumRemainingInputArgs();
            if (numdata < 1) {
                opserr << "WARNING: needs precision \n";
                return 0;
            }
            numdata = 1;
            if (OPS_GetIntInput(&numdata, &precision) < 0) return 0;
        } else if (type == "eleResponse") {
            numdata = OPS_GetNumRemainingInputArgs();
            if (numdata < 1) {
                opserr << "WANRING: elementResponse needs 'argc','argv'\n";
                return 0;
            }
            GmshRecorder::EleData edata;
            numdata = OPS_GetNumRemainingInputArgs();
            edata.resize(numdata);
            for (int i = 0; i < numdata; i++) {
                edata[i] = OPS_GetString();
                opserr << edata[i].c_str() << endln;
            }
            eledata.push_back(edata);
        }
        numdata = OPS_GetNumRemainingInputArgs();
    }

    // create recorder
    return new GmshRecorder(name, nodedata, eledata, indent, precision, write_graph_mesh, write_update_time, write_ele_updatetime);
}

GmshRecorder::GmshRecorder(const char *name, const NodeData& ndata,
                           const std::vector<EleData>& edata, int ind, int pre, int write_graph_mesh_, int write_update_time_, int write_ele_updatetime_)
    : Recorder(RECORDER_TAGS_GmshRecorder),  precision(pre),
      write_header_now(true), write_mesh_now(true), write_binary_mode(false), write_ele_updatetime(write_ele_updatetime_),
      filename(name), 
      timestep(), timeparts(), theFile(), nodedata(ndata), eledata(edata), theDomain(NULL), current_step(0),
      write_graph_mesh(write_graph_mesh_), write_update_time(write_update_time_)
{
    DEBUGSTREAM << "edata.size() = " << (int)edata.size() << endln;
    DEBUGSTREAM << "eledata.size() = " << (int)eledata.size() << endln;
}

GmshRecorder::GmshRecorder()
    : Recorder(RECORDER_TAGS_GmshRecorder),
    precision(0), 
    write_header_now(true), write_mesh_now(true), write_binary_mode(false), write_ele_updatetime(false),
    filename(),
    timestep(), timeparts(), 
    theFile(), theDomain(NULL), current_step(0),
    write_graph_mesh(0)
{
        GmshRecorder::setGMSHType();
}


GmshRecorder::~GmshRecorder()
{
    theFile.close();
}

int
GmshRecorder::record(int ctag, double timestamp)
{

    if (precision == 0)
        return 0;

    // get current time
    timestep.push_back(timestamp);

    if(write_update_time)
    {
        write_update_time_now();
        if(write_ele_updatetime)
        {
            write_mesh(); 
            write_eleupdatetime_now();
        }
        return 0;
    }


    // if(write_graph_mesh)
    if(false) // Disabled for now
    {
        write_element_graph();
    }
    else
    {
        write_mesh(); 
    
        if( nodedata.disp || 
            nodedata.vel || 
            nodedata.accel || 
            nodedata.incrdisp || 
            nodedata.reaction || 
            nodedata.rayleigh || 
            nodedata.pressure || 
            nodedata.unbalanced || 
            nodedata.mass)
            write_node_data();
        
        // opserr << "Write element output?  ";
        if (eledata.size() > 0)
        {
            // opserr << "  YES!\n";
            write_element_data();
        }
        else
        {
            // opserr << "  No!\n";
        }
    }

    current_step ++;

    return 0;
}

int
GmshRecorder::restart()
{
    timestep.clear();
    timeparts.clear();
    return 0;
}

int
GmshRecorder::domainChanged()
{
    return 0;
}

int
GmshRecorder::setDomain(Domain& domain)
{
    DEBUGSTREAM << "GmshRecorder::setDomain(Domain& domain)" << endln;
    theDomain = &domain;
    return 0;
}

int GmshRecorder::write_header()
{
    if (write_header_now)
    {
        theFile << "$MeshFormat\n"
                << "2.2 0 8\n";
        if(write_binary_mode)
        {
            int one = 1;
            theFile.write((const char *)&one, sizeof(int));
        }
        theFile << "$EndMeshFormat\n";
        write_header_now = false; // Don't do this again
    }

    return 0;
}



int
GmshRecorder::write_mesh()
{

    if(!write_mesh_now)
    {
        return 0;
    }
    int rank = 0;
    int nproc = 1;
#ifdef _PARALLEL_PROCESSING
    rank = theMachineBroker->getPID();  
    nproc = theMachineBroker->getNP();  
#endif

    DEBUGSTREAM << " Saving mesh for rank ---> " << rank << endln;

    if (theDomain == 0) {
        opserr << "WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << 0 << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;


    // open file
    theFile.close();
    std::stringstream ssmshname;
    ssmshname << filename << ".mesh." << rank << ".msh";
    std::string mshname = ssmshname.str();

    DEBUGSTREAM << "Opening " << mshname.c_str() << endln;


    theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out);
    if (theFile.fail()) {
        opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
        return -1;
    }
    theFile.precision(precision);
    theFile << std::scientific;

    write_header();

    theFile << "$Nodes\n";


    theFile << theDomain->getNumNodes() << "\n";
    
    // Write Node coordinates
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = 0;
    while ((theNode = theNodes()) != 0) {
        theFile << theNode->getTag() << ' ';
        const Vector& crds = theNode->getCrds();
        for (int j = 0; j < 3; j++) {
            if (j < crds.Size()) {
                theFile << crds(j) << ' ';
            } else {
                theFile << 0.0 << ' ';
            }
        }
        theFile << std::endl;        
    }
    theFile << "$EndNodes\n";


    theFile << "$Elements\n";




    // // connectivity
    int numel;

    if( rank == 0)
    {
        numel = theDomain->getNumElements() - (nproc - 1); // rank 0 holds pointers to subdomains as "elements"
    }
    else
    {
        numel = theDomain->getNumElements();
    }
    theFile <<  numel << "\n";

    ElementIter* eiter = &(theDomain->getElements());
    Element* theEle = 0;
    while ((theEle = (*eiter)()) != 0) {
        int classTag = theEle->getClassTag();
        int tag = theEle->getTag();
        int gmsheletype = gmshtypes[classTag];
        int ntags = 2;
        int phys = rank+1;
        int enti = rank+1;

        // DEBUGSTREAM << rank << ": " << tag << " " << classTag << " " << gmsheletype << endln;

        if(gmsheletype > 0)
        {
            theFile << tag << " "
                << gmsheletype << " "
                << ntags << " "
                << phys << " "
                << enti << " ";
        
                const ID& elenodes = theEle->getExternalNodes();
                int nnodes = elenodes.Size();
                for (int j = 0; j < nnodes; j++) 
                {
                    theFile << elenodes(j) << " ";
                }
                theFile << '\n';
        }
    }
    theFile << "$EndElements\n";

    theFile.close();
    write_header_now = true;
    write_mesh_now = false;

    return 0;
}

int
GmshRecorder::write_node_data()
{
    
    int rank = 0;

#ifdef _PARALLEL_PROCESSING
    rank = theMachineBroker->getPID();  
#endif

    DEBUGSTREAM << "GmshRecorder -- Saving mesh for rank ---> " << rank << endln;

    if (theDomain == 0) {
        opserr << "GmshRecorder::write_data() - WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    std::string viewname(" ");
    if(nodedata.disp) viewname = std::string("\"Displacement\"");
    if(nodedata.vel) viewname = std::string("\"Velocity\"");
    if(nodedata.accel) viewname = std::string("\"Acceleration\"");
    if(nodedata.incrdisp) viewname = std::string("\"Incremental Displacement\"");
    if(nodedata.reaction) viewname = std::string("\"Reaction Forces\"");
    if(nodedata.rayleigh) viewname = std::string("\"Rayleigh Damping Forces\"");
    if(nodedata.pressure) viewname = std::string("\"Pressure\"");
    if(nodedata.unbalanced) viewname = std::string("\"Unbalanced Forces\"");
    if(nodedata.mass) viewname = std::string("\"Mass\"");

    if(!theFile.is_open())
    {
        std::stringstream ss;
    
        // open file
        ss << filename << ".";

        if(nodedata.disp) ss << "disp";
        if(nodedata.vel) ss << "vel";
        if(nodedata.accel) ss << "accel";
        if(nodedata.incrdisp) ss << "incrdisp";
        if(nodedata.reaction) ss << "reaction";
        if(nodedata.rayleigh) ss << "rayleigh";
        if(nodedata.pressure) ss << "pressure";
        if(nodedata.unbalanced) ss << "unbalanced";
        if(nodedata.mass) ss << "mass";

        ss << "." << rank << ".msh";
      
        std::string mshname = ss.str();

        if(write_binary_mode)
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
        }
        else
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
        }


        if (theFile.fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theFile.precision(precision);
        theFile << std::scientific;

        write_header();
    }   

    // Data header
    theFile << "$NodeData\n";
    theFile << "1\n";
    theFile << viewname << std::endl;
    theFile << "1\n";
    theFile << timestep.back() << std::endl;
    theFile << "3\n";
    theFile << current_step << std::endl;
    theFile << "3\n";
    theFile << theDomain->getNumNodes() << "\n";
    
    // Write data
    NodeIter& theNodes = theDomain->getNodes();
    Node* theNode = 0;
    while ((theNode = theNodes()) != 0) {
        theFile << theNode->getTag() << ' ';
        if(nodedata.disp) write_data_line(theFile, theNode->getDisp());
        if(nodedata.vel) write_data_line(theFile, theNode->getVel());
        if(nodedata.accel) write_data_line(theFile, theNode->getAccel());
        if(nodedata.incrdisp) write_data_line(theFile, theNode->getIncrDisp());
        if(nodedata.reaction) write_data_line(theFile, theNode->getReaction());
        if(nodedata.unbalanced) write_data_line(theFile, theNode->getUnbalancedLoad());
    }
    theFile << "$EndNodeData\n";

    return 0;
}

int GmshRecorder::write_data_line(std::ofstream &s, const Vector & data, const int truncatesize)
{
        for (int j = 0; j < truncatesize; j++) {
            if (j < data.Size()) {
                s << data(j) << ' ';
            } else {
                s << 0.0 << ' ';
            }
        }
        s << std::endl;        
        return 0;
}


int
GmshRecorder::write_eleupdatetime_now()
{

#ifdef GMSHRECORDER_JAABELL_EXPERIMENTAL
  
    int rank = 0; 
    int nproc = 1;

#ifdef _PARALLEL_PROCESSING
    rank =theMachineBroker->getPID();  
    nproc =theMachineBroker->getNP();  
#endif 

    DEBUGSTREAM << "GmshRecorder -- write_eleupdatetime_now---> " << rank << endln;

    if (theDomain == 0) {
        opserr << "GmshRecorder::write_data() - WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << 0 << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;


    // open file
    std::stringstream ssmshname;
    ssmshname << filename << ".eleupdatetime." << rank << ".msh";
    std::string mshname = ssmshname.str();
  

    DEBUGSTREAM << "mshname = " << mshname.c_str() << endln;

    if(not theFile.is_open())
    {
        if(write_binary_mode)
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
        }
        else
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out );
        }


        if (theFile.fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theFile.precision(precision);
        theFile << std::scientific;

        write_header();
    }


    // Data header
    theFile << "$ElementData\n";
    theFile << "1\n";
    theFile << "\"Element Update Time (s)\"" << std::endl;
    theFile << "1\n";
    theFile << timestep.back() << std::endl;
    theFile << "3\n";
    theFile << current_step << std::endl;
    theFile << "1" << std::endl;

    int numel;
    if( rank == 0)
    {
        numel = theDomain->getNumElements() ; // rank 0 holds pointers to subdomains as "elements"
    }
    else
    {
        numel = theDomain->getNumElements();
    }
    theFile << numel << "\n";
    
    ElementIter &elementIter = theDomain->getElements();
    Element* theElement = 0; 
    while ( (theElement=elementIter()) != 0) 
    {
        int tag = theElement->getTag();
        double eleupdatetime =theElement->getTime();
        theFile << tag << ' ' << eleupdatetime << std::endl;
    } 
    theFile << "$EndElementData\n";



    // Data header
    theFile << "$ElementData\n";
    theFile << "1\n";
    theFile << "\"Element Partition\"" << std::endl;
    theFile << "1\n";
    theFile << timestep.back() << std::endl;
    theFile << "3\n";
    theFile << current_step << std::endl;
    theFile << "1" << std::endl;

    theFile << numel << "\n";
    
    ElementIter &elementIter2 = theDomain->getElements();
    theElement = 0; 
    while ( (theElement=elementIter2()) != 0) 
    {
        int tag = theElement->getTag();
        // double eleupdatetime =theElement->getTime();
        theFile << tag << ' ' << rank << std::endl;
    } 
    theFile << "$EndElementData\n";



    current_step++;

#endif

    return 0;
}


int
GmshRecorder::write_element_data()
{
    
    int rank = 0; 
    int nproc = 1;

#ifdef _PARALLEL_PROCESSING
    rank =theMachineBroker->getPID();  
    nproc =theMachineBroker->getNP();  
#endif 

    DEBUGSTREAM << "GmshRecorder -- Saving mesh for rank ---> " << rank << endln;

    if (theDomain == 0) {
        opserr << "GmshRecorder::write_data() - WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    // get time and part
    std::stringstream ss;
    ss.precision(precision);
    ss << std::scientific;
    ss << 0 << ' ' << timestep.back();
    std::string stime, spart;
    ss >> spart >> stime;


    // open file
    std::stringstream ssmshname;
    ssmshname << filename << ".eledata." << rank << ".msh";
    std::string mshname = ssmshname.str();
  

    DEBUGSTREAM << "mshname = " << mshname.c_str() << endln;

    if(!theFile.is_open())
    {
        if(write_binary_mode)
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
        }
        else
        {
            theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out | std::ios::binary);
        }


        if (theFile.fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theFile.precision(precision);
        theFile << std::scientific;

        write_header();
    }

  // element response
    DEBUGSTREAM << "eledata.size() = " << (int) eledata.size() << endln;
    // for(int i=0; i<(int)eledata.size(); i++) 
    for(int i=0; i<1; i++) 
    {
        // DEBUGSTREAM << "i = " << i << endln;
        // check data
        int argc = (int) eledata[i].size();
        if(argc == 0) 
            continue;

        std::vector<const char*> argv(argc);
        for(int j=0; j<argc; j++) 
        {
            // DEBUGSTREAM << "j = " << j << endln;
            argv[j] = eledata[i][j].c_str();
        }


        ElementIter &elementIter = theDomain->getElements();
        Element* theElement =  elementIter();
        const Vector* data =theDomain->getElementResponse(theElement->getTag(), &(argv[0]), argc);


        while(data == NULL)
        {
            theElement =  elementIter();
            data =theDomain->getElementResponse(theElement->getTag(), &(argv[0]), argc);
        }
        int datasize = data->Size();
        elementIter = theDomain->getElements();
        if(rank==0)
        {
            for(int d=0;d<nproc-1;d++)
                data =theDomain->getElementResponse(theElement->getTag(), &(argv[0]), argc);
        }

        // DEBUGSTREAM << "data = " << data << endln;
        // DEBUGSTREAM << "theElement = " << theElement << endln;

        if(datasize == 2)
        {
            datasize = 3;
        } 
        else if (datasize > 3)
        {
            datasize = 9;
        }

        // Data header
        theFile << "$ElementData\n";
        theFile << "1\n";
        theFile << "\"Element Data\"" << std::endl;
        theFile << "1\n";
        theFile << timestep.back() << std::endl;
        theFile << "3\n";
        theFile << current_step << std::endl;
        theFile << datasize << std::endl;

        int numel;
        if( rank == 0)
        {
            numel = theDomain->getNumElements() - (nproc - 1); // rank 0 holds pointers to subdomains as "elements"
        }
        else
        {
            numel = theDomain->getNumElements();
        }
        theFile << numel << "\n";
        
        while ( theElement != 0) 
        {
            int tag = theElement->getTag();
            const Vector* data =theDomain->getElementResponse(tag, &(argv[0]),argc);
            theFile << tag << ' ';
            if(data != NULL)
                write_data_line(theFile, *data, datasize);
            else
            {
                for(int d=0;d<datasize;d++)
                    theFile<<0<<" ";
                theFile << endln;
            }
            theElement =  elementIter();
        } 
        theFile << "$EndElementData\n";
    }

    current_step++;

    return 0;
}






int
GmshRecorder::write_element_graph()
{
    static bool do_only_once = true;

    int rank = 0; 
    int nproc = 1; 

#ifdef _PARALLEL_PROCESSING
    rank = theMachineBroker->getPID();  
    nproc = theMachineBroker->getNP();  

    // if (rank > 0)
    // {
    //     return 0;  // Only do this on proc0
    // }

    // DEBUGSTREAM << "Saving mesh for rank ---> " << rank << endln;

    if (theDomain == 0) {
        opserr << "WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    if(do_only_once)        
    {
        // get time and part
        std::stringstream ss;
        ss.precision(precision);
        ss << std::scientific;
        ss << 0 << ' ' << timestep.back();
        std::string stime, spart;
        ss >> spart >> stime;


        // open file
        theFile.close();
        std::stringstream ssmshname;
        ssmshname << filename << ".graphnodes." << rank << ".msh";
        std::string mshname = ssmshname.str();

        // DEBUGSTREAM << "Opening " << mshname.c_str() << endln;


        theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out);
        if (theFile.fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theFile.precision(precision);
        theFile << std::scientific;

        write_header();


        Graph& theGraph = theDomain->getElementGraph();
        int numvertex = theGraph.getNumVertex();
        int numedge = theGraph.getNumEdge();

        theFile << "$Nodes\n";

        theFile << numvertex << endln;

        //Iterate the graph, get the vertices (which correspond to the elements)
        VertexIter vertices = theGraph.getVertices();
        Vertex *v = vertices();
        while(v != 0)
        {
            int eletag = v->getRef();

            //Compute this element's centroid for graph placement
            Element *e = theDomain->getElement(eletag);
            const ID& thisElementNodes = e->getExternalNodes();

            int nelenodes = thisElementNodes.Size();

            double x0=0, y0=0, z0=0;
            for(int n = 0; n < nelenodes; n++)
            {
                int nodetag = thisElementNodes(n);
                Node* thisNode = theDomain->getNode(nodetag);
                const Vector &crds = thisNode->getCrds();
                if(crds.Size() >= 1) 
                    x0 += crds(0)/nelenodes;
                if(crds.Size() >= 2) 
                    y0 += crds(1)/nelenodes;
                if(crds.Size() >= 3) 
                    z0 += crds(2)/nelenodes;
            }

            theFile << eletag << " " << x0 << " " << y0 << " " << z0 << endln;

            v = vertices();
        }
        theFile << "$EndNodes\n";


        theFile << "$Elements\n";
        // theFile <<  numvertex  << "\n";
        theFile <<  numvertex  << "\n";

        vertices.reset();
        v = vertices();
        int maxvertextag = 0;
        while(v != 0)
        {
            int eletag = v->getRef();
            if (eletag > maxvertextag)
                maxvertextag = eletag;
            theFile << eletag << " " << GMSH_VERTEX << " " << "2 " << rank+1 << " " << rank+1 << " " << eletag << endln; 
            v = vertices();
        }

        theFile << "$EndElements\n";

        theFile.close();

        std::stringstream ssmshname2;
        ssmshname2 << filename << ".graphelements." << rank << ".msh";
        mshname = ssmshname2.str();

        theFile.open(mshname.c_str(), std::ios::trunc | std::ios::out);
        if (theFile.fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theFile.precision(precision);
        theFile << std::scientific;

        int numvertex_proc[nproc];


        int start_tag = maxvertextag;

#ifdef _PARALLEL_PROCESSING
        numvertex_proc[rank] = numvertex;
        // maxvertextag
        int maxvertextag_all;
        MPI_Allreduce(
            &maxvertextag,
            &maxvertextag_all,
            1,
            MPI_INT,
            MPI_MAX,
            MPI_COMM_WORLD);

        MPI_Allgather(
            &numvertex,
            1,
            MPI_INT,
            numvertex_proc,
            1,
            MPI_INT,
            MPI_COMM_WORLD);


        // DEBUGSTREAM << "(" << rank << ")";
        for(int r = 0; r<rank; r++)
        {   
            maxvertextag_all += numvertex_proc[r];
            // DEBUGSTREAM << " " << numvertex_proc[r] << "    " << maxvertextag_all << endln;
        }

        start_tag = maxvertextag_all;

        DEBUGSTREAM << "Rank " << rank << " gets tags from " << start_tag << " to " << start_tag + numvertex << endln;

#endif


        write_header_now = true;
        write_header();
        theFile << "$Elements\n";
        theFile <<  numvertex  << "\n";
        vertices.reset();
        v = vertices();
        int neweletag = start_tag;
        while(v != 0)
        {
            int eletag = v->getRef();
            const ID& adjacent_vertices = v->getAdjacency();
            for (int i = 0; i <= adjacent_vertices.Size(); i++)
            {
                theFile << neweletag << " " << GMSH_LINE << " " << "2 " << rank+1 << " " << rank+1 << " " << eletag << " " << adjacent_vertices(i) <<endln; 
                neweletag++;
            }
            v = vertices();
        }
        theFile << "$EndElements\n";
        
        do_only_once = false;
    }
#endif
    return 0;
}


int
GmshRecorder::write_update_time_now()
{

#ifdef GMSHRECORDER_JAABELL_EXPERIMENTAL

    static bool do_only_once = true;
    int rank = 0; 
    int nproc = 1; 
    static std::ofstream* theUpdateTimeFilePtr = 0;

#ifdef _PARALLEL_PROCESSING
    rank = theMachineBroker->getPID();  
    nproc = theMachineBroker->getNP();  
#endif
    if (theDomain == 0) {
        opserr << "WARNING: setDomain has not been called -- GmshRecorder\n";
        return -1;
    }

    if(do_only_once)        
    {


        // get time and part
        std::stringstream ss;
        ss.precision(precision);
        ss << std::scientific;
        ss << 0 << ' ' << timestep.back();
        std::string stime, spart;
        ss >> spart >> stime;


        // open file
        // theFile.close();
        std::stringstream ssmshname;
        ssmshname << filename << ".updatetime." << rank << "." << nproc << ".out";
        std::string mshname = ssmshname.str();

        // DEBUGSTREAM << "Opening " << mshname.c_str() << endln;
        theUpdateTimeFilePtr = new std::ofstream();

        theUpdateTimeFilePtr->open(mshname.c_str(), std::ios::trunc | std::ios::out);
        if (theUpdateTimeFilePtr->fail()) {
            opserr << "WARNING: Failed to open file " << mshname.c_str() << "\n";
            return -1;
        }
        theUpdateTimeFilePtr->precision(precision);
        *theUpdateTimeFilePtr << std::scientific;

        do_only_once = false;
    }

    double time = theDomain->getCurrentTime();
    double updatetime = theDomain->getUpdateTime();

    *theUpdateTimeFilePtr << time << " " << updatetime << endln;
#endif

    return 0;
}




int
GmshRecorder::sendSelf(int commitTag, Channel &theChannel)
{
    int length = filename.size();
    char * string = new char[length];
    strcpy(string, filename.c_str());
    int n_eledata = eledata.size();

    ID data(16);
    data(0) = length;
    data(1) = precision;
    data(2) = nodedata.disp;
    data(3) = nodedata.vel;
    data(4) = nodedata.accel;
    data(5) = nodedata.incrdisp;
    data(6) = nodedata.reaction;
    data(7) = nodedata.rayleigh;
    data(8) = nodedata.pressure;
    data(9) = nodedata.unbalanced;
    data(10) = nodedata.mass;
    data(11) = write_binary_mode;
    data(12) = n_eledata;
    data(13) = write_graph_mesh;
    data(14) = write_update_time;
    data(15) = write_ele_updatetime;

    ID sizesdata(n_eledata);
    for(int i = 0; i < n_eledata; i++)
    {
        sizesdata(i) = eledata[i].size();
    }


    if(theChannel.sendID(this->getDbTag(), commitTag, data) < 0)
    {
        opserr << "GmshRecorder::sendSelf - Problem sending data\n";
    }
    DEBUGSTREAM << "Sent " << data << endln;

    if(theChannel.sendID(this->getDbTag(), commitTag, sizesdata) < 0)
    {
        opserr << "GmshRecorder::sendSelf - Problem sending sizesdata\n";
    }

    Message msg(string, length);
    if( theChannel.sendMsg(this->getDbTag(), commitTag, msg) < 0)
    {
        opserr << "GmshRecorder::sendSelf - Problem sending filename\n";
    }
    DEBUGSTREAM << "Sent filename=" << filename.c_str() << endln;
    

    delete [] string;

    for(int i = 0; i < n_eledata; i++)
    {
        int send_how_many = sizesdata(i);
        for(int j = 0; j < send_how_many; j++)
        {
            std::string eledatastring = eledata[i][j];
            int string_length = eledatastring.size();
            char * charbuffer = new char[string_length];
            strcpy(charbuffer, eledatastring.c_str());

            static ID string_length_vec(1);
            string_length_vec(0) = string_length;

            if( theChannel.sendID(this->getDbTag(), commitTag, string_length_vec) < 0)
            {
                opserr << "GmshRecorder::sendSelf - Problem sending string_length_vec " << eledatastring.c_str() << endln;
            } 

            Message msg_string(charbuffer, string_length);
            if( theChannel.sendMsg(this->getDbTag(), commitTag, msg_string) < 0)
            {
                opserr << "GmshRecorder::sendSelf - Problem sending eledata " << eledatastring.c_str() << endln;
            }        

            DEBUGSTREAM << "Sent: eledatastring = " << eledatastring.c_str() << endln;

            delete [] charbuffer;
        }
    }
    

    return 0;
}

int
GmshRecorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    ID data(16);
    int length = 0;
    int n_eledata;
    if(theChannel.recvID(this->getDbTag(), commitTag, data) < 0)
    {
        opserr << "GmshRecorder::recvSelf - Problem sending filename string length\n";
    }    
    length = data(0);
    precision = data(1);
    nodedata.disp = data(2);
    nodedata.vel = data(3);
    nodedata.accel = data(4);
    nodedata.incrdisp = data(5);
    nodedata.reaction = data(6);
    nodedata.rayleigh = data(7);
    nodedata.pressure = data(8);
    nodedata.unbalanced = data(9);
    nodedata.mass = data(10);
    write_binary_mode = data(11);
    n_eledata = data(12);
    write_graph_mesh = data(13);
    write_update_time = data(14);
    write_ele_updatetime = data(15);

    DEBUGSTREAM << "Received data = " << data << endln;

    ID sizesdata(n_eledata);

    if(theChannel.recvID(this->getDbTag(), commitTag, sizesdata) < 0)
    {
        opserr << "GmshRecorder::recvSelf - Problem receiving sizesdata\n";
    } 
    DEBUGSTREAM << "Received sizesdata = " << sizesdata << endln;

    char * string = new char[length];

    Message msg(string, length );
    if( theChannel.recvMsg(this->getDbTag(), commitTag, msg) < 0)
    {
        opserr << "GmshRecorder::recvSelf - Problem receiving the filename\n";
    }
    filename = std::string(string, length);
    DEBUGSTREAM << "Received filename = " << filename.c_str() << endln;

    delete [] string;

    for(int i = 0; i < n_eledata; i++)
    {
        int recv_how_many = sizesdata(i);
        GmshRecorder::EleData edata;

        for(int j = 0; j < recv_how_many; j++)
        {

            static ID string_length_vec(1);

            if( theChannel.recvID(this->getDbTag(), commitTag, string_length_vec) < 0)
            {
                opserr << "GmshRecorder::recvSelf - Problem receiving string_length_vec " << endln;
            }        
 

            int string_length = string_length_vec(0);
            DEBUGSTREAM << "Received: string_length = " << string_length << endln;

            char * charbuffer = new char[string_length];
            Message msg_string(charbuffer, string_length);
            if( theChannel.recvMsg(this->getDbTag(), commitTag, msg_string) < 0)
            {
                opserr << "GmshRecorder::recvSelf - Problem receiving eledata "  << endln;
            }

            std::string eledatastring(charbuffer, string_length);

            DEBUGSTREAM << "Received: eledatastring = " << eledatastring.c_str()<< endln;

            edata.push_back(eledatastring);

            delete [] charbuffer ;
        }
        eledata.push_back(edata);
    }


    return 0;
}

void
GmshRecorder::setGMSHType()
{
    if (gmshtypes.empty() == false) {
        return;
    }
    gmshtypes[ELE_TAG_Subdomain] = GMSH_POLY_VERTEX;
    gmshtypes[ELEMENT_TAGS_WrapperElement] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_ElasticBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ModElasticBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElasticBeam3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_Beam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_beam2d02] = GMSH_LINE;
    gmshtypes[ELE_TAG_beam2d03] = GMSH_LINE;
    gmshtypes[ELE_TAG_beam2d04] = GMSH_LINE;
    gmshtypes[ELE_TAG_beam3d01] = GMSH_LINE;
    gmshtypes[ELE_TAG_beam3d02] = GMSH_LINE;
    gmshtypes[ELE_TAG_Truss] = GMSH_LINE;
    gmshtypes[ELE_TAG_TrussSection] = GMSH_LINE;
    gmshtypes[ELE_TAG_CorotTruss] = GMSH_LINE;
    gmshtypes[ELE_TAG_CorotTrussSection] = GMSH_LINE;
    gmshtypes[ELE_TAG_fElmt05] = GMSH_LINE;
    gmshtypes[ELE_TAG_fElmt02] = GMSH_LINE;
    gmshtypes[ELE_TAG_MyTruss] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLength] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthSection] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthND] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthContact2D] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthContact3D] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthContactASDimplex] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthContactNTS2D] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthInterface2D] = GMSH_LINE;
    gmshtypes[ELE_TAG_CoupledZeroLength] = GMSH_LINE;
    gmshtypes[ELE_TAG_ZeroLengthRocking] = GMSH_LINE;
    gmshtypes[ELE_TAG_NLBeamColumn2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_NLBeamColumn3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_LargeDispBeamColumn3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_FourNodeQuad] = GMSH_QUAD;
    gmshtypes[ELE_TAG_FourNodeQuad3d] = GMSH_QUAD;
    gmshtypes[ELE_TAG_Tri31] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_SixNodeTri] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_BeamWithHinges2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_BeamWithHinges3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_EightNodeBrick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_TwentyNodeBrick] = GMSH_QUADRATIC_HEXAHEDRON;
    gmshtypes[ELE_TAG_EightNodeBrick_u_p_U] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_TwentyNodeBrick_u_p_U] = GMSH_QUADRATIC_HEXAHEDRON;
    gmshtypes[ELE_TAG_FourNodeQuadUP] = GMSH_QUAD;
    gmshtypes[ELE_TAG_TotalLagrangianFD20NodeBrick] = GMSH_QUADRATIC_HEXAHEDRON;
    gmshtypes[ELE_TAG_TotalLagrangianFD8NodeBrick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_EightNode_LDBrick_u_p] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_EightNode_Brick_u_p] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_TwentySevenNodeBrick] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_BrickUP] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_Nine_Four_Node_QuadUP] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_Twenty_Eight_Node_BrickUP] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_Twenty_Node_Brick] = GMSH_QUADRATIC_HEXAHEDRON;
    gmshtypes[ELE_TAG_BBarFourNodeQuadUP] = GMSH_QUAD;
    gmshtypes[ELE_TAG_BBarBrickUP] = GMSH_QUAD;
    gmshtypes[ELE_TAG_PlateMITC4] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ShellMITC4] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ShellMITC9] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_ASDShellQ4] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ASDShellT3] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_Plate1] = GMSH_QUAD;
    gmshtypes[ELE_TAG_Brick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_BbarBrick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_FLBrick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_EnhancedQuad] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ConstantPressureVolumeQuad] = GMSH_QUAD;
    gmshtypes[ELE_TAG_NineNodeMixedQuad] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_NineNodeQuad] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_EightNodeQuad] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_DispBeamColumn2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TimoshenkoBeamColumn2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_DispBeamColumn3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_DispBeamColumnWarping3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_HingedBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_HingedBeam3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TwoPointHingedBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TwoPointHingedBeam3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_OnePointHingedBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_OnePointHingedBeam3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_BeamColumnJoint2d] = GMSH_QUAD;
    gmshtypes[ELE_TAG_BeamColumnJoint3d] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ForceBeamColumn2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ForceBeamColumnWarping2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ForceBeamColumn3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElasticForceBeamColumn2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElasticForceBeamColumnWarping2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElasticForceBeamColumn3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ForceBeamColumnCBDI2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ForceBeamColumnCBDI3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_DispBeamColumn2dInt] = GMSH_LINE;
    gmshtypes[ELE_TAG_InternalSpring] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SimpleJoint2D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_Joint2D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_Joint3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_ElastomericBearingPlasticity3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElastomericBearingPlasticity2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TwoNodeLink] = GMSH_LINE;
    gmshtypes[ELE_TAG_ActuatorCorot] = GMSH_LINE;
    gmshtypes[ELE_TAG_Actuator] = GMSH_LINE;
    gmshtypes[ELE_TAG_Adapter] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_ElastomericBearingBoucWen2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElastomericBearingBoucWen3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_FlatSliderSimple2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_FlatSliderSimple3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_FlatSlider2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_FlatSlider3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_SingleFPSimple2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_SingleFPSimple3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_SingleFP2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_SingleFP3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_DoubleFPSimple2d] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_DoubleFPSimple3d] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_DoubleFP2d] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_DoubleFP3d] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_TripleFPSimple2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TripleFPSimple3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TripleFP2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TripleFP3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_MultiFP2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_MultiFP3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_GenericClient] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_GenericCopy] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_PY_MACRO2D] = GMSH_LINE;
    gmshtypes[ELE_TAG_SimpleContact2D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SimpleContact3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_BeamContact3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SurfaceLoad] = GMSH_QUAD;
    gmshtypes[ELE_TAG_BeamContact2D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_BeamEndContact3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SSPquad] = GMSH_QUAD;
    gmshtypes[ELE_TAG_SSPquadUP] = GMSH_QUAD;
    gmshtypes[ELE_TAG_SSPbrick] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_SSPbrickUP] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_BeamContact2Dp] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_BeamContact3Dp] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_BeamEndContact3Dp] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_Quad4FiberOverlay] = GMSH_QUAD;
    gmshtypes[ELE_TAG_Brick8FiberOverlay] = GMSH_HEXAHEDRON;
    gmshtypes[ELE_TAG_QuadBeamEmbedContact] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_DispBeamColumn2dThermal] = GMSH_LINE;
    gmshtypes[ELE_TAG_TPB1D] = GMSH_LINE;
    gmshtypes[ELE_TAG_TFP_Bearing] = GMSH_LINE;
    gmshtypes[ELE_TAG_TFP_Bearing2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_TripleFrictionPendulum] = GMSH_LINE;
    gmshtypes[ELE_TAG_TripleFrictionPendulumX] = GMSH_LINE;
    gmshtypes[ELE_TAG_PFEMElement2D] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_FourNodeQuad02] = GMSH_QUAD;
    gmshtypes[ELE_TAG_cont2d01] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_cont2d02] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_CST] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_Truss2] = GMSH_LINE;
    gmshtypes[ELE_TAG_CorotTruss2] = GMSH_LINE;
    gmshtypes[ELE_Tag_ZeroLengthImpact3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_PFEMElement3D] = GMSH_TETRA;
    gmshtypes[ELE_TAG_PFEMElement2DCompressible] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_PFEMElement2DBubble] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_PFEMElement2Dmini] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_ElasticTimoshenkoBeam2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElasticTimoshenkoBeam3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElastomericBearingUFRP2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElastomericBearingUFRP3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_RJWatsonEQS2d] = GMSH_LINE;
    gmshtypes[ELE_TAG_RJWatsonEQS3d] = GMSH_LINE;
    gmshtypes[ELE_TAG_HDR] = GMSH_LINE;
    gmshtypes[ELE_TAG_ElastomericX] = GMSH_LINE;
    gmshtypes[ELE_TAG_LeadRubberX] = GMSH_LINE;
    gmshtypes[ELE_TAG_PileToe3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_N4BiaxialTruss] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_ShellDKGQ] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ShellNLDKGQ] = GMSH_QUAD;
    gmshtypes[ELE_TAG_MultipleShearSpring] = GMSH_LINE;
    gmshtypes[ELE_TAG_MultipleNormalSpring] = GMSH_LINE;
    gmshtypes[ELE_TAG_KikuchiBearing] = GMSH_LINE;
    gmshtypes[ELE_TAG_YamamotoBiaxialHDR] = GMSH_LINE;
    gmshtypes[ELE_TAG_MVLEM] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SFI_MVLEM] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_MVLEM_3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_SFI_MVLEM_3D] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_E_SFI_MVLEM_3D] = GMSH_POLY_VERTEX;
	gmshtypes[ELE_TAG_E_SFI] = GMSH_POLY_VERTEX;
	gmshtypes[ELE_TAG_MEFI] = GMSH_POLY_VERTEX;
    gmshtypes[ELE_TAG_PFEMElement2DFIC] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_CatenaryCable] = GMSH_LINE;
    gmshtypes[ELE_TAG_FourNodeTetrahedron] = GMSH_TETRA;
    gmshtypes[ELE_TAG_TriSurfaceLoad] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_ShellANDeS] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_ShellDKGT] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_ShellNLDKGT] = GMSH_TRIANGLE;
    gmshtypes[ELE_TAG_InertiaTruss] = GMSH_LINE;
    gmshtypes[ELE_TAG_ASDAbsorbingBoundary2D] = GMSH_QUAD;
    gmshtypes[ELE_TAG_ASDAbsorbingBoundary3D] = GMSH_HEXAHEDRON;
}

int GmshRecorder::flush() {
    if (theFile.is_open() && theFile.good()) {
        theFile.flush();
    }
    return 0;
}
