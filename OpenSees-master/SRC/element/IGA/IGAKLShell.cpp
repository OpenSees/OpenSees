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

// $Revision: 1.10 $
// $Date: 2011/03/10 22:51:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/IGAKLShell.cpp,v $

// Original implementation: Ed "C++" Love
// Reimplementation: Leopoldo Tesser, Diego A. Talledo, Véronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//STD library
#include <float.h>
#include <vector>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <IGAKLShell.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <MaterialResponse.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

#define min(a,b) ( (a)<(b) ? (a):(b) )

#include "gaussQuadrature.h"
#include "R3vectors.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// for debugging
#include <signal.h>

#include <set>
using namespace std;
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
  std::vector<T> values;
  for (T value = start; value < stop; value += step)
    values.push_back(value);
  return values;
}

static int numIGAKLShell = 0;


//static data
Matrix*  IGAKLShell::stiff = 0;
Matrix*  IGAKLShell::mass = 0;
Vector*  IGAKLShell::resid = 0;
// Vector*  IGAKLShell::load = 0;




//null constructor
IGAKLShell::IGAKLShell( ) :
  Element( 0, ELE_TAG_IGAKLShell ),
  connectedExternalNodes(4)
{

}


//*********************************************************************
//full constructor
IGAKLShell::IGAKLShell( int tag,
                        IGASurfacePatch *myPatch_,
                        const ID& nodes,
                        int ngauss_,
                        const Vector& xiE_,
                        const Vector& etaE_,
                        const ID& matTags):
  Element( tag, ELE_TAG_IGAKLShell ),
  ngauss(ngauss_),
  myPatch(myPatch_),
  xiE(xiE_),
  etaE(etaE_),
  connectedExternalNodes(nodes)
{
  if (numIGAKLShell == 0) {
    // opserr << "Using IGAKLShell - Developed by: Felipe Elgueta and Jose A. Abell (www.joseabell.com)\n";
    numIGAKLShell++;
  }
  // ngauss  = quadorder * quadorder;
  nLayers = myPatch->getNLayers();

  // ngauss=36;

  quadPoint = new Matrix(ngauss, 2);
  quadWeight = new Vector(ngauss);

  ID PQ = myPatch -> getOrders();
  int P = PQ(0);
  int Q = PQ(1);

  // P=2;
  // Q=2;

  gaussQuad2dNurbs(P + 1, Q + 1, quadPoint, quadWeight);

  // opserr << "quadPoint" << (*quadPoint) << endln;
  // opserr << "quadWeight" << (*quadWeight) << endln;
  // raise(SIGSEGV);

  // gaussQuad2dNurbs(P + 1, P + 1, quadPoint, quadWeight);
  // gaussQuad2dNurbs(6, 6, quadPoint, quadWeight);



  materialPointers = new NDMaterial** [ngauss]; //Double array
  for (int i = 0; i < ngauss; ++i)
  {
    materialPointers[i] = new NDMaterial* [nLayers];
  }


  for (int gp = 0 ;  gp < ngauss; gp++ )
  {
    for (int capa = 0; capa < nLayers; capa++)
    {
      NDMaterial* theReferenceMaterial = OPS_getNDMaterial(myPatch->getMatTag(capa)); // Pointer to NDMaterial
      NDMaterial* newmat = theReferenceMaterial->getCopy( );  // Copy of pointer to NDMaterial
      materialPointers[gp][capa] =  newmat;   // llama new

      if (materialPointers[gp][capa] == 0) {
        opserr << "ShellMITC4::constructor - failed to get a material of type: ShellSection\n";
      } //end if
    }
  } //end for gp

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  pressure = 0;

  load = 0; // null pointer

}

//******************************************************************

//destructor
IGAKLShell::~IGAKLShell( )
{

  // opserr << 'IGAKLShell::~IGAKLShell' << endln;
  // borrar todos los punteros. ARREGLAR PARA ARREGLO 2D
  int nLayers = myPatch->getNLayers();
  // opserr << "quadorder = " << quadorder << endln;
  // int i ;
  // for ( i = 0 ;  i < 4; i++ ) {

  //   delete materialPointers[i] ;
  //   materialPointers[i] = 0 ;

  // } //end for i
  for (int gp = 0; gp < ngauss; ++gp)
  {
    for (int capa = 0; capa < nLayers; ++capa)
    {
      delete materialPointers[gp][capa];
      materialPointers[gp][capa] = 0;
    }

    nodePointers[gp] = 0;

  }
  // if (Ki != 0)
  //   delete Ki;


  if (load != 0){
    delete load;
    load = 0;
  }

  // if (Ki != 0)
  //   delete Ki;

}
//**************************************************************************


//set domain
void  IGAKLShell::setDomain( Domain *theDomain )
{
  int i;

  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;

  // Initialiaze class-wide static data
  if (numIGAKLShell == 1)
  {
    stiff = new Matrix(NDOF, NDOF);
    mass = new Matrix(NDOF, NDOF);
    resid = new Vector(NDOF);
    // load = new Vector(NDOF);
  }

  //node pointers
  nodePointers = new Node* [connectedExternalNodes.Size()];

  for ( i = 0; i < connectedExternalNodes.Size(); i++ ) {
    nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
    if (nodePointers[i] == 0) {
      opserr << "IGAKLShell::setDomain - no node " << connectedExternalNodes(i);
      opserr << " exists in the model\n";
    }
  }

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  IGAKLShell::getNumExternalNodes( ) const
{
  return connectedExternalNodes.Size() ;
}


//return connected external nodes
const ID&  IGAKLShell::getExternalNodes( )
{
  return connectedExternalNodes ;
}


Node **
IGAKLShell::getNodePtrs(void)
{
  return nodePointers;
}

//return number of dofs
int  IGAKLShell::getNumDOF( )
{
  return 3 * connectedExternalNodes.Size() ;
}

//commit state
int  IGAKLShell::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "IGAKLShell::commitState () - failed in base class";
  }

  // bool nonLinearGeometry = myPatch->getAnalysisType();
  // if (nonLinearGeometry == false)
  // {
  //   this-> formStrainsLinear();
  // }

  // opserr << "IGAKLShell::commitState called!" << endln;
  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {

      // Should update the stresses here ?


      success += materialPointers[gp][capa]->commitState( ) ;
    }

  return success ;
}



//revert to last commit
int  IGAKLShell::revertToLastCommit( )
{
  int success = 0 ;

  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {
      success += materialPointers[gp][capa]->revertToLastCommit( ) ;
    }

  return success ;
}


//revert to start
int  IGAKLShell::revertToStart( )
{
  int success = 0 ;

  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {
      success += materialPointers[gp][capa]->revertToStart( ) ;
    }

  return success ;
}

//print out element data
void  IGAKLShell::Print( OPS_Stream &s, int flag )
{

  // s << "IGAKLShell # " << this->getTag() << " \n" 
  // << "    quadorder = " << quadorder << " \n" 
  // << "    ngauss = " << ngauss << " \n" 
  // << "    nLayers = " << nLayers << " \n" << endln; 

}

#define RESPONSETYPE_KLSHELL_FORCES 1
#define RESPONSETYPE_KLSHELL_STRESSES 2
#define RESPONSETYPE_KLSHELL_STRAINS 3
#define RESPONSETYPE_KLSHELL_FAKESECTION_DEFORMATIONS 101
#define RESPONSETYPE_KLSHELL_FAKESECTION_FORCES 102
#define RESPONSETYPE_KLSHELL_FAKESECTION_FORCEANDDEFORMATION 103
#define RESPONSETYPE_KLSHELL_FAKESECTION_FIBER 104

//For MPCO recorder:
// "IGAGauss1E" xi coordinates gauss [Vectror od doubles) all
// "IGAGauss2E" eta coord gauss [Vectror od doubles) volums ans surfaces
// "IGAGauss3E" mu coord gauss [Vectror od doubles) only pr volumes
// "IGAGauss1P" xi coordinates gauss [Vectror od doubles) all
// "IGAGauss2P" eta coord gauss [Vectror od doubles) volums ans surfaces
// "IGAGauss3P" mu coord gauss [Vectror od doubles) only pr volumes
// "IGAGaussWeight" (vector) for all

#define RESPONSETYPE_KLSHELL_IGAGauss1E 200
#define RESPONSETYPE_KLSHELL_IGAGauss2E 201
#define RESPONSETYPE_KLSHELL_IGAGauss3E 202
#define RESPONSETYPE_KLSHELL_IGAGauss1P 203
#define RESPONSETYPE_KLSHELL_IGAGauss2P 204
#define RESPONSETYPE_KLSHELL_IGAGauss3P 205
#define RESPONSETYPE_KLSHELL_IGAGaussWeight 206



Response*
IGAKLShell::setResponse(const char **argv, int argc, OPS_Stream &output)
{

  Response* theResponse = 0;

  // opserr << "IGAKLShell::setResponse - start argv list - argc = "  << argc << endln;
  // for (int i = 0; i < argc; ++i)
  // {
  //   opserr << "argv[" << i << "] = " << argv[i] << endln;
  // }
  // opserr << "IGAKLShell::setResponse - end argv list - argc = "  << argc << endln;

  // 1 - OPEN THE ELEMENT TAG
  output.tag("ElementOutput");
  output.attr("eleType", "IGAKLShell");
  output.attr("eleTag", this->getTag());
  int numNodes = this->getNumExternalNodes();
  const ID& nodes = this->getExternalNodes();
  static char nodeData[32];

  for (int i = 0; i < numNodes; i++) {
      sprintf(nodeData, "node%d", i + 1);
      output.attr(nodeData, nodes(i));
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ){
      //Not yet implemented
  }
  else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0) {
      if (argc < 2) {
          opserr << "KLShell::setResponse() - need to specify more data\n";
          return 0;
      }
      int pointNum = atoi(argv[1]);
      // opserr << "IGAKLShell::setResponse - requesting material at GP # " << pointNum << " argv[0] = " << argv[0] << " argv[1] = " << argv[1] << endln;
      if (pointNum > 0 && pointNum <= ngauss) {

          // 2 - OPEN THE GAUSS TAG
          // Pipe: assuming pointNum is the gauss point gp [1,ngauss]
          double wt = (*quadWeight)(pointNum-1);
          double ptU = (*quadPoint)(pointNum-1, 0);
          double ptV = (*quadPoint)(pointNum-1, 1);

          double xi = myPatch->parent2ParametricSpace(xiE, ptU);
          double eta = myPatch->parent2ParametricSpace(etaE, ptV);
          // double xi = xi_es_algo_que_el_pipe_tiene_que_calcular;
          // double eta = eta_es_algo_que_el_pipe_tiene_que_calcular;

          output.tag("GaussPoint");
          output.attr("number", pointNum);
          output.attr("eta", xi);
          output.attr("neta", eta);

          // REPLACE THIS LINE WITH AN EMULATION OF WHAT A LayeredShellFiberSection would do
          //theResponse = m_sections[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
          theResponse = this->emulateSectionSetResponse(&argv[2], argc - 2, output, pointNum, xi, eta);

          // 2* - CLOSE THE GAUSS TAG
          output.endTag();
      }
  }
  else if (strcmp(argv[0], "stresses") == 0) {
      
    for (int i = 0; i < ngauss; i++) {

      double wt = (*quadWeight)(i);
      double ptU = (*quadPoint)(i, 0); //[-1 , 1]
      double ptV = (*quadPoint)(i, 1); //[-1 , 1]

      double xi = myPatch->parent2ParametricSpace(xiE, ptU); //[min, max]
      double eta = myPatch->parent2ParametricSpace(etaE, ptV); //[min, max]
      // double xi = xi_es_algo_que_el_pipe_tiene_que_calcular;
      // double eta = eta_es_algo_que_el_pipe_tiene_que_calcular;


      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", xi);
      output.attr("neta", eta);

      output.tag("SectionOutput");
      output.attr("classType", this->getClassTag());
      output.attr("tag", this->getTag());

      output.tag("ResponseType", "p11");
      output.tag("ResponseType", "p22");
      output.tag("ResponseType", "p1212");
      output.tag("ResponseType", "m11");
      output.tag("ResponseType", "m22");
      output.tag("ResponseType", "m12");
      // output.tag("ResponseType", "q1");
      // output.tag("ResponseType", "q2");

      output.endTag(); // SectionOutput
      output.endTag(); // GaussPoint
    }

    const int stress_components_per_gauss_point = 6;
    int number_of_output_components = stress_components_per_gauss_point*ngauss;
    theResponse =  new ElementResponse(this, RESPONSETYPE_KLSHELL_STRESSES, Vector(number_of_output_components));
  }
  else if (strcmp(argv[0], "strains") == 0) {
      
    for (int i = 0; i < ngauss; i++) {

      double wt = (*quadWeight)(i);
      double ptU = (*quadPoint)(i, 0);
      double ptV = (*quadPoint)(i, 1);

      double xi = myPatch->parent2ParametricSpace(xiE, ptU);
      double eta = myPatch->parent2ParametricSpace(etaE, ptV);
      // double xi = xi_es_algo_que_el_pipe_tiene_que_calcular;
      // double eta = eta_es_algo_que_el_pipe_tiene_que_calcular;

      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", xi);
      output.attr("neta", eta);

      output.tag("SectionForceDeformation");
      output.attr("classType", this->getClassTag());
      output.attr("tag", this->getTag());

      output.tag("ResponseType", "e11");
      output.tag("ResponseType", "e22");
      output.tag("ResponseType", "e1212");
      output.tag("ResponseType", "k11");
      output.tag("ResponseType", "k22");
      output.tag("ResponseType", "k12");
      // output.tag("ResponseType", "q1");
      // output.tag("ResponseType", "q2");

      output.endTag(); // GaussPoint
      output.endTag(); // NdMaterialOutput
    }

    const int strain_components_per_gauss_point = 6;
    int number_of_output_components = strain_components_per_gauss_point*ngauss;
    theResponse =  new ElementResponse(this, RESPONSETYPE_KLSHELL_STRAINS, Vector(number_of_output_components));
  }
  else if (strcmp(argv[0], "material_and_layer_number") == 0 ) {
    if (argc < 3) {
      opserr << "IGAKLShell::setResponse() - need to specify more data\n";
      return 0;
    }
    
    int gp = atoi(argv[1]);
    int capa = atoi(argv[2]);
    // opserr << "IGAKLShell::setResponse eleTag = "<< this->getTag() << "called" << endln;

    if (gp >= 0 && gp < ngauss) { // was gp>0 and gp <= ngauss


      output.tag("GaussPoint");
      output.attr("number", gp);
      // output.attr("eta", sg[pointNum - 1]);
      // output.attr("neta", tg[pointNum - 1]);
      output.attr("layer", capa);

      theResponse =  materialPointers[gp][capa]->setResponse(&argv[3], argc - 3, output);

      output.endTag();
    }
  } else if(strcmp(argv[0], "IGAGauss1E") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss1E, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGauss2E") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss2E, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGauss3E") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss3E, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGauss1P") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss1P, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGauss2P") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss2P, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGauss3P") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGauss3P, Vector(ngauss));
  } else if(strcmp(argv[0], "IGAGaussWeight") == 0) {
        theResponse =  new ElementResponse(this, 
          RESPONSETYPE_KLSHELL_IGAGaussWeight, Vector(ngauss));
  }


  // 1* - CLOSE THE ELEMENT TAG
  output.endTag();
  return theResponse;

  // End of Massimo's code

}


// // Namespace added from Massimo's code
namespace {
    
    // A simple Singleton class used to generate unique IDs from
    // ordered integer sequences (not optimized... not necessary)
    class fakeSectionIdBuilder {
    private:
        int current_id = 0;
        std::map<std::vector<int>, int> ids;
        
    private:
        fakeSectionIdBuilder() = default;
        fakeSectionIdBuilder(const fakeSectionIdBuilder &) = delete;
        fakeSectionIdBuilder& operator = (const fakeSectionIdBuilder &) = delete;
        
    public:
        static fakeSectionIdBuilder& instance() {
            static fakeSectionIdBuilder _instance;
            return _instance;
        }
        
        int getId(const std::vector<int> &key) {
            auto it = ids.find(key);
            if (it == ids.end()) {
                ++current_id;
                ids[key] = current_id;
                return current_id;
            }
            else {
                return it->second;
            }
        }
    };
    
}

// THIS IS A COPY-PASTE OF THE LayeredShellFiberSection setResponse method.
// For the sake of readability, you can create a private method in your element
Response*
IGAKLShell::emulateSectionSetResponse(const char **argv, int argc,
                                      OPS_Stream &output, int gaussPointNum, double xi, double eta)
{
  Response *theResponse =0;

  static Vector vectorResponse(6);
  vectorResponse.Zero();

  // // opserr << "IGAKLShell::emulateSectionSetResponse - start argv list - argc = "  << argc << endln;
  // for (int i = 0; i < argc; ++i)
  // {
  //   opserr << "argv[" << i << "] = " << argv[i] << endln;
  // }
  // opserr << "IGAKLShell::emulateSectionSetResponse - end argv list - argc = "  << argc << endln;

  if (argc == 0)
  {
    return 0; 
  }

  output.tag("SectionForceDeformation");
  
  // You don't have these methods.. they belong to SectionForceDeformation.
  // the getClassType is not an issue... you can hardcode a name for it.
  // for the getTag.. you can either use an hard coded integer, or you can make it
  // smarter, generating a unique ID depending on integer sequences, that could be 
  // your sequence of materials and thicknesses (somehow encoded into integers...).
  // in this way all section with the same Layout (material and thickness stack order)
  // will have the same ID...
  
  // output.attr("secType", this->getClassType());
  // output.attr("secTag", this->getTag());
  
  output.attr("secType", "KLShellBuiltinSection");
  
  // std::vector<int> encoded;
  // for(int i = 0; i < nLayers; i++) {
  //   encoded.push_back(theMaterialIds[i]);
  //   encoded.push_back(float_to_int(theThicknesses[i])); // encode a float into an int the way you like
  //   encoded.push_back(float_to_int(theThetas[i])); // encode a float into an int the way you like encoded
  // }
  // output.attr("secTag", fakeSectionIdBuilder::instance().getId(encoded));
  output.attr("secTag", myPatch->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    output.tag("ResponseType","e11");
    output.tag("ResponseType","e22");
    output.tag("ResponseType","e1212");
    output.tag("ResponseType","k11");
    output.tag("ResponseType","k22");
    output.tag("ResponseType","k12");

    vectorResponse.resize(6);
    // theResponse =  new MaterialResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_DEFORMATIONS, vectorResponse);
    // Switching this type of response to elemental response... will this be a problem? For massimo. 
    theResponse =  new ElementResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_DEFORMATIONS, vectorResponse);
    
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    output.tag("ResponseType","p11");
    output.tag("ResponseType","p22");
    output.tag("ResponseType","p1212");
    output.tag("ResponseType","m11");
    output.tag("ResponseType","m22");
    output.tag("ResponseType","m12");
    vectorResponse.resize(6);
    // theResponse =  new MaterialResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_FORCES, vectorResponse);
    // Switching this type of response to elemental response... will this be a problem? For massimo. 
    theResponse =  new ElementResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_FORCES, vectorResponse);
  

  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    output.tag("ResponseType","e11");
    output.tag("ResponseType","e22");
    output.tag("ResponseType","e1212");
    output.tag("ResponseType","k11");
    output.tag("ResponseType","k22");
    output.tag("ResponseType","k12");
    output.tag("ResponseType","p11");
    output.tag("ResponseType","p22");
    output.tag("ResponseType","p1212");
    output.tag("ResponseType","m11");
    output.tag("ResponseType","m22");
    output.tag("ResponseType","m12");

    vectorResponse.resize(12);
    // theResponse =  new MaterialResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_FORCEANDDEFORMATION, vectorResponse);
    // Switching this type of response to elemental response... will this be a problem? For massimo. 
    theResponse =  new ElementResponse(this, RESPONSETYPE_KLSHELL_FAKESECTION_FORCEANDDEFORMATION, vectorResponse);
  }  
  else if (strcmp(argv[0],"fiber") == 0 || strcmp(argv[0],"Fiber") == 0) {
    if (argc < 3) {
      opserr << "LayeredShellFiberSection::setResponse() - need to specify more data\n";
      return 0;
    }
    int capa = atoi(argv[1]);

    int nCapas = myPatch -> getNLayers();

    if (capa > 0 && capa <= nCapas) {
      
      // double zLoc = Pipe
      double zLoc = myPatch->getZk(capa-1);
      // double thickness = PipePipe
      double thickness = myPatch->getThickness(capa-1);

      output.tag("FiberOutput");
      output.attr("number",capa);
      output.attr("zLoc",zLoc);
      output.attr("thickness",thickness);
      
      // theResponse =  theFibers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      // RESPONSETYPE_KLSHELL_FAKESECTION_FIBER --> in this case it is set by the fiber setResponse
      // opserr << "IGAKLShell::emulateSectionSetResponse called with argv[2] " << argv[2] << endln;
      theResponse =  materialPointers[gaussPointNum-1][capa-1]->setResponse(&argv[2], argc - 2, output);
      
      output.endTag();
    }
    else
    {
      theResponse = 0;
    }
  }
  output.endTag(); // SectionOutput
  return theResponse;
}



int
IGAKLShell::getResponse(int responseID, Information &eleInfo)
{
  int cnt = 0;
  const int stress_components_per_gauss_point = 6;
  int number_of_output_components = stress_components_per_gauss_point*ngauss;

  Vector stresses(number_of_output_components);
  Vector strains(number_of_output_components);

  // switch (responseID) {
  // case RESPONSETYPE_KLSHELL_FORCES: // global forces
  if (responseID == RESPONSETYPE_KLSHELL_FORCES || 
    responseID == RESPONSETYPE_KLSHELL_FAKESECTION_FORCES)
  {
    return eleInfo.setVector(this->getResistingForce());
  }
  // break;

  // case RESPONSETYPE_KLSHELL_STRESSES: // stresses
  if( responseID == RESPONSETYPE_KLSHELL_STRESSES || 
    responseID == RESPONSETYPE_KLSHELL_FAKESECTION_FORCES)
  {

    // Number of shape functions
    int noFuncs = myPatch->getNoFuncs();

    static Matrix dR(noFuncs, 3);
    static Matrix ddR(noFuncs, 3);

    static Matrix dr(3, noFuncs);
    static Matrix ddr(3, noFuncs);

    // First computing reference and deformed configuration

    bool nonLinearGeometry = myPatch->getAnalysisType();

    Vector point(3); // Vector for storing node coordinates
    Vector point_disp(3); // Vector for storing node displacement

    Matrix pts(connectedExternalNodes.Size(), 3); // Matrix for storing element nodes coordinates, noFuncs,3
    Matrix pts_d(connectedExternalNodes.Size(), 3); // Matrix for storing element nodes displaced coordinates

    // Get pts span matrix for jacobian
    for (int i = 0; i < connectedExternalNodes.Size(); ++i)
    {
      point = nodePointers[i]->getCrds();
      point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp

      // opserr << "   @ i = " << i << " point_disp = " << point_disp << endln;

      for (int j = 0; j < 3; ++j)
      {
        pts(i, j) = point(j);
        if (isnan(point_disp(j)))
        {
          opserr << "Nan found on formResidAndTangent = " << endln;
        }
        if (nonLinearGeometry)
        {
          pts_d(i, j) = pts(i, j) + point_disp(j);
        }
        else
        {
          pts_d(i, j) = pts(i, j);
        }
      }
    }

    for (int gp = 0; gp < ngauss; gp++) { // changed "i" for "gp"

      // opserr << "IGAKLShell::getResponse() - tag = " << this->getTag() << " computing stresses at gp = " << gp << endln;

      // Get material stress response
      // const Vector &sigma = materialPointers[i]->getStress();
      static Vector N_ca(3); 
      static Vector M_ca(3);


      // Computing curvature

      Vector R(noFuncs);
      Vector dRdxi(noFuncs);
      Vector dRdeta(noFuncs);
      Vector dR2dxi(noFuncs);
      Vector dR2deta(noFuncs);
      Vector dR2dxideta(noFuncs);

      R.Zero();
      dRdxi.Zero();
      dRdeta.Zero();
      dR2dxi.Zero();
      dR2deta.Zero();
      dR2dxideta.Zero();

      double wt = (*quadWeight)(gp);
      double ptU = (*quadPoint)(gp, 0);
      double ptV = (*quadPoint)(gp, 1);

      double xi = myPatch->parent2ParametricSpace(xiE, ptU);
      double eta = myPatch->parent2ParametricSpace(etaE, ptV);

      myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

      // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
      Matrix ddR(3, noFuncs);
      Matrix dR(2, noFuncs);
      // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
      for (int i = 0; i < noFuncs; ++i)
      {
        dr(0, i) = dRdxi(i);
        dr(1, i) = dRdeta(i);
        dr(2, i) = 0;
        ddr(0, i) = dR2dxi(i);
        ddr(1, i) = dR2deta(i);
        ddr(2, i) = dR2dxideta(i);
      }

      Matrix G(2, 3);
      G = dR  * pts;
      Matrix H(3, 3);
      H = ddR * pts;

      Matrix g(2, 3);
      g = dR  * pts_d;
      Matrix h(3, 3);
      h = ddR * pts_d;

      // Transposing because it is needed in cols drdxi, drdeta, etc
      G = transpose(2, 3, G);
      g = transpose(2, 3, g);

      H = transpose(3, 3, H);
      h = transpose(3, 3, h);

      // Shell geo temporal parameters
      Vector g3_tmp(3);
      Vector n_tmp(3);
      double lg3_tmp = 0; // = lG3
      Matrix T_E_G_tmp(3, 3); // = T_E_G
      Matrix T_G_E_tmp(3, 3); // = T_G_E
      Matrix Tb_tmp(3, 3); // = T_Gcon_E

      Matrix Gab_r(2, 2);// = Gab
      Gab_r.Zero();
      Matrix gab(2, 2); // = Gab
      gab.Zero();

      Vector Bv_r(3); // = Bv
      Bv_r.Zero();
      Vector bv(3); // = Bv
      bv.Zero();

      Matrix Tb(3, 3); // = T_Gcon_E
      Tb.Zero();

      // Call shell geometry functions
      shellGeo(G , H , g3_tmp , lg3_tmp  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);  // Reference configuration
      shellGeo(g , h , g3_tmp     , lg3_tmp , n_tmp     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp); // Deformed or actual configuration

      //
      // Compute the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
      //

      Vector K_cu(3);
      K_cu.Zero();
      Vector K_ca(3);
      K_ca.Zero();

      //curvature vector [K11,K22,K12] referred to curvilinear coor sys
      K_cu(0) = -(bv(0) - Bv_r(0));
      K_cu(1) = -(bv(1) - Bv_r(1));
      K_cu(2) = -(bv(2) - Bv_r(2));


      K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys, , Tb is T_Gcon_E




      // Looping through the layers

      double theta1 = 0; // In-plane rotation, assuming zero in the meantime
      Matrix T(3, 3);   // Rotating laminate matrix

      // opserr << "   @ gp = " << gp << " forming stresses with nLayers = " << nLayers << endln;
      for (int iPly = 0; iPly < nLayers; ++iPly)
      {
        const Matrix& Czig = materialPointers[gp][iPly]->getTangent();
        Vector Szeta = materialPointers[gp][iPly]->getStress();  // necesitas para plasticidad
        double iThickness = myPatch->getThickness(iPly);
        double iZ         = myPatch->getZk(iPly);// Mid center laminate position
        double iAngle     = myPatch->getAngle(iPly) - theta1;

        // opserr << "          iPly = " << iPly << " Szeta = " << Szeta << endln;

        T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
        T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
        T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;


        Matrix C(3, 3);
        C.addMatrixTripleProduct(0, T, Czig, 1.0);


        Szeta = transpose(3, 3, T) * Szeta; // Rotating stresses
        // Szeta = T * Szeta; // Rotating stresses

        // Integrating stres in every ply
        N_ca += iThickness * Szeta;
        M_ca += iThickness * iZ * Szeta + (C * K_ca) * pow(iThickness, 3) / 12.0 ;
      }

      // opserr << "   @ gp = " << gp << " M_ca = " << M_ca << endln;

      // N_ca = pipe_calcula_N_ca_para_este_punto_de_gauss();
      // M_ca = pipe_calcula_M_ca_para_este_punto_de_gauss();
      stresses(cnt) = N_ca(0);
      stresses(cnt + 1) = N_ca(1);
      stresses(cnt + 2) = N_ca(2);
      stresses(cnt + 3) = M_ca(0);
      stresses(cnt + 4) = M_ca(1);
      stresses(cnt + 5) = M_ca(2);
      cnt += 6;
    }

    return eleInfo.setVector(stresses);
  }
    // break;

  // case RESPONSETYPE_KLSHELL_STRAINS: //strain
  if(responseID == RESPONSETYPE_KLSHELL_STRAINS || 
    responseID == RESPONSETYPE_KLSHELL_FAKESECTION_DEFORMATIONS)
  {
    for (int gp = 0; gp < ngauss; gp++) { // Changed "i" to "gp"

      // Get section deformation
      static Vector E_ca(3); 
      static Vector K_ca(3);

      // Computing curvature


      // First computing reference and deformed configuration

      bool nonLinearGeometry = myPatch->getAnalysisType();

      Vector point(3); // Vector for storing node coordinates
      Vector point_disp(3); // Vector for storing node displacement

      Matrix pts(connectedExternalNodes.Size(), 3); // Matrix for storing element nodes coordinates, noFuncs,3
      Matrix pts_d(connectedExternalNodes.Size(), 3); // Matrix for storing element nodes displaced coordinates

      // Get pts span matrix for jacobian
      for (int i = 0; i < connectedExternalNodes.Size(); ++i)
      {
        point = nodePointers[i]->getCrds();
        point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
        for (int j = 0; j < 3; ++j)
        {
          pts(i, j) = point(j);
          if (isnan(point_disp(j)))
          {
            opserr << "Nan found on formResidAndTangent = " << endln;
          }
          if (nonLinearGeometry)
          {
            pts_d(i, j) = pts(i, j) + point_disp(j);
          }
          else
          {
            pts_d(i, j) = pts(i, j);
          }
        }
      }

      // Number of shape functions
      int noFuncs = myPatch->getNoFuncs();

      Vector R(noFuncs);
      Vector dRdxi(noFuncs);
      Vector dRdeta(noFuncs);
      Vector dR2dxi(noFuncs);
      Vector dR2deta(noFuncs);
      Vector dR2dxideta(noFuncs);

      R.Zero();
      dRdxi.Zero();
      dRdeta.Zero();
      dR2dxi.Zero();
      dR2deta.Zero();
      dR2dxideta.Zero();

      double wt = (*quadWeight)(gp);
      double ptU = (*quadPoint)(gp, 0);
      double ptV = (*quadPoint)(gp, 1);

      double xi = myPatch->parent2ParametricSpace(xiE, ptU);
      double eta = myPatch->parent2ParametricSpace(etaE, ptV);

      myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

      // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
      Matrix ddR(3, noFuncs);
      Matrix dR(2, noFuncs);
      for (int i = 0; i < noFuncs; ++i)
      {
        dR(0,i) = dRdxi(i);
        dR(1,i) = dRdeta(i);
        ddR(0, i) = dR2dxi(i);
        ddR(1, i) = dR2deta(i);
        ddR(2, i) = dR2dxideta(i);
      }

      Matrix G(2, 3);
      G = dR  * pts;
      Matrix H(3, 3);
      H = ddR * pts;

      Matrix g(2, 3);
      g = dR  * pts_d;
      Matrix h(3, 3);
      h = ddR * pts_d;

      // Transposing because it is needed in cols drdxi, drdeta, etc
      G = transpose(2, 3, G);
      g = transpose(2, 3, g);

      H = transpose(3, 3, H);
      h = transpose(3, 3, h);

      // Shell geo temporal parameters
      Vector g3_tmp(3);
      Vector n_tmp(3);
      double lg3_tmp = 0; // = lG3
      Matrix T_E_G_tmp(3, 3); // = T_E_G
      Matrix T_G_E_tmp(3, 3); // = T_G_E
      Matrix Tb_tmp(3, 3); // = T_Gcon_E

      Matrix Gab_r(2, 2);// = Gab
      Gab_r.Zero();
      Matrix gab(2, 2); // = Gab
      gab.Zero();

      Vector Bv_r(3); // = Bv
      Bv_r.Zero();
      Vector bv(3); // = Bv
      bv.Zero();

      Matrix Tb(3, 3); // = T_Gcon_E
      Tb.Zero();

      // Call shell geometry functions
      shellGeo(G , H , g3_tmp , lg3_tmp  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);  // Reference configuration
      shellGeo(g , h , g3_tmp     , lg3_tmp , n_tmp     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp); // Deformed or actual configuration

      //
      // Compute the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
      //

      Vector K_cu(3);
      Vector E_cu(3);
      K_cu.Zero();
      E_cu.Zero();


      //strain vector [E11,E22,E12] referred to curvilinear coor sys
      E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
      E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
      E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

      //curvature vector [K11,K22,K12] referred to curvilinear coor sys
      K_cu(0) = -(bv(0) - Bv_r(0));
      K_cu(1) = -(bv(1) - Bv_r(1));
      K_cu(2) = -(bv(2) - Bv_r(2));


      E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys, Tb is T_Gcon_E
      K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys, , Tb is T_Gcon_E

      // E_ca = pipe_calcula_E_ca_para_este_punto_de_gauss();
      // K_ca = pipe_calcula_K_ca_para_este_punto_de_gauss();
      strains(cnt) = E_ca(0);
      strains(cnt + 1) = E_ca(1);
      strains(cnt + 2) = E_ca(2);
      strains(cnt + 3) = K_ca(0);
      strains(cnt + 4) = K_ca(1);
      strains(cnt + 5) = K_ca(2);
      cnt += 6;
    }
    return eleInfo.setVector(strains);
  }


  // At this point we are requesting data about the gauss points
  // Here we build that data:
  Vector IGAGauss1E(ngauss);
  Vector IGAGauss2E(ngauss);
  Vector IGAGauss3E(0); // IGAKLshell only has 2 dimensions
  Vector IGAGauss1P(ngauss);
  Vector IGAGauss2P(ngauss);
  Vector IGAGauss3P(0); // IGAKLshell only has 2 dimensions
  Vector IGAGaussWeight(ngauss);

  for (int i = 0; i < ngauss; i++) {
    // to-do : tidy up.. optimize...
    IGAGaussWeight(i) = (*quadWeight)(i);
    IGAGauss1E(i) = (*quadPoint)(i, 0); //[-1 , 1]
    IGAGauss2E(i) = (*quadPoint)(i, 1); //[-1 , 1]

    double IGAGauss1P = myPatch->parent2ParametricSpace(xiE, IGAGauss1E(i)); //[min, max]
    double IGAGauss2P = myPatch->parent2ParametricSpace(etaE, IGAGauss2E(i));
  }

  if(responseID == RESPONSETYPE_KLSHELL_IGAGauss1E){
    eleInfo.setVector(IGAGauss1E);
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGauss2E){
    eleInfo.setVector(IGAGauss2E);
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGauss3E){
    eleInfo.setVector(IGAGauss3E); // IGAKLshell only has 2 dimensions;
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGauss1P){
    eleInfo.setVector(IGAGauss1P);
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGauss2P){
    eleInfo.setVector(IGAGauss2P);
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGauss3P){
    eleInfo.setVector(IGAGauss3P); // IGAKLshell only has 2 dimensions;
  } else if(responseID == RESPONSETYPE_KLSHELL_IGAGaussWeight){
    eleInfo.setVector(IGAGaussWeight);
  } 


    return -1; // Failure
}


//return stiffness matrix
const Matrix&  IGAKLShell::getTangentStiff( )
{
  // opserr << "IGAKLShell::getTangentStiff - eleTag" << this->getTag() << " called! " << endln;

  // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;

  // opserr << "Element number:  = " << this->getTag() << endln;

  // opserr << "xiE = " << xiE << endln;
  // opserr << "etaE = " << etaE << endln;

  int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;
  // formTangentNguyen();

  // opserr << "(*stiff) = " << (*stiff) << endln;

  return *stiff ;



  // stiff->Zero();

  // bool nonLinearGeometry = myPatch->getAnalysisType();

  // // opserr << "xiE = " << xiE << endln;
  // // opserr << "etaE = " << etaE << endln;

  // float wt;
  // float ptU;
  // float ptV;
  // double xi;
  // double eta;

  // int noFuncs = myPatch->getNoFuncs();

  // Vector R(noFuncs);
  // Vector dRdxi(noFuncs);
  // Vector dRdeta(noFuncs);
  // Vector dR2dxi(noFuncs);
  // Vector dR2deta(noFuncs);
  // Vector dR2dxideta(noFuncs);

  // double J2; //Jacobian for the parametric mapping

  // Vector point(3);
  // Vector point_disp(3);
  // Matrix pts(connectedExternalNodes.Size(), 3);
  // Matrix pts_d(connectedExternalNodes.Size(), 3);
  // // Get pts span matrix for jacobian
  // for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  // {
  //   point = nodePointers[i]->getCrds();
  //   point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
  //   // points_disp = nodePointers[i]->getDisp(); //getTrialDisp
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     pts(i, j) = point(j);
  //     if (isnan(point_disp(j)))
  //     {
  //       opserr << "Nan found on getTangentStiff = " << endln;
  //     }
  //     if (nonLinearGeometry)
  //     {
  //       pts_d(i, j) = pts(i, j) + point_disp(j);
  //     }
  //     else
  //     {
  //       pts_d(i, j) = pts(i, j);
  //     }
  //   }
  // }
  // // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;
  // // opserr << "pts = " << pts << endln;

  // // opserr << "pts = " << pts << endln;

  // // Loop over integrations points
  // for (int gp = 0; gp < ngauss; ++gp)
  // {
  //   wt = (*quadWeight)(gp);
  //   ptU = (*quadPoint)(gp, 0);
  //   ptV = (*quadPoint)(gp, 1);
  //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
  //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
  //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));
  //   // opserr << "xi = " << xi << endln;
  //   // opserr << "eta = " << eta << endln;
  //   R.Zero();
  //   dRdxi.Zero();
  //   dRdeta.Zero();
  //   dR2dxi.Zero();
  //   dR2deta.Zero();
  //   dR2dxideta.Zero();
  //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

  //   // opserr << "R = " << R << endln;
  //   // opserr << "dRdxi = " << dRdxi << endln;
  //   // opserr << "dRdeta = " << dRdeta << endln;
  //   // opserr << "dR2dxi = " << dR2dxi << endln;
  //   // opserr << "dR2deta = " << dR2deta << endln;
  //   // opserr << "dR2dxideta = " << dR2dxideta << endln;




  //   // Get first and second order jacobians
  //   Matrix dr(2, noFuncs);
  //   Matrix dr2(3, noFuncs);


  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     dr(0, i) = dRdxi(i);
  //     dr(1, i) = dRdeta(i);

  //     dr2(0, i) = dR2dxi(i);
  //     dr2(1, i) = dR2deta(i);
  //     dr2(2, i) = dR2dxideta(i);
  //   }

  //   //compute the jacobian of physical and parameter domain mapping
  //   // then the derivative w.r.t spatial physical coordinates, current configuration


  //   Matrix jacob = dr * pts;  //G: Covariant base vectors, needed in the reference configuration
  //   Matrix jacob2 = dr2 * pts; //H: Hessian matrix, needed in the reference configuration


  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //   // !    DESCRIPTION OF THE VARIABLES !
  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  //   // ! g: matrix with the 3 components of the in plane vectors g1, g2
  //   // ! dR: first derivatives of the shape functions
  //   // ! ddR: second derivatives of the shape functions
  //   // ! dn: first variation of the normal vector g3
  //   // ! lg3: length of the normal vector g3

  //   // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
  //   // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
  //   // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
  //   // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

  //   // ! dE_cu: first variation of the membrane strain in local coord. system
  //   // ! dE_ca: first variation of the membrane strain in global coord. system
  //   // ! dK_cu: first variation of the curvature in local coord. system
  //   // ! dK_ca: first variation of the curvature in global coord. system

  //   // ! ddE_cu: second variation of the membrane strain in local coord. system
  //   // ! ddE_ca: second variation of the membrane strain in global coord. system
  //   // ! ddK_cu: second variation of the curvature in local coord. system
  //   // ! ddK_ca: second variation of the curvature in global coord. system

  //   Vector E_cu(3);
  //   Vector K_cu(3);

  //   Vector dE_cu(3);
  //   Matrix dE_ca(3, 3 * noFuncs);
  //   Vector dK_cu(3);
  //   Matrix dK_ca(3, 3 * noFuncs);

  //   Matrix ddR(3, noFuncs);
  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     ddR(0, i) = dR2dxi(i);
  //     ddR(1, i) = dR2deta(i);
  //     ddR(2, i) = dR2dxideta(i);
  //   }


  //   Matrix G = dr  * pts;
  //   Matrix H = ddR * pts;

  //   Matrix g = dr  * pts_d;
  //   Matrix h = ddR * pts_d;

  //   G = transpose(2, 3, G); // Transposing because it is needed in cols drdxi, drdeta
  //   g = transpose(2, 3, g);
  //   H = transpose(3, 3, H);
  //   h = transpose(3, 3, h);

  //   // Shell geo reference parameters and temporal
  //   Vector g3_tmp(3);
  //   Vector n_tmp(3);
  //   double dA; // = lG3
  //   Matrix Gab_r(2, 2); // = Gab
  //   Vector Bv_r(3); // = Bv
  //   Matrix Tb(3, 3); // = T_Gcon_E
  //   Matrix T_E_G_tmp(3, 3); // = T_Gcon_E
  //   Matrix T_G_E_tmp(3, 3); // = T_Gcon_E

  //   // Shell geo deformed parameters and temporal
  //   Vector g3(3); // = g3
  //   Vector n(3); // = n
  //   double lg3; // = dA
  //   Matrix gab(2, 2); // = Gab
  //   Vector bv(3); // = Bv
  //   Matrix Tb_tmp(3, 3);

  //   // Call shell geometry functions
  //   shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);
  //   shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp);

  //   //
  //   // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
  //   //

  //   //strain vector [E11,E22,E12] referred to curvilinear coor sys
  //   E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
  //   E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
  //   E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

  //   //curvature vector [K11,K22,K12] referred to curvilinear coor sys
  //   K_cu(0) = -(bv(0) - Bv_r(0));
  //   K_cu(1) = -(bv(1) - Bv_r(1));
  //   K_cu(2) = -(bv(2) - Bv_r(2));


  //   Vector E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys
  //   Vector K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys

  //   // opserr << "E_ca = " << E_ca << endln;
  //   // opserr << "K_ca = " << K_ca << endln;

  //   //
  //   // Compute the First variation of strain and curvature w.r.t the d.o.f.
  //   //

  //   // Local variables
  //   int kr = 0;
  //   int dirr = 0;
  //   Matrix dg(3, 2);
  //   Matrix dR = transpose(2, noFuncs, dr);
  //   ddR = transpose(3, noFuncs, ddR);
  //   Matrix dg3(3, 3 * noFuncs);
  //   Vector g3dg3(3 * noFuncs);
  //   Vector g3dg3lg3_3(3 * noFuncs);
  //   Matrix dn(3, 3 * noFuncs);

  //   // For "slicing" matrices
  //   ID ur_only(1);
  //   ID us_only(1);
  //   auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

  //   double lg3_3 = pow(lg3, 3);



  //   for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     // opserr << "kr = " << kr << endln;
  //     // opserr << "dirr = " << dirr << endln;

  //     dg(dirr, 0) = dR(kr, 0);
  //     dg(dirr, 1) = dR(kr, 1);

  //     // Strain
  //     dE_cu(0) = dR(kr, 0) * g(dirr, 0);
  //     dE_cu(1) = dR(kr, 1) * g(dirr, 1);
  //     dE_cu(2) = 0.5 * (dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1));

  //     dE_ca(0, ur) = (Tb * dE_cu)(0);
  //     dE_ca(1, ur) = (Tb * dE_cu)(1);
  //     dE_ca(2, ur) = (Tb * dE_cu)(2);


  //     // Curvature
  //     dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
  //     dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
  //     dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);

  //     g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
  //     g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;


  //     dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
  //     dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
  //     dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);


  //     dK_cu(0) = -(ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     dK_cu(1) = -(ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     dK_cu(2) = -(ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // CHECK THIS OUT, I HAD TO SWAP THE INDICES FOR H

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(0, 1) * dn(1, ur) + h(0, 2) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(1, 0) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(1, 2) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(2, 0) * dn(0, ur) + h(2, 1) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     dK_ca(0, ur) = (Tb * dK_cu)(0);
  //     dK_ca(1, ur) = (Tb * dK_cu)(1);
  //     dK_ca(2, ur) = (Tb * dK_cu)(2);


  //   }


  //   // Second variation of strain and curvature w.r.t. dofs
  //   Vector ddE_cu(3);
  //   Matrix ddE_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declarin a 3D array, instead using 3 matrices
  //   Matrix ddE_ca_1(3 * noFuncs, 3 * noFuncs);
  //   Matrix ddE_ca_2(3 * noFuncs, 3 * noFuncs);

  //   Vector ddK_cu(3);
  //   Matrix ddK_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declarin a 3D array, instead using 3 matrices
  //   Matrix ddK_ca_1(3 * noFuncs, 3 * noFuncs);
  //   Matrix ddK_ca_2(3 * noFuncs, 3 * noFuncs);

  //   // Local variables
  //   int ks = 0;
  //   int dirs = 0;
  //   int dirt = 0;
  //   int ddir = 0;
  //   Vector ddg3(3);
  //   double tmp1 = 0.0;
  //   double tmp2 = 0.0;
  //   double lg3_5 = pow(lg3, 5);
  //   Vector ddn(3);



  //   if (nonLinearGeometry) // true if non-linear geometry
  //   {
  //     // opserr << "Using nonLinearGeometry!" << endln;
  //     for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //     {
  //       // Local node number kr and dof direction dirr
  //       kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //       dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //       for (int us = 0; us <= ur; ++us)
  //       {
  //         ks = (us + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //         dirs = ur - 3 * (ks); // Takes value 0 1 2 0 1 2 ...

  //         // Strain
  //         if (dirr == dirs)
  //         {
  //           ddE_cu(0) = dR(kr, 0) * dR(ks, 0);
  //           ddE_cu(1) = dR(kr, 1) * dR(ks, 1);
  //           ddE_cu(2) = 0.5 * (dR(kr, 0) * dR(ks, 1) + dR(kr, 1) * dR(ks, 0));
  //         }

  //         ddE_ca_0(ur, us) = (Tb * ddE_cu)(0);
  //         ddE_ca_1(ur, us) = (Tb * ddE_cu)(1);
  //         ddE_ca_2(ur, us) = (Tb * ddE_cu)(2);

  //         // Curvature

  //         dirt = 6 - dirr - dirs; // Check this
  //         ddir = dirr - dirs;

  //         if (ddir == -1 or ddir == 2)
  //         {
  //           ddg3(dirt) = dR(kr, 0) * dR(ks, 1) - dR(ks, 0) * dR(kr, 1);
  //         }
  //         else if (ddir == 1 or ddir == -2)
  //         {
  //           ddg3(dirt) = -dR(kr, 0) * dR(ks, 1) + dR(ks, 0) * dR(kr, 1);
  //         }

  //         ur_only.fill(ur);
  //         us_only.fill(us);

  //         tmp1 = -(ddg3 ^ g3) + (transpose(3, 1, dg3(xyz, ur_only)) * dg3(xyz, us_only) / lg3_3)(0, 0); // the (0,0) is to transform the 1x1 matrix into double
  //         tmp2 = 3 * g3dg3(ur) * g3dg3(us) / lg3_5;

  //         ddn(0) = ddg3(0) / lg3 - g3dg3lg3_3(us) * dg3(0, ur) - g3dg3lg3_3(ur) * dg3(0, us) + tmp1 * g3(0) + tmp2 * g3(0);
  //         ddn(1) = ddg3(1) / lg3 - g3dg3lg3_3(us) * dg3(1, ur) - g3dg3lg3_3(ur) * dg3(1, us) + tmp1 * g3(1) + tmp2 * g3(1);
  //         ddn(2) = ddg3(2) / lg3 - g3dg3lg3_3(us) * dg3(2, ur) - g3dg3lg3_3(ur) * dg3(2, us) + tmp1 * g3(2) + tmp2 * g3(2);

  //         // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(0, 1) * ddn(1) + h(0, 2) * ddn(2)));
  //         // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(1, 0) * ddn(0) + h(1, 1) * ddn(1) + h(1, 2) * ddn(2)));
  //         // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(2, 0) * ddn(0) + h(2, 1) * ddn(1) + h(2, 2) * ddn(2)));

  //         // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
  //         // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
  //         // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

  //         ddK_cu(0) = -(ddR(kr, 0) * dn(dirr, us) + ddR(ks, 0) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
  //         ddK_cu(1) = -(ddR(kr, 1) * dn(dirr, us) + ddR(ks, 1) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
  //         ddK_cu(2) = -(ddR(kr, 2) * dn(dirr, us) + ddR(ks, 2) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

  //         ddK_ca_0(ur, us) = (Tb * ddK_cu)(0);
  //         ddK_ca_1(ur, us) = (Tb * ddK_cu)(1);
  //         ddK_ca_2(ur, us) = (Tb * ddK_cu)(2);
  //       }
  //     }
  //   }

  //   // Leaving this here momentarily, further on it is done again, should delete and stick with this
  //   int nLayers = myPatch->getNLayers();

  //   Matrix A(3, 3);   // Membrane stiffness matrix
  //   Matrix B(3, 3);   // Coupling stifness matrix
  //   Matrix D(3, 3);   // Bending stiffness matrix
  //   Matrix T(3, 3);   // Rotating laminate matrix
  //   Matrix Cbar(3, 3); // Equivalent constitutive matrix
  //   for (int capa = 0; capa < nLayers; ++capa)
  //   {
  //     double iAngle     = myPatch->getAngle(capa);
  //     double iThickness = myPatch->getThickness(capa);
  //     double iZ         = myPatch->getZk(capa);// Mid center laminate position, 1 in the meantime

  //     const Matrix& Czig = materialPointers[gp][capa]->getTangent();

  //     T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
  //     T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
  //     T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

  //     Cbar.addMatrixTripleProduct(0, T, Czig, 1.0);

  //     A.addMatrix(1.0, Cbar, iThickness);
  //     B.addMatrix(1.0, Cbar, iThickness * iZ);
  //     D.addMatrix(1.0, Cbar, (iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12));
  //   }

  //   Vector N_ca = A * E_ca + B * K_ca; // Membrane forces
  //   Vector M_ca = B * E_ca + D * K_ca; // Bending moments

  //   // opserr << "A = " << A << endln;
  //   // opserr << "B = " << B << endln;
  //   // opserr << "D = " << D << endln;


  //   Matrix ke(3 * noFuncs, 3 * noFuncs);
  //   Matrix kem(3 * noFuncs, 3 * noFuncs);
  //   Matrix keb(3 * noFuncs, 3 * noFuncs);



  //   for (int ur = 0; ur < 3 * noFuncs ; ++ur)
  //   {
  //     ur_only.fill(ur);


  //     Matrix dN_ca = A * dE_ca(xyz, ur_only) + B * dK_ca(xyz, ur_only);
  //     Matrix dM_ca = B * dE_ca(xyz, ur_only) + D * dK_ca(xyz, ur_only);

  //     for (int us = 0; us <= ur; ++us)
  //     {
  //       us_only.fill(us);

  //       // Vectors for non-linear part
  //       Vector ddE_ca(3);
  //       ddE_ca(0) = ddE_ca_0(ur, us);
  //       ddE_ca(1) = ddE_ca_1(ur, us);
  //       ddE_ca(2) = ddE_ca_2(ur, us);

  //       Vector ddK_ca(3);
  //       ddK_ca(0) = ddK_ca_0(ur, us);
  //       ddK_ca(1) = ddK_ca_1(ur, us);
  //       ddK_ca(2) = ddK_ca_2(ur, us);

  //       // Membrane stiffness
  //       kem(ur, us) = (transpose(3, 1, dN_ca) * dE_ca(xyz, us_only))(0, 0); // Linear part

  //       // Bending stiffness
  //       keb(ur, us) = (transpose(3, 1, dM_ca) * dK_ca(xyz, us_only))(0, 0); // Linear part


  //       if (nonLinearGeometry) // 1 if nonLinear, 0 if Linear
  //       {
  //         kem(ur, us) += N_ca ^ ddE_ca; // Non-linear part
  //         keb(ur, us) += M_ca ^ ddK_ca; // Non-linear part
  //       }
  //       // Symmetric parts
  //       kem(us, ur) = kem(ur, us);
  //       keb(us, ur) = keb(ur, us);
  //     }
  //   }

  //   ke = (kem + keb) * (dA * J2 * wt); //Stiffness matrix on gauss point
  //   *stiff += ke; // Kiendl stiffness

  // }




  // // opserr << "K.noRows() = " << K.noRows() << endln;
  // // opserr << "K.noCols() = " << K.noCols() << endln;
  // // opserr << "ke.noRows() = " << ke.noRows() << endln;
  // // opserr << "ke.noCols() = " << ke.noCols() << endln;

  // // opserr << "Ke = " << memStiff * transpose(Bmem.noRows(), Bmem.noCols(), Bmem) * C * Bmem * J1 * J2 * wt + benStiff * transpose(Bben.noRows(), Bben.noCols(), Bben) * C * Bben * J1 * J2 * wt << endln;


  // // opserr << "Finished making Ke!!" << endln << endln;
  // // opserr << "K = " << K << endln;

  // // opserr << "K = " << K << endln;
  // // return K ;
  // // opserr << "*stiff = " << *stiff << endln;
  // return *stiff;
}

//return secant matrix
const Matrix&  IGAKLShell::getInitialStiff( )
{
  opserr << "IGAKLShell::getInitialStiff - eleTag" << this->getTag() << " called! " << endln;

  Matrix& K = *stiff;

  // K = getTangentStiff();

  // Rellenar K(i,j)  (inicial) aqui (no es necesario por ahora)

  return K ;
}


// //return mass matrix
// const Matrix&  IGAKLShell::getMass( )
// {
//   // opserr << "IGAKLShell::getMass - eleTag" << this->getTag() << " called! " << endln;

//   double t; // thickness [m]
//   double rho; //kg/m^3

//   double W; //Total laminate weight
//   W = 0.0;
//   int nLayers = myPatch->getNLayers();

//   for (int capa = 0; capa < nLayers; ++capa)
//   {
//     rho = (OPS_getNDMaterial(myPatch->getMatTag(capa)))->getRho(); // Density of material
//     t = myPatch->getThickness(capa);
//     W += rho * t;
//   }

//   mass->Zero();

//   int noFuncs = myPatch->getNoFuncs();

//   static Vector R(noFuncs);
//   static Vector dRdxi(noFuncs);
//   static Vector dRdeta(noFuncs);
//   static Vector dR2dxi(noFuncs);
//   static Vector dR2deta(noFuncs);
//   static Vector dR2dxideta(noFuncs);

//   R.resize(noFuncs);
//   dRdxi.resize(noFuncs);
//   dRdeta.resize(noFuncs);
//   dR2dxi.resize(noFuncs);
//   dR2deta.resize(noFuncs);
//   dR2dxideta.resize(noFuncs);

//   static Matrix N(3, 3 * noFuncs); N.resize(3, 3 * noFuncs); N.Zero();

//   float wt;
//   float ptU;
//   float ptV;
//   double xi;
//   double eta;

//   double J2;

//   Vector point(3);
//   static Matrix pts(connectedExternalNodes.Size(), 3);
//   pts.resize(connectedExternalNodes.Size(), 3);

//   // Get pts span matrix for jacobian
//   for (int i = 0; i < connectedExternalNodes.Size(); ++i)
//   {
//     // point = nodePointers[i]->getCrds();
//     point = nodePointers[i]->getTrialDisp();
//     for (int j = 0; j < 3; ++j)
//     {
//       pts(i, j) = point(j);
//     }
//   }

//   // Loop over integrations points
//   for (int gp = 0; gp < ngauss; ++gp)
//   {
//     wt = (*quadWeight)(gp);
//     ptU = (*quadPoint)(gp, 0);
//     ptV = (*quadPoint)(gp, 1);
//     xi = myPatch->parent2ParametricSpace(xiE, ptU);
//     eta = myPatch->parent2ParametricSpace(etaE, ptV);
//     J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

//     R.Zero(); dRdxi.Zero(); dRdeta.Zero(); dR2dxi.Zero(); dR2deta.Zero(); dR2dxideta.Zero();
//     myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

//     // Get first order jacobian
//     static Matrix dr(2, noFuncs); dr.resize(2, noFuncs); dr.Zero();

//     for (int i = 0; i < noFuncs; ++i)
//     {
//       dr(0, i) = dRdxi(i);
//       dr(1, i) = dRdeta(i);
//     }

//     //compute the jacobian of physical and parameter domain mapping
//     // then the derivative w.r.t spatial physical coordinates

//     Matrix jacob = dr * pts; 

//     // a1, a2 and a3 vectors (surface basis vectors)
//     // and its derivatives

//     static Vector a1(3);
//     static Vector a2(3);

//     for (int i = 0; i < 3; ++i)
//     {
//       a1(i) = jacob(0, i);
//       a2(i) = jacob(1, i);
//     }
//     Vector a3 = LovelyCrossProduct(a1, a2);
//     double norma = a3.Norm();
//     a3 /= norma;
//     double J1 = norma;

//     // Forming N matrix
//     N.Zero();
//     for (int i = 0; i < noFuncs; ++i)
//     {
//       N(0, 3 * i) = R(i);
//       N(1, 3 * i + 1) = R(i);
//       N(2, 3 * i + 2) = R(i);
//     }

//     // M += transpose(3, 3 * noFuncs, N) * N * rho * J1 * J2 * wt * t;
//     // M.addMatrixTransposeProduct(1.0, N, N, rho * J1 * J2 * wt * t);
//     mass -> addMatrixTransposeProduct(1.0, N, N, W * J1 * J2 * wt );
//     // (*mass).addMatrixTransposeProduct(1.0, N, N, W * J1 * wt );

//   }

//   // opserr << "Finished making M!!! " << endln << endln;
//   // opserr << "M = " << M << endln;

//   // Rellenar M(i,j)   aqui

//   return *mass ;
// }

//return mass matrix
const Matrix&  IGAKLShell::getMass( )
{
  // opserr << "IGAKLShell::getMass - eleTag" << this->getTag() << " called! " << endln;

  double t; // thickness [m]
  double rho; //kg/m^3

  double W; //Total laminate weight
  W = 0.0;
  int nLayers = myPatch->getNLayers();

  for (int capa = 0; capa < nLayers; ++capa)
  {
    rho = (OPS_getNDMaterial(myPatch->getMatTag(capa)))->getRho(); // Density of material
    t = myPatch->getThickness(capa);
    W += rho * t;
  }

  mass->Zero();

  int noFuncs = myPatch->getNoFuncs();

  static Vector R(noFuncs);
  static Vector dRdxi(noFuncs);
  static Vector dRdeta(noFuncs);
  static Vector dR2dxi(noFuncs);
  static Vector dR2deta(noFuncs);
  static Vector dR2dxideta(noFuncs);

  R.resize(noFuncs);
  dRdxi.resize(noFuncs);
  dRdeta.resize(noFuncs);
  dR2dxi.resize(noFuncs);
  dR2deta.resize(noFuncs);
  dR2dxideta.resize(noFuncs);

  static Matrix N(3, 3 * noFuncs); N.resize(3, 3 * noFuncs);

  float wt;
  float ptU;
  float ptV;
  double xi;
  double eta;

  double J2;

  Vector point(3);
  static Matrix pts(connectedExternalNodes.Size(), 3); pts.resize(connectedExternalNodes.Size(), 3); pts.Zero();

  // Get pts span matrix for jacobian
  for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  {
    point = nodePointers[i]->getCrds();

    // point = nodePointers[i]->getTrialDisp();

    // point = nodePointers[i]->getDisp();
    // point += nodePointers[i]->getIncrDisp();

    for (int j = 0; j < 3; ++j)
    {
      pts(i, j) = point(j);
    }
  }

  // Loop over integrations points
  for (int gp = 0; gp < ngauss; ++gp)
  {
    wt = (*quadWeight)(gp);
    ptU = (*quadPoint)(gp, 0);
    ptV = (*quadPoint)(gp, 1);
    xi = myPatch->parent2ParametricSpace(xiE, ptU);
    eta = myPatch->parent2ParametricSpace(etaE, ptV);
    J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

    R.Zero(); dRdxi.Zero(); dRdeta.Zero(); dR2dxi.Zero(); dR2deta.Zero(); dR2dxideta.Zero();
    myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

    // Get first order jacobian
    Matrix dr(2, noFuncs);

    for (int i = 0; i < noFuncs; ++i)
    {
      dr(0, i) = dRdxi(i);
      dr(1, i) = dRdeta(i);
    }

    //compute the jacobian of physical and parameter domain mapping
    // then the derivative w.r.t spatial physical coordinates

    Matrix jacob = dr * pts;

    // a1, a2 and a3 vectors (surface basis vectors)
    // and its derivatives

    Vector a1(3);
    Vector a2(3);

    for (int i = 0; i < 3; ++i)
    {
      a1(i) = jacob(0, i);
      a2(i) = jacob(1, i);
    }
    Vector a3 = LovelyCrossProduct(a1, a2);
    double norma = a3.Norm();
    a3 /= norma;
    double J1 = norma;

    // Forming N matrix
    N.Zero();
    for (int i = 0; i < noFuncs; ++i)
    {
      N(0, 3 * i) = R(i);
      N(1, 3 * i + 1) = R(i);
      N(2, 3 * i + 2) = R(i);
    }

    // M += transpose(3, 3 * noFuncs, N) * N * rho * J1 * J2 * wt * t;
    // M.addMatrixTransposeProduct(1.0, N, N, rho * J1 * J2 * wt * t);
    mass->addMatrixTransposeProduct(1.0, N, N, W * J1 * J2 * wt );
    // (*mass).addMatrixTransposeProduct(1.0, N, N, W * J1 * wt );

  }

  // opserr << "Finished making M!!! " << endln << endln;
  // opserr << "M = " << M << endln;

  // Rellenar M(i,j)   aqui

  return *mass ;
}

void IGAKLShell::shellGeo(Matrix G, Matrix H, Vector& G3, double& dA, Vector& N, Matrix& Gab, Vector& Bv, Matrix& T_Gcon_E, Matrix& T_E_G, Matrix& T_G_E) // Get geometric quantities

{
  // opserr << "Starting shellGeo!" << endln;


  //G: Covariant base vectors, needed in the reference configuration
  //H: Hessian matrix, needed in the reference configuration


  // Basis vector G3 (covariant) on reference configuration
  G3(0) = G(1, 0) * G(2, 1) - G(2, 0) * G(1, 1);
  G3(1) = G(2, 0) * G(0, 1) - G(0, 0) * G(2, 1);
  G3(2) = G(0, 0) * G(1, 1) - G(1, 0) * G(0, 1);

  // differential area dA = length of G3
  dA = G3.Norm();

  // Normal vector N
  N = G3 / dA;

  // Curvature coefficients in vector notation
  Bv.Zero();
  for (int i = 0; i < 3; ++i)
  {
    Bv(i) = H(0, i) * N(0) + H(1, i) * N(1) + H(2, i) * N(2) ;
  }

  // Covariant metric coefficients Gab
  Gab.addMatrixTransposeProduct(0, G, G, 1);

  // Contravariant metric coefficients Gab_con and contravariant base vectors G_con
  double invdetGab = 1.0 / ( Gab(0, 0) * Gab(1, 1) - Gab(0, 1) * Gab(0, 1) );
  // Testing this
  // double invdetGab = 1.0 / ( Gab(0, 0) * Gab(1, 1) - Gab(0, 1) * Gab(1, 0) );
  Matrix Gab_con(2, 2);
  Gab_con(0, 0) =  invdetGab * Gab(1, 1);
  Gab_con(0, 1) = -invdetGab * Gab(0, 1);
  Gab_con(1, 1) =  invdetGab * Gab(0, 0);
  Gab_con(1, 0) =  Gab_con(0, 1);
  // Gab_con(1, 0) =  -invdetGab*Gab(1, 0); // It should be this

  Matrix G_con(3, 2);
  G_con.Zero();
  G_con = G * transpose(2, 2, Gab_con);


  // Local cartesian coordinates
  Vector g1(3);
  Vector g2_con(3);
  for (int i = 0; i < 3; ++i)
  {
    g1(i) = G(i, 0);
    g2_con(i) = G_con(i, 1);
  }
  double lg1 = g1.Norm();
  double lg2_con = g2_con.Norm();

  Matrix E(3, 2);
  for (int i = 0; i < 3; ++i)
  {
    E(i, 0) = g1(i) / lg1;
    E(i, 1) = g2_con(i) / lg2_con;
  }



  // Transformation matrix T_Gcon_E from G_con to E
  // with 2 in last row for strain in Voigt notation
  Matrix EG(2, 2);
  EG.addMatrixTransposeProduct(0, E, G_con, 1);

  T_Gcon_E(0, 0) = pow(EG(0, 0), 2)        ; T_Gcon_E(0, 1) = pow(EG(0, 1), 2)        ; T_Gcon_E(0, 2) = 2 * EG(0, 0) * EG(0, 1)                       ;
  T_Gcon_E(1, 0) = pow(EG(1, 0), 2)        ; T_Gcon_E(1, 1) = pow(EG(1, 1), 2)        ; T_Gcon_E(1, 2) = 2 * EG(1, 0) * EG(1, 1)                       ;
  T_Gcon_E(2, 0) = 2 * EG(0, 0) * EG(1, 0) ; T_Gcon_E(2, 1) = 2 * EG(0, 1) * EG(1, 1) ; T_Gcon_E(2, 2) = 2 * EG(0, 0) * EG(1, 1) + EG(0, 1) * EG(1, 0) ;


  // Transformation matrix T_E_G from E to G (for PK2 stress)

  T_E_G(0, 0) = pow(EG(0, 0), 2)  ; T_E_G(0, 1) = pow(EG(1, 0), 2) ;  T_E_G(0, 2) = 2 * EG(0, 0) * EG(1, 0)               ;
  T_E_G(1, 0) = pow(EG(0, 1), 2)  ; T_E_G(1, 1) = pow(EG(1, 1), 2) ;  T_E_G(1, 2) = 2 * EG(0, 1) * EG(1, 1)               ;
  T_E_G(2, 0) = EG(0, 0) * EG(0, 1) ; T_E_G(2, 1) = EG(1, 0) * EG(1, 1);  T_E_G(2, 2) = EG(0, 0) * EG(1, 1) + EG(1, 0) * EG(0, 1) ;


  // Transformation matrix T_g_e from g to e (for Cauchy stress)
  EG.addMatrixTransposeProduct(0, E, G, 1);

  T_G_E(0, 0) = pow(EG(0, 0), 2)  ; T_G_E(0, 1) = pow(EG(0, 1), 2) ;  T_G_E(0, 2) = 2 * EG(0, 0) * EG(0, 1)               ;
  T_G_E(1, 0) = pow(EG(1, 0), 2)  ; T_G_E(1, 1) = pow(EG(1, 1), 2) ;  T_G_E(1, 2) = 2 * EG(1, 0) * EG(1, 1)               ;
  T_G_E(2, 0) = EG(0, 0) * EG(1, 0) ; T_G_E(2, 1) = EG(0, 1) * EG(1, 1);  T_G_E(2, 2) = EG(0, 0) * EG(1, 1) + EG(0, 1) * EG(1, 0) ;


// opserr << "Finished shellGeo!" << endln;
  // End shell geo function
}


void  IGAKLShell::zeroLoad( )
{
  // opserr << "IGAKLShell::zeroLoad - eleTag" << this->getTag() << " called! " << endln;
  if (load != 0)
  {
    load->Zero(); // Uncomment this for twisted shell
  }

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return ;
}


int
IGAKLShell::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // opserr << "IGAKLShell::addLoad - tag = " << this->getTag() << endln;
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    // opserr << "IGAKLShell::addLoad - type = LOAD_TAG_SelfWeight - tag = " << this->getTag() << endln;
    // added compatability with selfWeight class implemented for all continuum elements, C.McGann, U.W.
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0);
    appliedB[1] += loadFactor * data(1);
    appliedB[2] += loadFactor * data(2);
    return 0;
  } else if (type == LOAD_TAG_SurfaceLoader) {
    pressure = data(0);
    return 0;
  } else if (type == LOAD_TAG_IGAFollowerLoad) {
    // opserr << "IGAKLShell::addLoad - type = LOAD_TAG_IGAFollowerLoad - tag = " << this->getTag() << endln;

    double xi = data(0);
    double eta = data(1);

    int Nnodes = connectedExternalNodes.Size();
    int NDOF = 3 * Nnodes;

    if (load == 0)
    {
      load = new Vector(NDOF);
    }


    if (pointInElement(xi, eta))
    {
      // opserr << "IGAKLShell::addLoad - type = LOAD_TAG_IGAFollowerLoad -tag =  = " << this->getTag() << " called" << endln;
      Vector followerforce(3);

      followerforce(0) = data(2);
      followerforce(1) = data(3);
      followerforce(2) = data(4);

      // Obtain information from current configuration
      int noFuncs = myPatch->getNoFuncs();
      Vector R(noFuncs);
      Vector dRdxi(noFuncs);
      Vector dRdeta(noFuncs);
      Vector dR2dxi(noFuncs);
      Vector dR2deta(noFuncs);
      Vector dR2dxideta(noFuncs);
      R.Zero();
      dRdxi.Zero();
      dRdeta.Zero();
      dR2dxi.Zero();
      dR2deta.Zero();
      dR2dxideta.Zero();

      // opserr << "xi = " << xi << endln;
      // opserr << "eta = " << eta << endln;
      // // if (eta == 0)
      // {
      //   eta = 1;
      // }
      // myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      // myPatch->Nurbs2DBasis2ndDers(xi, 1.0, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // I need this for twisted shell
      myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // I need this for bending plate circle

      // dRdxi /= 16; //I need this for twisted shell
      // dRdeta /= 2; //I need this for twisted shell

      // opserr << "R = " << R << endln;
      // opserr << "dRdxi = " << dRdxi << endln;
      // opserr << "dRdeta = " << dRdeta << endln;

      // myPatch->Nurbs2DBasis2ndDers(1.0, 1.0, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // Hardcoding this to get the normal on the midplane

      // opserr << "R = " << R << endln;
      // opserr << "dRdxi = " << dRdxi << endln;
      // opserr << "dRdeta = " << dRdeta << endln;

      // myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

      // double wt;
      // double ptU;
      // double ptV;
      // double xi;
      // double eta;
      // double J2;

      // opserr << "R = " << R << endln;
      // opserr << "dRdxi = " << dRdxi << endln;

      // for (int gp = 0; gp < ngauss; ++gp)
      // {
      //   wt = (*quadWeight)(gp);
      //   ptU = (*quadPoint)(gp, 0);
      //   ptV = (*quadPoint)(gp, 1);
      //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
      //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
      //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

      //   opserr << "wt = " << wt << endln;
      //   opserr << "ptU = " << ptU << endln;
      //   opserr << "ptV = " << ptV << endln;

      //   R.Zero();
      //   dRdxi.Zero();
      //   dRdeta.Zero();
      //   dR2dxi.Zero();
      //   dR2deta.Zero();
      //   dR2dxideta.Zero();
      //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      // }

      Vector point(3);
      Vector point_disp(3);
      Matrix pts(connectedExternalNodes.Size(), 3);
      Matrix pts_d(connectedExternalNodes.Size(), 3);

      bool nonLinearGeometry = myPatch->getAnalysisType();

      // Get pts span matrix for jacobian
      for (int i = 0; i < connectedExternalNodes.Size(); ++i)
      {
        point = nodePointers[i]->getCrds();
        point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
        for (int j = 0; j < 3; ++j)
        {
          pts(i, j) = point(j);
          if (isnan(point_disp(j)))
          {
            opserr << "Nan found on addLoad = " << endln;
          }
          if (nonLinearGeometry)
          {
            pts_d(i, j) = pts(i, j) + point_disp(j);
          }
          else
          {
            pts_d(i, j) = pts(i, j);
          }
        }
      }

      // Get first and second order jacobians
      Matrix dr(2, noFuncs);
      Matrix dr2(3, noFuncs);


      for (int i = 0; i < noFuncs; ++i)
      {
        dr(0, i) = dRdxi(i);
        dr(1, i) = dRdeta(i);

        dr2(0, i) = dR2dxi(i);
        dr2(1, i) = dR2deta(i);
        dr2(2, i) = dR2dxideta(i);
      }

      Matrix ddR(3, noFuncs);
      for (int i = 0; i < noFuncs; ++i)
      {
        ddR(0, i) = dR2dxi(i);
        ddR(1, i) = dR2deta(i);
        ddR(2, i) = dR2dxideta(i);
      }


      Matrix g = dr  * pts_d; //G: Covariant base vectors, needed in the current configuration
      Matrix h = ddR * pts_d; //H: Hessian matrix, needed in the current configuration

      // Transposing because it is needed in cols drdxi, drdeta
      g = transpose(2, 3, g);
      h = transpose(3, 3, h); // Shouldn't need this


      // Shell geo deformed parameters and temporal
      Vector g3(3); // = g3
      Vector n(3); // = n
      double lg3_tmp; // = dA
      Matrix gab_tmp(2, 2); // = Gab
      Vector bv_tmp(3); // = Bv
      Matrix T_Gcon_E(3, 3);
      Matrix T_E_G_tmp(3, 3);
      Matrix T_G_E(3, 3);


      // Get T_G_E for base transformation and g3
      shellGeo(g , h , g3     , lg3_tmp , n     , gab_tmp   , bv_tmp   , T_Gcon_E, T_E_G_tmp, T_G_E);

      // Local covariant base vectors
      Vector g1(3);
      Vector g2(3);
      for (int i = 0; i < 3; ++i)
      {
        g1(i) = g(i, 0);
        g2(i) = g(i, 1);
      }

      // opserr << "g3 = " << g3 << endln;
      // To transform from local covariant g to local cartesian e
      Vector e1 = T_G_E * g1;
      Vector e2 = T_G_E * g2;

      Vector e3 = n;

      // for (int i = 0; i < 3; ++i)
      // {
      //   if (e3(i)<0)
      //   {
      //     e3(i)*=-1;
      //   }
      // }

      // Vector e1 = T_Gcon_E * g1; // Trying
      // Vector e2 = T_Gcon_E * g2;

      // Normalizing vectors
      e1 /= e1.Norm();
      e2 /= e2.Norm();
      e3 /= e3.Norm();

      opserr << "e1 = " << e1 << endln;
      opserr << "e2 = " << e2 << endln;
      opserr << "e3 = " << e3 << endln;
      opserr << "data(2) = " << data(2) << endln;
      opserr << "data(3) = " << data(3) << endln;
      opserr << "data(4) = " << data(4) << endln;

      followerforce = loadFactor * (data(2) * e1 + data(3) * e2 + data(4) * e3);
      // followerforce = (data(2) * e1 + data(3) * e2 + data(4) * e3);


      // followerforce(0) = 0;

      // opserr << "e3 = " << e3 << endln;

      // opserr << "followerforce = " << followerforce << endln;

      // N debe ser evaluado en xi, eta
      myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      Matrix N(3, 3 * noFuncs);
      for (int i = 0; i < noFuncs; ++i)
      {
        N(0, 3 * i) = R(i);
        N(1, 3 * i + 1) = R(i);
        N(2, 3 * i + 2) = R(i);
      }

      // opserr << "old load = " << *load << endln;

      // Error found here
      load->addMatrixTransposeVector(1.0, N, followerforce, 1.0);

      opserr << "load = " << *load << endln;

      opserr << "e3 = " << e3 << endln;
      opserr << "followerforce = " << followerforce << endln;
      opserr << "loadFactor = " << loadFactor << endln;


      // resid->addMatrixTransposeVector(1.0, N, followerforce, -1.0);

      // opserr << "*load = " << *load << endln;

      // opserr << "*load = " << *load << endln;
      // The error is that the shape functions are not interpolatory, so the above multiplication does not apply the force in the desired DOF, as the shape functions are greater elsewhere.

      // opserr << "xi = " << xi << endln;
      // opserr << "eta = " << eta << endln;
      // opserr << "xiE = " << xiE << endln;
      // opserr << "etaE = " << etaE << endln;
      // opserr << "R = " << R << endln;
      // opserr << "-1*followerforce = " << -1 * followerforce << endln;
      // opserr << "*load = " << *load << endln;
    }

    return 0;

  } else {
    opserr << "ShellMITC4::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
    return -1;
  }
}



int
IGAKLShell::addInertiaLoadToUnbalance(const Vector &accel)
{
  // opserr << "IGAKLShell::addInertiaLoadToUnbalance - eleTag" << this->getTag() << " called! " << endln;

  // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;
  static Vector res(NDOF);   /// ngdl
  res.resize(NDOF);
  res.Zero();
  // Calculo aceleraciones en nodos y los guardo en res o accel
  static Vector NodalAccelerations(NDOF);
  static Vector NodalDisplacements(NDOF);
  NodalAccelerations.resize(NDOF);
  NodalDisplacements.resize(NDOF);
  NodalAccelerations.Zero();
  NodalDisplacements.Zero();
  static Vector accel_i(3);
  static Vector disp_i(3);
  accel_i.Zero();
  disp_i.Zero();

  // Loop through nodes
  int k = 0;
  for (int i = 0; i < Nnodes; ++i)
  {
    Node *node_i = nodePointers[i];
    accel_i = node_i->getTrialAccel();
    disp_i = node_i->getTrialDisp();
    const Vector &Raccel = nodePointers[i]->getRV(accel);

    for (int j = 0; j < 3; ++j)
    {
      NodalAccelerations(k) = Raccel(j);
      NodalDisplacements(k) = disp_i(j);
      k += 1;
      // opserr << "NodalAccelerations(Nnodes * j + i) = " << NodalAccelerations(Nnodes * j + i) << endln;
    }
  }

  int tang_flag = 1 ; //get the tangent
  bool nonLinearGeometry = myPatch->getAnalysisType();

  // formResidAndTangent( tang_flag ) ;

  // Calculo masa M
  //
  // res = this->getMass() * NodalAccelerations;
  if (load == 0) 
    load = new Vector(NDOF);

  load->addMatrixVector(1.0, this->getMass(), NodalAccelerations, -1.0);

  return 0;
}



//get residual
const Vector&  IGAKLShell::getResistingForce() 
{
  // opserr << "IGAKLShell::getResistingForce - eleTag" << this->getTag() << " called! " << endln;

  int tang_flag = 1 ; //get the tangent
  bool nonLinearGeometry = myPatch->getAnalysisType();

  formResidAndTangent( tang_flag ) ;

  // Adding K*NodalDisplacements to resid

  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;
  static Vector res(NDOF);   /// ngdl
  res.resize(NDOF);
  res.Zero();

  // res = this->getResistingForce();

  // Calculo desplazamientos en nodos
  static Vector NodalDisplacements(NDOF);
  NodalDisplacements.resize(NDOF);
  NodalDisplacements.Zero();
  static Vector disp_i(3);
  disp_i.Zero();


  if (nonLinearGeometry == false)
  {
    // Loop through nodes
    int k = 0;
    for (int i = 0; i < Nnodes; ++i)
    {
      Node *node_i = nodePointers[i];
      disp_i = node_i->getTrialDisp();
      for (int j = 0; j < 3; ++j)
      {
        NodalDisplacements(k) = disp_i(j);
        k += 1;
      }
    }

    // Adding K*nodalDisplacements;
    // opserr << "Adding K*nodalDisplacements" << endln;
    *resid += (*stiff) * NodalDisplacements;
    // *resid += 1*(this->getTangentStiff()) * NodalDisplacements;

  }

  // *resid += this->getMass() * NodalAccelerations;


  // subtract external loads
  if (load != 0)
  {
    *resid -= *load;
    // opserr << "this->getTag() = " << this->getTag() << endln;
    // opserr << "*resid = " << *resid << endln;
  }

  // opserr << "*resid = " << *resid << endln;

  return *resid ;

  // resid->Zero();
  // // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;

  // // opserr << "Element number:  = " << this->getTag() << endln;

  // // opserr << "xiE = " << xiE << endln;
  // // opserr << "etaE = " << etaE << endln;

  // Vector gFact = myPatch->getGravityFactors();
  // // opserr << "gFact = " << gFact << endln;


  // double t = 0.0; // thickness [m]
  // double rho = 0.0; //kg/m^3

  // double W = 0.0; //Total laminate weight (kg/m^2)
  // int nLayers = myPatch->getNLayers();

  // for (int capa = 0; capa < nLayers; ++capa)
  // {
  //   rho = (OPS_getNDMaterial(myPatch->getMatTag(capa)))->getRho(); // Density of material
  //   t = myPatch->getThickness(capa); // Thickness
  //   W += rho * t;
  // }



  // bool nonLinearGeometry = myPatch->getAnalysisType();


  // float wt;
  // float ptU;
  // float ptV;
  // double xi;
  // double eta;

  // int noFuncs = myPatch->getNoFuncs();


  // Vector R(noFuncs);
  // Vector dRdxi(noFuncs);
  // Vector dRdeta(noFuncs);
  // Vector dR2dxi(noFuncs);
  // Vector dR2deta(noFuncs);
  // Vector dR2dxideta(noFuncs);

  // double J2; //Jacobian for the parametric mapping

  // Vector point(3);
  // Vector point_disp(3);
  // Matrix pts(connectedExternalNodes.Size(), 3);
  // Matrix pts_d(connectedExternalNodes.Size(), 3);
  // // Get pts span matrix for jacobian
  // for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  // {
  //   point = nodePointers[i]->getCrds();
  //   point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     pts(i, j) = point(j);
  //     if (isnan(point_disp(j)))
  //     {
  //       opserr << "Nan found on formResidAndTangent = " << endln;
  //     }
  //     if (nonLinearGeometry)
  //     {
  //       pts_d(i, j) = pts(i, j) + point_disp(j);
  //     }
  //     else
  //     {
  //       pts_d(i, j) = pts(i, j);
  //     }
  //   }

  // }

  // // opserr << "pts = " << pts << endln;

  // // Loop over integrations points
  // for (int gp = 0; gp < ngauss; ++gp)
  // {
  //   wt = (*quadWeight)(gp);
  //   ptU = (*quadPoint)(gp, 0);
  //   ptV = (*quadPoint)(gp, 1);
  //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
  //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
  //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));
  //   // opserr << "xi = " << xi << endln;
  //   // opserr << "eta = " << eta << endln;
  //   R.Zero();
  //   dRdxi.Zero();
  //   dRdeta.Zero();
  //   dR2dxi.Zero();
  //   dR2deta.Zero();
  //   dR2dxideta.Zero();
  //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
  //   // opserr << "Got derivatives! " <<  endln;

  //   // Get first and second order jacobians
  //   Matrix dr(2, noFuncs);
  //   Matrix dr2(3, noFuncs);


  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     dr(0, i) = dRdxi(i);
  //     dr(1, i) = dRdeta(i);

  //     dr2(0, i) = dR2dxi(i);
  //     dr2(1, i) = dR2deta(i);
  //     dr2(2, i) = dR2dxideta(i);
  //   }

  //   //compute the jacobian of physical and parameter domain mapping
  //   // then the derivative w.r.t spatial physical coordinates, current configuration


  //   Matrix jacob = dr * pts;  //G: Covariant base vectors, needed in the reference configuration
  //   Matrix jacob2 = dr2 * pts; //H: Hessian matrix, needed in the reference configuration


  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //   // !    DESCRIPTION OF THE VARIABLES !
  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  //   // ! g: matrix with the 3 components of the in plane vectors g1, g2
  //   // ! dR: first derivatives of the shape functions
  //   // ! ddR: second derivatives of the shape functions
  //   // ! dn: first variation of the normal vector g3
  //   // ! lg3: length of the normal vector g3

  //   // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
  //   // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
  //   // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
  //   // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

  //   // ! dE_cu: first variation of the membrane strain in local coord. system
  //   // ! dE_ca: first variation of the membrane strain in global coord. system
  //   // ! dK_cu: first variation of the curvature in local coord. system
  //   // ! dK_ca: first variation of the curvature in global coord. system

  //   // ! ddE_cu: second variation of the membrane strain in local coord. system
  //   // ! ddE_ca: second variation of the membrane strain in global coord. system
  //   // ! ddK_cu: second variation of the curvature in local coord. system
  //   // ! ddK_ca: second variation of the curvature in global coord. system

  //   Vector E_cu(3);
  //   Vector K_cu(3);

  //   Vector dE_cu(3);
  //   Matrix dE_ca(3, 3 * noFuncs);
  //   Vector dK_cu(3);
  //   Matrix dK_ca(3, 3 * noFuncs);

  //   Vector fiem(3 * noFuncs);
  //   Vector fieb(3 * noFuncs);
  //   Vector fie(3 * noFuncs);

  //   Matrix ddR(3, noFuncs);
  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     ddR(0, i) = dR2dxi(i);
  //     ddR(1, i) = dR2deta(i);
  //     ddR(2, i) = dR2dxideta(i);
  //   }


  //   Matrix G = dr  * pts;
  //   Matrix H = ddR * pts;

  //   Matrix g = dr  * pts_d;
  //   Matrix h = ddR * pts_d;


  //   G = transpose(2, 3, G); // Transposing because it is needed in cols drdxi, drdeta
  //   g = transpose(2, 3, g);
  //   H = transpose(3, 3, H);
  //   h = transpose(3, 3, h);


  //   // Shell geo reference parameters and temporal
  //   Vector g3_tmp(3);
  //   Vector n_tmp(3);
  //   double dA; // = lG3
  //   Matrix Gab_r(2, 2); // = Gab
  //   Vector Bv_r(3); // = Bv
  //   Matrix Tb(3, 3); // = T_Gcon_E
  //   Matrix T_E_G_tmp(3, 3); // = T_Gcon_E
  //   Matrix T_G_E_tmp(3, 3); // = T_Gcon_E

  //   // Shell geo deformed parameters and temporal
  //   Vector g3(3); // = g3
  //   Vector n(3); // = n
  //   double lg3; // = dA
  //   Matrix gab(2, 2); // = Gab
  //   Vector bv(3); // = Bv
  //   Matrix Tb_tmp(3, 3);

  //   // Call shell geometry functions
  //   shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);
  //   shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp);

  //   //
  //   // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
  //   //

  //   //strain vector [E11,E22,E12] referred to curvilinear coor sys
  //   E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
  //   E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
  //   E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

  //   //curvature vector [K11,K22,K12] referred to curvilinear coor sys
  //   K_cu(0) = -(bv(0) - Bv_r(0));
  //   K_cu(1) = -(bv(1) - Bv_r(1));
  //   K_cu(2) = -(bv(2) - Bv_r(2));


  //   Vector E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys
  //   Vector K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys

  //   //
  //   // Compute the First variation of strain and curvature w.r.t the d.o.f.
  //   //

  //   // Local variables
  //   int kr = 0;
  //   int dirr = 0;
  //   Matrix dg(3, 2);
  //   Matrix dR = transpose(2, noFuncs, dr);
  //   ddR = transpose(3, noFuncs, ddR);
  //   Matrix dg3(3, 3 * noFuncs);
  //   Vector g3dg3(3 * noFuncs);
  //   Vector g3dg3lg3_3(3 * noFuncs);
  //   Matrix dn(3, 3 * noFuncs);

  //   // For "slicing" matrices
  //   ID ur_only(1);
  //   ID us_only(1);
  //   auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

  //   double lg3_3 = pow(lg3, 3);



  //   for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     // opserr << "kr = " << kr << endln;
  //     // opserr << "dirr = " << dirr << endln;

  //     dg(dirr, 0) = dR(kr, 0);
  //     dg(dirr, 1) = dR(kr, 1);

  //     // Strain
  //     dE_cu(0) = dR(kr, 0) * g(dirr, 0);
  //     dE_cu(1) = dR(kr, 1) * g(dirr, 1);
  //     dE_cu(2) = 0.5 * (dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1));

  //     dE_ca(0, ur) = (Tb * dE_cu)(0);
  //     dE_ca(1, ur) = (Tb * dE_cu)(1);
  //     dE_ca(2, ur) = (Tb * dE_cu)(2);


  //     // Curvature
  //     dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
  //     dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
  //     dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);

  //     g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
  //     g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;


  //     dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
  //     dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
  //     dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);

  //     dK_cu(0) = -(ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     dK_cu(1) = -(ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     dK_cu(2) = -(ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // CHECK THIS OUT, I HAD TO SWAP THE INDICES FOR H

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(0, 1) * dn(1, ur) + h(0, 2) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(1, 0) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(1, 2) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(2, 0) * dn(0, ur) + h(2, 1) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     dK_ca(0, ur) = (Tb * dK_cu)(0);
  //     dK_ca(1, ur) = (Tb * dK_cu)(1);
  //     dK_ca(2, ur) = (Tb * dK_cu)(2);
  //   }

  //   int nLayers = myPatch->getNLayers();

  //   Matrix A(3, 3);   // Membrane stiffness matrix
  //   Matrix B(3, 3);   // Coupling stifness matrix
  //   Matrix D(3, 3);   // Bending stiffness matrix
  //   Matrix T(3, 3);   // Rotating laminate matrix
  //   Matrix Cbar(3, 3); // Equivalent constitutive matrix
  //   for (int capa = 0; capa < nLayers; ++capa)
  //   {
  //     double iAngle     = myPatch->getAngle(capa);
  //     double iThickness = myPatch->getThickness(capa);
  //     double iZ         = myPatch->getZk(capa);// Mid center laminate position, 1 in the meantime

  //     const Matrix& Czig = materialPointers[gp][capa]->getTangent();

  //     T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
  //     T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
  //     T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

  //     Cbar.addMatrixTripleProduct(0, T, Czig, 1.0);

  //     A.addMatrix(1.0, Cbar, iThickness);
  //     B.addMatrix(1.0, Cbar, iThickness * iZ);
  //     D.addMatrix(1.0, Cbar, (iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12));
  //   }

  //   Vector N_ca = A * E_ca + B * K_ca; // Membrane forces
  //   Vector M_ca = B * E_ca + D * K_ca; // Bending moments


  //   for (int ur = 0; ur < 3 * noFuncs ; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     ur_only.fill(ur);

  //     Matrix dN_ca = A * dE_ca(xyz, ur_only) + B * dK_ca(xyz, ur_only);
  //     Matrix dM_ca = B * dE_ca(xyz, ur_only) + D * dK_ca(xyz, ur_only);

  //     Vector dE_ca_here(3);
  //     dE_ca_here(0) = dE_ca(xyz, ur_only)(0, 0);
  //     dE_ca_here(1) = dE_ca(xyz, ur_only)(1, 0);
  //     dE_ca_here(2) = dE_ca(xyz, ur_only)(2, 0);

  //     Vector dK_ca_here(3);
  //     dK_ca_here(0) = dK_ca(xyz, ur_only)(0, 0);
  //     dK_ca_here(1) = dK_ca(xyz, ur_only)(1, 0);
  //     dK_ca_here(2) = dK_ca(xyz, ur_only)(2, 0);

  //     // Force vectors
  //     fiem(ur) = (N_ca ^ dE_ca_here);
  //     fieb(ur) = (M_ca ^ dK_ca_here);

  //   }

  //   Matrix N(3, 3 * noFuncs);
  //   Vector pressure_vector = pressure * g3;

  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     N(0, 3 * i) = R(i);
  //     N(1, 3 * i + 1) = R(i);
  //     N(2, 3 * i + 2) = R(i);
  //   }


  //   fie = (fiem + fieb) * (dA * J2 * wt); //Residual vector on gauss point

  //   fie.addMatrixTransposeVector(1.0, N, pressure_vector, dA * J2 * wt);

  //   // Body forces
  //   Vector b(appliedB, 3);

  //   if (applyLoad == 1)
  //   {
  //     fie.addMatrixTransposeVector(1.0, N, b, W * dA * J2 * wt);
  //   } else {
  //     fie.addMatrixTransposeVector(1.0, N, gFact, W * dA * J2 * wt);
  //   }

  //   *resid += fie; // Kiendl stiffness


  // }

  // // Aqui se rellena el vector de fuerzas residuales...

  // // // subtract external loads
  // // if (load != 0)
  // //   *resid -= *load;

  // // Adding follower load vector to resid (substracting actually to make it positive)
  // // opserr << "*resid = " << *resid << endln;
  // // opserr << "*load = " << -1*(*load) << endln;
  // *resid-=(*load);
  // // load->Zero();

  // return *resid ;


}


void IGAKLShell::formStrainsLinear()
{
  // Gauss integration parameters
  float wt;
  float ptU;
  float ptV;
  double xi;
  double eta;

  // Number of shape functions and layers
  int noFuncs = myPatch->getNoFuncs();
  int nLayers = myPatch->getNLayers();

  static Vector E_cu(3);
  static Vector E_ca(3);
  static Vector K_cu(3);
  static Vector K_ca(3);

  static Matrix pts(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  static Matrix pts_d(noFuncs, 3); // Matrix for storing element nodes displaced coordinates

  static Matrix G_t(2, 3);
  static Matrix H_t(3, 3);
  static Matrix g_t(2, 3);
  static Matrix h_t(3, 3);
  // Shell geo reference parameters and temporal
  static Vector g3_tmp(3);
  static Vector n_tmp(3);
  static double dA = 0; // = lG3
  static Matrix Gab_r(2, 2);// = Gab
  static Vector Bv_r(3); // = Bv
  static Matrix Tb(3, 3); // = T_Gcon_E
  static Matrix T_E_G_tmp(3, 3); // = T_E_G
  static Matrix T_G_E_tmp(3, 3); // = T_G_E
  static Matrix T_G_E(3, 3); // = T_G_E for testing
  // Shell geo deformed parameters and temporal
  static Vector g3(3); // = g3
  static Vector n(3); // = n
  static double lg3 = 0; // = dA
  static Matrix gab(2, 2); // = Gab
  static Vector bv(3); // = Bv
  static Matrix Tb_tmp(3, 3);

  static Matrix dR(noFuncs, 3);
  static Matrix ddR(noFuncs, 3);

  static Matrix dr(3, noFuncs);
  static Matrix ddr(3, noFuncs);

  // Getting ply matrices
  static Matrix A(3, 3);   // Membrane stiffness matrix
  static Matrix B(3, 3);   // Coupling stifness matrix
  static Matrix D(3, 3);   // Bending stiffness matrix
  static Matrix T(3, 3);   // Rotating laminate matrix
  static Matrix C(3, 3); // Equivalent constitutive matrix
  static double W; //Total laminate weight (kg/m^2)
  static Vector Ezeta(3);
  static double iAngle;    
  static double iThickness;
  static double iZ;        
  static double iRho;      

  static Vector R(noFuncs);  // 
  static Vector dRdxi(noFuncs);
  static Vector dRdeta(noFuncs);
  static Vector dR2dxi(noFuncs);
  static Vector dR2deta(noFuncs);
  static Vector dR2dxideta(noFuncs);

  //Jacobian for the parametric mapping
  static double J2;


  // Resizing
  R.resize(noFuncs);  // 
  dRdxi.resize(noFuncs);
  dRdeta.resize(noFuncs);
  dR2dxi.resize(noFuncs);
  dR2deta.resize(noFuncs);
  dR2dxideta.resize(noFuncs);
  pts.resize(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  pts_d.resize(noFuncs, 3); // Matrix for storing element nodes displaced coordinates


  // START OF CALCULATION
  pts.Zero();
  pts_d.Zero();

  // Get pts span matrix for jacobian
  for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  {
    const Vector& point  = nodePointers[i]->getCrds() ; // Vector for storing node coordinates
    const Vector& point_disp = nodePointers[i]->getTrialDisp(); // Vector for storing node displacement
    const Vector& point_disp_commited = nodePointers[i]->getDisp(); // Vector for storing node displacement

    // opserr << "nodePointers[i]->getDisp() = " << (nodePointers[i]->getDisp()) << endln;
    // opserr << "nodePointers[i]->getTrialDisp() = " << (nodePointers[i]->getTrialDisp()) << endln;

    for (int j = 0; j < 3; ++j)
    {
      pts(i, j) = point(j);
      if (isnan(point_disp(j)))
      {
        opserr << "Nan found on formResidAndTangent = " << endln;
      }
      // I need to add the last commited disp to get the current configuration while using linear algorithm
      // pts_d(i, j) = pts(i, j);
      pts_d(i, j) = pts(i, j) + point_disp_commited(j);;
    }
  }

  // Loop over integrations points
  for (int gp = 0; gp < ngauss; ++gp)
  {
    wt = (*quadWeight)(gp);
    ptU = (*quadPoint)(gp, 0);
    ptV = (*quadPoint)(gp, 1);

    xi = myPatch->parent2ParametricSpace(xiE, ptU);
    eta = myPatch->parent2ParametricSpace(etaE, ptV);
    J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

    R.Zero();
    dRdxi.Zero();
    dRdeta.Zero();
    dR2dxi.Zero();
    dR2deta.Zero();
    dR2dxideta.Zero();
    myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    // !    DESCRIPTION OF THE VARIABLES !
    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    // ! g: matrix with the 3 components of the in plane vectors g1, g2
    // ! dR: first derivatives of the shape functions
    // ! ddR: second derivatives of the shape functions
    // ! dn: first variation of the normal vector g3
    // ! lg3: length of the normal vector g3

    // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
    // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
    // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
    // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys


    E_cu.Zero();
    E_ca.Zero();
    K_cu.Zero();
    K_ca.Zero();

    // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
    for (int i = 0; i < noFuncs; ++i)
    {
      dr(0, i) = dRdxi(i);
      dr(1, i) = dRdeta(i);
      dr(2, i) = 0;
      ddr(0, i) = dR2dxi(i);
      ddr(1, i) = dR2deta(i);
      ddr(2, i) = dR2dxideta(i);
    }

    G_t = dr  * pts;
    H_t = ddr * pts;

    g_t = dr  * pts_d;
    h_t = ddr * pts_d;

    // Transposing because it is needed in cols drdxi, drdeta, etc
    // G is finally in shape (3,2) and H(3,3)
    Matrix G(3,2); G.Zero();
    Matrix g(3,2); g.Zero();
    Matrix H(3,3); H.Zero();
    Matrix h(3,3); h.Zero();

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        G(i,j) = G_t(j,i);
        g(i,j) = g_t(j,i);
      }
    }

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        H(i,j) = H_t(j,i);
        h(i,j) = h_t(j,i);
      }
    }
    
    // G = transpose(2, 3, G);
    // g = transpose(2, 3, g);

    // H = transpose(3, 3, H);
    // h = transpose(3, 3, h);



    for (int i = 0; i < noFuncs; ++i)
    {
      dR(i, 0) = dRdxi(i);
      dR(i, 1) = dRdeta(i);
      dR(i, 2) = 0;

      ddR(i, 0) = dR2dxi(i);
      ddR(i, 1) = dR2deta(i);
      ddR(i, 2) = dR2dxideta(i);
    }


    Gab_r.Zero(); Bv_r.Zero(); Tb.Zero(); T_G_E.Zero(); g3.Zero(); n.Zero(); gab.Zero(); bv.Zero();
    // Call shell geometry functions
    shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E);  // Reference configuration
    shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp); // Deformed or actual configuration


    //
    // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
    //

    E_cu.Zero();
    K_cu.Zero();

    //strain vector [E11,E22,E12] referred to curvilinear coor sys
    E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
    E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
    E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

    //curvature vector [K11,K22,K12] referred to curvilinear coor sys
    K_cu(0) = -(bv(0) - Bv_r(0));
    K_cu(1) = -(bv(1) - Bv_r(1));
    K_cu(2) = -(bv(2) - Bv_r(2));


    E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys, Tb is T_Gcon_E
    K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys, , Tb is T_Gcon_E

    A.Zero();         // Membrane stiffness matrix
    B.Zero();         // Coupling stifness matrix
    D.Zero();         // Bending stiffness matrix
    W = 0.0; //Total laminate weight (kg/m^2)

    // opserr << "Getting Tangent" << endln;

    for (int capa = 0; capa < nLayers; ++capa)
    {
      iAngle     = myPatch -> getAngle(capa);
      iThickness = myPatch -> getThickness(capa);
      iZ         = myPatch -> getZk(capa);// Mid center laminate position
      iRho       = (OPS_getNDMaterial(myPatch -> getMatTag(capa))) -> getRho(); // Density of material kg/m^3;

      W += iRho * iThickness;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      // Aqui le pedimos al material su C

      // calcular Ezeta
      Ezeta.Zero();
      Ezeta = T * (E_ca + iZ * K_ca); // Cartesian

      // Giving the trail strain to the material, rotated and in cartesian coordinates
      materialPointers[gp][capa] -> setTrialStrain(Ezeta); 
    }

  }




}

void IGAKLShell::formResidAndTangent( int tang_flag )
{

  // opserr << "IGAKLShell::formResidAndTangent - eleTag" << this->getTag() << " called! " << endln;

  //zero stiffness and residual
  stiff->Zero();
  resid->Zero();

  bool nonLinearGeometry = myPatch->getAnalysisType();

  Vector gFact = myPatch->getGravityFactors();

  // Gauss integration parameters
  float wt;
  float ptU;
  float ptV;
  double xi;
  double eta;

  // Number of shape functions and layers
  int noFuncs = myPatch->getNoFuncs();
  int nLayers = myPatch->getNLayers();


  // Start of static declarations
  static Matrix ke(3 * noFuncs, 3 * noFuncs);
  static Matrix kem(3 * noFuncs, 3 * noFuncs);
  static Matrix keb(3 * noFuncs, 3 * noFuncs);

  static Vector E_cu(3);
  static Vector E_ca(3);
  static Vector K_cu(3);
  static Vector K_ca(3);
  static Vector dE_cu(3);
  static Matrix dE_ca(3, 3 * noFuncs);
  static Vector dK_cu(3);
  static Matrix dK_ca(3, 3 * noFuncs);
  static Vector N_ca(3);
  static Vector M_ca(3);
  static Vector fiem(3 * noFuncs);
  static Vector fieb(3 * noFuncs);
  static Vector fie(3 * noFuncs);

  static Matrix pts(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  static Matrix pts_d(noFuncs, 3); // Matrix for storing element nodes displaced coordinates

  static Matrix G_t(2, 3);
  static Matrix H_t(3, 3);
  static Matrix g_t(2, 3);
  static Matrix h_t(3, 3);
  // Shell geo reference parameters and temporal
  static Vector g3_tmp(3);
  static Vector n_tmp(3);
  static double dA = 0; // = lG3
  static Matrix Gab_r(2, 2);// = Gab
  static Vector Bv_r(3); // = Bv
  static Matrix Tb(3, 3); // = T_Gcon_E
  static Matrix T_E_G_tmp(3, 3); // = T_E_G
  static Matrix T_G_E_tmp(3, 3); // = T_G_E
  static Matrix T_G_E(3, 3); // = T_G_E for testing
  // Shell geo deformed parameters and temporal
  static Vector g3(3); // = g3
  static Vector n(3); // = n
  static double lg3 = 0; // = dA
  static Matrix gab(2, 2); // = Gab
  static Vector bv(3); // = Bv
  static Matrix Tb_tmp(3, 3);

  static Matrix dR(noFuncs, 3);
  static Matrix ddR(noFuncs, 3);

  static Matrix dr(3, noFuncs);
  static Matrix ddr(3, noFuncs);

  // Getting ply matrices
  static Matrix A(3, 3);   // Membrane stiffness matrix
  static Matrix B(3, 3);   // Coupling stifness matrix
  static Matrix D(3, 3);   // Bending stiffness matrix
  static Matrix T(3, 3);   // Rotating laminate matrix
  static Matrix C(3, 3); // Equivalent constitutive matrix
  static double W; //Total laminate weight (kg/m^2)
  static Vector Ezeta(3);
  static double iAngle;    
  static double iThickness;
  static double iZ;        
  static double iRho;      

  // Variables for stiffness assembly
  static double lg3_3;
  static double lg3_5;
  static int kr;
  static int dirr;
  static Matrix dg(3, 2);
  static Matrix dg3(3, 3 * noFuncs);
  static Vector g3dg3(3 * noFuncs);
  static Vector g3dg3lg3_3(3 * noFuncs);
  static Matrix dn(3, 3 * noFuncs);

  static Vector ddE_cu(3);
  static Matrix ddE_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  static Matrix ddE_ca_1(3 * noFuncs, 3 * noFuncs);
  static Matrix ddE_ca_2(3 * noFuncs, 3 * noFuncs);
 
  static Vector ddK_cu(3);
  static Matrix ddK_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  static Matrix ddK_ca_1(3 * noFuncs, 3 * noFuncs);
  static Matrix ddK_ca_2(3 * noFuncs, 3 * noFuncs);

  static int ks;
  static int dirs;
  static int dirt;
  static int ddir;
  static double tmp1;
  static double tmp2;
  static Vector ddn(3);
  static Vector ddg3(3);

  static Vector Szeta(3);

  static Vector dE_ca_ur(3);
  static Vector dK_ca_ur(3);
  static Vector dN_ca(3);
  static Vector dM_ca(3);
  static Vector ddE_ca(3);
  static Vector ddK_ca(3);
  static Vector dE_ca_us(3);
  static Vector dK_ca_us(3);

  static Vector R(noFuncs);  // 
  static Vector dRdxi(noFuncs);
  static Vector dRdeta(noFuncs);
  static Vector dR2dxi(noFuncs);
  static Vector dR2deta(noFuncs);
  static Vector dR2dxideta(noFuncs);



  // Start of resizings
  ke.resize(3 * noFuncs, 3 * noFuncs);
  kem.resize(3 * noFuncs, 3 * noFuncs);
  keb.resize(3 * noFuncs, 3 * noFuncs);
  dE_ca.resize(3, 3 * noFuncs);
  dK_ca.resize(3, 3 * noFuncs);
  fiem.resize(3 * noFuncs);
  fieb.resize(3 * noFuncs);
  fie.resize(3 * noFuncs);
  pts.resize(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  pts_d.resize(noFuncs, 3); // Matrix for storing element nodes displaced coordinates
  dg3.resize(3, 3 * noFuncs);
  g3dg3.resize(3 * noFuncs);
  g3dg3lg3_3.resize(3 * noFuncs);
  dn.resize(3, 3 * noFuncs);
  ddE_ca_0.resize(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  ddE_ca_1.resize(3 * noFuncs, 3 * noFuncs);
  ddE_ca_2.resize(3 * noFuncs, 3 * noFuncs);
  ddK_ca_0.resize(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  ddK_ca_1.resize(3 * noFuncs, 3 * noFuncs);
  ddK_ca_2.resize(3 * noFuncs, 3 * noFuncs);
  R.resize(noFuncs);  // 
  dRdxi.resize(noFuncs);
  dRdeta.resize(noFuncs);
  dR2dxi.resize(noFuncs);
  dR2deta.resize(noFuncs);
  dR2dxideta.resize(noFuncs);

  //Jacobian for the parametric mapping
  static double J2;


  // START OF CALCULATION
  pts.Zero();
  pts_d.Zero();


  // Get pts span matrix for jacobian
  for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  {
    const Vector& point  = nodePointers[i]->getCrds() ; // Vector for storing node coordinates
    const Vector& point_disp = nodePointers[i]->getTrialDisp(); // Vector for storing node displacement
    const Vector& point_disp_commited = nodePointers[i]->getDisp(); // Vector for storing node displacement

    // opserr << "nodePointers[i]->getDisp() = " << (nodePointers[i]->getDisp()) << endln;
    // opserr << "nodePointers[i]->getTrialDisp() = " << (nodePointers[i]->getTrialDisp()) << endln;

    for (int j = 0; j < 3; ++j)
    {
      pts(i, j) = point(j);
      if (isnan(point_disp(j)))
      {
        opserr << "Nan found on formResidAndTangent = " << endln;
      }
      if (nonLinearGeometry)
      {
        pts_d(i, j) = pts(i, j) + point_disp(j);
      }
      else
      {
        // Pipe : should I add getDisp() here? The strains are always zero otherwise, and maybe 
        // I need to add the last commited disp to get the current configuration while using linear algorithm
        // pts_d(i, j) = pts(i, j);
        pts_d(i, j) = pts(i, j) + point_disp_commited(j);;
      }
    }
  }


  // Loop over integrations points
  for (int gp = 0; gp < ngauss; ++gp)
  {
    wt = (*quadWeight)(gp);
    ptU = (*quadPoint)(gp, 0);
    ptV = (*quadPoint)(gp, 1);

    // opserr << " ptU = " << ptU << endln;
    // opserr << " ptV = " << ptV << endln;
    // opserr << " wt = " << wt << endln;

    xi = myPatch->parent2ParametricSpace(xiE, ptU);
    eta = myPatch->parent2ParametricSpace(etaE, ptV);
    J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

    R.Zero();
    dRdxi.Zero();
    dRdeta.Zero();
    dR2dxi.Zero();
    dR2deta.Zero();
    dR2dxideta.Zero();
    myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    // !    DESCRIPTION OF THE VARIABLES !
    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    // ! g: matrix with the 3 components of the in plane vectors g1, g2
    // ! dR: first derivatives of the shape functions
    // ! ddR: second derivatives of the shape functions
    // ! dn: first variation of the normal vector g3
    // ! lg3: length of the normal vector g3

    // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
    // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
    // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
    // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

    // ! dE_cu: first variation of the membrane strain in local coord. system
    // ! dE_ca: first variation of the membrane strain in global coord. system
    // ! dK_cu: first variation of the curvature in local coord. system
    // ! dK_ca: first variation of the curvature in global coord. system

    // ! ddE_cu: second variation of the membrane strain in local coord. system
    // ! ddE_ca: second variation of the membrane strain in global coord. system
    // ! ddK_cu: second variation of the curvature in local coord. system
    // ! ddK_ca: second variation of the curvature in global coord. system

    ke.Zero(); 

    kem.Zero(); E_cu.Zero(); dE_cu.Zero(); E_ca.Zero(); dE_ca.Zero();
    N_ca.Zero();

    keb.Zero(); K_cu.Zero(); dK_cu.Zero(); K_ca.Zero();dK_ca.Zero();
    M_ca.Zero();


    // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
    for (int i = 0; i < noFuncs; ++i)
    {
      dr(0, i) = dRdxi(i);
      dr(1, i) = dRdeta(i);
      dr(2, i) = 0;
      ddr(0, i) = dR2dxi(i);
      ddr(1, i) = dR2deta(i);
      ddr(2, i) = dR2dxideta(i);
    }

    G_t = dr  * pts;
    H_t = ddr * pts;

    g_t = dr  * pts_d;
    h_t = ddr * pts_d;

    // Transposing because it is needed in cols drdxi, drdeta, etc
    // G is finally in shape (3,2) and H(3,3)
    Matrix G(3,2); G.Zero();
    Matrix g(3,2); g.Zero();
    Matrix H(3,3); H.Zero();
    Matrix h(3,3); h.Zero();

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        G(i,j) = G_t(j,i);
        g(i,j) = g_t(j,i);
      }
    }

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        H(i,j) = H_t(j,i);
        h(i,j) = h_t(j,i);
      }
    }

    // opserr << 'G = ' << G << endln;
    // opserr << 'H = ' << H << endln;
    
    // G = transpose(2, 3, G);
    // g = transpose(2, 3, g);

    // H = transpose(3, 3, H);
    // h = transpose(3, 3, h);



    for (int i = 0; i < noFuncs; ++i)
    {
      dR(i, 0) = dRdxi(i);
      dR(i, 1) = dRdeta(i);
      dR(i, 2) = 0;

      ddR(i, 0) = dR2dxi(i);
      ddR(i, 1) = dR2deta(i);
      ddR(i, 2) = dR2dxideta(i);
    }


    Gab_r.Zero(); Bv_r.Zero(); Tb.Zero(); T_G_E.Zero(); g3.Zero(); n.Zero(); gab.Zero(); bv.Zero();
    // Call shell geometry functions
    shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E);  // Reference configuration
    shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp); // Deformed or actual configuration

    // Compute the in-plane rotation of the plies
    double theta1 = 0.0;
    Vector bG1(3);
    bG1(0) = G(0, 0);
    bG1(1) = G(1, 0);
    bG1(2) = G(2, 0);

    if (bG1(0) > 0 && bG1(1) > 0) {
      theta1 = atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) > 0 && bG1(1) < 0) {
      theta1 = 2 * M_PI - atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0)<0 && bG1(1)>0) {
      theta1 = M_PI - atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) < 0 && bG1(1) < 0) {
      theta1 = M_PI + atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) > 0 && bG1(1) == 0) {
      theta1 = 0;

    } else if (bG1(0) < 0 && bG1(1) == 0) {
      theta1 = M_PI;

    } else if (bG1(0) == 0 && bG1(1) > 0) {
      theta1 = 90.0 / 180.0 * M_PI;

    } else if (bG1(0) == 0 && bG1(1) < 0) {
      theta1 = 270.0 / 180.0 * M_PI;
    }

    theta1 = 0; // Not working, making zero in-plane rotation
    //

    //
    // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
    //

    E_cu.Zero();
    K_cu.Zero();

    //strain vector [E11,E22,E12] referred to curvilinear coor sys
    E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
    E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
    E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

    //curvature vector [K11,K22,K12] referred to curvilinear coor sys
    K_cu(0) = -(bv(0) - Bv_r(0));
    K_cu(1) = -(bv(1) - Bv_r(1));
    K_cu(2) = -(bv(2) - Bv_r(2));


    E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys, Tb is T_Gcon_E
    K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys, , Tb is T_Gcon_E

    A.Zero();         // Membrane stiffness matrix
    B.Zero();         // Coupling stifness matrix
    D.Zero();         // Bending stiffness matrix
    W = 0.0; //Total laminate weight (kg/m^2)

    // opserr << "Getting Tangent" << endln;

    for (int capa = 0; capa < nLayers; ++capa)
    {
      iAngle     = myPatch -> getAngle(capa) - theta1;
      iThickness = myPatch -> getThickness(capa);
      iZ         = myPatch -> getZk(capa);// Mid center laminate position
      iRho       = (OPS_getNDMaterial(myPatch -> getMatTag(capa))) -> getRho(); // Density of material kg/m^3;

      W += iRho * iThickness;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      // Aqui le pedimos al material su C

      // calcular Ezeta
      Ezeta.Zero();
      Ezeta = T * (E_ca + iZ * K_ca); // Cartesian

      // Giving the trail strain to the material, rotated and in cartesian coordinates
      materialPointers[gp][capa] -> setTrialStrain(Ezeta); 

      const Matrix& Czig = materialPointers[gp][capa]->getTangent();

      C.addMatrixTripleProduct(0, T, Czig, 1.0);

      A.addMatrix(1.0, C, iThickness);
      B.addMatrix(1.0, C, iThickness * iZ);
      D.addMatrix(1.0, C, iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12.0 );
    }

    lg3_3 = pow(lg3, 3);
    lg3_5 = pow(lg3, 5);


    //
    // Compute the First variation of strain and curvature w.r.t the d.o.f.
    //

    // Local variables
    kr = 0; dirr = 0; dg.Zero(); dg3.Zero(); g3dg3.Zero(); g3dg3lg3_3.Zero(); dn.Zero();

    // For "slicing" matrices
    ID ur_only(1);
    ID us_only(1);
    auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

    for (int ur = 0; ur < 3 * noFuncs; ++ur)
    {
      // Local node number kr and dof direction dirr
      kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
      dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

      dg(dirr, 0) = dR(kr, 0);
      dg(dirr, 1) = dR(kr, 1);


      // ADDING THIS, WEIRD
      // if (dirr == 2)
      // {
      //   dg.Zero();
      // }
      // ADDING THIS, WEIRD


      // First variation of the membrane strain in local coord. system: 1st term on RHS eqn. 1.174 of Computational FSI
      dE_cu(0) = dR(kr, 0) * g(dirr, 0);
      dE_cu(1) = dR(kr, 1) * g(dirr, 1);
      dE_cu(2) = 0.5 * ( dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1) );

      // First variation of the membrane strain in global coord. system: RHS eqn. 1.174 of Computational FSI
      dE_ca(0, ur) = (Tb * dE_cu)(0);
      dE_ca(1, ur) = (Tb * dE_cu)(1);
      dE_ca(2, ur) = (Tb * dE_cu)(2);

      //
      // dE_ca(0, ur) = Tb(0, 0) * dE_cu(0) + Tb(0, 1) * dE_cu(1) + Tb(0, 2) * dE_cu(2);
      // dE_ca(1, ur) = Tb(1, 0) * dE_cu(0) + Tb(1, 1) * dE_cu(1) + Tb(1, 2) * dE_cu(2);
      // dE_ca(2, ur) = Tb(2, 0) * dE_cu(0) + Tb(2, 1) * dE_cu(1) + Tb(2, 2) * dE_cu(2);


      // First variation of the curvature in local coord. system: 1st term on RHS eqn. 1.175 of Computational FSI
      dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
      dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
      dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);


      // opserr << "Artem's" << endln;
      // opserr << 'dg3(0, ur) = ' << dg3(0, ur)*1e-10 << endln;
      // opserr << 'dg3(1, ur) = ' << dg3(1, ur)*1e-10 << endln;
      // opserr << 'dg3(2, ur) = ' << dg3(2, ur)*1e-10 << endln;


      // Kratos code THIS GETS CLOSE!
      if (dirr == 0)
      {
        dg3(0, ur) = 0;
        dg3(1, ur) = -dR(kr,0) * g(2,1) + dR(kr,1) * g(2,0);
        dg3(2, ur) = dR(kr,0) * g(1,1) - dR(kr,1) * g(1,0);
      }
      else if (dirr == 1)
      {
        dg3(0, ur) = dR(kr,0) * g(2,1) - dR(kr,1) * g(2,0); 
        dg3(1, ur) = 0;
        dg3(2, ur) = -dR(kr,0) * g(0,1) + dR(kr,1) * g(0,0);
      }
      else if (dirr == 2)
      {
        dg3(0, ur) = -dR(kr,0) * g(1,1) + dR(kr,1) * g(1,0); 
        dg3(1, ur) = dR(kr,0) * g(0,1) - dR(kr,1) * g(0,0);
        dg3(2, ur) = 0;
      }

      // Kratos code
      // if (dirr == 0)
      // {
      //   dg3(0, ur) = 0;
      //   dg3(1, ur) = dR(kr,0) * g(2,2) - dR(kr,1) * g(1,2); 
      //   dg3(2, ur) = -dR(kr,0) * g(2,1) + dR(kr,1) * g(1,1); 
      // }
      // else if (dirr == 1)
      // {
      //   dg3(0, ur) = -dR(kr,0) * g(2,2) + dR(kr,1) * g(1,2);
      //   dg3(1, ur) = 0;
      //   dg3(2, ur) = dR(kr,0) * g(2,0) - dR(kr,1) * g(1,0);
      // }
      // else if (dirr == 2)
      // {
      //   dg3(0, ur) = dR(kr,0) * g(2,1) - dR(kr,1) * g(1,1);
      //   dg3(1, ur) = -dR(kr,0) * g(2,0) + dR(kr,1) * g(1,0);
      //   dg3(2, ur) = 0;
      // }

      // opserr << "Kratos" << endln;
      // opserr << 'dg3(0, ur) = ' << dg3(0, ur)*1e-10 << endln;
      // opserr << 'dg3(1, ur) = ' << dg3(1, ur)*1e-10 << endln;
      // opserr << 'dg3(2, ur) = ' << dg3(2, ur)*1e-10 << endln;


      g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
      g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;

      dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
      dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
      dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);

      // if (dirr==0 or dirr==1)
      // {
      //   dn(0,ur)=0;
      //   dn(1,ur)=0;
      //   dn(2,ur)=0;
      // }


      // opserr << "dn(0,ur) = " << dn(0,ur) << endln;
      // opserr << "dn(1,ur) = " << dn(1,ur) << endln;
      // opserr << "dn(2,ur) = " << dn(2,ur) << endln << endln;
      // opserr << "h = " << h << endln << endln;

      // opserr << "ddR(kr, 0) * n(dirr) = " << (ddR(kr, 0) * n(dirr)) << endln;
      // opserr << "1*(h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur))) = " << 1*(h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)) << endln;
      // opserr << endln;
      // opserr << "ddR(kr, 1) * n(dirr) = " << (ddR(kr, 1) * n(dirr)) << endln;
      // opserr << "1*(h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)) = " << 1*(h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)) << endln;
      // opserr << endln;
      // opserr << "ddR(kr, 2) * n(dirr) = " << (ddR(kr, 2) * n(dirr)) << endln;
      // opserr << "1*(h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)) = " << 1*(h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)) << endln;
      // opserr << endln;

      // Original , extremely weird behaviour with the 0.00001 part
      // dK_cu(0) = - (ddR(kr, 0) * n(dirr) + 1 * (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur))) ;
      // dK_cu(1) = - (ddR(kr, 1) * n(dirr) + 1 * (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur))) ;
      // dK_cu(2) = - (ddR(kr, 2) * n(dirr) + 1 * (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur))) ;


      // testing
      dK_cu(0) = 0 - (ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur))) ;
      dK_cu(1) = 0 - (ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur))) ;
      dK_cu(2) = 0 - (ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur))) ;


      // Kratos code
      // if (dirr == 0)
      // {
      //   dK_cu(0) = 0 - (ddR(kr,0)*n(0) + h(0,0)*dn(0,ur) + h(1,0)*dn(0,ur)+h(2,0)*dn(0,ur));
      //   dK_cu(1) = 0 - (ddR(kr,0)*n(1) + h(0,0)*dn(1,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(1,ur));
      //   dK_cu(2) = 0 - (ddR(kr,0)*n(2) + h(0,0)*dn(2,ur) + h(1,0)*dn(2,ur)+h(2,0)*dn(2,ur));
      // }
      // else if (dirr == 1)
      // {
      //   dK_cu(0) = 0 - (ddR(kr,2)*n(0) + h(0,1)*dn(0,ur) + h(1,1)*dn(0,ur)+h(2,1)*dn(0,ur));
      //   dK_cu(1) = 0 - (ddR(kr,2)*n(1) + h(0,1)*dn(1,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(1,ur));
      //   dK_cu(2) = 0 - (ddR(kr,2)*n(2) + h(0,1)*dn(2,ur) + h(1,1)*dn(2,ur)+h(2,1)*dn(2,ur));
      // }
      // else if (dirr == 2)
      // {
      //   dK_cu(0) = 0 - (ddR(kr,1)*n(0) + h(0,2)*dn(0,ur) + h(1,2)*dn(0,ur)+h(2,2)*dn(0,ur));
      //   dK_cu(1) = 0 - (ddR(kr,1)*n(1) + h(0,2)*dn(1,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(1,ur));
      //   dK_cu(2) = 0 - (ddR(kr,1)*n(2) + h(0,2)*dn(2,ur) + h(1,2)*dn(2,ur)+h(2,2)*dn(2,ur));
      // }

      // // Kratos code
      if (dirr == 0)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(0) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(0) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(0) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }
      else if (dirr == 1)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(1) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(1) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(1) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }
      else if (dirr == 2)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(2) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(2) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(2) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }



      // First variation of the curvature in global coord. system: RHS eqn. 1.175 of Computational FSI
      dK_ca(0, ur) = (Tb * dK_cu)(0);
      dK_ca(1, ur) = (Tb * dK_cu)(1);
      dK_ca(2, ur) = (Tb * dK_cu)(2);

      // dK_ca(0, ur) = Tb(0, 0) * dK_cu(0) + Tb(0, 1) * dK_cu(1) + Tb(0, 2) * dK_cu(2);
      // dK_ca(1, ur) = Tb(1, 0) * dK_cu(0) + Tb(1, 1) * dK_cu(1) + Tb(1, 2) * dK_cu(2);
      // dK_ca(2, ur) = Tb(2, 0) * dK_cu(0) + Tb(2, 1) * dK_cu(1) + Tb(2, 2) * dK_cu(2);
    }


    // Second variation of strain and curvature w.r.t. dofs
    ddE_cu.Zero(); ddE_ca_0.Zero(); ddE_ca_1.Zero(); ddE_ca_2.Zero(); ddK_cu.Zero(); ddK_ca_0.Zero(); ddK_ca_1.Zero(); ddK_ca_2.Zero();

    // Local variables
    ks = 0; dirs = 0; dirt = 0; ddir = 0; tmp1 = 0.0; tmp2 = 0.0; ddg3.Zero(); ddn.Zero();

    if (nonLinearGeometry) // true if non-linear geometry
    {
      // opserr << "Using nonLinearGeometry!" << endln;
      for (int ur = 0; ur < 3 * noFuncs; ++ur)
      {
        // Local node number kr and dof direction dirr
        kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
        dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

        for (int us = 0; us <= ur; ++us)
        {
          ks = (us + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
          // ks = (us + 2) / 3 ; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
          dirs = us - 3 * (ks); // Takes value 0 1 2 0 1 2 ...

          // Strain
          ddE_cu.Zero();
          if (dirr == dirs)
          {
            ddE_cu(0) = dR(kr, 0) * dR(ks, 0);
            ddE_cu(1) = dR(kr, 1) * dR(ks, 1);
            ddE_cu(2) = 0.5 * ( dR(kr, 0) * dR(ks, 1) + dR(kr, 1) * dR(ks, 0) );
          }

          ddE_ca_0(ur, us) = (Tb * ddE_cu)(0);
          ddE_ca_1(ur, us) = (Tb * ddE_cu)(1);
          ddE_ca_2(ur, us) = (Tb * ddE_cu)(2);

          // Curvature
          ddg3.Zero();

          dirt = 6 - dirr - dirs; // Check this
          ddir = dirr - dirs;

          if (ddir == -1 || ddir == 2)
          {
            // opserr << "dirt - 3 = " << dirt-3 << endln;
            ddg3(dirt - 3) = dR(kr, 0) * dR(ks, 1) - dR(ks, 0) * dR(kr, 1);
          }
          else if (ddir == 1 || ddir == -2)
          {
            // opserr << "dirt - 3 = " << dirt-3 << endln;
            ddg3(dirt - 3) = -dR(kr, 0) * dR(ks, 1) + dR(ks, 0) * dR(kr, 1);
          }

          // ur_only.fill(ur);
          // us_only.fill(us);
          ur_only(0) = ur;
          us_only(0) = us;

          tmp1 = -( (ddg3 ^ g3) + (transpose(3, 1, dg3(xyz, ur_only)) * dg3(xyz, us_only)) (0, 0) ) / lg3_3 ; // the (0,0) is to transform the 1x1 matrix into double
          tmp2 = 3 * g3dg3(ur) * g3dg3(us) / lg3_5;



          ddn(0) = ddg3(0) / lg3 - g3dg3lg3_3(us) * dg3(0, ur) - g3dg3lg3_3(ur) * dg3(0, us) + tmp1 * g3(0) + tmp2 * g3(0);
          ddn(1) = ddg3(1) / lg3 - g3dg3lg3_3(us) * dg3(1, ur) - g3dg3lg3_3(ur) * dg3(1, us) + tmp1 * g3(1) + tmp2 * g3(1);
          ddn(2) = ddg3(2) / lg3 - g3dg3lg3_3(us) * dg3(2, ur) - g3dg3lg3_3(ur) * dg3(2, us) + tmp1 * g3(2) + tmp2 * g3(2);

          // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(0, 1) * ddn(1) + h(0, 2) * ddn(2)));
          // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(1, 0) * ddn(0) + h(1, 1) * ddn(1) + h(1, 2) * ddn(2)));
          // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(2, 0) * ddn(0) + h(2, 1) * ddn(1) + h(2, 2) * ddn(2)));

          // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
          // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
          // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

          ddK_cu(0) = -( ddR(kr, 0) * dn(dirr, us) + ddR(ks, 0) * dn(dirs, ur) + h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2) );
          ddK_cu(1) = -( ddR(kr, 1) * dn(dirr, us) + ddR(ks, 1) * dn(dirs, ur) + h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2) );
          ddK_cu(2) = -( ddR(kr, 2) * dn(dirr, us) + ddR(ks, 2) * dn(dirs, ur) + h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2) );

          ddK_ca_0(ur, us) = (Tb * ddK_cu)(0);
          ddK_ca_1(ur, us) = (Tb * ddK_cu)(1);
          ddK_ca_2(ur, us) = (Tb * ddK_cu)(2);
        }
      }
    }


    // // Computing layer forces, membrane and bending  (asume material lineal)
    // N_ca = A * E_ca + B * K_ca; // Membrane forces
    // M_ca = B * E_ca + D * K_ca; // Bending moments

    // opserr << "N_ca old = " << N_ca << endln;
    // opserr << "M_ca old = " << M_ca << endln;

    N_ca.Zero();
    M_ca.Zero();

    //alternativamente // necesitas para plasticidad
    for (int capa = 0; capa < nLayers; ++capa)
    {
      iThickness = myPatch->getThickness(capa);
      iZ         = myPatch->getZk(capa);// Mid center laminate position
      iAngle     = myPatch->getAngle(capa) - theta1;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      const Matrix& Czig = materialPointers[gp][capa]->getTangent();
      C.addMatrixTripleProduct(0, T, Czig, 1.0);

      Szeta.Zero();
      Szeta = transpose(3, 3, T) * materialPointers[gp][capa]->getStress();  // Rotating stresses

      // Integrating stress in every ply
      N_ca += iThickness * Szeta;
      M_ca += iThickness * iZ * Szeta + (C * K_ca) * pow(iThickness, 3) / 12.0 ;

    }

    // opserr << "N_ca new = " << N_ca << endln;
    // opserr << "M_ca new = " << M_ca << endln;

    kem.Zero();
    keb.Zero();

    // Looping over the dofs ur and us
    for (int ur = 0; ur < 3 * noFuncs ; ++ur)
    {

      dE_ca_ur(0) = dE_ca(0, ur);
      dE_ca_ur(1) = dE_ca(1, ur);
      dE_ca_ur(2) = dE_ca(2, ur);

      dK_ca_ur(0) = dK_ca(0, ur);
      dK_ca_ur(1) = dK_ca(1, ur);
      dK_ca_ur(2) = dK_ca(2, ur);

      dN_ca = A * dE_ca_ur + B * dK_ca_ur;
      dM_ca = B * dE_ca_ur + D * dK_ca_ur;

      for (int us = 0; us <= ur; ++us)
      {
        // Vectors for non-linear part

        ddE_ca(0) = ddE_ca_0(ur, us);
        ddE_ca(1) = ddE_ca_1(ur, us);
        ddE_ca(2) = ddE_ca_2(ur, us);

        ddK_ca(0) = ddK_ca_0(ur, us);
        ddK_ca(1) = ddK_ca_1(ur, us);
        ddK_ca(2) = ddK_ca_2(ur, us);

        // Membrane stiffness
        // kem(ur, us) = (transpose(3, 1, dN_ca) * dE_ca(xyz, us_only)) (0, 0); // Linear part
        dE_ca_us(0) = dE_ca(0, us);
        dE_ca_us(1) = dE_ca(1, us);
        dE_ca_us(2) = dE_ca(2, us);
        // opserr << "dN_ca = " << dN_ca << endln;
        // opserr << "dE_ca_here = " << dE_ca_here << endln;
        // opserr << "dN_ca ^ dE_ca_here = " << (dN_ca ^ dE_ca_here) << endln;
        kem(ur, us) = dN_ca ^ dE_ca_us; // Linear part
        // kem(ur, us) = dN_ca(0) * dE_ca_us(0) + dN_ca(1) * dE_ca_us(1) + dN_ca(2) * dE_ca_us(2); // Linear part
        // opserr << "kem(ur,us) = " << kem(ur,us) << endln;


        // Bending stiffness
        // keb(ur, us) = (transpose(3, 1, dM_ca) * dK_ca(xyz, us_only)) (0, 0); // Linear part
        dK_ca_us(0) = dK_ca(0, us);
        dK_ca_us(1) = dK_ca(1, us);
        dK_ca_us(2) = dK_ca(2, us);
        // opserr << "dM_ca = " << dM_ca << endln;
        // opserr << "dK_ca_here = " << dK_ca_here << endln;
        // opserr << "dM_ca ^ dK_ca_here = " << (dM_ca ^ dK_ca_here) << endln;
        keb(ur, us) = dM_ca ^ dK_ca_us; // Linear part
        // keb(ur, us) = dM_ca(0) * dK_ca_us(0) + dM_ca(1) * dK_ca_us(1) + dM_ca(2) * dK_ca_us(2); // Linear part
        // opserr << "keb(ur,us) = " << keb(ur,us) << endln;

        if (nonLinearGeometry) // 1 if nonLinear, 0 if Linear
        {
          kem(ur, us) += N_ca ^ ddE_ca; // Non-linear part
          keb(ur, us) += M_ca ^ ddK_ca; // Non-linear part
        }
      }

      // Force vectors (Residual)
      if (nonLinearGeometry)
      {
        fiem(ur) = (N_ca ^ dE_ca_ur);
        fieb(ur) = (M_ca ^ dK_ca_ur); 
      }

    }

    // Debugging
    // opserr << "kem = " << (kem * dA * J2 * wt) << endln;
    // opserr << "keb = " << (keb * dA * J2 * wt) << endln;
    // // raise(SIGSEGV);

    // Adding stiffness matrices
    ke.addMatrix(0, kem + keb, dA) ; //Stiffness matrix on gauss point
    for (int ur = 0; ur < 3 * noFuncs; ++ur)
    {
      for (int us = 0; us < 3 * noFuncs; ++us)
      {
        if (us >= ur) {
          ke(ur, us) = ke(us, ur);
        }
      }
    }

    Matrix N(3, 3 * noFuncs);
    Vector pressure_vector = pressure * g3;

    for (int i = 0; i < noFuncs; ++i)
    {
      N(0, 3 * i) = R(i);
      N(1, 3 * i + 1) = R(i);
      N(2, 3 * i + 2) = R(i);
    }

    fie = (fiem + fieb) * (dA * J2 * wt); //Residual vector on gauss point
    fie.addMatrixTransposeVector(1.0, N, pressure_vector, dA * J2 * wt);

    // opserr << "fie = " << fie << endln;

    // Body forces
    Vector b(appliedB, 3);

    if (applyLoad == 1)
    {
      fie.addMatrixTransposeVector(1.0, N, b, -1 * W * dA * J2 * wt);
    } else {
      fie.addMatrixTransposeVector(1.0, N, gFact, -1 * W * dA * J2 * wt);
    }

    // Add contribution to residual vector
    *resid += fie;

    // Add contribution to stiffness matrix if tang_flag == 1
    if (tang_flag == 1 ) {
      // *stiff += ke;
      stiff->addMatrix(1, ke, J2 * wt);
    }

    // opserr << "ke = " << ke << endln;

  }
  // opserr << "*stiff = " << *stiff << endln;
  // raise(SIGSEGV);
}

void IGAKLShell::formTangentNguyen()
{
  // opserr << "IGAKLShell::getTangentStiff - eleTag" << this->getTag() << " called! " << endln;

  // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;

  // opserr << "Element number:  = " << this->getTag() << endln;

  // opserr << "xiE = " << xiE << endln;
  // opserr << "etaE = " << etaE << endln;

  // Matrix& K = *stiff;
  // K.Zero();
  stiff->Zero();

  // opserr << "xiE = " << xiE << endln;
  // opserr << "etaE = " << etaE << endln;

  float wt;
  float ptU;
  float ptV;
  double xi;
  double eta;

  int noFuncs = myPatch->getNoFuncs();

  Vector R(noFuncs);
  Vector dRdxi(noFuncs);
  Vector dRdeta(noFuncs);
  Vector dR2dxi(noFuncs);
  Vector dR2deta(noFuncs);
  Vector dR2dxideta(noFuncs);

  double J2; //Jacobian for the parametric mapping

  int ders; //Only because Nurbs2dBasis2ndDers returns int

  Vector point(3);
  Matrix pts(connectedExternalNodes.Size(), 3);
  // Get pts span matrix for jacobian
  for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  {
    point = nodePointers[i]->getCrds();
    for (int j = 0; j < 3; ++j)
    {
      pts(i, j) = point(j);
    }
  }

  // opserr << "pts = " << pts << endln;

  // Loop over integrations points
  for (int gp = 0; gp < ngauss; ++gp)
  {
    wt = (*quadWeight)(gp);
    ptU = (*quadPoint)(gp, 0);
    ptV = (*quadPoint)(gp, 1);
    xi = myPatch->parent2ParametricSpace(xiE, ptU);
    eta = myPatch->parent2ParametricSpace(etaE, ptV);
    J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));
    // opserr << "xi = " << xi << endln;
    // opserr << "eta = " << eta << endln;
    R.Zero();
    dRdxi.Zero();
    dRdeta.Zero();
    dR2dxi.Zero();
    dR2deta.Zero();
    dR2dxideta.Zero();
    ders = myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
    // opserr << "Got derivatives! " <<  endln;

    // Get first and second order jacobians
    Matrix dr(2, noFuncs);
    Matrix dr2(3, noFuncs);

    for (int i = 0; i < noFuncs; ++i)
    {
      dr(0, i) = dRdxi(i);
      dr(1, i) = dRdeta(i);

      dr2(0, i) = dR2dxi(i);
      dr2(1, i) = dR2deta(i);
      dr2(2, i) = dR2dxideta(i);
    }

    //compute the jacobian of physical and parameter domain mapping
    // then the derivative w.r.t spatial physical coordinates


    Matrix jacob = dr * pts;
    Matrix jacob2 = dr2 * pts;

    double dxdxi = jacob(0, 0);
    double dydxi = jacob(0, 1);
    double dxdet = jacob(1, 0);
    double dydet = jacob(1, 1);

    Matrix j33(3, 3);

    j33(0, 0) = dxdxi * dxdxi; j33(1, 0) = dydxi * dydxi; j33(2, 0) = 2 * dxdxi * dydxi;
    j33(0, 1) = dxdet * dxdet; j33(1, 1) = dydet * dydet; j33(2, 1) = 2 * dxdet * dydet;
    j33(0, 2) = dxdxi * dxdet; j33(1, 2) = dydxi * dydet; j33(2, 2) = dxdxi * dydet + dxdet * dydxi;

    // a1, a2 and a3 vectors (surface basis vectors)
    // and its derivatives

    Vector a1(3);
    Vector a2(3);

    for (int i = 0; i < 3; ++i)
    {
      a1(i) = jacob(0, i);
      a2(i) = jacob(1, i);
    }
    Vector a3 = LovelyCrossProduct(a1, a2);
    double norma = a3.Norm();
    a3 /= norma;
    double J1 = norma;

    Vector a11(3);
    Vector a22(3);
    Vector a12(3);

    for (int i = 0; i < 3; ++i)
    {
      a11(i) = jacob2(0, i);
      a22(i) = jacob2(1, i);
      a12(i) = jacob2(2, i);
    }

    //dot products of ai and ei

    double a1e1 = a1(0);
    double a1e2 = a1(1);
    double a1e3 = a1(2);
    double a2e1 = a2(0);
    double a2e2 = a2(1);
    double a2e3 = a2(2);

    //R_I,2*a1 + R_I,1*a2 for all shape functions
    Matrix dRIa(3, noFuncs);
    for (int i = 0; i < noFuncs; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        dRIa(j, i) = dRdeta(i) * a1(j) + dRdxi(i) * a2(j);
      }
    }

    //compute the constitutive matrix C
    double a_11 = a1 ^ a1;
    double a_12 = a1 ^ a2;
    double a_21 = a2 ^ a1;
    double a_22 = a2 ^ a2;

    // making ihat,jhat and khat
    Vector i1(3);
    Vector i2(3);
    i1(0) = 1;
    i2(1) = 1;
    Vector i3 = LovelyCrossProduct(i1, i2); //Testing cross product

    Vector aa1(2);
    Vector aa2(2);

    Matrix aa1_mat(2, 2);
    Matrix aa2_mat(2, 2);

    aa1_mat(0, 0) = a_11; aa1_mat(0, 1) = a_21;
    aa1_mat(1, 0) = a_12; aa1_mat(1, 1) = a_22;

    aa2_mat(0, 0) = a_11; aa2_mat(0, 1) = a_21;
    aa2_mat(1, 0) = a_12; aa2_mat(1, 1) = a_22;


    // A*x=b
    aa1_mat.Solve(i1, aa1);
    aa2_mat.Solve(i2, aa2);

    double au11 = aa1(0);
    double au12 = aa1(1);
    double au22 = aa2(1);

    // HARDCODING THIS
    double nu;// = 0.3; // -
    double E;// = 2.1e11; //Young's modulus N/m^2

    // double memStiff = E * t / (1 - nu * nu);
    // double benStiff = E * t * t * t / 12 / (1 - nu * nu);

    Matrix C(3, 3);

    // C (0, 0) = au11 * au11;                               C (0, 1) = nu * au11 * au22 + (1 - nu) * au12 * au12; C (0, 2) = au11 * au12;
    // C (1, 0) = nu * au11 * au22 + (1 - nu) * au12 * au12; C (1, 1) = au22 * au22;                               C (1, 2) = au22 * au12;
    // C (2, 0) = au11 * au12;                               C (2, 1) = au22 * au12;                               C (2, 2) = 0.5 * ((1 - nu) * au11 * au22 + (1 + nu) * au12 * au12);

    // C --> punteros a plane-stress materials   getTangent  (por punto gauss y capa) //Armamos A, B y D

    // opserr << "\nNow integrating through thickness" << endln;
    int nLayers = myPatch->getNLayers();

    Matrix A(3, 3);   // Membrane stiffness matrix
    Matrix B(3, 3);   // Coupling stifness matrix
    Matrix D(3, 3);   // Bending stiffness matrix
    Matrix T(3, 3);   // Rotating laminate matrix
    Matrix Cbar(3, 3); // Equivalent constitutive matrix
    double stiffMat; // Factor for constituvie composite matrix

    // opserr << "nLayers = " << nLayers << endln;
    for (int capa = 0; capa < nLayers; ++capa)
    {
      double iAngle     = myPatch->getAngle(capa);
      double iThickness = myPatch->getThickness(capa);
      double iZ         = myPatch->getZk(capa);// Mid center laminate position, 1 in the meantime

      // opserr << "iAngle = " << iAngle << endln;
      // opserr << "iThickness = " << iThickness << endln;
      // opserr << "iZ = " << iZ << endln;


      const Matrix& Czig = materialPointers[gp][capa]->getTangent();

      nu = Czig(0, 1) / Czig(0, 0);
      E = Czig(0, 0) * (1 - nu * nu);
      stiffMat = E / (1 - nu * nu);

      C(0, 0) = au11 * au11;                               C(0, 1) = nu * au11 * au22 + (1 - nu) * au12 * au12; C(0, 2) = au11 * au12;
      C(1, 0) = nu * au11 * au22 + (1 - nu) * au12 * au12; C(1, 1) = au22 * au22;                               C(1, 2) = au22 * au12;
      C(2, 0) = au11 * au12;                               C(2, 1) = au22 * au12;                               C(2, 2) = 0.5 * ((1 - nu) * au11 * au22 + (1 + nu) * au12 * au12);

      C *= stiffMat;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      Cbar.addMatrixTripleProduct(0, T, C, 1.0);

      A.addMatrix(1.0, Cbar, iThickness);
      B.addMatrix(1.0, Cbar, iThickness * iZ);
      D.addMatrix(1.0, Cbar, iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12 );

      // opserr << "A = " << A << endln;
      // opserr << "B = " << B << endln;
      // opserr << "D = " << D << endln;
    }
    // opserr << "Finished making composite laminate " << endln;
    // opserr << "A = " << A << endln;
    // opserr << "B = " << B << endln;
    // opserr << "D = " << D << endln;

    // opserr << "A = " << A << endln;
    // opserr << "B = " << B << endln;
    // opserr << "D = " << D << endln;


    // for (int capa = 0; capa < nLayers; capa++)
    // {
    //   opserr << "capa = " << matTags(capa) << endln;
    // }

    // Vector theta = myPatch->theta;




    // Membrane and bending B matrices
    Matrix Bmem(3, noFuncs * 3);
    Bmem.Zero();
    Matrix Bben(3, noFuncs * 3);
    Bben.Zero();

    float dRIdx;
    float dRIdy;
    Vector BI1(3);
    Vector BI2(3);
    Vector BI3(3);
    ID id_(3);

    for (int i = 0; i < noFuncs; ++i)
    {
      dRIdx = dRdxi(i);
      dRIdy = dRdeta(i);

      for (int j = 0; j < 3; ++j)
      {
        id_(j) = 3 * i + j;
      }

      Bmem(0, id_(0)) = dRIdx * a1e1; Bmem(0, id_(1)) = dRIdx * a1e2; Bmem(0, id_(2)) = dRIdx * a1e3;
      Bmem(1, id_(0)) = dRIdy * a2e1; Bmem(1, id_(1)) = dRIdy * a2e2; Bmem(1, id_(2)) = dRIdy * a2e3;
      Bmem(2, id_(0)) = dRIa(0, i);   Bmem(2, id_(1)) = dRIa(1, i);   Bmem(2, id_(2)) = dRIa(2, i);


      BI1 = -dR2dxi(i) * a3 + 1 / norma * (dRIdx * LovelyCrossProduct(a11, a2) + dRIdy * LovelyCrossProduct(a1, a11)\
                                           +(a3 ^ a11) * (dRIdx * LovelyCrossProduct(a2, a3) + dRIdy * LovelyCrossProduct(a3, a1)));

      BI2 = -dR2deta(i) * a3 + 1 / norma * (dRIdx * LovelyCrossProduct(a22, a2) + dRIdy * LovelyCrossProduct(a1, a22)\
                                            +(a3 ^ a22) * (dRIdx * LovelyCrossProduct(a2, a3) + dRIdy * LovelyCrossProduct(a3, a1)));

      BI3 = -dR2dxideta(i) * a3 + 1 / norma * (dRIdx * LovelyCrossProduct(a12, a2) + dRIdy * LovelyCrossProduct(a1, a12)\
            +(a3 ^ a12) * (dRIdx * LovelyCrossProduct(a2, a3) + dRIdy * LovelyCrossProduct(a3, a1)));

      Bben(0, id_(0)) = BI1(0);     Bben(0, id_(1)) = BI1(1);     Bben(0, id_(2)) = BI1(2);
      Bben(1, id_(0)) = BI2(0);     Bben(1, id_(1)) = BI2(1);     Bben(1, id_(2)) = BI2(2);
      Bben(2, id_(0)) = 2 * BI3(0); Bben(2, id_(1)) = 2 * BI3(1); Bben(2, id_(2)) = 2 * BI3(2);

    }

    // Filling stiffness matrix
    // K.addMatrixTripleProduct(1.0, Bmem, C, memStiff * J1 * J2 * wt);
    // K.addMatrixTripleProduct(1.0, Bben, C, benStiff * J1 * J2 * wt);

    stiff->addMatrixTripleProduct(1.0, Bmem, A, J1 * J2 * wt);
    stiff->addMatrixTripleProduct(1.0, Bmem, B, Bben, J1 * J2 * wt);
    stiff->addMatrixTripleProduct(1.0, Bben, B, Bmem, J1 * J2 * wt);
    stiff->addMatrixTripleProduct(1.0, Bben, D, J1 * J2 * wt);

    // opserr << "Ke = " <<  transpose(Bben.noRows(), Bben.noCols(), Bben) * D * Bben * J1 * J2 * wt << endln;
  }

  // opserr << "Finished making Ke!!" << endln << endln;
  // opserr << "K = " << K << endln;
  // return *stiff;
  // *stiff+=K;
  // return K ;
}

//get residual with inertia terms
const Vector&  IGAKLShell::getResistingForceIncInertia( )
{
  // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;
  static Vector res(NDOF);   /// ngdl
  res.resize(NDOF);
  res.Zero();
  // Calculo aceleraciones en nodos y los guardo en res o accel
  static Vector NodalAccelerations(NDOF);
  static Vector NodalDisplacements(NDOF);
  NodalAccelerations.resize(NDOF);
  NodalDisplacements.resize(NDOF);
  NodalAccelerations.Zero();
  NodalDisplacements.Zero();
  static Vector accel_i(3);
  static Vector disp_i(3);
  accel_i.Zero();
  disp_i.Zero();


  // Loop through nodes
  int k = 0;
  for (int i = 0; i < Nnodes; ++i)
  {
    Node *node_i = nodePointers[i];
    accel_i = node_i->getTrialAccel();
    disp_i = node_i->getTrialDisp();

    for (int j = 0; j < 3; ++j)
    {
      NodalAccelerations(k) = accel_i(j);
      NodalDisplacements(k) = disp_i(j);
      k += 1;
      // opserr << "NodalAccelerations(Nnodes * j + i) = " << NodalAccelerations(Nnodes * j + i) << endln;
    }
  }

  int tang_flag = 1 ; //get the tangent
  bool nonLinearGeometry = myPatch->getAnalysisType();

  // formResidAndTangent( tang_flag ) ;

  // Calculo masa M
  //
  res = this->getResistingForce() + this->getMass() * NodalAccelerations;
  *resid=res;
  // res = this->getMass() * NodalAccelerations;

  // res = *resid;
  // res += this->getMass() * NodalAccelerations;
  // *resid += this->getMass() * NodalAccelerations;

  // if (nonLinearGeometry == false)
  // {
  //   // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  //   res += 1 * (this->getTangentStiff() * NodalDisplacements);
  //   // res += (*stiff) * NodalDisplacements;

  //   // *resid += 1 * (this->getTangentStiff() * NodalDisplacements);
  //   // *resid += (*stiff) * NodalDisplacements;
  // }

  // add the damping forces if rayleigh damping
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0){
    
    // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called Rayleigh! " << endln;
    res += this->getRayleighDampingForces();
    // *resid += this->getRayleighDampingForces();
  }

  // // subtract external loads
  // if (load != 0)
  // {
  //   // *resid -= *load;
  //   res -= *load;
  //   opserr << "this->getTag() = " << this->getTag() << endln;
  //   // opserr << "*resid = " << *resid << endln;
  // }

  // res = *resid;
  // *resid = res;

  // opserr << "M = " << this->getMass() << endln;

  // opserr << "NodalAccelerations = " << NodalAccelerations << endln;

  // opserr << "this->getMass() * NodalAccelerations = " << this->getMass() * NodalAccelerations << endln;

  // *resid += this->getMass() * NodalAccelerations;
  // res = getResistingForce() + this->getMass() * NodalAccelerations;

  // res = this->resid

  // opserr << "Finished IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  // opserr << "res = " << res;
  *resid=res;
  // return res;
  return *resid;
}


//**********************************************************************

Matrix
IGAKLShell::transpose( int dim1,
                       int dim2,
                       const Matrix &M )
{
  int i ;
  int j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {
    for ( j = 0; j < dim2; j++ )
      Mtran(j, i) = M(i, j) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  IGAKLShell::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  return res;
}

int  IGAKLShell::recvSelf (int commitTag,
                           Channel &theChannel,
                           FEM_ObjectBroker &theBroker)
{
  int res = 0;



  return res;
}
//**************************************************************************

int
IGAKLShell::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

  int error = 0;

  return error;
}


bool IGAKLShell::pointInElement(double xi, double eta) const
{
  bool in_xi = xiE(0) <= xi && xi <= xiE(1);
  bool in_eta = etaE(0) <= eta && eta <= etaE(1);
  return in_xi && in_eta;
}

#define PARAM_ADVANCESTATE 1
#define PARAM_RESETSTRESSES 2

int IGAKLShell :: setParameter(const char **argv, int argc, Parameter &param)
{
  // opserr << "IGAKLShell :: setParameter element with tag " << this->getTag() << "called" << endln;
  if (argc < 1)
    return -1;

  int res = -1;
  int matRes = res;

  for (int gp = 0; gp < ngauss; gp++ ) {
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa) {
      // opserr << "Material calling setParameter" << endln;
      matRes = materialPointers[gp][capa]->setParameter(argv, argc, param);
      if (matRes != -1) {
        res = matRes;
      } else {
        opserr << "IGAKLShell :: setParameter - failed" << endln;
        res = -1;
      }
    }
  }
  return res;
}

int IGAKLShell :: updateParameter(int parameterID, Information &info)
{
  int res = -1;
  int matRes = res;

  if (parameterID == res)
  {
    return -1;
  } else {
    for (int gp = 0; gp < ngauss; gp++ ) {
      for (int capa = 0; capa < myPatch->getNLayers(); ++capa) {
        matRes = materialPointers[gp][capa]->updateParameter(parameterID, info);
        if (matRes != -1) {
          res = matRes;
        } else {
          opserr << "IGAKLShell :: updateParameter failed " << endln;
          res = -1;
        }
      }
    }
  }
  return res;
}
