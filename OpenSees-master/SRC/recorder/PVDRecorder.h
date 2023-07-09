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
                                                                        
// $Revision: 1.0 $
// $Date: 2015-11-12 $
                                                                        
#ifndef PVDRecorder_h
#define PVDRecorder_h

// Written: Minjie Zhu
//
// Description: This file contains the class definition for 
// PVDRecorder. A PVDRecorder is used to store all responses in pvd format.


#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <ID.h>
#include <Recorder.h>

class Node;
class Element;

class PVDRecorder: public Recorder
{
public:
    struct NodeData {
	NodeData():disp(false),vel(false),accel(false),incrdisp(false),reaction(false),
		   rayleigh(false),pressure(false),unbalanced(false),
		   mass(false),numeigen(0){}
	bool disp,vel,accel,incrdisp,reaction,rayleigh,pressure,unbalanced,mass;
	int numeigen;
    };
    typedef std::vector<std::string> EleData;
    
public:
    PVDRecorder(const char *filename, const NodeData& ndata,
		const std::vector<EleData>& edata, int ind=2, int pre=10, double dt=0, double relDeltaTTol = 0.00001);
    PVDRecorder();
    ~PVDRecorder();

    int record(int commitTag, double timeStamp);
    int restart();
    int flush();
    int domainChanged();    
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

    virtual int vtu();
    virtual int pvd();
    virtual void addEleData(const EleData& edata) {eledata.push_back(edata);}

private:
    virtual void indent();
    virtual void incrLevel() {indentlevel++;}
    virtual void decrLevel() {indentlevel--;}
    virtual void getParts();
    virtual int savePart(int partno, int ctag, int ndf);
    virtual int savePart0(int ndf);
    virtual int savePartParticle(int partno, int gtag, int ndf);
    void getfilename(const char* name);
    
private:
    int indentsize, precision, indentlevel;
    std::string pathname,basename;
    std::vector<double> timestep;
    std::vector<ID> timeparts;
    std::ofstream theFile;
    char quota;
    std::map<int,ID> parts;
    NodeData nodedata;
    std::vector<EleData> eledata;
    Domain* theDomain;
    std::map<int,int> partnum;
    double dT, nextTime;
    double relDeltaTTol;

public:
    enum VtkType {
	VTK_VERTEX=1,VTK_POLY_VERTEX=2,VTK_LINE=3,VTK_POLY_LINE=4,
	VTK_TRIANGLE=5,VTK_TRIANGLE_STRIP=6,VTK_POLYGON=7,
	VTK_PIXEL=8,VTK_QUAD=9,VTK_TETRA=10,VTK_VOXEL=11,
	VTK_HEXAHEDRON=12,VTK_WEDGE=12,VTK_PYRAMID=14,
	VTK_QUADRATIC_EDGE=21,VTK_QUADRATIC_TRIANGLE=22,
	VTK_QUADRATIC_QUAD=23,VTK_QUADRATIC_TETRA=24,
	VTK_QUADRATIC_HEXAHEDRON=25
    };
    static void setVTKType();

private:
    static std::map<int,VtkType> vtktypes;

};

#endif
