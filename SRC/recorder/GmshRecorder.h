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
                                                                        
                                                                        
#ifndef GmshRecorder_h
#define GmshRecorder_h

// Written: Jose Abell
//
// Description: This file contains the class definition for 
// GmshRecorder. A GmshRecorder is used to store mesh and all responses in gmsh format.


#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <ID.h>
#include <Recorder.h>


#define GMSHRECORDER_MAX_FILENAME_SIZE 100

// #define GMSHRECORDER_JAABELL_EXPERIMENTAL 1

class Node;
class Element;

class GmshRecorder: public Recorder
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
    GmshRecorder(const char *filename, const NodeData& ndata,
		const std::vector<EleData>& edata, int ind=2, int pre=10, int write_graph_mesh_=0, int write_update_time_=0, int write_ele_updatetime=0);
    GmshRecorder();
    ~GmshRecorder();

    int record(int commitTag, double timeStamp);
    int restart();
    int flush();
    int domainChanged();    
    int setDomain(Domain &theDomain);
    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

protected:

    virtual int write_mesh();
    virtual int write_header();
    virtual int write_node_data();
    virtual int write_element_data();
    virtual int write_element_graph();
    virtual int write_update_time_now();
    virtual int write_eleupdatetime_now();
    virtual int write_data_line(std::ofstream& s, const Vector & data, int truncatesize=3);

private:
    int  precision;
    bool write_header_now;
    bool write_mesh_now;
    bool write_binary_mode;
    bool write_ele_updatetime;
    std::string filename;
    std::vector<double> timestep;
    std::vector<ID> timeparts;
    std::ofstream theFile;
    NodeData nodedata;
    std::vector<EleData> eledata;
    Domain* theDomain;
    int current_step;
    int write_graph_mesh;
    int write_update_time;

public:
    enum GmshType {
	GMSH_VERTEX=15,GMSH_POLY_VERTEX=-1,GMSH_LINE=1,GMSH_POLY_LINE=-1,
	GMSH_TRIANGLE=2,GMSH_TRIANGLE_STRIP=-1,GMSH_POLYGON=-1,
	GMSH_PIXEL=-1,GMSH_QUAD=3,GMSH_TETRA=4,GMSH_VOXEL=-1,
	GMSH_HEXAHEDRON=5,GMSH_WEDGE=-1,GMSH_PYRAMID=7,
	GMSH_QUADRATIC_EDGE=8,GMSH_QUADRATIC_TRIANGLE=9,
	GMSH_QUADRATIC_QUAD=10,GMSH_QUADRATIC_TETRA=11,
	GMSH_QUADRATIC_HEXAHEDRON=12
    };
    static void setGMSHType();

private:
    static std::map<int,GmshType> gmshtypes;

};

// GMSH element types
// line_2_node                 = 1
// triangle_3_node             = 2
// quadrangle_4_node           = 3
// tetrahedron_4_node          = 4
// hexahedron_8_node           = 5
// prism_6_node                = 6
// pyramid_5_node              = 7
// line_3_node                 = 8
// triangle_6_node             = 9
// quadrangle_9_node           = 10
// tetrahedron_10_node         = 11
// hexahedron_27_node          = 12
// prism_18_node               = 13
// pyramid_14_node             = 14
// point_1_node                = 15
// quadrangle_8_node           = 16
// hexahedron_20_node          = 17
// prism_15_node               = 18
// pyramid_13_node             = 19
// triangle_9_node_incomplete  = 20
// triangle_10_node            = 21
// triangle_12_node_incomplete = 22
// triangle_15_node            = 23
// triangle_15_node_incomplete = 24
// triangle_21_node            = 25
// edge_4_node                 = 26
// edge_5_node                 = 27
// edge_6_node                 = 28
// tetrahedron_20_node         = 29
// tetrahedron_35_node         = 30
// tetrahedron_56_node         = 31
// hexahedron_64_node          = 92
// hexahedron_125_node         = 93

#endif
