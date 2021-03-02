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


// ============================================================================
// 2019 By Jose Abell @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// 
// H5DRM provides a load pattern that implements reading H5DRM datasets
// (based on HDF5) containing ground motions useful in the context of
// soil-structure interaction modeling using the Domain Reduction method.
//
// Important assumptions and information:
//
//  + A brute force method for matching FEM domiain nodes to those of the 
//    DRM dataset is used. There is space for improvements here (A scale 
//    factor can be applied to the DRM dataset point coordinates to deal 
//    with non-matching units, DRM datasets should be consistently 
//    created using base SI units and the coordinate system convention).
//    Nodal matching is sufficiently fast for problems with less than
//    1M nodes, it only occurs once during analysis.
//  + Coordinate system for H5DRM datasets:   Z is up, X is east and Y is
//    north. Units should be SI base units always. The coordinate system
//    of the FEM model should match that of the DRM except, maybe, in
//    in the units of the point coordinates. 
//  + This H5DRM pattern currently only works for 3-DOF 8 node bricks in
//    the DRM layer the rest of the domain can be anything.  
//  + Althounh the H5DRM format can store ground motions based on displacements, 
//    velocities or accelerations. The present implementation can use 
//    displacement-onlyH5DRM datasets (will double-differentiate internaly) or 
//    displacement and acceleration datasets (where the data is used directly.
//    use this if filterning of motions is handled independently for displacements 
//    and accelerations, although datasets are heavier). There is an untested
//    implementation that is based on velocities and will integrate and 
//    differentiate motions based on that to derive displacements and 
//    accelerations needed in classic DRM. 
//  + Linear interpolation is done on displacements and accelerations if
//    current analysis timestep does not "fall" in any of the motion time 
//    samples. This can be improved since a linear interpolation in 
//    accelerations is a cubic interpolation in displacements.  
//  + Based on the original papers by Jacobo Bielak [1] and [2]. 
//  + The previous remark means that velocities are not directly used. 
//    This might be important if damping in the DRM layer is thought to be 
//    needed. Please let me (jaabell) know if you think this is needed in
//    some case, I would be happy to discuss and implement this and also
//    know why it might be needed. (coupled poroelasticity?)
//  + It has been tested to work in parallel. 
//  + The H5DRM data format specification is documented in [3].
//  + Based on and and extends the work done for my PhD thesis [4].
// 
//  References: 
//   [1] Bielak, J., Loukakis, K., Hisada, Y., & Yoshimura, C. (2003). 
//       Domain Reduction Method for Three--Dimensional Earthquake 
//       Modeling in Localized Regions. Part {I}: Theory. Bulletin of the 
//       Seismological Society of America, 93(2), 817–824.
//   [2] Yoshimura, C., Bielak, J., & Hisada, Y. (2003). Domain Reduction 
//       Method for Three--Dimensional Earthquake Modeling in Localized 
//       Regions. Part {II}: Verif\/ication and Examples. Bulletin of the 
//       Seismological Society of America, 93(2), 825–840.
//   [3] Abell J.A., Crempien J.G.F and Recabarren M. (2019) ShakerMaker|H5DRM: 
//       a library for the simulation ofseismic wave-field data enabling 
//       high-fidelity earthquake safety modeling, and a format to store them. 
//       (To be submitted 2019 to SoftwareX)
//   [4] Abell, J. A. (2016). Earthquake-Soil-Structure Interaction Modeling 
//       of Nuclear Power Plants for Near-Field Events. 
//       PhD Thesis. University of California, Davis.
// 
// ============================================================================

#ifndef H5DRM_h
#define H5DRM_h

#include <LoadPattern.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>
#include <Domain.h>
#include <NodeIter.h>
#include <Node.h>
#include <ElementIter.h>
#include <Element.h>
#include <Channel.h>
#include <ErrorHandler.h>
#include <hdf5.h>
#include <map>
#include <vector>
#include <algorithm>  // For std::min and std::max functions
#include <string>

#define H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS 50
#define H5DRM_MAX_RETURN_OPEN_OBJS 100
#define H5DRM_MAX_FILENAME 200
#define H5DRM_MAX_STRINGSIZE 80
#define H5DRMerror std::cerr  << "H5DRM(" << myrank << ") - ERROR: "
#define H5DRMwarning std::cout  <<  "H5DRM(" << myrank << ") - Warning: "
#define H5DRMout std::cout  << "H5DRM(" << myrank << ") : "

class Vector;
class Matrix;

class Plane
{
public:
    Plane(const hid_t& id_h5drm_file, int plane_number, double crd_scale = 1);
    ~Plane();

    bool locate_point(const Vector& x, double& xi1, double& xi2, double& distance) const ;
    
    bool get_ij_coordinates_to_point(double xi1, double xi2, int& i, int& j) const;

    int getNumber() const;
    
    int getStationNumber(int i, int j) const;

    void print(FILE * fptr) const;
    
private:
    int number;
    bool internal;
    int *stations;
    Vector v0;
    Vector v1;
    Vector v2;
    Vector xi1;
    Vector xi2;
};



class H5DRM : public LoadPattern
{
public:
    H5DRM();
    H5DRM(int tag, std::string HDF5filename_, double cFactor_ = 1.0, double crd_scale_ = 1, double distance_tolerance_ = 1e-3);
    ~H5DRM();
    void clean_all_data(); // Called by destructor and if domain changes

    void setDomain(Domain *theDomain);
    void applyLoad(double time);
    void Print(ostream &s, int flag = 0);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
                    FEM_ObjectBroker &theBroker);
    LoadPattern *getCopy(void);
    virtual std::string getName() const
    {
        return "H5DRM: " + HDF5filename;
    }

protected:

    bool  ComputeDRMLoads(double t);
    bool  ComputeDRMMotions(double next_integration_time);
    bool  drm_differentiate_displacements(double next_integration_time);
    bool  drm_integrate_velocity(double next_integration_time);
    bool  drm_direct_read(double next_integration_time);
    Vector *getNodalLoad(int node, double time);

    void intitialize();

    void node_matching_BruteForce(double d_tol, const ID& internal, const Matrix& xyz, const Vector& drmbox_x0, double& d_err, int & n_nodes_found);
    void node_matching_UsePlaneInfo(double d_tol, const ID& internal, const Matrix& xyz, const Vector& drmbox_x0, double& d_err, int & n_nodes_found);
    void node_matching_OctTree(double d_tol, const ID& internal, const Matrix& xyz, const Vector& drmbox_x0, double& d_err, int & n_nodes_found);

private:

    std::string HDF5filename;   // Name of the HDF5 dataset containing the DRM motions  
                
    ID Elements;                // integer list containing the DRM elements
    ID Nodes;                   // integer list containing the nodes belonging to the DRM elements
    ID IsBoundary;              // integer list containing a flag indicating whether that node is a boundary node

    Vector DRMForces;           // vector containing the DRM equivalent-loads
    Vector DRMDisplacements;    // vector containing the displacements at DRM boundary nodes
    Vector DRMDisplacements0;   // vector containing the initial displacements at DRM boundary nodes
    Vector DRMAccelerations;    // vector containing the accelerations at DRM boundary nodes
    
    double last_integration_time;
    double t1, t2, tstart, tend, dt;    // specifies the time increment used in load path vector
    double cFactor;                     // additional factor on the returned load factor
    double crd_scale;                   // Scaling for the point coordinates of the DRM dataset
    double distance_tolerance;          // Distance tolerance for node-matching algorithm
    int step, step1, step2;
    double drmbox_xmax, drmbox_xmin, drmbox_ymax, drmbox_ymin, drmbox_zmax, drmbox_zmin;

    int number_of_DRM_elements;
    int number_of_DRM_nodes;
    int number_of_DRM_nodes_b;
    int number_of_DRM_nodes_e;
    int number_of_timesteps;
    int thetimeSteps;

    int maxnodetag;

    bool is_initialized;

    int number_of_local_elements;

    std::map<int, int> nodetag2station_id;
    std::map<int, int> nodetag2local_pos;
    ID station_id2data_pos;

    //HDF5 handles
    hid_t id_drm_file;
    hid_t id_velocity;
    hid_t id_velocity_dataspace;
    hid_t id_displacement;
    hid_t id_displacement_dataspace;
    hid_t id_acceleration;
    hid_t id_acceleration_dataspace;
    hid_t id_one_node_memspace;
    hid_t id_xfer_plist;

    int myrank;         // MPI Process-id (rank) in the case of parallel processing

    std::vector<Plane*> planes;
};

#endif
