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
// Please read detailed description in H5DRM.h.
// ============================================================================

#ifdef _H5DRM

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>
#include <time.h>
#include <hdf5.h>
#include <H5DRMLoadPattern.h>

#include <Channel.h>
#include <elementAPI.h>

#include <OPS_Globals.h>
#include <map>
#include <vector>
#include <Message.h>

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
#include <mpi.h>
#endif

#include <limits>   //For std::numeric_limits<double>::epsilon()

#include <Node.h>
#include <NodeIter.h>

#include <Element.h>
#include <ElementIter.h>


#define numNodeDOF 3  // Only nodes with 3-dofs per node can be used in DRM... :/

using namespace std;

#define DEBUG_NODE_MATCHING false
#define DEBUG_DRM_INTEGRATION false
#define DEBUG_DRM_FORCES false
#define DEBUG_WITH_GMSH false


Vector H5DRM_calculate_cross_product(const Vector& a, const Vector& b);
bool read_int_dataset_into_id(const hid_t& h5drm_dataset, std::string dataset_name, ID& result);
bool read_scalar_double_dataset_into_double(const hid_t& h5drm_dataset, std::string dataset_name, double& result);
bool read_double_dataset_into_vector(const hid_t& h5drm_dataset, std::string dataset_name, Vector& result);
bool read_double_dataset_into_matrix(const hid_t& h5drm_dataset, std::string dataset_name, Matrix& result);
bool read_int_dataset_into_array(const hid_t& h5drm_dataset, std::string dataset_name, int *& result);


#define HDF5_CLOSE_AND_REPORT(id, close_func, msg)  \
    if ((id) > 0 && (close_func(id)) < 0)           \
    {                                               \
        H5DRMerror << "Error closing " << msg << "." << endl; \
    }


static int numH5DRMpatterns = 0;

void* OPS_H5DRMLoadPattern()
{

    if (numH5DRMpatterns == 0)
    {
        opserr << "H5DRM - An HDF5 based DRM implementation. \n"
               << "          By:  Jose Abell (Prof. Universidad de los Andes, Chile) \n";
    }
    numH5DRMpatterns++;

    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "insufficient number of args\n";
        return 0;
    }

    LoadPattern *thePattern = 0;

    int num = 1;
    int tag = 0;
    OPS_GetIntInput(&num, &tag);

    std::string filename = OPS_GetString();

    double factor = 1.0;
    OPS_GetDoubleInput(&num, &factor);

    double crd_scale = 1.0;
    double distance_tolerance = 1e-3;
    bool do_coordinate_transformation = false;
    double T00 = 1.0; double T01 = 0.0; double T02 = 0.0;
    double T10 = 0.0; double T11 = 1.0; double T12 = 0.0;
    double T20 = 0.0; double T21 = 0.0; double T22 = 1.0;
    double x00 = 0.0; double x01 = 0.0; double x02 = 0.0;

    if (OPS_GetNumRemainingInputArgs() > 1)
    {
        OPS_GetDoubleInput(&num, &crd_scale);
        opserr << "crd_scale = " << crd_scale << endln;
    }

    if (OPS_GetNumRemainingInputArgs() > 1)
    {
        OPS_GetDoubleInput(&num, &distance_tolerance);
        opserr << "distance_tolerance = " << distance_tolerance << endln;
    }

    int int_do_coordinate_transformation;
    if (OPS_GetNumRemainingInputArgs() > 1)
    {
        OPS_GetIntInput(&num, &int_do_coordinate_transformation);
        do_coordinate_transformation = int_do_coordinate_transformation;
        opserr << "do_coordinate_transformation = " << do_coordinate_transformation << endln;
    }

    if (OPS_GetNumRemainingInputArgs() == 12)
    {
        OPS_GetDoubleInput(&num, &T00);
        OPS_GetDoubleInput(&num, &T01);
        OPS_GetDoubleInput(&num, &T02);
        OPS_GetDoubleInput(&num, &T10);
        OPS_GetDoubleInput(&num, &T11);
        OPS_GetDoubleInput(&num, &T12);
        OPS_GetDoubleInput(&num, &T20);
        OPS_GetDoubleInput(&num, &T21);
        OPS_GetDoubleInput(&num, &T22);
        OPS_GetDoubleInput(&num, &x00);
        OPS_GetDoubleInput(&num, &x01);
        OPS_GetDoubleInput(&num, &x02);
        opserr << "T = " <<  endln;
        opserr << "    " << T00 << " " << T01 << " " << T02 << endln;
        opserr << "    " << T10 << " " << T11 << " " << T12 << endln;
        opserr << "    " << T20 << " " << T21 << " " << T22 << endln;
        opserr << "x0 = " <<  endln;
        opserr << "     " << x00 << " " << x01 << " " << x02 << endln;
    }



    thePattern = new H5DRMLoadPattern(tag, filename,
                           factor,
                           crd_scale,
                           distance_tolerance,
                           do_coordinate_transformation,
                           T00, T01, T02,
                           T10, T11, T12,
                           T20, T21, T22,
                           x00, x01, x02);



    return thePattern;
}


















H5DRMLoadPattern::H5DRMLoadPattern()
    : LoadPattern(0, PATTERN_TAG_H5DRM),
      DRM_F(100),
      DRM_D(100),
      DRM_A(100),
      last_integration_time(0),
      t1(0), t2(0), tend(0),
      cFactor(1),
      crd_scale(0),
      distance_tolerance(0),
      flag_initialized(false),
      station_id2data_pos(100),
      ih5_fname(0), ih5_vel(0), ih5_vel_ds(0), ih5_dis(0), ih5_dis_ds(0), ih5_acc(0), ih5_acc_ds(0), ih5_one_node_ms(0), ih5_xfer_plist(0),       
      do_coordinate_transformation(true),
      T(3, 3),
      x0(3)
{

    T(0, 0) = 1;
    T(0, 1) = 0;
    T(0, 2) = 0;
    T(1, 0) = 0;
    T(1, 1) = 1;
    T(1, 2) = 0;
    T(2, 0) = 0;
    T(2, 1) = 0;
    T(2, 2) = 1;
    x0(0) = 0;
    x0(1) = 0;
    x0(2) = 0;


    MPI_local_rank = 0;
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_local_rank);
#endif
    H5DRMout << "H5DRMLoadPattern - empty constructor\n";
}



H5DRMLoadPattern::H5DRMLoadPattern(
    int tag,
    std::string dataset_fname_,
    double cFactor_,
    double crd_scale_,
    double distance_tolerance_,
    bool do_coordinate_transformation_,
    double T00, double T01, double T02,
    double T10, double T11, double T12,
    double T20, double T21, double T22,
    double x00, double x01, double x02
)
    : LoadPattern(tag, PATTERN_TAG_H5DRM),
      dataset_fname(dataset_fname_),
      DRM_F(100),
      DRM_D(100),
      DRM_A(100),
      last_integration_time(0),
      t1(0), t2(0), tend(0),
      cFactor(cFactor_),
      crd_scale(crd_scale_),
      distance_tolerance(distance_tolerance_),
      flag_initialized(false),
      station_id2data_pos(100),
      ih5_fname(0), ih5_vel(0), ih5_vel_ds(0), ih5_dis(0), ih5_dis_ds(0), ih5_acc(0), ih5_acc_ds(0), ih5_one_node_ms(0), ih5_xfer_plist(0),
      do_coordinate_transformation(do_coordinate_transformation_),
      T(3, 3),
      x0(3)
{

    MPI_local_rank = 0;
#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    MPI_Comm_rank(MPI_COMM_WORLD, &MPI_local_rank);
#endif


    T(0, 0) = T00;
    T(0, 1) = T01;
    T(0, 2) = T02;
    T(1, 0) = T10;
    T(1, 1) = T11;
    T(1, 2) = T12;
    T(2, 0) = T20;
    T(2, 1) = T21;
    T(2, 2) = T22;
    x0(0) = x00;
    x0(1) = x01;
    x0(2) = x02;


    if (numH5DRMpatterns == 0 && MPI_local_rank == 0)
    {
        opserr << " \n"
               << "---------------------------------------------------------------------------\n"
               << "H5DRM - Domain reduction method load pattern. \n"
               << "       \n"
               << "      By:  Jose Abell (Prof. Universidad de los Andes, Chile) \n"
               << "       \n"
               << "      Bug reports to OpenSees github repo and use @jaabell to get my\n"
               << "      attention. Send input files and DRM datasets through dropbox or\n"
               << "      similar service.  \n"
               << "       \n"
               << "      www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl \n"
               << "---------------------------------------------------------------------------\n";
    }
    numH5DRMpatterns++;
}

std::string vector_to_string(Vector v)
{
    std::stringstream ss;
    ss << "( ";
    for (int i = 0; i < v.Size(); ++i)
    {
        ss << v[i] << " ";
    }
    ss << ")";
    string s = ss.str();

    return s;
}


void H5DRMLoadPattern::do_intitialization()
{
    Domain *theDomain = this->getDomain();

    if (theDomain == 0)
    {
        return;
    }

    if (MPI_local_rank == 0)
    {
        H5DRMout << "initializing - filename : " << dataset_fname << "\n";
    }

    //===========================================================================
    // Open the specified file for read only
    //===========================================================================
    hid_t file_access_plist   = H5Pcreate(H5P_FILE_ACCESS);


    ih5_fname = H5Fopen(dataset_fname.c_str(), H5F_ACC_RDONLY, file_access_plist);

    if (ih5_fname < 0)
    {
        H5DRMerror << "Error opening HDF5 file: " << dataset_fname << endl;
        return;
    }


    //===========================================================================
    // MAP DRM data to local domain
    //===========================================================================

    H5DRMout << "Initialize" << endln;

    ID internal;
    Matrix xyz;
    Vector drmbox_x0;
    read_double_dataset_into_vector(ih5_fname, "DRM_Metadata/drmbox_x0", drmbox_x0);
    read_double_dataset_into_matrix(ih5_fname, "DRM_Data/xyz", xyz);
    read_int_dataset_into_id(ih5_fname, "DRM_Data/internal", internal);
    read_int_dataset_into_id(ih5_fname, "DRM_Data/data_location", station_id2data_pos);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_xmax", drmbox_xmax);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_xmin", drmbox_xmin);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_ymax", drmbox_ymax);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_ymin", drmbox_ymin);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_zmax", drmbox_zmax);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/drmbox_zmin", drmbox_zmin);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/dt", dt);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/tstart", tstart);
    read_scalar_double_dataset_into_double(ih5_fname, "DRM_Metadata/tend", tend);

    N_timesteps = int((tend - tstart) / dt);

    last_integration_time = theDomain->getCurrentTime();

    int NDRM_points = xyz.noRows();

    std::string transformation_status("(not transformed)");

    if (do_coordinate_transformation)
    {
        transformation_status = "(transformed)";

        if (MPI_local_rank == 0)
        {
            H5DRMout << "Applying coordinate transformation  xyz (new) = T xyz (old) + x0  with" << endl;
            H5DRMout << "T = " << "(" << T(0, 0) << " " << T(0, 1) << " " << T(0, 2) << ") (" << T(1, 0) << " " << T(1, 1) << " " << T(1, 2) << ") (" << T(2, 0) << " " << T(2, 1) << " " << T(2, 2) << ")" << endln;
            H5DRMout << "x0 = " << vector_to_string(x0) << endln;
            H5DRMout << "crd_scale = " << crd_scale << endln;
            H5DRMout << "NDRM_points = " << NDRM_points << endln;

            H5DRMout << "drmbox_x0 = " << vector_to_string(drmbox_x0) << " (before transformation)" << endln;
        }
        // utility for point transformation
        auto lam_transform = [this, &drmbox_x0](Vector & point) {
            static Vector copy(3);
            copy = point;
            copy -= drmbox_x0;
            copy *= crd_scale;
            point.addMatrixVector(0.0, T, copy, 1.0);
            point += x0;
        };

        // transform all points
        for (int i = 0; i < NDRM_points; ++i)
        {
            static Vector row(3);
            row(0) = xyz(i, 0);
            row(1) = xyz(i, 1);
            row(2) = xyz(i, 2);
            lam_transform(row);
            xyz(i, 0) = row(0);
            xyz(i, 1) = row(1);
            xyz(i, 2) = row(2);
        }

        // transform bounding box
        static Vector pmax(3);
        static Vector pmin(3);
        pmax(0) = drmbox_xmax;
        pmax(1) = drmbox_ymax;
        pmax(2) = drmbox_zmax;
        pmin(0) = drmbox_xmin;
        pmin(1) = drmbox_ymin;
        pmin(2) = drmbox_zmin;
        lam_transform(pmax);
        lam_transform(pmin);
        drmbox_xmax = pmax(0);
        drmbox_ymax = pmax(1);
        drmbox_zmax = pmax(2);
        drmbox_xmin = pmin(0);
        drmbox_ymin = pmin(1);
        drmbox_zmin = pmin(2);

        // transform drm x0, which is not the user-defined x0
        drmbox_x0 = x0;

        if (MPI_local_rank == 0)
        {
            H5DRMout << "drmbox_x0 = " << vector_to_string(drmbox_x0) << " (after transformation)" << endln;
        }
    }

    double d_tol = distance_tolerance;
    double d_err = 0;
    int n_nodes_found = 0;

    node_matching_BruteForce(d_tol, internal, xyz, drmbox_x0, d_err, n_nodes_found );


    if (MPI_local_rank == 0)
    {
        H5DRMout << "Dataset has " << NDRM_points << " data-points\n";
        H5DRMout << "drmbox_x0   =  " << vector_to_string(drmbox_x0) << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_xmax =  " << drmbox_xmax << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_xmin =  " << drmbox_xmin << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_ymax =  " << drmbox_ymax << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_ymin =  " << drmbox_ymin << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_zmax =  " << drmbox_zmax << " " << transformation_status.c_str() << "\n";
        H5DRMout << "drmbox_zmin =  " << drmbox_zmin << " " << transformation_status.c_str() << "\n";
        H5DRMout << "DRM BBox X  =  " << drmbox_xmax - drmbox_xmin << " " << transformation_status.c_str() << "\n";
        H5DRMout << "DRM BBox Y  =  " << drmbox_ymax - drmbox_ymin << " " << transformation_status.c_str() << "\n";
        H5DRMout << "DRM BBox Z  =  " << drmbox_zmax - drmbox_zmin << " " << transformation_status.c_str() << "\n";
    }

    if (DEBUG_WITH_GMSH)
    {
        if (MPI_local_rank == 0 )
        {
            char viewstations[100];
            sprintf(viewstations, "stations.msh");
            FILE * fptr = fopen(viewstations, "w");
            fprintf(fptr, "$MeshFormat\n");
            fprintf(fptr, "2.2 0 8\n");
            fprintf(fptr, "$EndMeshFormat\n");
            fprintf(fptr, "$Nodes\n");
            fprintf(fptr, "%d\n", NDRM_points);
            for (int i = 0; i < NDRM_points; ++i)
            {
                fprintf(fptr, "%d %f %f %f\n", i + 1, xyz(i, 0), xyz(i, 1), xyz(i, 2));
            }
            fprintf(fptr, "$EndNodes\n");
            fprintf(fptr, "$Elements\n");
            fprintf(fptr, "%d\n", NDRM_points);
            for (int i = 0; i < NDRM_points; ++i)
            {
                fprintf(fptr, "%d 15 2 1 1 %d\n", i + 1, i + 1);
            }
            fprintf(fptr, "$EndElements\n");
            fclose(fptr);
        }

        char viewdrm[100];
        sprintf(viewdrm, "drm.%d.msh", MPI_local_rank);
        FILE * fptrdrm = fopen(viewdrm, "w");

        fprintf(fptrdrm, "$MeshFormat\n");
        fprintf(fptrdrm, "2.2 0 8\n");
        fprintf(fptrdrm, "$EndMeshFormat\n");
        fprintf(fptrdrm, "$Nodes\n");
        fprintf(fptrdrm, "%d\n", theDomain->getNumNodes());

        NodeIter& node_iter = theDomain->getNodes();
        Node* node_ptr = 0;
        int drmtag = NDRM_points;
        while ((node_ptr = node_iter()) != 0)
        {
            const Vector& xyz_node = node_ptr->getCrds();
            fprintf(fptrdrm, "%d %f %f %f\n", ++drmtag, xyz_node(0), xyz_node(1), xyz_node(2));
        }
        fprintf(fptrdrm, "$EndNodes\n");


        fprintf(fptrdrm, "$Elements\n");
        fprintf(fptrdrm, "%d\n", theDomain->getNumNodes());
        int ndrmpoint = drmtag - 1;
        drmtag = n_nodes_found + 1;
        for (int i = 0; i < ndrmpoint; ++i)
        {
            fprintf(fptrdrm, "%d 15 2 1 %d %d\n", drmtag, 5 + MPI_local_rank, drmtag);
            ++drmtag;
        }
        fprintf(fptrdrm, "$EndElements\n");


        fclose(fptrdrm);
    }

    H5DRMout << "Found and connected " << n_nodes_found << " of " << NDRM_points << " nodes (d_err = " << d_err << ")\n";
    DRM_F.resize(3 * n_nodes_found);
    DRM_D.resize(3 * n_nodes_found);
    DRM_A.resize(3 * n_nodes_found);

    DRM_F.Zero();
    DRM_D.Zero();
    DRM_A.Zero();

    char mshfilename[100];
    FILE * fptr = 0;
    if (DEBUG_WITH_GMSH)
    {
        sprintf(mshfilename, "drmnodes.%d.msh", MPI_local_rank);
        fptr = fopen(mshfilename, "w");
        fprintf(fptr, "$MeshFormat\n");
        fprintf(fptr, "2.2 0 8\n");
        fprintf(fptr, "$EndMeshFormat\n");
        fprintf(fptr, "$NodeData\n");
        fprintf(fptr, "    1\n");
        fprintf(fptr, "\"DRM Nodes\"\n");
        fprintf(fptr, "1\n");
        fprintf(fptr, "0.0000000000e+00\n");
        fprintf(fptr, "3\n");
        fprintf(fptr, "0\n");
        fprintf(fptr, "1\n");
        fprintf(fptr, "%d\n", DRM_Nodes.Size());
        for (int i = 0; i < DRM_Nodes.Size(); ++i)
        {
            fprintf(fptr, "%d %d\n", DRM_Nodes(i), 1 * DRM_Boundary_Flag(i) + 2 * (1 - DRM_Boundary_Flag(i)));
        }
        fprintf(fptr, "$EndNodeData\n");
        fclose(fptr);
    }


    // Identify Elements on DRM boundary
    Element * element_ptr = 0;
    ElementIter& element_iter = theDomain->getElements();
    while ((element_ptr = element_iter()) != 0)
    {
        int element_tag = element_ptr->getTag();
        const ID& element_nodelist = element_ptr->getExternalNodes();
        int nnodes = element_nodelist.Size();
        int count = 0;
        for (int nodenum = 0; nodenum < nnodes; ++nodenum)
        {
            int node_tag = element_nodelist(nodenum);
            if (DRM_Nodes.getLocation(node_tag) >= 0)
                count++;
        }
        if (count == nnodes)
            DRM_Elements.insert(element_tag);
    }

    N_local_elements = DRM_Elements.Size();
    H5DRMout << "Found  " << DRM_Elements.Size() << " elements on DRM boundary. \n";


    if (DEBUG_WITH_GMSH)
    {
        sprintf(mshfilename, "drmelements.%d.msh", MPI_local_rank);
        fptr = fopen(mshfilename, "w");
        fprintf(fptr, "$MeshFormat\n");
        fprintf(fptr, "2.2 0 8\n");
        fprintf(fptr, "$EndMeshFormat\n");
        fprintf(fptr, "$ElementData\n");
        fprintf(fptr, "    1\n");
        fprintf(fptr, "\"DRM Elements\"\n");
        fprintf(fptr, "1\n");
        fprintf(fptr, "0.0000000000e+00\n");
        fprintf(fptr, "3\n");
        fprintf(fptr, "0\n");
        fprintf(fptr, "1\n");
        fprintf(fptr, "%d\n", DRM_Elements.Size());
        for (int i = 0; i < DRM_Elements.Size(); ++i)
        {
            fprintf(fptr, "%d %d\n", DRM_Elements(i), 1);
        }
        fprintf(fptr, "$EndElementData\n");
        fclose(fptr);

        sprintf(mshfilename, "drmforces.%d.msh", MPI_local_rank);
        fptr = fopen(mshfilename, "w");
        fprintf(fptr, "$MeshFormat\n");
        fprintf(fptr, "2.2 0 8\n");
        fprintf(fptr, "$EndMeshFormat\n");
        fclose(fptr);

        sprintf(mshfilename, "drmdisplacements.%d.msh", MPI_local_rank);
        fptr = fopen(mshfilename, "w");
        fprintf(fptr, "$MeshFormat\n");
        fprintf(fptr, "2.2 0 8\n");
        fprintf(fptr, "$EndMeshFormat\n");
        fclose(fptr);

        sprintf(mshfilename, "drmaccelerations.%d.msh", MPI_local_rank);
        fptr = fopen(mshfilename, "w");
        fprintf(fptr, "$MeshFormat\n");
        fprintf(fptr, "2.2 0 8\n");
        fprintf(fptr, "$EndMeshFormat\n");
        fclose(fptr);
    }


    //===========================================================================
    // Open Displacements and Accelerations arrays for reading
    //===========================================================================
    htri_t found_it = 0;


    found_it = H5Lexists(ih5_fname, "DRM_Data/velocity", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening velocity dataset.\n";


        ih5_vel = H5Dopen(ih5_fname, "DRM_Data/velocity", H5P_DEFAULT);

        if (ih5_vel < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/velocity\" array!\n ";
            return;
        }

        ih5_vel_ds = H5Dget_space(ih5_vel);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        ih5_one_node_ms  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    found_it = H5Lexists(ih5_fname, "DRM_Data/acceleration", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening acceleration dataset.\n";


        ih5_acc = H5Dopen(ih5_fname, "DRM_Data/acceleration", H5P_DEFAULT);

        if (ih5_acc < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/acceleration\" array!\n ";
            return;
        }

        ih5_acc_ds = H5Dget_space(ih5_acc);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        ih5_one_node_ms  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    found_it = H5Lexists(ih5_fname, "DRM_Data/displacement", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening displacement dataset.\n";


        ih5_dis = H5Dopen(ih5_fname, "DRM_Data/displacement", H5P_DEFAULT);

        if (ih5_dis < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/displacement\" array!\n ";
            return;
        }

        ih5_dis_ds = H5Dget_space(ih5_dis);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_PREALLOC_TSTEPS};
        ih5_one_node_ms  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    ih5_xfer_plist = H5Pcreate( H5P_DATASET_XFER);


    //===========================================================================
    // Set status to initialized and ready to compute loads
    //===========================================================================
    if (MPI_local_rank == 0)
    {
        H5DRMout << "dt =   " << dt << " \n";
        H5DRMout << "N_timesteps =   " << N_timesteps << " \n";
        H5DRMout << "tstart =   " << tstart << " \n";
        H5DRMout << "tend =   " << tend << " \n";
        H5DRMout << "DRM initialization done!\n";
    }

    flag_initialized = true;
}



H5DRMLoadPattern::~H5DRMLoadPattern()
{
    cleanup();
}

void H5DRMLoadPattern::cleanup()
{
    // Clear mappings
    nodetag2station_id.clear();
    nodetag2local_pos.clear();

    // Close individual HDF5 resources
    HDF5_CLOSE_AND_REPORT(ih5_vel, H5Dclose, "motion dataset");
    HDF5_CLOSE_AND_REPORT(ih5_vel_ds, H5Sclose, "motion dataspace");
    HDF5_CLOSE_AND_REPORT(ih5_one_node_ms, H5Sclose, "one node memspace");

    // Get the number of open objects associated with the file
    ssize_t num_open_objs = H5Fget_obj_count(ih5_fname, H5F_OBJ_ALL);
    if (num_open_objs < 0) {
        printf("H5DRMLoadPattern::cleanup - Error retrieving open HDF5 objects.\n");
        return;
    }

    // If there are open objects, retrieve their IDs and close them
    if (num_open_objs > 0) {
        // Allocate memory to store object IDs
        hid_t *obj_ids = (hid_t *)malloc(num_open_objs * sizeof(hid_t));
        if (!obj_ids) {
            printf("H5DRMLoadPattern::cleanup - Memory allocation error.\n");
            return;
        }

        // Get the object IDs
        ssize_t num_retrieved = H5Fget_obj_ids(ih5_fname, H5F_OBJ_ALL, num_open_objs, obj_ids);
        if (num_retrieved < 0) {
            printf("H5DRMLoadPattern::cleanup - Error retrieving object IDs.\n");
            free(obj_ids);
            return;
        }

        // Close each retrieved object ID
        for (ssize_t i = 0; i < num_retrieved; i++) {
            H5Oclose(obj_ids[i]);
        }

        // Free allocated memory
        free(obj_ids);
    }

    // Close the HDF5 file itself
    HDF5_CLOSE_AND_REPORT(ih5_fname, H5Fclose, "H5 file: \"" + dataset_fname + "\"");

    flag_initialized = false;
}

void
H5DRMLoadPattern::setDomain(Domain * theDomain)
{
    this->LoadPattern::setDomain(theDomain);
}


void
H5DRMLoadPattern::applyLoad(double time)
{
    if (!flag_initialized)
    {
        do_intitialization();
    }

    //If we don't have elements, no point in doing anything :/
    if (N_local_elements <= 0)
    {
        return;
    }
    else
    {
        Domain *theDomain = this->getDomain();

        bool computedLoads = CalculateBoundaryForces(time);

        if (theDomain == 0 || DRM_Elements.Size() == 0 || !computedLoads)
        {
            opserr << "H5DRMLoadPattern::applyLoad -- Error! Failed to compute DRM loads at time = " << time << endln;
            return;
        }

        Node *theNode = 0;
        static Vector load(3);
        load.Zero();

        FILE* fptr_forces = 0;
        FILE* fptr_displ = 0;
        FILE* fptr_accel = 0;


        if (DEBUG_WITH_GMSH)
        {
            static int step = 0;
            char mshfilename[100];

            // drmforces
            sprintf(mshfilename, "drmforces.%d.msh", MPI_local_rank);
            fptr_forces = fopen(mshfilename, "a");
            // drmdisplacements
            sprintf(mshfilename, "drmdisplacements.%d.msh", MPI_local_rank);
            fptr_displ = fopen(mshfilename, "a");
            // drmaccelerations
            sprintf(mshfilename, "drmaccelerations.%d.msh", MPI_local_rank);
            fptr_accel = fopen(mshfilename, "a");

            fprintf(fptr_forces, "$NodeData\n");
            fprintf(fptr_forces, "    1\n");
            fprintf(fptr_forces, "\"DRM Forces\"\n");
            fprintf(fptr_forces, "1\n");
            fprintf(fptr_forces, "%f\n", time);
            fprintf(fptr_forces, "3\n");
            fprintf(fptr_forces, "%d\n", step);
            fprintf(fptr_forces, "3\n");
            fprintf(fptr_forces, "%d\n", DRM_Nodes.Size());

            fprintf(fptr_displ, "$NodeData\n");
            fprintf(fptr_displ, "    1\n");
            fprintf(fptr_displ, "\"DRM Displacements\"\n");
            fprintf(fptr_displ, "1\n");
            fprintf(fptr_displ, "%f\n", time);
            fprintf(fptr_displ, "3\n");
            fprintf(fptr_displ, "%d\n", step);
            fprintf(fptr_displ, "3\n");
            fprintf(fptr_displ, "%d\n", DRM_Nodes.Size());

            fprintf(fptr_accel, "$NodeData\n");
            fprintf(fptr_accel, "    1\n");
            fprintf(fptr_accel, "\"DRM Acceleration\"\n");
            fprintf(fptr_accel, "1\n");
            fprintf(fptr_accel, "%f\n", time);
            fprintf(fptr_accel, "3\n");
            fprintf(fptr_accel, "%d\n", step);
            fprintf(fptr_accel, "3\n");
            fprintf(fptr_accel, "%d\n", DRM_Nodes.Size());


            ++step;
        }

        //Apply DRM forces to nodes
        for (int i = 0; i < DRM_Nodes.Size(); i++)
        {
            int nodeTag = DRM_Nodes[i] ;
            int local_pos = nodetag2local_pos[nodeTag];
            theNode = theDomain->getNode( nodeTag );
            if ( theNode == 0 )
                continue;

            load(0) = DRM_F(3 * local_pos + 0);
            load(1) = DRM_F(3 * local_pos + 1);
            load(2) = DRM_F(3 * local_pos + 2);


            if (DEBUG_WITH_GMSH)
            {
                fprintf(fptr_forces, "%d %f %f %f\n", nodeTag, load(0), load(1), load(2));
                fprintf(fptr_displ, "%d %f %f %f\n", nodeTag, DRM_D(3 * local_pos + 0), DRM_D(3 * local_pos + 1), DRM_D(3 * local_pos + 2));
                fprintf(fptr_accel, "%d %f %f %f\n", nodeTag, DRM_A(3 * local_pos + 0), DRM_A(3 * local_pos + 1), DRM_A(3 * local_pos + 2));
            }
            //Add to current nodal unbalanced load
            theNode->addUnbalancedLoad(load);
        }


        if (DEBUG_WITH_GMSH)
        {
            fprintf(fptr_forces, "$EndNodeData\n");
            fclose(fptr_forces);
            fprintf(fptr_displ, "$EndNodeData\n");
            fclose(fptr_displ);
            fprintf(fptr_accel, "$EndNodeData\n");
            fclose(fptr_accel);
        }
    }
}


bool H5DRMLoadPattern::ComputeDRMMotions(double next_integration_time)
{
    bool have_displacement = ih5_dis > 0;
    bool have_acceleration = ih5_acc > 0;

    if (!have_displacement)
    {
        opserr << "Error - H5DRM file " << dataset_fname.c_str() << " does not have a displacements dataset. FAILING" << endln;
        return false;
    }

    if (!have_acceleration)
    {
        opserr << "Error - H5DRM file " << dataset_fname.c_str() << " does not have an acceleration dataset. FAILING" << endln;
        return false;
    }

    if (have_displacement && have_acceleration)
    {
        return drm_direct_read(next_integration_time);
    }

    return false;
}

bool H5DRMLoadPattern::drm_direct_read(double t)
{

    if (DRM_Nodes.Size() == 0)
    {
        H5DRMout << " This process has no DRM nodes. Nothing to be done by H5DRM" << endln;
        return false;
    }

    DRM_D.Zero();
    DRM_A.Zero();

    if (t < tstart || t > tend)
    {
        H5DRMout << "t = " << t << " tstart = " << tstart << " tend = " << tend << " DRM Not computing forces (t < tstart or t > tend)"  << endln;
        return true;
    }

    int i1 = (int) floor( (t - tstart) / dt);
    int i2 = (int) floor( (t - tstart) / dt) + 1;
    double t1 = i1 * dt + tstart;
    double t2 = i2 * dt + tstart;
    double dtau = (t - t1) / (t2 - t1);

    if (MPI_local_rank == 0)
    {
        H5DRMout << "t = " << t << " dt = " << dt << " i1 = " << i1 << " i2 = " << i2 << " t1 = " << t1 << " t2 = " << t2 << " dtau = " << dtau << endln;
    }

    ih5_acc_ds = H5Dget_space(ih5_acc);
    ih5_dis_ds = H5Dget_space(ih5_dis);

    double umax = -std::numeric_limits<double>::infinity();
    double amax = -std::numeric_limits<double>::infinity();
    double umin =  std::numeric_limits<double>::infinity();
    double amin =  std::numeric_limits<double>::infinity();

    for (int n = 0; n < DRM_Nodes.Size(); ++n)
    {
        int nodeTag = DRM_Nodes(n);
        int station_id = nodetag2station_id[nodeTag];
        int data_pos = station_id2data_pos[station_id];
        int local_pos = nodetag2local_pos[nodeTag];

        double d1[3], d2[3];
        double a1[3], a2[3];

        d1[0] = d1[1] = d1[2] = 0.;
        d2[0] = d2[1] = d2[2] = 0.;
        a1[0] = a1[1] = a1[2] = 0.;
        a2[0] = a2[1] = a2[2] = 0.;

        hsize_t start[2]  = {(hsize_t) data_pos , (hsize_t) i1};
        hsize_t stride[2] = {1                  , 1};
        hsize_t count[2]  = {3                  , 1};
        hsize_t block[2]  = {1                  , 1};

        //Selection in dataspace
        H5Sselect_hyperslab(
            ih5_dis_ds,
            H5S_SELECT_SET, start, stride, count, block );
        H5Sselect_hyperslab(
            ih5_acc_ds,
            H5S_SELECT_SET, start, stride, count, block );


        //Selection in memspace
        hsize_t rank_one_array = 1;
        hsize_t one_node_data_dims[1] = {3};
        hsize_t one_node_data_maxdims[1] = {3};
        hid_t memspace  = H5Screate_simple(rank_one_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory


        hsize_t mem_start[1] = {0};
        hsize_t mem_stride[1] = {1};
        hsize_t mem_count[1] = {3};
        hsize_t mem_block[1] = {1};
        H5Sselect_hyperslab(
            memspace,
            H5S_SELECT_SET, mem_start, mem_stride, mem_count, mem_block );

        //Read data
        herr_t errorflag1 = H5Dread( ih5_dis, H5T_NATIVE_DOUBLE, memspace,
                                     ih5_dis_ds, ih5_xfer_plist,  d1 );
        herr_t errorflag2 = H5Dread(ih5_acc, H5T_NATIVE_DOUBLE, memspace,
                                    ih5_acc_ds, ih5_xfer_plist,  a1 );
        start[1] = (hsize_t) i2;
        H5Sselect_hyperslab(
            ih5_dis_ds,
            H5S_SELECT_SET, start, stride, count, block );
        H5Sselect_hyperslab(
            ih5_acc_ds,
            H5S_SELECT_SET, start, stride, count, block );


        herr_t errorflag3 = H5Dread( ih5_dis, H5T_NATIVE_DOUBLE, memspace,
                                     ih5_dis_ds, ih5_xfer_plist,  d2 );
        herr_t errorflag4 = H5Dread(ih5_acc, H5T_NATIVE_DOUBLE, memspace,
                                    ih5_acc_ds, ih5_xfer_plist,  a2 );
        H5Sclose(memspace);

        bool nanfound = false;
        for (int i = 0; i < 3; ++i)
        {
            if ( isnan(d1[i])  ||
                    isnan(a1[i])  ||
                    isnan(d2[i])  ||
                    isnan(a2[i]) )
            {
                nanfound = true;
            }
            umax = d1[i] > umax ? d1[i] : umax;
            umax = d2[i] > umax ? d2[i] : umax;
            amax = a1[i] > amax ? a1[i] : amax;
            amax = a2[i] > amax ? a2[i] : amax;
            umin = d1[i] < umin ? d1[i] : umin;
            umin = d2[i] < umin ? d2[i] : umin;
            amin = a1[i] < amin ? a1[i] : amin;
            amin = a2[i] < amin ? a2[i] : amin;
        }


        if (errorflag1 < 0 || errorflag2 < 0 || errorflag3 < 0 || errorflag4 < 0 || nanfound)
        {
            H5DRMerror << "H5DRMLoadPattern::drm_direct_read - Failed to read displacement or acceleration array!!\n" <<
                       " n = " << n << endln <<
                       " nodeTag = " << nodeTag << endln <<
                       " station_id = " << station_id << endln <<
                       " i1 = " << i1 << endln <<
                       " data_pos = " << data_pos << endln <<
                       " local_pos = " << local_pos << endln <<
                       " Nt = " << 1 << endln <<
                       " last_integration_time = " << last_integration_time << endln <<
                       " start   = [" << start[0] << ", " << start[1] << "]" << endln <<
                       " mem_start  = [" << mem_start[0] << "]" << endln <<
                       " stride  = [" << stride[0] << ", " << stride[1] << "]" << endln <<
                       " count   = [" << count[0] << ", " << count[1] << "]" << endln <<
                       " block   = [" << block[0] << ", " << block[1] << "]" << endln <<
                       " errorflag1 = [" << errorflag1 << ", " << block[1] << "]" << endln <<
                       " errorflag2 = [" << errorflag2 << ", " << block[1] << "]" << endln <<
                       " errorflag3 = [" << errorflag3 << ", " << block[1] << "]" << endln <<
                       " errorflag4 = [" << errorflag4 << ", " << block[1] << "]" << endln <<
                       " nanfound   = [" << nanfound << ", " << block[1] << "]" << endln ;

            exit(-1);
        }

        d1[2] = -d1[2];
        d2[2] = -d2[2];
        a1[2] = -a1[2];
        a2[2] = -a2[2];


        DRM_D(3 * local_pos + 0) = d1[0] * (1 - dtau) + d2[0] * (dtau);
        DRM_D(3 * local_pos + 1) = d1[1] * (1 - dtau) + d2[1] * (dtau);
        DRM_D(3 * local_pos + 2) = d1[2] * (1 - dtau) + d2[2] * (dtau);

        DRM_A(3 * local_pos + 0) = a1[0] * (1 - dtau) + a2[0] * (dtau);
        DRM_A(3 * local_pos + 1) = a1[1] * (1 - dtau) + a2[1] * (dtau);
        DRM_A(3 * local_pos + 2) = a1[2] * (1 - dtau) + a2[2] * (dtau);
    }


    H5DRMout << "t = " << t << " u = (" << umin << ", " << umax << ") a = (" << amin << ", " << amax << ")" << endln;

    return true;
}


bool H5DRMLoadPattern::drm_differentiate_displacements(double t)
{

    if (DRM_Nodes.Size() == 0)
        return false;

    DRM_D.Zero();
    DRM_A.Zero();

    if (t < tstart || t > tend)
    {
        return true;
    }

    int i1 = (int) floor( (t - tstart) / dt);
    int i2 = (int) floor( (t - tstart) / dt) + 1;
    double t1 = i1 * dt + tstart;
    double t2 = i2 * dt + tstart;

    double dtau = (t - t1) / (t2 - t1);

    ih5_dis_ds = H5Dget_space(ih5_dis);
    hsize_t i_first = i1 - 1;
    hsize_t i_last  = i1 + 2;

    i_first = i_first < 0 ? 0 : i_first;
    i_first = i_first > (hsize_t) N_timesteps - 1 ? N_timesteps - 1 : i_first;

    i_last = i_last < 0 ? 0 : i_last;
    i_last = i_last > (hsize_t) N_timesteps - 1 ? N_timesteps - 1 : i_last;

    hsize_t i_len = (i_last - i_first) + 1;
    H5DRMout << "t = " << t
             << " dt = " << dt
             << " i1 = " << i1
             << " i2 = " << i2
             << " i_first = " << i_first
             << " i_last = " << i_last
             << " i_len = " << i_len
             << " t1 = " << t1
             << " t2 = " << t2
             << " dtau = " << dtau << endln;

    double umax = -std::numeric_limits<double>::infinity();
    double amax = -std::numeric_limits<double>::infinity();
    double umin =  std::numeric_limits<double>::infinity();
    double amin =  std::numeric_limits<double>::infinity();
    double dt2 = dt * dt;



    for (int n = 0; n < DRM_Nodes.Size(); ++n)
    {
        int nodeTag = DRM_Nodes(n);
        int station_id = nodetag2station_id[nodeTag];
        int data_pos = station_id2data_pos[station_id];
        int local_pos = nodetag2local_pos[nodeTag];

        double d0[3][4];
        double d1[3], d2[3];
        double a1[3], a2[3];

        d0[0][0] = d0[1][0] = d0[2][0] = 0.;
        d0[0][1] = d0[1][1] = d0[2][1] = 0.;
        d0[0][2] = d0[1][2] = d0[2][2] = 0.;
        d0[0][3] = d0[1][3] = d0[2][3] = 0.;
        d1[0] = d1[1] = d1[2] = 0.;
        d2[0] = d2[1] = d2[2] = 0.;
        a1[0] = a1[1] = a1[2] = 0.;
        a2[0] = a2[1] = a2[2] = 0.;

        hsize_t start[2]  = {(hsize_t) data_pos , (hsize_t) i_first};
        hsize_t stride[2] = {1                  , 1};
        hsize_t count[2]  = {3                  , i_len};
        hsize_t block[2]  = {1                  , 1};

        //Selection in dataspace
        H5Sselect_hyperslab(
            ih5_dis_ds,
            H5S_SELECT_SET, start, stride, count, block );

        //Selection in memspace
        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, i_len};
        hsize_t one_node_data_maxdims[2] = {3, i_len};
        hid_t memspace  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory


        hsize_t mem_start[2] = {0, 0};
        hsize_t mem_stride[2] = {1, 1};
        hsize_t mem_count[2] = {3, i_len};
        hsize_t mem_block[2] = {1, 1};
        H5Sselect_hyperslab(
            memspace,
            H5S_SELECT_SET, mem_start, mem_stride, mem_count, mem_block );

        //Read data
        herr_t errorflag1 = 0;
        if (i1 > 2)
        {
            errorflag1 = H5Dread( ih5_dis, H5T_NATIVE_DOUBLE, memspace,
                                  ih5_dis_ds, ih5_xfer_plist,  d0 );
        }


        for (int dof = 0; dof < 3; ++dof)
        {
            d1[dof] = d0[dof][1];
            d2[dof] = d0[dof][2];
            a1[dof] = (d0[dof][0] - 2 * d0[dof][1] + d0[dof][2]) / dt2;
            a2[dof] = (d0[dof][1] - 2 * d0[dof][2] + d0[dof][3]) / dt2;
        }


        H5Sclose(memspace);

        bool nanfound = false;
        for (int i = 0; i < 3; ++i)
        {

            if ( isnan(d1[i])  ||
                    isnan(a1[i])  ||
                    isnan(d2[i])  ||
                    isnan(a2[i]) )
            {
                nanfound = true;
            }
            umax = d1[i] > umax ? d1[i] : umax;
            umax = d2[i] > umax ? d2[i] : umax;
            amax = a1[i] > amax ? a1[i] : amax;
            amax = a2[i] > amax ? a2[i] : amax;
            umin = d1[i] < umin ? d1[i] : umin;
            umin = d2[i] < umin ? d2[i] : umin;
            amin = a1[i] < amin ? a1[i] : amin;
            amin = a2[i] < amin ? a2[i] : amin;
        }


        if (errorflag1 < 0 ||
                // errorflag2 < 0 || errorflag3 < 0 || errorflag4 < 0 ||
                nanfound)
        {
            H5DRMerror << "Failed to read displacement or acceleration array!!\n" <<
                       " n = " << n << endln <<
                       " nodeTag = " << nodeTag << endln <<
                       " station_id = " << station_id << endln <<
                       " i1 = " << i1 << endln <<
                       " data_pos = " << data_pos << endln <<
                       " local_pos = " << local_pos << endln <<
                       " Nt = " << 1 << endln <<
                       " last_integration_time = " << last_integration_time << endln <<
                       " start   = [" << start[0] << ", " << start[1] << "]" << endln <<
                       " mem_start  = [" << mem_start[0] << ", " << mem_start[1] << "]" << endln <<
                       " stride  = [" << stride[0] << ", " << stride[1] << "]" << endln <<
                       " count   = [" << count[0] << ", " << count[1] << "]" << endln <<
                       " block   = [" << block[0] << ", " << block[1] << "]" << endln;
            exit(-1);
        }

        d1[2] = -d1[2];
        d2[2] = -d2[2];
        a1[2] = -a1[2];
        a2[2] = -a2[2];


        DRM_D(3 * local_pos + 0) = d1[0] * (1 - dtau) + d2[0] * (dtau);
        DRM_D(3 * local_pos + 1) = d1[1] * (1 - dtau) + d2[1] * (dtau);
        DRM_D(3 * local_pos + 2) = d1[2] * (1 - dtau) + d2[2] * (dtau);

        DRM_A(3 * local_pos + 0) = a1[0] * (1 - dtau) + a2[0] * (dtau);
        DRM_A(3 * local_pos + 1) = a1[1] * (1 - dtau) + a2[1] * (dtau);
        DRM_A(3 * local_pos + 2) = a1[2] * (1 - dtau) + a2[2] * (dtau);
    }

    H5DRMout << "t = " << t << " u = (" << umin << ", " << umax << ") a = (" << amin << ", " << amax << ")" << endln;

    return true;
}

bool H5DRMLoadPattern::drm_integrate_velocity(double next_integration_time)
{
    // TODO: the idea here is that it should be enough to just store
    // DRM velocities then integrate to get displacement and differentiate
    // to get accelerations.
    exit(-1);
    return false;
}

bool H5DRMLoadPattern::CalculateBoundaryForces(double currentTime)
{
    if (MPI_local_rank == 0)
    {
        H5DRMout << "CalculateBoundaryForces @ t = " << currentTime << "\n";
    }

    constexpr int NDF = 3;
    constexpr int MaxNodes = 8;
    
    static Vector Peff_b(MaxNodes * NDF); 
    static Vector Peff_e(MaxNodes * NDF); 
    
    static Vector u_b(MaxNodes * NDF); 
    static Vector a_b(MaxNodes * NDF); 
    static Vector u_e(MaxNodes * NDF); 
    static Vector a_e(MaxNodes * NDF); 

    static Matrix M_be(MaxNodes * NDF, MaxNodes * NDF); 
    static Matrix K_be(MaxNodes * NDF, MaxNodes * NDF); 

    static ID BoundaryNodes(MaxNodes); 
    BoundaryNodes.Zero();
    
    static ID ExteriorNodes(MaxNodes); 
    ExteriorNodes.Zero();

    DRM_F.Zero();

    if (currentTime < tstart || currentTime > tend)
    {
        // Assume that DRM forces are null outside of the defined range...
        // If there are residual displacements, this will produce a 
        // sudden return to zero response at the end of analysis
        // better not to exceed the end of the analysis
        return true;
    }

    if (N_local_elements > 0)
    {
        if (!ComputeDRMMotions(currentTime))
        {
            return false;
        }

        Domain* theDomain = this->getDomain();
        Element* theElement = theDomain->getElement(DRM_Elements[0]);
        int numElementNodes = 8;
        
        BoundaryNodes.resize(MaxNodes);
        ExteriorNodes.resize(MaxNodes);

        for (int elemIndex = 0; elemIndex < DRM_Elements.Size(); ++elemIndex)
        {
            int elementTag = DRM_Elements[elemIndex];
            theElement = theDomain->getElement(elementTag);
            const ID& elementNodeIDs = theElement->getExternalNodes();
            numElementNodes = elementNodeIDs.Size();

            // Identify boundary and exterior nodes
            int boundaryCount = 0, exteriorCount = 0;
            for (int nodeIndex = 0; nodeIndex < numElementNodes; ++nodeIndex)
            {
                int nodeTag = elementNodeIDs(nodeIndex);
                int localPosition = nodetag2local_pos[nodeTag];

                if (DRM_Boundary_Flag[localPosition])
                {
                    BoundaryNodes(boundaryCount++) = nodeIndex;
                }
                else
                {
                    ExteriorNodes(exteriorCount++) = nodeIndex;
                }
            }

            // Compute DRM forces if there are identified nonzero b and e nodes
            if (boundaryCount != 0 && exteriorCount != 0)
            {
                // Get K_be and M_be matrices
                M_be.resize(NDF*boundaryCount, NDF*exteriorCount);
                M_be.Zero();
                K_be.resize(NDF*boundaryCount, NDF*exteriorCount);
                K_be.Zero();

                const Matrix& M_ele = theElement->getMass();
                const Matrix& K_ele = theElement->getTangentStiff();

                for (int i = 0; i < boundaryCount; ++i)
                {
                    int b = BoundaryNodes(i);
                    for (int j = 0; j < exteriorCount; ++j)
                    {
                        int e = ExteriorNodes(j);
                        for (int ci = 0; ci < NDF; ++ci)
                        {
                            for (int cj = 0; cj < NDF; ++cj)
                            {
                                M_be(NDF*i+ci,NDF*j+cj) = M_ele(NDF*b+ci,NDF*e+cj);
                                K_be(NDF*i+ci,NDF*j+cj) = K_ele(NDF*b+ci,NDF*e+cj);
                            }
                        }
                    }
                }

                // Get the nodal displacements and accelerations at b and e nodes
                u_b.resize(NDF*boundaryCount);
                u_b.Zero();
                a_b.resize(NDF*boundaryCount);
                a_b.Zero();
                u_e.resize(NDF*exteriorCount);
                u_e.Zero();
                a_e.resize(NDF*exteriorCount);
                a_e.Zero();

                for (int i = 0; i < boundaryCount; ++i)
                {
                    int b = BoundaryNodes(i);
                    int nodeTag = elementNodeIDs(b);
                    int localPosition = nodetag2local_pos[nodeTag];
                    for (int ci = 0; ci < 3; ++ci)
                    {
                        u_b(3 * i + ci) = DRM_D[3 * localPosition + ci];
                        a_b(3 * i + ci) = DRM_A[3 * localPosition + ci];
                    }
                }

                for (int j = 0; j < exteriorCount; ++j)
                {
                    int e = ExteriorNodes(j);
                    int nodeTag = elementNodeIDs(e);
                    int localPosition = nodetag2local_pos[nodeTag];
                    for (int cj = 0; cj < 3; ++cj)
                    {
                        u_e(3 * j + cj) = DRM_D[3 * localPosition + cj];
                        a_e(3 * j + cj) = DRM_A[3 * localPosition + cj];
                    }
                }

                // Peff_b = -Mbe ae - Kbe ue
                Peff_b.resize(NDF*boundaryCount);
                Peff_b.addMatrixVector(0, M_be, a_e, -1.0);
                Peff_b.addMatrixVector(1, K_be, u_e, -1.0);

                // Peff_e = Meb ab + Keb ub
                Peff_e.resize(NDF*exteriorCount);
                Peff_e.addMatrixTransposeVector(0, M_be, a_b, 1.0);
                Peff_e.addMatrixTransposeVector(1, K_be, u_b, 1.0);

                // Now add the DRM forces for this element to the global DRM forces
                for (int i = 0; i < boundaryCount; ++i)
                {
                    int b = BoundaryNodes(i);
                    int nodeTag = elementNodeIDs(b);
                    int localPosition = nodetag2local_pos[nodeTag];
                    for (int ci = 0; ci < 3; ++ci)
                    {
                        DRM_F[3 * localPosition + ci] += Peff_b[3*i+ci];
                    }
                }

                for (int j = 0; j < exteriorCount; ++j)
                {
                    int e = ExteriorNodes(j);
                    int nodeTag = elementNodeIDs(e);
                    int localPosition = nodetag2local_pos[nodeTag];
                    for (int cj = 0; cj < 3; ++cj)
                    {
                        DRM_F[3 * localPosition + cj] += Peff_e[3*j+cj];
                    }
                }

                // Check for NaN
                for (int k = 0; k < numElementNodes; ++k)
                {
                    int nodeTag = elementNodeIDs(k);
                    int localPosition = nodetag2local_pos[nodeTag];

                    if (isnan(DRM_F(3 * localPosition)) || isnan(DRM_F(3 * localPosition + 1)) || isnan(DRM_F(3 * localPosition + 2)))
                    {
                        H5DRMerror << "NAN Detected!!! \n";
                        H5DRMerror << "    nodeTag = " << nodeTag << endln;
                        H5DRMerror << "    localPosition = " << localPosition << endln;
                        
                        exit(-1);
                    }
                }
            }
        }
    }

    if (MPI_local_rank == 0)
    {
        H5DRMout << "CalculateBoundaryForces. END\n";
    }

    if (DEBUG_DRM_FORCES)
    {
        char debugFilename[100];
        sprintf(debugFilename, "drmforces.%d.txt", MPI_local_rank);
        FILE* filePointer = fopen(debugFilename, "w");

        fprintf(filePointer, "%f ", currentTime);

        for (int i = 0; i < DRM_F.Size(); ++i)
        {
            fprintf(filePointer, "%f ", DRM_F(i));
        }
        fprintf(filePointer, "\n");
    }

    return true;
}



int
H5DRMLoadPattern::sendSelf(int commitTag, Channel & theChannel)
{

    H5DRMout << "sending filename: " << dataset_fname << endl;

    static Vector data(3);
    data(0) = cFactor;
    data(1) = crd_scale;
    data(2) = distance_tolerance;

    char drmfilename[H5DRM_MAX_FILENAME];
    strcpy(drmfilename, dataset_fname.c_str());
    Message filename_msg(drmfilename, H5DRM_MAX_FILENAME);

    if (theChannel.sendMsg(0, 0, filename_msg) < 0)
    {
        cerr << "H5DRMLoadPattern::sendSelf - error sending filename_msg\n";
        return -1;
    }

    if (theChannel.sendVector(0, 0, data) < 0)
    {
        cerr << "H5DRMLoadPattern::sendSelf -- error sending data\n";
        return -1;
    }


    return 0;
}

int
H5DRMLoadPattern::recvSelf(int commitTag, Channel & theChannel,
                FEM_ObjectBroker & theBroker)
{
    static Vector data(3);
    char drmfilename[H5DRM_MAX_FILENAME];
    Message filename_msg(drmfilename, H5DRM_MAX_FILENAME);

    if (theChannel.recvMsg(0, 0, filename_msg) < 0)
    {
        cerr << "H5DRMLoadPattern::receiveSelf - error receiving filename_msg\n";
        return -1;
    }

    if (theChannel.recvVector(0, 0, data) < 0)
    {
        cerr << "H5DRMLoadPattern::receiveSelf - error receiving data\n";
        return -1;
    }
    cFactor = data(0);
    crd_scale = data(1);
    distance_tolerance = data(2);

    dataset_fname = drmfilename;

    return 0;
}

void
H5DRMLoadPattern::Print(ostream & s, int flag)
{

}


LoadPattern *
H5DRMLoadPattern::getCopy(void)
{
    return new H5DRMLoadPattern(this->getTag(), dataset_fname);
}



Vector *
H5DRMLoadPattern::getNodalLoad(int nodeTag, double time)
{


    return 0;
}








void H5DRMLoadPattern::node_matching_BruteForce(double d_tol, const ID & internal, const Matrix & xyz, const Vector & drmbox_x0, double & d_err, int & n_nodes_found)
{

    if (MPI_local_rank == 0)
    {
        H5DRMout << "node_matching_BruteForce - Begin! d_tol = " << d_tol << "\n";
        H5DRMout << "                                  drmbox_x0(0) = " << drmbox_x0(0) << "\n";
        H5DRMout << "                                  drmbox_x0(1) = " << drmbox_x0(1) << "\n";
        H5DRMout << "                                  drmbox_x0(2) = " << drmbox_x0(2) << "\n";
    }

    char debugfilename[100];
    sprintf(debugfilename, "debugdrmbruteforce.%d.txt", MPI_local_rank);
    FILE * fptrdrm;
    if (DEBUG_NODE_MATCHING)
    {
        fptrdrm = fopen(debugfilename, "w");
    }

    int Nstations = xyz.noRows();
    bool *accounted_for = new bool[Nstations+1];

    for (int i = 0; i <= Nstations; ++i)
    {
        accounted_for[i] = false;
    }

    Domain *theDomain = this->getDomain();
    NodeIter& node_iter = theDomain->getNodes();
    int NDRM_points = theDomain->getNumNodes();
    Node* node_ptr = 0;
    int drmtag = NDRM_points;
    int local_pos = 0;
    while ((node_ptr = node_iter()) != 0)
    {
        int tag = node_ptr->getTag();
        const Vector& node_xyz  =  node_ptr->getCrds();
        double dmin = std::numeric_limits<double>::infinity();
        int ii_station_min = 0;

        if (DEBUG_NODE_MATCHING)
        {
            fprintf(fptrdrm, "%d %f %f %f\n", ++drmtag, node_xyz[0] , node_xyz[1] , node_xyz[2] );
        }

        Vector station_xyz(3);
        for (int ii = 0; ii < xyz.noRows(); ++ii)
        {
            station_xyz(0) = xyz(ii, 0);
            station_xyz(1) = xyz(ii, 1);
            station_xyz(2) = xyz(ii, 2);
            double d = (node_xyz - station_xyz).Norm();
            if (d < dmin)
            {
                dmin = d;
                ii_station_min = ii;
            }
        }

        if (fabs(dmin) < d_tol)
        {
            int station_id = ii_station_min;
            static Vector station_xyz(3);
            for (int dir = 0; dir < 3; ++dir)
            {
                station_xyz(dir) = xyz(station_id, dir);
            }
            double this_d_err = (station_xyz - node_xyz ).Norm();
            if (this_d_err > d_err)
            {
                d_err = this_d_err;
            }
            nodetag2station_id.insert ( std::pair<int, int>(tag, station_id) );
            nodetag2local_pos.insert (  std::pair<int, int>(tag, local_pos) );
            DRM_Nodes[local_pos] = tag;
            DRM_Boundary_Flag[local_pos] = internal(station_id);

            accounted_for[station_id] = true;

            ++n_nodes_found;
            ++local_pos;
        }
        else
        {
            if (DEBUG_NODE_MATCHING)
            {
                fprintf(fptrdrm, "Node # %05d @ (%4.2f, %4.2f, %4.2f) rejected \n", tag, node_xyz(0) , node_xyz(1) , node_xyz(2) );
                fprintf(fptrdrm, "closest pal @ (%4.2f, %4.2f, %4.2f) ii =  %05d  with dmin = %f \n", xyz(ii_station_min, 0) , xyz(ii_station_min, 1) , xyz(ii_station_min, 1), ii_station_min , dmin);
            }
        }
    }

    if (DEBUG_NODE_MATCHING)
    {
        fclose(fptrdrm);
    }

    int n_accounted_for = 0;

#if defined(_PARALLEL_PROCESSING) || defined(_PARALLEL_INTERPRETERS)
    bool * accounted_for0 = new bool[Nstations];
    MPI_Reduce(
    accounted_for,
    accounted_for0,
    Nstations,
    MPI_C_BOOL,
    MPI_LOR,
    0,
    MPI_COMM_WORLD);


    for (int i = 0; i < Nstations; ++i)
    {
        if (accounted_for0[i])
        {
            n_accounted_for++;
        }
    }

    delete [] accounted_for;
    delete [] accounted_for0;

#else

    for (int i = 0; i < Nstations; ++i)
    {
        if (accounted_for[i])
        {
            n_accounted_for++;
        }
    }

    delete [] accounted_for;

#endif


    if (MPI_local_rank == 0)
    {
        H5DRMout << "node_matching_BruteForce - End!\n";
        H5DRMout << "Accounted for " << n_accounted_for << " out of " << Nstations << " stations\n";
    }

    return;
}














//================================================================================
//Utility functions and classes
//================================================================================



Vector H5DRM_calculate_cross_product(const Vector & a, const Vector & b)
{
    Vector a_cross_b(3); // Store the result here

    if ( (a.Size() != 3) || (b.Size() != 3) )
    {
        opserr << "Error: H5DRM_calculate_cross_product only defined for 3x1 vectors.\n";
        exit(-1);
    }

    a_cross_b(0) =   a(1) * b(2) - b(1) * a(2);
    a_cross_b(1) = - a(0) * b(2) + b(0) * a(2);
    a_cross_b(2) =   a(0) * b(1) - b(0) * a(1);

    return a_cross_b;
}

bool read_int_dataset_into_id(const hid_t & h5drm_dataset, std::string dataset_name, ID & result)
{
    hid_t ih5_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t ih5_ds = H5Dget_space(ih5_dataset);
    int ndims = H5Sget_simple_extent_ndims( ih5_ds);
    if (ndims != 1)
    {
        opserr << "H5DRM !! failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "H5DRM !! read_int_dataset_into_id - array dimension should be 1\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim;
    hsize_t maxdim;
    H5Sget_simple_extent_dims(ih5_ds, &dim, &maxdim);
    H5Sselect_all( ih5_ds );
    hid_t ih5_ms  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t ih5_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    int *d = new int[dim];

    H5Dread( ih5_dataset, H5T_NATIVE_INT, ih5_ms, ih5_ds, ih5_xfer_plist,  d);
    result.resize(dim);
    for (hsize_t i = 0; i < dim; ++i)
    {
        result(i) = d[i];
    }

    H5Sclose(ih5_ds);
    H5Sclose(ih5_ms);
    H5Dclose(ih5_dataset);

    delete [] d;

    return true;
}

bool read_double_dataset_into_vector(const hid_t & h5drm_dataset, std::string dataset_name, Vector & result)
{
    hid_t ih5_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t ih5_ds = H5Dget_space(ih5_dataset);
    int ndims = H5Sget_simple_extent_ndims( ih5_ds);
    if (ndims != 1)
    {
        opserr << "H5DRM !! failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "H5DRM !! read_double_dataset_into_vector - array dimension should be 1.\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim;
    hsize_t maxdim;
    H5Sget_simple_extent_dims(ih5_ds, &dim, &maxdim);
    H5Sselect_all( ih5_ds );
    hid_t ih5_ms  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t ih5_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    double *d = new double[dim];

    H5Dread( ih5_dataset, H5T_NATIVE_DOUBLE, ih5_ms, ih5_ds, ih5_xfer_plist,  d);
    result.resize(dim);
    for (hsize_t i = 0; i < dim; ++i)
    {
        result(i) = d[i];
    }

    H5Sclose(ih5_ds);
    H5Sclose(ih5_ms);
    H5Dclose(ih5_dataset);
    delete [] d;
    return true;
}


bool read_scalar_double_dataset_into_double(const hid_t & h5drm_dataset, std::string dataset_name, double & result)
{
    hid_t ih5_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t ih5_ds = H5Dget_space(ih5_dataset);
    int ndims = H5Sget_simple_extent_ndims( ih5_ds);
    if (ndims != 0)
    {
        opserr << "H5DRM !! failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "H5DRM !! read_scalar_double_dataset_into_double - array dimension should be 0\n";
        Vector error(-1);
        return false;
    }
    H5Sselect_all( ih5_ds );
    hsize_t dim = 1;
    hid_t ih5_ms  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t ih5_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    H5Dread( ih5_dataset, H5T_NATIVE_DOUBLE, ih5_ms, ih5_ds, ih5_xfer_plist,  &result);

    H5Sclose(ih5_ds);
    H5Sclose(ih5_ms);
    H5Dclose(ih5_dataset);

    return true;
}



bool read_double_dataset_into_matrix(const hid_t & h5drm_dataset, std::string dataset_name, Matrix & result)
{
    hid_t ih5_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t ih5_ds = H5Dget_space(ih5_dataset);
    int ndims = H5Sget_simple_extent_ndims( ih5_ds);
    if (ndims != 2)
    {
        opserr << "H5DRM !! failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "H5DRM !! read_double_dataset_into_matrix - array dimension should be 2\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim[2];
    hsize_t maxdim[2];
    H5Sget_simple_extent_dims(ih5_ds, dim, maxdim);
    H5Sselect_all( ih5_ds );
    hid_t ih5_ms  = H5Screate_simple(2, dim, 0);       // create dataspace of memory
    hid_t ih5_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    double *d = new double [dim[0]*dim[1]];

    if (d == 0)
    {
        opserr << "Fatal! no memory! 2163" << endln;
        exit(-1);
    }

    H5Dread( ih5_dataset, H5T_NATIVE_DOUBLE, ih5_ms, ih5_ds, ih5_xfer_plist,  d);
    result.resize(dim[0], dim[1]);

    for (int i = 0; i < (int) dim[0]; ++i)
    {
        for (int j = 0; j < (int) dim[1]; ++j)
        {
            result(i, j) = d[i * dim[1] + j];
        }
    }

    H5Sclose(ih5_ds);
    H5Sclose(ih5_ms);
    H5Dclose(ih5_dataset);

    delete [] d;
    return true;
}


bool read_int_dataset_into_array(const hid_t & h5drm_dataset, std::string dataset_name, int *& result)
{
    hid_t ih5_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t ih5_ds = H5Dget_space(ih5_dataset);
    int ndims = H5Sget_simple_extent_ndims( ih5_ds);
    if (ndims != 2)
    {
        opserr << "H5DRM !! failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "H5DRM !! read_int_dataset_into_array - array dimension should be 2\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim[2];
    hsize_t maxdim[2];
    H5Sget_simple_extent_dims(ih5_ds, dim, maxdim);
    H5Sselect_all( ih5_ds );
    hid_t ih5_ms  = H5Screate_simple(2, dim, 0);       // create dataspace of memory
    hid_t ih5_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    result = new int [dim[0]*dim[1]];

    H5Dread( ih5_dataset, H5T_NATIVE_INT, ih5_ms, ih5_ds, ih5_xfer_plist,  result);

    H5Sclose(ih5_ds);
    H5Sclose(ih5_ms);
    H5Dclose(ih5_dataset);

    return true;
}



#endif // _H5DRM
