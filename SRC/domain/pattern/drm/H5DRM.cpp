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
#include <iostream>
#include <string>
#include <math.h>
#include <time.h>
#include <hdf5.h>
#include <H5DRM.h>

#include <Channel.h>
#include <elementAPI.h>

#include <OPS_Globals.h>
#include <map>
#include <vector>
#include <Message.h>

#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#endif


#include <limits>   //For std::numeric_limits<double>::epsilon()

#include <Node.h>
#include <NodeIter.h>

#include <Element.h>
#include <ElementIter.h>


#define numNodeDOF 3  // Only nodes with 3-dofs per node canbe used in DRM... :/

using namespace std;

#define DEBUG_NODE_MATCHING false
#define DEBUG_DRM_INTEGRATION false
#define DEBUG_DRM_FORCES false


Vector H5DRM_calculate_cross_product(const Vector& a, const Vector& b);
bool read_int_dataset_into_id(const hid_t& h5drm_dataset, std::string dataset_name, ID& result);
bool read_scalar_double_dataset_into_double(const hid_t& h5drm_dataset, std::string dataset_name, double& result);
bool read_double_dataset_into_vector(const hid_t& h5drm_dataset, std::string dataset_name, Vector& result);
bool read_double_dataset_into_matrix(const hid_t& h5drm_dataset, std::string dataset_name, Matrix& result);
bool read_int_dataset_into_array(const hid_t& h5drm_dataset, std::string dataset_name, int *& result);
inline void convert_h5drmcrd_to_ops_crd(Vector&v );
inline void convert_h5drmcrd_to_ops_crd(Matrix&xyz );


static int numH5DRMpatterns = 0;

void *
OPS_ADD_RUNTIME_VPV(OPS_H5DRM)
{

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

    if (OPS_GetNumRemainingInputArgs() < 3)
    {
        OPS_GetDoubleInput(&num, &crd_scale);
        opserr << "crd_scale = " << crd_scale << endln;
    }

    thePattern = new H5DRM(tag, filename, factor, crd_scale);

    return thePattern;
}


















H5DRM::H5DRM()
    : LoadPattern(0, PATTERN_TAG_H5DRM),
      DRMForces(100),
      DRMDisplacements(100),
      DRMAccelerations(100),
      last_integration_time(0),
      crd_scale(0),
      distance_tolerance(0),
      maxnodetag(0),
      station_id2data_pos(100)
{
    is_initialized = false;
    t1 =  t2 =  tend = 0;
    cFactor = 1;
    step = step1 = step2 = 0;

    id_velocity = id_displacement = id_acceleration = 0;
    id_velocity_dataspace = id_displacement_dataspace = id_acceleration_dataspace = 0;


    myrank = 0;
#ifdef _PARALLEL_PROCESSING
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif
    H5DRMout << "H5DRM - empty constructor\n";
}

H5DRM::H5DRM(
    int tag,
    std::string HDF5filename_,
    double cFactor_,
    double crd_scale_,
    double distance_tolerance_)
    : LoadPattern(tag, PATTERN_TAG_H5DRM),
      HDF5filename(HDF5filename_),
      DRMForces(100),
      DRMDisplacements(100),
      DRMAccelerations(100),
      last_integration_time(0),
      crd_scale(crd_scale_),
      distance_tolerance(distance_tolerance_),
      maxnodetag(0),
      station_id2data_pos(100)
{

    id_velocity = id_displacement = id_acceleration = 0;
    id_velocity_dataspace = id_displacement_dataspace = id_acceleration_dataspace = 0;

    is_initialized = false;
    t1 =  t2 =  tend = 0;
    cFactor = cFactor_;
    step = step1 = step2 = 0;
    myrank = 0;
#ifdef _PARALLEL_PROCESSING
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#endif

    if (numH5DRMpatterns == 0)
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

void H5DRM::intitialize()
{
    Domain *theDomain = this->getDomain();

    if (theDomain == 0)
    {
        return;
    }

    if (myrank == 0)
        H5DRMout << "intializing - filename : " << HDF5filename << "\n";

    //===========================================================================
    // Open the specified file for read only
    //===========================================================================

    hid_t file_access_plist   = H5Pcreate(H5P_FILE_ACCESS);


    unsigned flags = H5F_ACC_RDONLY;
    id_drm_file = H5Fopen(HDF5filename.c_str(), flags, file_access_plist);

    if (id_drm_file < 0)
    {
        H5DRMerror << "Unable to open file: " << HDF5filename << endl;
        return;
    }


    //===========================================================================
    // MAP DRM data to local domain
    //===========================================================================

    for (int i = 0; i < 10; ++i)
    {
        // Plane * newplane = new Plane(id_drm_file, i);
        // planes.push_back(newplane);
    }


    H5DRMout << "Initialize" << endln;

    ID internal;
    Matrix xyz;
    Vector drmbox_x0;
    read_double_dataset_into_vector(id_drm_file, "DRM_Metadata/drmbox_x0", drmbox_x0);
    read_double_dataset_into_matrix(id_drm_file, "DRM_Data/xyz", xyz);
    read_int_dataset_into_id(id_drm_file, "DRM_Data/internal", internal);
    read_int_dataset_into_id(id_drm_file, "DRM_Data/data_location", station_id2data_pos);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_xmax", drmbox_xmax);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_xmin", drmbox_xmin);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_ymax", drmbox_ymax);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_ymin", drmbox_ymin);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_zmax", drmbox_zmax);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/drmbox_zmin", drmbox_zmin);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/dt", dt);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/tstart", tstart);
    read_scalar_double_dataset_into_double(id_drm_file, "DRM_Metadata/tend", tend);

    number_of_timesteps = int((tend - tstart)/dt);

    // last_integration_time = tstart;
    last_integration_time = theDomain->getCurrentTime();

    xyz *= crd_scale;
    drmbox_x0 *= crd_scale;
    drmbox_xmax *= crd_scale;
    drmbox_xmin *= crd_scale;
    drmbox_ymax *= crd_scale;
    drmbox_ymin *= crd_scale;
    drmbox_zmax *= crd_scale;
    drmbox_zmin *= crd_scale;

    convert_h5drmcrd_to_ops_crd(drmbox_x0);
    convert_h5drmcrd_to_ops_crd(xyz);

    int NDRM_points = xyz.noRows();
    double d_tol = distance_tolerance;
    double d_err = 0;
    int n_nodes_found = 0;

    node_matching_BruteForce(d_tol, internal, xyz, drmbox_x0, d_err, n_nodes_found );


    if (myrank == 0)
    {
        H5DRMout << "Dataset has " << NDRM_points << " data-points\n";
        opserr << "drmbox_x0   =  " << drmbox_x0 << " \n";
        opserr << "drmbox_xmax =  " << drmbox_xmax << " \n";
        opserr << "drmbox_xmin =  " << drmbox_xmin << " \n";
        opserr << "drmbox_ymax =  " << drmbox_ymax << " \n";
        opserr << "drmbox_ymin =  " << drmbox_ymin << " \n";
        opserr << "drmbox_zmax =  " << drmbox_zmax << " \n";
        opserr << "drmbox_zmin =  " << drmbox_zmin << " \n";
        opserr << "dx =  " << drmbox_xmax - drmbox_xmin << " \n";
        opserr << "dy =  " << drmbox_ymax - drmbox_ymin << " \n";
        opserr << "dz =  " << drmbox_zmax - drmbox_zmin << " \n";
    }

    if (myrank == 0)
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
    sprintf(viewdrm, "drm.%d.msh", myrank);
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
        fprintf(fptrdrm, "%d 15 2 1 %d %d\n", drmtag, 5 + myrank, drmtag);
        ++drmtag;
    }
    fprintf(fptrdrm, "$EndElements\n");


    fclose(fptrdrm);

    H5DRMout << "Found and connected " << n_nodes_found << " nodes (d_err = " << d_err << ")\n";
    DRMForces.resize(3 * n_nodes_found);
    DRMDisplacements.resize(3 * n_nodes_found);
    DRMAccelerations.resize(3 * n_nodes_found);

    DRMForces.Zero();
    DRMDisplacements.Zero();
    DRMAccelerations.Zero();

    char mshfilename[100];
    sprintf(mshfilename, "drmnodes.%d.msh", myrank);
    FILE * fptr = fopen(mshfilename, "w");
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
    fprintf(fptr, "%d\n", Nodes.Size());
    for (int i = 0; i < Nodes.Size(); ++i)
    {
        fprintf(fptr, "%d %d\n", Nodes(i), 1 * IsBoundary(i) + 2 * (1 - IsBoundary(i)));
    }
    fprintf(fptr, "$EndNodeData\n");
    fclose(fptr);







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
            if (Nodes.getLocation(node_tag) >= 0)
                count++;
        }
        if (count == nnodes)
            Elements.insert(element_tag);
    }

    number_of_local_elements = Elements.Size();
    H5DRMout << "Found  " << Elements.Size() << " elements on DRM boundary. \n";


    // char mshfilename[100];
    sprintf(mshfilename, "drmelements.%d.msh", myrank);
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
    fprintf(fptr, "%d\n", Elements.Size());
    for (int i = 0; i < Elements.Size(); ++i)
    {
        fprintf(fptr, "%d %d\n", Elements(i), 1);
    }
    fprintf(fptr, "$EndElementData\n");
    fclose(fptr);

    sprintf(mshfilename, "drmforces.%d.msh", myrank);
    fptr = fopen(mshfilename, "w");
    fprintf(fptr, "$MeshFormat\n");
    fprintf(fptr, "2.2 0 8\n");
    fprintf(fptr, "$EndMeshFormat\n");
    fclose(fptr);

    sprintf(mshfilename, "drmdisplacements.%d.msh", myrank);
    fptr = fopen(mshfilename, "w");
    fprintf(fptr, "$MeshFormat\n");
    fprintf(fptr, "2.2 0 8\n");
    fprintf(fptr, "$EndMeshFormat\n");
    fclose(fptr);

    sprintf(mshfilename, "drmaccelerations.%d.msh", myrank);
    fptr = fopen(mshfilename, "w");
    fprintf(fptr, "$MeshFormat\n");
    fprintf(fptr, "2.2 0 8\n");
    fprintf(fptr, "$EndMeshFormat\n");
    fclose(fptr);



    //===========================================================================
    // Open Displacements and Accelerations arrays for reading
    //===========================================================================
    htri_t found_it = 0;


    found_it = H5Lexists(id_drm_file, "DRM_Data/velocity", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening velocity dataset.\n";


        id_velocity = H5Dopen(id_drm_file, "DRM_Data/velocity", H5P_DEFAULT);

        if (id_velocity < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/velocity\" array!\n ";
            return;
        }

        id_velocity_dataspace = H5Dget_space(id_velocity);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        id_one_node_memspace  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    found_it = H5Lexists(id_drm_file, "DRM_Data/acceleration", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening acceleration dataset.\n";


        id_acceleration = H5Dopen(id_drm_file, "DRM_Data/acceleration", H5P_DEFAULT);

        if (id_acceleration < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/acceleration\" array!\n ";
            return;
        }

        id_acceleration_dataspace = H5Dget_space(id_acceleration);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        id_one_node_memspace  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    found_it = H5Lexists(id_drm_file, "DRM_Data/displacement", H5P_DEFAULT );
    if (found_it)
    {
        H5DRMout << "Opening displacement dataset.\n";


        id_displacement = H5Dopen(id_drm_file, "DRM_Data/displacement", H5P_DEFAULT);

        if (id_displacement < 0)
        {
            H5DRMerror << "Could not open \"DRM_Data/displacement\" array!\n ";
            return;
        }

        id_displacement_dataspace = H5Dget_space(id_displacement);

        hsize_t rank_two_array = 2;
        hsize_t one_node_data_dims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        hsize_t one_node_data_maxdims[2] = {3, H5DRM_NUM_OF_PRECOMPUTED_TIMESTEPS};
        id_one_node_memspace  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory
    }

    id_xfer_plist = H5Pcreate( H5P_DATASET_XFER);


//===========================================================================
// Set status to initialized and ready to compute loads
//===========================================================================
    if (myrank == 0)
    {
        H5DRMout << "dt =   " << dt << " \n";
        H5DRMout << "number_of_timesteps =   " << number_of_timesteps << " \n";
        H5DRMout << "tstart =   " << tstart << " \n";
        H5DRMout << "tend =   " << tend << " \n";
        H5DRMout << "DRM initialization done!\n";
    }

    is_initialized = true;
}



H5DRM::~H5DRM()
{
    clean_all_data();
    // for (int i = 0; i < 10; ++i)
    // for (auto iter = planes.begin(); iter != planes.end(); ++iter)
    // {
    //     delete *iter;// planes[i];
    // }
    // planes.clear();
}

void H5DRM::clean_all_data()
{

    nodetag2station_id.clear();
    nodetag2local_pos.clear();

    if (id_velocity > 0 && H5Dclose(id_velocity) < 0)
    {
        H5DRMerror << "Unable to close motion dataset. " << endl;
    }

    if (id_velocity_dataspace > 0 && H5Sclose(id_velocity_dataspace) < 0)
    {
        H5DRMerror << "Unable to close motion dataspace. " << endl;
    }

    if (id_one_node_memspace > 0 && H5Sclose(id_one_node_memspace) < 0)
    {
        H5DRMerror << "Unable to close one node memspace. " << endl;
    }

    hid_t obj_id_list[H5DRM_MAX_RETURN_OPEN_OBJS];
    hsize_t n_obj_open = 10;
    while (id_drm_file > 0 and n_obj_open > 0)
    {
        int n_objects_closed = 0;
        n_obj_open = H5Fget_obj_count(id_drm_file, H5F_OBJ_DATASET | H5F_OBJ_GROUP | H5F_OBJ_ATTR | H5F_OBJ_LOCAL );
        cout << "H5DRM -- N of HDF5 objects open = " << n_obj_open << endl;

        if (n_obj_open <= 0)
        {
            break;
        }

        int number_of_open_objects;

        //Close datasets
        number_of_open_objects = H5Fget_obj_ids( id_drm_file, H5F_OBJ_DATASET | H5F_OBJ_LOCAL, H5DRM_MAX_RETURN_OPEN_OBJS, obj_id_list );
        if (number_of_open_objects > 0)
        {
            cout << "                                   - Closing " <<  number_of_open_objects << " datasets.\n";
        }
        for (int i = 0; i < std::min(number_of_open_objects, H5DRM_MAX_RETURN_OPEN_OBJS) ; i++ )
        {
            H5Dclose(obj_id_list[i]);
            n_objects_closed++;
        }

        //Close groups
        number_of_open_objects = H5Fget_obj_ids( id_drm_file, H5F_OBJ_GROUP | H5F_OBJ_LOCAL, H5DRM_MAX_RETURN_OPEN_OBJS, obj_id_list );
        if (number_of_open_objects > 0)
        {
            cout << "                                   - Closing " <<  number_of_open_objects << " groups.\n";
        }
        for (int i = 0; i < std::min(number_of_open_objects, H5DRM_MAX_RETURN_OPEN_OBJS) ; i++ )
        {
            H5Gclose(obj_id_list[i]);
            n_objects_closed++;
        }

        //Close groups
        number_of_open_objects = H5Fget_obj_ids( id_drm_file, H5F_OBJ_ATTR | H5F_OBJ_LOCAL, H5DRM_MAX_RETURN_OPEN_OBJS, obj_id_list );
        if (number_of_open_objects > 0)
        {
            cout << "                                   - Closing " <<  number_of_open_objects << " attributes.\n";
        }
        for (int i = 0; i < std::min(number_of_open_objects, H5DRM_MAX_RETURN_OPEN_OBJS) ; i++ )
        {
            H5Aclose(obj_id_list[i]);
            n_objects_closed++;
        }

        if (n_objects_closed == 0)
        {
            //This guards against the possibility of H5Fget_obj_count misbehaving.
            // ot that the open objects are not datasets, groups or attributes.
            break;
        }

    }

    if (id_drm_file > 0 && H5Fclose(id_drm_file) < 0)
    {
        H5DRMerror << "Unable to close HDF5 file \"" << HDF5filename << "\"\n";
    }

    is_initialized = false;
}

void
H5DRM::setDomain(Domain *theDomain)
{
    this->LoadPattern::setDomain(theDomain);
}


void
H5DRM::applyLoad(double time)
{
    if (not is_initialized)
    {
        intitialize();
    }

    //If we don't have elements, no point in doing anything :/
    if (number_of_local_elements <= 0)
    {
        return;
    }
    else
    {
        Domain *theDomain = this->getDomain();

        bool computedLoads = ComputeDRMLoads(time);

        if (theDomain == 0 || Elements.Size() == 0 || !computedLoads)
        {
            opserr << "H5DRM::applyLoad -- Error! Failed to compute DRM loads at time = " << time << endln;
            return;
        }

        Node *theNode = 0;
        static Vector load(3);
        load.Zero();

        FILE* fptr_forces = 0;
        FILE* fptr_displ = 0;
        FILE* fptr_accel = 0;

        static int step = 0;
        char mshfilename[100];

        // drmforces
        sprintf(mshfilename, "drmforces.%d.msh", myrank);
        fptr_forces = fopen(mshfilename, "a");
        // drmdisplacements
        sprintf(mshfilename, "drmdisplacements.%d.msh", myrank);
        fptr_displ = fopen(mshfilename, "a");
        // drmaccelerations
        sprintf(mshfilename, "drmaccelerations.%d.msh", myrank);
        fptr_accel = fopen(mshfilename, "a");
        
        fprintf(fptr_forces, "$NodeData\n");
        fprintf(fptr_forces, "    1\n");
        fprintf(fptr_forces, "\"DRM Forces\"\n");
        fprintf(fptr_forces, "1\n");
        fprintf(fptr_forces, "%f\n", time);
        fprintf(fptr_forces, "3\n");
        fprintf(fptr_forces, "%d\n", step);
        fprintf(fptr_forces, "3\n");
        fprintf(fptr_forces, "%d\n", Nodes.Size());

        fprintf(fptr_displ, "$NodeData\n");
        fprintf(fptr_displ, "    1\n");
        fprintf(fptr_displ, "\"DRM Displacements\"\n");
        fprintf(fptr_displ, "1\n");
        fprintf(fptr_displ, "%f\n", time);
        fprintf(fptr_displ, "3\n");
        fprintf(fptr_displ, "%d\n", step);
        fprintf(fptr_displ, "3\n");
        fprintf(fptr_displ, "%d\n", Nodes.Size());

        fprintf(fptr_accel, "$NodeData\n");
        fprintf(fptr_accel, "    1\n");
        fprintf(fptr_accel, "\"DRM Acceleration\"\n");
        fprintf(fptr_accel, "1\n");
        fprintf(fptr_accel, "%f\n", time);
        fprintf(fptr_accel, "3\n");
        fprintf(fptr_accel, "%d\n", step);
        fprintf(fptr_accel, "3\n");
        fprintf(fptr_accel, "%d\n", Nodes.Size());


        ++step;

        for (int i = 0; i < Nodes.Size(); i++)
        {
            int nodeTag = Nodes[i] ;
            int local_pos = nodetag2local_pos[nodeTag];
            theNode = theDomain->getNode( nodeTag );
            if ( theNode == 0 )
                continue;

            if (IsBoundary[local_pos])
            {
                load(0) = -DRMForces(3 * local_pos + 0);
                load(1) = -DRMForces(3 * local_pos + 1);
                load(2) = -DRMForces(3 * local_pos + 2);
            }
            else
            {
                load(0) = DRMForces(3 * local_pos + 0);
                load(1) = DRMForces(3 * local_pos + 1);
                load(2) = DRMForces(3 * local_pos + 2);
            }

            fprintf(fptr_forces, "%d %f %f %f\n", nodeTag, load(0), load(1), load(2));
            fprintf(fptr_displ, "%d %f %f %f\n", nodeTag, DRMDisplacements(3 * local_pos + 0), DRMDisplacements(3 * local_pos + 1), DRMDisplacements(3 * local_pos + 2));
            fprintf(fptr_accel, "%d %f %f %f\n", nodeTag, DRMAccelerations(3 * local_pos + 0), DRMAccelerations(3 * local_pos + 1), DRMAccelerations(3 * local_pos + 2));

            //Add to current nodal unbalanced load
            theNode->addUnbalancedLoad(load);
        }
        fprintf(fptr_forces, "$EndNodeData\n");
        fclose(fptr_forces);
        fprintf(fptr_displ, "$EndNodeData\n");
        fclose(fptr_displ);
        fprintf(fptr_accel, "$EndNodeData\n");
        fclose(fptr_accel);
    }
}


bool H5DRM::ComputeDRMMotions(double next_integration_time)
{
    bool have_displacement = id_displacement > 0;
    // bool have_velocity = id_velocity > 0;
    // bool have_acceleration = id_acceleration > 0;

    if(have_displacement )//&& !have_acceleration)
        return drm_differentiate_displacements(next_integration_time);

    //JAA Disable these modes for now

    // if(have_displacement && have_acceleration) // Note: disabled!
    //     return drm_direct_read(next_integration_time);

    // if(!have_displacement && !have_acceleration && have_velocity) 
    //     return drm_integrate_velocity(next_integration_time);

    return false;
}

bool H5DRM::drm_direct_read(double t)
{

    if (Nodes.Size() == 0)
        return false;

    DRMDisplacements.Zero();
    DRMAccelerations.Zero();

    if (t < tstart || t > tend)
    {
        return true;
    }
    
    int i1 = (int) floor( (t - tstart) / dt);
    int i2 = (int) floor( (t - tstart) / dt) + 1;
    double t1 = i1*dt + tstart;
    double t2 = i2*dt + tstart;

    double dtau = (t - t1)/(t2-t1);

    id_acceleration_dataspace = H5Dget_space(id_acceleration);
    id_displacement_dataspace = H5Dget_space(id_displacement);

    H5DRMout << "t = " << t << " dt = " << dt << " i1 = " << i1 << " i2 = " << i2 << " t1 = " << t1 << " t2 = " << t2 << " dtau = " << dtau << endln;

    double umax = -std::numeric_limits<double>::infinity();
    double amax = -std::numeric_limits<double>::infinity();
    double umin =  std::numeric_limits<double>::infinity();
    double amin =  std::numeric_limits<double>::infinity();

    for (int n = 0; n < Nodes.Size(); ++n)
    {
        int nodeTag = Nodes(n);
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
            id_displacement_dataspace,
            H5S_SELECT_SET, start, stride, count, block );
        H5Sselect_hyperslab(
            id_acceleration_dataspace,
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
        herr_t errorflag1 = H5Dread( id_displacement, H5T_NATIVE_DOUBLE, memspace,
                                    id_displacement_dataspace, id_xfer_plist,  d1 );
        herr_t errorflag2 = H5Dread(id_acceleration, H5T_NATIVE_DOUBLE, memspace,
                                    id_acceleration_dataspace, id_xfer_plist,  a1 );
        start[1] = (hsize_t) i2;
        H5Sselect_hyperslab(
            id_displacement_dataspace,
            H5S_SELECT_SET, start, stride, count, block );
        H5Sselect_hyperslab(
            id_acceleration_dataspace,
            H5S_SELECT_SET, start, stride, count, block );


        herr_t errorflag3 = H5Dread( id_displacement, H5T_NATIVE_DOUBLE, memspace,
                                    id_displacement_dataspace, id_xfer_plist,  d2 );
        herr_t errorflag4 = H5Dread(id_acceleration, H5T_NATIVE_DOUBLE, memspace,
                                    id_acceleration_dataspace, id_xfer_plist,  a2 );
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
            H5DRMerror << "H5DRM::drm_direct_read - Failed to read displacement or acceleration array!!\n" <<
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

        
        DRMDisplacements(3 * local_pos + 0) = d1[0]*(1-dtau) + d2[0]*(dtau);
        DRMDisplacements(3 * local_pos + 1) = d1[1]*(1-dtau) + d2[1]*(dtau);
        DRMDisplacements(3 * local_pos + 2) = d1[2]*(1-dtau) + d2[2]*(dtau);

        DRMAccelerations(3 * local_pos + 0) = a1[0]*(1-dtau) + a2[0]*(dtau);
        DRMAccelerations(3 * local_pos + 1) = a1[1]*(1-dtau) + a2[1]*(dtau);
        DRMAccelerations(3 * local_pos + 2) = a1[2]*(1-dtau) + a2[2]*(dtau);
    }

    H5DRMout << "t = " << t << " u = (" << umin << ", " << umax << ") a = (" << amin << ", " << amax << ")" << endln;

    return true;
}


bool H5DRM::drm_differentiate_displacements(double t)
{

    if (Nodes.Size() == 0)
        return false;

    DRMDisplacements.Zero();
    DRMAccelerations.Zero();

    if (t < tstart || t > tend)
    {
        return true;
    }
    
    int i1 = (int) floor( (t - tstart) / dt);
    int i2 = (int) floor( (t - tstart) / dt) + 1;
    double t1 = i1*dt + tstart;
    double t2 = i2*dt + tstart;

    double dtau = (t - t1)/(t2-t1);

    id_displacement_dataspace = H5Dget_space(id_displacement);
    hsize_t i_first = i1 - 1;
    hsize_t i_last  = i1 + 2;

    i_first = i_first < 0 ? 0 : i_first;
    i_first = i_first > (hsize_t) number_of_timesteps-1 ? number_of_timesteps-1 : i_first;

    i_last = i_last < 0 ? 0 : i_last;
    i_last = i_last > (hsize_t) number_of_timesteps-1 ? number_of_timesteps-1 : i_last;

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
    double dt2 = dt*dt;



    for (int n = 0; n < Nodes.Size(); ++n)
    {
        int nodeTag = Nodes(n);
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
            id_displacement_dataspace,
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
            errorflag1 = H5Dread( id_displacement, H5T_NATIVE_DOUBLE, memspace,
                                    id_displacement_dataspace, id_xfer_plist,  d0 );
        }
        

        for (int dof = 0; dof < 3; ++dof)
        {
            d1[dof] = d0[dof][1];
            d2[dof] = d0[dof][2];
            a1[dof] = (d0[dof][0] -2*d0[dof][1] +d0[dof][2])/dt2;
            a2[dof] = (d0[dof][1] -2*d0[dof][2] +d0[dof][3])/dt2;
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

        
        DRMDisplacements(3 * local_pos + 0) = d1[0]*(1-dtau) + d2[0]*(dtau);
        DRMDisplacements(3 * local_pos + 1) = d1[1]*(1-dtau) + d2[1]*(dtau);
        DRMDisplacements(3 * local_pos + 2) = d1[2]*(1-dtau) + d2[2]*(dtau);

        DRMAccelerations(3 * local_pos + 0) = a1[0]*(1-dtau) + a2[0]*(dtau);
        DRMAccelerations(3 * local_pos + 1) = a1[1]*(1-dtau) + a2[1]*(dtau);
        DRMAccelerations(3 * local_pos + 2) = a1[2]*(1-dtau) + a2[2]*(dtau);
    }

    H5DRMout << "t = " << t << " u = (" << umin << ", " << umax << ") a = (" << amin << ", " << amax << ")" << endln;

    return true;
}

bool H5DRM::drm_integrate_velocity(double next_integration_time)
{

    exit(-1);
    return false;

//     FILE * fptr = 0;
//     if (DEBUG_DRM_INTEGRATION)
//     {
//         char debugfilename[100];
//         sprintf(debugfilename, "drmintegration.%d.txt", myrank);
//         fptr = fopen(debugfilename, "w");
//     }


//     if (Nodes.Size() == 0)
//         return false;

//     if (next_integration_time < tstart || next_integration_time > tend)
//     {
//         DRMDisplacements.Zero();
//         DRMAccelerations.Zero();
//         return false;
//     }

//     int dir = 1;
//     double t1 = last_integration_time;
//     double t2 = next_integration_time;
//     if (t1 > t2)
//     {
//         // Integrating backwards
//         t1 = next_integration_time;
//         t2 = last_integration_time;
//         dir = -1;
//     }
//     hsize_t i1 = (hsize_t) (t1 - tstart) / dt;
//     hsize_t i2 = (hsize_t) (t2 - tstart) / dt + 1;
//     hsize_t Nt = i2 - i1;

//     hsize_t motions_dims[2];
//     id_velocity_dataspace = H5Dget_space(id_velocity);

//     hsize_t imax_t = motions_dims[1] - 1;


//     i1 = i1 < 0 ? 0 : i1 ;
//     i1 = i1 > imax_t ? imax_t : i1;
//     i2 = i2 < 0 ? 0 : i2 ;
//     i2 = i2 > imax_t ? imax_t : i2;
//     Nt = i1 == i2 ? 1 : Nt;

//     for (int n = 0; n < Nodes.Size(); ++n)
//     {
//         int nodeTag = Nodes(n);
//         int station_id = nodetag2station_id[nodeTag];
//         int data_pos = station_id2data_pos[station_id];
//         int local_pos = nodetag2local_pos[nodeTag];

//         double v[3][Nt];

//         hsize_t start[2]  = {(hsize_t) data_pos , (hsize_t)i1};
//         hsize_t stride[2] = {1                  , 1};
//         hsize_t count[2]  = {1                  , 1};
//         hsize_t block[2]  = {3                  , (hsize_t) Nt};

//         //Selection in dataspace
//         H5Sselect_hyperslab(
//             id_velocity_dataspace,
//             H5S_SELECT_SET, start, stride, count, block );

//         //Selection in memspace
//         hsize_t rank_two_array = 2;
//         hsize_t one_node_data_dims[2] = {3, Nt};
//         hsize_t one_node_data_maxdims[2] = {3, Nt};
//         hid_t memspace  = H5Screate_simple(rank_two_array, one_node_data_dims, one_node_data_maxdims);       // create dataspace of memory


//         hsize_t mem_start[2] = {0 , 0};
//         H5Sselect_hyperslab(
//             memspace,
//             H5S_SELECT_SET, mem_start, stride, count, block );

//         //Read data
//         herr_t errorflag = H5Dread( id_velocity, H5T_NATIVE_DOUBLE, memspace,
//                                     id_velocity_dataspace, id_xfer_plist,  v );

//         H5Sclose(memspace);

//         if (errorflag < 0)
//         {
//             H5DRMerror << "Failed to read velocity array!!\n" <<
//                        " n = " << n << endln <<
//                        " nodeTag = " << nodeTag << endln <<
//                        " station_id = " << station_id << endln <<
//                        " i1 = " << i1 << endln <<
//                        " data_pos = " << data_pos << endln <<
//                        " local_pos = " << local_pos << endln <<
//                        " Nt = " << Nt << endln <<
//                        " last_integration_time = " << last_integration_time << endln <<
//                        " start   = [" << start[0] << ", " << start[1] << "]" << endln <<
//                        " mem_start  = [" << mem_start[0] << ", " << mem_start[1] << "]" << endln <<
//                        " stride  = [" << stride[0] << ", " << stride[1] << "]" << endln <<
//                        " count   = [" << count[0] << ", " << count[1] << "]" << endln <<
//                        " block   = [" << block[0] << ", " << block[1] << "]" << endln;
//             exit(-1);
//         }

//         for (hsize_t i = 0; i <  Nt; ++i)
//         {
//             double v0 = v[0][i];
//             double v1 = v[1][i];
//             v[0][i] = v0;
//             v[1][i] = v1;
//             v[2][i] = -v[2][i];
//         }


// <<<<<<< HEAD
// 		for (hsize_t i = 0; i < Nt; ++i)
// 		{
// 			double dtau = 0;
// 			double tau_1 = tstart + i * dt;
// 			double tau_2 = tstart + (i + 1) * dt;
// 			tau_1 = tau_1 > t1 ? tau_1 : t1;
// 			tau_2 = tau_2 < t2 ? tau_2 : t2;
// 			dtau = tau_2 - tau_1;

// 			if (dtau <= 0)
// 				continue;



// 			if (DEBUG_DRM_INTEGRATION)
// 			{
// 				/* FMK
// 				fprintf(fptr, "i = %d dtau=%f tau_1=%f tau_2=%f dt=%f local_pos=%d dir=%d\n", (int) i, dtau, tau_1, tau_2, dt, local_pos, dir );
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0][i], v[0][i + 1] );
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1][i], v[1][i + 1] );
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2][i], v[2][i + 1] );
// 				******/
// 				fprintf(fptr, "i = %d dtau=%f tau_1=%f tau_2=%f dt=%f local_pos=%d dir=%d\n", (int)i, dtau, tau_1, tau_2, dt, local_pos, dir);
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0 + i * Nt], v[0 + (i + 1)*3]);
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1 + i * Nt], v[1 + (i + 1)*3]);
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2 + i * Nt], v[2 + (i + 1)*3]);

// 			}

// 			double u1 = DRMDisplacements(3 * local_pos + 0);
// 			double u2 = DRMDisplacements(3 * local_pos + 1);
// 			double u3 = DRMDisplacements(3 * local_pos + 2);
// 			/* FMK
// 			double du1 = (v[0][i] + v[0][i + 1]) * (dir * dtau / 2);
// 			double du2 = (v[1][i] + v[1][i + 1]) * (dir * dtau / 2);
// 			double du3 = (v[2][i] + v[2][i + 1]) * (dir * dtau / 2);
// 			******/
// 			double du1 = (v[0 + i * 3] + v[0 + (i + 1)*3]) * (dir * dtau / 2);
// 			double du2 = (v[1 + i * 3] + v[1 + (i + 1)*3]) * (dir * dtau / 2);
// 			double du3 = (v[2 + i * 3] + v[2 + (i + 1)*3]) * (dir * dtau / 2);

// 			DRMDisplacements(3 * local_pos + 0) += du1;
// 			DRMDisplacements(3 * local_pos + 1) += du2;
// 			DRMDisplacements(3 * local_pos + 2) += du3;

// 			if (DEBUG_DRM_INTEGRATION)
// 			{
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 0) = %f \n", DRMDisplacements(3 * local_pos + 0));
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 1) = %f \n", DRMDisplacements(3 * local_pos + 1));
// 				fprintf(fptr, "    DRMDisplacements(3 * local_pos + 2) = %f \n", DRMDisplacements(3 * local_pos + 2));
// 			}

// 			// bool found_nan = false;
// 			if (isnan(u1) || isnan(du1) ||
// 				isnan(u2) || isnan(du2) ||
// 				isnan(u3) || isnan(du3) ||
// 				isnan(dt) || isnan(dtau))
// 			{
// 				H5DRMerror << "NAN Detected!!! \n";
// 				H5DRMerror << "    nodeTag = " << nodeTag << endln;
// 				H5DRMerror << "    local_pos = " << local_pos << endln;
// 				printf("    i = %d dtau=%f tau_1=%f tau_2=%f dt=%f local_pos=%d dir=%d\n", (int)i, dtau, tau_1, tau_2, dt, local_pos, dir);
// 				printf("        u1 = %f du1 = %f \n", u1, du1);
// 				printf("        u2 = %f du2 = %f \n", u2, du2);
// 				printf("        u3 = %f du3 = %f \n", u3, du3);
// 				/* FMK
// 				printf("        DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0][i], v[0][i + 1] );
// 				printf("        DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1][i], v[1][i + 1] );
// 				printf("        DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2][i], v[2][i + 1] );
// 				*/
// 				printf("        DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0 + i * 3], v[0 + (i + 1)*3]);
// 				printf("        DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1 + i * 3], v[1 + (i + 1)*3]);
// 				printf("        DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2 + i * 3], v[2 + (i + 1)*3]);
// 				exit(-1);
// 			}

// 		}
// =======
//         for (hsize_t i = 0; i <  Nt; ++i)
//         {
//             double dtau = 0;
//             double tau_1 = tstart + i * dt;
//             double tau_2 = tstart + (i + 1) * dt;
//             tau_1 = tau_1 > t1 ? tau_1 : t1;
//             tau_2 = tau_2 < t2 ? tau_2 : t2;
//             dtau = tau_2 - tau_1;

//             if (dtau <= 0)
//                 continue;



//             if (DEBUG_DRM_INTEGRATION)
//             {
//                 fprintf(fptr, "i = %d dtau=%f tau_1=%f tau_2=%f dt=%f local_pos=%d dir=%d\n", (int) i, dtau, tau_1, tau_2, dt, local_pos, dir );
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0][i], v[0][i + 1] );
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1][i], v[1][i + 1] );
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2][i], v[2][i + 1] );
//             }

//             double u1 = DRMDisplacements(3 * local_pos + 0);
//             double u2 = DRMDisplacements(3 * local_pos + 1);
//             double u3 = DRMDisplacements(3 * local_pos + 2);
//             double du1 = (v[0][i] + v[0][i + 1]) * (dir * dtau / 2);
//             double du2 = (v[1][i] + v[1][i + 1]) * (dir * dtau / 2);
//             double du3 = (v[2][i] + v[2][i + 1]) * (dir * dtau / 2);

//             DRMDisplacements(3 * local_pos + 0) += du1;
//             DRMDisplacements(3 * local_pos + 1) += du2;
//             DRMDisplacements(3 * local_pos + 2) += du3;

//             if (DEBUG_DRM_INTEGRATION)
//             {
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 0) = %f \n", DRMDisplacements(3 * local_pos + 0) );
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 1) = %f \n", DRMDisplacements(3 * local_pos + 1) );
//                 fprintf(fptr, "    DRMDisplacements(3 * local_pos + 2) = %f \n", DRMDisplacements(3 * local_pos + 2) );
//             }

//             // bool found_nan = false;
//             if (isnan(u1) || isnan(du1) ||
//                     isnan(u2) || isnan(du2) ||
//                     isnan(u3) || isnan(du3) ||
//                     isnan(dt) || isnan(dtau) )
//             {
//                 H5DRMerror << "NAN Detected!!! \n";
//                 H5DRMerror << "    nodeTag = " << nodeTag << endln;
//                 H5DRMerror << "    local_pos = " << local_pos << endln;
//                 printf("    i = %d dtau=%f tau_1=%f tau_2=%f dt=%f local_pos=%d dir=%d\n", (int)i,  dtau, tau_1, tau_2, dt, local_pos, dir );
//                 printf("        u1 = %f du1 = %f \n", u1, du1);
//                 printf("        u2 = %f du2 = %f \n", u2, du2);
//                 printf("        u3 = %f du3 = %f \n", u3, du3);
//                 printf("        DRMDisplacements(3 * local_pos + 0) = %f -->  v[0][i] = %f  v[0][i+1] = %f\n", DRMDisplacements(3 * local_pos + 0), v[0][i], v[0][i + 1] );
//                 printf("        DRMDisplacements(3 * local_pos + 1) = %f -->  v[1][i] = %f  v[1][i+1] = %f\n", DRMDisplacements(3 * local_pos + 1), v[1][i], v[1][i + 1] );
//                 printf("        DRMDisplacements(3 * local_pos + 2) = %f -->  v[2][i] = %f  v[2][i+1] = %f\n", DRMDisplacements(3 * local_pos + 2), v[2][i], v[2][i + 1] );
//                 exit(-1);
//             }

//         }
// >>>>>>> 732044f5b1846fa0047e4516e0717dc0d22bffc3

//         int eval_i = (int) (t2 - t1) / dt;


//         DRMAccelerations(3 * local_pos + 0) = (v[0][eval_i] - v[0][eval_i + 1]) / dt;
//         DRMAccelerations(3 * local_pos + 1) = (v[1][eval_i] - v[1][eval_i + 1]) / dt;
//         DRMAccelerations(3 * local_pos + 2) = (v[2][eval_i] - v[2][eval_i + 1]) / dt;

//         if (DEBUG_DRM_INTEGRATION)
//         {
//             fprintf(fptr, "eval_i = %d \n", eval_i );
//             fprintf(fptr, "     v[0][eval_i] = %f  v[0][eval_i+1] = %f\n",  v[0][eval_i], v[0][eval_i + 1] );
//             fprintf(fptr, "     v[1][eval_i] = %f  v[1][eval_i+1] = %f\n",  v[1][eval_i], v[1][eval_i + 1] );
//             fprintf(fptr, "     v[2][eval_i] = %f  v[2][eval_i+1] = %f\n",  v[2][eval_i], v[2][eval_i + 1] );
//             fprintf(fptr, "     DRMAccelerations(3 * local_pos + 0) = %f\n",  DRMAccelerations(3 * local_pos + 0) );
//             fprintf(fptr, "     DRMAccelerations(3 * local_pos + 1) = %f\n",  DRMAccelerations(3 * local_pos + 1) );
//             fprintf(fptr, "     DRMAccelerations(3 * local_pos + 2) = %f\n",  DRMAccelerations(3 * local_pos + 2) );

//         }
//     }



//     last_integration_time = t2;

//     return true;
}


bool
H5DRM::ComputeDRMLoads(double t)
{
    if (myrank == 0)
        H5DRMout << "ComputeDRMLoads.... Begin (t = " <<  t << ").\n";


    int NDOF = 3;
    int NIEMAX = 8;
    static Vector Fm(NIEMAX * NDOF); Fm.Zero();
    static Vector Fk(NIEMAX * NDOF); Fk.Zero();
    static Vector u_e(NIEMAX * NDOF); u_e.Zero();
    static Vector udd_e(NIEMAX * NDOF); udd_e.Zero();
    static Vector v1(NDOF); v1.Zero();
    static Vector v2(NDOF); v2.Zero();
    static ID B_node(NIEMAX); B_node.Zero();
    static ID E_node(NIEMAX); E_node.Zero();

    DRMForces.Zero();
    
    if ( t < tstart || tend < t)
    {
        // Will assume that DRM forces are null outside of the defined range...
        return true;
    }


    if (number_of_local_elements > 0)
    {
        bool computed = ComputeDRMMotions(t);

        if (!computed)
            return false;

        Domain *theDomain = this->getDomain();

        Element *theElement = theDomain->getElement( Elements[0] );
        int NIE = 8;
        B_node.resize(NIEMAX);
        E_node.resize(NIEMAX);

        for (int e = 0; e < Elements.Size(); e++)
        {
            int eleTag = Elements[e];

            //Pointer to current element
            theElement = theDomain->getElement(eleTag);

            //List of current element nodes
            const ID &elementNodes = theElement->getExternalNodes();

            //Identify boundary and exterior nodes for this element
            NIE = elementNodes.Size();
            
            int nB = 0, nE = 0;
            for ( int ii = 0; ii < NIE; ii++)
            {
                int nodeTag = elementNodes(ii);
                int local_pos = nodetag2local_pos[nodeTag];

                if ( IsBoundary[local_pos] == 1 )
                {
                    B_node(nB) = ii;
                    nB ++;
                }
                else
                {
                    E_node(nE) = ii;
                    nE ++;
                }
            }

            if ( nB != 0 and nE != 0 )
            {
                //Mass and stiffness matrices
                Matrix Me = theElement->getMass();
                Matrix Ke = theElement->getTangentStiff();

                //Zero out the diagonal block of Boundary nodes
                for (int m = 0; m < nB; m++)
                    for (int n = 0; n < nB; n++)
                        for (int d = 0; d < NDOF; d++)
                            for (int e = 0; e < NDOF; e++)
                            {
                                Me( B_node(m)*NDOF + d, B_node(n)*NDOF + e ) = 0.0;
                                Ke( B_node(m)*NDOF + d, B_node(n)*NDOF + e ) = 0.0;
                            }


                //Zero out the diagonal block of Exterior nodes
                for (int m = 0; m < nE; m++)
                    for (int n = 0; n < nE; n++)
                        for (int d = 0; d < NDOF; d++)
                            for (int e = 0; e < NDOF; e++)
                            {
                                Me( E_node(m)*NDOF + d, E_node(n)*NDOF + e ) = 0.0;
                                Ke( E_node(m)*NDOF + d, E_node(n)*NDOF + e ) = 0.0;
                            }

                u_e.resize(3 * NIE);
                udd_e.resize(3 * NIE);
                Fm.resize(3 * NIE);
                Fk.resize(3 * NIE);

                u_e.Zero();
                udd_e.Zero();
                Fm.Zero();
                Fk.Zero();

                for (int k = 0; k < NIE; ++k)
                {
                    int nodeTag = elementNodes(k);
                    int local_pos = nodetag2local_pos[nodeTag];
                    u_e(3 * k      ) = DRMDisplacements[3 * local_pos];// - DRMDisplacements0[3 * local_pos] ;
                    u_e(3 * k + 1  ) = DRMDisplacements[3 * local_pos + 1];// - DRMDisplacements0[3 * local_pos + 1] ;
                    u_e(3 * k + 2  ) = DRMDisplacements[3 * local_pos + 2];// - DRMDisplacements0[3 * local_pos + 2] ;
                    udd_e(3 * k    ) = DRMAccelerations[3 * local_pos];
                    udd_e(3 * k + 1) = DRMAccelerations[3 * local_pos + 1];
                    udd_e(3 * k + 2) = DRMAccelerations[3 * local_pos + 2];
                }
                // cout << endl;

                Fm.addMatrixVector(0.0, Me, udd_e, 1.0);
                Fk.addMatrixVector(0.0, Ke, u_e, 1.0);

                for (int k = 0; k < NIE; k++)
                {
                    int nodeTag = elementNodes(k);
                    int local_pos = nodetag2local_pos[nodeTag];
                    DRMForces( 3 * local_pos  + 0) +=  Fk(3 * k  + 0) + Fm(3 * k  + 0);
                    DRMForces( 3 * local_pos  + 1) +=  Fk(3 * k  + 1) + Fm(3 * k  + 1);
                    DRMForces( 3 * local_pos  + 2) +=  Fk(3 * k  + 2) + Fm(3 * k  + 2);


                    if (isnan(DRMForces( 3 * local_pos  + 0) ) ||
                            isnan(DRMForces( 3 * local_pos  + 1) ) ||
                            isnan(DRMForces( 3 * local_pos  + 2) ) )
                    {
                        H5DRMerror << "NAN Detected!!! \n";
                        H5DRMerror << "    nodeTag = " << nodeTag << endln;
                        H5DRMerror << "    local_pos = " << local_pos << endln;
                        for (int dd = 0; dd < NIE; dd++)
                        {
                            H5DRMerror << "    Fk(" << dd << ") = " << Fk(dd) <<
                                       "  u_e(" << dd << ") = " << u_e(dd) <<
                                       "  Fm(" << dd << ") = " << Fm(dd) <<
                                       "  udd_e(" << dd << ") = " << udd_e(dd) << endln;
                        }
                        exit(-1);

                    }
                }
            }
        }
    }
    if (myrank == 0)
        H5DRMout << "ComputeDRMLoads.... Done.\n";



    if (DEBUG_DRM_FORCES)
    {
        char debugfilename[100];
        sprintf(debugfilename, "drmforces.%d.txt", myrank);
        FILE * fptr = fopen(debugfilename, "w");

        fprintf(fptr, "%f ", t);

        for (int i = 0; i < DRMForces.Size(); ++i)
        {
            fprintf(fptr, "%f ", DRMForces(i));
        }
        fprintf(fptr, "\n");

    }



    return true;
}




int
H5DRM::sendSelf(int commitTag, Channel & theChannel)
{

    H5DRMout << "sending filename: " << HDF5filename << endl;

    static Vector data(3);
    data(0) = cFactor;
    data(1) = crd_scale;
    data(2) = distance_tolerance;

    char drmfilename[H5DRM_MAX_FILENAME];
    strcpy(drmfilename, HDF5filename.c_str());
    Message filename_msg(drmfilename, H5DRM_MAX_FILENAME);

    if (theChannel.sendMsg(0, 0, filename_msg) < 0)
    {
        cerr << "H5DRM::sendSelf -- failed to send HDF5filename\n";
        return -1;
    }

    if (theChannel.sendVector(0, 0, data) < 0)
    {
        cerr << "H5DRM::sendSelf -- failed to send cFactor\n";
        return -1;
    }


    return 0;
}

int
H5DRM::recvSelf(int commitTag, Channel & theChannel,
                FEM_ObjectBroker & theBroker)
{
    H5DRMout << "receiving...\n";
    static Vector data(3);
    char drmfilename[H5DRM_MAX_FILENAME];
    Message filename_msg(drmfilename, H5DRM_MAX_FILENAME);

    // strcpy(drmfilename, HDF5filename.c_str());

    if (theChannel.recvMsg(0, 0, filename_msg) < 0)
    {
        cerr << "H5DRM::receiveSelf -- failed to receive HDF5filename\n";
        return -1;
    }

    if (theChannel.recvVector(0, 0, data) < 0)
    {
        cerr << "H5DRM::receiveSelf -- failed to receive cFactor\n";
        return -1;
    }
    cFactor = data(0);
    crd_scale = data(1);
    distance_tolerance = data(2);

    HDF5filename = drmfilename;
    H5DRMout << "received filename is " <<  drmfilename << "\n";

    return 0;
}

void
H5DRM::Print(ostream & s, int flag)
{

}



// method to obtain a blank copy of the LoadPattern
LoadPattern *
H5DRM::getCopy(void)
{
    return new H5DRM(this->getTag(), HDF5filename);
}



Vector *
H5DRM::getNodalLoad(int nodeTag, double time)
{


    return 0;
}




void H5DRM::node_matching_UsePlaneInfo(double d_tol, const ID & internal, const Matrix & xyz, const Vector & drmbox_x0, double & d_err, int & n_nodes_found)
{

    // char debugfilename[100];
    // sprintf(debugfilename, "debuguseplaneinfo.%d.txt", myrank);
    // FILE * fptr = fopen(debugfilename, "w");

    // Domain *theDomain = this->getDomain();
    // NodeIter& node_iter = theDomain->getNodes();
    // Node* node_ptr = 0;
    // int drmtag = NDRM_points;
    // while ((node_ptr = node_iter()) != 0)
    // {
    //     int tag = node_ptr->getTag();
    //     const Vector& node_xyz  =  node_ptr->getCrds();
    //     int i_min_plane = -1;
    //     double dmin = std::numeric_limits<double>::infinity();
    //     double xi1_min, xi2_min;
    //     fprintf(fptrdrm, "%d %f %f %f\n", ++drmtag, node_xyz[0] + drmbox_x0[0], node_xyz[1] + drmbox_x0[1], node_xyz[2] + drmbox_x0[2]);

    //     for (auto plane_it = std::begin(planes); plane_it != std::end(planes); ++plane_it)
    //     {
    //         Plane* plane = *plane_it;
    //         int plane_number = plane->getNumber();
    //         double xi1, xi2, d;
    //         plane->locate_point(node_xyz + drmbox_x0, xi1, xi2, d);
    //         if (fabs(d) < dmin)
    //         {
    //             dmin = d;
    //             i_min_plane = plane_number;
    //             xi1_min = xi1;
    //             xi2_min = xi2;
    //         }
    //     }
    //     if (i_min_plane == -1)
    //     {
    //         printf( "node_matching_UsePlaneInfo - Node # %05d @ (%4.2f, %4.2f, %4.2f) was not matched to any plane.\n", tag, node_xyz(0) + drmbox_x0(0), node_xyz(1) + drmbox_x0(1), node_xyz(2) + drmbox_x0(2));
    //     }
    //     if (fabs(dmin) < d_tol)
    //     {
    //         int i, j;
    //         Plane* plane = planes[i_min_plane];
    //         plane->get_ij_coordinates_to_point(xi1_min, xi2_min, i, j);
    //         if (i < 0 || j < 0)
    //         {
    //             printf( "node_matching_UsePlaneInfo - Node # %05d @ (%4.2f, %4.2f, %4.2f) \n", tag, node_xyz(0) + drmbox_x0(0), node_xyz(1) + drmbox_x0(1), node_xyz(2) + drmbox_x0(2));
    //             printf( "     at min distance d = %4.2f from plane %02d \n", dmin, i_min_plane);
    //             printf( "     xi1 = %4.2f xi2 = %4.2f \n", xi1_min, xi2_min);
    //             printf( "     i   = %4d xi2 = %4d \n", i, j);
    //         }
    //         else
    //         {
    //             int station_id = plane->getStationNumber(i, j);
    //             static Vector station_xyz(3);
    //             for (int dir = 0; dir < 3; ++dir)
    //                 station_xyz(dir) = xyz(station_id, dir);
    //             double this_d_err = (station_xyz - node_xyz - drmbox_x0).Norm();
    //             if (this_d_err > d_err)
    //                 d_err = this_d_err;
    //             nodetag2hdf5_pos.insert ( std::pair<int, int>(tag, station_id) );
    //             Nodes.insert(tag);
    //             ++n_nodes_found;
    //         }
    //     }
    //     else
    //     {
    //         Plane* plane = planes[i_min_plane];
    //         // int station_id = plane->getStationNumber(i, j);
    //         fprintf(fptr, "Node # %05d @ (%4.2f, %4.2f, %4.2f) rejected \n", tag, node_xyz(0) + drmbox_x0(0), node_xyz(1) + drmbox_x0(1), node_xyz(2) + drmbox_x0(2));
    //         fprintf(fptr, "     at distance d = %4.2f from plane %02d (xi1 = %5.3f, xi2 = %5.3f)\n", dmin, i_min_plane, xi1_min, xi2_min);
    //         plane->print(fptr);

    //     }
    // }
    // fclose(fptr);

    return;
}




void H5DRM::node_matching_OctTree(double d_tol, const ID & internal, const Matrix & xyz, const Vector & drmbox_x0, double & d_err, int & n_nodes_found)
{
    return;
}



void H5DRM::node_matching_BruteForce(double d_tol, const ID & internal, const Matrix & xyz, const Vector & drmbox_x0, double & d_err, int & n_nodes_found)
{

    H5DRMout << "node_matching_BruteForce - Begin!\n";

    char debugfilename[100];
    sprintf(debugfilename, "debugdrmbruteforce.%d.txt", myrank);
    FILE * fptrdrm;
    if (DEBUG_NODE_MATCHING)
        fptrdrm = fopen(debugfilename, "w");

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
            fprintf(fptrdrm, "%d %f %f %f\n", ++drmtag, node_xyz[0] + drmbox_x0[0], node_xyz[1] + drmbox_x0[1], node_xyz[2] + drmbox_x0[2]);
        Vector station_xyz(3);
        for (int ii = 0; ii < xyz.noRows(); ++ii)
        {
            station_xyz(0) = xyz(ii, 0);
            station_xyz(1) = xyz(ii, 1);
            station_xyz(2) = xyz(ii, 2);
            double d = (node_xyz + drmbox_x0 - station_xyz).Norm();
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
                station_xyz(dir) = xyz(station_id, dir);
            double this_d_err = (station_xyz - node_xyz - drmbox_x0).Norm();
            if (this_d_err > d_err)
                d_err = this_d_err;
            nodetag2station_id.insert ( std::pair<int, int>(tag, station_id) );
            nodetag2local_pos.insert (  std::pair<int, int>(tag, local_pos) );
            Nodes[local_pos] = tag;
            IsBoundary[local_pos] = internal(station_id);
            ++n_nodes_found;
            ++local_pos;

        }
        else
        {
            if (DEBUG_NODE_MATCHING)
                fprintf(fptrdrm, "Node # %05d @ (%4.2f, %4.2f, %4.2f) rejected \n", tag, node_xyz(0) + drmbox_x0(0), node_xyz(1) + drmbox_x0(1), node_xyz(2) + drmbox_x0(2));
        }
    }
    if (DEBUG_NODE_MATCHING)
        fclose(fptrdrm);

    H5DRMout << "node_matching_BruteForce - End!\n";

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
    hid_t id_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t id_dataspace = H5Dget_space(id_dataset);
    int ndims = H5Sget_simple_extent_ndims( id_dataspace);
    if (ndims != 1)
    {
        opserr << "failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "read_int_dataset_into_id - array dimension should be 1\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim;
    hsize_t maxdim;
    H5Sget_simple_extent_dims(id_dataspace, &dim, &maxdim);
    H5Sselect_all( id_dataspace );
    hid_t id_memspace  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t id_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    int d[dim];

    H5Dread( id_dataset, H5T_NATIVE_INT, id_memspace, id_dataspace, id_xfer_plist,  d);
    result.resize(dim);
    for (hsize_t i = 0; i < dim; ++i)
    {
        result(i) = d[i];
    }

    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_dataset);

    return true;
}

bool read_double_dataset_into_vector(const hid_t & h5drm_dataset, std::string dataset_name, Vector & result)
{
    hid_t id_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t id_dataspace = H5Dget_space(id_dataset);
    int ndims = H5Sget_simple_extent_ndims( id_dataspace);
    if (ndims != 1)
    {
        opserr << "failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "read_double_dataset_into_vector - array dimension should be 1.\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim;
    hsize_t maxdim;
    H5Sget_simple_extent_dims(id_dataspace, &dim, &maxdim);
    H5Sselect_all( id_dataspace );
    hid_t id_memspace  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t id_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    double d[dim];

    H5Dread( id_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, id_xfer_plist,  d);
    result.resize(dim);
    for (hsize_t i = 0; i < dim; ++i)
    {
        result(i) = d[i];
    }

    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_dataset);

    return true;
}


bool read_scalar_double_dataset_into_double(const hid_t & h5drm_dataset, std::string dataset_name, double & result)
{
    hid_t id_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t id_dataspace = H5Dget_space(id_dataset);
    int ndims = H5Sget_simple_extent_ndims( id_dataspace);
    if (ndims != 0)
    {
        opserr << "failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "read_scalar_double_dataset_into_double - array dimension should be 0\n";
        Vector error(-1);
        return false;
    }
    H5Sselect_all( id_dataspace );
    hsize_t dim = 1;
    hid_t id_memspace  = H5Screate_simple(1, &dim, 0);       // create dataspace of memory
    hid_t id_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    H5Dread( id_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, id_xfer_plist,  &result);

    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_dataset);

    return true;
}



bool read_double_dataset_into_matrix(const hid_t & h5drm_dataset, std::string dataset_name, Matrix & result)
{
    hid_t id_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t id_dataspace = H5Dget_space(id_dataset);
    int ndims = H5Sget_simple_extent_ndims( id_dataspace);
    if (ndims != 2)
    {
        opserr << "failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "read_double_dataset_into_matrix - array dimension should be 2\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim[2];
    hsize_t maxdim[2];
    H5Sget_simple_extent_dims(id_dataspace, dim, maxdim);
    H5Sselect_all( id_dataspace );
    hid_t id_memspace  = H5Screate_simple(2, dim, 0);       // create dataspace of memory
    hid_t id_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    double d[dim[0]][dim[1]];

    H5Dread( id_dataset, H5T_NATIVE_DOUBLE, id_memspace, id_dataspace, id_xfer_plist,  d);
    result.resize(dim[0], dim[1]);
    for (hsize_t i = 0; i < dim[0]; ++i)
    {
        for (hsize_t j = 0; j < dim[1]; ++j)
        {
            result(i, j) = d[i][j];
        }
    }

    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_dataset);

    return true;
}


bool read_int_dataset_into_array(const hid_t & h5drm_dataset, std::string dataset_name, int *& result)
{
    hid_t id_dataset = H5Dopen(h5drm_dataset, dataset_name.c_str(), H5P_DEFAULT);
    hid_t id_dataspace = H5Dget_space(id_dataset);
    int ndims = H5Sget_simple_extent_ndims( id_dataspace);
    if (ndims != 2)
    {
        opserr << "failure trying to open dataset: " << dataset_name.c_str() << endln;
        opserr << "read_int_dataset_into_array - array dimension should be 2\n";
        Vector error(-1);
        return false;
    }
    hsize_t dim[2];
    hsize_t maxdim[2];
    H5Sget_simple_extent_dims(id_dataspace, dim, maxdim);
    H5Sselect_all( id_dataspace );
    hid_t id_memspace  = H5Screate_simple(2, dim, 0);       // create dataspace of memory
    hid_t id_xfer_plist = H5Pcreate(H5P_DATASET_XFER);

    result = new int [dim[0]*dim[1]];

    H5Dread( id_dataset, H5T_NATIVE_INT, id_memspace, id_dataspace, id_xfer_plist,  result);

    H5Sclose(id_dataspace);
    H5Sclose(id_memspace);
    H5Dclose(id_dataset);

    return true;
}














Plane::Plane(const hid_t & id_h5drm_file, int plane_number, double crd_scale) :
    number(plane_number),
    internal(false),
    stations(0),
    v0(3), v1(3), v2(3)
{
    char plane_group_name[H5DRM_MAX_STRINGSIZE];
    sprintf(plane_group_name, "/DRM_Planes/plane_%02d", plane_number);
    hid_t id_plane_grp = H5Gopen(id_h5drm_file, plane_group_name, H5P_DEFAULT);

    if ( !read_double_dataset_into_vector(id_plane_grp, "v0", v0));
    if ( !read_double_dataset_into_vector(id_plane_grp, "v1", v1));
    if ( !read_double_dataset_into_vector(id_plane_grp, "v2", v2));
    if ( !read_double_dataset_into_vector(id_plane_grp, "xi1", xi1));
    if ( !read_double_dataset_into_vector(id_plane_grp, "xi2", xi2));
    if ( !read_int_dataset_into_array(id_plane_grp, "stations", stations));

    v0 *= crd_scale;  //Convert scales km -> m
    v1 *= crd_scale;  //Convert scales km -> m
    v2 *= crd_scale;  //Convert scales km -> m

    convert_h5drmcrd_to_ops_crd(v0);
    convert_h5drmcrd_to_ops_crd(v1);
    convert_h5drmcrd_to_ops_crd(v2);
}


Plane::~Plane()
{
    if (stations != 0)
    {
        delete stations;
    }
}


bool Plane::locate_point(const Vector & x, double & xi1_, double & xi2_, double & distance) const
{
    static Matrix A(2, 2);
    static Vector b(2);
    static Vector xi(2);

    A(0, 0) = v1 ^ v1;
    A(0, 1) = v1 ^ v2;
    A(1, 0) = v2 ^ v1;
    A(1, 1) = v1 ^ v1;

    Vector xperp = x - v0;
    b(0) = xperp ^ v1;
    b(1) = xperp ^ v2;

    A.Solve(b, xi);

    xi1_ = xi(0);
    xi2_ = xi(1);

    int i, j;
    get_ij_coordinates_to_point(xi1_, xi2_, i, j);

    xi1_ = xi1(i);
    xi2_ = xi2(j);
    distance = (x - (v0 + xi1(i) * v1 + xi2(j) * v2 )).Norm();

    return true;
}


bool Plane::get_ij_coordinates_to_point(double x1, double x2, int& i, int& j) const
{
    int N1 = xi1.Size();
    int N2 = xi2.Size();
    double xi1start = xi1(0);
    double xi1end = xi1(N1 - 1);
    double xi2start = xi2(0);
    double xi2end = xi2(N2 - 1);
    i = (int) round( (x1 - xi1start) / (xi1end - xi1start) * (N1 - 1) );
    j = (int) round( (x2 - xi2start) / (xi2end - xi2start) * (N2 - 1) );

    i = i < 0     ? 0    : i;
    i = i > N1 - 1  ? N1 - 1 : i;
    j = j < 0     ? 0    : j;
    j = j > N2 - 1  ? N2 - 1 : j;

    return true;
}


int Plane::getNumber() const
{
    return number;
}


int Plane::getStationNumber(int i, int j) const
{
    int N1 = xi1.Size();
    int N2 = xi2.Size();

    if ( i < 0 || i >= N1 || j < 0 || j >= N2)
    {
        return -1;
    }

    if ( stations == 0)
    {
        opserr << "Plane::getStationNumber - stations uninitialized at Plane # "  << number << "\n";
        return -1;
    }
    else
    {
        return stations[N2 * i + j];
    }
}

void Plane::print(FILE * fptr) const
{
    fprintf(fptr, "    Plane # %d v0 = (%5.3f, %5.3f, %5.3f) v1 = (%5.3f, %5.3f, %5.3f) v2 = (%5.3f, %5.3f, %5.3f)\n", number, v0(0), v0(1), v0(2), v1(0), v1(1), v1(2), v2(0), v2(1), v2(2));
}



void convert_h5drmcrd_to_ops_crd(Vector & v )
{
    // static Vector v_tmp(3);
    // v_tmp = v;
    // v(0) = v_tmp(1);
    // v(1) = v_tmp(0);
    // v(2) = -v_tmp(2);
}

void convert_h5drmcrd_to_ops_crd(Matrix & xyz )
{
    // int nrows = xyz.noRows();
    // Matrix tmp_xyz(xyz);
    // for (int i = 0; i < nrows; ++i)
    // {
    //     xyz(i, 0) = tmp_xyz(i, 1);
    //     xyz(i, 1) = tmp_xyz(i, 0);
    //     xyz(i, 2) = -tmp_xyz(i, 2);
    // }
}

