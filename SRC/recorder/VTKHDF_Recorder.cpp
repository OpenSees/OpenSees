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
**   Frank McKenna (fmckenna@ce.berkeley.edu)                        **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                      **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                    **
**                                                                    **
** ****************************************************************** */

#include "VTKHDF_Recorder.h"
#include <sstream>
#include <elementAPI.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <NodeIter.h>
#include <Element.h>
#include <ElementIter.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Message.h>
#include <classTags.h>
#include <iostream>
#include "hdf5.h"

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif


// Static variable initialization
std::map<int,VTKHDF_Recorder::VtkType> VTKHDF_Recorder::vtktypes;



// Constructor
SilentStream::SilentStream() : OPS_Stream(0) {}

// Destructor
SilentStream::~SilentStream() {}




OutputDataHDF::OutputDataHDF()
{
    disp     = false; vel      = false; accel    = false; 
    reaction = false; mass     = false; unbalancedLoad = false;
    stress3D6 = false; strain3D6 = false;
    stress2D3 = false; strain2D3 = false;
}

OutputDataHDF &OutputDataHDF::operator=(const OutputDataHDF &other) 
{
    // Guard against self-assignment
    if (this == &other) {
        return *this;
    }
    // Copy all members
    disp       = other.disp;     
    vel        = other.vel;    
    accel      = other.accel;  
    reaction   = other.reaction; 
    mass       = other.mass;     
    unbalancedLoad = other.unbalancedLoad;
    stress3D6 = other.stress3D6;
    strain3D6 = other.strain3D6;
    stress2D3 = other.stress2D3;
    strain2D3 = other.strain2D3;

    
    return *this;
}



void* OPS_VTKHDF_Recorder()
{
    // Check minimum arguments
    int numdata = OPS_GetNumRemainingInputArgs();
    if(numdata < 1) {
        opserr << "WARNING: insufficient number of arguments\n";
        return 0;
    }

    // Get filename
    const char* name = OPS_GetString();

    // Initialize default values
    OutputDataHDF outputData;
    std::vector<VTKHDF_Recorder::EleData> eledata;
    double dT = 0.0;
    double rTolDt = 0.00001;

    // Parse remaining arguments
    numdata = OPS_GetNumRemainingInputArgs();
    while(numdata > 0) {
        const char* type = OPS_GetString();

        // Displacement options
        if(strcmp(type, "disp") == 0) {
            outputData.disp = true;
        }
        // Velocity options
        else if(strcmp(type, "vel") == 0) {
            outputData.vel = true;
        }
        // Acceleration options
        else if(strcmp(type, "accel") == 0) {
            outputData.accel = true;
        }
        // Reaction options
        else if(strcmp(type, "reaction") == 0) {
            outputData.reaction = true;
        }
        // Other response options
        else if(strcmp(type, "mass") == 0) {
            outputData.mass = true;
        } else if(strcmp(type, "unbalancedLoad") == 0) {
            outputData.unbalancedLoad = true;
        }

        else if (strcmp(type, "stress3D6") == 0) {
            outputData.stress3D6 = true;
        } else if (strcmp(type, "strain3D6") == 0) {
            outputData.strain3D6 = true;
        }
        else if (strcmp(type, "stress2D3") == 0) {
            outputData.stress2D3 = true;
        } else if(strcmp(type, "strain2D3") == 0) {
            outputData.strain2D3 = true;
        }
        // Time step options
        else if(strcmp(type, "-dT") == 0) {
            if(OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: needs dT\n";
                return 0;
            }
            numdata = 1;
            if(OPS_GetDoubleInput(&numdata, &dT) < 0) {
                opserr << "WARNING: failed to read dT\n";
                return 0;
            }
            dT = (dT < 0) ? 0 : dT;
        }
        else if(strcmp(type, "-rTolDt") == 0) {
            if(OPS_GetNumRemainingInputArgs() < 1) {
                opserr << "WARNING: needs rTolDt\n";
                return 0;
            }
            numdata = 1;
            if(OPS_GetDoubleInput(&numdata, &rTolDt) < 0) {
                opserr << "WARNING: failed to read rTolDt\n";
                return 0;
            }
            rTolDt = (rTolDt < 0) ? 0 : rTolDt;
        }

        numdata = OPS_GetNumRemainingInputArgs();
    }
    // // print summary of parmeters for the user
    // opserr << "VTKHDF_Recorder: " << name << " -dT " << dT << " -rTolDt " << rTolDt << endln;
    // opserr << "disp " << (outputData.disp ? "true" : "false") << endln;
    // opserr << "vel " << (outputData.vel ? "true" : "false") << endln;
    // opserr << "accel " << (outputData.accel ? "true" : "false") << endln;
    // opserr << "reaction " << (outputData.reaction ? "true" : "false") << endln;
    // opserr << "mass " << (outputData.mass ? "true" : "false") << endln; 
    // opserr << "unbalancedLoad " << (outputData.unbalancedLoad ? "true" : "false") << endln;
    // opserr << "stress3D6 " << (outputData.stress3D6 ? "true" : "false") << endln;
    // opserr << "strain3D6 " << (outputData.strain3D6 ? "true" : "false") << endln;
    // opserr << "stress2D3 " << (outputData.stress2D3 ? "true" : "false") << endln;
    // opserr << "strain2D3 " << (outputData.strain2D3 ? "true" : "false") << endln;

    // Create the recorder
    return new VTKHDF_Recorder(name, outputData, eledata, dT, rTolDt);
}

VTKHDF_Recorder::VTKHDF_Recorder(const char *inputName, 
                                const OutputDataHDF& outData,
                                const std::vector<EleData>& edata, 
                                double dt, double rTolDt)
    :Recorder(RECORDER_TAGS_VTKHDF_Recorder), 
     // Initialize member variables
     name(nullptr),
     outputData(outData),
     theDomain(0),
     nextTimeStampToRecord(0.0),
     deltaT(dt),
     relDeltaTTol(rTolDt),
     counter(0),
     initializationDone(false),
     initDone(false)
{

   
    name = new char[strlen(inputName) + 1];
    strcpy(name, inputName);

    
    VTKHDF_Recorder::setVTKType();    // -----------------------
    // 1) Create (or overwrite) the HDF5 file using the C API with SWMR support
    // -----------------------
    // Create file access property list for SWMR (Single Writer Multiple Reader)
    hid_t fapl_id = H5Pcreate(H5P_FILE_ACCESS);
    if (fapl_id < 0) {
        opserr << "Error: Could not create file access property list for " << name << endln;
        return;
    }
    
    // Enable latest version of the library format (required for SWMR)
    herr_t ret = H5Pset_libver_bounds(fapl_id, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
    if (ret < 0) {
        opserr << "Error: Could not set library version bounds for " << name << endln;
        H5Pclose(fapl_id);
        return;
    }

    // Create the file with SWMR-compatible settings
    file_id = H5Fcreate(name, H5F_ACC_TRUNC, H5P_DEFAULT, fapl_id);
    if (file_id < 0) {
        opserr << "Error: Could not create HDF5 file " << name << endln;
        H5Pclose(fapl_id);
        return;
    }
    
    // Close the file temporarily to reopen it in SWMR write mode
    H5Fclose(file_id);
    
    // Reopen the file in SWMR write mode
    file_id = H5Fopen(name, H5F_ACC_RDWR | H5F_ACC_SWMR_WRITE, fapl_id);
    if (file_id < 0) {
        opserr << "Error: Could not reopen HDF5 file in SWMR write mode " << name << endln;
        H5Pclose(fapl_id);
        return;
    }
    
    // Close the file access property list
    H5Pclose(fapl_id);

    // 4) Create the group "/VTKHDF"
    group_id = H5Gcreate(file_id, "/VTKHDF",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (group_id < 0) {
        opserr << "Error: Could not create group '/VTKHDF' in HDF5 file " << name << endln;
        H5Fclose(file_id);
        return;
    }

    // -----------------------
    // 2) Write the "Version" attribute (two integers [2, 3])
    // -----------------------
    // Create a simple 1D dataspace for two integers
    hsize_t dims[1] = {2};
    hid_t version_dataspace_id = H5Screate_simple(1, dims, nullptr);
    if (version_dataspace_id < 0) {
        opserr << "Error: Could not create dataspace for 'Version' attribute in " 
               << name << endln;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    hid_t version_attr_id = H5Acreate(group_id,
                                      "Version",
                                      H5T_NATIVE_INT,
                                      version_dataspace_id,
                                      H5P_DEFAULT,
                                      H5P_DEFAULT);
    if (version_attr_id < 0) {
        opserr << "Error: Could not create attribute 'Version' in group '/VTKHDF' in file "
               << name << endln;
        H5Sclose(version_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    int version_data[2] = {2, 3};
    herr_t status = H5Awrite(version_attr_id, H5T_NATIVE_INT, version_data);
    if (status < 0) {
        opserr << "Error: Could not write 'Version' attribute in file " << name << endln;
        H5Aclose(version_attr_id);
        H5Sclose(version_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // Close version attribute and dataspace
    H5Aclose(version_attr_id);
    H5Sclose(version_dataspace_id);

    // -----------------------
    // 3) Write the "Type" attribute: "UnstructuredGrid"
    // -----------------------
    // Create a dataspace for a single (scalar) string
    hid_t type_dataspace_id = H5Screate(H5S_SCALAR);
    if (type_dataspace_id < 0) {
        opserr << "Error: Could not create dataspace for 'Type' attribute in file "
               << name << endln;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // Create a string datatype (fixed-length or variable-length).
    // Here we use a fixed-length approach with the length of the string.
    const char *typeStr = "UnstructuredGrid";
    size_t typeLen = strlen(typeStr);

    hid_t str_type_id = H5Tcopy(H5T_C_S1);
    if (str_type_id < 0) {
        opserr << "Error: Could not copy H5T_C_S1 for attribute 'Type' in file "
               << name << endln;
        H5Sclose(type_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }
    // Set the size of the string type (including space for '\0', if desired)
    status = H5Tset_size(str_type_id, typeLen);
    if (status < 0) {
        opserr << "Error: Could not set size for 'Type' string in file " << name << endln;
        H5Tclose(str_type_id);
        H5Sclose(type_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }
    // Decide how to handle strings shorter than typeLen (null-pad, etc.)
    // For example:
    H5Tset_strpad(str_type_id, H5T_STR_NULLTERM);

    // Create the attribute
    hid_t type_attr_id = H5Acreate(group_id,
                                   "Type",
                                   str_type_id,
                                   type_dataspace_id,
                                   H5P_DEFAULT,
                                   H5P_DEFAULT);
    if (type_attr_id < 0) {
        opserr << "Error: Could not create attribute 'Type' in group '/VTKHDF' in file "
               << name << endln;
        H5Tclose(str_type_id);
        H5Sclose(type_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // Write the string
    status = H5Awrite(type_attr_id, str_type_id, typeStr);
    if (status < 0) {
        opserr << "Error: Could not write attribute 'Type' in file " << name << endln;
        H5Aclose(type_attr_id);
        H5Tclose(str_type_id);
        H5Sclose(type_dataspace_id);
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // Close all resources
    H5Aclose(type_attr_id);
    H5Tclose(str_type_id);
    H5Sclose(type_dataspace_id);

    // -----------------------
    // 4) create pointData group "/VTKHDF/pointData"
    // -----------------------

    point_data_group = H5Gcreate(group_id, "/VTKHDF/PointData",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (point_data_group < 0) {
        opserr << "Error: Could not create group '/VTKHDF/pointData' in HDF5 file " << name << endln;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // -----------------------
    // 5) create Cells group "/VTKHDF/CellData"
    // -----------------------
    cell_data_group = H5Gcreate(group_id, "/VTKHDF/CellData",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (cell_data_group < 0) {
        opserr << "Error: Could not create group '/VTKHDF/CellData' in HDF5 file " << name << endln;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // -----------------------
    // 6) create Points dataset "/VTKHDF/Points" with shape (0,3) and resizable (unlimited,3)
    // -----------------------
    {
        hsize_t dims[2] = {0, 3};
        hsize_t maxdims[2] = {H5S_UNLIMITED, 3};
        hid_t space_id = H5Screate_simple(2, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'Points'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[2] = {100, 3};
        H5Pset_chunk(plist_id, 2, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "Points",
            H5T_NATIVE_DOUBLE,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'Points'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);

    }

    // -----------------------
    // 7) create Connectivity dataset "/VTKHDF/Connectivity" with shape (0,) and resizable (unlimited)
    // -----------------------
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'Connectivity'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {100};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "Connectivity",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'Connectivity'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);

    }


    // -----------------------
    // 8) create Offset dataset "/VTKHDF/Offsets" with shape (0,) and resizable (unlimited)
    // -----------------------
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'Offset'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {100};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "Offsets",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'Offset'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);

    }

    // --------------------------------
    // 9) create the dataset VTKHDF/Types with shape (0,) and resizable (unlimited)
    // --------------------------------
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'Types'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {100};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "Types",
            H5T_NATIVE_UCHAR,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'Types'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);

    }


    // --------------------------------
    // 10) create a steps group "/VTKHDF/steps"
    // --------------------------------
    steps_group = H5Gcreate(group_id, "/VTKHDF/Steps",
                               H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (steps_group < 0) {
        opserr << "Error: Could not create group '/VTKHDF/steps' in HDF5 file " << name << endln;
        H5Gclose(group_id);
        H5Fclose(file_id);
        return;
    }

    // --------------------------------
    // 11) create NSteps attribute for steps group and set it zero
    // --------------------------------
    // Create a dataspace for a single integer
    {
        // Create NSteps attribute
        const char* attr_name = "NSteps";
        hsize_t dims[1] = {1};
        hid_t nsteps_dataspace_id = H5Screate_simple(1, dims, nullptr);
        hid_t nsteps_attr_id = H5Acreate(steps_group,
                                        attr_name,
                                        H5T_NATIVE_INT,
                                        nsteps_dataspace_id,
                                        H5P_DEFAULT,
                                        H5P_DEFAULT);
        
        int nsteps_data = 0;
        herr_t status = H5Awrite(nsteps_attr_id, H5T_NATIVE_INT, &nsteps_data);

        // Check overall success and cleanup
        if (nsteps_dataspace_id < 0 || nsteps_attr_id < 0 || status < 0) {
            opserr << "Error: Failed to initialize 'NSteps' attribute in file " << name << endln;
            if (nsteps_attr_id >= 0) H5Aclose(nsteps_attr_id);
            if (nsteps_dataspace_id >= 0) H5Sclose(nsteps_dataspace_id);
            if (steps_group >= 0) H5Gclose(steps_group);
            if (group_id >= 0) H5Gclose(group_id);
            if (file_id >= 0) H5Fclose(file_id);
            return;
        }

        // Clean up handles
        H5Aclose(nsteps_attr_id);
        H5Sclose(nsteps_dataspace_id);
    }

    // --------------------------------
    // 12) 
    // create groups NumberOfPoints and NumberOfCells 
    // and NumberOfConnectivityIds with
    // resizable length zero shape
    // Create NumberOfPoints dataset
    // --------------------------------
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'NumberOfPoints'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "NumberOfPoints",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'NumberOfPoints'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // Create NumberOfCells dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'NumberOfCells'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "NumberOfCells",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'NumberOfCells'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // Create NumberOfConnectivityIds dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'NumberOfConnectivityIds'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            group_id, "NumberOfConnectivityIds",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'NumberOfConnectivityIds'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }


    // -----------------------
    // 13) 
    // create the dataset VTKHDF/Steps/PointOffsets  with shape (0) and resizable
    // create the dataset VTKHDF/Steps/Values  with shape (0) and resizable
    // create the dataset VTKHDF/Steps/PartOffsets with shape (0) and resizable
    // create the dataset VTKHDF/Steps/NumberOfParts with shape (0) and resizable
    // -----------------------
    // Create PointOffsets dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'PointOffsets'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "PointOffsets",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'PointOffsets'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // Create Values dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'Values'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "Values",
            H5T_NATIVE_DOUBLE,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'Values'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    // Create PartOffsets dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'PartOffsets'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {1};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "PartOffsets",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'PartOffsets'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    // Create NumberOfParts dataset
    {
        hsize_t dims[1] = {0};
        hsize_t maxdims[1] = {H5S_UNLIMITED};
        hid_t space_id = H5Screate_simple(1, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'NumberOfParts'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[1] = {100};
        H5Pset_chunk(plist_id, 1, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "NumberOfParts",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'NumberOfParts'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    // --------------------------------
    // 14)
    // create dataset VTKHDF/Steps/CellOffsets with shape (0,1) and resizable (unlimited,1)
    // create dataset VTKHDF/Steps/ConnectivityIdOffsets  with shape (0,1) and resizable (unlimited,1)
    // --------------------------------
    // Create CellOffsets dataset
    {
        hsize_t dims[2] = {0, 1};
        hsize_t maxdims[2] = {H5S_UNLIMITED, 1};
        hid_t space_id = H5Screate_simple(2, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'CellOffsets'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[2] = {1, 1};
        H5Pset_chunk(plist_id, 2, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "CellOffsets",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'CellOffsets'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    
    // Create ConnectivityIdOffsets dataset
    {
        hsize_t dims[2] = {0, 1};
        hsize_t maxdims[2] = {H5S_UNLIMITED, 1};
        hid_t space_id = H5Screate_simple(2, dims, maxdims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for 'ConnectivityIdOffsets'\n";
            return;
        }

        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        hsize_t chunk_dims[2] = {1, 1};
        H5Pset_chunk(plist_id, 2, chunk_dims);

        hid_t dset_id = H5Dcreate(
            steps_group, "ConnectivityIdOffsets",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );
        if (dset_id < 0) {
            opserr << "Error creating dataset 'ConnectivityIdOffsets'\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        H5Pclose(plist_id);
        H5Dclose(dset_id);
    }

    // --------------------------------
    // 15)
    // create group VTKHDF/Steps/PointDataOffsets 
    // create group VTKHDF/Steps/CellDataOffsets
    // create group VTKHDF/Steps/FieldDataOffsets
    // create group VTKHDF/Steps/FieldDataSizes
    // --------------------------------
    // Create PointDataOffsets group
    {
        point_data_offsets_group = H5Gcreate(steps_group, "PointDataOffsets",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (point_data_offsets_group < 0) {
            opserr << "Error: Could not create group '/VTKHDF/Steps/PointDataOffsets' in HDF5 file " << name << endln;
            H5Gclose(steps_group);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }
    }

    // Create CellDataOffsets group
    {
        cell_data_offsets_group = H5Gcreate(steps_group, "CellDataOffsets",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (cell_data_offsets_group < 0) {
            opserr << "Error: Could not create group '/VTKHDF/Steps/CellDataOffsets' in HDF5 file " << name << endln;
            H5Gclose(steps_group);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }
    }

    // Create FieldDataOffsets group
    {
        hid_t field_data_offsets_group = H5Gcreate(steps_group, "FieldDataOffsets",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (field_data_offsets_group < 0) {
            opserr << "Error: Could not create group '/VTKHDF/Steps/FieldDataOffsets' in HDF5 file " << name << endln;
            H5Gclose(steps_group);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }
    }

    // Create FieldDataSizes group
    {
        hid_t field_data_sizes_group = H5Gcreate(steps_group, "FieldDataSizes",
                                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        if (field_data_sizes_group < 0) {
            opserr << "Error: Could not create group '/VTKHDF/Steps/FieldDataSizes' in HDF5 file " << name << endln;
            H5Gclose(steps_group);
            H5Gclose(group_id);
            H5Fclose(file_id);
            return;
        }
    } 

    // --------------------------------
    // 16)
    // check if displacement is requested create displacement dataset under PointDataOffsets
    // Create dataset for displacement if requested
    // --------------------------------
    if (outputData.disp) {
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[1] = {0}; // Initial size (0 rows)
        hsize_t max_dims[1] = {H5S_UNLIMITED}; // Maximum size (unlimited rows)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(1, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for displacement data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for displacement dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[1] = {100};
        if (H5Pset_chunk(plist_id, 1, chunk_dims) < 0) {
            opserr << "Error setting chunk size for displacement dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_offsets_group, "displacement",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for displacement\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // check if displacement is requested create pointData group for displacement
    // Create group for displacement if requested
    if (outputData.disp) {
        
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[2] = {0, 3}; // Initial size (0 rows, 3 columns)
        hsize_t max_dims[2] = {H5S_UNLIMITED, 3}; // Maximum size (unlimited rows, fixed 3 columns)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(2, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for displacement data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for displacement dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[2] = {100, 3};
        if (H5Pset_chunk(plist_id, 2, chunk_dims) < 0) {
            opserr << "Error setting chunk size for displacement dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_group, "displacement",
            H5T_NATIVE_DOUBLE,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for displacement\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // check if the acceleration is requested create offset dataset under PointDataOffsets
    // Create dataset for acceleration if requested
    if (outputData.accel) {
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[1] = {0}; // Initial size (0 rows)
        hsize_t max_dims[1] = {H5S_UNLIMITED}; // Maximum size (unlimited rows)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(1, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for acceleration data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for acceleration dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[1] = {100};
        if (H5Pset_chunk(plist_id, 1, chunk_dims) < 0) {
            opserr << "Error setting chunk size for acceleration dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_offsets_group, "acceleration",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for acceleration\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // check if the acceleration is requested create pointData group for acceleration
    // Create group for acceleration if requested
    if (outputData.accel) {
        
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[2] = {0, 3}; // Initial size (0 rows, 3 columns)
        hsize_t max_dims[2] = {H5S_UNLIMITED, 3}; // Maximum size (unlimited rows, fixed 3 columns)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(2, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for acceleration data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for acceleration dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[2] = {100, 3};
        if (H5Pset_chunk(plist_id, 2, chunk_dims) < 0) {
            opserr << "Error setting chunk size for acceleration dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_group, "acceleration",
            H5T_NATIVE_DOUBLE,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for acceleration\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // check if the velocity is requested create offset dataset under PointDataOffsets
    // Create dataset for velocity if requested
    if (outputData.vel) {
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[1] = {0}; // Initial size (0 rows)
        hsize_t max_dims[1] = {H5S_UNLIMITED}; // Maximum size (unlimited rows)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(1, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for velocity data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for velocity dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[1] = {100};
        if (H5Pset_chunk(plist_id, 1, chunk_dims) < 0) {
            opserr << "Error setting chunk size for velocity dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_offsets_group, "velocity",
            H5T_NATIVE_INT,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for velocity\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }

    // check if the velocity is requested create pointData group for velocity
    // Create group for velocity if requested
    if (outputData.vel) {
        
        // Define initial dimensions as 0 since we don't know the number of nodes yet
        hsize_t dims[2] = {0, 3}; // Initial size (0 rows, 3 columns)
        hsize_t max_dims[2] = {H5S_UNLIMITED, 3}; // Maximum size (unlimited rows, fixed 3 columns)

        // Create a dataspace with these dimensions
        hid_t space_id = H5Screate_simple(2, dims, max_dims);
        if (space_id < 0) {
            opserr << "Error creating dataspace for velocity data\n";
            return;
        }

        // Create a property list for chunked storage
        hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
        if (plist_id < 0) {
            opserr << "Error creating property list for velocity dataset\n";
            H5Sclose(space_id);
            return;
        }

        // Define the chunk size (e.g., 100 rows per chunk)
        hsize_t chunk_dims[2] = {100, 3};
        if (H5Pset_chunk(plist_id, 2, chunk_dims) < 0) {
            opserr << "Error setting chunk size for velocity dataset\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Create the dataset with chunked storage
        hid_t dset_id = H5Dcreate(
            point_data_group, "velocity",
            H5T_NATIVE_DOUBLE,
            space_id,
            H5P_DEFAULT, plist_id, H5P_DEFAULT
        );

        if (dset_id < 0) {
            opserr << "Error creating dataset for velocity\n";
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return;
        }

        // Clean up
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        H5Sclose(space_id);
    }



    if (outputData.stress3D6) {
        int res = 0;
        // Create 3DStress6 dataset
        const char* dataset_name = "stress3D6";
        hsize_t ncols = 6;
        hsize_t chunk_rows = 500;
        hid_t type_id = H5T_NATIVE_DOUBLE;
        hsize_t max_rows = H5S_UNLIMITED;
        res = this->createDataset(cell_data_group, dataset_name, ncols, chunk_rows, type_id, max_rows);


        // Create 3DStress6 offset dataset
        dataset_name = "stress3D6";
        ncols = 1;
        chunk_rows = 500;
        type_id = H5T_NATIVE_INT;
        max_rows = H5S_UNLIMITED;
        // res += this->createDataset(cell_data_offsets_group, dataset_name, ncols, chunk_rows, type_id, max_rows);
        res += this->createOffsetDataset(cell_data_offsets_group, dataset_name, chunk_rows, type_id, max_rows);

        current_Stress3D6Offset = 0;
        if (res < 0) {
            opserr << "Error creating 3DStress6 datasets\n";
            return;
        }
    }


    if (outputData.strain3D6) {
        int res = 0;
        // Create 3DStrain6 dataset
        const char* dataset_name = "strain3D6";
        hsize_t ncols = 6;
        hsize_t chunk_rows = 500;
        hid_t type_id = H5T_NATIVE_DOUBLE;
        hsize_t max_rows = H5S_UNLIMITED;
        res = this->createDataset(cell_data_group, dataset_name, ncols, chunk_rows, type_id, max_rows);

        // Create 3DStrain6 offset dataset
        dataset_name = "strain3D6";
        ncols = 1;
        chunk_rows = 500;
        type_id = H5T_NATIVE_INT;
        max_rows = H5S_UNLIMITED;
        // res += this->createDataset(cell_data_offsets_group, dataset_name, ncols, chunk_rows, type_id, max_rows);
        res += this->createOffsetDataset(cell_data_offsets_group, dataset_name, chunk_rows, type_id, max_rows);

        current_Strain3D6Offset = 0;
        if (res < 0) {
            opserr << "Error creating 3DStrain6 datasets\n";
            return;
        }
    }
    if (outputData.stress2D3) {
        int res = 0;
        // Create 2DStress3 dataset
        const char* dataset_name = "stress2D3";
        hsize_t ncols = 3;
        hsize_t chunk_rows = 500;
        hid_t type_id = H5T_NATIVE_DOUBLE;
        hsize_t max_rows = H5S_UNLIMITED;
        res = this->createDataset(cell_data_group, dataset_name, ncols, chunk_rows, type_id, max_rows);

        // Create 2DStress3 offset dataset
        dataset_name = "stress2D3";
        ncols = 1;
        chunk_rows = 500;
        type_id = H5T_NATIVE_INT;
        max_rows = H5S_UNLIMITED;
        res += this->createOffsetDataset(cell_data_offsets_group, dataset_name, chunk_rows, type_id, max_rows);


        current_Stress2D3Offset = 0;
        if (res < 0) {
            opserr << "Error creating 2DStress3 datasets\n";
            return;
        }
    }

    if (outputData.strain2D3) {
        int res = 0;
        // Create 2DStrain3 dataset
        const char* dataset_name = "strain2D3";
        hsize_t ncols = 3;
        hsize_t chunk_rows = 500;
        hid_t type_id = H5T_NATIVE_DOUBLE;
        hsize_t max_rows = H5S_UNLIMITED;
        res = this->createDataset(cell_data_group, dataset_name, ncols, chunk_rows, type_id, max_rows);

        // Create 2DStrain3 offset dataset
        dataset_name = "strain2D3";
        ncols = 1;
        chunk_rows = 500;
        type_id = H5T_NATIVE_INT;
        max_rows = H5S_UNLIMITED;
        res += this->createOffsetDataset(cell_data_offsets_group, dataset_name, chunk_rows, type_id, max_rows);

        current_Strain2D3Offset = 0;
        if (res < 0) {
            opserr << "Error creating 2DStrain3 datasets\n";
            return;
        }
    }





    current_PointOffset = 0;
    current_CellOffset = 0;
    current_ConnectivityIdOffset = 0;
    current_PartOffset = 0;

    next_PointOffset = 0;
    next_CellOffset = 0;
    next_ConnectivityIdOffset = 0;
    next_PartOffset = 0;


    initDone = false;
    numSteps = 0;    // opserr << "initilization done\n";
    CurrentDispOffset = 0;
    CurrentVelOffset = 0;
    CurrentAccelOffset = 0;

    // Enable SWMR mode for concurrent reading
    if (H5Fstart_swmr_write(file_id) < 0) {
        opserr << "Warning: Could not enable SWMR write mode for " << name << endln;
        opserr << "         File will still be created but concurrent reading may not work optimally" << endln;
    }

}





VTKHDF_Recorder::VTKHDF_Recorder()
  :Recorder(RECORDER_TAGS_VTKHDF_Recorder),
   theDomain(0),
   nextTimeStampToRecord(0.0),
   deltaT(0.0),
   relDeltaTTol(0.00001),
   counter(0),
   initializationDone(false)
{
  name = NULL;
  VTKHDF_Recorder::setVTKType();
  initDone = false;

}


int VTKHDF_Recorder::writeMesh() {
    // 1) Check domain
    if (theDomain == nullptr) {
        opserr << "WARNING: No domain found -- VTKHDF_Recorder::writeMesh\n";
        return -1;
    }

    // Clear old data structures (if you store them as class members)
    theNodeMapping.clear();
    theEleMapping.clear();
    theNodeTags.clear();
    theEleTags.clear();
    theEleClassTags.clear();
    theEleVtkTags.clear();
    theEleVtkOffsets.clear();

    // 2) Gather node info
    NodeIter &theNodes = theDomain->getNodes();
    Node *theNode;
    numNode = 0;  
    maxNDM = 0;   // max spatial dimension
    maxNDF = 0;   // max degrees of freedom

    while ((theNode = theNodes()) != nullptr) {
        int nodeTag = theNode->getTag();

        const Vector &crd = theNode->getCrds();
        if (crd.Size() > maxNDM) {
            maxNDM = crd.Size();
        }

        const Vector &disp = theNode->getTrialDisp();
        if (disp.Size() > maxNDF) {
            maxNDF = disp.Size();
        }

        // Build mapping (nodeTag -> index)
        theNodeMapping[nodeTag] = numNode;
        theNodeTags.push_back(nodeTag);
        numNode++;
    }

    // 3) Gather element info
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    numElement = 0;
    numConnectivityIds = 0;
    
    theEleVtkOffsets.push_back(numConnectivityIds);
    while ((theElement = theElements()) != nullptr) {
        int eleTag   = theElement->getTag();
        int classTag = theElement->getClassTag();

        // Check if recognized by our VTK type map
        auto it = vtktypes.find(classTag);
        if (it != vtktypes.end()) {
            int vtkType = it->second;

            theEleMapping[eleTag]   = numElement;
            theEleTags.push_back(eleTag);
            theEleClassTags.push_back(classTag);
            theEleVtkTags.push_back(vtkType);

            const ID &extNodes = theElement->getExternalNodes();
            int nN = extNodes.Size();

            // Offsets array is cumulative
            numConnectivityIds += nN; 
            theEleVtkOffsets.push_back(numConnectivityIds);

            numElement++;
        } else {
            // Possibly skip subdomains or warn for unrecognized element
            if (classTag != ELE_TAG_Subdomain) {
                opserr << "VTKHDF_Recorder::writeMesh - unknown vtk type for element "
                       << eleTag << " (classTag=" << classTag << ")" << endln;
            }
        }
    }

    // -------------------------------------------------------------
    // 4) Build arrays for HDF5
    // -------------------------------------------------------------

    // (A) Points => shape (numNode,3)
    std::vector<double> pointsData(numNode * 3, 0.0);
    {
        NodeIter &theNodes2 = theDomain->getNodes();
        while ((theNode = theNodes2()) != nullptr) {
            int nodeTag         = theNode->getTag();
            int mappedIndex     = theNodeMapping[nodeTag];
            const Vector &crd   = theNode->getCrds();

            pointsData[3 * mappedIndex + 0] = (crd.Size() > 0) ? crd(0) : 0.0;
            pointsData[3 * mappedIndex + 1] = (crd.Size() > 1) ? crd(1) : 0.0;
            pointsData[3 * mappedIndex + 2] = (crd.Size() > 2) ? crd(2) : 0.0;
        }
    }

    // (B) Connectivity => shape (numConnectivityIds)
    std::vector<int> connectivity(numConnectivityIds);
    {
        ElementIter &theElements2 = theDomain->getElements();
        int idx = 0;
        while ((theElement = theElements2()) != nullptr) {
            int classTag = theElement->getClassTag();
            if (vtktypes.find(classTag) != vtktypes.end()) {
                const ID &extNodes = theElement->getExternalNodes();
                for (int i = 0; i < extNodes.Size(); i++) {
                    int nodeTag = extNodes(i);
                    connectivity[idx++] = static_cast<int>(theNodeMapping[nodeTag]);
                }
            }
        }
    }

    // (C) Offsets => shape (numElement+1,1) so that ParaView can read it
    std::vector<int> offsets(theEleVtkOffsets.begin(), theEleVtkOffsets.end()); 


    // (D) Types => shape (numElement)
    std::vector<unsigned char> types(theEleVtkTags.begin(), theEleVtkTags.end());

    // opserr << " Gathering data done\n";
    // -------------------------------------------------------------
    // 5) Open group "/VTKHDF" from the already opened file_id
    // -------------------------------------------------------------
    // Make sure file_id is still valid (i.e., not closed in constructor).
    hid_t vtkhdfGroup = H5Gopen(file_id, "/VTKHDF", H5P_DEFAULT);
    if (vtkhdfGroup < 0) {
        opserr << "Error: Could not open group '/VTKHDF' in file to write mesh.\n";
        return -1;
    }

    // -------------------------------------------------------------
    // 6) resize the Point and add the points data
    // -------------------------------------------------------------
    {
        // Open the existing dataset
        hid_t dset_id = H5Dopen(vtkhdfGroup, "Points", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'Points'\n";
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'Points'\n";
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[2];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[2] = {current_dims[0] + static_cast<hsize_t>(numNode), 3};
        
        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'Points'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'Points'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[2] = {current_dims[0], 0};  // Start from the end of existing data
        hsize_t count[2] = {static_cast<hsize_t>(numNode), 3};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(2, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'Points'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, new_space_id, H5P_DEFAULT, pointsData.data()) < 0) {
            opserr << "Error writing data to 'Points'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);

    }
    // opserr << "Points done\n";
    

    // -------------------------------------------------------------
    // 7) resize the Connectivity and add the connectivity data
    // -------------------------------------------------------------
    {
        // Open the existing dataset
        hid_t dset_id = H5Dopen(vtkhdfGroup, "Connectivity", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'Connectivity'\n";
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'Connectivity'\n";
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + static_cast<hsize_t>(numConnectivityIds)};
        
        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'Connectivity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'Connectivity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {static_cast<hsize_t>(numConnectivityIds)};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'Connectivity'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space_id, new_space_id, H5P_DEFAULT, connectivity.data()) < 0) {
            opserr << "Error writing data to 'Connectivity'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }
    // opserr << "Connectivity done\n";

    // -------------------------------------------------------------
    // 8) resize the Offsets and add the offsets data
    // -------------------------------------------------------------
    {
        // Open the existing dataset
        hid_t dset_id = H5Dopen(vtkhdfGroup, "Offsets", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'Offsets'\n";
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'Offsets'\n";
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + static_cast<hsize_t>(numElement+1)};
        
        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'Offsets'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'Offsets'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {static_cast<hsize_t>(numElement+1)};  // Size of new data

        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'Offsets'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space_id, new_space_id, H5P_DEFAULT, offsets.data()) < 0) {
            opserr << "Error writing data to 'Offsets'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }
    // opserr << "Offsets done\n";


    // -------------------------------------------------------------
    // 9) resize the Types and add the types data
    // -------------------------------------------------------------
    {
        // Open the existing dataset
        hid_t dset_id = H5Dopen(vtkhdfGroup, "Types", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'Types'\n";
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'Types'\n";
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + static_cast<hsize_t>(numElement)};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'Types'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'Types'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {static_cast<hsize_t>(numElement)};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'Types'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            H5Gclose(vtkhdfGroup);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_UCHAR, mem_space_id, new_space_id, H5P_DEFAULT, types.data()) < 0) {
            opserr << "Error writing data to 'Types'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }
    // opserr << "Types done\n";
    H5Gclose(vtkhdfGroup);

    // -------------------------------------------------------------
    // 10 ) update the current offsets
    // -------------------------------------------------------------
    // for having a changing mesh you may need to adjust this ofsset
    // espically current_ConnectivityIdOffset and current_CellOffset
    // right now we are assuming that the mesh is not changing


    current_PointOffset = next_PointOffset;
    current_CellOffset = next_CellOffset;
    current_ConnectivityIdOffset = next_ConnectivityIdOffset;
    current_PartOffset = next_PartOffset;

    next_PointOffset = current_PointOffset + numNode;
    next_CellOffset = current_CellOffset + 0;
    next_ConnectivityIdOffset = current_ConnectivityIdOffset + 0;
    next_PartOffset = current_PartOffset + 0;


    // -------------------------------------------------------------
    // 11 ) Edit specific recorders information
    // -------------------------------------------------------------    
    // 3DStress6

    if (outputData.stress3D6) {
        opserr << "Setting up stress3D6 responses\n";
        stress3D6Responses.clear();
        stress3D6Responses.reserve(numElement);
        for (auto i: theEleTags) {
            Element *theElement = theDomain->getElement(i);
            SilentStream theOutputHandler;
            const char* args[] = {"stress3D6"};
            Response *theResponse = theElement->setResponse(args, 1, theOutputHandler);
            stress3D6Responses.push_back(theResponse);
        }
    }

    // 3DStrain6
    if (outputData.strain3D6) {
        opserr << "Setting up strain3D6 responses\n";
        strain3D6Responses.clear();
        strain3D6Responses.reserve(numElement);
        for (auto i: theEleTags) {
            Element *theElement = theDomain->getElement(i);
            SilentStream theOutputHandler;
            const char* args[] = {"strain3D6"};
            Response *theResponse = theElement->setResponse(args, 1, theOutputHandler);
            strain3D6Responses.push_back(theResponse);
        }
    }

    // 2DStress3
    if (outputData.stress2D3) {
        opserr << "Setting up stress2D3 responses\n";
        stress2D3Responses.clear();
        stress2D3Responses.reserve(numElement);
        for (auto i: theEleTags) {
            // const char* istring = std::to_string(i).c_str();
            // opserr << istring << endln;
            Element *theElement = theDomain->getElement(i);
            SilentStream theOutputHandler;
            const char* args[] = {"stress2D3"};
            Response *theResponse = theElement->setResponse(args, 1, theOutputHandler);
            stress2D3Responses.push_back(theResponse);
        }
    }

    // 2DStrain3
    if (outputData.strain2D3) {
        opserr << "Setting up strain2D3 responses\n";
        strain2D3Responses.clear();
        strain2D3Responses.reserve(numElement);
        for (auto i: theEleTags) {
            Element *theElement = theDomain->getElement(i);
            SilentStream theOutputHandler;
            const char* args[] = {"strain2D3"};
            Response *theResponse = theElement->setResponse(args, 1, theOutputHandler);
            strain2D3Responses.push_back(theResponse);
        }
    }



    // -------------------------------------------------------------
    // Mark that we've(re)written the mesh
    // -------------------------------------------------------------
    initDone = true;
    return 0;
}


int VTKHDF_Recorder::record(int commitTag, double timeStamp)
{
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    if (!initDone) {
        this->writeMesh();
    }

    if (deltaT == 0.0 || timeStamp - nextTimeStampToRecord >= -deltaT * relDeltaTTol) {
        if (this->writeStep(timeStamp) < 0) {
            opserr << "VTKHDF_Recorder::record() - writeStep() failed\n";
            return -1;
        }

        // // displacment
        if (outputData.disp) {
            // opserr<<"writeDisp"<<endln;
            if (this->writeDisp() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeDisp() failed\n";
                return -1;
            }
        }

        // // velocity
        if (outputData.vel) {
            // opserr<<"writeVel"<<endln;
            if (this->writeVel() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeVel() failed\n";
                return -1;
            }
        }

        // // acceleration
        if (outputData.accel) {
            // opserr<<"writeAccel"<<endln;
            if (this->writeAccel() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeAccel() failed\n";
                return -1;
            }
        }

        if (outputData.stress3D6) {
            // opserr<<"writeStress3D"<<endln;
            if (this->writeStrees3D6() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeStress3D() failed\n";
                return -1;
            }
        }

        if (outputData.strain3D6) {
            // opserr<<"writeStrain3D"<<endln;
            if (this->writeStrain3D6() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeStrain3D() failed\n";
                return -1;
            }
        }

        if (outputData.stress2D3) {
            // opserr<<"writeStress2D"<<endln;
            if (this->writeStress2D3() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeStress2D() failed\n";
                return -1;
            }
        }

        if (outputData.strain2D3) {
            // opserr<<"writeStrain2D"<<endln;
            if (this->writeStrain2D3() < 0) {
                opserr << "VTKHDF_Recorder::record() - writeStrain2D() failed\n";
                return -1;
            }
        }

        if (deltaT != 0.0) {
            nextTimeStampToRecord = timeStamp + deltaT;
        }
    }
    H5Fflush(file_id, H5F_SCOPE_GLOBAL);
    return 0;
}


int VTKHDF_Recorder::writeStress2D3(void) {

    std::vector<double> stressData(numElement * 3, 0.0);
    for (auto i: theEleTags) {
        Element *theEle = theDomain->getElement(i);
        int MappedIndex = theEleMapping[i];
        if (stress2D3Responses[MappedIndex] == 0) {
            continue;
        }
        stress2D3Responses[MappedIndex]->getResponse();
        Information &info = stress2D3Responses[MappedIndex]->getInformation();
        const Vector &stress = info.getData();
        for (int j = 0; j < 3; j++) {
            stressData[3 * MappedIndex + j] = stress(j);
        }
    }

    int res = 0;
    res += this->extendDataset(cell_data_group, "stress2D3", stressData.data(), H5T_NATIVE_DOUBLE, numElement, 3);
    res += this->extendOffsetDataset(cell_data_offsets_group, "stress2D3", &current_Stress2D3Offset, H5T_NATIVE_INT, 1);
    current_Stress2D3Offset += numElement;
    return res;
}

int VTKHDF_Recorder::writeStrain2D3(void) {

    std::vector<double> strainData(numElement * 3, 0.0);
    for (auto i: theEleTags) {
        Element *theEle = theDomain->getElement(i);
        int MappedIndex = theEleMapping[i];
        if (strain2D3Responses[MappedIndex] == 0) {
            continue;
        }
        strain2D3Responses[MappedIndex]->getResponse();
        Information &info = strain2D3Responses[MappedIndex]->getInformation();
        const Vector &strain = info.getData();
        for (int j = 0; j < 3; j++) {
            strainData[3 * MappedIndex + j] = strain(j);
        }
    }

    int res = 0;
    res += this->extendDataset(cell_data_group, "strain2D3", strainData.data(), H5T_NATIVE_DOUBLE, numElement, 3);
    res += this->extendOffsetDataset(cell_data_offsets_group, "strain2D3", &current_Strain2D3Offset, H5T_NATIVE_INT, 1);
    current_Strain2D3Offset += numElement;
    return res;
}

int VTKHDF_Recorder::writeStrees3D6(void) {

    std::vector<double> stressData(numElement * 6, 0.0);
    for (auto i: theEleTags) {
        Element *theEle = theDomain->getElement(i);
        int MappedIndex = theEleMapping[i];
        if (stress3D6Responses[MappedIndex]==0) {
            continue;
        }
        stress3D6Responses[MappedIndex]->getResponse();
        Information &info = stress3D6Responses[MappedIndex]->getInformation();
        const Vector &stress = info.getData();
        for (int j = 0; j < 6; j++) {
            stressData[6 * MappedIndex + j] = stress(j);
        }
    }

    int res = 0;
    res += this->extendDataset(cell_data_group, "stress3D6", stressData.data(), H5T_NATIVE_DOUBLE, numElement, 6);
    res += this->extendOffsetDataset(cell_data_offsets_group, "stress3D6", &current_Stress3D6Offset, H5T_NATIVE_INT, 1);
    current_Stress3D6Offset += numElement;
    return res;
}


int VTKHDF_Recorder::writeStrain3D6(void) {

    std::vector<double> strainData(numElement * 6, 0.0);
    for (auto i: theEleTags) {
        Element *theEle = theDomain->getElement(i);
        int MappedIndex = theEleMapping[i];
        if (strain3D6Responses[MappedIndex] == 0) {
            continue;
        }
        strain3D6Responses[MappedIndex]->getResponse();
        Information &info = strain3D6Responses[MappedIndex]->getInformation();
        const Vector &strain = info.getData();
        for (int j = 0; j < 6; j++) {
            strainData[6 * MappedIndex + j] = strain(j);
        }
    }

    int res = 0;
    res += this->extendDataset(cell_data_group, "strain3D6", strainData.data(), H5T_NATIVE_DOUBLE, numElement, 6);
    res += this->extendOffsetDataset(cell_data_offsets_group, "strain3D6", &current_Strain3D6Offset, H5T_NATIVE_INT, 1);
    current_Strain3D6Offset += numElement;
    return res;
}

int VTKHDF_Recorder::writeVel(void) 
{
    // write vel
    std::vector<double> velData(numNode * 3, 0.0);
    for (auto i : theNodeTags) {
        Node *theNode=theDomain->getNode(i);
        const Vector &vel=theNode->getVel();
        int mappedIndex     = theNodeMapping[i];
        int size = theNode->getCrds().Size();

        velData[3 * mappedIndex + 0] = (size > 0) ? vel(0) : 0.0;
        velData[3 * mappedIndex + 1] = (size > 1) ? vel(1) : 0.0;
        velData[3 * mappedIndex + 2] = (size > 2) ? vel(2) : 0.0;
    }

    {
        // Open dataset in point_data_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_group, "velocity", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'velocity'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'velocity'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[2];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[2] = {current_dims[0] + static_cast<hsize_t>(numNode), 3};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'velocity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'velocity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[2] = {current_dims[0], 0};  // Start from the end of existing data
        hsize_t count[2] = {static_cast<hsize_t>(numNode), 3};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(2, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'velocity'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, new_space_id, H5P_DEFAULT, velData.data()) < 0) {
            opserr << "Error writing data to 'velocity'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }

    {
        // Open dataset in point_data_offsets_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_offsets_group, "velocity", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'velocity'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'velocity'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + 1};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'velocity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'velocity'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {1};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'velocity'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space_id, new_space_id, H5P_DEFAULT, &CurrentVelOffset) < 0) {
            opserr << "Error writing data to 'velocity'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }

    CurrentVelOffset += numNode;
    return 0;
}


int VTKHDF_Recorder::writeAccel(void) 
{
    // write accel
    std::vector<double> accelData(numNode * 3, 0.0);
    for (auto i : theNodeTags) {
        Node *theNode=theDomain->getNode(i);
        const Vector &accel=theNode->getAccel();
        int mappedIndex     = theNodeMapping[i];
        int size = theNode->getCrds().Size();

        accelData[3 * mappedIndex + 0] = (size > 0) ? accel(0) : 0.0;
        accelData[3 * mappedIndex + 1] = (size > 1) ? accel(1) : 0.0;
        accelData[3 * mappedIndex + 2] = (size > 2) ? accel(2) : 0.0;
    }

    {
        // Open dataset in point_data_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_group, "acceleration", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'acceleration'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'acceleration'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[2];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[2] = {current_dims[0] + static_cast<hsize_t>(numNode), 3};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'acceleration'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'acceleration'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[2] = {current_dims[0], 0};  // Start from the end of existing data
        hsize_t count[2] = {static_cast<hsize_t>(numNode), 3};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(2, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'acceleration'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, new_space_id, H5P_DEFAULT, accelData.data()) < 0) {
            opserr << "Error writing data to 'acceleration'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }

    {
        // Open dataset in point_data_offsets_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_offsets_group, "acceleration", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'acceleration'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'acceleration'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + 1};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'acceleration'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'acceleration'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {1};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'acceleration'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space_id, new_space_id, H5P_DEFAULT, &CurrentAccelOffset) < 0) {
            opserr << "Error writing data to 'acceleration'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }

    CurrentAccelOffset += numNode;
    return 0;
    
}


int VTKHDF_Recorder::writeDisp(void)
{
    // write disp 
    std::vector<double> dispData(numNode * 3, 0.0);
    for (auto i : theNodeTags) {
        Node *theNode=theDomain->getNode(i);
        const Vector &disp=theNode->getDisp();
        int mappedIndex     = theNodeMapping[i];
        int size = theNode->getCrds().Size();

        dispData[3 * mappedIndex + 0] = (size > 0) ? disp(0) : 0.0;
        dispData[3 * mappedIndex + 1] = (size > 1) ? disp(1) : 0.0;
        dispData[3 * mappedIndex + 2] = (size > 2) ? disp(2) : 0.0;

    }

    {
        // Open dataset in point_data_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_group, "displacement", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'displacement'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'displacement'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[2];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[2] = {current_dims[0] + static_cast<hsize_t>(numNode), 3};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'displacement'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'displacement'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[2] = {current_dims[0], 0};  // Start from the end of existing data
        hsize_t count[2] = {static_cast<hsize_t>(numNode), 3};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(2, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'displacement'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space_id, new_space_id, H5P_DEFAULT, dispData.data()) < 0) {
            opserr << "Error writing data to 'displacement'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);        
    }


    // write dispofsset into 
    {
        // Open dataset in point_data_offsets_group and resize it to the new number of steps
        hid_t dset_id = H5Dopen(point_data_offsets_group, "displacement", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error opening dataset 'displacement'\n";
            return -1;
        }

        // Get the current dataspace and dimensions
        hid_t space_id = H5Dget_space(dset_id);
        if (space_id < 0) {
            opserr << "Error getting dataspace for 'displacement'\n";
            H5Dclose(dset_id);
            return -1;
        }

        // Get the current dimensions
        hsize_t current_dims[1];
        if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
            opserr << "Error getting current dimensions\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Calculate new dimensions (current rows + new rows, same columns)
        hsize_t new_dims[1] = {current_dims[0] + 1};

        // Extend the dataset
        if (H5Dset_extent(dset_id, new_dims) < 0) {
            opserr << "Error setting new extent for 'displacement'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Get the new dataspace
        hid_t new_space_id = H5Dget_space(dset_id);
        if (new_space_id < 0) {
            opserr << "Error getting new dataspace for 'displacement'\n";
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Select the hyperslab where we want to write new data
        hsize_t offset[1] = {current_dims[0]};  // Start from the end of existing data
        hsize_t count[1] = {1};  // Size of new data
        if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
            opserr << "Error selecting hyperslab\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }

        // Create memory space for new data
        hid_t mem_space_id = H5Screate_simple(1, count, NULL);
        if (mem_space_id < 0) {
            opserr << "Error creating memory space for 'displacement'\n";
            H5Sclose(new_space_id);
            H5Sclose(space_id);
            H5Dclose(dset_id);
            return -1;
        }


        // Write only the new data
        if (H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space_id, new_space_id, H5P_DEFAULT, &CurrentDispOffset) < 0) {
            opserr << "Error writing data to 'displacement'\n";
        }

        // Close all resources
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
    }

    // update the current offset
    CurrentDispOffset += numNode;
    return 0;
}



int VTKHDF_Recorder::writeStep(double timeStamp)
{
    // Update numSteps and attribute
    numSteps++;
    const char* attr_name = "NSteps";
    
    // Open and update NSteps attribute
    hid_t attr_id = H5Aopen(steps_group, attr_name, H5P_DEFAULT);
    if (attr_id < 0) {
        opserr << "Error: Failed to open 'NSteps' attribute\n";
        return -1;
    }
    
    H5Awrite(attr_id, H5T_NATIVE_INT, &numSteps);
    H5Aclose(attr_id);

    // Values for each dataset
    const int newNumberOfParts = 1;  // Always one rank per recorder
    const int newPointOffsets = current_PointOffset;
    const int newCellOffsets = current_CellOffset;
    const int newPartOffsets = current_PartOffset;
    const int newConnectivityIdOffsets = current_ConnectivityIdOffset;
    const double value = timeStamp;
    const int newNumPoints = numNode;
    const int newNumCells = numElement;
    const int newNumberOfConnectivityIds = numConnectivityIds;
    

    // Debug output (can be removed for production)
    // opserr << " Number of parts: " << newNumberOfParts << endln;
    // opserr << " Point offsets: " << newPointOffsets << endln;
    // opserr << " Cell offsets: " << newCellOffsets << endln;
    // opserr << " Part offsets: " << newPartOffsets << endln;
    // opserr << " Value: " << value << endln;
    // opserr << " Number of points: " << newNumPoints << endln;
    // opserr << " Number of cells: " << newNumCells << endln;
    // opserr << " Number of connectivity ids: " << newNumberOfConnectivityIds << endln;

    // Set up dimensions for dataset operations
    hsize_t newSize[1] = {static_cast<hsize_t>(numSteps)};
    hsize_t start[1] = {static_cast<hsize_t>(numSteps - 1)};
    hsize_t count[1] = {1};

    // NumberOfParts
    {
        hid_t dset_id = H5Dopen(steps_group, "NumberOfParts", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'NumberOfParts' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newNumberOfParts);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }
    // PointOffsets
    {
        hid_t dset_id = H5Dopen(steps_group, "PointOffsets", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'PointOffsets' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newPointOffsets);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }

    // PartOffsets
    {
        hid_t dset_id = H5Dopen(steps_group, "PartOffsets", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'PartOffsets' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newPartOffsets);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }

    // Values
    {
        hid_t dset_id = H5Dopen(steps_group, "Values", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'Values' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, mem_space, file_space, H5P_DEFAULT, &value);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }

    // NumberOfPoints
    {
        hid_t dset_id = H5Dopen(group_id, "NumberOfPoints", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'NumberOfPoints' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newNumPoints);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }
    // NumberOfCells
    {
        hid_t dset_id = H5Dopen(group_id, "NumberOfCells", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'NumberOfCells' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newNumCells);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }

    // NumberOfConnectivityIds
    {
        hid_t dset_id = H5Dopen(group_id, "NumberOfConnectivityIds", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'NumberOfConnectivityIds' dataset\n";
            return -1;
        }

        H5Dset_extent(dset_id, newSize);
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start, NULL, count, NULL);

        // Create memory space and write data
        hid_t mem_space = H5Screate_simple(1, count, NULL);
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newNumberOfConnectivityIds);

        // Clean up resources
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }

    // Handle CellOffsets dataset
    {
        // Set new dimensions for the current step
        hsize_t newSize2[2] = {static_cast<hsize_t>(numSteps), 1};
        hsize_t start2[2] = {static_cast<hsize_t>(numSteps - 1), 0};
        hsize_t count2[2] = {1, 1};

        // Open existing dataset
        hid_t dset_id = H5Dopen(steps_group, "CellOffsets", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'CellOffsets' dataset\n";
            return -1;
        }

        // Extend the dataset
        H5Dset_extent(dset_id, newSize2);

        // Get the file space and select the hyperslab
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start2, NULL, count2, NULL);

        // Create memory space for single value
        hid_t mem_space = H5Screate_simple(2, count2, NULL);

        // Write the cell offset value
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newCellOffsets);

        // Clean up
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }
    // Handle ConnectivityIdOffsets
    {
        // Set new dimensions for the current step
        hsize_t newSize2[2] = {static_cast<hsize_t>(numSteps), 1};
        hsize_t start2[2] = {static_cast<hsize_t>(numSteps - 1), 0};
        hsize_t count2[2] = {1, 1};

        // Open existing dataset
        hid_t dset_id = H5Dopen(steps_group, "ConnectivityIdOffsets", H5P_DEFAULT);
        if (dset_id < 0) {
            opserr << "Error: Failed to open 'ConnectivityIdOffsets' dataset\n";
            return -1;
        }

        // Extend the dataset
        H5Dset_extent(dset_id, newSize2);

        // Get the file space and select the hyperslab
        hid_t file_space = H5Dget_space(dset_id);
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, start2, NULL, count2, NULL);

        // Create memory space for single value
        hid_t mem_space = H5Screate_simple(2, count2, NULL);

        // Write the connectivity id offset value
        H5Dwrite(dset_id, H5T_NATIVE_INT, mem_space, file_space, H5P_DEFAULT, &newConnectivityIdOffsets);

        // Clean up
        H5Sclose(mem_space);
        H5Sclose(file_space);
        H5Dclose(dset_id);
    }
    
    
    return 0;
}


int VTKHDF_Recorder::restart(void)
{
    return 0;
}

int VTKHDF_Recorder::domainChanged(void)
{
    opserr << "Warning: VTKHDF_Recorder::domainChanged() - This recorder does not support domain changes yet.\n";
    opserr << "if you changed the domain, you need to create a new recorder.\n";
    opserr << "if you just used the domainChanged() mesthod without changing the domain, you can ignore this warning.\n";
    this->writeMesh();
    return -1;
}

int VTKHDF_Recorder::setDomain(Domain &domain)
{
  theDomain = &domain;
  return 0;
}

int VTKHDF_Recorder::sendSelf(int commitTag, Channel &theChannel)
{
    // Implementation for sending data
    return 0;
}

int VTKHDF_Recorder::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // Implementation for receiving data
    return 0;
}





// int VTKHDF_Recorder::writeStep_stress3D6() {
//     /*
//     This function is used to write the 3D stress data to the HDF5 file. The function writes the following:

//     1) collect the 3D stress data from the elements
//     2) Write the 3D stress data to the dataset "3DStress6" under group "/VTKHDF/CellData" in the HDF5 file
//     3) Write the offset to the dataset "3DStress6" under group "/VTKHDF/Steps/CellDataOffsets" in the HDF5 file
//     */

//     // 1) Collect the 3D stress data from the elements
//     std::vector<double> stress(numElement * 6, 0.0);
    
//     // std::vector<double> velData(numNode * 3, 0.0);
//     {
//         for (auto i : theEleTags) {
//             Element *theElement = theDomain->getElement(i);
//             int mappedIndex = theEleMapping[i];

//             int classTag = theElement->getClassTag();

//             if (classTag == ELE_TAG_STDQUADUP) {
//                 // stdQuadUP
//                 Vector stressVec(6);
//                 Information &eleInfo = 
//                 theElement->getResponse(1, stressVec);
//                 theElement->getResponse(1,)

//             }
        
//         }
        
// }






VTKHDF_Recorder::~VTKHDF_Recorder()
{
    // close the HDF5 file
    if (name) {
        delete [] name;
        name = nullptr;
    }

    if (file_id >= 0) { H5Fclose(file_id);}

    if (group_id >= 0) {H5Gclose(group_id);}

    if (point_data_group >= 0) {H5Gclose(point_data_group);}

    if (cell_data_group >= 0) {H5Gclose(cell_data_group);}

    if (steps_group >= 0) {H5Gclose(steps_group);}

    if (point_data_offsets_group >= 0) {H5Gclose(point_data_offsets_group);}

    if (cell_data_offsets_group >= 0) {H5Gclose(cell_data_offsets_group);}

}





void VTKHDF_Recorder::setVTKType()
{
    if (!vtktypes.empty()) {
	    return;
    }
    //    vtktypes[ELE_TAG_Subdomain] = VTK_POLY_VERTEX;
    vtktypes[ELEMENT_TAGS_WrapperElement] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElasticBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ModElasticBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_Beam2d] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d02] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d03] = VTK_LINE;
    vtktypes[ELE_TAG_beam2d04] = VTK_LINE;
    vtktypes[ELE_TAG_beam3d01] = VTK_LINE;
    vtktypes[ELE_TAG_beam3d02] = VTK_LINE;
    vtktypes[ELE_TAG_Truss] = VTK_LINE;
    vtktypes[ELE_TAG_TrussSection] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTruss] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTrussSection] = VTK_LINE;
    vtktypes[ELE_TAG_fElmt05] = VTK_LINE;
    vtktypes[ELE_TAG_fElmt02] = VTK_LINE;
    vtktypes[ELE_TAG_MyTruss] = VTK_LINE;
    vtktypes[ELE_TAG_ZeroLength] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthSection] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthND] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContactASDimplex] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthContactNTS2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthInterface2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_CoupledZeroLength] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ZeroLengthRocking] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_NLBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_NLBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_LargeDispBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_FourNodeQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_FourNodeQuad3d] = VTK_QUAD;
    vtktypes[ELE_TAG_Tri31] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_SixNodeTri] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_BeamWithHinges2d] = VTK_LINE;
    vtktypes[ELE_TAG_BeamWithHinges3d] = VTK_LINE;
    vtktypes[ELE_TAG_EightNodeBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentyNodeBrick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNodeBrick_u_p_U] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentyNodeBrick_u_p_U] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_FourNodeQuadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_TotalLagrangianFD20NodeBrick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_TotalLagrangianFD8NodeBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNode_LDBrick_u_p] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EightNode_Brick_u_p] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_TwentySevenNodeBrick] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BrickUP] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_Nine_Four_Node_QuadUP] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Twenty_Eight_Node_BrickUP] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Twenty_Node_Brick] = VTK_QUADRATIC_HEXAHEDRON;
    vtktypes[ELE_TAG_BBarFourNodeQuadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_BBarBrickUP] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_PlateMITC4] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellMITC4] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellMITC9] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ASDShellQ4] = VTK_QUAD;
    vtktypes[ELE_TAG_ASDShellT3] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_Plate1] = VTK_QUAD;
    vtktypes[ELE_TAG_Brick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_BbarBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_FLBrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_EnhancedQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_ConstantPressureVolumeQuad] = VTK_QUAD;
    vtktypes[ELE_TAG_NineNodeMixedQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_NineNodeQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_EightNodeQuad] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DispBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_TimoshenkoBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumnWarping3d] = VTK_LINE;
    vtktypes[ELE_TAG_HingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_HingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoPointHingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoPointHingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_OnePointHingedBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_OnePointHingedBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_BeamColumnJoint2d] = VTK_QUAD;
    vtktypes[ELE_TAG_BeamColumnJoint3d] = VTK_QUAD;
    vtktypes[ELE_TAG_ForceBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnWarping2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumn2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumnWarping2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticForceBeamColumn3d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnCBDI2d] = VTK_LINE;
    vtktypes[ELE_TAG_ForceBeamColumnCBDI3d] = VTK_LINE;
    vtktypes[ELE_TAG_DispBeamColumn2dInt] = VTK_LINE;
    vtktypes[ELE_TAG_InternalSpring] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SimpleJoint2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Joint2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Joint3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElastomericBearingPlasticity3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingPlasticity2d] = VTK_LINE;
    vtktypes[ELE_TAG_TwoNodeLink] = VTK_LINE;
    vtktypes[ELE_TAG_ActuatorCorot] = VTK_LINE;
    vtktypes[ELE_TAG_Actuator] = VTK_LINE;
    vtktypes[ELE_TAG_Adapter] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ElastomericBearingBoucWen2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingBoucWen3d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSliderSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSliderSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSlider2d] = VTK_LINE;
    vtktypes[ELE_TAG_FlatSlider3d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFPSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFPSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_SingleFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_DoubleFPSimple2d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFPSimple3d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFP2d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DoubleFP3d] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_TripleFPSimple2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFPSimple3d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_MultiFP2d] = VTK_LINE;
    vtktypes[ELE_TAG_MultiFP3d] = VTK_LINE;
    vtktypes[ELE_TAG_GenericClient] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_GenericCopy] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PY_MACRO2D] = VTK_LINE;
    vtktypes[ELE_TAG_SimpleContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SimpleContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SurfaceLoad] = VTK_QUAD;
    vtktypes[ELE_TAG_BeamContact2D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamEndContact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SSPquad] = VTK_QUAD;
    vtktypes[ELE_TAG_SSPquadUP] = VTK_QUAD;
    vtktypes[ELE_TAG_SSPbrick] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_SSPbrickUP] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_BeamContact2Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamContact3Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_BeamEndContact3Dp] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_Quad4FiberOverlay] = VTK_QUAD;
    vtktypes[ELE_TAG_Brick8FiberOverlay] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_QuadBeamEmbedContact] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_DispBeamColumn2dThermal] = VTK_LINE;
    vtktypes[ELE_TAG_TPB1D] = VTK_LINE;
    vtktypes[ELE_TAG_TFP_Bearing] = VTK_LINE;
    vtktypes[ELE_TAG_TFP_Bearing2d] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFrictionPendulum] = VTK_LINE;
    vtktypes[ELE_TAG_TripleFrictionPendulumX] = VTK_LINE;
    vtktypes[ELE_TAG_PFEMElement2D] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_FourNodeQuad02] = VTK_QUAD;
    vtktypes[ELE_TAG_cont2d01] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_cont2d02] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_CST] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_Truss2] = VTK_LINE;
    vtktypes[ELE_TAG_CorotTruss2] = VTK_LINE;
    vtktypes[ELE_Tag_ZeroLengthImpact3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PFEMElement3D] = VTK_TETRA;
    vtktypes[ELE_TAG_PFEMElement2DCompressible] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2DBubble] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2Dmini] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ElasticTimoshenkoBeam2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElasticTimoshenkoBeam3d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingUFRP2d] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericBearingUFRP3d] = VTK_LINE;
    vtktypes[ELE_TAG_RJWatsonEQS2d] = VTK_LINE;
    vtktypes[ELE_TAG_RJWatsonEQS3d] = VTK_LINE;
    vtktypes[ELE_TAG_HDR] = VTK_LINE;
    vtktypes[ELE_TAG_ElastomericX] = VTK_LINE;
    vtktypes[ELE_TAG_LeadRubberX] = VTK_LINE;
    vtktypes[ELE_TAG_PileToe3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_N4BiaxialTruss] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_ShellDKGQ] = VTK_QUAD;
    vtktypes[ELE_TAG_ShellNLDKGQ] = VTK_QUAD;
    vtktypes[ELE_TAG_MultipleShearSpring] = VTK_LINE;
    vtktypes[ELE_TAG_MultipleNormalSpring] = VTK_LINE;
    vtktypes[ELE_TAG_KikuchiBearing] = VTK_LINE;
    vtktypes[ELE_TAG_ComponentElement2d] = VTK_LINE;
    vtktypes[ELE_TAG_YamamotoBiaxialHDR] = VTK_LINE;
    vtktypes[ELE_TAG_MVLEM] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SFI_MVLEM] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_MVLEM_3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_SFI_MVLEM_3D] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_E_SFI_MVLEM_3D] = VTK_POLY_VERTEX;
	vtktypes[ELE_TAG_E_SFI] = VTK_POLY_VERTEX;
	vtktypes[ELE_TAG_MEFI] = VTK_POLY_VERTEX;
    vtktypes[ELE_TAG_PFEMElement2DFIC] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_TaylorHood2D] = VTK_QUADRATIC_TRIANGLE;
    vtktypes[ELE_TAG_PFEMElement2DQuasi] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_MINI] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_CatenaryCable] = VTK_LINE;
    vtktypes[ELE_TAG_FourNodeTetrahedron] = VTK_TETRA;
    vtktypes[ELE_TAG_PFEMElement3DBubble] = VTK_TETRA;
    vtktypes[ELE_TAG_TriSurfaceLoad] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ShellDKGT] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_ShellNLDKGT] = VTK_TRIANGLE;
    vtktypes[ELE_TAG_InertiaTruss] = VTK_LINE;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary2D] = VTK_QUAD;
    vtktypes[ELE_TAG_ASDAbsorbingBoundary3D] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_FSIFluidElement2D] = VTK_QUAD;
    vtktypes[ELE_TAG_FSIInterfaceElement2D] = VTK_LINE;
    vtktypes[ELE_TAG_FSIFluidBoundaryElement2D] = VTK_LINE;
    vtktypes[ELE_TAG_PML2D_3] = VTK_QUAD;
    vtktypes[ELE_TAG_PML2D_5] = VTK_QUAD;
    vtktypes[ELE_TAG_PML2D_12] = VTK_QUAD;
    vtktypes[ELE_TAG_PML2DVISCOUS] = VTK_QUAD;
    vtktypes[ELE_TAG_PML2D] = VTK_QUAD;
    vtktypes[ELE_TAG_PML3D] = VTK_HEXAHEDRON;
    vtktypes[ELE_TAG_PML3DVISCOUS] = VTK_HEXAHEDRON;
}



/**
 * Creates an HDF5 dataset with automatic dimension detection based on columns
 * 
 * param group_id Parent group ID where dataset will be created
 * param dataset_name Name of the dataset to create
 * param ncols Number of columns (0 for 1D dataset, >0 for 2D dataset)
 * param chunk_rows Number of rows per chunk (default 100)
 * param type_id HDF5 datatype for the dataset (default H5T_NATIVE_DOUBLE)
 * param max_rows Maximum number of rows (H5S_UNLIMITED for unlimited, 0 for fixed size)
 * return hid_t Dataset ID if successful, negative value if failed
 */
int VTKHDF_Recorder::createDataset(
    hid_t group_id,
    const char* dataset_name,
    hsize_t ncols,
    hsize_t chunk_rows,
    hid_t type_id,
    hsize_t max_rows
) {
   
    // Determine if it's 1D or 2D based on ncols
    int ndims = (ncols == 0) ? 1 : 2;
    
    // Declare fixed-size arrays
    hsize_t dims[2];
    hsize_t max_dims[2];
    hsize_t chunk_dims[2];

    // Set dimensions based on whether it's 1D or 2D
    if (ndims == 1) {
        dims[0] = 0;
        max_dims[0] = max_rows;
        chunk_dims[0] = chunk_rows;
    } else {
        dims[0] = 0;
        dims[1] = ncols;
        max_dims[0] = max_rows;
        max_dims[1] = ncols;
        chunk_dims[0] = chunk_rows;
        chunk_dims[1] = ncols;
    }

    // Create dataspace
    hid_t space_id = H5Screate_simple(ndims, dims, max_dims);
    if (space_id < 0) {
        fprintf(stderr, "Error creating dataspace for %s\n", dataset_name);
        return -1;
    }

    // Create property list for chunked storage
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    if (plist_id < 0) {
        fprintf(stderr, "Error creating property list for %s\n", dataset_name);
        H5Sclose(space_id);
        return -1;
    }

    // Set chunk size
    if (H5Pset_chunk(plist_id, ndims, chunk_dims) < 0) {
        fprintf(stderr, "Error setting chunk size for %s\n", dataset_name);
        H5Pclose(plist_id);
        H5Sclose(space_id);
        return -1;
    }

    // Create the dataset
    hid_t dset_id = H5Dcreate(
        group_id,
        dataset_name,
        type_id,
        space_id,
        H5P_DEFAULT,
        plist_id,
        H5P_DEFAULT
    );

    // Check if dataset was created successfully
    if (dset_id < 0) {
        fprintf(stderr, "Error creating dataset %s\n", dataset_name);
        H5Pclose(plist_id);
        H5Sclose(space_id);
        return -1;
    }

    // Clean up
    H5Pclose(plist_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);

    return 0;
}



int VTKHDF_Recorder::extendDataset(hid_t dataGroup,
                     const char* datasetName, 
                     const void* newData,
                     hid_t datatype,
                     size_t numNewElements,
                     size_t numColumns) {
    // Open dataset
    hid_t dset_id = H5Dopen(dataGroup, datasetName, H5P_DEFAULT);
    if (dset_id < 0) {
        opserr << "Error opening dataset '" << datasetName << "'\n";
        return -1;
    }

    // Get the current dataspace
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        opserr << "Error getting dataspace for '" << datasetName << "'\n";
        H5Dclose(dset_id);
        return -1;
    }

    // Get rank and dimensions
    int rank = H5Sget_simple_extent_ndims(space_id);
    std::vector<hsize_t> current_dims(rank);
    if (H5Sget_simple_extent_dims(space_id, current_dims.data(), NULL) < 0) {
        opserr << "Error getting current dimensions\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Prepare new dimensions
    std::vector<hsize_t> new_dims = current_dims;
    new_dims[0] += static_cast<hsize_t>(numNewElements);  // Extend first dimension

    // Extend the dataset
    if (H5Dset_extent(dset_id, new_dims.data()) < 0) {
        opserr << "Error setting new extent for '" << datasetName << "'\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Get the new dataspace
    hid_t new_space_id = H5Dget_space(dset_id);
    if (new_space_id < 0) {
        opserr << "Error getting new dataspace for '" << datasetName << "'\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Prepare hyperslab selection
    std::vector<hsize_t> offset(rank, 0);
    std::vector<hsize_t> count(rank);
    
    offset[0] = current_dims[0];  // Start from end of existing data
    count[0] = numNewElements;    // Number of new elements
    
    if (rank == 2) {
        count[1] = numColumns;    // Use specified number of columns for 2D
    }

    // Select hyperslab for writing
    if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset.data(), 
                            NULL, count.data(), NULL) < 0) {
        opserr << "Error selecting hyperslab\n";
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Create memory space
    hid_t mem_space_id = H5Screate_simple(rank, count.data(), NULL);
    if (mem_space_id < 0) {
        opserr << "Error creating memory space for '" << datasetName << "'\n";
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Write data with the specified datatype
    if (H5Dwrite(dset_id, datatype, mem_space_id, new_space_id, 
                    H5P_DEFAULT, newData) < 0) {
        opserr << "Error writing data to '" << datasetName << "'\n";
    }

    // Clean up
    H5Sclose(mem_space_id);
    H5Sclose(new_space_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);

    return 0;
}




int VTKHDF_Recorder::createOffsetDataset(hid_t group_id, 
                                        const char* dataset_name, 
                                        hsize_t chunk_rows,
                                        hid_t type_id,
                                        hsize_t max_rows) {
    // Define initial dimensions as 0 since we don't know the number of nodes yet
    hsize_t dims[1] = {0}; // Initial size (0 rows)
    hsize_t max_dims[1] = {max_rows}; // Maximum size (unlimited rows)

    // Create a dataspace with these dimensions
    hid_t space_id = H5Screate_simple(1, dims, max_dims);
    if (space_id < 0) {
        fprintf(stderr, "Error creating dataspace for %s\n", dataset_name);
        return -1;
    }

    // Create a property list for chunked storage
    hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
    if (plist_id < 0) {
        fprintf(stderr, "Error creating property list for %s\n", dataset_name);
        H5Sclose(space_id);
        return -1;
    }

    // Define the chunk size (e.g., 100 rows per chunk)
    hsize_t chunk_dims[1] = {chunk_rows};
    if (H5Pset_chunk(plist_id, 1, chunk_dims) < 0) {
            fprintf(stderr, "Error setting chunk size for %s\n", dataset_name);
            H5Pclose(plist_id);
            H5Sclose(space_id);
            return -1;
        }

    // Create the dataset with chunked storage
    hid_t dset_id = H5Dcreate(
        group_id,
        dataset_name,
        type_id,
        space_id,
        H5P_DEFAULT, 
        plist_id, 
        H5P_DEFAULT
        );

    if (dset_id < 0) {
        fprintf(stderr, "Error creating dataset for %s\n", dataset_name);
        H5Pclose(plist_id);
        H5Sclose(space_id);
        return -1;
    }

    // Clean up
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Sclose(space_id);// Create dataset with 1 column

    return 0;
}

int VTKHDF_Recorder::extendOffsetDataset(hid_t group_id, 
                  const char* dataset_name, 
                  const void* new_value, 
                  hid_t datatype,
                  size_t extra_rows) {
    // Open dataset
    hid_t dset_id = H5Dopen(group_id, dataset_name, H5P_DEFAULT);
    if (dset_id < 0) {
        opserr << "Error opening dataset '" << dataset_name << "'\n";
        return -1;
    }

    // Get current dataspace and dimensions
    hid_t space_id = H5Dget_space(dset_id);
    if (space_id < 0) {
        opserr << "Error getting dataspace for '" << dataset_name << "'\n";
        H5Dclose(dset_id);
        return -1;
    }

    // Get current dimensions
    hsize_t current_dims[1];
    if (H5Sget_simple_extent_dims(space_id, current_dims, NULL) < 0) {
        opserr << "Error getting current dimensions\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Calculate new dimensions
    hsize_t new_dims[1] = {current_dims[0] + static_cast<hsize_t>(extra_rows)};

    // Extend dataset
    if (H5Dset_extent(dset_id, new_dims) < 0) {
        opserr << "Error setting new extent for '" << dataset_name << "'\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Get new dataspace
    hid_t new_space_id = H5Dget_space(dset_id);
    if (new_space_id < 0) {
        opserr << "Error getting new dataspace for '" << dataset_name << "'\n";
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Select hyperslab for new data
    hsize_t offset[1] = {current_dims[0]};
    hsize_t count[1] = {1};
    if (H5Sselect_hyperslab(new_space_id, H5S_SELECT_SET, offset, NULL, count, NULL) < 0) {
        opserr << "Error selecting hyperslab\n";
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Create memory space for new data
    hid_t mem_space_id = H5Screate_simple(1, count, NULL);
    if (mem_space_id < 0) {
        opserr << "Error creating memory space for '" << dataset_name << "'\n";
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Write new data using the provided datatype
    if (H5Dwrite(dset_id, datatype, mem_space_id, new_space_id, H5P_DEFAULT, new_value) < 0) {
        opserr << "Error writing data to '" << dataset_name << "'\n";
        H5Sclose(mem_space_id);
        H5Sclose(new_space_id);
        H5Sclose(space_id);
        H5Dclose(dset_id);
        return -1;
    }

    // Clean up
    H5Sclose(mem_space_id);
    H5Sclose(new_space_id);
    H5Sclose(space_id);
    H5Dclose(dset_id);

    return 0;
}
