/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT' in the main directory for information on usage and **
** redistribution, and for a DISCLAIMER OF ALL WARRANTIES.            **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
#ifndef VTKHDF_RECORDER_H
#define VTKHDF_RECORDER_H

// Written: Amin Pakzad, 2025
// Purpose: Class definition for VTKHDF_Recorder to store responses in HDF5 format.

// Include existing headers plus HDF5

#include <string>
#include <fstream>
#include <vector>
#include <map>
#include <ID.h>
#include "Response.h"
#include <Recorder.h>
#include <vector>
#include <cstdint>
#include "hdf5.h"


// Define the OutputData class
class OutputDataHDF {
public:
    OutputDataHDF();
    OutputDataHDF &operator=(const OutputDataHDF &other);

    bool disp;
    bool vel;
    bool accel; 
    bool reaction; 
    bool mass, unbalancedLoad;
    bool stress3D6, strain3D6;
    bool stress2D3, strain2D3;
};


class SilentStream : public OPS_Stream {
public:
    SilentStream();
    ~SilentStream() override;

    // Implement all pure virtual methods with no-op behavior
    int tag(const char* tagName) override { return 0; }
    int tag(const char* tagName, const char* value) override { return 0; }
    int endTag() override { return 0; }
    int attr(const char* name, int value) override { return 0; }
    int attr(const char* name, double value) override { return 0; }
    int attr(const char* name, const char* value) override { return 0; }
    int write(Vector& data) override { return 0; }
    int flush() override { return 0; }
    int open() override { return 0; }
    int close(openMode nextOpen = APPEND) override { return 0; }

    // Implement the write methods with no-op behavior
    OPS_Stream& write(const char* s, int n) override { return *this; }
    OPS_Stream& write(const unsigned char* s, int n) override { return *this; }
    OPS_Stream& write(const signed char* s, int n) override { return *this; }
    OPS_Stream& write(const void* s, int n) override { return *this; }
    OPS_Stream& write(const double* s, int n) override { return *this; }

    // Implement output operators with no-op behavior
    OPS_Stream& operator<<(char c) override { return *this; }
    OPS_Stream& operator<<(unsigned char c) override { return *this; }
    OPS_Stream& operator<<(signed char c) override { return *this; }
    OPS_Stream& operator<<(const char* s) override { return *this; }
    OPS_Stream& operator<<(const unsigned char* s) override { return *this; }
    OPS_Stream& operator<<(const signed char* s) override { return *this; }
    OPS_Stream& operator<<(const void* p) override { return *this; }
    OPS_Stream& operator<<(int n) override { return *this; }
    OPS_Stream& operator<<(unsigned int n) override { return *this; }
    OPS_Stream& operator<<(long n) override { return *this; }
    OPS_Stream& operator<<(unsigned long n) override { return *this; }
    OPS_Stream& operator<<(short n) override { return *this; }
    OPS_Stream& operator<<(unsigned short n) override { return *this; }
    OPS_Stream& operator<<(bool b) override { return *this; }
    OPS_Stream& operator<<(double n) override { return *this; }
    OPS_Stream& operator<<(float n) override { return *this; }

    // Implement parallel methods with no-op behavior
    void setAddCommon(int flag) override {}
    int setOrder(const ID& order) override { return 0; }
    int sendSelf(int commitTag, Channel& theChannel) override { return 0; }
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) override { return 0; }
};

// Define the VTKHDF_Recorder class
class VTKHDF_Recorder : public Recorder {
public:
    typedef std::vector<std::string> EleData;

    VTKHDF_Recorder(const char *filename, const OutputDataHDF &ndata,
                    const std::vector<EleData> &edata, double dt = 0, double rTolDt = 0.00001);
    VTKHDF_Recorder();
    ~VTKHDF_Recorder();

    int record(int commitTag, double timeStamp) override;
    int restart() override;
    int domainChanged() override;
    int setDomain(Domain &domain) override;
    int sendSelf(int commitTag, Channel &theChannel) override;
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker) override;

private:
    int createDataset(hid_t group_id, const char* dataset_name, 
                          hsize_t ncols, hsize_t chunk_rows = 100,
                          hid_t type_id = H5T_NATIVE_DOUBLE, 
                          hsize_t max_rows = H5S_UNLIMITED);

    int createOffsetDataset(hid_t group_id, const char* dataset_name, 
                                hsize_t chunk_rows = 100,
                                hid_t type_id = H5T_NATIVE_DOUBLE,
                                hsize_t max_rows = H5S_UNLIMITED);

    int extendDataset(hid_t dataGroup,
                     const char* datasetName, 
                     const void* newData,
                     hid_t datatype,
                     size_t numNewElements,
                     size_t numColumns = 0);
    

    int extendOffsetDataset(hid_t group_id, 
                  const char* dataset_name, 
                  const void* new_value, 
                  hid_t datatype,
                  size_t extra_rows = 1);
    

    OutputDataHDF outputData;

protected:
    bool initDone;
    virtual void addEleData(const EleData &edata) { eledata.push_back(edata); }
    std::vector<EleData> eledata;

private:
    int writeMesh();
    int writeStep(double timeStamp);
    int writeDisp();
    int writeVel();
    int writeAccel();
    int writeStrees3D6();
    int writeStrain3D6();
    int writeStress2D3();
    int writeStrain2D3();



    Domain *theDomain;

    double nextTimeStampToRecord;
    double deltaT;
    double relDeltaTTol;

    char *name;
    int counter;

    int ndm;
    int ndf;
    
    int numNode;
    int numElement;
    int numConnectivityIds;
    int maxNDM;
    int maxNDF;

    hid_t group_id;
    hid_t file_id;

    bool initializationDone;

    hid_t dispOffsetDataset;
    hid_t dispDataset;
    

    int current_PointOffset;
    int current_CellOffset;
    int current_ConnectivityIdOffset;
    int current_PartOffset;
    int next_PointOffset;
    int next_CellOffset;
    int next_ConnectivityIdOffset;
    int next_PartOffset;



    // data specif offsets
    int CurrentDispOffset;
    int CurrentVelOffset;
    int CurrentAccelOffset;

    // 3Dstrees6
    int current_Stress3D6Offset;
    std::vector<Response *> stress3D6Responses; // specific to 3DStress6 responses

    // 3DStrain6
    int current_Strain3D6Offset;
    std::vector<Response *> strain3D6Responses; // specific to 3DStrain6 responses

    // 2Dstrees3
    int current_Stress2D3Offset;
    std::vector<Response *> stress2D3Responses; // specific to 2DStress3 responses

    // 2DStrain3
    int current_Strain2D3Offset;
    std::vector<Response *> strain2D3Responses; // specific to 2DStrain3 responses





    unsigned int numSteps;



    std::map<int, int> theNodeMapping;
    std::map<int, int> theEleMapping;

    std::vector<int> theNodeTags;
    std::vector<int> theEleTags;
    std::vector<int> theEleClassTags;
    std::vector<unsigned char> theEleVtkTags;
    std::vector<int> theEleVtkOffsets;

    // Create groups for temporal data
    hid_t steps_group;
    hid_t point_data_group;
    hid_t cell_data_group;
    hid_t point_data_offsets_group;
    hid_t cell_data_offsets_group;

    hid_t nsteps_attr_id;



    // Keep VtkType enum and related static members
    enum VtkType {
        VTK_VERTEX = 1,
        VTK_POLY_VERTEX = 2,
        VTK_LINE = 3,
        VTK_POLY_LINE = 4,
        VTK_TRIANGLE = 5,
        VTK_TRIANGLE_STRIP = 6,
        VTK_POLYGON = 7,
        VTK_PIXEL = 8,
        VTK_QUAD = 9,
        VTK_TETRA = 10,
        VTK_VOXEL = 11,
        VTK_HEXAHEDRON = 12,
        VTK_WEDGE = 12,
        VTK_PYRAMID = 14,
        VTK_QUADRATIC_EDGE = 21,
        VTK_QUADRATIC_TRIANGLE = 22,
        VTK_QUADRATIC_QUAD = 23,
        VTK_QUADRATIC_TETRA = 24,
        VTK_QUADRATIC_HEXAHEDRON = 25
    };

    static std::map<int, VtkType> vtktypes;
    static void setVTKType();
};

#endif
