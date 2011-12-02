/* 
 * @author: gnp <petropoulos@gmail.com>
 *
 * @Description: DRM Load Pattern Wrapper Header File
 *               purpose is to ease the implementation
 *               of the send/recv self methods
 *
 * @Date: 7/7/08
 *
 */


#ifndef DRMLoadPatternWrapper_h
#define DRMLoadPatternWrapper_h

#include <LoadPattern.h>
#include <Domain.h>
#include <LoadPattern.h>
#include "DRMBoundaryLayerDecorator.h"
#include "DRMLoadPattern.h"
#include "DRMInputHandler.h"
#include <Vector.h>
#include <ID.h>
#include <Message.h>
#include <Matrix.h>
#include <Element.h>
#include <map>
#include <set>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string>
#include <classTags.h>

class DRMLoadPatternWrapper : public LoadPattern
{
  public :
    
    DRMLoadPatternWrapper();
  
    DRMLoadPatternWrapper(int _tag, double _cfact, char** _in_files, int _files,  double _dt,
			  int _num_steps, int* _file_data, int _fileData_size,
			  int __nd1, int __nd2,
			  double* _drm_box_crds, double* _eleD,
			  int _steps_to_cache);

    virtual ~DRMLoadPatternWrapper();

    void initialize();

    void applyLoad(double time);

    int sendSelf(int commitTag, Channel& theChannel); 

    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker); 

    void Print(OPS_Stream &s, int flag);

    LoadPattern* getCopy();
    
 private:

    char** infiles;
    int files;
    double dt;
    int num_steps;
    int* file_data;
    int fileData_size;
    int nd1;
    int nd2;
    double* drm_box_crds;
    double* eleD;
    int steps_to_cache;

    DRMLoadPattern* myPattern;
    
    double factor;
    
    bool initialized;
    bool cleanUpAfterMySelf;
};
#endif

