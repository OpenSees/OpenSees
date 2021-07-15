                                                                        
// Written: Minjie

// Description: command to create CrdTransf

#include <string.h>
#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>
#include <elementAPI.h>

void* OPS_LinearCrdTransf2d();
void* OPS_LinearCrdTransf3d();
void* OPS_PDeltaCrdTransf2d();
void* OPS_PDeltaCrdTransf3d();
void* OPS_CorotCrdTransf2d();
void* OPS_CorotCrdTransf3d();
void* OPS_CorotCrdTransfWarping3d();  

namespace {
    static void* OPS_LinearCrdTransf()
    {
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if (ndm == 2 && ndf == 3) {
	    return OPS_LinearCrdTransf2d();
	} else if (ndm == 3 && ndf == 6) {
	    return OPS_LinearCrdTransf3d();
	} else {
	    opserr<<"current NDM and NDF is incompatible with frame elements\n";
	    return 0;
	}
    }

    void* OPS_PDeltaCrdTransf()
    {
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if(ndm == 2 && ndf == 3) {
	    return OPS_PDeltaCrdTransf2d();
	} else if(ndm == 3 && ndf == 6) {
	    return OPS_PDeltaCrdTransf3d();
	} else {
	    opserr<<"current NDM and NDF is incompatible with frame elements\n";
	    return 0;
	}
    }
    void* OPS_CorotCrdTransf()
    {
	int ndm = OPS_GetNDM();
	int ndf = OPS_GetNDF();
	if(ndm == 2 && ndf == 3) {
	    return OPS_CorotCrdTransf2d();
	} else if(ndm == 3 && ndf == 6) {
	    return OPS_CorotCrdTransf3d();
	} else if(ndm == 3 && ndf == 7) {
	    return OPS_CorotCrdTransfWarping3d();
	} else {	  
	    opserr<<"current NDM and NDF is incompatible with frame elements\n";
	    return 0;
	}
    }
    
}

int OPS_CrdTransf()
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"WARNING too few arguments: geomTransf type? tag? ...\n";
	return -1;
    }

    const char* type = OPS_GetString();

    // create it
    CrdTransf* transf = 0;
    if (strcmp(type, "Linear") == 0) {
    	transf = (CrdTransf*)OPS_LinearCrdTransf();
    } else if(strcmp(type, "PDelta") == 0) {
    	transf = (CrdTransf*)OPS_PDeltaCrdTransf();
    } else if(strcmp(type, "Corotational") == 0) {
    	transf = (CrdTransf*)OPS_CorotCrdTransf();
    } else {
    	opserr << "WARNING geomTransf type "<< type <<" is unknown\n";
    	return -1;
    }

    // check
    if (transf == 0) {
	opserr << "WARNING failed to create geomTransf object\n";
	return -1;
    }

    // add it
    if (OPS_addCrdTransf(transf) == false) {
	opserr << "WARNING failed to add CrdTransf\n";
	delete transf;
	return -1;
    }

    return 0;
}
