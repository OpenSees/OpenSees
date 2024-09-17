
// ============================================================================
// 2021 By Jose Abell @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Implements simple routines to test nDMaterials using direct strains.
//
// ============================================================================

#include <TimeSeries.h>
#include <map>
#include <string.h>
#include <elementAPI.h>
#include <vector>
#include <PathTimeSeries.h>
#include <PathSeries.h>
#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>
#include <Information.h>
#include <Response.h>
#include <MaterialResponse.h>

//Forward declaration for indexing purposes
Tcl_CmdProc OPS_NDSetStrain;               // To set a trial strain on the material
Tcl_CmdProc OPS_NDPrintStrain;             // Printing strains to the screen
Tcl_CmdProc OPS_NDPrintStress;             // Printing strains to the screen
Tcl_CmdProc OPS_NDGetStrain;               // Return current strains to a list variable
Tcl_CmdProc OPS_NDGetStress;               // Return current stresses to a list variable
Tcl_CmdProc OPS_NDGetTangentStiffness;     // Return current tangent stiffness to a list variable
Tcl_CmdProc OPS_NDCommitState;             // To commit the state of the material after last trial strain
Tcl_CmdProc OPS_NDUpdateIntegerParameter;  // To change an integer parameter in the material
Tcl_CmdProc OPS_NDUpdateDoubleParameter;   // To change a double parameter in the material
Tcl_CmdProc OPS_NDGetResponse;             // To get other responses

int OPS_NDSetStrain(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    double strains[6];

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        opserr << "OPS_NDSetStrain - got incorrect integer tag for material" << "\n";
        return 0;
    }


    numdata = 6;
    if (OPS_GetDoubleInput(&numdata, strains) < 0) {
        opserr << "OPS_NDSetStrain - need 6 components of floating-point strains" << "\n";
        return 0;
    }

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDSetStrain - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }


    Vector new_strain(6);

    new_strain(0) = strains[0];
    new_strain(1) = strains[1];
    new_strain(2) = strains[2];
    new_strain(3) = strains[3];
    new_strain(4) = strains[4];
    new_strain(5) = strains[5];

    mat->setTrialStrain(new_strain);

    return 0;
}

int OPS_NDPrintStrain(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDPrintStrain - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    const Vector &strain = mat->getStrain();

    return 0;
}

int OPS_NDPrintStress(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDPrintStress - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    const Vector &stress = mat->getStress();

    return 0;
}


int OPS_NDGetStrain(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    int size = 6;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDGetStrain - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    const Vector &strain = mat->getStrain();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        values[i] = strain(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



int OPS_NDGetStress(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    int size = 6;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        return 0;
    }


    NDMaterial* mat = OPS_getNDMaterial(tag);
    if (mat == nullptr)
    {
        opserr << "OPS_NDGetStress() - Material tag " << tag << " not declared" << "\n";
    }

    const Vector &stress = mat->getStress();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        values[i] = stress(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}

int OPS_NDGetResponse(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0)
    {
        return 0;
    }


    NDMaterial* mat = OPS_getNDMaterial(tag);
    if (mat == nullptr)
    {
        opserr << "OPS_NDGetStress() - Material tag " << tag << " not declared" << "\n";
    }


    const char* response_type = OPS_GetString();

    Response* theResponse = mat->setResponse(&response_type, 1, opserr);
    theResponse->getResponse();
    Information& theInformation = theResponse->getInformation();

    const Vector &respose = *theInformation.theVector;

    int size = respose.Size();
    std::vector<double> values(size);
    for (int i = 0; i < size; i++) {
        values[i] = respose(i);
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



int OPS_NDGetTangentStiffness(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    int size = 36;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;

    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDGetTangentStiffness - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    const Matrix &stiffness = mat->getTangent();

    std::vector<double> values(size);
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++)
        {
            double one_value = stiffness(i, j);
            values[6 * i + j] = one_value;
        }
    }
    if (OPS_SetDoubleOutput(&size, &values[0], false) < 0) {
        opserr << "WARNING OPS_NDGetStress - failed to set double inputs\n";
        return 0;
    }


    return 0;
}



int OPS_NDCommitState(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    double strains[6];
    int size = 6;
    double stressdata[6] = {1 , 2 , 3, 4, 5, 6};

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;


    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDCommitState - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }


    mat->commitState();

    return 0;
}


int OPS_NDUpdateIntegerParameter(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{

    opserr << "OPS_NDUpdateIntegerParameter" << "\n";

    int tag = 0;
    int responseID = 0;
    int theNewIntegerParameterValue = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
    opserr << "tag = " << tag  << "\n";
    if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
    opserr << "responseID = " << responseID  << "\n";
    if (OPS_GetIntInput(&numdata, &theNewIntegerParameterValue) < 0) return 0;
    opserr << "theNewIntegerParameterValue = " << theNewIntegerParameterValue  << "\n";


    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_getNDMaterial - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    Information info;

    info.theInt = theNewIntegerParameterValue;

    mat->updateParameter(responseID, info);

    return 0;
}


int OPS_NDUpdateDoubleParameter(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{
    int tag = 0;
    int responseID = 0;
    double theNewDoubleParameterValue = 0;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) return 0;
    if (OPS_GetIntInput(&numdata, &responseID) < 0) return 0;
    if (OPS_GetDoubleInput(&numdata, &theNewDoubleParameterValue) < 0) return 0;


    NDMaterial* mat = OPS_getNDMaterial(tag);

    if (mat == nullptr)
    {
        opserr << "OPS_NDUpdateDoubleParameter - material with tag " << tag << " does not exist" << "\n";
        return 0;
    }

    Information info;

    info.theDouble = theNewDoubleParameterValue;

    mat->updateParameter(responseID, info);

    return 0;
}

typedef std::unordered_map<std::string, Tcl_CmdProc> invoke_material_dispatch {
    {"SetStrain",              &OPS_NDSetStrain},
    {"CommitState",            &OPS_NDCommitState},
    {"PrintStress",            &OPS_NDPrintStress},
    {"PrintStrain",            &OPS_NDPrintStrain},
    {"GetStrain",              &OPS_NDGetStrain},
    {"GetStress",              &OPS_NDGetStress},
    {"GetTangentStiffness",    &OPS_NDGetTangentStiffness},
    {"UpdateIntegerParameter", &OPS_NDUpdateIntegerParameter},
    {"UpdateDoubleParameter",  &OPS_NDUpdateDoubleParameter},
    {"GetResponse",            &OPS_NDGetResponse}
};


int
OPS_NDTest(ClientData clientData, Tcl_Interp *interp, int argc, const char** const argv)
{

    // Identify what specific command of Patch we're calling
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING too few arguments: NDTest cmd? \n";
        opserr << " available commands: SetStrain|CommitState|GetStrain|GetStress \n";
        return -1;
    }

    const char* type = OPS_GetString();

    OPS_ParsingFunctionMap::const_iterator iter = functionMap.find(type);
    if (iter == functionMap.end()) {
        opserr << "WARNING NDTest type " << type << " is unknown\n";
        return -1;
    }

    // Call the function
    (*iter->second)();

    return 0;

}
