
// Description: command to create damping

#include <math.h>
#include <Vector.h>
#include <Matrix.h>
#include <string.h>
#include <elementAPI.h>

#include <UniformDamping.h>
#include <SecStifDamping.h>
#include <URDDamping.h>
#include <URDDampingbeta.h>


void* OPS_UniformDamping()
{
    
    Damping* damping = 0;
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 3) {
        opserr << "insufficient arguments: uniformdanping\n";
        return 0;
    }


    int dampingTag = 0;
    double dData[6];
    dData[3] = 0.0;
    dData[4] = 1e20;
    TimeSeries* facSeries = 0;
    int numData = 0;

    //get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &dampingTag) != 0) {
            opserr << "WARNING invalid dampingtag in uniform tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    numData = 3;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in uniformdamping with tag: " << dampingTag << endln;
        return 0;
    }
    numRemainingArgs -= 3;


    //options
    while (numRemainingArgs > 1) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-activateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[3]) != 0) {
                opserr << "WARNING invalid gamma in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-deactivateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[4]) != 0) {
                opserr << "WARNING invalid nu in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-factor") == 0) {
            numData = 1;
            int tsTag;
            if (OPS_GetIntInput(&numData, &tsTag) != 0) {
                opserr << "WARNING invalid amplitude in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
            facSeries = OPS_getTimeSeries(tsTag);
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in uniformdamping with tag?" << dampingTag << endln;
            return 0;
        }
        numRemainingArgs -= 2;
    }
    damping = new UniformDamping(dampingTag, dData[0]*2.0, dData[1], dData[2], dData[3], dData[4], facSeries);
    return damping;
}

void* OPS_SecStifDamping()
{

    Damping* damping = 0;
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 1) {
        opserr << "insufficient arguments: uniformdanping\n";
        return 0;
    }


    int dampingTag = 0;
    double dData[4];
    dData[1] = 0.0;
    dData[2] = 1e20;
    TimeSeries* facSeries = 0;
    int numData = 0;

    //get tag if provided
    if (numRemainingArgs == 2 || numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &dampingTag) != 0) {
            opserr << "WARNING invalid dampingtag in uniform tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }
    opserr << "dampingTag = " << dampingTag << endln;//test

    numData = 1;
    if (OPS_GetDouble(&numData, dData) != 0) {
        opserr << "WARNING invalid double data in SecStifdamping with tag: " << dampingTag << endln;
        return 0;
    }
    numRemainingArgs -= 1;

    opserr << "numRemainingArgs = " << numRemainingArgs << endln;//test

    //options
    while (numRemainingArgs > 1) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-activateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[1]) != 0) {
                opserr << "WARNING invalid gamma in SecStifdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-deactivateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[2]) != 0) {
                opserr << "WARNING invalid nu in SecStifdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-factor") == 0) {
            numData = 1;
            int tsTag;
            if (OPS_GetIntInput(&numData, &tsTag) != 0) {
                opserr << "WARNING invalid amplitude in SecStifdamping with tag?" << dampingTag << endln;
                return 0;
            }
            facSeries = OPS_getTimeSeries(tsTag);
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in SecStifdamping with tag?" << dampingTag << endln;
            return 0;
        }
        numRemainingArgs -= 2;
    }
    damping = new SecStifDamping(dampingTag, dData[0], dData[1], dData[2], facSeries);
    return damping;

}

void* OPS_URDDamping()
{

    Damping* damping = 0;
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 2) {
        opserr << "insufficient arguments: URDdamping\n";
        return 0;
    }


    int dampingTag = 0;
    int numfreq;
    int iData[3];     //numfreq prttag maxiter
    double dData[3];  //dptol ta td
    double tmpetafeq;
    dData[0] = 0.05;
    dData[1] = 0.0;
    dData[2] = 1e20;
    TimeSeries* facSeries = 0;
    int numData = 0;

    //get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10 || numRemainingArgs == 12 || numRemainingArgs == 14 || numRemainingArgs == 16 || numRemainingArgs == 18 || numRemainingArgs == 20 || numRemainingArgs == 22 || numRemainingArgs == 24 || numRemainingArgs == 26 || numRemainingArgs == 28 || numRemainingArgs == 30) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &dampingTag) != 0) {
            opserr << "WARNING invalid dampingtag in URD tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }
    opserr << "dampingTag = " << dampingTag << endln;//test

    numData = 1;
    if (OPS_GetIntInput(&numData, &iData[0]) != 0) {
        opserr << "WARNING invalid gamma in URDdamping with tag?" << dampingTag << endln;
        return 0;
    }
    numRemainingArgs -= 1;

    numfreq = iData[0];

    opserr << "numfreq = " << numfreq << endln;//test
    
    if (numfreq < 2) {
        opserr << "WARNING - n needs to be larger than 1\n ";
        return  0;
    }
    
    // Write the obtained number into the matrix
    Matrix* etaFreq = new Matrix(numfreq, 2);
    for (int i = 0; i < numfreq; i++) {
        for (int j = 0; j < 2; j++) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &tmpetafeq) != 0) {
                opserr << "WARNING invalid factor series - want: damping URD tag n freq1 zeta1 ... freqn zetan\n ";
                return 0;
            }
            numRemainingArgs -= 1;
            if (tmpetafeq <= 0.0) opserr << "URDDamping::URDDamping:  Invalid frequency or damping ratio\n";
            (*etaFreq)(i, j) = tmpetafeq;
        }
        //(*etaFreq)(i, 1) *= (2.0);
    }

    opserr << "numRemainingArgs = " << numRemainingArgs << endln;//test

    //options
    while (numRemainingArgs > 1) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-activateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[1]) != 0) {
                opserr << "WARNING invalid gamma in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-deactivateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[2]) != 0) {
                opserr << "WARNING invalid nu in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-factor") == 0) {
            numData = 1;
            int tsTag;
            if (OPS_GetIntInput(&numData, &tsTag) != 0) {
                opserr << "WARNING invalid amplitude in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
            facSeries = OPS_getTimeSeries(tsTag);
        }
        else if (strcmp(argvS, "-tol") == 0 || strcmp(argvS, "-tolerence") == 0) {
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &dData[0]) != 0) {
                opserr << "WARNING invalid amplitude in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-maxiter") == 0 || strcmp(argvS, "-iter") == 0) {
            numData = 1;
            if (OPS_GetIntInput(&numData, &iData[2]) != 0) {
                opserr << "WARNING invalid amplitude in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-prttag") == 0 || strcmp(argvS, "-print") == 0) {
            iData[1] = 1;
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in uniformdamping with tag?" << dampingTag << endln;
            return 0;
        }
        numRemainingArgs -= 2;
    }
    damping = new URDDamping(dampingTag, iData[0], etaFreq, dData[0], dData[1], dData[2], facSeries, iData[1], iData[2]);
    return damping;

}

void* OPS_URDbetaDamping()
{

    Damping* damping = 0;
    int numRemainingArgs = OPS_GetNumRemainingInputArgs();

    if (numRemainingArgs < 1) {
        opserr << "insufficient arguments: uniformdanping\n";
        return 0;
    }

    int dampingTag = 0;
    int iData[1];
    int numfreq;
    double dData[3];
    dData[0] = 0.0;
    dData[1] = 1e20;
    TimeSeries* facSeries = 0;
    int numData = 0;

    //get tag if provided
    if (numRemainingArgs == 4 || numRemainingArgs == 6 || numRemainingArgs == 8 || numRemainingArgs == 10 || numRemainingArgs == 12 || numRemainingArgs == 14 || numRemainingArgs == 16 || numRemainingArgs == 18 || numRemainingArgs == 20 || numRemainingArgs == 22 || numRemainingArgs == 24 || numRemainingArgs == 26 || numRemainingArgs == 28 || numRemainingArgs == 30) {
        numData = 1;
        if (OPS_GetIntInput(&numData, &dampingTag) != 0) {
            opserr << "WARNING invalid URDbetadampingtag in uniform tag?" << endln;
            return 0;
        }
        numRemainingArgs -= 1;
    }

    opserr << "dampingTag = " << dampingTag << endln;//test

    numData = 1;
    if (OPS_GetIntInput(&numData, iData) != 0) {
        opserr << "WARNING invalid gamma in URDbetadamping with tag?" << dampingTag << endln;
        return 0;
    }
    numRemainingArgs -= 1;

    numfreq = iData[0];

    if (numfreq < 2) {
        opserr << "WARNING - n needs to be larger than 1\n ";
        return  0;
    }

    opserr << "numfreq = " << numfreq << endln;//test

    Vector* tmpbeta = new Vector(numfreq);
    Vector* tmpomegac = new Vector(numfreq);
    for (int i = 0; i < numfreq; i++) {
        numData = 1;
        if (OPS_GetDoubleInput(&numData, &(*tmpomegac)(i)) != 0) {
            opserr << "WARNING invalid factor series - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan\n ";
            return 0;
        }
        numRemainingArgs -= 1;
        if ((*tmpomegac)(i) <= 0.0) opserr << "URDDamping::URDDamping:  Invalid frequency\n";
        (*tmpomegac)(i) *= (6.28318530718);
        if (OPS_GetDoubleInput(&numData, &(*tmpbeta)(i)) != 0) {
            opserr << "WARNING invalid factor series - want: damping URDbeta tag n freqc1 beta1 ... freqcn betan\n ";
            return 0;
        }
        numRemainingArgs -= 1;
    }
    opserr << "numRemainingArgs = " << numRemainingArgs << endln;//test
    //options
    while (numRemainingArgs > 1) {
        const char* argvS = OPS_GetString();

        if (strcmp(argvS, "-activateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[0]) != 0) {
                opserr << "WARNING invalid gamma in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-deactivateTime") == 0) {
            numData = 1;
            if (OPS_GetDouble(&numData, &dData[1]) != 0) {
                opserr << "WARNING invalid nu in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
        }
        else if (strcmp(argvS, "-factor") == 0) {
            numData = 1;
            int tsTag;
            if (OPS_GetIntInput(&numData, &tsTag) != 0) {
                opserr << "WARNING invalid amplitude in uniformdamping with tag?" << dampingTag << endln;
                return 0;
            }
            facSeries = OPS_getTimeSeries(tsTag);
        }
        else {
            opserr << "WARNING unknown option: " << argvS << "  in uniformdamping with tag?" << dampingTag << endln;
            return 0;
        }
        numRemainingArgs -= 2;
    }
    damping = new URDDampingbeta(dampingTag, iData[0], tmpomegac, tmpbeta, dData[0], dData[1], facSeries);
    return damping;

}

int OPS_Damping()    
{
    if (OPS_GetNumRemainingInputArgs() < 2) {
        opserr << "WARNING too few arguments: damping type? tag? ...\n";
        return -1;
    }

    const char* type = OPS_GetString();


    //create it
    Damping* damping = 0;
    if (strcmp(type, "uniform") == 0) {
        damping = (Damping*)OPS_UniformDamping();
    }
    else if (strcmp(type, "SecStif") == 0) {
        damping = (Damping*)OPS_SecStifDamping();
    }
    else if (strcmp(type, "URD") == 0) {
        damping = (Damping*)OPS_URDDamping();
    }
    else if (strcmp(type, "URDbeta") == 0) {
        damping = (Damping*)OPS_URDbetaDamping();
    }
    else {
        opserr << "warning damping type " << type << " is unknown\n";
        return -1;
    }

    // check
    if (damping == 0) {
	opserr << "warning failed to create damping object\n";
	return -1;
    }

    // add the damping to the modelbuilder
    if (OPS_addDamping(damping) != true) {
	opserr << "warning  could not add damping to model builder\n";
	delete damping;
	return -1;
    }

    return 0;
}
