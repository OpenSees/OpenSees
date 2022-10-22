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
                                                                        
// Original implementation: José Abell (UANDES), Massimo Petracca (ASDEA)
//
// ASDPlasticMaterial
//
// Fully general templated material class for plasticity modeling


#include "ArmstrongFrederickTensor_EV.h"
#include <Vector.h>

VoigtVector ArmstrongFrederickTensor_EV::derivative(3, 3, 0.0);

ArmstrongFrederickTensor_EV::ArmstrongFrederickTensor_EV( double ha_, double cr_) : EvolvingVariable(VoigtVector(3, 3, 0.0)), ha(ha_), cr(cr_)
{}

ArmstrongFrederickTensor_EV::ArmstrongFrederickTensor_EV( double ha_, double cr_, VoigtVector& alpha0) : EvolvingVariable(alpha0), ha(ha_), cr(cr_)
{}

const VoigtVector& ArmstrongFrederickTensor_EV::getDerivative(const VoigtVector &depsilon,
        const VoigtVector &m,
        const VoigtVector& stress) const
{
    using namespace ASDPlasticMaterialGlobals;
    //Zero de static variable
    derivative *= 0;



    //Compute the derivative (hardening function)
    const VoigtVector &alpha = this->getVariableConstReference();
    static VoigtVector mdev(3, 3, 0);
    mdev *= 0;
    mdev(i, j) = m(i, j) - m(k, k) / 3 * kronecker_delta(i, j);
    derivative(i, j) =  (2. / 3.) * ha * m(i, j) - cr * sqrt((2. / 3.) * m(k, l) * m(k, l)) * alpha(i, j);




    return derivative;
}

void ArmstrongFrederickTensor_EV::check_hardening_saturation_limit(VoigtVector& backstress, VoigtVector const& plasticFlow_m ){
    using namespace ASDPlasticMaterialGlobals;
    double limit_length = SQRT_2_over_3 * ha / cr ;

    // Limit direction is unit vector in the plastic flow direction.
    static VoigtVector limit_direction(3,3,0.0);
    limit_direction(i,j) = plasticFlow_m(i,j) / sqrt (plasticFlow_m(k,l) * plasticFlow_m(k,l)) ;

    // Limit each component
    for (int i_ = 0; i_ < 3; ++i_)
        for (int j_ = 0; j_ < 3; ++j_){
            if(limit_direction(i_,j_) < 0) limit_direction(i_,j_) = - limit_direction(i_,j_);
            if(backstress(i_,j_) > +limit_length*limit_direction(i_,j_) ) backstress(i_,j_) = +limit_length*limit_direction(i_,j_);
            if(backstress(i_,j_) < -limit_length*limit_direction(i_,j_) ) backstress(i_,j_) = -limit_length*limit_direction(i_,j_);
        }
    // // // ====================================================
    // // // debug purpose printing
    // // // ====================================================
    // cout<<"backstress(0,1)"<<backstress(0,1)<<endl;
    // cout<<"new"<<endl;
    // cout<<backstress<<endl;
    // cout<<"---------------------------------\n";
    // // // ====================================================
}


int ArmstrongFrederickTensor_EV::sendSelf(int commitTag, Channel &theChannel)
{
    //Shove all data into single vector for sending
    static Vector data(9 + 9 + 2);
    const VoigtVector &a = this->getVariableConstReference();
    const VoigtVector &a_committed = this->getVariableConstReference();
    int pos = 0;

    data(pos++) = ha;
    data(pos++) = cr;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            data(pos++) = a(i, j);
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            data(pos++) = a_committed(i, j);
        }
    if (theChannel.sendVector(0, commitTag, data) != 0)
    {
        cerr << "ArmstrongFrederickTensor_EV::sendSelf() - Failed sending data" << endl;
        return -1;
    }

    return 0;
}

int ArmstrongFrederickTensor_EV::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static Vector data(9 + 9 + 2);
    if (theChannel.receiveVector(0, commitTag, data) != 0)
    {
        cerr << "ArmstrongFrederickTensor_EV::recvSelf() - Failed recieving data" << endl;
        return -1;
    }

    //Extract data from vector
    int pos = 0;
    static VoigtVector tmp_a(3, 3, 0);
    static VoigtVector tmp_a_committed(3, 3, 0);
    ha = data(pos++);
    cr = data(pos++);

    // TODO: This loop should get taken care of by the variables themselves.
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp_a(i, j) = data(pos++);
        }
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {
            tmp_a_committed(i, j) = data(pos++) ;
        }

    this->setVar(tmp_a);
    this->setCommittedVar(tmp_a_committed);

    return 0;
}