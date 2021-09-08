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
                                                                        
// $Revision: 1.0 $
// $Date: 2020-11-27 18:07:11 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/domain/DomainModalProperties.h,v $
                                                                        
// Written: Massimo Petracca (ASDEA Software) 
// Created: Fri Nov 27 18:07:11: 2020
// Revision: A
//
// Description: This file contains the class definition for DomainModalProperties.
// It contains all data optionally computed after an Eigenvalue analysis.
//
// What: "@(#) DomainModalProperties.h, revA"

#ifndef DomainModalProperties_h
#define DomainModalProperties_h

#include <Vector.h>
#include <Matrix.h>
#include <string>

class Domain;

class DomainModalProperties
{
public:
    DomainModalProperties(bool unorm = false);
    DomainModalProperties(const DomainModalProperties&) = default;
    DomainModalProperties& operator = (const DomainModalProperties&) = default;

public:
    bool compute(Domain* domain);
    void print();
    void print(const std::string& file_name);

public:
    inline bool isUnormalized() const { return m_unorm; }
    inline const Vector& eigenVectorScaleFactors() const { return m_unorm_scale_factors; }
    inline const Vector& centerOfMass() const { return m_center_of_mass; }
    inline const Vector& totalMass() const { return m_total_mass; }
    inline const Vector& totalFreeMass() const { return m_total_free_mass; }
    inline const Vector& eigenvalues() const { return m_eigenvalues; }
    inline const Vector& generalizedMasses() const { return m_generalized_mass_matrix; }
    inline const Matrix& modalParticipationFactors() const { return m_modal_participation_factors; }
    inline const Matrix& modalParticipationMasses() const { return m_modal_participation_masses; }
    inline const Matrix& modalParticipationMassesCumulative() const { return m_modal_participation_masses_cumulative; }
    inline const Matrix& modalParticipationMassRatios() const { return m_modal_participation_mass_ratios; }
    inline const Matrix& modalParticipationMassRatiosCumulative() const { return m_modal_participation_mass_ratios_cumulative; }

private:
    // optional flag to force eigenvector displacement-normalization
    bool m_unorm;
    // the scale factors for each eigenvector (equal to 1 if unorm==false) (size = NUM_MODES)
    Vector m_unorm_scale_factors;
    // the center of mass (size = NDM)
    Vector m_center_of_mass;
    // the total mass of the domain including the mass of fixed DOFs (size = NDF : 3 if NDM==2, 6 if NDM==3)
    Vector m_total_mass;
    // the total mass of the domain excluding the mass of fixed DOFs (size = NDF : 3 if NDM==2, 6 if NDM==3)
    Vector m_total_free_mass;
    // eigenvalues (size = NUM_MODES)
    Vector m_eigenvalues;
    // diag(V'*M*V)
    // the diagonal of the generalized mass matrix (size = NUM_MODES)
    // we store this because it the eigenvectors are not properly normalized
    // the diagonal entries won't be 1
    Vector m_generalized_mass_matrix;
    // (V'*M*R)/diag(V'*M*V)
    // the modal participation factors (size = NUM_MODES x NDF)
    Matrix m_modal_participation_factors;
    // ((V'*M*R)^2)/diag(V'*M*V)
    // the modal participation masses (size = NUM_MODES x NDF)
    Matrix m_modal_participation_masses;
    Matrix m_modal_participation_masses_cumulative;
    // ((V'*M*R)^2)/diag(V'*M*V)/m_total_free_mass*100
    // the modal participation mass ratio (%) (size = NUM_MODES x NDF)
    Matrix m_modal_participation_mass_ratios;
    Matrix m_modal_participation_mass_ratios_cumulative;
};

#endif