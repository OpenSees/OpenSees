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

#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <OPS_Globals.h>
#include <DomainModalProperties.h>
#include <Domain.h>
#include <analysis/model/AnalysisModel.h>
#include <DOF_Group.h>
#include <FE_Element.h>
#include <Element.h>
#include <ElementIter.h>
#include <Node.h>
#include <NodeIter.h>

#define DMP_VERBOSE
#define DMP_DEBUG

#define DMP_ERR_INFO "( function: " << __func__ << ", file: \"" << __FILE__ << "\", line: " << __LINE__ << " )\n"
#define DMP_ERR(X) opserr << "FATAL ERROR: " << X << DMP_ERR_INFO, exit(-1)

#define DMP_DBL_LARGE 1.0e200

// Anonymous namespace for utilities
namespace
{
    // a triplet for the sparse mass matrix
    struct triplet_t
    {
        int i;
        int j;
        double val;
        triplet_t() = default;
        triplet_t(int r, int c, double v)
            : i(r), j(c), val(v)
        {}
        inline bool operator < (const triplet_t& b) const {
            if (i < b.i) return true;
            if (i > b.i) return false;
            if (j < b.j) return true;
            if (j > b.j) return false;
            return (val < b.val);
        }
        inline bool accum(triplet_t& b) {
            if (b.i == i && b.j == j) {
                val += b.val;
                b.val = 0.0;
                b.i = b.j = -1;
                return true;
            }
            return false;
        }
    };

    // a simple structure for the sparse mass matrix
    struct sparse_matrix_t
    {
        struct row_t {
            size_t pos = 0;
            size_t count = 0;
        };
        // all triplets
        std::vector<triplet_t> triplets;
        // row access data
        std::vector<row_t> rows;

        // append a new triplet
        inline void append(int i, int j, double value) {
            triplets.push_back({ i, j, value });
        }

        // get a value, returns 0.0 if not found
        inline double get(int i, int j) const {
            if (i < 0 || i >= static_cast<int>(rows.size()))
                return 0.0;
            const auto& row_data = rows[static_cast<size_t>(i)];
            if (row_data.count > 0) {
                for (size_t c = 0; c < row_data.count; ++c) {
                    const auto& it = triplets[c + row_data.pos];
                    if (it.j == j)
                        return it.val; // found
                    if (it.j > j)
                        break; // no need to search further
                }
            }
            return 0.0;
        }

        // call this after all append calls
        inline void finish() {
            if (triplets.size() == 0)
                return;
            // sort triplets first by row, then by column
            std::sort(triplets.begin(), triplets.end());
            size_t i_current = 0;
            size_t count = 1;
            for (size_t i = 1; i < triplets.size(); ++i) {
                auto& it = triplets[i];
                if (!triplets[i_current].accum(it)) {
                    i_current = i;
                    ++count;
                }
            }
            // exclude duplicates
            std::vector<triplet_t> aux;
            aux.swap(triplets);
            triplets.resize(count);
            count = 0;
            for (const auto& it : aux)
                if (it.i > -1)
                    triplets[count++] = it;
            // make row structure for fast access
            size_t max_row = static_cast<size_t>(triplets.back().i);
            rows.resize(max_row + 1);
            for (size_t i = 0; i < triplets.size(); ++i) {
                const auto& it = triplets[i];
                auto& row_data = rows[static_cast<size_t>(it.i)];
                if (row_data.count == 0)
                    row_data.pos = i; // on first access, save the position
                ++(row_data.count); // increase counter
            }
        }
    };

    // find domain size
    int domainSize(Domain *domain) {
        int ndm = 0;
        Node* node;
        NodeIter& theNodes = domain->getNodes();
        while ((node = theNodes()) != 0) {
            int trial = node->getCrds().Size();
            if (ndm == 0)
                ndm = trial;
            if (ndm != trial)
                DMP_ERR("Cannot mix nodes with different dimensions\n");
        }
        if ((ndm != 2) && (ndm != 3))
            DMP_ERR("DomainModalProperties can be calculated only when NDM is 2 or 3, not " << ndm << "\n");
        return ndm;
    }

    // a simple structure to collect nodes in a contiguos array
    // and with a map that maps the node tag to its position in the array
    struct node_map_t
    {
        enum DOF_TYPE { DOF_EXCLUDED = -1, DOF_FIXED = -2};

        // the contiguous list of nodes
        std::vector<Node*> nodes;
        // for each node, save a re-mapped ID to the contiguous list of dofs
        // (-1 means "exluded", -2 means "fixed")
        std::vector<ID> node_ids;
        // for each node, save a list of local DOF id (0 to ndf)
        std::vector<std::vector<int> > node_u_flags;
        // for each node, save the first position in the contiguos list of dofs, ndf for each node (3 in 2D and 6 in 3D)
        std::map<int, size_t> pos;

        node_map_t(Domain* domain, int ndm, int ndf) {
            size_t n = static_cast<size_t>(domain->getNumNodes());
            if (n > 0) {
                nodes.resize(n);
                node_ids.resize(n);
                node_u_flags.resize(n);
                Node* nodePtr;
                NodeIter& theNodes = domain->getNodes();
                size_t counter = 0;
                while ((nodePtr = theNodes()) != 0) {
                    // node
                    nodes[counter] = nodePtr;
                    // map
                    pos.emplace(std::make_pair(nodePtr->getTag(), counter));
                    // id (original), needed to see if this dof is fixed
                    const ID& id_source = nodePtr->getDOF_GroupPtr()->getID();
                    // id (remapped)
                    auto& id = node_ids[counter];
                    int node_ndf = nodePtr->getNumberDOF();
                    id.resize(node_ndf);
                    int offset = static_cast<int>(counter)* ndf; // each node has ndf DOFs in this context
                    for (int j = 0; j < node_ndf; ++j) {
                        if (j < ndf) {
                            // exclude pressure dofs
                            // we can detect a pressure node (node_ndf=4) in 3d.. but not in 2d!
                            if (node_ndf == 4 && ndf == 6 && j == 3) {
                                id(j) = DOF_EXCLUDED;
                            }
                            else {
                                if (id_source(j) == -1)
                                    id(j) = DOF_FIXED;
                                else
                                    id(j) = offset;
                            }
                            ++offset;
                        }
                        else {
                            // exclude any dof >= ndf
                            id(j) = DOF_EXCLUDED;
                        }
                    }
                    // local dof flags
                    std::vector<int>& flags = node_u_flags[counter];
                    flags.resize(node_ndf);
                    for (int j = 0; j < node_ndf; ++j)
                        flags[static_cast<size_t>(j)] = id(j) == DOF_EXCLUDED ? -1 : j;
                    // go on
                    ++counter;
                }
            }

#ifdef DMP_VERBOSE
            opserr << "NODE RE-MAPPING:\n";
            opserr << "num eq: " << (int)(n * ndf) << "\n";
            for (size_t i = 0; i < n; ++i) {
                auto node = nodes[i];
                const auto& id = node_ids[i];
                const auto& uf = node_u_flags[i];
                opserr << "[" << (int)i << "] Node " << node->getTag() << " -> (";
                for (int j = 0; j < id.Size(); ++j)
                    opserr << " " << id(j);
                opserr << " ) {";
                for (size_t j = 0; j < uf.size(); ++j)
                    opserr << " " << static_cast<int>(uf[j]);
                opserr << " }\n";
            }
#endif // DMP_VERBOSE
        }

        inline size_t getPosition(int tag) const {
            auto it = pos.find(tag);
            if (it == pos.end())
                DMP_ERR("Cannot find node " << tag << "\n");
            return it->second;
        }

    };

    // a simple structure to collect element data
    struct ele_map_t
    {
        // the contiguous list of elements
        std::vector<Element*> elements;
        // for each element, save a re-mapped ID to the contiguous list of dofs
        std::vector<ID> element_ids;
        // for each element, save a list of nodal position pointing to the contiguos position of the node
        // in the node map
        std::vector<std::vector<size_t> > element_node_pos;
        // for each element, save a list of flags. true if the dof is translational, false otherwise
        std::vector<std::vector<int> > element_u_flags;

        ele_map_t(Domain* domain, const node_map_t& nm) {
            size_t count = static_cast<size_t>(domain->getNumElements());
            elements.resize(count);
            element_ids.resize(count);
            element_u_flags.resize(count);
            element_node_pos.resize(count);
            count = 0;
            {
                Element* elePtr;
                ElementIter& theEles = domain->getElements();
                while ((elePtr = theEles()) != 0) {
                    // save element
                    elements[count] = elePtr;
                    // element node tags
                    const ID& ele_nodes = elePtr->getExternalNodes();
                    // count number of dofs
                    size_t dof_count = 0;
                    for (int i = 0; i < ele_nodes.Size(); ++i) {
                        size_t pos = nm.getPosition(ele_nodes(i));
                        const ID& node_id = nm.node_ids[pos];
                        dof_count += static_cast<size_t>(node_id.Size());
                    }
                    // initialize element data
                    ID& id = element_ids[count];
                    std::vector<int>& uf = element_u_flags[count];
                    std::vector<size_t>& positions = element_node_pos[count];
                    id.resize(static_cast<int>(dof_count));
                    uf.resize(dof_count);
                    positions.resize(dof_count);
                    // remap element data
                    dof_count = 0;
                    for (int i = 0; i < ele_nodes.Size(); ++i) {
                        size_t pos = nm.getPosition(ele_nodes(i));
                        const ID& node_id = nm.node_ids[pos];
                        const std::vector<int>& node_uf = nm.node_u_flags[pos];
                        if (static_cast<int>(dof_count) + node_id.Size() > id.Size())
                            DMP_ERR("FE_Element::getID() Size < sum(size(mapped node IDs))");
                        for (int j = 0; j < node_id.Size(); ++j) {
                            id(dof_count) = node_id(j);
                            uf[dof_count] = node_uf[static_cast<size_t>(j)];
                            positions[dof_count] = pos;
                            ++dof_count;
                        }
                    }
                    ++count;
#ifdef DMP_VERBOSE
                    int dof_count_print = 0;
                    size_t uf_count_print = 0;
                    size_t pos_count_print = 0;
                    opserr << "ELEMENT " << elePtr->getTag() << " RE-MAPPING:\n";
                    for (size_t i = 0; i < ele_nodes.Size(); ++i) {
                        size_t pos = nm.getPosition(ele_nodes(i));
                        const ID& node_id = nm.node_ids[pos];
                        opserr << "   [" << ele_nodes(i) << "] (";
                        for (int j = 0; j < node_id.Size(); ++j)
                            opserr << " " << id(dof_count_print++);
                        opserr << " ) {";
                        for (size_t j = 0; j < static_cast<size_t>(node_id.Size()); ++j)
                            opserr << " " << uf[uf_count_print++];
                        opserr << " } <";
                        for (size_t j = 0; j < static_cast<size_t>(node_id.Size()); ++j)
                            opserr << " " << static_cast<int>(positions[pos_count_print++]);
                        opserr << " >\n";
                    }
#endif // DMP_VERBOSE
                }
            }
        }
    };

}

bool DomainModalProperties::compute(Domain* domain)
{
    /*
    Notes:
    1 - we assemble the global mass matrix and eigenvectors with our own (plain) numbering
    2 - we make room for rotational DOFs even if some nodes are translational only, because we
        also include the rotational effect given by the translational masses gyrating about the
        center of mass
    3 - the mass matrix may be consistent. however we need also a diagonalized version of it
        to compute the center of mass, the total mass of the structure, and the rotational masses
        due to the translational masses gyrating about the center of mass
    */

    // number of eigen-modes
    int num_eigen = domain->getEigenvalues().Size();
    if (num_eigen < 1)
        DMP_ERR("No Eigenvalue provided.\n");
    // eigenvalues
    m_eigenvalues = domain->getEigenvalues();
    // number of dimensions
    int ndm = domainSize(domain);
    // max number of DOFs per node, translational and rotational only
    int ndf = ndm == 2 ? 3 : 6;
    // number of nodes
    int num_nodes = domain->getNumNodes();
    // number of equations (not the real one, include rotational dofs even if not present)
    int num_eq = num_nodes * ndf;

    // initialize members
    m_center_of_mass.resize(ndm);
    m_total_mass.resize(ndf);
    m_total_free_mass.resize(ndf);
    m_generalized_mass_matrix.resize(num_eigen);
    m_modal_participation_factors.resize(num_eigen, ndf);
    m_modal_participation_masses.resize(num_eigen, ndf);
    m_modal_participation_masses_cumulative.resize(num_eigen, ndf);
    m_modal_participation_mass_ratios.resize(num_eigen, ndf);
    m_modal_participation_mass_ratios_cumulative.resize(num_eigen, ndf);

    // map all nodes
    node_map_t nodemap(domain, ndm, ndf);

    // map all elements
    ele_map_t elemap(domain, nodemap);

    // the sparse mass matrix
    sparse_matrix_t M;
    // the equivalent diagonalized masses for each node (total)
    Matrix ML(num_nodes, ndf);
    // the equivalent diagonalized masses for each node (only at free DOFs)
    Matrix MLfree(num_nodes, ndf);
    // the array of eigenvectors
    std::vector<Vector> V(num_eigen);
    for (auto& iV : V) {
        iV.resize(num_eq);
        iV.Zero();
    }

    // asseble in the sparse mass matrix
    Vector aux_MD; // auxiliary: sum of each row of iM
    Vector aux_ML(ndf); // auxiliary: sum of MD for each DOF
    Vector aux_MC(ndf); // auxiliary: sum of consistent mass iM diagonal terms only
    Vector aux_C(ndf); // auxiliary: scale factors
    Vector aux_Lumped; // auxiliary: lumped diagonalized mass matrix
    auto addM = [&M, &ML, &MLfree, &aux_MD, &aux_ML, &aux_MC, &aux_C, &aux_Lumped](
        const Matrix& iM, const ID& iD, 
        const std::vector<int> &uflags, const std::vector<size_t> &npos) {
        int n = iD.Size();
        if (iM.noRows() != n || iM.noCols() != n)
            DMP_ERR("Error: inconsistent mass matrix and ID\n");
        for (int i = 0; i < n; ++i) {
            int iloc = iD(i);
            if (iloc >= 0) {
                for (int j = 0; j < n; ++j) {
                    int jloc = iD(j);
                    if (jloc >= 0) {
                        double value = iM(i, j);
                        if (value != 0.0)
                            M.append(iloc, jloc, value);
                    }
                }
            }
        }
        // compute diagonal terms summing each row (it can have negative terms, so
        // it cannot be directly used as a lumping approach)
        aux_MD.resize(n);
        aux_MD.Zero();
        for (int i = 0; i < n; ++i) {
            if (uflags[static_cast<size_t>(i)] >= 0) {
                for (int j = 0; j < n; ++j) {
                    aux_MD(i) += iM(i, j);
                }
            }
        }
        // sum terms in MD DOF by DOF. This will give the total mass of the element
        // for each DOF
        aux_ML.Zero();
        for (int i = 0; i < n; i++) {
            int local_dof = uflags[static_cast<size_t>(i)];
            if (local_dof >= 0)
                aux_ML(local_dof) += aux_MD(i);
        }
        // sum consistent diagonal terms of iM DOF by DOF
        aux_MC.Zero();
        for (int i = 0; i < n; i++) {
            int local_dof = uflags[static_cast<size_t>(i)];
            if (local_dof >= 0)
                aux_MC(local_dof) += iM(i, i);
        }
        // obtain scale factors
        for (int i = 0; i < aux_C.Size(); ++i) {
            double ml = aux_ML(i);
            double mc = aux_MC(i);
            aux_C(i) = std::abs(mc) > 0.0 ? ml / mc : 0.0;
        }
        // obtain lumped mass matrix
        aux_Lumped.resize(n);
        aux_Lumped.Zero();
        for (int i = 0; i < n; ++i) {
            int local_dof = uflags[static_cast<size_t>(i)];
            if (local_dof >= 0) {
                aux_Lumped(i) = iM(i, i) * aux_C(local_dof);
                int node_pos = npos[static_cast<size_t>(i)];
                ML(node_pos, local_dof) += aux_Lumped(i);
                if (iD(i) >= 0)
                    MLfree(node_pos, local_dof) += aux_Lumped(i);
            }
        }
    };

    // assemble in the array of eigenvectors
    auto addV = [&V, num_eigen](const Matrix& iV, const ID& iD) {
        int n = iD.Size();
        if (iV.noRows() != n || iV.noCols() != static_cast<int>(V.size()))
            DMP_ERR("Error: inconsistent eigenvector matrix and ID\n");
        for (int j = 0; j < static_cast<int>(V.size()); ++j) {
            auto& jV = V[static_cast<size_t>(j)];
            for (int i = 0; i < n; ++i) {
                int iloc = iD(i);
                if (iloc >= 0)
                    jV(iloc) = iV(i, j);
            }
        }
    };

    // assemble all element contributions
    {
        for (size_t i = 0; i < elemap.elements.size(); ++i) {
            Element* element = elemap.elements[i];
            const ID& iD = elemap.element_ids[i];
            const std::vector<int>& uflags = elemap.element_u_flags[i];
            const std::vector<size_t>& npos = elemap.element_node_pos[i];
            const Matrix& iM = element->getMass();
            addM(iM, iD, uflags, npos);
        }
    }
    // assemble all nodal contributions
    {
        std::vector<size_t> npos;
        for (size_t i = 0; i < nodemap.nodes.size(); ++i) {
            Node* node = nodemap.nodes[i];
            const ID& iD = nodemap.node_ids[i];
            const std::vector<int>& uflags = nodemap.node_u_flags[i];
            npos.resize(static_cast<size_t>(iD.Size()));
            std::fill(npos.begin(), npos.end(), i);
            const Matrix& iV = node->getEigenvectors();
            const Matrix& iM = node->getMass();
            addM(iM, iD, uflags, npos);
            addV(iV, iD);
        }
    }

    // done assemling M
    M.finish();
#ifdef DMP_DEBUG
    for (const auto& it : M.triplets) {
        if (it.i < 0 || it.i >= num_eq || it.j < 0 || it.j >= num_eq)
            DMP_ERR("Triplet indices are out of bounds\n");
    }
#endif // DMP_DEBUG

    // compute the center of mass.
    // note that if the total mass in one of the global directions is 0
    // we sould take the geometric center, as otherwise it will be 0.
    // TODO: check wether we need the total or free mass.
    // accordingly, for the geometric center... should we use free nodes only?
    // let's start with free nodes only... by the way they are the only ones
    // accounted for by the eigenvalue analysis...
    {
        Vector geometric_center(ndm);
        Vector com_weight(ndm); // weight for center of mass
        Vector cog_weight(ndm); // weight for geometric center
        m_center_of_mass.Zero();
        for (int i = 0; i < num_nodes; ++i) {
            Node* node = nodemap.nodes[static_cast<size_t>(i)];
            const ID& ids = nodemap.node_ids[static_cast<size_t>(i)];
            const Vector& pos = node->getCrds();
            for (int j = 0; j < ndm; ++j) {
                if ((j < ids.Size()) && (ids(j) >= 0)) {
                    double Mij = MLfree(i, j); // use the free masses only
                    double Pi = pos(j);
                    geometric_center(j) += Pi;
                    cog_weight(j) += 1.0;
                    m_center_of_mass(j) += Pi * Mij;
                    com_weight(j) += Mij;
                }
            }
        }
        for (int j = 0; j < ndm; ++j) {
            // finish computing geometric center
            if (cog_weight(j) > 0.0)
                geometric_center(j) /= cog_weight(j);
            // finish computing center of mass.
            // select the geometric center if not mass is provided in this component
            if (com_weight(j) > 0.0)
                m_center_of_mass(j) /= com_weight(j);
            else
                m_center_of_mass(j) = geometric_center(j);
        }
    }

    // Now the sparse matrix M contains the orignal mass (consistent or lumped).
    // ML is a lumped version of M. for rotational DOFs it contains only the 
    // rotary masses directly input by the user.
    // However we also need to compute, for the modal properties, the rotational
    // masses due to the translational masses gyrating about the center of mass.
    // Of course this goes only in ML, not in M.
    auto compute_extra_rotary_mass = [&nodemap, num_nodes, ndf, this](Matrix& iML) {
        for (int i = 0; i < num_nodes; ++i) {
            Node* node = nodemap.nodes[static_cast<size_t>(i)];
            const ID& ids = nodemap.node_ids[static_cast<size_t>(i)];
            const Vector& pos = node->getCrds();
            double dx = pos(0) - m_center_of_mass(0);
            double dy = pos(1) - m_center_of_mass(1);
            double mx = iML(i, 0);
            double my = iML(i, 1);
            if (ndf == 3) { // 2D case
                iML(i, 2) += dx * dx * my + dy * dy * mx;
            }
            else { // 3D case
                double dz = pos(2) - m_center_of_mass(2);
                double mz = iML(i, 2);
                iML(i, 3) += dy * dy * mz + dz * dz * my;
                iML(i, 4) += dx * dx * mz + dz * dz * mx;
                iML(i, 5) += dx * dx * my + dy * dy * mx;
            }
        }
    };
    compute_extra_rotary_mass(ML);
    compute_extra_rotary_mass(MLfree);

    // compute the total mass of the domain (total and free-only)
    m_total_mass.Zero();
    m_total_free_mass.Zero();
    for (int j = 0; j < ndf; ++j) {
        double msum = 0.0;
        double msum_free = 0.0;
        for (int i = 0; i < num_nodes; ++i) {
            msum += ML(i, j);
            msum_free += MLfree(i, j);
        }
        m_total_mass(j) = msum;
        m_total_free_mass(j) = msum_free;
    }

    // now we can compute all the modal masses and participation factors.
    // we can do it mode by mode due to their orthogonality.
    
    // a temporary to store the V'*M product
    Vector jVTM(num_eq);
    // defines the magnitude of the rigid body response of a DOF 
    // to imposed rigid body motion (displacement or infinitesimal rotation) in the i-direction.
    // each (ndf x 1) block correspond to a node (ordered sequentially as in nodemap), and it is defined as:
    // | 1   0   0   0  dz -dz | |e1|
    // | 0   1   0 -dz   0  dx | |e2|
    // | 0   0   1  dy -dx   0 | |e3|
    // | 0   0   0   1   0   0 | |e4|
    // | 0   0   0   0   1   0 | |e5|
    // | 0   0   0   0   0   1 | |e6|
    // where ei is 1, and all the other are 0.
    Vector R(num_eq);
    // for each j mode...
    for (int j = 0; j < num_eigen; ++j) {

        // current mode eigenvector
        const Vector& jV = V[j];

        // compute [V' * M] : a temporary to avoid extra calculations
        jVTM.Zero();
        for (const auto& it : M.triplets) 
            jVTM[it.j] += it.val * jV(it.i);

        // compute [V' * M * V] : generalized mass matrix
        double GM = jVTM ^ jV;
        m_generalized_mass_matrix(j) = GM;
        double invGM = GM == 0.0 ? DMP_DBL_LARGE : 1.0 / GM;

        // for each DOF ...
        for (int i = 0; i < ndf; ++i) {

            // compute the rigid-body-mode vector R
            R.Zero();
            for (int inode = 0; inode < num_nodes; ++inode) {
                int index = inode * ndf;
                R(index + i) = 1.0; // one at the current DOF for direct translational or rotational effects
                if (i >= ndm) {
                    const Vector& pos = nodemap.nodes[static_cast<size_t>(inode)]->getCrds();
                    double dx = pos(0) - m_center_of_mass(0);
                    double dy = pos(1) - m_center_of_mass(1);
                    if (ndf == 3) { // 2D case
                        if(i == 2) {
                            R(index + 0) = -dy;
                            R(index + 1) = dx;
                        }
                    }
                    else { // 3D case
                        double dz = pos(2) - m_center_of_mass(2);
                        if (i == 3) {
                            R(index + 1) = -dz;
                            R(index + 2) = dy;
                        }
                        else if (i == 4) {
                            R(index + 0) = dz;
                            R(index + 2) = -dx;
                        }
                        else if (i == 5) {
                            R(index + 0) = -dy;
                            R(index + 1) = dx;
                        }
                    }
                }
            }

            // compute [V' * M * R] : generalized load vector
            double L = jVTM ^ R;

            // compute [(V' * M * R) / diag(V' * M * V)] : modal participation factors
            m_modal_participation_factors(j, i) = L * invGM;

            // compute [(V' * M * R)^2 / diag(V' * M * V)] : effective modal masses
            m_modal_participation_masses(j, i) = L * L * invGM;
        }
    }

    // now we can compute the cumulative masses and ratios
    for (int j = 0; j < ndf; ++j) {
        double tot_mass = m_total_free_mass(j);
        double inv_tot_mass = tot_mass == 0.0 ? DMP_DBL_LARGE : 1.0 / tot_mass;
        double accum_1 = 0.0;
        double accum_2 = 0.0;
        for (int i = 0; i < num_eigen; ++i) {
            double mpm = m_modal_participation_masses(i, j);
            double mpmr = mpm * inv_tot_mass;
            accum_1 += mpm;
            accum_2 += mpmr;
            m_modal_participation_mass_ratios(i, j) = mpmr;
            m_modal_participation_masses_cumulative(i, j) = accum_1;
            m_modal_participation_mass_ratios_cumulative(i, j) = accum_2;
        }
    }

    // done
    return true;
}

namespace
{

#define DMP_OUT_COMMENT "#"
#define DMP_OUT_RECORD "*"
#define DMP_OUT_BLANK " "
#define DMP_OUT_FLOAT(X) std::setw(12) << std::setprecision(6) << X
#define DMP_OUT_GEN(X) std::setw(12) << X

    template<class TStream>
    void print_internal(TStream& out, const DomainModalProperties& dmp)
    {
        // header
        out << DMP_OUT_COMMENT << " MODAL ANALYSIS REPORT\n\n";

        // problem size
        out << DMP_OUT_RECORD << " 1. DOMAIN SIZE:\n" << dmp.centerOfMass().Size() << "\n\n";
        
        // eigenvalues and derived quantities
        out << DMP_OUT_RECORD << " 2. EIGENVALUE ANALYSIS:\n"
            << DMP_OUT_COMMENT
            << DMP_OUT_GEN("MODE") << DMP_OUT_GEN("LAMBDA") << DMP_OUT_GEN("OMEGA")
            << DMP_OUT_GEN("FREQUENCY") << DMP_OUT_GEN("PERIOD") << "\n";
    }
}

void DomainModalProperties::print()
{
    std::stringstream ss;
    print_internal(ss, *this);
    std::string s = ss.str();
    opserr << s.c_str();
}

void DomainModalProperties::print(const std::string& file_name)
{
    std::ofstream ss(file_name);
    if (!ss.is_open())
        DMP_ERR("Cannot open file \"" << file_name.c_str() << "\"\n");
    print_internal(ss, *this);
    ss.close();
}


