#pragma once

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

const double pi = 3.14159265358979323846;
const double twopi = 2.0 * pi;

void eraseDoubleArrayIndex(int len, double** doubleArr, int index) {
    for (int i = 0; i < len; i++) {
        for (int j = index + 1; j < len; j++) {
            doubleArr[i][j-1] = doubleArr[i][j];
        }
        doubleArr[i][len-1] = -1.0;
    }
}

void eraseDoubleArrayRow(int len, double** doubleArr, int index) {
    for (int i = index + 1; i < len; i++) {
        for (int j = 0; j < len; j++) {
            doubleArr[i-1][j] = doubleArr[i][j];
        }
    }

    for (int j = 0; j < len; j++) doubleArr[len-1][j] = -1.0;
}

void printArr(double** arr, int n, int m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {

            // printing the value of an element
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }
}

void orderJetsByDecreasingPt(vector<double>*& in_etas,
                     vector<double>*& in_phis,
                     vector<double>*& in_pts,
                     vector<double>*& in_ms) {
    vector<double> working_etas(*in_etas);
    vector<double> working_phis(*in_phis);
    vector<double> working_pts(*in_pts);
    vector<double> working_ms(*in_ms);

    vector<double>* final_jet_etas = new vector<double>();
    vector<double>* final_jet_phis = new vector<double>();
    vector<double>* final_jet_pts = new vector<double>();
    vector<double>* final_jet_ms = new vector<double>();

    while (working_etas.size() > 0) {
        double largest_pt = 0;
        int largest_pt_index = -1;
        for (int i = 0; i < working_etas.size(); i++) {
            if (working_pts[i] > largest_pt) {
                largest_pt = working_pts[i];
                largest_pt_index = i;
            }
        }

        final_jet_etas->emplace_back(working_etas[largest_pt_index]);
        final_jet_phis->emplace_back(working_phis[largest_pt_index]);
        final_jet_pts->emplace_back(working_pts[largest_pt_index]);
        final_jet_ms->emplace_back(working_ms[largest_pt_index]);

        working_etas.erase(working_etas.begin() + largest_pt_index);
        working_phis.erase(working_phis.begin() + largest_pt_index);
        working_pts.erase(working_pts.begin() + largest_pt_index);
        working_ms.erase(working_ms.begin() + largest_pt_index);
    }

    delete in_etas;
    delete in_phis;
    delete in_pts;
    delete in_ms;

    in_etas = final_jet_etas;
    in_phis = final_jet_phis;
    in_pts = final_jet_pts;
    in_ms = final_jet_ms;
}

// -----------------------------------------------------------------------------

double rapidityFromPseudoRapidity(double eta, double phi, double pt, double m) {
    double numerator = sqrt((m * m) + (pt * pt * cosh(eta) * cosh(eta))) + (pt * sinh(eta));
    double denominator = sqrt((m * m) + (pt * pt));
    return log(numerator/denominator);
}

// https://en.wikipedia.org/wiki/Pseudorapidity
// https://en.wikipedia.org/wiki/Energy-momentum_relation
void etaPhiPtToFourMomentum(double eta, double phi, double pt, double m, double& px, double& py, double& pz, double& E) {
    px = pt * cos(phi);
    py = pt * sin(phi);
    pz = pt * sinh(eta);

    double p_mag_squared = px * px + py * py + pz * pz;
    E = sqrt(p_mag_squared + m * m);
}

void etaPhiPtToFourMomentumNonPseudo(double eta, double phi, double pt, double m, double& px, double& py, double& pz, double& E) {
    px = pt * cos(phi);
    py = pt * sin(phi);
    
    double y = rapidityFromPseudoRapidity(eta, phi, pt, m);
    double mT = sqrt((pt * pt) + (m * m));
    pz = mT * sinh(y);
    E = mT * cosh(y);
}

void fourMomentumToEtaPhiPt(double px, double py, double pz, double E, double& eta, double& phi, double& pt, double& m) {
    pt = sqrt(px * px + py * py);
    phi = atan2(py, px);

    double p_mag_squared = px * px + py * py + pz * pz;
    eta = atanh(pz / sqrt(p_mag_squared));
    m = sqrt(E * E - p_mag_squared);
}

// Verified with FastJet library
void addFourMomenta(double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2, double& eta_result, double& phi_result, double& pt_result, double& m_result) {
    // 4-Momentum #1 Cartesian Conversion
    double px1, py1, pz1, E1;
    etaPhiPtToFourMomentum(eta1, phi1, pt1, m1, px1, py1, pz1, E1);

    // 4-Momentum #2 Cartesian Conversion
    double px2, py2, pz2, E2;
    etaPhiPtToFourMomentum(eta2, phi2, pt2, m2, px2, py2, pz2, E2);

    // 4-Momentum Addition
    double px_result = px1 + px2;
    double py_result = py1 + py2;
    double pz_result = pz1 + pz2;
    double E_result = E1 + E2;

    // 4-Momentum conversion to (eta, phi) space
    fourMomentumToEtaPhiPt(px_result, py_result, pz_result, E_result, eta_result, phi_result, pt_result, m_result);
    // https://en.wikipedia.org/wiki/Invariant_mass#Collider_experiments
    //  m_result = sqrt(2 * pt1 * pt2 * (cosh(eta1 - eta2) - cos(phi1 - phi2)));
}

// -----------------------------------------------------------------------------

// returns the distance squared in (eta, phi) space
double distanceBetweenJetsSquaredPseudo(double eta1, double phi1, double eta2, double phi2) {
    double differenceEta = eta2 - eta1;
    double differencePhi = abs(phi2 - phi1); if (differencePhi > pi) differencePhi = twopi - differencePhi;

    return (differenceEta * differenceEta) + (differencePhi * differencePhi);
}

// returns the distance squared in (y, phi) space
double distanceBetweenJetsSquared(double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2) {
    double differenceY = rapidityFromPseudoRapidity(eta2, phi2, pt2, m2) - rapidityFromPseudoRapidity(eta1, phi1, pt1, m1);
    double differencePhi = abs(phi2 - phi1); if (differencePhi > pi) differencePhi = twopi - differencePhi;
    return (differenceY * differenceY) + (differencePhi * differencePhi);
}

// returns the distance in (eta, phi) space
double distanceBetweenJets(double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2) {
    return sqrt(distanceBetweenJetsSquared(eta1, phi1, pt1, m1, eta2, phi2, pt2, m2));
}

// based off of this paper https://arxiv.org/pdf/1407.2922

// anti-kT has n = -1
double jetReclusteringMetric_anti_kT(double R, double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2) {
    return min(1.0/(pt1 * pt1), 1.0/(pt2 * pt2)) * distanceBetweenJetsSquared(eta1, phi1, pt1, m1, eta2, phi2, pt2, m2) / (R * R);
}

// Cambridge/Aachen (C/A) has n = 0
double jetReclusteringMetric_CA(double R, double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2) {
    return distanceBetweenJetsSquared(eta1, phi1, pt1, m1, eta2, phi2, pt2, m2) / (R * R);
}
// kT has n = 1
double jetReclusteringMetric_kT(double R, double eta1, double phi1, double pt1, double m1, double eta2, double phi2, double pt2, double m2) {
    return min((pt1 * pt1), (pt2 * pt2)) * distanceBetweenJetsSquared(eta1, phi1, pt1, m1, eta2, phi2, pt2, m2) / (R * R);
}

// -----------------------------------------------------------------------------

// https://fastjet.fr/repo/fastjet-doc-3.4.3.pdf
void antiKtAlgorithm(double R,
                     vector<double>* in_jet_etas,
                     vector<double>* in_jet_phis,
                     vector<double>* in_jet_pts,
                     vector<double>* in_jet_ms,
                     vector<double>*& result_jet_etas,
                     vector<double>*& result_jet_phis,
                     vector<double>*& result_jet_pts,
                     vector<double>*& result_jet_ms) {
    // R__ASSERT(in_jet_etas); R__ASSERT(in_jet_phis); R__ASSERT(in_jet_pts); R__ASSERT(in_jet_ms);

    // R__ASSERT(in_jet_etas->size() == in_jet_phis->size());
    // R__ASSERT(in_jet_etas->size() == in_jet_pts->size());
    // R__ASSERT(in_jet_etas->size() == in_jet_ms->size());
    const int indexB = -2;

    vector<double> working_jet_etas(*in_jet_etas);
    vector<double> working_jet_phis(*in_jet_phis);
    vector<double> working_jet_pts(*in_jet_pts);
    vector<double> working_jet_ms(*in_jet_ms);

    int numParticles = working_jet_etas.size();

    vector<double>* final_jet_etas = new vector<double>();
    vector<double>* final_jet_phis = new vector<double>();
    vector<double>* final_jet_pts = new vector<double>();
    vector<double>* final_jet_ms = new vector<double>();

    int stepNum = 0;
    while (numParticles > 0) {
        // STEP 1: Find the minimum distance between two jets / jet and beam
        //cout << "Step " << stepNum++ << ": numParticles=" << numParticles << " | final_jet_etas.size()=" << final_jet_etas->size() << endl;

        double minDistance = HUGE_VALF; int index_i_min = -1; int index_j_min = -1; 
        for (int i = 0; i < numParticles; i++) {
            double eta_i = working_jet_etas.at(i);
            double phi_i = working_jet_phis.at(i);
            double pt_i = working_jet_pts.at(i);
            double m_i = working_jet_ms.at(i);

            for (int j = i + 1; j < numParticles; j++) {
                double eta_j = working_jet_etas.at(j);
                double phi_j = working_jet_phis.at(j);
                double pt_j = working_jet_pts.at(j);
                double m_j = working_jet_ms.at(j);

                double d_ij = jetReclusteringMetric_anti_kT(R, eta_i, phi_i, pt_i, m_i, eta_j, phi_j, pt_j, m_j);
                // cout << "Anti-kt Metric between jet " << i << " (" << eta_i << " " << phi_i << " " << pt_i << " " << m_i << ") and jet " << j << " (" << eta_j << " " << phi_j << " " << pt_j << " " << m_j << "): " << d_ij << endl;

                if (d_ij < minDistance) {
                    minDistance = d_ij;
                    index_i_min = i;
                    index_j_min = j;
                }
            }

            double d_iB = 1.0 / (pt_i * pt_i);
            // cout << "Anti-kt Metric between jet " << i << " (" << eta_i << " " << phi_i << " " << pt_i << " " << m_i << ") and beam: " << d_iB << endl;
            if (d_iB < minDistance) {
                minDistance = d_iB;
                index_i_min = i;
                index_j_min = indexB;
            }
        }

        // STEP 2: Combine Closests Jets
        if (index_j_min == indexB) {

            // Take out jet i and place it in the final jet array
            double eta_i = working_jet_etas.at(index_i_min);
            double phi_i = working_jet_phis.at(index_i_min);
            double pt_i = working_jet_pts.at(index_i_min);
            double m_i = working_jet_ms.at(index_i_min);

            final_jet_etas->emplace_back(eta_i);
            final_jet_phis->emplace_back(phi_i);
            final_jet_pts->emplace_back(pt_i);
            final_jet_ms->emplace_back(m_i);

            // Erase the data for jet i
            working_jet_etas.erase(working_jet_etas.begin() + index_i_min);
            working_jet_phis.erase(working_jet_phis.begin() + index_i_min);
            working_jet_pts.erase(working_jet_pts.begin() + index_i_min);
            working_jet_ms.erase(working_jet_ms.begin() + index_i_min);
        }
        else {

            // Combine jets i and j
            double eta_i = working_jet_etas.at(index_i_min);
            double phi_i = working_jet_phis.at(index_i_min);
            double pt_i = working_jet_pts.at(index_i_min);
            double m_i = working_jet_ms.at(index_i_min);

            double eta_j = working_jet_etas.at(index_j_min);
            double phi_j = working_jet_phis.at(index_j_min);
            double pt_j = working_jet_pts.at(index_j_min);
            double m_j = working_jet_ms.at(index_j_min);

            double new_eta_i, new_phi_i, new_pt_i, new_m_i;
            addFourMomenta(eta_i, phi_i, pt_i, m_i, eta_j, phi_j, pt_j, m_j, new_eta_i, new_phi_i, new_pt_i, new_m_i);

            // Put the combination of jets i and j in index i
            working_jet_etas.at(index_i_min) = new_eta_i;
            working_jet_phis.at(index_i_min) = new_phi_i;
            working_jet_pts.at(index_i_min) = new_pt_i;
            working_jet_ms.at(index_i_min) = new_m_i;

            // Erase the data for jet j
            working_jet_etas.erase(working_jet_etas.begin() + index_j_min);
            working_jet_phis.erase(working_jet_phis.begin() + index_j_min);
            working_jet_pts.erase(working_jet_pts.begin() + index_j_min);
            working_jet_ms.erase(working_jet_ms.begin() + index_j_min);
        }

        numParticles = working_jet_etas.size();
    }

    orderJetsByDecreasingPt(final_jet_etas, final_jet_phis, final_jet_pts, final_jet_ms);

    result_jet_etas = final_jet_etas;
    result_jet_phis = final_jet_phis;
    result_jet_pts = final_jet_pts;
    result_jet_ms = final_jet_ms;
}


//#define DEBUG_FAST_ALGO

// Adding in optimization where distances are stored
void antiKtAlgorithm_fast(double R,
                     vector<double>* in_jet_etas,
                     vector<double>* in_jet_phis,
                     vector<double>* in_jet_pts,
                     vector<double>* in_jet_ms,
                     vector<double>*& result_jet_etas,
                     vector<double>*& result_jet_phis,
                     vector<double>*& result_jet_pts,
                     vector<double>*& result_jet_ms) {
    const int indexB = -2;

    vector<double> working_jet_etas(*in_jet_etas);
    vector<double> working_jet_phis(*in_jet_phis);
    vector<double> working_jet_pts(*in_jet_pts);
    vector<double> working_jet_ms(*in_jet_ms);

    int numParticles = working_jet_etas.size();
    const int originalNumParticles = numParticles;

    vector<double>* final_jet_etas = new vector<double>();
    vector<double>* final_jet_phis = new vector<double>();
    vector<double>* final_jet_pts = new vector<double>();
    vector<double>* final_jet_ms = new vector<double>();

    double** fastAlgoStoredDistances = new double*[originalNumParticles];
    for (int i = 0; i < originalNumParticles; i++) {
        fastAlgoStoredDistances[i] = new double[originalNumParticles];
        for (int j = 0; j < originalNumParticles; j++) fastAlgoStoredDistances[i][j] = -1.0;
    }

    int stepNum = 0;
    while (numParticles > 0) {
        // STEP 1: Find the minimum distance between two jets / jet and beam
        
#ifdef DEBUG_FAST_ALGO
        cout << "Step " << stepNum++ << ": numParticles=" << numParticles << " | final_jet_etas.size()=" << final_jet_etas->size() << endl;
        cout << "Jet Pts: ";
        for (int i = 0; i < working_jet_pts.size(); i++) cout << working_jet_pts[i] << " ";
        cout << endl << endl;
        printArr(fastAlgoStoredDistances, working_jet_etas.size(), working_jet_etas.size());
        cout << endl;
#endif

        double minDistance = HUGE_VAL; int index_i_min = -1; int index_j_min = -1; 
        for (int i = 0; i < numParticles; i++) {
            double pt_i = working_jet_pts.at(i);

            for (int j = i + 1; j < numParticles; j++) {
                double d_ij = fastAlgoStoredDistances[i][j];
                if (d_ij == -1.0) {
                    double eta_i = working_jet_etas.at(i);
                    double phi_i = working_jet_phis.at(i);
                    double m_i = working_jet_ms.at(i);

                    double eta_j = working_jet_etas.at(j);
                    double phi_j = working_jet_phis.at(j);
                    double pt_j = working_jet_pts.at(j);
                    double m_j = working_jet_ms.at(j);

                    d_ij = jetReclusteringMetric_anti_kT(R, eta_i, phi_i, pt_i, m_i, eta_j, phi_j, pt_j, m_j);
                    fastAlgoStoredDistances[i][j] = d_ij;
                }

                if (d_ij < minDistance) {
                    minDistance = d_ij;
                    index_i_min = i;
                    index_j_min = j;
                }
            }

            double d_iB = 1.0 / (pt_i * pt_i);
            if (d_iB < minDistance) {
                minDistance = d_iB;
                index_i_min = i;
                index_j_min = indexB;
            }
        }

#ifdef DEBUG_FAST_ALGO
        cout << "After computing new distances" << endl;
        cout << endl;
        printArr(fastAlgoStoredDistances, working_jet_etas.size(), working_jet_etas.size());
        cout << endl;

        cout << "index_i_min = " << index_i_min << " | index_j_min = " << index_j_min << endl;
#endif

        // STEP 2: Combine / Finalize Smallest Dist. Jet
        if (index_j_min == indexB) { // Take one jet out as final

            // Take out jet i and place it in the final jet array
            double eta_i = working_jet_etas.at(index_i_min);
            double phi_i = working_jet_phis.at(index_i_min);
            double pt_i = working_jet_pts.at(index_i_min);
            double m_i = working_jet_ms.at(index_i_min);

            final_jet_etas->emplace_back(eta_i);
            final_jet_phis->emplace_back(phi_i);
            final_jet_pts->emplace_back(pt_i);
            final_jet_ms->emplace_back(m_i);

            working_jet_etas.erase(working_jet_etas.begin() + index_i_min);
            working_jet_phis.erase(working_jet_phis.begin() + index_i_min);
            working_jet_pts.erase(working_jet_pts.begin() + index_i_min);
            working_jet_ms.erase(working_jet_ms.begin() + index_i_min);

            // Erase the metrics at the row and column of the ith jet and move all greater indexes back by 1
            eraseDoubleArrayIndex(numParticles, fastAlgoStoredDistances, index_i_min);
            eraseDoubleArrayRow(numParticles, fastAlgoStoredDistances, index_i_min);
        }
        else { // Combine Closests Jets

            // Combine jets i and j
            double eta_i = working_jet_etas.at(index_i_min);
            double phi_i = working_jet_phis.at(index_i_min);
            double pt_i = working_jet_pts.at(index_i_min);
            double m_i = working_jet_ms.at(index_i_min);

            double eta_j = working_jet_etas.at(index_j_min);
            double phi_j = working_jet_phis.at(index_j_min);
            double pt_j = working_jet_pts.at(index_j_min);
            double m_j = working_jet_ms.at(index_j_min);

            double new_eta_i, new_phi_i, new_pt_i, new_m_i;
            addFourMomenta(eta_i, phi_i, pt_i, m_i, eta_j, phi_j, pt_j, m_j, new_eta_i, new_phi_i, new_pt_i, new_m_i);

            // Put the combination of jets i and j in index i
            working_jet_etas.at(index_i_min) = new_eta_i;
            working_jet_phis.at(index_i_min) = new_phi_i;
            working_jet_pts.at(index_i_min) = new_pt_i;
            working_jet_ms.at(index_i_min) = new_m_i;

            // Erase the data for jet j
            working_jet_etas.erase(working_jet_etas.begin() + index_j_min);
            working_jet_phis.erase(working_jet_phis.begin() + index_j_min);
            working_jet_pts.erase(working_jet_pts.begin() + index_j_min);
            working_jet_ms.erase(working_jet_ms.begin() + index_j_min);

            // Erase the metrics at the row and column of the ith jet and move all greater indexes back by 1
            eraseDoubleArrayIndex(numParticles, fastAlgoStoredDistances, index_j_min);
            eraseDoubleArrayRow(numParticles, fastAlgoStoredDistances, index_j_min);
            
            // Erase the stored for values for the old metrics of other jets with the jet at position i
            for (int i = 0; i < numParticles; i++) fastAlgoStoredDistances[i][index_i_min] = -1;
            for (int j = 0; j < numParticles; j++) fastAlgoStoredDistances[index_i_min][j] = -1;
        }

        numParticles = working_jet_etas.size();
    }

    orderJetsByDecreasingPt(final_jet_etas, final_jet_phis, final_jet_pts, final_jet_ms);

    result_jet_etas = final_jet_etas;
    result_jet_phis = final_jet_phis;
    result_jet_pts = final_jet_pts;
    result_jet_ms = final_jet_ms;
}

// -----------------------------------------------------------------------------

void JetMath() {
    cout << "Math for Jet Reclustering" << endl;
}

