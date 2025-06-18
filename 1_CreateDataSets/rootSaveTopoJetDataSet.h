#include <iostream>
#include <vector>
#include "TChain.h"

#include "../JetMath.h"

using namespace std;

void WriteDataSet(string name, vector<double>* etas,
                     vector<double>* phis,
                     vector<double>* pts,
                     vector<double>* ms) {
    ofstream testdatafile;
    testdatafile.open(("../../SavedDataSets/" + name + ".bin").c_str());
    
    int numMembers = etas->size();
    testdatafile.write(reinterpret_cast<const char*>(&numMembers), sizeof(int));
    testdatafile.write(reinterpret_cast<const char*>(etas->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(phis->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(pts->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(ms->data()), numMembers*sizeof(double));
    testdatafile.close();
}

void save_info_topoclusters(TChain& inputDataChain, std::string& data_prefix) {
    std::cout << "TODO" << std::endl;
}


void save_info_jets(TChain& inputDataChain, std::string& data_prefix) {
    vector<float> *eta_pseudorapidity = nullptr, *phi_azimuth_angle = nullptr, *mass_MeV = nullptr, *pt_transverse_momentum = nullptr;

    std::string branch_pt_name = data_prefix + "_pt";
    std::string branch_eta_name = data_prefix + "_eta";
    std::string branch_phi_name = data_prefix + "_phi";
    std::string branch_m_name = data_prefix + "_m";

    inputDataChain.SetBranchAddress(branch_pt_name.c_str(), &pt_transverse_momentum);
    inputDataChain.SetBranchAddress(branch_eta_name.c_str(), &eta_pseudorapidity);
    inputDataChain.SetBranchAddress(branch_phi_name.c_str(), &phi_azimuth_angle);
    inputDataChain.SetBranchAddress(branch_m_name.c_str(), &mass_MeV);

    const int n = inputDataChain.GetEntries();
    vector<double> masses_GeV;
    vector<double> colors_transverse_momentum;
    vector<double> x_eta_pseudorapidities;
    // vector<double> x_theta_polar_angles = vector<double>(n * 100);
    vector<double> y_phi_azimuthal_angles;

    std::cout << "How Many Jets Should Be Included in this Dataset? ";
    long numJetsLimit = 0; cin >> numJetsLimit;
    std::cout << "Give the DataSet a Name: ";
    std::string dataSetName = ""; cin >> dataSetName;

    long totalNumJets = 0;
    for (size_t i = 0; totalNumJets < numJetsLimit; i++){
        inputDataChain.GetEntry(i);

        assert(mass_MeV->size() == phi_azimuth_angle->size());
        assert(mass_MeV->size() == eta_pseudorapidity->size());
        assert(mass_MeV->size() == pt_transverse_momentum->size());

        for (size_t j = 0; totalNumJets < numJetsLimit && j < mass_MeV->size(); j++, totalNumJets++) {
            auto element_mass = mass_MeV->at(j);
            masses_GeV.emplace_back(element_mass / 1000.0);

            auto element_phi = phi_azimuth_angle->at(j);
            y_phi_azimuthal_angles.emplace_back(element_phi);

            auto element_eta = eta_pseudorapidity->at(j);
            x_eta_pseudorapidities.emplace_back(element_eta);

            auto element_pt = pt_transverse_momentum->at(j);
            colors_transverse_momentum.emplace_back(element_pt / 1000.0);
        }
    }
    cout << "# Jets in input: " << x_eta_pseudorapidities.size() << endl;
    WriteDataSet(dataSetName, &x_eta_pseudorapidities, &y_phi_azimuthal_angles, &colors_transverse_momentum, &masses_GeV);
}