#include <string>
#include <iostream>
#include <filesystem>
#include <iostream>
#include <vector>
#include "TChain.h"

#include "../JetMath.h"

void draw_jet_scatter_plot(const std::string datasetname,
                            std::vector<double>* x_eta_pseudorapidities,
                            std::vector<double>* y_phi_azimuthal_angles,
                            std::vector<double>* colors_transverse_momentum,
                            std::vector<double>* masses_GeV) {
    gStyle->SetPalette(kBird, 0, 0.6);

    TCanvas* angle_plot_canvas = new TCanvas(("angle_plot_canvas_" + datasetname).c_str(), ("Angle Locations in " + datasetname).c_str(), 710, 10, 700, 500);

    std::cout << x_eta_pseudorapidities->size() << std::endl;

    TScatter* angle_plot = new TScatter(x_eta_pseudorapidities->size(), &(*x_eta_pseudorapidities)[0], &(*y_phi_azimuthal_angles)[0], &(*colors_transverse_momentum)[0], &(*colors_transverse_momentum)[0]);
    angle_plot->SetTitle(("Angle Locations of Jets in " + datasetname + ";#eta [pseudorapidity];#phi [azimuthal angle] (rad);t (GeV)").c_str());
    angle_plot->SetMarkerStyle(20);
    angle_plot->Draw("AP");
}

void draw_jet_scatter_plot_3d(const std::string datasetname,
                            std::vector<double>* x_eta_pseudorapidities,
                            std::vector<double>* y_phi_azimuthal_angles,
                            std::vector<double>* colors_transverse_momentum,
                            std::vector<double>* masses_GeV) {

    const double R = 1;     // radius

    TCanvas *c1 = new TCanvas(("c1_" + datasetname).c_str(), ("3D Cylinder Scatter of " + datasetname).c_str(), 800, 600);
    TPolyMarker3D *pm3d = new TPolyMarker3D(masses_GeV->size(), 20);  // size 20 marker

    gStyle->SetPalette(kBird, 0, 0.6);
    pm3d->SetMarkerStyle(20);  // full dot

    double maxMomentum = *std::max_element(colors_transverse_momentum->begin(), colors_transverse_momentum->end());
    for (int i = 0; i < masses_GeV->size(); ++i) {
        double phi = y_phi_azimuthal_angles->at(i);
        double theta = 2*atan(exp(-x_eta_pseudorapidities->at(i)));
        
        double normMomentum = colors_transverse_momentum->at(i) / maxMomentum; // 0 to 1
        int colorIndex = TColor::GetColorPalette(static_cast<int>(normMomentum * 99));  // scale to palette index 0–99

        double z = -log(normMomentum) * R * cos(theta); //x_eta_pseudorapidities->at(i)
        double x = -log(normMomentum) * R * sin(phi) * cos(theta);// R * sin(phi);
        double y = -log(normMomentum) * R * sin(phi) * sin(theta);

        TPolyMarker3D *pm = new TPolyMarker3D(1);
        pm->SetPoint(i, x, y, z);
        pm->SetMarkerStyle(20);
        pm->SetMarkerSize(1);
        pm->SetMarkerColor(colorIndex);
        pm->Draw("same");
    }

    gPad->Update();
    gPad->GetView()->ShowAxis();
}

void draw_jet_mass_histogram(const std::string datasetname, std::vector<double>* masses_GeV) {
    TH1F* mass_historgram = new TH1F(("mass_dist_" + datasetname).c_str(),
                                 "Mass Distributions;Mass (GeV);# Occurences",
                                 100, // Number of Bins
                                 0.0, // Lower X Boundary
            /*max is 500.0*/     50.0); // Upper X Boundary

    mass_historgram->SetLineWidth(3);
    mass_historgram->SetLineColor(kBlue);

    mass_historgram->FillN(masses_GeV->size(), &(*masses_GeV)[0], nullptr);
    
    TCanvas* histogram_canvas = new TCanvas(("histogram_canvas_" + datasetname).c_str(), ("Mass Distributions in " + datasetname).c_str());
    mass_historgram->Draw();
}

void visualizeJets() {
    // Collect User Input on which Data Set to Visualize
    cout << " -- Present Datasets -- " << endl << endl;

    int i = 0;
    std::vector<std::string> datasetNames;
    for (const auto & entry : std::filesystem::directory_iterator("../SavedDataSets/")) {
        std::string fileName = entry.path().stem();
        if (fileName[0] == '.') continue;

        datasetNames.push_back(fileName);
        std::cout << "\t" << i++ << ") " << fileName << std::endl;
    }

    std::cout << std::endl << "Choose Dataset: ";
    int datasetNameIndex = 0; cin >> datasetNameIndex;
    std::string datasetName = datasetNames[datasetNameIndex];

    // Collect User Input on which R size to use in Anti-kt Algorithm
    std::cout << std::endl << "What R should be used in the Anti-kt Algorithm? ";
    double R = 0; cin >> R;

    // Load Data Set
    int numMembersPositive = 0;

    ifstream testdatafile;
    testdatafile.open("../SavedDataSets/" + datasetName + ".bin");
    testdatafile.read(reinterpret_cast<char*>(&numMembersPositive), sizeof(int));
    
    vector<double> test_etas(numMembersPositive, 0);
    vector<double> test_phis(numMembersPositive, 0);
    vector<double> test_pts(numMembersPositive, 0);
    vector<double> test_ms(numMembersPositive, 0);

    testdatafile.read(reinterpret_cast<char*>(test_etas.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_phis.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_pts.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_ms.data()), numMembersPositive*sizeof(double));
    testdatafile.close();

    vector<double>* result_jets_etas;
    vector<double>* result_jets_phis;
    vector<double>* result_jets_pts;
    vector<double>* result_jets_m;
    antiKtAlgorithm_fast(R,
                    &test_etas, &test_phis, &test_pts, &test_ms,
                    result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);


    // Collect User Input on which Data Set Visualions to Use
    cout << " -- Choose the Following Visualizations -- " << endl << endl;

    cout << "\t0) Mass Histogram" << endl;
    cout << "\t1) 2D Scatter Plot" << endl;
    cout << "\t2) 3D Scatter Plot" << endl;

    cout << "\t3) After Anti-kt Mass Histogram" << endl;
    cout << "\t4) After Anti-kt 2D Scatter Plot" << endl;
    cout << "\t5) After Anti-kt 3D Scatter Plot" << endl;

    std::cout << std::endl << "Choose Visulation(s), comma separated without spaces: ";
    std::string inputNumStr; std::cin >> inputNumStr;

    std::istringstream ss(inputNumStr);
    std::string token;
    while(std::getline(ss, token, ',')) {
        int chosenVisulaization = std::stoi(token);

        switch (chosenVisulaization) {
            case 0:
                draw_jet_mass_histogram(datasetName, &test_ms);
                break;
            case 1:
                draw_jet_scatter_plot(datasetName, &test_etas, &test_phis, &test_pts, &test_ms);
                break;
            case 2:
                draw_jet_scatter_plot_3d(datasetName, &test_etas, &test_phis, &test_pts, &test_ms);
                break;
            case 3:
                draw_jet_mass_histogram(datasetName + " after the Anti-kt Algorithm", result_jets_m);
                break;
            case 4:
                draw_jet_scatter_plot(datasetName + " after the Anti-kt Algorithm", result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
                break;
            case 5:
                draw_jet_scatter_plot_3d(datasetName + " after the Anti-kt Algorithm", result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
                break;
        }
    }
}

