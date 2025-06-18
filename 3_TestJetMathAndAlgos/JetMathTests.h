#include "../JetMath.h"

#include <chrono>

void WriteDataSet(string name, vector<double>* etas,
                     vector<double>* phis,
                     vector<double>* pts,
                     vector<double>* ms) {
    ofstream testdatafile;
    testdatafile.open(("../SavedDataSets/" + name + ".bin").c_str());
    
    int numMembers = etas->size();
    testdatafile.write(reinterpret_cast<const char*>(&numMembers), sizeof(int));
    testdatafile.write(reinterpret_cast<const char*>(etas->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(phis->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(pts->data()), numMembers*sizeof(double));
    testdatafile.write(reinterpret_cast<const char*>(ms->data()), numMembers*sizeof(double));
    testdatafile.close();
}

void TestRealJetAlgoOnDataSet(string dataset) {
    system(("fastjetExample/short-example ../SavedDataSets/" + dataset + ".bin").c_str());
    cout << endl;
}

void WriteTestDataSet() {
    vector<double> test_etas { -1.270057e+00, 1.281192e+00, 1.248796e+00, -1.204515e+00, -1.170003e+00, 1.136300e+00, 1.115080e+00, -1.252079e+00, -1.284597e+00, 1.148225e+00,  };
    vector<double> test_phis { -2.249272e+00, 1.029662e+00, 9.470263e-01, -2.162736e+00, -2.185710e+00, 8.551231e-01, 8.703419e-01, -2.191254e+00, 3.539034e-01, 9.599912e-01,  };
    vector<double> test_pts  { 1.020695e+06, 9.089080e+05, 5.859013e+05, 5.470122e+05, 2.649277e+05, 2.075110e+05, 1.835946e+05, 1.362937e+05, 6.147102e+04, 7.637645e+04,  };
    vector<double> test_ms   { 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,  };

    WriteDataSet("TestDataSet", &test_etas, &test_phis, &test_pts, &test_ms);
}

void WriteTestDataSet2() {
    vector<double> test_etas { 1.2355, 1.61391, 1.79862 };
    vector<double> test_phis { -2.80263, 2.85216, -2.84391 };
    vector<double> test_pts  { 20.2006, 22.9769, 16.2807 };
    vector<double> test_ms   { 4.79154, 4.20558, 2.93298 };

    WriteDataSet("TestDataSet2", &test_etas, &test_phis, &test_pts, &test_ms);
}

// Good test case for wrapping around phi coordinate
void WriteTestDataSet3() {
    vector<double> test_etas { 1.52282,  1.61391};
    vector<double> test_phis { -2.82105, 2.85216 };
    vector<double> test_pts  { 36.4736, 22.9769 };
    vector<double> test_ms   { 13.2785, 4.20558 };

    WriteDataSet("TestDataSet3", &test_etas, &test_phis, &test_pts, &test_ms);
}

void CompareAlgorithmResults(string dataset) {
    auto begin_time = std::chrono::steady_clock::now();
    TestRealJetAlgoOnDataSet(dataset);
    auto end_time = std::chrono::steady_clock::now();
    cout << "Took " <<  std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() / 1000.0 << " ms" << endl << endl;

    // Read data from Test File
    int numMembers = 0;

    ifstream testdatafile;
    testdatafile.open(("../SavedDataSets/" + dataset + ".bin").c_str());
    testdatafile.read(reinterpret_cast<char*>(&numMembers), sizeof(int));

    vector<double> test_etas(numMembers, 0);
    vector<double> test_phis(numMembers, 0);
    vector<double> test_pts(numMembers, 0);
    vector<double> test_ms(numMembers, 0);

    testdatafile.read(reinterpret_cast<char*>(test_etas.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_phis.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_pts.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_ms.data()), numMembers*sizeof(double));
    testdatafile.close();

    // double r_eta, r_phi, r_pt, r_m;
    // addFourMomenta(test_etas[0], test_phis[0], test_pts[0], test_ms[0], test_etas[1], test_phis[1], test_pts[1], test_ms[1], r_eta, r_phi, r_pt, r_m);
    // cout << r_eta << " | " << r_phi << " | " << r_pt << " | " << r_m << endl;

    const double R = 0.7;
    cout << endl << "IMPLEMENTED ALGORITHM:" << endl << endl;
    cout << "# Jets in input: " << numMembers << endl;
    cout << "Running anti-kT algorithm at R=" << R << endl;

    // for (int i = 0; i < numMembers; i++)
    //     cout << "Jet " << i << "'s pt: " << test_pts[i] << endl;

    vector<double>* result_jets_etas;
    vector<double>* result_jets_phis;
    vector<double>* result_jets_pts;
    vector<double>* result_jets_m;

    begin_time = std::chrono::steady_clock::now();
    antiKtAlgorithm(R,
                    &test_etas, &test_phis, &test_pts, &test_ms,
                    result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
    end_time = std::chrono::steady_clock::now();
    cout << "Took " <<  std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() / 1000.0 << " ms" << endl << endl;

    for (unsigned i = 0; i < result_jets_etas->size(); i++) {
        cout << "jet " << i << ": "<< result_jets_etas->at(i) << " " 
                   << result_jets_phis->at(i) << " " << result_jets_pts->at(i) << " " << result_jets_m->at(i) << endl; 
    }

    cout << endl << "IMPLEMENTED FAST ALGORITHM:" << endl << endl;
    cout << "Running fast anti-kT algorithm at R=" << R << endl;

    begin_time = std::chrono::steady_clock::now();
    antiKtAlgorithm_fast(R,
                    &test_etas, &test_phis, &test_pts, &test_ms,
                    result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
    end_time = std::chrono::steady_clock::now();
    cout << "Took " <<  std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() / 1000.0 << " ms" << endl << endl;

    for (unsigned i = 0; i < result_jets_etas->size(); i++) {
        cout << "jet " << i << ": "<< result_jets_etas->at(i) << " " 
                   << result_jets_phis->at(i) << " " << result_jets_pts->at(i) << " " << result_jets_m->at(i) << endl; 
    }
}

void CompareAlgorithmSpeeds(string dataset) {
    // Read data from Test File
    int numMembers = 0;

    ifstream testdatafile;
    testdatafile.open(("../SavedDataSets/" + dataset + ".bin").c_str());
    testdatafile.read(reinterpret_cast<char*>(&numMembers), sizeof(int));
    
    vector<double> test_etas(numMembers, 0);
    vector<double> test_phis(numMembers, 0);
    vector<double> test_pts(numMembers, 0);
    vector<double> test_ms(numMembers, 0);

    testdatafile.read(reinterpret_cast<char*>(test_etas.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_phis.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_pts.data()), numMembers*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(test_ms.data()), numMembers*sizeof(double));
    testdatafile.close();

    cout << numMembers << endl;

    const double R = 0.7;
    const int deltaI = 20;

    vector<double> speed_curve_points_x;
    vector<double> speed_curve_normalAlgo_points_y;
    vector<double> speed_curve_fastAlgo_points_y;
    for (int i = numMembers-1; i >= 0; i -= deltaI) {
        vector<double>* result_jets_etas;
        vector<double>* result_jets_phis;
        vector<double>* result_jets_pts;
        vector<double>* result_jets_m;

        auto begin_time = std::chrono::steady_clock::now();
        antiKtAlgorithm(R,
                        &test_etas, &test_phis, &test_pts, &test_ms,
                        result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
        auto end_time = std::chrono::steady_clock::now();
        double normalAlgo_deltaTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() / 1000.0;

        begin_time = std::chrono::steady_clock::now();
        antiKtAlgorithm_fast(R,
                        &test_etas, &test_phis, &test_pts, &test_ms,
                        result_jets_etas, result_jets_phis, result_jets_pts, result_jets_m);
        end_time = std::chrono::steady_clock::now();
        double fastAlgo_deltaTime = std::chrono::duration_cast<std::chrono::microseconds>(end_time - begin_time).count() / 1000.0;

        speed_curve_points_x.emplace_back(i);
        speed_curve_normalAlgo_points_y.emplace_back(normalAlgo_deltaTime / 1000.0);
        speed_curve_fastAlgo_points_y.emplace_back(fastAlgo_deltaTime  / 1000.0);

        for (int j = 0; j < deltaI; j++) {
            test_etas.erase(test_etas.begin() + (i - j));
            test_phis.erase(test_phis.begin() + (i - j));
            test_pts.erase(test_pts.begin() + (i - j));
            test_ms.erase(test_ms.begin() + (i - j));
        }

        cout << i << endl;
    }

    cout << speed_curve_points_x.size() << endl;

    TCanvas* speed_plot_canvas = new TCanvas("speed_plot_canvas", "Algorithm Speed Comparison (R=0.7)", 710, 10, 700, 500);
    TMultiGraph *mg = new TMultiGraph();

    TGraph* gr = new TGraph(speed_curve_points_x.size(), &speed_curve_points_x[0], &speed_curve_normalAlgo_points_y[0]);
    gr->SetTitle("Normal Algorithm");
    gr->SetMarkerColor(kSpring); // https://root.cern.ch/doc/master/classTColor.html
    gr->SetMarkerStyle(20);
    gr->GetXaxis()->SetRangeUser(0, numMembers);
    mg->Add(gr);

    TGraph* gr2 = new TGraph(speed_curve_points_x.size(), &speed_curve_points_x[0], &speed_curve_fastAlgo_points_y[0]);
    gr2->SetTitle("Optimized Algorithm");
    gr2->SetMarkerColor(kBlue); // https://root.cern.ch/doc/master/classTColor.html
    gr2->SetMarkerStyle(20);
    gr->GetXaxis()->SetRangeUser(0, numMembers);
    mg->Add(gr2);

    speed_plot_canvas->BuildLegend();

    mg->SetTitle("Algorithm Speed Comparison (R=0.7);Number of Input Jets;Time for Algorithm (s)");
    mg->GetXaxis()->SetRangeUser(0, numMembers);
    mg->Draw("APL");
}

// Library Testing
void JetMathTests() {
    //CompareAlgorithmResults("TestDataSet");

    //CompareAlgorithmResults("TestData60Jets");

    CompareAlgorithmResults("PositiveDataSet0");

    //CompareAlgorithmSpeeds("PositiveDataSet0");
}

