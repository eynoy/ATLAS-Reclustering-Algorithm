#include "../JetMath.h"

#include <TApplication.h>
#include <TGraph.h>
#include <TColor.h>

bool AcceptanceAlgorithm(double eta, double phi, double pt, double m, double pt_cutoff) {
    return pt >= pt_cutoff;
}

const double dpt_cutoff = 1.0;
void GenerateRocCurve(double R) {
    // Positive Data Set: 1116 Jets AntiKt4lctopoCaloCalJets from Zprime_tthad_mu200
    int numMembersPositive = 0;

    ifstream testdatafile;
    testdatafile.open("../SavedDataSets/PositiveDataSet.bin");
    testdatafile.read(reinterpret_cast<char*>(&numMembersPositive), sizeof(int));
    
    vector<double> p_test_etas(numMembersPositive, 0);
    vector<double> p_test_phis(numMembersPositive, 0);
    vector<double> p_test_pts(numMembersPositive, 0);
    vector<double> p_test_ms(numMembersPositive, 0);

    testdatafile.read(reinterpret_cast<char*>(p_test_etas.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_phis.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_pts.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_ms.data()), numMembersPositive*sizeof(double));
    testdatafile.close();

    vector<double>* positive_result_jets_etas;
    vector<double>* positive_result_jets_phis;
    vector<double>* positive_result_jets_pts;
    vector<double>* positive_result_jets_m;
    antiKtAlgorithm_fast(R,
                    &p_test_etas, &p_test_phis, &p_test_pts, &p_test_ms,
                    positive_result_jets_etas, positive_result_jets_phis, positive_result_jets_pts, positive_result_jets_m);


    // Negative Data Set: 1116 Jets AntiKt4lctopoCaloCalJets from JZ3WithSW_mu200
    int numMembersNegative = 0;

    testdatafile.open("../SavedDataSets/NegativeDataSet.bin");
    testdatafile.read(reinterpret_cast<char*>(&numMembersNegative), sizeof(int));
    
    vector<double> n_test_etas(numMembersNegative, 0);
    vector<double> n_test_phis(numMembersNegative, 0);
    vector<double> n_test_pts(numMembersNegative, 0);
    vector<double> n_test_ms(numMembersNegative, 0);

    testdatafile.read(reinterpret_cast<char*>(n_test_etas.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_phis.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_pts.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_ms.data()), numMembersNegative*sizeof(double));
    testdatafile.close();

    vector<double>* negative_result_jets_etas;
    vector<double>* negative_result_jets_phis;
    vector<double>* negative_result_jets_pts;
    vector<double>* negative_result_jets_m;
    antiKtAlgorithm_fast(R,
                    &n_test_etas, &n_test_phis, &n_test_pts, &n_test_ms,
                    negative_result_jets_etas, negative_result_jets_phis, negative_result_jets_pts, negative_result_jets_m);

    int numClusteredPositive = positive_result_jets_etas->size();
    int numClusteredNegative = negative_result_jets_etas->size();

    cout << "pos: " << numMembersPositive << " -> " << numClusteredPositive << endl;
    cout << "neg: " << numMembersNegative << " -> " << numClusteredNegative << endl;

    vector<double> roc_curve_points_x;
    vector<double> roc_curve_points_y;

    vector<double> efficiency_curve_points_x;
    vector<double> efficiency_curve_points_y;
    for (double pt_cutoff = dpt_cutoff; pt_cutoff < 5000.0; pt_cutoff += dpt_cutoff) {

        int numTruePositive = 0;
        for (int i = 0; i < numClusteredPositive; i++) {
            if (AcceptanceAlgorithm(positive_result_jets_etas->at(i),
                                    positive_result_jets_phis->at(i),
                                    positive_result_jets_pts->at(i),
                                    positive_result_jets_m->at(i),
                                    pt_cutoff)) {
                numTruePositive++;
            }
        }
        double truePositiveRate = ((double)numTruePositive)/((double)numClusteredPositive);
        roc_curve_points_x.emplace_back(truePositiveRate);

        efficiency_curve_points_x.emplace_back(pt_cutoff);
        efficiency_curve_points_y.emplace_back(truePositiveRate);

        int numFalsePositive = 0;
        for (int i = 0; i < numClusteredNegative; i++) {
            if (AcceptanceAlgorithm(negative_result_jets_etas->at(i),
                                    negative_result_jets_phis->at(i),
                                    negative_result_jets_pts->at(i),
                                    negative_result_jets_m->at(i),
                                    pt_cutoff)) {
                numFalsePositive++;
            }
        }
        double falsePositiveRate = ((double)numFalsePositive)/((double)numClusteredNegative);
        double backgroundRejection = 1.0 / falsePositiveRate;

        roc_curve_points_y.emplace_back(backgroundRejection);
    }

    // Create ROC Cruve Graph
    TGraph* gr = new TGraph(roc_curve_points_x.size(), &roc_curve_points_x[0], &roc_curve_points_y[0]);
    gr->SetTitle("ROC Curve;True Positive Rate;Background Rejection (1/FPR)");
    
    TCanvas* roc_plot_canvas = new TCanvas("roc_plot_canvas", "ROC Curve", 710, 10, 700, 500);
    
    gr->SetMarkerColor(kBlue);
    gr->SetMarkerStyle(20);

    gr->GetYaxis()->SetRangeUser(0, 10);

    gr->Draw("APL");

    // Create Eff Cruve Graph
    TGraph* gr2 = new TGraph(efficiency_curve_points_x.size(), &efficiency_curve_points_x[0], &efficiency_curve_points_y[0]);
    gr2->SetTitle("Efficiency Curve;P_t Cutoff;True Positive Rate");
    
    TCanvas* eff_plot_canvas = new TCanvas("eff_plot_canvas", "Efficiency Curve", 0, 10, 700, 500);
    
    gr2->SetMarkerColor(kBlue);
    gr2->SetMarkerStyle(20);

    gr2->GetYaxis()->SetRangeUser(0, 1);

    gr2->Draw("APL");
}

int colorPalleteRocCurves[6] = {8, 9, 20, 46, 5, 6};
void GenerateRocCurves(double R_min, double R_max, double dR) {
    // Positive Data Set: 1116 Jets AntiKt4lctopoCaloCalJets from Zprime_tthad_mu200
    int numMembersPositive = 0;

    ifstream testdatafile;
    testdatafile.open("../SavedDataSets/PositiveDataSet.bin");
    testdatafile.read(reinterpret_cast<char*>(&numMembersPositive), sizeof(int));
    cout << numMembersPositive << " postive data set jets." << endl;
    
    vector<double> p_test_etas(numMembersPositive, 0);
    vector<double> p_test_phis(numMembersPositive, 0);
    vector<double> p_test_pts(numMembersPositive, 0);
    vector<double> p_test_ms(numMembersPositive, 0);

    testdatafile.read(reinterpret_cast<char*>(p_test_etas.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_phis.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_pts.data()), numMembersPositive*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(p_test_ms.data()), numMembersPositive*sizeof(double));
    testdatafile.close();

    vector<double>* positive_result_jets_etas;
    vector<double>* positive_result_jets_phis;
    vector<double>* positive_result_jets_pts;
    vector<double>* positive_result_jets_m;

    // Negative Data Set: 1116 Jets AntiKt4lctopoCaloCalJets from JZ3WithSW_mu200
    int numMembersNegative = 0;

    testdatafile.open("../SavedDataSets/NegativeDataSet.bin");
    testdatafile.read(reinterpret_cast<char*>(&numMembersNegative), sizeof(int));
    cout << numMembersNegative << " negative data set jets." << endl;

    vector<double> n_test_etas(numMembersNegative, 0);
    vector<double> n_test_phis(numMembersNegative, 0);
    vector<double> n_test_pts(numMembersNegative, 0);
    vector<double> n_test_ms(numMembersNegative, 0);

    testdatafile.read(reinterpret_cast<char*>(n_test_etas.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_phis.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_pts.data()), numMembersNegative*sizeof(double));
    testdatafile.read(reinterpret_cast<char*>(n_test_ms.data()), numMembersNegative*sizeof(double));
    testdatafile.close();

    vector<double>* negative_result_jets_etas;
    vector<double>* negative_result_jets_phis;
    vector<double>* negative_result_jets_pts;
    vector<double>* negative_result_jets_m;

    TCanvas* roc_plot_canvas = new TCanvas("roc_plot_canvas", "ROC Curve", 710, 10, 700, 500); int numColor = 0;
    TMultiGraph *mg = new TMultiGraph();
    for (double R = R_min; R <= R_max; R += dR, numColor++) {
        cout << "\tPass R=" << R << endl;
        antiKtAlgorithm_fast(R,
                        &p_test_etas, &p_test_phis, &p_test_pts, &p_test_ms,
                        positive_result_jets_etas, positive_result_jets_phis, positive_result_jets_pts, positive_result_jets_m); cout << "\t\tPassed Pos Antikt" << endl;

        antiKtAlgorithm_fast(R,
                        &n_test_etas, &n_test_phis, &n_test_pts, &n_test_ms,
                        negative_result_jets_etas, negative_result_jets_phis, negative_result_jets_pts, negative_result_jets_m); cout << "\t\tPassed Neg Antikt" << endl;

        int numClusteredPositive = positive_result_jets_etas->size();
        int numClusteredNegative = negative_result_jets_etas->size();

        cout << "\t\tBefore loop" << endl;
        vector<double> roc_curve_points_x;
        vector<double> roc_curve_points_y;
        for (double pt_cutoff = dpt_cutoff; pt_cutoff < 5000.0; pt_cutoff += dpt_cutoff) {

            int numTruePositive = 0;
            for (int i = 0; i < numClusteredPositive; i++) {
                if (AcceptanceAlgorithm(positive_result_jets_etas->at(i),
                                        positive_result_jets_phis->at(i),
                                        positive_result_jets_pts->at(i),
                                        positive_result_jets_m->at(i),
                                        pt_cutoff)) {
                    numTruePositive++;
                }
            }
            double truePositiveRate = ((double)numTruePositive)/((double)numClusteredPositive);
            roc_curve_points_x.emplace_back(truePositiveRate);

            int numFalsePositive = 0;
            for (int i = 0; i < numClusteredNegative; i++) {
                if (AcceptanceAlgorithm(negative_result_jets_etas->at(i),
                                        negative_result_jets_phis->at(i),
                                        negative_result_jets_pts->at(i),
                                        negative_result_jets_m->at(i),
                                        pt_cutoff)) {
                    numFalsePositive++;
                }
            }
            double falsePositiveRate = ((double)numFalsePositive)/((double)numClusteredNegative);
            double backgroundRejection = 1.0 / falsePositiveRate;

            roc_curve_points_y.emplace_back(backgroundRejection);
        }

        cout << "Begin gen ROC Curves" << endl;

        // Create ROC Curve Graph
        TGraph* gr = new TGraph(roc_curve_points_x.size(), &roc_curve_points_x[0], &roc_curve_points_y[0]);
        gr->SetTitle(("R = " + std::to_string(R).substr(0,4)).c_str());
        
        gr->SetMarkerColor(colorPalleteRocCurves[numColor]); // https://root.cern.ch/doc/master/classTColor.html
        gr->SetMarkerStyle(20);

        gr->GetYaxis()->SetRangeUser(0, 10);

        mg->Add(gr);
    }

    mg->SetTitle("ROC Curve;True Positive Rate;Background Rejection (1/FPR)");

    mg->GetYaxis()->SetRangeUser(0, 10);

    mg->Draw("APL");

    roc_plot_canvas->BuildLegend();
}

void JetRocCurves() {
    //GenerateRocCurve(0.7);
    GenerateRocCurves(0.25, 1.00, 0.25);
}

