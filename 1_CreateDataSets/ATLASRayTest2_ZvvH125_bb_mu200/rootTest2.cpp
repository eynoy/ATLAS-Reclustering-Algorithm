#include "../rootSaveTopoJetDataSet.h"

TChain in_chain("ntuple");
int rootTest2() {
    //in_chain.Add("TestDataATLAS/user.mmazza.23179734.STREAM_TREE._*.root");
    in_chain.Add("ZvvH125_bb_mu200_*.root");

    cout << endl << " -- Tree has " << in_chain.GetEntries() << " entries. -- " << endl << endl;

    cout << " -- Each entry has the following datapoints: -- " << endl << endl;

    cout << "       CaloCalTopoClusters_et  (vector<float>*)\n \
        CaloCalTopoClusters_eta (vector<float>*)\n \
        CaloCalTopoClusters_phi (vector<float>*)\n \
        CaloCalTopoClusters_m  (vector<float>*)\n \
        Calo420TopoClusters_et  (vector<float>*)\n \
        Calo420TopoClusters_eta  (vector<float>*)\n \
        Calo420TopoClusters_phi  (vector<float>*)\n \
        Calo420TopoClusters_m  (vector<float>*)\n \
        Calo422TopoClusters_et  (vector<float>*)\n \
        Calo422TopoClusters_eta  (vector<float>*)\n \
        Calo422TopoClusters_phi  (vector<float>*)\n \
        Calo422TopoClusters_m  (vector<float>*)\n \
        AntiKt4lctopoCaloCalJets_pt  (vector<float>*)\n \
        AntiKt4lctopoCaloCalJets_eta  (vector<float>*)\n \
        AntiKt4lctopoCaloCalJets_phi  (vector<float>*)\n \
        AntiKt4lctopoCaloCalJets_m  (vector<float>*)\n \
        AntiKt4lctopoCaloCalJets_nConstituents  (vector<int>*)\n \
        AntiKt4emtopoCalo420Jets_pt  (vector<float>*)\n \
        AntiKt4emtopoCalo420Jets_eta  (vector<float>*)\n \
        AntiKt4emtopoCalo420Jets_phi  (vector<float>*)\n \
        AntiKt4emtopoCalo420Jets_m  (vector<float>*)\n \
        AntiKt4emtopoCalo420Jets_nConstituents  (vector<int>*)\n \
        AntiKt4emtopoCalo422Jets_pt  (vector<float>*)\n \
        AntiKt4emtopoCalo422Jets_eta  (vector<float>*)\n \
        AntiKt4emtopoCalo422Jets_phi  (vector<float>*)\n \
        AntiKt4emtopoCalo422Jets_m  (vector<float>*)\n \
        AntiKt4emtopoCalo422Jets_nConstituents  (vector<int>*)\n \
        eventNumber (int?)\n \
        runNumber (int?)\n \
        weight (int?)\n \
        distFrontBunchTrain (int?)\n \
        BCID (int?)\n \
        averageInteractionsPerCrossing (int?)\n \
        cells_e  (vector<float>*)\n \
        cells_time       (vector<float>*)\n \
        cells_quality    (vector<unsigned int>*)\n \
        cells_provenance  (vector<unsigned int>*)\n" << endl;

    cout << " -- Please choose from the following datasets : -- " << endl << endl;

    cout << "\t1) CaloCalTopoClusters [\"standard topological clusters\"]" << endl;
    cout << "\t2) Calo420TopoClusters [\"alternative algorithm parameter clusters\"]" << endl;
    cout << "\t3) Calo422TopoClusters [\"alternative algorithm parameter clusters\"]" << endl << endl;
    cout << "\t4) AntiKt4lctopoCaloCalJets [\"standard anti-kT R=0.4 constructed jets from (1)\"]" << endl;
    cout << "\t5) AntiKt4emtopoCalo420Jets [\"anti-kT R=0.4 constructed jets from (2)\"]" << endl;
    cout << "\t6) AntiKt4emtopoCalo422Jets [\"anti-kT R=0.4 constructed jets from (3)\"]" << endl;

    cout << endl << "Please type the number you choose: ";

    int choice = 0; cin >> choice;

    std::string prefixes[6] = {"CaloCalTopoClusters",
                               "Calo420TopoClusters",
                               "Calo422TopoClusters",
                               "AntiKt4lctopoCaloCalJets",
                               "AntiKt4emtopoCalo420Jets",
                               "AntiKt4emtopoCalo422Jets"};
    std:string choice_prefix = prefixes[choice - 1];

    switch (choice) {
        case 1:
        case 2:
        case 3:
            cout << " -- Focusing on the Following Cluster Data Points: --" << endl;

            cout << "       " << choice_prefix << "_et  (vector<float>*)\n" <<
                choice_prefix << "_eta (vector<float>*)\n" <<
                choice_prefix << "_phi (vector<float>*)\n" <<
                choice_prefix << "_m  (vector<float>*)\n" << endl << endl;

            save_info_topoclusters(in_chain, choice_prefix);
            break;
        case 4:
        case 5:
        case 6:
            cout << " -- Focusing on the Following Jet Data Points: -- " << endl << endl;


            cout << "       " << choice_prefix << "_pt  (vector<float>*)\n" <<
                choice_prefix << "_eta (vector<float>*)\n" <<
                choice_prefix << "_phi (vector<float>*)\n" <<
                choice_prefix << "_m  (vector<float>*)\n" << endl << endl;

            save_info_jets(in_chain, choice_prefix);
            break;
        default:
            cout << "Not a valid number" << endl;
            break;
    }

    return 0;
}