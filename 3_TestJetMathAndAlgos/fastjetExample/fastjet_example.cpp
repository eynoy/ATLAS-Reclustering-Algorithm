/*
 * g++ -std=c++20 fastjet_example.cpp -o short-example `fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`
 * ./short-example
 */

#include "fastjet/ClusterSequence.hh"
#include "../../JetMath.h"
#include <iostream>
#include <fstream>
using namespace fastjet;
using namespace std;

int main(int argc, char *argv[]) {

  if (argc == 1) {
    cout << "No file name given" << endl;
    return -1;
  }

  int numMembers = 0;

  ifstream testdatafile;
  testdatafile.open(argv[1]);
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

  vector<PseudoJet> particles;
  for (int i = 0; i < numMembers; i++) {
    double eta = test_etas[i];
    double phi = test_phis[i];
    double pt = test_pts[i];
    double m = test_ms[i];

    double px, py, pz, E;
    etaPhiPtToFourMomentum(eta, phi, pt, m, px, py, pz, E);

    particles.push_back(PseudoJet(px, py, pz, E));
  }
  
  /*
  cout << "\nCurrently storing " << particles.size() << " particles. They are: " << endl;
  for (int i = 0; i < particles.size(); i++) {
    cout << "\tParticle " << i << ": { eta = " << particles.at(i).eta() << \
                                    ", phi = " << particles.at(i).phi_std() << \
                                    ", pt = " << particles.at(i).pt() << \
                                    ", m = " << particles.at(i).m() << \
                                    ", Et = " << particles.at(i).Et() << \
                                    ", E = " << particles.at(i).E() <<" }" << endl;
  }
  cout << endl;

  PseudoJet sumParticle = particles.at(0) + particles.at(1);
  cout << "\tParticle Sum is: { eta = " << sumParticle.eta() << \
                                    ", phi = " << sumParticle.phi_std() << \
                                    ", pt = " << sumParticle.pt() << \
                                    ", m = " << sumParticle.m() << \
                                    ", E = " << sumParticle.E() << " }" << endl;
  */

  //choose a jet definition
  double R = 0.7;
  JetDefinition jet_def(antikt_algorithm, R, E_scheme, N2Plain);

  // run the clustering, extract the jets
  ClusterSequence cs(particles, jet_def);
  vector<PseudoJet> jets = sorted_by_pt(cs.inclusive_jets());

  // print out some infos
  cout << "FASTJET ALGORITHM:" << endl << endl;
  cout << "# Jets in input: " << numMembers << endl;
  cout << "Clustering with " << jet_def.description() << ": " << endl;

  // print the jets
  cout <<   "        eta phi pt m" << endl;
  for (unsigned i = 0; i < jets.size(); i++) {
    cout << "jet " << i << ": "<< jets[i].eta() << " " 
                   << jets[i].phi_std() << " " << jets[i].pt() << " " << jets[i].m() << endl;
    // vector<PseudoJet> constituents = jets[i].constituents();
    // for (unsigned j = 0; j < constituents.size(); j++) {
    //   cout << "    constituent " << j << ": " << constituents[j].eta() << " " 
    //                << constituents[j].phi_std() << " " << constituents[j].pt() << " " << constituents[j].m() << " w/ dij=" << (j == 1 ? constituents[j].delta_R(constituents[j-1]) : 0.0) << endl;
    // }
  }
} 