#include <iostream>

#include "classes/HQClass.h"

using namespace std;

int analysis(){


    TChain *tree = new TChain("tree/HexaQAnalysis");
    tree->Add("../test/MC_tree.root");

    HQClass hqhand;
    hqhand.Init(tree);


    Long64_t nEntries = tree->GetEntries();    
    cout<<nEntries<<endl;

    for (int i=0; i<nEntries; i++){
        hqhand.GetEntry(i);

        cout<<hqhand.gen_mass->size()<<endl;
        for(int genp=0; genp<hqhand.gen_mass->size(); genp++){
            cout<<hqhand.gen_mass->at(genp)<<endl;
        }

        cout<<hqhand.nTrack<<endl;
        for(int trk=0; trk<hqhand.nTrack; trk++){
            cout<<hqhand.track_pt->at(trk)<<endl;
        }


    }

    return 0;

}
