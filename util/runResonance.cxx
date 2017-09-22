#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

#include "TChain.h"
#include "TSystem.h"

#include "xAODRootAccess/Init.h"

#include "xAODTruthAna/TruthSelectorBase.h"
#include "xAODTruthAna/DiHiggsResonance.h"

#include "xAODTruthAna/ChainHelper.h"
#include "xAODTruthAna/string_utils.h"

void help()
{
    cout << "runAcceptance" << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "  -f                file list, file, or dir" << endl;
    cout << "  -n|--nevents      number of events to process" << endl;
    cout << "  -d                debug levent (int)" << endl;
    cout << "  -x                resonant x-mass" << endl;
    cout << "  -l|--lumi         integrated luminosity assumption (in fb-1) (default = 100 fb-1)" << endl;
    cout << "  -h                print this help message" << endl;
}

int main(int argc, char** argv)
{
    int nevents = -1;
    int dbg = 0;
    string filelist = "";
    int x_mass = -1;
    float lumi = 100.;

    int optin(1);
    while(optin < argc) {
        string opt = argv[optin];
        if          (opt == "-f") { filelist = argv[++optin]; }
        else if     (opt == "-n" || opt == "--nevents") { nevents = atoi(argv[++optin]); }
        else if     (opt == "-d") { dbg = atoi(argv[++optin]); }
        else if     (opt == "-x") { x_mass = atoi(argv[++optin]); }
        else if     (opt == "-l" || opt == "--lumi") { lumi = atof(argv[++optin]); }
        else if     (opt == "-h") { help(); return 0; }
        else {
            cout << "Unknown command line argument : '" << opt << "'" << endl;
            help();
            return 1;
        }
        optin++;
    } // while

    if(filelist == "") {
        cout << "You did not provide an input, exiting" << endl;
        return 1;
    }

    if(lumi < 0) {
        cout << "You provided a negative lumi, exiting" << endl;
        return 1;
    }

    cout << "DiHiggsResonance" << endl;
    xAOD::Init("DiHiggsResonance");
    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, filelist, dbg);
    if(file_err) return 1;
    Long64_t n_entries = chain->GetEntries();
    chain->ls();

    // build the selector
    res::DiHiggsResonance* ana = new res::DiHiggsResonance();
    ana->set_lumi(lumi);
    ana->set_debug(dbg);
    ana->set_input_samplename(filelist);
    ana->set_x_mass(x_mass);

    if(nevents < 0) nevents = n_entries;
    cout << "----------------------------------------" << endl;
    cout << " Total entries in input   : " << n_entries << endl;
    cout << " Total entries to process : " << nevents << endl;
    cout << "----------------------------------------" << endl;

    chain->Process(ana, filelist.c_str(), nevents);
    cout << "TruthSelector job done" << endl;
    delete chain;
    return 0;

}
