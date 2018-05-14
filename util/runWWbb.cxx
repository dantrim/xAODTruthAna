#include <cstdlib>
#include <string>
#include <iostream>
using namespace std;

#include "TChain.h"
#include "TSystem.h"

#include "xAODRootAccess/Init.h"

#include "xAODTruthAna/TruthSelectorBase.h"
#include "xAODTruthAna/WWbbTruth.h"

#include "xAODTruthAna/ChainHelper.h"
#include "xAODTruthAna/string_utils.h"

void help()
{
    cout << "runWWbb" << endl;
    cout << endl;
    cout << " Options:" << endl;
    cout << "  -f                   file list, file, or dir" << endl;
    cout << "  -n|--nevents         number of events to process" << endl;
    cout << "  -d                   debug level (int)" << endl;
    cout << "  --dsid               sample dsid [REQUIRED]" << endl;
    cout << "  --sumw               do sumw calc only" << endl;
    cout << "  -l|--lumi            set lumi (in pb) [default=1000 pb]" << endl;
    cout << "  -s|--suffix          add a suffix to output root file" << endl;
    cout << "  --hh                 set for running over hh signals" << endl;
    cout << "  --top-sample         calculate top variables" << endl;
    cout << "  --bjet-eff           emulate b-tagging efficiency" << endl;
    cout << "  --skip-maps          do not use xsec or sumw map files" << endl;
    cout << "  -h                   print this help message" << endl;
}

int main(int argc, char** argv)
{
    int nevents = -1;
    string filelist = "";
    int dbg = 0;
    bool do_sumw = false;
    double lumi = 1.0; /* fb-1 */
    int dsid = 0;
    string suffix = "";
    bool is_hh_signal = false;
    bool use_bjet_eff = false;
    bool top_sample = false;
    bool skip_maps = false;

    int optin(1);
    while(optin < argc) {
        string opt = argv[optin];
        if          (opt == "-f") { filelist = argv[++optin]; }
        else if     (opt == "-n" || opt == "--nevents") { nevents = atoi(argv[++optin]); }
        else if     (opt == "-d") { dbg = atoi(argv[++optin]); }
        else if     (opt == "--sumw") { do_sumw = true; }
        else if     (opt == "-l" || opt == "--lumi") { lumi = atof(argv[++optin]); }
        else if     (opt == "-h") { help(); return 0; }
        else if     (opt == "-s" || opt == "--suffix") { suffix = argv[++optin]; }
        else if     (opt == "--dsid") { dsid = atoi(argv[++optin]); }
        else if     (opt == "--bjet-eff") { use_bjet_eff = true; }
        else if     (opt == "--hh") { is_hh_signal = true; }
        else if     (opt == "--top-sample") { top_sample = true; }
        else if     (opt == "--skip-maps") { skip_maps = true; }
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

    if(dsid == 0) {
        cout << "You did not provide an input sample dsid, exiting" << endl;
        return 1;
    }

    cout << "WWbb" << endl;
    xAOD::Init("WWbb");
    TChain* chain = new TChain("CollectionTree");
    int file_err = ChainHelper::addInput(chain, filelist, true);
    if(file_err) return 1;
    Long64_t n_entries = chain->GetEntries();
    chain->ls();

    // build the selector
    WWbbTruth* ana = new WWbbTruth();
    ana->set_debug(dbg);
    ana->set_input_samplename(filelist);
    ana->do_sumw_calc(do_sumw);
    ana->set_lumi(lumi);
    ana->set_dsid(dsid);
    ana->set_suffix(suffix);
    ana->set_use_bjet_eff(use_bjet_eff);
    ana->set_is_hh_signal(is_hh_signal);
    ana->set_skip_maps(skip_maps);
    ana->set_top_sample(top_sample);

    if(nevents < 0) nevents = n_entries;
    cout << "-----------------------------------------------" << endl;
    cout << " Total entries in input    : " << n_entries << endl;
    cout << " Total entries to process  : " << nevents << endl;
    cout << "-----------------------------------------------" << endl;

    chain->Process(ana, filelist.c_str(), nevents);
    cout << "WWbb job done" << endl;
    delete chain;
    return 0;
}
