#include <cstdlib>
#include <string>
#include <iostream>

#include "TChain.h"
#include "TSystem.h"

#include "TChain.h"
#include "xAODRootAccess/Init.h"


#include "xAODTruthAna/TruthSelector.h"
#include "xAODTruthAna/ChainHelper.h"
#include "xAODTruthAna/string_utils.h"

using namespace std;


void help()
{
    cout << "Options:" << endl
         << "   -i|--input        name of input container"          << endl
         << "   -n|--num-events   default: -1 (all events)"         << endl
         << "   -f|--filelist     default: fileList.txt"            << endl
         << "   -d|--dbg-level    default: 0 (quiet)"               << endl
         << "   --outfile-name    default: 'superTruth.root'"       << endl
    << endl;
}

int main(int argc, char** argv)
{

    int nEvt        = -1;
    int dbg         = 0;
    string inputSampleName = "";
    string fileList = "fileList.txt";
    string outputFileName = "superTruth.root";

    cout << "runTruthSelector" << endl;
    cout << endl;

    int optin(1);
    while (optin < argc) {
        std::string opt = argv[optin];
        if      (opt=="-i" || opt=="--input")   { inputSampleName = argv[++optin]; }
        else if (opt=="-n" || opt=="--num-events") { nEvt = atoi(argv[++optin]); }
        else if (opt=="-f" || opt=="--filelist") { fileList = argv[++optin]; }
        else if (opt=="-d" || opt=="--dbg-level") { dbg = atoi(argv[++optin]); }
        else if (opt=="--outfile-name") { outputFileName = argv[++optin]; }
        else {
            cout << "Unknown command line argument : '" << opt << "'" << endl;
            help();
            return 1;
        }
        optin++;
    } // while

    if(inputSampleName == "") {
        cout << "You must provide an input sample name!" << endl;
        help();
        return 1;
    }

    xAOD::Init("TruthSelector");

    TChain* chain = new TChain("CollectionTree");
    int fileErr = ChainHelper::addInput(chain, fileList, dbg);
    if(fileErr) return 1;
    Long64_t nEntries = chain->GetEntries();
    chain->ls();

    // Build the TSelector
    TruthSelector* susyAna = new TruthSelector();
    susyAna->setDebug(dbg);
    susyAna->m_inputSampleName = inputSampleName;
    susyAna->m_outputFileName = outputFileName;

    if(nEvt<0) nEvt = nEntries;
    cout << endl;
    cout << "Total entries   : " << nEntries << endl;
    cout << "Process entries : " << nEvt << endl;
    chain->Process(susyAna, inputSampleName.c_str(), nEvt, 0); 

    cout << endl;
    cout << "TruthSelector job done" << endl;
    delete chain;

    return 0;
}

