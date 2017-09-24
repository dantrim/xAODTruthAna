#include "xAODTruthAna/WWbbTruth.h"

//std/stl
#include <iostream>
#include <cstdio> // printf
#include <sstream>
#include <math.h> // sqrt, cos
#include <fstream>
using namespace std;

//ROOT
#include "TH1F.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"

//xAOD/EDM
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/ReturnCheck.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODTruth/TruthEventContainer.h"

//ANA
#include "xAODTruthAna/xaod_utils.h"
#include "xAODTruthAna/MT2_ROOT.h"

std::string ANANAME = "WWbbTruth";


//////////////////////////////////////////////////////////////////////////////
WWbbTruth::WWbbTruth() :
    m_outfilename(""),
    m_outfile_setup(false),
    m_do_sumw_calc(false),
    m_total_sumw(0.0),
    m_dsid(0),
    m_weight(0.0),
    m_sumw(0.0),
    m_xsec(0.0),
    m_lumi(1.)
{
    cout << "WWbbTruth()" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::SlaveBegin(TTree* /*tree*/)
{
    timer()->Start();
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::initialize_xsec_map()
{
    xsec_map.clear();
    string filename = "xsec_wwbb.txt";
    std::ifstream dfs(filename.c_str());
    int d;
    double x;
    while(dfs >> d >> x) {
        xsec_map[d] = x;
    }

    cout << "------------------------------------------" << endl;
    cout << " XSEC LOADED " << endl;
    for(auto x : xsec_map) {
        cout << " > dsid " << x.first << "  xsec [pb] = " << x.second << endl;
    }
    cout << "------------------------------------------" << endl;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::initialize_sumw_map()
{
    sumw_map.clear();

    string filename = "sumw_wwbb.txt";
    std::ifstream dfs(filename.c_str());
    if(!dfs.good()) {
        cout << "ERROR sumw map file could not be found!" << endl;
        exit(1);
    }
    int d;
    double s;
    while(dfs >> d >> s) {
        sumw_map[d] = s;
    }
    cout << "------------------------------------------" << endl;
    cout << " SUMW LOADED " << endl;
    for(auto x : sumw_map) {
        cout << " > dsid " << x.first << "  sumw = " << x.second << endl;
    }
    cout << "------------------------------------------" << endl;

    
}
//////////////////////////////////////////////////////////////////////////////
Bool_t WWbbTruth::Process(Long64_t entry)
{
    static Long64_t chain_entry = -1;
    chain_entry++;
    event()->getEntry(chain_entry);

    const xAOD::EventInfo* ei = 0;
    RETURN_CHECK(GetName(), event()->retrieve(ei, "EventInfo") );

    if(!m_outfile_setup && !m_do_sumw_calc) {
        cout << "setting up output ROOT file" << endl;

        initialize_sumw_map();
        initialize_xsec_map();

        m_xsec = xsec_map[m_dsid];
        m_sumw = sumw_map[m_dsid];

        m_outfile_setup = true;
    }

    string fn = "WWbbTruth::Process    ";
    if(dbg() || chain_entry%5000==0) {
        cout << fn << " **** Processing entry " << setw(6) << chain_entry << " **** " << endl;
    }

    if(m_do_sumw_calc) {
        m_total_sumw += ei->mcEventWeight();
    }
    else {
        cout << fn << "Non-sumw calculation loop not yet set!" << endl;

        m_weight = ei->mcEventWeight();

        exit(1);
    }

    return true;
}
//////////////////////////////////////////////////////////////////////////////
void WWbbTruth::Terminate()
{
    cout << "WWbbTruth::Terminate" << endl;

    if(m_do_sumw_calc) {
        cout << "WWbbTruth::Terminate    TOTAL SUMW = " << m_total_sumw << endl;
    }
}


