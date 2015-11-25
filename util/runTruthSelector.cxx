#include <cstdlib>
#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Algorithm.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>
#include "SampleHandler/ScanDir.h"
#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoop/OutputStream.h>


#include "xAODTruthAna/TruthSelector.h"
#include "xAODTruthAna/utils.h"

#include <iostream>
#include <string>
using namespace std;


void help() {
    cout << "Usage: "
         <<"runTruthSelector [options]"<<endl
         <<"-o --output        "       <<"\t"<<" output directory"<<endl
         <<"-f --filelist      "       <<"\t"<<" input filelist"  <<endl
         <<"-i --input-dir     "       <<"\t"<<" input directory" << endl
         <<"        to be scanned if you do not provide a filelist" << endl
         <<"-p --file-pattern  "       <<"\t"<<" pattern used to scan directory"<<endl
         <<"-s --sample-name   "       <<"\t"<<" sample name (used for output)"<<endl
         <<"-t --tag           "       <<"\t"<<" tag used in outputname if submitting to the grid"<<endl
         <<"-n --num-events    "       <<"\t"<<" number of events to process"<<endl
         <<"-v --verbose       "       <<"\t"<<" print more"<<endl
         <<"-g --grid          "       <<"\t"<<" submit jobs to the grid with prun"<<endl
         <<endl;
}

int main(int argc, char* argv[]) {

    string sampleName;
    string outputDir;
    string inputFileList;
    string inputDirectory;
    string inputFilePattern = "*.pool.root*"; // used when scanning input directory
    string gridTag; // if submitted to grid use this in the output name for tagging
    int numEvents=-1;
    bool verbose = false;
    bool runGrid = false;

    for(int i = 1; i < argc; i++) {
        string opt = argv[i];
        if     (opt=="-s" || opt=="--sample-name") sampleName = argv[++i];
        else if(opt=="-o" || opt=="--output") outputDir = argv[++i];
        else if(opt=="-f" || opt=="--filelist") inputFileList = argv[++i];
        else if(opt=="-i" || opt=="--input-dir") inputDirectory = argv[++i];
        else if(opt=="-p" || opt=="--file-pattern") inputFilePattern = argv[++i];
        else if(opt=="-t" || opt=="--tag") gridTag = argv[++i];
        else if(opt=="-n" || opt=="--num-events") numEvents = atoi(argv[++i]);
        else if(opt=="-v" || opt=="--verbose") verbose = true;
        else if(opt=="-g" || opt=="--grid") runGrid = true;
        else {
            cout << "unknown option '"<<opt<<"'"<<endl;
            help();
            return 1;
        }
    } // i


    xAOD::Init().ignore();  // setup the job for xAOD access

    SH::SampleHandler sh;

    const char* inputFilePath = gSystem->ExpandPathName(inputDirectory.c_str());
    if(not runGrid) {
        SH::ScanDir().sampleDepth(0).samplePattern(inputFilePattern).scan(sh, inputFilePath);
    } else {
        SH::scanDQ2(sh, inputFilePath);
    }

   // if(not sampleName.size()) {
   //     cout << "required sample name"<<endl;
   //     return 1;
   // }

    //if(not outputDir.size()) {
    //    cout << "required output directory" << endl;
    //    return 1;
    //}
    //mkdirIfNeeded(basedir(outputDir));

    //if(inputFileList.size())
    //    SH::readFileList(sh, sampleName, inputFileList);
    //else if(inputDirectory.size() && inputFilePattern.size())
    //    SH::ScanDir().sampleDepth(0).samplePattern(inputFilePattern).scan(sh, inputDirectory);
    //else {
    //    cout <<"required input (filelist or directory)"<<endl;
    //    return 1;
    //}

    sh.setMetaString("nc_tree", "CollectionTree");
    if(verbose)
        sh.print();

    EL::Job job;
    job.sampleHandler(sh);
    if(numEvents>0)
        job.options()->setDouble(EL::Job::optMaxEvents, numEvents);
//    EL::OutputStream output ("myOutput");
//    job.outputAdd (output);
//    EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
//    job.algsAdd (ntuple);

    TruthSelector* alg = new TruthSelector();
    //alg->verbose = verbose;
    job.algsAdd(alg);

    cout << "outputDir " << outputDir <<endl;

 //   if(!runGrid) {
        EL::DirectDriver driver;
        driver.submit( job, outputDir );
 //   } else {
 //       EL::PrunDriver driver;
 //       driver.options()->setString("nc_outputSampleName", "user.dantrim."+gridTag+".%in:name[2]%.%in:name[6]%");
 //       driver.submitOnly(job, outputDir);
 //   }

    return 0;

}
