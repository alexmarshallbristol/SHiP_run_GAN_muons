#include <math.h>
#include "TSystem.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TFile.h"
#include "FairPrimaryGenerator.h"
#include "MuonBackGenerator.h"
#include "TDatabasePDG.h"               // for TDatabasePDG
#include "TMath.h"                      // for Sqrt
#include "vetoPoint.h"
#include "ShipMCTrack.h"
#include "TMCProcess.h"
#include <algorithm>
#include <fstream>

// read events from Pythia8/Geant4 base simulation (only target + hadron absorber

// -----   Default constructor   -------------------------------------------
MuonBackGenerator::MuonBackGenerator() {}
// -------------------------------------------------------------------------
// -----   Default constructor   -------------------------------------------
Bool_t MuonBackGenerator::Init(const char* fileName) {
  return Init(fileName, 0, false);
}
// -----   Default constructor   -------------------------------------------
Bool_t MuonBackGenerator::Init(const char* fileName, const int firstEvent, const Bool_t fl = false ) {
  fLogger = FairLogger::GetLogger();
  fLogger->Info(MESSAGE_ORIGIN,"Opening input file %s",fileName);
  if (0 == strncmp("/eos",fileName,4) ) {
     TString tmp = gSystem->Getenv("EOSSHIP");
     tmp+=fileName;
     fInputFile  = TFile::Open(tmp);
  }else{
  fInputFile  = new TFile(fileName);
  }
  if (fInputFile->IsZombie()) {
    fLogger->Fatal(MESSAGE_ORIGIN, "Error opening the Signal file:",fInputFile);
  }
  fn = firstEvent;
  fPhiRandomize = fl;
  fSameSeed = 0;
  fsmearBeam = 0; // default no beam smearing, use SetSmearBeam(sb) if different, sb [cm]

   fTree = (TTree *)fInputFile->Get("pythia8-Geant4");
   fNevents   = fTree->GetEntries();

   std::cout << "fNevents " << fNevents <<"\n";
   Double_t fUniquieID,fPdgCode,fMotherId,fPx,fPy,fPz,fM,fStartX,fStartY,fStartZ,fW,fProcID;

  fTree->SetBranchAddress("pythiaid",&fPdgCode);
  fTree->SetBranchAddress("mother_id",&fMotherId);
  fTree->SetBranchAddress("px",&fPx);
  fTree->SetBranchAddress("py",&fPy);
  fTree->SetBranchAddress("pz",&fPz);
  fTree->SetBranchAddress("x",&fStartX);
  fTree->SetBranchAddress("y",&fStartY);
  fTree->SetBranchAddress("z",&fStartZ);
  fTree->SetBranchAddress("w",&fW);
  fTree->SetBranchAddress("process_id",&fProcID);


  return kTRUE;
}
// -----   Destructor   ----------------------------------------------------
MuonBackGenerator::~MuonBackGenerator()
{
}
// -------------------------------------------------------------------------

// -----   Passing the event   ---------------------------------------------
Bool_t MuonBackGenerator::ReadEvent(FairPrimaryGenerator* cpg)
{
TDatabasePDG* pdgBase = TDatabasePDG::Instance();

Double_t mass,e,tof,phi;
Double_t dx = 0, dy = 0;

std::vector<int> muList;

//define variables to pull from tree
Double_t fUniquieID,fPdgCode,fMotherId,fPx,fPy,fPz,fM,fStartX,fStartY,fStartZ,fW,fProcID;

while (fn<fNevents)
{

  if (fn>fNevents)
  {
    std::cout << "3 " << fn <<"\n";
    fLogger->Info(MESSAGE_ORIGIN,"End of file reached %i",fNevents);
    return kFALSE;
  }


  fTree->SetBranchAddress("pythiaid",&fPdgCode);
  fTree->SetBranchAddress("mother_id",&fMotherId);
  fTree->SetBranchAddress("px",&fPx);
  fTree->SetBranchAddress("py",&fPy);
  fTree->SetBranchAddress("pz",&fPz);
  fTree->SetBranchAddress("x",&fStartX);
  fTree->SetBranchAddress("y",&fStartY);
  fTree->SetBranchAddress("z",&fStartZ);
  fTree->SetBranchAddress("w",&fW);
  fTree->SetBranchAddress("process_id",&fProcID);



  fTree->GetEntry(fn);

  if (fn % 1000 == 0)
  {
    std::cout << fn <<"\n";
  }
  //add track to simulation

Double_t r = 5 + 0.8 * gRandom->Gaus();
     Double_t phi = gRandom->Uniform(0., 2.) * TMath::Pi();
     dx = r * TMath::Cos(phi);
     dy = r * TMath::Sin(phi);

  // cpg->AddTrack(13,px2,py2,pz2,x2,y2,z2+2084.5,-1,true,ecut2,0,1);
  // cpg->AddTrack(13,px2,py2,pz2,x2,y2,z2,-1,true,ecut2,0,1);
  // std::cout << 13 <<" "<<px2<<" "<<py2<<" "<<pz2<<" "<<x2<<" "<<y2<<" "<<-6542<<" "<<-1<<" "<<true<<" "<<ecut2<<" "<<0<<" "<<1 <<"\n";
  // std::cout << "fStartZ" << fStartZ<<"\n";
  cpg->AddTrack(fPdgCode,fPx,fPy,fPz,fStartX+dx,fStartY+dy,fStartZ+2228.5-144.0,-1,true,0,0,fW);

  fn++;

  return kTRUE;
}
}

// -------------------------------------------------------------------------
Int_t MuonBackGenerator::GetNevents()
{
 return fNevents;
}
void MuonBackGenerator::CloseFile()
{
 fInputFile->Close();
 fInputFile->Delete();
 delete fInputFile;
}

ClassImp(MuonBackGenerator)


