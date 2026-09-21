//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/ClusterSequence.hh"
//#include "fastjet/Selector.hh"
#include "llp_MuonSystem_CA_TrigEff_mdsnano.h"
#include "RazorHelper.h"
#include "RazorAnalyzer.h"
#include "TreeMuonSystemCombination_TnP.h"
#include "TreeMuonSystem_Skim_Merge_TnP.h"

#include "CACluster.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "Math/Vector3D.h"
#include <iostream>
#include <random>
#include <bitset>
//C++ includes
#include "assert.h"

//ROOT includes
#include "TH1F.h"

//using namespace fastjet;
using namespace std::chrono;
using namespace std;
using namespace ROOT::Math;

struct greater_than_pt
{
  inline bool operator() (const TLorentzVector& p1, const TLorentzVector& p2){return p1.Pt() > p2.Pt();}
};

struct leptons
{
  TLorentzVector lepton;
  int pdgId;
  float dZ;
  // bool passLooseId;
  // bool passMediumId;
  bool passId;
  bool passVetoId;
  bool passLooseIso;
  bool passTightIso;
  bool passVTightIso;
  bool passVVTightIso;
  bool passLooseId;
  bool isGlobal;

};


struct jets
{
  TLorentzVector jet;
  float time;
  bool passId;
  // bool passLooseId;
  // bool passMediumId;
  // bool passTightId;
  bool isCSVL;
  int ecalNRechits;
  float ecalRechitE;
  float jetChargedEMEnergyFraction;
  float jetNeutralEMEnergyFraction;
  float jetChargedHadronEnergyFraction;
  float jetNeutralHadronEnergyFraction;
  bool jetPassMuFrac;
  float jetPtJESUp;
  float jetPtJESDown;
  float jetEJESUp;
  float jetEJESDown;
  float JecUnc;

  float electronEnergyFraction;
  float neutralEmEnergyFraction;
  float chargedHadronEnergyFraction;
  float neutralHadronEnergyFraction;
  float muonEnergyFraction;
  float highPtJet; //added by Alex

};

//lepton highest pt comparator
struct largest_pt
{
  inline bool operator() (const leptons& p1, const leptons& p2){return p1.lepton.Pt() > p2.lepton.Pt();}
} my_largest_pt;

//jet highest pt comparator
struct largest_pt_jet
{
  inline bool operator() (const jets& p1, const jets& p2){return p1.jet.Pt() > p2.jet.Pt();}
} my_largest_pt_jet;


void llp_MuonSystem_CA_TrigEff_mdsnano::Analyze(bool isData, int options, string outputfilename, string analysisTag)
{
  //initialization: create one TTree for each analysis box
  cout << "Initializing..." << endl;
  cout << "IsData = " << isData << "\n";
  cout << "options = " << options << "\n";

  //---------------------------
  //options format: MH/MX/ctau/condor: 1000/300/0/1
  // mh can be 3-4 digits, mx is always 3 digits, ctau is one digit(number of zeros), last digit is condor option
  // mh can be 3-4 digits, mx is always 3 digits, ctau is 2 digit(number of zeros), last digit is condor option
  //
  //
  // int mx = int(options/1000)%1000;
  // int mh = options/1000000;
  // int ctau = pow(10, int(options/10)%10) * int(int(options/100)%10);
  //
  // cout<<"mh "<<mh<<", mx "<<mx<<", ctau "<<ctau<<endl;

  bool signalScan = int(options/10) == 1;
  int option = options%10;
  // if (options % 1){
  //   option = 1; // used when running condor
  // }
  // else{
  //   option = 0;// used when running locally
  // }

  if( isData )
  {
    std::cout << "[INFO]: running on data with option: " << option << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running on MC with option: " << option << std::endl;
  }
  if( signalScan )
  {
    std::cout << "[INFO]: running with Signal scan" << std::endl;
  }
  else
  {
    std::cout << "[INFO]: running without Signal scan " << option << std::endl;
  }



  const float ELE_MASS = 0.000511;
  const float MU_MASS  = 0.105658;
  const float Z_MASS   = 91.2;

  if (analysisTag == ""){
    analysisTag = "Razor2016_80X";

  }


  const int zh_lepton0_cut = 15;
  const int zh_lepton1_cut = 15;

  const int wh_Muon_pt_cut = 25;
  const int wh_elePt_cut = 35;



  //-----------------------------------------------
  //Set up Output File
  //-----------------------------------------------
  string outfilename = outputfilename;
  if (outfilename == "") outfilename = "MuonSystem_Tree.root";
  TFile *outFile;
  if (isData || !signalScan) outFile = new TFile(outfilename.c_str(), "RECREATE");


  TreeMuonSystem_Skim_Merge_TnP *MuonSystem = new TreeMuonSystem_Skim_Merge_TnP();
  MuonSystem->CreateTree();
  MuonSystem->tree_->SetAutoFlush(0);
  MuonSystem->InitTree();

  // for signals, need one output file for each signal point
  map<pair<int,int>, TFile*> Files2D;
  map<pair<int,int>, TTree*>Trees2D;
  map<pair<int,int>, TH1F*> NEvents2D;
  map<pair<int,int>, TH1F*> accep2D;
  map<pair<int,int>, TH1F*> accep_met2D;
  map<pair<int,int>, TH1F*> Total2D;



  //histogram containing total number of processed events (for normalization)
  TH1F *NEvents = new TH1F("NEvents", "NEvents", 1, 1, 2);
  TH1F *Total = new TH1F("Total", "Total", 1, 1, 2);

  TH1F *accep = new TH1F("accep", "acceptance", 1, 1, 2);
  TH1F *accep_met = new TH1F("accep_met", "acceptance_met", 1, 1, 2);

  TH1F *Nmet200 = new TH1F("Nmet200", "Nmet200", 1, 1, 2);
  TH1F *NmetFilter = new TH1F("NmetFilter", "NmetFilter", 1, 1, 2);
  TH1F *Nlep0 = new TH1F("Nlep0", "Nlep0", 1, 1, 2);
  TH1F *Njet1 = new TH1F("Njet1", "Njet1", 1, 1, 2);
  TH1F *NcosmicVeto = new TH1F("NcosmicVeto", "NcosmicVeto", 1, 1, 2);

  //JetDefinition jet_def( antikt_algorithm, .4 );
  //fastjet::JetDefinition jet_def(fastjet::cambridge_algorithm, 0.4);

  //vector<fastjet::PseudoJet> input_particles;


  char* cmsswPath;
  cmsswPath = getenv("CMSSW_BASE");
  string pathname;
  if(cmsswPath != NULL) pathname = string(cmsswPath) + "/src/llp_analyzer/data/JEC/";
  if(cmsswPath != NULL and option == 1) pathname = "JEC/"; //run on condor if option == 1




  //--------------------------------
  //Initialize helper
  //--------------------------------
  cout << "Pre Razor Helper" << endl;
  RazorHelper *helper = 0;
  helper = new RazorHelper(analysisTag, isData);
  cout << "Post Razor Helper" << endl;


  //*************************************************************************
  //Look over Input File Events
  //*************************************************************************
  if (fChain == 0) return;
  cout << "Total Events: " << fChain->GetEntries() << "\n";
  Long64_t nbytes = 0, nb = 0;
  clock_t start, end;
  start = clock();
  //float maxEvents = 505000;
  
  //variables for checking acceptance on cut-by-cut basis
  float greater_than_one_tag = 0;
  float two_satsifying_probe = 0;
  float leading_muon_pt = 0;
  float subleading_muon_pt = 0;
  float Zmumu_events_passed = 0;
  float total_events_passed = 0; 

  bool checkProbeMuons = false; //this only counts the number of probe muons and skips rest of code!

  
  float maxEvents = fChain->GetEntries();
  for (Long64_t jentry=0; jentry<maxEvents; jentry++) {
    //begin event
    if(jentry % 1000 == 0)
    {
      end = clock();
      double time_taken = double(end - start) / double(CLOCKS_PER_SEC);
      cout << "Processing entry " << jentry << endl;
      cout << "Time taken by program is : " << time_taken << endl;
      start = clock();
    }
    /*
    if ((jentry>6086000 && jentry<6087000)||(jentry>6743000 && jentry<6744000)||(jentry>4044000 && jentry<4045000)||(jentry>5008000 && jentry<5009000)||(jentry>7767000 && jentry<7768000) || (jentry>8023000 && jentry<8024000)){
	    continue;
    }
    */    
    // if (jentry<4000)continue;
    // cout<<jentry<<endl;
    //cout<<"About to load tree"<<endl;
    Long64_t ientry = LoadTree(jentry);
    //cout<<"Loaded tree"<<endl;
    if (ientry < 0) break;
    //GetEntry(ientry);
    //cout<<"Getting entry from fChain"<<endl;
    nb = fChain->GetEntry(jentry); nbytes += nb;
    //cout<<nbytes<<endl;
    //cout<<"Trying to Get Entry"<<endl;
    fChain->GetEntry(jentry);
    //cout<<"Got Entry"<<endl;
    //cout<<"RunNum: "<<runNum<<endl;
    //cout<<"eventNum: "<<eventNum<<endl;
    //cout<<"Got entry from fChain, about to check HLT_IsoMu24"<<endl;
    
    //cout<<"Passed HLT_IsoMu24, about to init variables in MuonSystem"<<endl;
    //cout<<"Past Cut"<<endl;
    //fill normalization histogram
    //cout<<"About to Init Variables"<<endl;
    MuonSystem->InitVariables();
    auto nCscRechits = ncscRechits;
    //cout<<"nCscRechits: "<<nCscRechits<<endl;
    
    
    //branch mapping for MDS NANO variables
    auto nCscSeg = ncscSegments;
    auto nDtSeg = ndtSegments;
    auto nDtRechits = ndtRecHits;
    auto nRpc = nrpcRecHits;
    auto cscRechitsPhi = cscRechits_Phi;
    auto cscRechitsEta = cscRechits_Eta;
    auto cscRechitsX = cscRechits_X;
    auto cscRechitsY = cscRechits_Y;
    auto cscRechitsZ = cscRechits_Z;
    auto cscRechitsTpeak = cscRechits_Tpeak;
    auto cscRechitsTwire = cscRechits_Twire;
    auto cscRechitsChamber = cscRechits_Chamber;
    auto cscRechitsStation = cscRechits_Station;

    auto dtRechitCorrectX = dtRecHits_X;
    auto dtRechitCorrectY = dtRecHits_Y;
    auto dtRechitCorrectZ = dtRecHits_Z;
    auto dtRechitCorrectEta = dtRecHits_Eta;
    auto dtRechitCorrectPhi = dtRecHits_Phi;
    auto dtRechitStation = dtRecHits_Station;
    auto dtRechitWheel = dtRecHits_Wheel;
    auto dtRechitSuperLayer = dtRecHits_SuperLayer;

    auto dtSegStation = dtSegments_Station;
    auto dtSegPhi = dtSegments_Phi;
    auto dtSegEta = dtSegments_Eta;
    auto rpcPhi = rpcRecHits_Phi;
    auto rpcEta = rpcRecHits_Eta;
    auto rpcX = rpcRecHits_X;
    auto rpcY = rpcRecHits_Y;
    auto rpcZ = rpcRecHits_Z;
    
    auto rpcBx = rpcRecHits_Bx;
    auto rpcRegion = rpcRecHits_Region;
    auto rpcRing = rpcRecHits_Ring;
    //cout<<"Inited Variables"<<endl;
    // std::cout << "deb1 " << jentry << std::endl;
    //cout<<"Initialized Variables in Muon System"<<endl;


    /*
    if (!isData && signalScan)
    {

      string mh_substring = lheComments->substr(lheComments->find("MH-")+3);
      int mh = stoi(mh_substring.substr(0,mh_substring.find('_')));
      string mx_substring = lheComments->substr(lheComments->find("MS-")+3);
      int mx = stoi(mx_substring.substr(0,mx_substring.find('_')));
      string ctau_substring = lheComments->substr(lheComments->find("ctauS-")+6);
      int ctau = stoi(ctau_substring.substr(0,ctau_substring.find('_')));
      MuonSystem->mH = mh;
      MuonSystem->mX = mx;
      MuonSystem->ctau = ctau;

      // if (mh2 != mh || mx2!=mx || ctau2!=ctau) continue;
      // cout<<*lheComments<<endl;

      pair<int,int> signalPair = make_pair(mx, ctau);

      if (Files2D.count(signalPair) == 0){ //create file and tree
        //format file name
        string thisFileName = outfilename;
        thisFileName.erase(thisFileName.end()-5, thisFileName.end());
        thisFileName += "_" + to_string(mx) + "_" + to_string(ctau) + ".root";

        Files2D[signalPair] = new TFile(thisFileName.c_str(), "recreate");
        Trees2D[signalPair] =  MuonSystem->tree_->CloneTree(0);
        NEvents2D[signalPair] = new TH1F(Form("NEvents%d%d", mx, ctau), "NEvents", 1,0.5,1.5);
        Total2D[signalPair] = new TH1F(Form("Total%d%d", mx, ctau), "Total", 1,0.5,1.5);
        accep2D[signalPair] = new TH1F(Form("accep2D%d%d", mx, ctau), "acceptance", 1,0.5,1.5);
        accep_met2D[signalPair] = new TH1F(Form("accep_met2D%d%d", mx, ctau), "acceptance_met", 1,0.5,1.5);



        cout << "Created new output file " << thisFileName << endl;
      }
          //Fill NEvents hist
      NEvents2D[signalPair]->Fill(1.0, genWeight);


    }
    */
    //event info
    if (isData)
    {
      NEvents->Fill(1);
      MuonSystem->weight = 1;
    }
    
    else
    {
      MuonSystem->weight = genWeight;
      NEvents->Fill(1, genWeight);
    }
    if (!HLT_IsoMu24) continue;
    MuonSystem->runNum = run;
    MuonSystem->lumiSec = luminosityBlock;
    MuonSystem->evtNum = event;

    //check event flags
    MuonSystem->Flag_goodVertices = Flag_goodVertices;
      MuonSystem->Flag_globalSuperTightHalo2016Filter = Flag_globalSuperTightHalo2016Filter;
      MuonSystem->Flag_EcalDeadCellTriggerPrimitiveFilter = Flag_EcalDeadCellTriggerPrimitiveFilter;
      MuonSystem->Flag_BadPFMuonFilter = Flag_BadPFMuonFilter;
      MuonSystem->Flag_BadPFMuonDzFilter = Flag_BadPFMuonDzFilter;
      MuonSystem->Flag_hfNoisyHitsFilter = Flag_hfNoisyHitsFilter;
      MuonSystem->Flag_eeBadScFilter = Flag_eeBadScFilter;
      MuonSystem->Flag_all = Flag_eeBadScFilter && Flag_hfNoisyHitsFilter && Flag_BadPFMuonDzFilter && Flag_BadPFMuonFilter && Flag_EcalDeadCellTriggerPrimitiveFilter
                              && Flag_globalSuperTightHalo2016Filter && Flag_goodVertices;
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;

      // Flag_ecalBadCalibFilter for nanoAOD: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#ECal_BadCalibration_Filter_Flag
      if (analysisTag == "Summer24") MuonSystem->Flag_ecalBadCalibFilter = Flag_ecalBadCalibFilter;
      else{
        MuonSystem->Flag_ecalBadCalibFilter = true;
        if (isData && run >= 362433 && run<=367144)
        {
          if (PuppiMET_pt > 100)
          {
            for(int i = 0; i < nJet; i++)
            {
              if (Jet_pt[i]<50) continue;
              if (!(Jet_eta[i] <= -0.1 && Jet_eta[i]>=-0.5 && Jet_phi[i] <-1.8 && Jet_phi[i]> -2.1)) continue;
              if (!(Jet_neEmEF[i] >0.9 || Jet_chEmEF[i]>0.9)) continue;
              if (deltaPhi(PuppiMET_phi, Jet_phi[i])<2.9) continue;
              Flag_ecalBadCalibFilter = false;
            }
          }
        }
      }
     // jet veto map, following selections here: https://cms-jerc.web.cern.ch/Recommendations/#jet-veto-maps
      MuonSystem->jetVeto = true;
      
      for(int i = 0; i < nJet; i++)
      {
        if (Jet_pt[i] <= 15) continue;
        if (Jet_neEmEF[i] + Jet_chEmEF[i] >= 0.9) continue;
        std::bitset<sizeof(int) * 8> jetID(Jet_jetId[i]);
        //cout<<"Filter bit in base 10: "<<binaryNumber<<endl;
        if (!jetID.test(1)) continue;
        //if (!jetPassIDTight[i]) continue;
        //remove overlaps
        bool overlap = false;
        for(int j = 0; j < nMuon; j++)
        {
          if (!Muon_isPFcand[j])continue;
          if (RazorAnalyzerMerged::deltaR(Jet_eta[i],Jet_phi[i],Muon_eta[j], Muon_phi[j]) < 0.2) overlap = true;
        }
        if (overlap) continue;
        //cout<<"about to call jet veto map"<<endl;
        helper->getJetVetoMap(0,1);
        if (helper->getJetVetoMap(Jet_eta[i],Jet_phi[i])>0.0) MuonSystem->jetVeto = false;
        if (analysisTag == "Summer24" && helper->getJetVetoFpixMap(Jet_eta[i],Jet_phi[i])>0.0) MuonSystem->jetVeto = false;
      } 
    if (!isData){
      MuonSystem->npu = Pileup_nTrueInt;
      MuonSystem->pileupWeight = helper->getPileupWeight(Pileup_nTrueInt);
      MuonSystem->pileupWeightUp = helper->getPileupWeightUp(Pileup_nTrueInt) / MuonSystem->pileupWeight;
      MuonSystem->pileupWeightDown = helper->getPileupWeightDown(Pileup_nTrueInt) / MuonSystem->pileupWeight;
    }



    /*commented out for new make class
    if (!isData)
    {
       for(int i = 0; i < 2;i++)
       {
	       MuonSystem->gLLP_eta[MuonSystem->nGLLP] = gLLP_eta[i];
         MuonSystem->gLLP_phi[MuonSystem->nGLLP] = gLLP_phi[i];
         MuonSystem->gLLP_e[MuonSystem->nGLLP] = gLLP_e[i];
         MuonSystem->gLLP_pt[MuonSystem->nGLLP] = gLLP_pt[i];

         MuonSystem->gLLP_decay_vertex_r[MuonSystem->nGLLP] = sqrt(gLLP_decay_vertex_x[i]*gLLP_decay_vertex_x[i]+gLLP_decay_vertex_y[i]*gLLP_decay_vertex_y[i]);
         MuonSystem->gLLP_decay_vertex_x[MuonSystem->nGLLP] = gLLP_decay_vertex_x[i];
         MuonSystem->gLLP_decay_vertex_y[MuonSystem->nGLLP] = gLLP_decay_vertex_y[i];
         MuonSystem->gLLP_decay_vertex_z[MuonSystem->nGLLP] = gLLP_decay_vertex_z[i];
         float beta = gLLP_beta[i];
         float gLLP_decay_vertex = sqrt(pow(MuonSystem->gLLP_decay_vertex_r[i], 2) + pow(MuonSystem->gLLP_decay_vertex_z[i],2));
         float gamma = 1.0/sqrt(1-beta*beta);
         MuonSystem->gLLP_ctau[MuonSystem->nGLLP] = gLLP_decay_vertex/(beta * gamma);
         MuonSystem->gLLP_beta[MuonSystem->nGLLP] = gLLP_beta[i];


           if (abs(MuonSystem->gLLP_eta[i]) < 2.4
             && abs(MuonSystem->gLLP_decay_vertex_z[i])<1100 && abs(MuonSystem->gLLP_decay_vertex_z[i])>400
             && MuonSystem->gLLP_decay_vertex_r[i] < 695.5) MuonSystem->gLLP_csc[MuonSystem->nGLLP] = true;
           if (abs(MuonSystem->gLLP_decay_vertex_z[i])< 661.0
             && MuonSystem->gLLP_decay_vertex_r[i] < 800
              && MuonSystem->gLLP_decay_vertex_r[i] > 200.0) MuonSystem->gLLP_dt[MuonSystem->nGLLP] = true;

            MuonSystem->nGLLP++;
       }
      

    }//end of isData
    */
      //get NPU
      //MuonSystem->npv = nPV; commented out for new make class
      //MuonSystem->rho = fixedGridRhoFastjetAll; commented out for new make class
      MuonSystem->met = MET_pt; //modified
      MuonSystem->metPhi = MET_phi; //modified
      /* commented out for new make class
      if(signalScan && !isData)Total2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      if(signalScan && !isData)
      {
        accep2D[make_pair(MuonSystem->mX, MuonSystem->ctau)]->Fill(1.0, genWeight*MuonSystem->pileupWeight);
      }
      else if (!isData)
      {
        accep->Fill(1.0, genWeight*MuonSystem->pileupWeight);

      }
      */

      //Triggers
      //for(int i = 0; i < NTriggersMAX; i++){
        //MuonSystem->HLTDecision[i] = HLTDecision[i];
        //cout<<"Loading triggers"<<endl;
        MuonSystem->HLT_CscCluster_Loose = HLT_CscCluster_Loose;
        MuonSystem->HLT_CscCluster_Medium = HLT_CscCluster_Medium;
        MuonSystem->HLT_CscCluster_Tight = HLT_CscCluster_Tight;
        MuonSystem->HLT_IsoMu20 = HLT_IsoMu20;
        MuonSystem->HLT_IsoMu24 = HLT_IsoMu24;
        MuonSystem->L1_SingleMuShower_Nominal = L1_SingleMuShower_Nominal;
        MuonSystem->L1_SingleMuShower_Tight = L1_SingleMuShower_Tight;
      //}

      //*************************************************************************
      //Start Object Selection
      //*************************************************************************

      std::vector<leptons> Leptons_Probe;
      std::vector<leptons> Leptons_notProbe;
      //-------------------------------
      //Muons
      //-------------------------------
      
      //cout<<"nMuon "<<nMuon<<endl;
      bool one_good_muon = false; 
      for( int i = 0; i < nMuon; i++ )  
      {
       	/*commented out for clusters not matched with muons
	if(!Muon_looseId[i]) continue;
        //cout<<"Found not loose muon"<<endl;
        //cout<<"Muon_pt: "<<Muon_pt[i]<<endl;
        if(Muon_pt[i] < 20.0) continue; //do we want this cut, only hard muons?
        //cout<<"Found hard, not loose muon"<<endl;
        //cout<<"Muon_eta"<<Muon_eta[i]<<endl;
        if(fabs(Muon_eta[i]) > 2.4) continue;
	*/
        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons_Probe)
        {
          if (RazorAnalyzerMerged::deltaR(Muon_eta[i],Muon_phi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;

        leptons tmpMuon;
        tmpMuon.lepton.SetPtEtaPhiM(Muon_pt[i],Muon_eta[i], Muon_phi[i], MU_MASS);
        tmpMuon.pdgId = 13 * -1 * Muon_charge[i];
        tmpMuon.dZ = Muon_dz[i];
        tmpMuon.passId = Muon_tightId[i];
        //ask about iso variables
        //float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / Muon_pt[i];
        float muonIso = Muon_pfRelIso04_all[i]; //added by alex, address set to Muon_pfRelIso04_all branch orignally fron nanoAODs
        //tmpMuon.passLooseIso = muonIso<0.25;
        tmpMuon.passLooseIso = muonIso<0.40; //loosen iso for probe muon
        tmpMuon.passTightIso = muonIso<0.15;
        tmpMuon.passVTightIso = muonIso<0.10;
        tmpMuon.passVVTightIso = muonIso<0.05;
        tmpMuon.passLooseId = Muon_looseId[i];
        tmpMuon.isGlobal = Muon_isGlobal[i];


        tmpMuon.passVetoId = false;
	    if (Muon_looseId[i] && Muon_pt[i] >= 20.0 && fabs(Muon_eta[i] < 2.4) && tmpMuon.passLooseIso){
		    one_good_muon = true;
        Leptons_Probe.push_back(tmpMuon);
	      }
        else{
          Leptons_notProbe.push_back(tmpMuon);
        }
        //if (!tmpMuon.passLooseIso) continue; Commented out for non-muon clusters
        
      }
      if(!one_good_muon) continue;
      MuonSystem->numProbeMuons = Leptons_Probe.size();
      
      if (checkProbeMuons){
         MuonSystem->tree_->Fill();
        continue;
      }
      //-------------------------------
      //Electrons
      //-------------------------------
      
      /*COMMENTED OUT BECAUSE MERGED NTUPLES DO NOT CONTAIN ELECTRONS
      for( int i = 0; i < nElectrons; i++ )
      {
        if (!ele_passCutBasedIDVeto[i]) continue;
        if(elePt[i] < 35) continue;
        if(fabs(eleEta[i]) > 2.5) continue;

        //remove overlaps
        bool overlap = false;
        for(auto& lep : Leptons)
        {
          if (RazorAnalyzer_TnP::deltaR(eleEta[i],elePhi[i],lep.lepton.Eta(),lep.lepton.Phi()) < 0.3) overlap = true;
        }
        if(overlap) continue;
        leptons tmpElectron;
        tmpElectron.lepton.SetPtEtaPhiM(elePt[i],eleEta[i], elePhi[i], ELE_MASS);
        tmpElectron.pdgId = 11 * -1 * eleCharge[i];
        tmpElectron.dZ = ele_dZ[i];
        tmpElectron.passId = ele_passCutBasedIDTight[i];
        Leptons.push_back(tmpElectron);
      }
      */
      
      sort(Leptons_Probe.begin(), Leptons_Probe.end(), my_largest_pt);
      sort(Leptons_notProbe.begin(), Leptons_notProbe.end(), my_largest_pt);
      //cout<< "length of leptons: " << Leptons.size() << endl;
      
      for ( auto &tmp : Leptons_Probe ) //note that these lepton variables will only include muon info
      {
        MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
        MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
        MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
        MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
        MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
        MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
        MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
        MuonSystem->lepPassProbe[MuonSystem->nLeptons] = true;
        MuonSystem->nLeptons++;
      }

      for ( auto &tmp : Leptons_notProbe ) //note that these lepton variables will only include muon info
      {
        MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
        MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
        MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
        MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
        MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
        MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
        MuonSystem->lepTightId[MuonSystem->nLeptons] = tmp.passId;
        MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
        MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
        MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
        MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
        MuonSystem->lepPassProbe[MuonSystem->nLeptons] = false;
        MuonSystem->nLeptons++;
      }
      
    // the following code , until jets, is for reconsturction of the Z-peak. It is taken from https://github.com/cms-lpc-llp/llp_analyzer/blob/master/analyzers/llp_MuonSystem_TnP.cc#L489-L526
    // now, we will comment out the Z-reconstruction code, since we are no longer doing Tag-N-Probe for the trigger efficiency measurement
    /*
    double ZMass = -999;
    double ZPt = -999;
    double tmpDistToZPole = 9999;
    pair<uint,uint> ZCandidateLeptonIndex;
    bool foundZ = false;
    TLorentzVector ZCandidate;
    double leadingLepPt = 0.0;
    int tagCount=0;
    for( uint i = 0; i < Leptons.size(); i++ )
    {
      for( uint j = i+1; j < Leptons.size(); j++ )
      {
        if (!( Leptons[i].pdgId == -1*Leptons[j].pdgId )) continue;// same flavor opposite charge
        double tmpMass = (Leptons[i].lepton+Leptons[j].lepton).M();
        //select the pair closest to Z pole mass
        if ( fabs( tmpMass - Z_MASS) < tmpDistToZPole)
        {
          tmpDistToZPole = tmpMass;
          if (Leptons[i].pdgId > 0)
          {
            ZCandidateLeptonIndex = pair<int,int>(i,j);
          }
          else
          {
            ZCandidateLeptonIndex = pair<int,int>(j,i);
          }
          ZMass = tmpMass;
          ZPt = (Leptons[i].lepton+Leptons[j].lepton).Pt();
          ZCandidate = Leptons[i].lepton+Leptons[j].lepton;
          leadingLepPt = max(Leptons[i].lepton.Pt(),Leptons[j].lepton.Pt());
          foundZ = true;
        }
      }
    }

    // if (foundZ  && Leptons.size() == 2 && leadingLepPt > zh_lepton0_cut)
    if (foundZ  && Leptons.size() == 2 )
    {
      MuonSystem->ZMass = ZMass;
      MuonSystem->ZPt   = ZPt;
      MuonSystem->ZEta  = ZCandidate.Eta();
      MuonSystem->ZPhi  = ZCandidate.Phi();
      MuonSystem->ZleptonIndex1 = ZCandidateLeptonIndex.first;
      MuonSystem->ZleptonIndex2 = ZCandidateLeptonIndex.second;
      MuonSystem->category = 2;
    } // endif foundZ
    else{
      for ( unsigned int i = Leptons.size(); i>0; --i )
      {
        int index = i-1;
        if (abs(Leptons[index].pdgId) == 13 && Leptons[index].lepton.Pt() < wh_Muon_pt_cut)
        {
           Leptons.erase(Leptons.begin() + index);
        }
        else if (abs(Leptons[index].pdgId) == 11 && Leptons[index].lepton.Pt() < wh_elePt_cut)
        {
          Leptons.erase(Leptons.begin() + index);
        }
      }
      if (Leptons.size() == 1) MuonSystem->category = 1;
      else MuonSystem->category = 0;
    }
    bool tag = false;
    for ( auto &tmp : Leptons )
    {
      MuonSystem->lepE[MuonSystem->nLeptons]      = tmp.lepton.E();
      MuonSystem->lepPt[MuonSystem->nLeptons]     = tmp.lepton.Pt();
      MuonSystem->lepEta[MuonSystem->nLeptons]    = tmp.lepton.Eta();
      MuonSystem->lepPhi[MuonSystem->nLeptons]    = tmp.lepton.Phi();
      MuonSystem->lepPdgId[MuonSystem->nLeptons]  = tmp.pdgId;
      MuonSystem->lepDZ[MuonSystem->nLeptons]     = tmp.dZ;
      MuonSystem->lepPassId[MuonSystem->nLeptons] = tmp.passId;
      MuonSystem->lepPassLooseIso[MuonSystem->nLeptons] = tmp.passLooseIso;
      MuonSystem->lepPassTightIso[MuonSystem->nLeptons] = tmp.passTightIso;
      MuonSystem->lepPassVTightIso[MuonSystem->nLeptons] = tmp.passVTightIso;
      MuonSystem->lepPassVVTightIso[MuonSystem->nLeptons] = tmp.passVVTightIso;
      
      if (!isData)
      {
        Commenting out all razor helper (corrections specific to Run 2?)
        MuonSystem->lepTriggerSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons] = helper->getMuonScaleFactor(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), true, true, true);
        MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, false);
        MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, true, true);
        MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, false);
        MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons] = helper->getMuonMCEfficiency(tmp.lepton.Pt(), tmp.lepton.Eta(), false, false, true);
        
      }


      
      if (MuonSystem->lepPassTightIso[MuonSystem->nLeptons] && MuonSystem->lepPassId[MuonSystem->nLeptons])
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = true;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepTightIdSF[MuonSystem->nLeptons] * MuonSystem->lepTightIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepTightIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
        tag = true;
        tagCount++;
      }
      else
      {
        MuonSystem->lepTag[MuonSystem->nLeptons] = false;
        if(!isData)
        {
          MuonSystem->lepSF[MuonSystem->nLeptons] = MuonSystem->lepTriggerSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIdSF[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoSF[MuonSystem->nLeptons];
          MuonSystem->lepEff[MuonSystem->nLeptons] = MuonSystem->lepTriggerMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIdMCEfficiency[MuonSystem->nLeptons] * MuonSystem->lepLooseIsoMCEfficiency[MuonSystem->nLeptons];
        }
        else{
          MuonSystem->lepSF[MuonSystem->nLeptons] = 1.0;
          MuonSystem->lepEff[MuonSystem->nLeptons] = 1.0;
        }
      }



      if(!isData)
      {
        for (int i=0; i < nGenParticle; i++)
        { float tmpDR = deltaR(gParticleEta[i],gParticlePhi[i],MuonSystem->lepEta[MuonSystem->nLeptons],MuonSystem->lepPhi[MuonSystem->nLeptons]);
          if ((abs(gParticleId[i]) == 13) && abs(gParticleMotherId[i]) == 23 && tmpDR<0.4) MuonSystem->lepFromZ[MuonSystem->nLeptons] = true;
        }
      }
      MuonSystem->nLeptons++;
    }
    
    //compute acceptances
    if(MuonSystem->nLeptons==2)two_satsifying_probe++;
    //if (abs(MuonSystem->lepPdgId[0])!=13)continue; not relevant, all muons
    //if (abs(MuonSystem->ZMass)<50)continue; removed to see whole Z-peak
    if (abs(MuonSystem->lepPt[0])>35)leading_muon_pt++;
    if (abs(MuonSystem->lepPt[1])>20)subleading_muon_pt++;
    // if(MuonSystem->ZMass<120)continue; remove to see whole Z-peak
    if (tagCount>0)greater_than_one_tag++;


    if(MuonSystem->nLeptons!=2)continue;
    if(MuonSystem->category!=2)continue;
    if (abs(MuonSystem->lepPdgId[0])!=13)continue;
    //if (abs(MuonSystem->ZMass)<50)continue; removed to see whole Z-peak
    if (abs(MuonSystem->lepPt[0])<35)continue;
    if (abs(MuonSystem->lepPt[1])<20)continue;
    // if(MuonSystem->ZMass<120)continue; remove to see whole Z-peak
    if (tagCount==0) continue;
    //cout<<"at lep tag"<<endl;
    //if(MuonSystem->lepTag[0] == MuonSystem->lepTag[1]) continue;
    //cout<<"made it past lepTag cuts" << endl;
      //require one tag one probe
    if(!isData)
    {
      MuonSystem->lepOverallSF = 1.0 - (1.0 - MuonSystem->lepSF[0] * MuonSystem->lepEff[0]) * (1 - MuonSystem->lepSF[1] * MuonSystem->lepEff[1]);
      MuonSystem->lepOverallSF = MuonSystem->lepOverallSF / (1.0 - (1.0 - MuonSystem->lepEff[0]) * (1 - MuonSystem->lepEff[1]));
    }
    else{
      MuonSystem->lepOverallSF = 1.0;
    }

    //cout<<"made it past muon cuts!" << endl;
  
  TLorentzVector met;
  met.SetPtEtaPhiE(metType1Pt,0,metType1Phi,metType1Pt);
  if ( Leptons.size() > 0 )
  {
    TLorentzVector visible = Leptons[0].lepton;
    MuonSystem->MT = GetMT(visible,met);
  }
  */
   Zmumu_events_passed++;   
    //cout<<"past z stuff"<<endl;
    //-----------------------------------------------
    //Select Jets
    //-----------------------------------------------
    
    std::vector<jets> Jets;
    //cout<<"nJet: "<<nJet<<endl;
    for(int i = 0; i < nJet; i++)
    {
      
      if (fabs(Jet_eta[i]) >= 3.0)continue;
      if (i<4){
        //cout<<Jet_pt[i]<<endl;
      }
      //if( Jet_pt[i] < 30 ) continue; //changed from 20 to 30
      //cout<<"Found Jet with pT > 30"<<endl;
      //if (!jetPassIDLoose[i]) continue;
      //------------------------------------------------------------
      //exclude selected muons and electrons from the jet collection
      //------------------------------------------------------------
      double deltaR = -1;
      for(auto& lep : Leptons_Probe){
        double thisDR = RazorAnalyzerMerged::deltaR(Jet_eta[i],Jet_phi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton

      deltaR = -1;
      for(auto& lep : Leptons_notProbe){
        double thisDR = RazorAnalyzerMerged::deltaR(Jet_eta[i],Jet_phi[i],lep.lepton.Eta(),lep.lepton.Phi());
        if(deltaR < 0 || thisDR < deltaR) deltaR = thisDR;
      }
      //if(deltaR > 0 && deltaR < 0.4) continue; //jet matches a selected lepton


      TLorentzVector thisJet = makeTLorentzVector( Jet_pt[i], Jet_eta[i], Jet_phi[i], Jet_mass[i] ); //was energy before, not sure about mass

      jets tmpJet;
      tmpJet.jet    = thisJet;
      //tmpJet.passId = jetPassIDTight[i];
      tmpJet.passId = true; //not sure of analog in Nano, changed to true for now
      tmpJet.highPtJet = Jet_pt[i] >= 30;
      Jets.push_back(tmpJet);

      }

      sort(Jets.begin(), Jets.end(), my_largest_pt_jet);



      for ( auto &tmp : Jets )
      {
        if(tmp.jet.Pt()<30)continue;

        MuonSystem->jetE[MuonSystem->nJets] = tmp.jet.E();
        MuonSystem->jetPt[MuonSystem->nJets] = tmp.jet.Pt();
        MuonSystem->jetEta[MuonSystem->nJets] = tmp.jet.Eta();
        MuonSystem->jetPhi[MuonSystem->nJets] = tmp.jet.Phi();
        MuonSystem->jetTightPassId[MuonSystem->nJets] = tmp.passId;


        MuonSystem->nJets++;
      }
      

     /*
      MuonSystem->nDTRechits  = 0;

      int nDTRechitsChamberMinus12 = 0;
      int nDTRechitsChamberMinus11 = 0;
      int nDTRechitsChamber10 = 0;
      int nDTRechitsChamberPlus11 = 0;
      int nDTRechitsChamberPlus12 = 0;
      int nDTRechitsChamberMinus22 = 0;
      int nDTRechitsChamberMinus21 = 0;
      int nDTRechitsChamber20 = 0;
      int nDTRechitsChamberPlus21 = 0;
      int nDTRechitsChamberPlus22 = 0;
      int nDTRechitsChamberMinus32 = 0;
      int nDTRechitsChamberMinus31 = 0;
      int nDTRechitsChamber30 = 0;
      int nDTRechitsChamberPlus31 = 0;
      int nDTRechitsChamberPlus32 = 0;
      int nDTRechitsChamberMinus42 = 0;
      int nDTRechitsChamberMinus41 = 0;
      int nDTRechitsChamber40 = 0;
      int nDTRechitsChamberPlus41 = 0;
      int nDTRechitsChamberPlus42 = 0;

      for (int i = 0; i < nDtRechits; i++) {
        //cout<<"check readings"<<endl;
        //cout<<dtRechitStation[i]<<endl;
        //cout<<dtRechitWheel[i]<<endl;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus12++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 0) nDTRechitsChamber10++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus11++;
        if (dtRechitStation[i] == 1 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus12++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus22++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 0) nDTRechitsChamber20++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus21++;
        if (dtRechitStation[i] == 2 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus22++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus32++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 0) nDTRechitsChamber30++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus31++;
        if (dtRechitStation[i] == 3 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus32++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -2) nDTRechitsChamberMinus42++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == -1) nDTRechitsChamberMinus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 0) nDTRechitsChamber40++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 1) nDTRechitsChamberPlus41++;
        if (dtRechitStation[i] == 4 && dtRechitWheel[i] == 2) nDTRechitsChamberPlus42++;




      }


      if ( nDTRechitsChamberMinus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber10 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus11 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus12 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber20 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus21 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus22 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber30 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus31 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus32 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus42 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberMinus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamber40 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus41 > 50) MuonSystem->nDtRings++;
      if ( nDTRechitsChamberPlus42 > 50) MuonSystem->nDtRings++;

      */

      vector<Rechits> points;
      vector<int> cscRechitsClusterId;
      points.clear();
      int nCscRechitsChamberPlus11 = 0;
      int nCscRechitsChamberPlus12 = 0;
      int nCscRechitsChamberPlus13 = 0;
      int nCscRechitsChamberPlus21 = 0;
      int nCscRechitsChamberPlus22 = 0;
      int nCscRechitsChamberPlus31 = 0;
      int nCscRechitsChamberPlus32 = 0;
      int nCscRechitsChamberPlus41 = 0;
      int nCscRechitsChamberPlus42 = 0;

      int nCscRechitsChamberMinus11 = 0;
      int nCscRechitsChamberMinus12 = 0;
      int nCscRechitsChamberMinus13 = 0;
      int nCscRechitsChamberMinus21 = 0;
      int nCscRechitsChamberMinus22 = 0;
      int nCscRechitsChamberMinus31 = 0;
      int nCscRechitsChamberMinus32 = 0;
      int nCscRechitsChamberMinus41 = 0;
      int nCscRechitsChamberMinus42 = 0;
      for (int i = 0; i < nCscRechits; i++) {
        //cout<<"entering rechits code"<<endl;
        //pick out the right bits for chamber
        int chamber = ((cscRechitsDetId[i] >> 3) & 077); //https://github.com/cms-sw/cmssw/blob/master/DataFormats/MuonDetId/interface/CSCDetId.h#L147

        int layer = (cscRechitsDetId[i] & 07);
        Rechits p;
        p.phi = cscRechitsPhi[i];
        p.eta = cscRechitsEta[i];
        p.x = cscRechitsX[i];
        p.y = cscRechitsY[i];
        p.z = cscRechitsZ[i];
        p.t = cscRechitsTpeak[i];
        p.twire = cscRechitsTwire[i];
        p.station = cscRechitsStation[i];
        p.chamber = cscRechitsChamber[i];
        p.layer = layer;
        p.superlayer = 0;
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);
        cscRechitsClusterId.push_back(-1);

        if (cscRechitsChamber[i] == 11)  nCscRechitsChamberPlus11++;
        if (cscRechitsChamber[i] == 12)  nCscRechitsChamberPlus12++;
        if (cscRechitsChamber[i] == 13)  nCscRechitsChamberPlus13++;
        if (cscRechitsChamber[i] == 21)  nCscRechitsChamberPlus21++;
        if (cscRechitsChamber[i] == 22)  nCscRechitsChamberPlus22++;
        if (cscRechitsChamber[i] == 31)  nCscRechitsChamberPlus31++;
        if (cscRechitsChamber[i] == 32)  nCscRechitsChamberPlus32++;
        if (cscRechitsChamber[i] == 41)  nCscRechitsChamberPlus41++;
        if (cscRechitsChamber[i] == 42)  nCscRechitsChamberPlus42++;
        if (cscRechitsChamber[i] == -11)  nCscRechitsChamberMinus11++;
        if (cscRechitsChamber[i] == -12)  nCscRechitsChamberMinus12++;
        if (cscRechitsChamber[i] == -13)  nCscRechitsChamberMinus13++;
        if (cscRechitsChamber[i] == -21)  nCscRechitsChamberMinus21++;
        if (cscRechitsChamber[i] == -22)  nCscRechitsChamberMinus22++;
        if (cscRechitsChamber[i] == -31)  nCscRechitsChamberMinus31++;
        if (cscRechitsChamber[i] == -32)  nCscRechitsChamberMinus32++;
        if (cscRechitsChamber[i] == -41)  nCscRechitsChamberMinus41++;
        if (cscRechitsChamber[i] == -42)  nCscRechitsChamberMinus42++;
      }
      MuonSystem->nCscRings = 0;
      if ( nCscRechitsChamberPlus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberPlus42 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus11 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus12 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus13 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus21 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus22 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus31 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus32 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus41 > 50) MuonSystem->nCscRings++;
      if ( nCscRechitsChamberMinus42 > 50) MuonSystem->nCscRings++;
      //Do DBSCAN Clustering

      int min_point = 50;  //minimum number of Rechitss to call it a cluster
      float epsilon = 0.4; //cluster radius parameter
      //cout<<"cluster points size: "<<points.size()<<endl;
      CACluster ds(min_point, epsilon, points);
      
      ds.run();
      
      //ds.result();

      ds.clusterProperties();

      //cout<<"Num clusters per merge: "<<ds.clusters.size()<<endl;
      //ds.merge_clusters();
      //cout<<"Num clusters post merge: "<<ds.clusters.size()<<endl;

      //ds.sort_clusters();
      
      /*
      ds.clusterProperties();
      ds.merge_clusters();
      ds.sort_clusters();
      */

      //cout<<ds.clusters.size()<<endl;

      MuonSystem->nCscRechitClusters = 0;
      bool event_cluster_matched_timed = false;
      float minClusterTime = -5.0;
      float maxClusterTime=  12.5;
      if (ds.clusters.size()!=1){
	      continue;
      }
      for ( auto &tmp : ds.clusters  ) {
        //cout<<"entering cluster variables"<<endl;
          MuonSystem->cscRechitClusterX[MuonSystem->nCscRechitClusters] =tmp.x;
          MuonSystem->cscRechitClusterY[MuonSystem->nCscRechitClusters] =tmp.y;
          MuonSystem->cscRechitClusterZ[MuonSystem->nCscRechitClusters] =tmp.z;
          MuonSystem->cscRechitClusterTimeWeighted[MuonSystem->nCscRechitClusters] = tmp.tWeighted;
          MuonSystem->cscRechitClusterTimeSpreadWeightedAll[MuonSystem->nCscRechitClusters] = tmp.TSpreadWeightedAll;
          
          MuonSystem->cscRechitClusterTime[MuonSystem->nCscRechitClusters] = tmp.tTotal;
          MuonSystem->cscRechitClusterTimeSpread[MuonSystem->nCscRechitClusters] = tmp.TSpread;

          MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters] =tmp.eta;
          MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters] = tmp.phi;

          MuonSystem->cscRechitClusterSize[MuonSystem->nCscRechitClusters] = tmp.nhits;

          MuonSystem->cscRechitClusterNRechitChamberPlus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus11;
          MuonSystem->cscRechitClusterNRechitChamberPlus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus12;
          MuonSystem->cscRechitClusterNRechitChamberPlus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus13;
          MuonSystem->cscRechitClusterNRechitChamberPlus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus21;
          MuonSystem->cscRechitClusterNRechitChamberPlus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus22;
          MuonSystem->cscRechitClusterNRechitChamberPlus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus31;
          MuonSystem->cscRechitClusterNRechitChamberPlus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus32;
          MuonSystem->cscRechitClusterNRechitChamberPlus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus41;
          MuonSystem->cscRechitClusterNRechitChamberPlus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberPlus42;
          MuonSystem->cscRechitClusterNRechitChamberMinus11[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus11;
          MuonSystem->cscRechitClusterNRechitChamberMinus12[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus12;
          MuonSystem->cscRechitClusterNRechitChamberMinus13[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus13;
          MuonSystem->cscRechitClusterNRechitChamberMinus21[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus21;
          MuonSystem->cscRechitClusterNRechitChamberMinus22[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus22;
          MuonSystem->cscRechitClusterNRechitChamberMinus31[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus31;
          MuonSystem->cscRechitClusterNRechitChamberMinus32[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus32;
          MuonSystem->cscRechitClusterNRechitChamberMinus41[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus41;
          MuonSystem->cscRechitClusterNRechitChamberMinus42[MuonSystem->nCscRechitClusters] = tmp.nCscRechitsChamberMinus42;
          MuonSystem->cscRechitClusterMaxChamber[MuonSystem->nCscRechitClusters] = tmp.maxChamber;
          MuonSystem->cscRechitClusterMaxChamberRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxChamberRechits/tmp.nhits;
          MuonSystem->cscRechitClusterNChamber[MuonSystem->nCscRechitClusters] = tmp.nChamber;
          MuonSystem->cscRechitClusterMaxStation[MuonSystem->nCscRechitClusters] = tmp.maxStation;
          MuonSystem->cscRechitClusterMaxStationRatio[MuonSystem->nCscRechitClusters] = 1.0*tmp.maxStationRechits/tmp.nhits;

          MuonSystem->cscRechitClusterNStation10[MuonSystem->nCscRechitClusters] = tmp.nStation10;
          MuonSystem->cscRechitClusterAvgStation10[MuonSystem->nCscRechitClusters] = tmp.avgStation10;




          /*Comment out jet veto - no jets!
          //Jet veto/ muon veto
          MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] = 0.0;
          MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer_TnP::deltaR(jetEta[i], jetPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] ) {
              MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = jetPt[i];
              MuonSystem->cscRechitClusterJetVetoE[MuonSystem->nCscRechitClusters]  = jetE[i];
              MuonSystem->cscRechitClusterJetVetoTightId[MuonSystem->nCscRechitClusters]  = jetPassIDTight[i];
              MuonSystem->cscRechitClusterJetVetoLooseId[MuonSystem->nCscRechitClusters]  = jetPassIDLoose[i];

            }

          }
          */
         
          float min_deltaR = 15.;
          int index = 999;
          bool matchedSingleCluster = false;
          bool matchedSingleCluster_notProbe = false;
          bool matchedSingleClusterHighPtJet = false;
          bool matchedSingleClusterLowPtJet = false;
          bool passTimeSingleCluster = false;
          bool noHits_Me1112_SingleCluster = false;
          //cout<<"here"<<endl;

          //flag clusters with no hits in ME11/12
          if (tmp.nCscRechitsChamberMinus11 == 0 && tmp.nCscRechitsChamberMinus12 == 0 && tmp.nCscRechitsChamberPlus11 == 0 && tmp.nCscRechitsChamberPlus12 == 0) noHits_Me1112_SingleCluster = true;
          MuonSystem->cscRechitCluster_passME1112Veto[MuonSystem->nCscRechitClusters] = noHits_Me1112_SingleCluster;
          float minDeltaR = 100000; float currentDeltaR = 0;
          for(auto &tmp : Leptons_Probe )
          { 
            //if (!MuonSystem->lepTag[i]) continue; //first find tagged muon
            //for (int j = 0; j < Leptons.size(); j++){
              //if (i==j) continue; //skip if same muon
              currentDeltaR = RazorAnalyzerMerged::deltaR(tmp.lepton.Eta(), tmp.lepton.Phi(), MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
              if (currentDeltaR < 0.4 && tmp.lepton.Pt() > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] && currentDeltaR < minDeltaR) {
                minDeltaR = currentDeltaR; //this way, we always match to muon with smallest delta R less than 0.4
                MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = tmp.lepton.Pt();
                //MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[j];
                MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = tmp.isGlobal;
                MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = tmp.passLooseId;
                matchedSingleCluster = true;
                }
            }
          
          MuonSystem->cscRechitCluster_matchToProbeMuon[MuonSystem->nCscRechitClusters] = matchedSingleCluster;
          if (!matchedSingleCluster){
          for(auto &tmp : Leptons_notProbe )
          { 
            //if (!MuonSystem->lepTag[i]) continue; //first find tagged muon
            //for (int j = 0; j < Leptons.size(); j++){
              //if (i==j) continue; //skip if same muon
              currentDeltaR = RazorAnalyzerMerged::deltaR(tmp.lepton.Eta(), tmp.lepton.Phi(), MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
              if (currentDeltaR < 0.4 && tmp.lepton.Pt() > MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters] && currentDeltaR < minDeltaR) {
                minDeltaR = currentDeltaR; //this way, we always match to muon with smallest delta R less than 0.4
                MuonSystem->cscRechitClusterMuonVetoPt[MuonSystem->nCscRechitClusters]  = tmp.lepton.Pt();
                //MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[j];
                MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = tmp.isGlobal;
                MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = tmp.passLooseId;
                matchedSingleCluster_notProbe = true;
                }
            }
          }
            MuonSystem->cscRechitCluster_matchToNotProbeMuon[MuonSystem->nCscRechitClusters] = matchedSingleCluster_notProbe;
          MuonSystem->cscRechitCluster_matchToProbeAndJet[MuonSystem->nCscRechitClusters] = false;
          float minDeltaRJet = 100000; float currentDeltaRJet = 0;
          for(auto &tmp : Jets )
          { 
            //if (!MuonSystem->lepTag[i]) continue; //first find tagged muon
            //for (int j = 0; j < Leptons.size(); j++){
              //if (i==j) continue; //skip if same muon
              currentDeltaR = RazorAnalyzerMerged::deltaR(tmp.jet.Eta(), tmp.jet.Phi(), MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]);
              if (currentDeltaRJet < 0.4 && tmp.jet.Pt() > MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters] && currentDeltaRJet < minDeltaRJet) {
                if (matchedSingleCluster||matchedSingleCluster_notProbe) {
                  MuonSystem->cscRechitCluster_matchToProbeAndJet[MuonSystem->nCscRechitClusters] = true;
                  continue;
                }
                
                minDeltaRJet = currentDeltaRJet; //this way, we always match to muon with smallest delta R less than 0.4
                MuonSystem->cscRechitClusterJetVetoPt[MuonSystem->nCscRechitClusters]  = tmp.jet.Pt();
                //MuonSystem->cscRechitClusterMuonVetoE[MuonSystem->nCscRechitClusters]  = muonE[j];
                //MuonSystem->cscRechitClusterMuonVetoGlobal[MuonSystem->nCscRechitClusters]  = Muon_isGlobal[i];
                //MuonSystem->cscRechitClusterMuonVetoLooseId[MuonSystem->nCscRechitClusters]  = Muon_looseId[i];
                if (tmp.highPtJet) matchedSingleClusterHighPtJet = true;
                else matchedSingleClusterLowPtJet = true;
                
            }
          }
            MuonSystem->cscRechitCluster_matchToHighPtJet[MuonSystem->nCscRechitClusters] = matchedSingleClusterHighPtJet;
            MuonSystem->cscRechitCluster_matchToLowPtJet[MuonSystem->nCscRechitClusters] = matchedSingleClusterLowPtJet;
            MuonSystem->cscRechitCluster_notMatched[MuonSystem->nCscRechitClusters] = !(matchedSingleCluster||matchedSingleCluster_notProbe||matchedSingleClusterHighPtJet||matchedSingleClusterLowPtJet);
          if(!isData)
          {
            // match to gen level LLP
            min_deltaR = 15.;
            index = 999;
            for(int j = 0; j < MuonSystem->nGLLP;j++)
            {

              double current_delta_r = RazorAnalyzerMerged::deltaR(MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters], MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters], MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = true;
            else MuonSystem->cscRechitCluster_match_gLLP[MuonSystem->nCscRechitClusters] = false;


            MuonSystem->cscRechitCluster_match_gLLP_minDeltaR[MuonSystem->nCscRechitClusters] = min_deltaR;
            MuonSystem->cscRechitCluster_match_gLLP_index[MuonSystem->nCscRechitClusters] = index;
            MuonSystem->cscRechitCluster_match_gLLP_eta[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_eta[index];
            MuonSystem->cscRechitCluster_match_gLLP_phi[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_phi[index];
            MuonSystem->cscRechitCluster_match_gLLP_decay_r[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
            MuonSystem->cscRechitCluster_match_gLLP_decay_z[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
            MuonSystem->cscRechitCluster_match_gLLP_csc[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_csc[index];
            MuonSystem->cscRechitCluster_match_gLLP_dt[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_dt[index];
            MuonSystem->cscRechitCluster_match_gLLP_e[MuonSystem->nCscRechitClusters] = MuonSystem->gLLP_e[index];
          }
          
          //match to MB1 DT segments
          /*
          for (int i = 0; i < nDtSeg; i++) {
            if (RazorAnalyzer_trigEff::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              MuonSystem->cscRechitCluster_match_dtSeg_0p4[MuonSystem->nCscRechitClusters] ++;
              if (dtSegStation[i] == 1) MuonSystem->cscRechitCluster_match_MB1Seg_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }

          //match to RPC hits in RE1/2
          for (int i = 0; i < nRpc; i++) {
            float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
            if (RazorAnalyzer_trigEff::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->cscRechitClusterEta[MuonSystem->nCscRechitClusters],MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters]) < 0.4 )
            {
              if (rpcR < 461.0 && rpcR > 275 && abs(rpcZ[i]) > 663 && abs(rpcZ[i]) < 730) MuonSystem->cscRechitCluster_match_RE12_0p4[MuonSystem->nCscRechitClusters] ++;
              if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661) MuonSystem->cscRechitCluster_match_RB1_0p4[MuonSystem->nCscRechitClusters] ++;
            }
          }
          */
          MuonSystem->cscRechitClusterMet_dPhi[MuonSystem->nCscRechitClusters] =  RazorAnalyzerMerged::deltaPhi(MuonSystem->cscRechitClusterPhi[MuonSystem->nCscRechitClusters],MuonSystem->metPhi);

          //check if there is at least one cluster that passes the time veto
          if (tmp.tTotal <= maxClusterTime && tmp.tTotal >= minClusterTime) {
             passTimeSingleCluster = true;
          }
          MuonSystem->cscRechitCluster_PassTimeVeto[MuonSystem->nCscRechitClusters] = passTimeSingleCluster;

          if (passTimeSingleCluster && matchedSingleCluster && noHits_Me1112_SingleCluster) event_cluster_matched_timed = true;

      //Now, compute HLT Decision manually
          //We'll use this to compute the expected HLT decision, compare with the actual trigger decision, and determine efficiency for
          //events with multiple clusters
          int eta_nStation_bin = 0;
          if (fabs(tmp.eta)<1.9 && tmp.nStation10 <= 1){
            eta_nStation_bin = 0;
          }
          if (fabs(tmp.eta)<1.9 && tmp.nStation10 > 1){
            eta_nStation_bin = 1;
          }
          if (fabs(tmp.eta)>=1.9 && tmp.nStation10 <= 1){
            eta_nStation_bin = 2;
          }
          if (fabs(tmp.eta)>=1.9 && tmp.nStation10 > 1){
            eta_nStation_bin = 3;
          }

          int loose_rechit_thresh = 10000; int medium_rechit_thresh = 10000; int tight_rechit_thresh = 10000;
          switch (eta_nStation_bin){
            case 0:
              loose_rechit_thresh = 200;
              medium_rechit_thresh = 300;
              tight_rechit_thresh = 500;
              break;

            case 1:
              loose_rechit_thresh = 100;
              medium_rechit_thresh = 100;
              tight_rechit_thresh = 100;
              break;

            case 2:
              loose_rechit_thresh = 500;
              medium_rechit_thresh = 800;
              tight_rechit_thresh = 800;
              break;

            case 3:
              loose_rechit_thresh = 500;
              medium_rechit_thresh = 500;
              tight_rechit_thresh = 500;
              break;
          }
          
          MuonSystem->cscRechitCluster_HLTCscCluster_Loose_Decision[MuonSystem->nCscRechitClusters] = (tmp.nhits >= loose_rechit_thresh && passTimeSingleCluster && noHits_Me1112_SingleCluster);
          MuonSystem->cscRechitCluster_HLTCscCluster_Medium_Decision[MuonSystem->nCscRechitClusters] = (tmp.nhits >= medium_rechit_thresh && passTimeSingleCluster && noHits_Me1112_SingleCluster);
          MuonSystem->cscRechitCluster_HLTCscCluster_Tight_Decision[MuonSystem->nCscRechitClusters] = (tmp.nhits >= tight_rechit_thresh && passTimeSingleCluster && noHits_Me1112_SingleCluster);
          MuonSystem->nCscRechitClusters++;
      }
      //if (!event_cluster_matched_timed) continue;
      total_events_passed++; 
      // DT cluster

      points.clear();
      // cout<<"here"<<endl;
      /*
      for (int i = 0; i < nDtRechits; i++) {
        Rechits p;

        p.phi = dtRechitCorrectPhi[i];
        p.eta = dtRechitCorrectEta[i];
        p.x = dtRechitCorrectX[i];
        p.y = dtRechitCorrectY[i];
        p.z = dtRechitCorrectZ[i];
        p.t = dtRechitTime[i];
        p.twire = dtRechitTime[i];
        p.station = dtRechitStation[i];
        p.chamber = dtRechitWheel[i];
        p.superlayer = dtRechitSuperLayer[i];
        p.clusterID = UNCLASSIFIED;
        points.push_back(p);

      }
      // cout<<"here"<<endl;

      //Do DBSCAN Clustering
      int min_point_dt = 50;  //minimum number of segments to call it a cluster
      float epsilon_dt = 0.2; //cluster radius parameter
      CACluster ds_dtRechit(min_point_dt, epsilon_dt, points);
      
      ds_dtRechit.run();
      
      //ds.result();

      ds_dtRechit.clusterProperties();

      //cout<<"Num clusters per merge: "<<ds.clusters.size()<<endl;
      //ds_dtRechit.merge_clusters();
      //cout<<"Num clusters post merge: "<<ds.clusters.size()<<endl;

      //ds_dtRechit.sort_clusters();

      MuonSystem->nDtRechitClusters = 0;

      for ( auto &tmp : ds_dtRechit.clusters  ) {
         //cout<<"entering dt cluster variables"<<endl;
        //remove overlaps
        bool overlap = false;
        for(int i = 0; i < MuonSystem->nCscRechitClusters; i++)
        {
          if (RazorAnalyzer_trigEff::deltaR(MuonSystem->cscRechitClusterEta[i],MuonSystem->cscRechitClusterPhi[i],tmp.eta, tmp.phi)<0.4) overlap = true;
        }
        if (overlap) continue;


          MuonSystem->dtRechitClusterX[MuonSystem->nDtRechitClusters] =tmp.x;
          MuonSystem->dtRechitClusterY[MuonSystem->nDtRechitClusters] =tmp.y;
          MuonSystem->dtRechitClusterZ[MuonSystem->nDtRechitClusters] =tmp.z;
          if (abs(tmp.z) < 126.8) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 0;
          else if (tmp.z > 126.8 && tmp.z < 395.4) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 1;
          else if (tmp.z < -126.8 && tmp.z > -395.4)MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -1;
          else if (tmp.z<0) MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = -2;
          else MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] = 2;
          MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters] =tmp.eta;
          MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters] =tmp.phi;
          MuonSystem->dtRechitClusterSize[MuonSystem->nDtRechitClusters] = tmp.nhits;

          unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
          default_random_engine generator (seed);

	        // default_random_engine generator;
    	   uniform_real_distribution<double> distribution(0.0,1.0);
	       float prob = 0.03;
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters]++;
         }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<12; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters]++;
      	 }
         for (int i=0; i<8; ++i) {
           if ( distribution(generator) < prob) MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters]++;
      	 }

         MuonSystem->dtRechitClusterNoiseHit[MuonSystem->nDtRechitClusters] = MuonSystem->dtRechitClusterNoiseHitStation1[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation2[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation3[MuonSystem->nDtRechitClusters] +
                                                                              MuonSystem->dtRechitClusterNoiseHitStation4[MuonSystem->nDtRechitClusters];


          MuonSystem->dtRechitClusterNHitStation1[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation1;
        	MuonSystem->dtRechitClusterNHitStation2[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation2;
        	MuonSystem->dtRechitClusterNHitStation3[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation3;
        	MuonSystem->dtRechitClusterNHitStation4[MuonSystem->nDtRechitClusters] = tmp.nDtRechitsStation4;

        	MuonSystem->dtRechitClusterMaxChamber[MuonSystem->nDtRechitClusters] = tmp.maxChamber;
        	MuonSystem->dtRechitClusterMaxChamberRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxChamberRechits/tmp.nhits;
        	MuonSystem->dtRechitClusterNChamber[MuonSystem->nDtRechitClusters] = tmp.nChamber;
        	MuonSystem->dtRechitClusterMaxStation[MuonSystem->nDtRechitClusters] = tmp.maxStation;
        	MuonSystem->dtRechitClusterMaxStationRatio[MuonSystem->nDtRechitClusters] = 1.0*tmp.maxStationRechits/tmp.nhits;
          MuonSystem->dtRechitClusterNStation10[MuonSystem->nDtRechitClusters] = tmp.nStation10;
          MuonSystem->dtRechitClusterAvgStation10[MuonSystem->nDtRechitClusters] = tmp.avgStation10;
          */

          /*Comment out jet veto - no jets!
          //Jet veto/ muon veto
          MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] = 0.0;
          MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters] = 0.0;


          // jet veto
          for(int i = 0; i < nJets; i++)
          {
            if (fabs(jetEta[i]>3.0)) continue;
            if (RazorAnalyzer_TnP::deltaR(jetEta[i], jetPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && jetPt[i] > MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterJetVetoPt[MuonSystem->nDtRechitClusters]  = jetPt[i];
              MuonSystem->dtRechitClusterJetVetoE[MuonSystem->nDtRechitClusters]  = jetE[i];
              MuonSystem->dtRechitClusterJetVetoLooseId[MuonSystem->nDtRechitClusters]  = jetPassIDLoose[i];
              MuonSystem->dtRechitClusterJetVetoTightId[MuonSystem->nDtRechitClusters]  = jetPassIDTight[i];

            }


          }
          */
         /*
          for(int i = 0; i < nMuon; i++)
          {
            if (fabs(Muon_eta[i]>3.0)) continue;
            //float muonIso = (muon_chargedIso[i] + fmax(0.0,  muon_photonIso[i] + muon_neutralHadIso[i] - 0.5*muon_pileupIso[i])) / Muon_pt[i];
            float muonIso = Muon_pfRelIso04_all[i]; //added by alex, address set to Muon_pfRelIso04_all branch orignally fron nanoAODs
            if (RazorAnalyzer_trigEff::deltaR(Muon_eta[i], Muon_phi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 && Muon_pt[i] > MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters] ) {
              MuonSystem->dtRechitClusterMuonVetoPt[MuonSystem->nDtRechitClusters]  = Muon_pt[i];
              //MuonSystem->dtRechitClusterMuonVetoE[MuonSystem->nDtRechitClusters]  = muonE[i];
              MuonSystem->dtRechitClusterMuonVetoGlobal[MuonSystem->nDtRechitClusters]  = Muon_isGlobal[i];
              MuonSystem->dtRechitClusterMuonVetoLooseId[MuonSystem->nDtRechitClusters]  = Muon_looseId[i];


            }
          }
          
          for (int i = 0; i < nDtSeg; i++) {
              if (RazorAnalyzer_trigEff::deltaR(dtSegEta[i], dtSegPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
              if (abs(RazorAnalyzer_trigEff::deltaPhi(dtSegPhi[i],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]))>2) {
                if (dtSegStation[i] == 1) MuonSystem->dtRechitClusterNOppositeSegStation1[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 2) MuonSystem->dtRechitClusterNOppositeSegStation2[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 3) MuonSystem->dtRechitClusterNOppositeSegStation3[MuonSystem->nDtRechitClusters]  +=1;
                if (dtSegStation[i] == 4) MuonSystem->dtRechitClusterNOppositeSegStation4[MuonSystem->nDtRechitClusters]  +=1;
              }
         }
          


          // match to gen-level LLP
          float min_deltaR = 15.;
          int index = 999;
          if (!isData)
          {
            for(int j = 0; j < MuonSystem->nGLLP;j++)
            {
              double current_delta_r = RazorAnalyzer_trigEff::deltaR(tmp.eta, tmp.phi, MuonSystem->gLLP_eta[j], MuonSystem->gLLP_phi[j]);
              if (current_delta_r < min_deltaR)
              {
                min_deltaR = current_delta_r;
                index = j;
              }
            }
            if (min_deltaR < 0.4)MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = true;
            else MuonSystem->dtRechitCluster_match_gLLP[MuonSystem->nDtRechitClusters] = false;

             MuonSystem->dtRechitCluster_match_gLLP_minDeltaR[MuonSystem->nDtRechitClusters] = min_deltaR;
             MuonSystem->dtRechitCluster_match_gLLP_index[MuonSystem->nDtRechitClusters] = index;
             MuonSystem->dtRechitCluster_match_gLLP_eta[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_eta[index];
             MuonSystem->dtRechitCluster_match_gLLP_phi[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_phi[index];
             MuonSystem->dtRechitCluster_match_gLLP_decay_r[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_r[index];
             MuonSystem->dtRechitCluster_match_gLLP_decay_z[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_decay_vertex_z[index];
             MuonSystem->dtRechitCluster_match_gLLP_csc[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_csc[index];
             MuonSystem->dtRechitCluster_match_gLLP_dt[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_dt[index];
             MuonSystem->dtRechitCluster_match_gLLP_e[MuonSystem->nDtRechitClusters] = MuonSystem->gLLP_e[index];

          }


          //match to MB1 DT segments
          MuonSystem->nCscRechits = ncscRechits;

          for (int i = 0; i < nDtRechits; i++) {
            if (RazorAnalyzer_trigEff::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.5 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p5[MuonSystem->nDtRechitClusters] ++;
            }
            if (RazorAnalyzer_trigEff::deltaR(dtRechitCorrectEta[i], dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters],MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
            {
              if (dtRechitStation[i] == 1) MuonSystem->dtRechitCluster_match_MB1hits_0p4[MuonSystem->nDtRechitClusters] ++;
            }
            if(abs(dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])==1 && dtRechitStation[i] == 1)
            {
              if (abs(RazorAnalyzer_trigEff::deltaPhi(dtRechitCorrectPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < TMath::Pi()/4.0 )
              {
                if (dtRechitWheel[i]-MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters] == 1) MuonSystem->dtRechitCluster_match_MB1hits_cosmics_plus[MuonSystem->nDtRechitClusters] ++;
                else MuonSystem->dtRechitCluster_match_MB1hits_cosmics_minus[MuonSystem->nDtRechitClusters] ++;
              }
            }

          }

         std::vector<int> dtRechitCluster_match_rpcBx;

         //match to RPC hits with dPhi<0.5 and same wheel in DT
         for (int i = 0; i < nRpc; i++) {
           float rpcR = sqrt(rpcX[i]*rpcX[i] + rpcY[i]*rpcY[i]);
           if (rpcRegion[i]!=0) continue;
           if (abs(RazorAnalyzer_trigEff::deltaPhi(rpcPhi[i], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters])) < 0.5 )
           {
             if (rpcRing[i] == MuonSystem->dtRechitClusterWheel[MuonSystem->nDtRechitClusters])
             {
               dtRechitCluster_match_rpcBx.push_back(rpcBx[i]);
               MuonSystem->dtRechitCluster_match_RPChits_dPhi0p5[MuonSystem->nDtRechitClusters]++;
               if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_dPhi0p5[MuonSystem->nDtRechitClusters] ++;

             }
           }
           if(RazorAnalyzer_trigEff::deltaR(rpcEta[i], rpcPhi[i], MuonSystem->dtRechitClusterEta[MuonSystem->nDtRechitClusters], MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters]) < 0.4 )
           {
             if (rpcR < 470 && rpcR > 380 && abs(rpcZ[i]) < 661)MuonSystem->dtRechitCluster_match_RB1_0p4[MuonSystem->nDtRechitClusters] ++;
           }
         }
         int max_occurence = 0;
         int max_bx = -999;
         for (unsigned int l = 0; l < dtRechitCluster_match_rpcBx.size(); l++)
         {
           int counter = 0;
           for(unsigned int j = 0; j < dtRechitCluster_match_rpcBx.size(); j ++)
           {
             if (dtRechitCluster_match_rpcBx[j] == dtRechitCluster_match_rpcBx[l]) counter++;
           }
           if (counter>max_occurence)
           {
             max_occurence = counter;
             max_bx = dtRechitCluster_match_rpcBx[l];
           }
         }
          MuonSystem->dtRechitCluster_match_RPCBx_dPhi0p5[MuonSystem->nDtRechitClusters] = max_bx;

          MuonSystem->dtRechitClusterMet_dPhi[MuonSystem->nDtRechitClusters] =  RazorAnalyzer_trigEff::deltaPhi(MuonSystem->dtRechitClusterPhi[MuonSystem->nDtRechitClusters],MuonSystem->metPhi);

          MuonSystem->nDtRechitClusters++;
          
        }
    
      // if (isData && MuonSystem->nDtRechitClusters + MuonSystem->nCscRechitClusters < 2) continue;
      */
      if(!isData && signalScan)
      {
        pair<int,int> smsPair = make_pair(MuonSystem->mX, MuonSystem->ctau);
        Trees2D[smsPair]->Fill();
      }
      else
      {
        //cout << "filling" << endl;
        MuonSystem->tree_->Fill();
      }


    }
      if(!isData && signalScan)
      {
        for(auto &filePtr : Files2D)
         {
           cout << "Writing output tree (" << filePtr.second->GetName() << ")" << endl;
           filePtr.second->cd();
           Trees2D[filePtr.first]->Write();
           NEvents2D[filePtr.first]->Write("NEvents");
           Total2D[filePtr.first]->Write("Total");
           accep2D[filePtr.first]->Write("acceptance");
           accep_met2D[filePtr.first]->Write("acceptance_met");
           filePtr.second->Close();

         }
      }
      else if (!isData)
      {
         cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
         cout << "Writing output trees..." << endl;
         outFile->cd();
         MuonSystem->tree_->Write();
         NEvents->Write();
         accep->Write("acceptance");
         accep_met->Write("acceptance_met");
         outFile->Close();
      }


      else
      {
        cout << "Filled Total of " << NEvents->GetBinContent(1) << " Events\n";
        cout << "Writing output trees..." << endl;
        outFile->cd();
        MuonSystem->tree_->Write();
        Nmet200->Write();
        NmetFilter->Write();
        Nlep0->Write();
        Njet1->Write();
        NcosmicVeto->Write();
        NEvents->Write();
        // outFile->Write();
        outFile->Close();
      }
      cout<<"Acceptances: "<<endl;
      cout<<"Total: "<<total_events_passed/maxEvents<<endl;
      cout<<"Two Muons Passing Probe with Oppose Charges: "<<two_satsifying_probe/maxEvents*100<<"%"<<endl;
      cout<<"At Least One Tag: "<<greater_than_one_tag/maxEvents*100<<"%"<<endl;
      cout<<"Leading Muon pT>35: "<<leading_muon_pt/maxEvents*100<<"%"<<endl;
      cout<<"Subleading Muon pT>20: "<<subleading_muon_pt/maxEvents*100<<"%"<<endl;
      cout<<"Total Efficiency - TrigEff Events: "<<Zmumu_events_passed/maxEvents*100<<"%"<<endl;


}
