// Author:  Ben Tannenwald
// Date:    Oct 30, 2018
// Purpose: Macro to analyze ttbar sample

#include <fstream>
#include <TFile.h>
#include <TDirectory.h>
#include <TTree.h>
#include <TBranch.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TH1.h>
#include <TLorentzVector.h>

#include "YggdrasilEventVars.h"

using namespace std;

void printProgBar( int percent );

void ttbar3()
{
  // *** 0. Set input and output files, set histogram output style
  // ** A. Set input file
  TChain* fChain = new TChain("ttHTreeMaker/worldTree");
  fChain->AddFile("$PWD/yggdrasil_treeMaker_ttbar_2l2nu.root");
  // ** B. Set output file
  TFile *outfile = new TFile( "$PWD/ttbar3.root", "RECREATE");
  //fstream output;
  //output.open("Event Categories.txt",ios_base::in|ios_base::out|ios_base::trunc);
    
  
  // *** 1. Define canvasses, histograms, and 'global' variables
  // ** A. Canvas
  //TCanvas *c0 = new TCanvas("c0", "c0", 50, 50, 800, 600);
  // ** B. Histograms
  TH1D* h_jet_n = new TH1D("h_jet_n", "Number of jets with cuts", 12, 0, 12);//jet distribution
  TH1D *h_Muon_n = new TH1D("h_Muon_n","Number of Muons with cuts",10,0,10);//Muon distribution
  TH1D *h_Electron_n = new TH1D("h_Electron_n", "Number of electrons with cuts", 10,0,10);//electron distribution
  TH1D *h_cat1_njets = new TH1D("h_cat1_njets","Distribution of jets from category1 (ee)",10,0,10);
  TH1D *h_cat2_njets = new TH1D("h_cat2_njets","Distribution of jets from category2 (emu)",10,0,10);
  TH1D *h_cat3_njets = new TH1D("h_cat3_njets","Distribution of jets from category3 (mumu)",10,0,10);
  TH1D *h_cat1_mll = new TH1D("h_cat1_mll","m_ll from category1 (ee)",500,0,500);
  TH1D *h_cat2_mll = new TH1D("h_cat2_mll","m_ll from category2 (emu)",500,0,500);
  TH1D *h_cat3_mll = new TH1D("h_cat3_mll","m_ll from category3 (mumu)",500,0,500);
  TH1D *h_sum_njets = new TH1D("h_sum_njets","Total distribution of jets from all 3 categories",10,0,10);
  TH1D *h_sum_mll = new TH1D("h_sum_mll","Total distribution of dilepton invariant mass from all 3 categories",500,0,500);
  // ** C. Variables
  int nJets;
  int nMuons;
  int nMuPt,nMuEta;
  int neCut;
  int jetCount;
  int lepCount;
  int nE = 0; //Number of electrons without cuts
  int nM = 0; //Number of Muons without cuts
  int nJ = 0; //Number of jets with pT>30GeV
  int nCat1 = 0;
  int nCat2 = 0;
  int nCat3 = 0;
  double jetCom = 0;
  double lepCom = 0;
  




  // *** 3. Set tree structure
  yggdrasilEventVars *eve;
  eve=0;
  //eventTree->SetBranchAddress("eve.", &eve );
  fChain->SetBranchAddress("eve.", &eve );
  
  // *** 4. Start looping! 
  long t_entries = fChain->GetEntries();
  //long t_entries = 2500; // for tests using small subset of events
  cout << "Chain entries: " << t_entries << endl;
    
  for(int i = 0; i < t_entries; i++) {
    // Keep track of progress through file
    if ( t_entries > 100) {
      if ((i+1)%(5*t_entries/100)==0)  printProgBar(100*i/t_entries +1); }
    if (i == t_entries-1) {printProgBar(100); cout << endl;}
    
    fChain->GetEntry(i);

    // re-iniialize event-level variables
    nJets = 0;
    nMuons = 0;
    nMuPt = nMuEta = 0;
    jetCom = lepCom = jetCount = lepCount = 0;
    neCut = 0;
    nE = nM = nJ = 0;


    double Echarge[20] = {0};
    double Mcharge[20] = {0};
    double EpT[20] = {0};
    double MpT[20] = {0};
    double jetpT[20] = {0};
    double Eeta[20] = {0};
    double Ephi[20] = {0};
    double Ee[20] = {0};
    double Meta[20] = {0};
    double Mphi[20] = {0};
    double Me[20] = {0};
    TLorentzVector p1;                  //lorentz vector for lepton 1 p1
    TLorentzVector p2;                  //lorentz vector for lepton 2 p2
    TLorentzVector ptotal;              //sum of the p1 and p2: (p1+p2);

    double mll = 0;                     //invariant mass for dilepton
    
    // ** I. Loop over jets
    //std::cout << "Evt: " << eve->evt_ << ", njets: " << eve->jet_pt_[0].size() << ", nLep: " << eve->lepton_pt_.size() << std::endl;
    for (unsigned int j = 0; j < eve->jet_pt_[0].size(); j++){
      bool JetSatisfyCut=0; 

      //calc nJets
      if ( (eve->jet_pt_[0][j] > 20)&&(TMath::Abs(eve->jet_eta_[0][j])<2.4)&&(eve->jet_puid_[0][j]==7) ){
        JetSatisfyCut = 1;
        if(eve->jet_pt_[0][j]>30)   nJ++;
          for(int k=0;k<eve->lepton_pt_.size();k++)
          {
              if((eve->lepton_isMuon_[k])&&((eve->lepton_pt_[k])>15)&&(TMath::Abs(eve->lepton_eta_[k])<2.4)&&(eve->lepton_relIso_[k]<0.25)&&(eve->lepton_isTight_[k]==1))
              {
                  if(sqrt((eve->jet_eta_[0][j]-eve->lepton_eta_[k])*(eve->jet_eta_[0][j]-eve->lepton_eta_[k])+(eve->jet_phi_[0][j]-eve->lepton_phi_[k])*(eve->jet_phi_[0][j]-eve->lepton_phi_[k]))<=0.4)
                  JetSatisfyCut = 0;
                  break;
              }
              else if((eve->lepton_isMuon_[k]==0)&&((eve->lepton_pt_[k])>15)&&(TMath::Abs(eve->lepton_eta_[k])<2.4)&&((TMath::Abs(eve->lepton_scEta_[k])<1.4442)||(TMath::Abs(eve->lepton_scEta_[k])>1.5660))&&(eve->lepton_isTight_[k]==1))
              {
                  if(sqrt((eve->jet_eta_[0][j]-eve->lepton_eta_[k])*(eve->jet_eta_[0][j]-eve->lepton_eta_[k])+(eve->jet_phi_[0][j]-eve->lepton_phi_[k])*(eve->jet_phi_[0][j]-eve->lepton_phi_[k]))<=0.4)
                  JetSatisfyCut = 0;
                  break;
              }
          }
	    //nJets++;
      }
      if(JetSatisfyCut) nJets++;
    }
    h_jet_n->Fill(nJets);

    // ** I. Loop over leptons
    for (unsigned int l=0; l < eve->lepton_pt_.size(); l++){
      if ( (eve->lepton_isMuon_[l]) ){                      //lepton is Muon
        nM++;
        Mcharge[nM-1] = eve->lepton_charge_[l];
        MpT[nM-1] = eve->lepton_pt_[l];
        Meta[nM-1] = eve->lepton_eta_[l];
        Mphi[nM-1] = eve->lepton_phi_[l];
        Me[nM-1] = eve->lepton_e_[l];
        if(((eve->lepton_pt_[l])>15)&&(TMath::Abs(eve->lepton_eta_[l])<2.4)&&(eve->lepton_relIso_[l]<0.25)&&(eve->lepton_isTight_[l]==1))    nMuons++;
      } 
      if(eve->lepton_isMuon_[l]==0){//electrons
        nE++;
        Echarge[nE-1] = eve->lepton_charge_[l];
        EpT[nE-1] = eve->lepton_pt_[l];
        Eeta[nE-1] = eve->lepton_eta_[l];
        Ephi[nE-1] = eve->lepton_phi_[l];
        Ee[nE-1] = eve->lepton_e_[l];
        if(((eve->lepton_pt_[l])>15)&&(TMath::Abs(eve->lepton_eta_[l])<2.4)&&((TMath::Abs(eve->lepton_scEta_[l])<1.4442)||(TMath::Abs(eve->lepton_scEta_[l])>1.5660))&&(eve->lepton_isTight_[l]==1))    neCut++;
      }
    }
    h_Muon_n->Fill(nMuons);
    h_Electron_n->Fill(neCut);

    //fChain->Draw("jet_pt_[0]");
    /*TH1F *h_jet_pt = (TH1F*)gPad->GetPrimitive("h_jet_pt");
    fChain->Draw("h_jet_n:jet_eta_[0]","jet_pt_[0]>20");
    TH1F *h_jet_eta = (TH1F*)gPad->GetPrimitive("h_jet_eta");
    fChain->Draw("h_jet_n:jet_phi_[0]","jet_pt_[0]>20");
    TH1F *h_jet_phi = (TH1F*)gPad->GetPrimitive("h_jet_phi");*/



    //**** do the categories!
    if((nM==0)&&(nE==2)&&(Echarge[0]+Echarge[1]==0)&&((EpT[0]>25)||(EpT[1]>25))&&(nJ>=2)&&(eve->MET_Type1xy_[0]>40)){
        h_cat1_njets->Fill(nJets);
        p1.SetPtEtaPhiE(EpT[0],Eeta[0],Ephi[0],Ee[0]);
        p2.SetPtEtaPhiE(EpT[1],Eeta[1],Ephi[1],Ee[1]);
        ptotal = p1+p2;
        mll = ptotal.M();
        h_cat1_mll->Fill(mll);
    //Cat1Index[nCat1] = i;
    //nCat1++;
    }
    else if((nM==1)&&(nE==1)&&(Echarge[0]+Mcharge[0]==0)&&((EpT[0]>25)||(MpT[0]>25))&&(nJ>=2)){
        h_cat2_njets->Fill(nJets);
        p1.SetPtEtaPhiE(EpT[0],Eeta[0],Ephi[0],Ee[0]);
        p2.SetPtEtaPhiE(MpT[0],Meta[0],Mphi[0],Me[0]);
        ptotal = p1+p2;
        mll = ptotal.M();
        h_cat2_mll->Fill(mll);
    //Cat2Index[nCat2] = i;
    //nCat2++;
    }
    else if((nM==2)&&(nE==0)&&(Mcharge[0]+Mcharge[1]==0)&&((MpT[0]>25)||(MpT[1]>25))&&(nJ>=2)&&(eve->MET_Type1xy_[0]>40)){
        h_cat3_njets->Fill(nJets);
        p1.SetPtEtaPhiE(MpT[0],Meta[0],Mphi[0],Me[0]);
        p2.SetPtEtaPhiE(MpT[1],Meta[1],Mphi[1],Me[1]);
        ptotal = p1+p2;
        mll = ptotal.M();
        h_cat3_mll->Fill(mll);
    //Cat3Index[nCat3] = i;
    //nCat3++;
    }
    //std::cout << "Now we are in entry:"<<i<<std::endl;
    //std::cout << "Now the numbers are: "<<nCat1<<"   "<<nCat2<<"   "<<nCat3<<std::endl;
  } // end of event loop

  //output << "Cat1" << "  " << "Cat2" << "   " << "Cat3" << endl;
  h_sum_njets->Add(h_cat1_njets);
  h_sum_njets->Add(h_cat2_njets);
  h_sum_njets->Add(h_cat3_njets);

  h_sum_mll->Add(h_cat1_mll);
  h_sum_mll->Add(h_cat2_mll);
  h_sum_mll->Add(h_cat3_mll);




  // *** 4. Save plots to output file
  //output.close();
  outfile->cd();
  outfile->Write();
  
  return 0;
}
    

void printProgBar( int percent )
{
  string bar;
  
  for(int i = 0; i < 50; i++){
    if( i < (percent/2)){
      bar.replace(i,1,"=");
    }else if( i == (percent/2)){
      bar.replace(i,1,">");
    }else{
      bar.replace(i,1," ");
    }
  }
  
  cout<< "\r" "[" << bar << "] ";
  cout.width( 3 );
  cout<< percent << "%     " << std::flush;
}

