#define _USE_MATH_DEFINES
#include "TEfficiency.h"
#include "TClonesArray.h"
#include "math.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TStyle.h"
#include "TCanvas.h"
#include <TMVA/ROCCurve.h>
#include "TBranch.h"
#include "TBasket.h"
#include "TLorentzVector.h"
//#include <TDatabasePDG.h>
#include <TGraphAsymmErrors.h>
#include "TH1.h"
#include "TH2.h"
#include <vector>
#include <iostream>
#include <string>
#include "CB.C"
#include "tdrstyle.C"

using namespace std;

 std::ostringstream NameString;
  std::ostringstream NameString2;
 std::unordered_map<std::string, TH1F*> hist;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////analysis for AK15 AK8 AK4 Jets//////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void draw1d3hist(TH1F* h1, string nc="NEvents"){

        TCanvas *canva = new TCanvas("canva","",900,780);
        gStyle->SetPadBorderMode(0);
        canva->cd();

        h1->GetYaxis()->SetTitle(nc.c_str());
        h1->GetXaxis()->SetTitle("score");
        //rapidity_data->GetXaxis()->SetTitle("M #mu#mu, GeV");
        //h1->GetYaxis()->SetTitleSize(0.05);
        //h1->GetYaxis()->SetLabelSize(0.05);
        h1->GetYaxis()->SetLabelOffset(0.005);
        h1->GetYaxis()->SetTitleOffset(1.95);

        gPad->SetTopMargin(0.07);
        gPad->SetRightMargin(0.1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.1);

        //h1->Rebin(40);
        //hist["genimuonnrapidityqg"]->Rebin(40);
        //hist["genimuonnrapidityqq"]->Rebin(40);

        //h1->GetXaxis()->SetRangeUser(1,600);

        h1->SetStats(0);
                                                  
  h1->SetMarkerStyle(37);
        h1->SetMarkerColor(kBlack);
        h1->SetMarkerSize(2.0);

      


        h1->Draw("hist p");
       
        //gPad->SetLogx(1);
        //gPad->SetLogy(1);

        TLegend *h_leg_tmp = new TLegend(0.6,0.6,0.8,0.8, "13 TeV");
        h_leg_tmp->Clear();
        h_leg_tmp->SetHeader(nc.c_str());
        h_leg_tmp->SetFillStyle(0);
        h_leg_tmp->SetTextSize(0.03);
        h_leg_tmp->SetBorderSize(0);
        h_leg_tmp->AddEntry(h1,"Total","lp");
       
        h_leg_tmp->Draw();
        string name = nc+".pdf";
        canva->SaveAs(name.c_str());
     canva->Close();

}



void roccurves_2l_channel() {


Float_t  LeadAK8_pt,LeadAK8_eta,LeadAK8_phi,LeadAK8_mass,LeadAK8_msoftdrop;

Float_t  HJ1_eta[50],HJ1_phi[50],HJ1_pt[50],HJ2_eta[50],HJ2_phi[50], HJ2_pt[50],HJ1_HJ2_dR_noFSR[50],Jet_eta[50],Jet_mass[50],Jet_phi[50],Jet_pt[50],H_pt,H_eta,H_phi;

Float_t GenJet_eta,GenJet_phi,GenJet_pt,GenJet_mass;
Float_t SubJet_eta,SubJet_phi,SubJet_pt,SubJet_mass;

Float_t GenJetAK8_eta,GenJetAK8_pt,GenJetAK8_phi,GenJetAK8_mass;


Float_t GenPart_eta[150],GenPart_mass[150],GenPart_phi[150],GenPart_pt[150];
 int  GenPart_genPartIdxMother[150], GenPart_pdgId[150], GenPart_status[150], GenPart_statusFlags[150];
UInt_t nGenPart;

Float_t LeadAK15_pt,LeadAK15_eta,LeadAK15_phi,LeadAK15_mass,LeadAK15_msoftdrop;

Float_t  nLeadAK15;

Float_t Electron_pt[50],Electron_phi[50],Electron_eta[50],Electron_mass[50], Muon_pt[50], Muon_eta[50], Muon_phi[50], Muon_mass[50]; 




Float_t JetPt_1,JetPt_2;
Float_t DeepJet_CvsL_1, DeepJet_CvsB_1,DeepJet_CvsL_2,DeepJet_CvsB_2,H_mass,V_mass,V_pt,HVdPhi;

Float_t  nLeadAK8;
UInt_t nJet;

int isZmm,isZee;
int controlSample;
int hJetInd1,hJetInd2;
bool twoResolvedJets;
UInt_t  nElectron, nMuon;



Float_t FatJet_ParticleNetMD_probQCDb,       
FatJet_ParticleNetMD_probQCDbb,      
FatJet_ParticleNetMD_probQCDc,       
FatJet_ParticleNetMD_probQCDcc,      
FatJet_ParticleNetMD_probQCDothers,  
FatJet_ParticleNetMD_probXbb,        
FatJet_ParticleNetMD_probXcc,        
FatJet_ParticleNetMD_probXqq,        
FatJet_btagCMVA,                     
FatJet_btagCSVV2,                    
FatJet_btagDDBvL,                    
FatJet_btagDDBvL_noMD ,              
FatJet_btagDDCvB,               
FatJet_btagDDCvB_noMD,               
FatJet_btagDDCvL,               
FatJet_btagDDCvL_noMD,               
FatJet_btagDeepB,                    
FatJet_btagHbb,                      
FatJet_deepTagMD_H4qvsQCD,           
FatJet_deepTagMD_HbbvsQCD,           
FatJet_deepTagMD_TvsQCD,             
FatJet_deepTagMD_WvsQCD,             
FatJet_deepTagMD_ZHbbvsQCD,          
FatJet_deepTagMD_ZHccvsQCD,          
FatJet_deepTagMD_ZbbvsQCD,           
FatJet_deepTagMD_ZvsQCD,             
FatJet_deepTagMD_bbvsLight,          
FatJet_deepTagMD_ccvsLight,          
FatJet_deepTag_H,                    
FatJet_deepTag_QCD,                
FatJet_deepTag_QCDothers,          
FatJet_deepTag_TvsQCD,                
FatJet_deepTag_WvsQCD,              
FatJet_deepTag_ZvsQCD,             
AK15Puppi_ParticleNetMD_probQCD,      
AK15Puppi_ParticleNetMD_probXbb,      
AK15Puppi_ParticleNetMD_probXcc,      
AK15Puppi_ParticleNetMD_probXqq,      
AK15Puppi_btagCSVV2,                  
AK15Puppi_btagDeepB,                  
AK15Puppi_btagJP;                    










//TString dirC = "/nfs/dust/cms/user/fismoilo/VHccAnalysisNtuples/Test_VHcc_noJetId/ZH125ToCC_ZLL_powheg";

	TChain *t = new TChain("Events");
//t->Add(dirC+"*.root");

 TChain *chain = new TChain("tree");
   for (int fileNum=1;fileNum < 12;fileNum++) {
      t->AddFile(Form("output_ZH125ToCC_ZLL_powheg_%d.root", fileNum));
//std::cout<<"File number = "<<fileNum<<std::endl; 
}
  

t->SetBranchAddress("twoResolvedJets",&twoResolvedJets);
        t->SetBranchAddress("isZmm",&isZmm);
        t->SetBranchAddress("isZee",&isZee);
        t->SetBranchAddress("JetPt_1",&JetPt_1);
        t->SetBranchAddress("JetPt_2",&JetPt_2);
        t->SetBranchAddress("DeepJet_CvsL_1",&DeepJet_CvsL_1);
        t->SetBranchAddress("DeepJet_CvsL_2",&DeepJet_CvsL_2);
        t->SetBranchAddress("DeepJet_CvsB_1",&DeepJet_CvsB_1);
        t->SetBranchAddress("DeepJet_CvsB_2",&DeepJet_CvsB_2);
        t->SetBranchAddress("hJetInd1",&hJetInd1);
        t->SetBranchAddress("hJetInd2",&hJetInd2);

        t->SetBranchAddress( "GenPart_pt",GenPart_pt);
        t->SetBranchAddress( "GenPart_phi",GenPart_phi);
        t->SetBranchAddress( "GenPart_eta",GenPart_eta);
        t->SetBranchAddress( "nGenPart",&nGenPart);
        t->SetBranchAddress( "GenPart_mass",GenPart_mass);

 //       t->SetBranchAddress("GenJet_eta",&GenJet_eta);
  //      t->SetBranchAddress("GenJet_phi",&GenJet_phi);
 //       t->SetBranchAddress("GenJet_pt",&GenJet_pt);
 //       t->SetBranchAddress("GenJet_mass",&GenJet_mass);

 //       t->SetBranchAddress("GenJetAK8_eta",&GenJet_eta);
 //       t->SetBranchAddress("GenJetAK8_phi",&GenJet_phi);
 //       t->SetBranchAddress("GenJetAK8_pt",&GenJet_pt);
 //       t->SetBranchAddress("GenJet_mass",&GenJet_mass);

 //       t->SetBranchAddress("SubJet_eta",&SubJet_eta);
 //       t->SetBranchAddress("SubJet_phi",&SubJet_phi);
 //       t->SetBranchAddress("SubJet_pt",&SubJet_pt);
 //       t->SetBranchAddress("SubJet_mass",&SubJet_mass);






        t->SetBranchAddress( "GenPart_pdgId",GenPart_pdgId);
        t->SetBranchAddress( "GenPart_genPartIdxMother",GenPart_genPartIdxMother);
        t->SetBranchAddress( "GenPart_status",GenPart_status);
        t->SetBranchAddress( "GenPart_statusFlags",GenPart_statusFlags);

        t->SetBranchAddress( "Electron_pt",Electron_pt);
        t->SetBranchAddress( "Electron_phi",Electron_phi);
        t->SetBranchAddress( "Electron_eta",Electron_eta);
        t->SetBranchAddress( "nElectron",&nElectron);
        t->SetBranchAddress( "Electron_mass",Electron_mass);

        t->SetBranchAddress( "Muon_pt",Muon_pt);
        t->SetBranchAddress( "Muon_phi",Muon_phi);
        t->SetBranchAddress( "Muon_eta",Muon_eta);
        t->SetBranchAddress( "nMuon",&nMuon);
        t->SetBranchAddress( "Muon_mass",Muon_mass);


        t->SetBranchAddress("controlSample",&controlSample);
        t->SetBranchAddress("H_mass",&H_mass);
        t->SetBranchAddress("V_pt",&V_pt);
        t->SetBranchAddress("HVdPhi",&HVdPhi);
        t->SetBranchAddress("V_mass",&V_mass);

        t->SetBranchAddress("HJ1_HJ2_dR_noFSR",HJ1_HJ2_dR_noFSR);



	t->SetBranchAddress("LeadAK8_pt",&LeadAK8_pt);
	t->SetBranchAddress("LeadAK8_eta",&LeadAK8_eta);
	t->SetBranchAddress("LeadAK8_phi",&LeadAK8_phi);
	t->SetBranchAddress("LeadAK8_mass",&LeadAK8_mass);
	t->SetBranchAddress("nLeadAK8",&nLeadAK8);
        t->SetBranchAddress("LeadAK8_msoftdrop",&LeadAK8_msoftdrop);

        t->SetBranchAddress("LeadAK15_pt",&LeadAK15_pt);
        t->SetBranchAddress("LeadAK15_eta",&LeadAK15_eta);
        t->SetBranchAddress("LeadAK15_phi",&LeadAK15_phi);
        t->SetBranchAddress("LeadAK15_mass",&LeadAK15_mass);
	t->SetBranchAddress("LeadAK15_msoftdrop",&LeadAK15_msoftdrop);
        t->SetBranchAddress("nLeadAK15",&nLeadAK15);

        t->SetBranchAddress("nJet",&nJet);
        t->SetBranchAddress("Jet_eta",Jet_eta);
        t->SetBranchAddress("Jet_phi",Jet_phi);
        t->SetBranchAddress("Jet_pt",Jet_pt);
        t->SetBranchAddress("Jet_mass",Jet_mass);


        t->SetBranchAddress("H_eta",&H_eta);
        t->SetBranchAddress("H_phi",&H_phi);
        t->SetBranchAddress("H_pt",&H_pt);
        t->SetBranchAddress("H_mass",&H_mass);

        t->SetBranchAddress("HJ1_eta",HJ1_eta);
        t->SetBranchAddress("HJ1_phi",HJ1_phi);
        t->SetBranchAddress("HJ1_pt",HJ1_pt);

        t->SetBranchAddress("HJ2_eta",HJ2_eta);
        t->SetBranchAddress("HJ2_phi",HJ2_phi);
        t->SetBranchAddress("HJ2_pt",HJ2_pt);

 t->SetBranchAddress("FatJet_ParticleNetMD_probQCDb",&FatJet_ParticleNetMD_probQCDb);     
 t->SetBranchAddress("FatJet_ParticleNetMD_probQCDbb",&  FatJet_ParticleNetMD_probQCDbb    );
 t->SetBranchAddress("FatJet_ParticleNetMD_probQCDc",&FatJet_ParticleNetMD_probQCDc  );     
 t->SetBranchAddress("FatJet_ParticleNetMD_probQCDcc",& FatJet_ParticleNetMD_probQCDcc     );
 t->SetBranchAddress("FatJet_ParticleNetMD_probQCDothers",& FatJet_ParticleNetMD_probQCDothers );
 t->SetBranchAddress("FatJet_ParticleNetMD_probXbb", &FatJet_ParticleNetMD_probXbb       );
 t->SetBranchAddress("FatJet_ParticleNetMD_probXcc", & FatJet_ParticleNetMD_probXcc      );
 t->SetBranchAddress("FatJet_ParticleNetMD_probXqq", & FatJet_ParticleNetMD_probXqq      );
 t->SetBranchAddress("FatJet_btagCMVA",        &FatJet_btagCMVA             );
 t->SetBranchAddress("FatJet_btagCSVV2",      & FatJet_btagCSVV2             );
 t->SetBranchAddress("FatJet_btagDDBvL",       &FatJet_btagDDBvL             );
 t->SetBranchAddress("FatJet_btagDDBvL_noMD" ,  &FatJet_btagDDBvL_noMD            );
 t->SetBranchAddress("FatJet_btagDDCvB"     ,   &FatJet_btagDDCvB            );
 t->SetBranchAddress("FatJet_btagDDCvB_noMD",   &FatJet_btagDDCvB_noMD            );
 t->SetBranchAddress("FatJet_btagDDCvL"     ,   &FatJet_btagDDCvL            );
 t->SetBranchAddress("FatJet_btagDDCvL_noMD",   &FatJet_btagDDCvL_noMD            );
 t->SetBranchAddress("FatJet_btagDeepB",         &FatJet_btagDeepB           );
 t->SetBranchAddress("FatJet_btagHbb",           &FatJet_btagHbb           );
 t->SetBranchAddress("FatJet_deepTagMD_H4qvsQCD", & FatJet_deepTagMD_H4qvsQCD         );
 t->SetBranchAddress("FatJet_deepTagMD_HbbvsQCD", &FatJet_deepTagMD_HbbvsQCD          );
 t->SetBranchAddress("FatJet_deepTagMD_TvsQCD",   &FatJet_deepTagMD_TvsQCD          );
 t->SetBranchAddress("FatJet_deepTagMD_WvsQCD",   &FatJet_deepTagMD_WvsQCD          );
 t->SetBranchAddress("FatJet_deepTagMD_ZHbbvsQCD", &FatJet_deepTagMD_ZHbbvsQCD         );
 t->SetBranchAddress("FatJet_deepTagMD_ZHccvsQCD", &FatJet_deepTagMD_ZHccvsQCD         );
 t->SetBranchAddress("FatJet_deepTagMD_ZbbvsQCD",  & FatJet_deepTagMD_ZbbvsQCD        );
 t->SetBranchAddress("FatJet_deepTagMD_ZvsQCD",     &FatJet_deepTagMD_ZvsQCD        );
 t->SetBranchAddress("FatJet_deepTagMD_bbvsLight", &FatJet_deepTagMD_bbvsLight         );
 t->SetBranchAddress("FatJet_deepTagMD_ccvsLight",  &FatJet_deepTagMD_ccvsLight        );
 t->SetBranchAddress("FatJet_deepTag_H",        &FatJet_deepTag_H            );
 t->SetBranchAddress("FatJet_deepTag_QCD",      &FatJet_deepTag_QCD          );
 t->SetBranchAddress("FatJet_deepTag_QCDothers", &FatJet_deepTag_QCDothers         );
 t->SetBranchAddress("FatJet_deepTag_TvsQCD",    &FatJet_deepTag_TvsQCD           ); 
 t->SetBranchAddress("FatJet_deepTag_WvsQCD",  &FatJet_deepTag_WvsQCD            );
 t->SetBranchAddress("FatJet_deepTag_ZvsQCD",  &FatJet_deepTag_ZvsQCD           );
 t->SetBranchAddress("AK15Puppi_ParticleNetMD_probQCD",  &AK15Puppi_ParticleNetMD_probQCD   ); 
 t->SetBranchAddress("AK15Puppi_ParticleNetMD_probXbb", &AK15Puppi_ParticleNetMD_probXbb     );
 t->SetBranchAddress("AK15Puppi_ParticleNetMD_probXcc",  &AK15Puppi_ParticleNetMD_probXcc    );
 t->SetBranchAddress("AK15Puppi_ParticleNetMD_probXqq", &AK15Puppi_ParticleNetMD_probXqq     );
 t->SetBranchAddress("AK15Puppi_btagCSVV2", &AK15Puppi_btagCSVV2                 );
 t->SetBranchAddress("AK15Puppi_btagDeepB",  &AK15Puppi_btagDeepB                );
     t->SetBranchAddress("AK15Puppi_btagJP", &AK15Puppi_btagJP);

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////To plot histogramms read  .root file below////////////////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	TFile *fnew = new TFile("histogrammi2.root", "RECREATE");


	TH1F *ak_eta = new TH1F("LeadAK8_eta","",100,-3.5,3.5);
	TH1F *ak_phi = new TH1F("LeadAK8_phi","",100,-3.5,3.5);
	TH1F *ak_pt = new TH1F("LeadAK8_pt","",20,0,1000);
	TH1F *ak_mass = new TH1F("LeadAK8_mass","",100,20,300);
	TH1F *ak_msoftdrop= new TH1F("LeadAK8_msoftdrop","",50,0,300);
        TH1F *ak_num = new TH1F("nLeadAK8","",100,0,50);

        TH1F *ak15_eta = new TH1F("LeadAK15_eta","",100,-3.5,3.5);
        TH1F *ak15_phi = new TH1F("AK15_phi","",100,-3.5,3.5);
        TH1F *ak15_pt = new TH1F("LeadAK15_pt","",20,0,1000);
 
        TH1F *ak4_eta = new TH1F("ak4_eta","",100,-3.5,3.5);
        TH1F *ak4_phi = new TH1F("ak4_phi","",100,-3.5,3.5);
        TH1F *ak4_pt = new TH1F("ak4_pt","",20,0,1000);
        TH1F *ak4_mass = new TH1F("ak4_mass","",20,0,1000);
    
TH1F *fatJet_ParticleNetMD_probQCDb= new TH1F("tagger_FatJet_ParticleNetMD1","",50,0,1);      
TH1F *fatJet_ParticleNetMD_probQCDbb= new TH1F("tagger_FatJet_ParticleNetMD2","",50,0,1);        
TH1F *fatJet_ParticleNetMD_probQCDc= new TH1F("tagger_FatJet_ParticleNetMD3","",50,0,1);         
TH1F *fatJet_ParticleNetMD_probQCDcc= new TH1F("tagger_FatJet_ParticleNetMD4","",50,0,1);        
TH1F *fatJet_ParticleNetMD_probQCDothers= new TH1F("tagger_FatJet_ParticleNetMD5","",50,0,1);    
TH1F *fatJet_ParticleNetMD_probXbb= new TH1F("tagger6","",50,0,1);          
TH1F *fatJet_ParticleNetMD_probXcc= new TH1F("tagger7","",50,0,1);          
TH1F *fatJet_ParticleNetMD_probXqq= new TH1F("tagger8","",50,0,1);          
TH1F *fatJet_btagCMVA= new TH1F("tagger9","",50,0,1);                       
TH1F *fatJet_btagCSVV2= new TH1F("tagger10","",50,0,1);                      
TH1F *fatJet_btagDDBvL= new TH1F("tagger11","",50,0,1);                      
TH1F *fatJet_btagDDBvL_noMD = new TH1F("tagger12","",50,0,1);                
TH1F *fatJet_btagDDCvB= new TH1F("tagger13","",50,0,1);                      
TH1F *fatJet_btagDDCvB_noMD= new TH1F("tagger14","",50,0,1);                 
TH1F *fatJet_btagDDCvL = new TH1F("tagger15","",50,0,1);                     
TH1F *fatJet_btagDDCvL_noMD= new TH1F("tagger16","",50,0,1);                 
TH1F *fatJet_btagDeepB= new TH1F("tagger17","",50,0,1);                      
TH1F *fatJet_btagHbb= new TH1F("tagger18","",50,0,1);                        
TH1F *fatJet_deepTagMD_H4qvsQCD= new TH1F("tagger19","",50,0,1);             
TH1F *fatJet_deepTagMD_HbbvsQCD= new TH1F("tagger20","",50,0,1);             
TH1F *fatJet_deepTagMD_TvsQCD= new TH1F("tagger21","",50,0,1);               
TH1F *fatJet_deepTagMD_WvsQCD= new TH1F("tagger22","",50,0,1);               
TH1F *fatJet_deepTagMD_ZHbbvsQCD= new TH1F("tagger23","",50,0,1);            
TH1F *fatJet_deepTagMD_ZHccvsQCD= new TH1F("tagger24","",50,0,1);            
TH1F *fatJet_deepTagMD_ZbbvsQCD= new TH1F("tagger25","",50,0,1);             
TH1F *fatJet_deepTagMD_ZvsQCD= new TH1F("tagger26","",50,0,1);               
TH1F *fatJet_deepTagMD_bbvsLight= new TH1F("tagger27","",50,0,1);            
TH1F *fatJet_deepTagMD_ccvsLight= new TH1F("tagger28","",50,0,1);            
TH1F *fatJet_deepTag_H= new TH1F("tagger29","",50,0,1);                      
TH1F *fatJet_deepTag_QCD= new TH1F("tagger30","",50,0,1);                  
TH1F *fatJet_deepTag_QCDothers= new TH1F("tagger31","",50,0,1);            
TH1F *fatJet_deepTag_TvsQCD= new TH1F("tagger32","",50,0,1);                  
TH1F *fatJet_deepTag_WvsQCD= new TH1F("tagger33","",50,0,1);                
TH1F *fatJet_deepTag_ZvsQCD= new TH1F("tagger34","",50,0,1);               
TH1F *aK15Puppi_ParticleNetMD_probQCD= new TH1F("tagger35","",50,0,1);        
TH1F *aK15Puppi_ParticleNetMD_probXbb= new TH1F("tagger36","",50,0,1);        
TH1F *aK15Puppi_ParticleNetMD_probXcc= new TH1F("tagger37","",50,0,1);        
TH1F *aK15Puppi_ParticleNetMD_probXqq= new TH1F("tagger38","",50,0,1);        
TH1F *aK15Puppi_btagCSVV2= new TH1F("tagger39","",50,0,1);                    
TH1F *aK15Puppi_btagDeepB= new TH1F("tagger40","",50,0,1);                    
TH1F *aK15Puppi_btagJP= new TH1F("tagger41","",50,0,1);       
    
  

TH1F *fatJet_ParticleNetMD_probQCDb_01= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_01","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_015= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_015","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_02= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_02","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_025= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_025","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_03= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_03","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_035= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_035","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_04= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_04","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_045= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_045","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_05= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_05","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_055= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_055","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_06= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_06","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_065= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_065","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_07= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_07","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_075= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_075","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_08= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_08","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_085= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_085","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_09= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_09","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_095= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_095","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDb_1= new TH1F("tagger_FatJet_ParticleNetMD_probQCDb_1","",50,0,1);


TH1F *fatJet_ParticleNetMD_probQCDc_01= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_01","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_015= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_015","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_02= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_02","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_025= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_025","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_03= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_03","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_035= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_035","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_04= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_04","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_045= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_045","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_05= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_05","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_055= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_055","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_06= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_06","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_065= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_065","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_07= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_07","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_075= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_075","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_08= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_08","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_085= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_085","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_09= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_09","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_095= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_095","",50,0,1);
TH1F *fatJet_ParticleNetMD_probQCDc_1= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_1","",50,0,1);    



TH1F *fatJet_deepTagMD_bbvsLight_01= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_01","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_015= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_015","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_02= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_02","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_025= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_025","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_03= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_03","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_035= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_035","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_04= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_04","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_045= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_045","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_05= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_05","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_055= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_055","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_06= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_06","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_065= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_065","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_07= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_07","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_075= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_075","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_08= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_08","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_085= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_085","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_09= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_09","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_095= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_095","",50,0,1);
TH1F *fatJet_deepTagMD_bbvsLight_1= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_1","",50,0,1);


TH1F *fatJet_deepTagMD_ccvsLight_01= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_01","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_015= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_015","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_02= new TH1F("tagger_FatJet_deepTagMD_cvsLight_02","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_025= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_025","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_03= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_03","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_035= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_035","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_04= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_04","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_045= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_045","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_05= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_05","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_055= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_055","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_06= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_06","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_065= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_065","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_07= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_07","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_075= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_075","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_08= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_08","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_085= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_085","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_09= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_09","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_095= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_095","",50,0,1);
TH1F *fatJet_deepTagMD_ccvsLight_1= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_1","",50,0,1);

    
    
TH1F *fat_btagDDBvL_01= new TH1F("tagger_btagDDBvL_01","",50,0,1);
TH1F *fat_btagDDBvL_015= new TH1F("tagger_btagDDBvL_015","",50,0,1);
TH1F *fat_btagDDBvL_02= new TH1F("tagger_btagDDBvL_02","",50,0,1);
TH1F *fat_btagDDBvL_025= new TH1F("tagger_btagDDBvL_025","",50,0,1);
TH1F *fat_btagDDBvL_03= new TH1F("tagger_btagDDBvL_03","",50,0,1);
TH1F *fat_btagDDBvL_035= new TH1F("tagger_btagDDBvL_035","",50,0,1);
TH1F *fat_btagDDBvL_04= new TH1F("tagger_btagDDBvL_04","",50,0,1);
TH1F *fat_btagDDBvL_045= new TH1F("tagger_btagDDBvL_045","",50,0,1);
TH1F *fat_btagDDBvL_05= new TH1F("tagger_btagDDBvL_05","",50,0,1);
TH1F *fat_btagDDBvL_055= new TH1F("tagger_btagDDBvL_055","",50,0,1);
TH1F *fat_btagDDBvL_06= new TH1F("tagger_btagDDBvL_06","",50,0,1);
TH1F *fat_btagDDBvL_065= new TH1F("tagger_btagDDBvL_065","",50,0,1);
TH1F *fat_btagDDBvL_07= new TH1F("tagger_btagDDBvL_07","",50,0,1);
TH1F *fat_btagDDBvL_075= new TH1F("tagger_btagDDBvL_075","",50,0,1);
TH1F *fat_btagDDBvL_08= new TH1F("tagger_btagDDBvL_08","",50,0,1);
TH1F *fat_btagDDBvL_085= new TH1F("tagger_btagDDBvL_085","",50,0,1);
TH1F *fat_btagDDBvL_09= new TH1F("tagger_btagDDBvL_09","",50,0,1);
TH1F *fat_btagDDBvL_095= new TH1F("tagger_btagDDBvL_095","",50,0,1);
TH1F *fat_btagDDBvL_1= new TH1F("tagger_btagDDBvL_1","",50,0,1);

        
TH1F *fat_btagDDCvL_01= new TH1F("tagger_btagDDCvL_01","",50,0,1);
TH1F *fat_btagDDCvL_015= new TH1F("tagger_btagDDCvL_015","",50,0,1);
TH1F *fat_btagDDCvL_02= new TH1F("tagger_btagDDCvL_02","",50,0,1);
TH1F *fat_btagDDCvL_025= new TH1F("tagger_btagDDCvL_025","",50,0,1);
TH1F *fat_btagDDCvL_03= new TH1F("tagger_btagDDCvL_03","",50,0,1);
TH1F *fat_btagDDCvL_035= new TH1F("tagger_btagDDCvL_035","",50,0,1);
TH1F *fat_btagDDCvL_04= new TH1F("tagger_btagDDCvL_04","",50,0,1);
TH1F *fat_btagDDCvL_045= new TH1F("tagger_btagDDCvL_045","",50,0,1);
TH1F *fat_btagDDCvL_05= new TH1F("tagger_btagDDCvL_05","",50,0,1);
TH1F *fat_btagDDCvL_055= new TH1F("tagger_btagDDCvL_055","",50,0,1);
TH1F *fat_btagDDCvL_06= new TH1F("tagger_btagDDCvL_06","",50,0,1);
TH1F *fat_btagDDCvL_065= new TH1F("tagger_btagDDCvL_065","",50,0,1);
TH1F *fat_btagDDCvL_07= new TH1F("tagger_btagDDCvL_07","",50,0,1);
TH1F *fat_btagDDCvL_075= new TH1F("tagger_btagDDCvL_075","",50,0,1);
TH1F *fat_btagDDCvL_08= new TH1F("tagger_btagDDCvL_08","",50,0,1);
TH1F *fat_btagDDCvL_085= new TH1F("tagger_btagDDCvL_085","",50,0,1);
TH1F *fat_btagDDCvL_09= new TH1F("tagger_btagDDCvL_09","",50,0,1);
TH1F *fat_btagDDCvL_095= new TH1F("tagger_btagDDCvL_095","",50,0,1);
TH1F *fat_btagDDCvL_1= new TH1F("tagger_btagDDCvL_1","",50,0,1);
 
        
TH1F *fat_btagDDCvB_01= new TH1F("tagger_btagDDCvB_01","",50,0,1);
TH1F *fat_btagDDCvB_015= new TH1F("tagger_btagDDCvB_015","",50,0,1);
TH1F *fat_btagDDCvB_02= new TH1F("tagger_btagDDCvB_02","",50,0,1);
TH1F *fat_btagDDCvB_025= new TH1F("tagger_btagDDCvB_025","",50,0,1);
TH1F *fat_btagDDCvB_03= new TH1F("tagger_btagDDCvB_03","",50,0,1);
TH1F *fat_btagDDCvB_035= new TH1F("tagger_btagDDCvB_035","",50,0,1);
TH1F *fat_btagDDCvB_04= new TH1F("tagger_btagDDCvB_04","",50,0,1);
TH1F *fat_btagDDCvB_045= new TH1F("tagger_btagDDCvB_045","",50,0,1);
TH1F *fat_btagDDCvB_05= new TH1F("tagger_btagDDCvB_05","",50,0,1);
TH1F *fat_btagDDCvB_055= new TH1F("tagger_btagDDCvB_055","",50,0,1);
TH1F *fat_btagDDCvB_06= new TH1F("tagger_btagDDCvB_06","",50,0,1);
TH1F *fat_btagDDCvB_065= new TH1F("tagger_btagDDCvB_065","",50,0,1);
TH1F *fat_btagDDCvB_07= new TH1F("tagger_btagDDCvB_07","",50,0,1);
TH1F *fat_btagDDCvB_075= new TH1F("tagger_btagDDCvB_075","",50,0,1);
TH1F *fat_btagDDCvB_08= new TH1F("tagger_btagDDCvB_08","",50,0,1);
TH1F *fat_btagDDCvB_085= new TH1F("tagger_btagDDCvB_085","",50,0,1);
TH1F *fat_btagDDCvB_09= new TH1F("tagger_btagDDCvB_09","",50,0,1);
TH1F *fat_btagDDCvB_095= new TH1F("tagger_btagDDCvB_095","",50,0,1);
TH1F *fat_btagDDCvB_1= new TH1F("tagger_btagDDCvB_1","",50,0,1);    
    
 /*   
    FatJet_btagDDBvL,                    
FatJet_btagDDBvL_noMD ,              
FatJet_btagDDCvB,               
FatJet_btagDDCvB_noMD,               
FatJet_btagDDCvL,               
FatJet_btagDDCvL_noMD,               
FatJet_btagDeepB,                    
FatJet_btagHbb,     
*/
    
     
    
     
    
    
    
float ParticleNetMD_probQCDb_01 ;     
float ParticleNetMD_probQCDb_015;              
float ParticleNetMD_probQCDb_02;       
float ParticleNetMD_probQCDb_025;         
float ParticleNetMD_probQCDb_03;        
float ParticleNetMD_probQCDb_035;
float ParticleNetMD_probQCDb_04;
float ParticleNetMD_probQCDb_045;      
float ParticleNetMD_probQCDb_05;              
float ParticleNetMD_probQCDb_055;       
float ParticleNetMD_probQCDb_06;         
float ParticleNetMD_probQCDb_065;        
float ParticleNetMD_probQCDb_07;
float ParticleNetMD_probQCDb_075;        
float ParticleNetMD_probQCDb_08;       
float ParticleNetMD_probQCDb_085;         
float ParticleNetMD_probQCDb_09;        
float ParticleNetMD_probQCDb_095;
float ParticleNetMD_probQCDb_1;
    
float ParticleNetMD_probQCDc_01;      
float ParticleNetMD_probQCDc_015;              
float ParticleNetMD_probQCDc_02;       
float ParticleNetMD_probQCDc_025;         
float ParticleNetMD_probQCDc_03;        
float ParticleNetMD_probQCDc_035;
float ParticleNetMD_probQCDc_04;
float ParticleNetMD_probQCDc_045;      
float ParticleNetMD_probQCDc_05;              
float ParticleNetMD_probQCDc_055;       
float ParticleNetMD_probQCDc_06;         
float ParticleNetMD_probQCDc_065;        
float ParticleNetMD_probQCDc_07;
float ParticleNetMD_probQCDc_075;        
float ParticleNetMD_probQCDc_08;       
float ParticleNetMD_probQCDc_085;         
float ParticleNetMD_probQCDc_09 ;       
float ParticleNetMD_probQCDc_095;
float ParticleNetMD_probQCDc_1;    
    
float btagDDBvL_01 ;     
float btagDDBvL_015;              
float btagDDBvL_02;       
float btagDDBvL_025;         
float btagDDBvL_03;        
float btagDDBvL_035;
float btagDDBvL_04;
float btagDDBvL_045;      
float btagDDBvL_05;              
float btagDDBvL_055;       
float btagDDBvL_06;         
float btagDDBvL_065;        
float btagDDBvL_07;
float btagDDBvL_075;        
float btagDDBvL_08;       
float btagDDBvL_085;         
float btagDDBvL_09;        
float btagDDBvL_095;
float btagDDBvL_1;
    
float btagDDBvL_noMD_01 ;     
float btagDDBvL_noMD_015;              
float btagDDBvL_noMD_02;       
float btagDDBvL_noMD_025;         
float btagDDBvL_noMD_03;        
float btagDDBvL_noMD_035;
float btagDDBvL_noMD_04;
float btagDDBvL_noMD_045;      
float btagDDBvL_noMD_05;              
float btagDDBvL_noMD_055;       
float btagDDBvL_noMD_06;         
float btagDDBvL_noMD_065;        
float btagDDBvL_noMD_07;
float btagDDBvL_noMD_075;        
float btagDDBvL_noMD_08;       
float btagDDBvL_noMD_085;         
float btagDDBvL_noMD_09;        
float btagDDBvL_noMD_095;
float btagDDBvL_noMD_1;

    float btagDDCvL_01 ;     
float btagDDCvL_015;              
float btagDDCvL_02;       
float btagDDCvL_025;         
float btagDDCvL_03;        
float btagDDCvL_035;
float btagDDCvL_04;
float btagDDCvL_045;      
float btagDDCvL_05;              
float btagDDCvL_055;       
float btagDDCvL_06;         
float btagDDCvL_065;        
float btagDDCvL_07;
float btagDDCvL_075;        
float btagDDCvL_08;       
float btagDDCvL_085;         
float btagDDCvL_09;        
float btagDDCvL_095;
float btagDDCvL_1;
    
float btagDDCvL_noMD_01 ;     
float btagDDCvL_noMD_015;              
float btagDDCvL_noMD_02;       
float btagDDCvL_noMD_025;         
float btagDDCvL_noMD_03;        
float btagDDCvL_noMD_035;
float btagDDCvL_noMD_04;
float btagDDCvL_noMD_045;      
float btagDDCvL_noMD_05;              
float btagDDCvL_noMD_055;       
float btagDDCvL_noMD_06;         
float btagDDCvL_noMD_065;        
float btagDDCvL_noMD_07;
float btagDDCvL_noMD_075;        
float btagDDCvL_noMD_08;       
float btagDDCvL_noMD_085;         
float btagDDCvL_noMD_09;        
float btagDDCvL_noMD_095;
float btagDDCvL_noMD_1;

 float btagDDCvB_01 ;     
float btagDDCvB_015;              
float btagDDCvB_02;       
float btagDDCvB_025;         
float btagDDCvB_03;        
float btagDDCvB_035;
float btagDDCvB_04;
float btagDDCvB_045;      
float btagDDCvB_05;              
float btagDDCvB_055;       
float btagDDCvB_06;         
float btagDDCvB_065;        
float btagDDCvB_07;
float btagDDCvB_075;        
float btagDDCvB_08;       
float btagDDCvB_085;         
float btagDDCvB_09;        
float btagDDCvB_095;
float btagDDCvB_1;
    
float btagDDCvB_noMD_01 ;     
float btagDDCvB_noMD_015;              
float btagDDCvB_noMD_02;       
float btagDDCvB_noMD_025;         
float btagDDCvB_noMD_03;        
float btagDDCvB_noMD_035;
float btagDDCvB_noMD_04;
float btagDDCvB_noMD_045;      
float btagDDCvB_noMD_05;              
float btagDDCvB_noMD_055;       
float btagDDCvB_noMD_06;         
float btagDDCvB_noMD_065;        
float btagDDCvB_noMD_07;
float btagDDCvB_noMD_075;        
float btagDDCvB_noMD_08;       
float btagDDCvB_noMD_085;         
float btagDDCvB_noMD_09;        
float btagDDCvB_noMD_095;
float btagDDCvB_noMD_1;


        std::vector<Double_t> bbvsLight; 
        std::vector<Double_t> ccvsLight;

        std::vector<Double_t> bjet;
        std::vector<Double_t> cjet;
    
   std::vector<Double_t>  btagDDBvL;
   std::vector<Double_t>  btagDDBvL_noMD;
    
      std::vector<Double_t>  btagDDCvB;
   std::vector<Double_t>  btagDDCvB_noMD; 
    
    std::vector<Double_t>  btagDDCvL;
   std::vector<Double_t>  btagDDCvL_noMD; 
    
 int NumberOfEvents =66540;
//	int NumberOfEvents =69594;
	std::cout<< t->GetEntries()<<std::endl;

	for (int i=0; i<NumberOfEvents;i++) {

		if ((i % 10000) == 0)
		{
			std::cout<<"Event = "<<i<<std::endl;
		}
		t->GetEntry(i);
	




        TLorentzVector HJ1;
        TLorentzVector HJ2;		
        
        TLorentzVector Hcc;
        TLorentzVector AK8Higgs,AK15Higgs;

        TLorentzVector Zboson;
        TLorentzVector Zee,Zmumu,Jet,Electron1, Electron2, Muon1, Muon2; 

        TLorentzVector GenHiggs;
        TLorentzVector LeadAK15;
        TLorentzVector LeadAK8;

        TLorentzVector GenJet,GenJet1,GenJet2;



  Hcc.SetPtEtaPhiM(H_pt,H_eta,H_phi,H_mass);


fatJet_ParticleNetMD_probQCDb->Fill(FatJet_ParticleNetMD_probQCDb );              
fatJet_ParticleNetMD_probQCDbb->Fill( FatJet_ParticleNetMD_probQCDbb  );       
fatJet_ParticleNetMD_probQCDc->Fill(FatJet_ParticleNetMD_probQCDc    );       
fatJet_ParticleNetMD_probQCDcc->Fill(FatJet_ParticleNetMD_probQCDcc      );   
fatJet_ParticleNetMD_probQCDothers->Fill(FatJet_ParticleNetMD_probQCDothers  );    
fatJet_ParticleNetMD_probXbb->Fill(FatJet_ParticleNetMD_probXbb);            
fatJet_ParticleNetMD_probXcc->Fill(FatJet_ParticleNetMD_probXcc);            
fatJet_ParticleNetMD_probXqq->Fill(FatJet_ParticleNetMD_probXqq);            
fatJet_btagCMVA->Fill(FatJet_btagCMVA );                        
fatJet_btagCSVV2->Fill(FatJet_btagCSVV2);                        
fatJet_btagDDBvL->Fill(FatJet_btagDDBvL);                        
fatJet_btagDDBvL_noMD->Fill(FatJet_btagDDBvL_noMD);                  
fatJet_btagDDCvB->Fill(FatJet_btagDDCvB);                   
fatJet_btagDDCvB_noMD->Fill(FatJet_btagDDCvB_noMD);                   
fatJet_btagDDCvL->Fill(FatJet_btagDDCvL);                   
fatJet_btagDDCvL_noMD->Fill(FatJet_btagDDCvL_noMD);                   
fatJet_btagDeepB->Fill(FatJet_btagDeepB);                        
fatJet_btagHbb->Fill(FatJet_btagHbb);                          
fatJet_deepTagMD_H4qvsQCD->Fill(FatJet_deepTagMD_H4qvsQCD);               
fatJet_deepTagMD_HbbvsQCD->Fill(FatJet_deepTagMD_HbbvsQCD);               
fatJet_deepTagMD_TvsQCD->Fill(FatJet_deepTagMD_TvsQCD);                 
fatJet_deepTagMD_WvsQCD->Fill(FatJet_deepTagMD_WvsQCD);                 
fatJet_deepTagMD_ZHbbvsQCD->Fill(FatJet_deepTagMD_ZHbbvsQCD);              
fatJet_deepTagMD_ZHccvsQCD->Fill(FatJet_deepTagMD_ZHccvsQCD);              
fatJet_deepTagMD_ZbbvsQCD->Fill(FatJet_deepTagMD_ZbbvsQCD);               
fatJet_deepTagMD_ZvsQCD->Fill(FatJet_deepTagMD_ZvsQCD );                
fatJet_deepTagMD_bbvsLight->Fill(FatJet_deepTagMD_bbvsLight);              
fatJet_deepTagMD_ccvsLight->Fill(FatJet_deepTagMD_ccvsLight);              
fatJet_deepTag_H->Fill(FatJet_deepTag_H );                       
fatJet_deepTag_QCD->Fill(FatJet_deepTag_QCD );                   
fatJet_deepTag_QCDothers->Fill(FatJet_deepTag_QCDothers);              
fatJet_deepTag_TvsQCD->Fill(FatJet_deepTag_TvsQCD );                   
fatJet_deepTag_WvsQCD->Fill(FatJet_deepTag_WvsQCD);                  
fatJet_deepTag_ZvsQCD->Fill(FatJet_deepTag_ZvsQCD );                
aK15Puppi_ParticleNetMD_probQCD->Fill(AK15Puppi_ParticleNetMD_probQCD);          
aK15Puppi_ParticleNetMD_probXbb->Fill(AK15Puppi_ParticleNetMD_probXbb );         
aK15Puppi_ParticleNetMD_probXcc->Fill(AK15Puppi_ParticleNetMD_probXcc);          
aK15Puppi_ParticleNetMD_probXqq->Fill(AK15Puppi_ParticleNetMD_probXqq );         
aK15Puppi_btagCSVV2->Fill(AK15Puppi_btagCSVV2 );                     
aK15Puppi_btagDeepB->Fill(AK15Puppi_btagDeepB );                     
aK15Puppi_btagJP->Fill(AK15Puppi_btagJP );        
        
        
        
        
        ////////////////// FatJet_ParticleNetMD_probQCDb////////////////////
  /*      
         if (FatJet_ParticleNetMD_probQCDb>0.1){
fatJet_ParticleNetMD_probQCDb_01->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.15){
fatJet_ParticleNetMD_probQCDb_015->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.2){
fatJet_ParticleNetMD_probQCDb_02->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.25){
fatJet_ParticleNetMD_probQCDb_025->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.3){
fatJet_ParticleNetMD_probQCDb_03->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.35){
fatJet_ParticleNetMD_probQCDb_035->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.4){
fatJet_ParticleNetMD_probQCDb_04->Fill(FatJet_ParticleNetMD_probQCDb);}
        
         if (FatJet_ParticleNetMD_probQCDb>0.45){
fatJet_ParticleNetMD_probQCDb_045->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.5){
fatJet_ParticleNetMD_probQCDb_05->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.55){
fatJet_ParticleNetMD_probQCDb_055->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.6){
fatJet_ParticleNetMD_probQCDb_06->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.65){
fatJet_ParticleNetMD_probQCDb_065->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.7){
fatJet_ParticleNetMD_probQCDb_07->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.75){
fatJet_ParticleNetMD_probQCDb_075->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.8){
fatJet_ParticleNetMD_probQCDb_08->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.85){
fatJet_ParticleNetMD_probQCDb_085->Fill(FatJet_ParticleNetMD_probQCDb);}
         if (FatJet_ParticleNetMD_probQCDb>0.9){
fatJet_ParticleNetMD_probQCDb_09->Fill(FatJet_ParticleNetMD_probQCDb);}
   if (FatJet_ParticleNetMD_probQCDb>0.95){
fatJet_ParticleNetMD_probQCDb_095->Fill(FatJet_ParticleNetMD_probQCDb);}
  
         if (FatJet_ParticleNetMD_probQCDb==1){
fatJet_ParticleNetMD_probQCDb_1->Fill(FatJet_ParticleNetMD_probQCDb);}
        
       ////////////////// FatJet_ParticleNetMD_probQCDc////////////////////
        
        
      if (FatJet_ParticleNetMD_probQCDc>0.1){
fatJet_ParticleNetMD_probQCDc_01->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.15){
fatJet_ParticleNetMD_probQCDc_015->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.2){
fatJet_ParticleNetMD_probQCDc_02->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.25){
fatJet_ParticleNetMD_probQCDc_025->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.3){
fatJet_ParticleNetMD_probQCDc_03->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.35){
fatJet_ParticleNetMD_probQCDc_035->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.4){
fatJet_ParticleNetMD_probQCDc_04->Fill(FatJet_ParticleNetMD_probQCDc);}
        
         if (FatJet_ParticleNetMD_probQCDc>0.45){
fatJet_ParticleNetMD_probQCDc_045->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.5){
fatJet_ParticleNetMD_probQCDc_05->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.55){
fatJet_ParticleNetMD_probQCDc_055->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.6){
fatJet_ParticleNetMD_probQCDc_06->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.65){
fatJet_ParticleNetMD_probQCDc_065->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.7){
fatJet_ParticleNetMD_probQCDc_07->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.75){
fatJet_ParticleNetMD_probQCDc_075->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.8){
fatJet_ParticleNetMD_probQCDc_08->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.85){
fatJet_ParticleNetMD_probQCDc_085->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc>0.9){
fatJet_ParticleNetMD_probQCDc_09->Fill(FatJet_ParticleNetMD_probQCDc);}
   if (FatJet_ParticleNetMD_probQCDc>0.95){
fatJet_ParticleNetMD_probQCDc_095->Fill(FatJet_ParticleNetMD_probQCDc);}
         if (FatJet_ParticleNetMD_probQCDc==1){
fatJet_ParticleNetMD_probQCDc_1->Fill(FatJet_ParticleNetMD_probQCDc);}
 */
        
        
        
        
        
          if (FatJet_btagDDBvL>0.1){
fat_btagDDBvL_01->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.15){
fat_btagDDBvL_015->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.2){
fat_btagDDBvL_02->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.25){
fat_btagDDBvL_025->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.3){
fat_btagDDBvL_03->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.35){
fat_btagDDBvL_035->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.4){
fat_btagDDBvL_04->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.45){
fat_btagDDBvL_045->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.5){
fat_btagDDBvL_05->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.55){
fat_btagDDBvL_055->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.6){
fat_btagDDBvL_06->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.65){
fat_btagDDBvL_065->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.7){
fat_btagDDBvL_07->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.75){
fat_btagDDBvL_075->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.8){
fat_btagDDBvL_08->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.85){
fat_btagDDBvL_085->Fill(FatJet_btagDDBvL);}
        
         if (FatJet_btagDDBvL>0.9){
fat_btagDDBvL_09->Fill(FatJet_btagDDBvL);}
        
   if (FatJet_btagDDBvL>0.95){
fat_btagDDBvL_095->Fill(FatJet_btagDDBvL);}
  
         if (FatJet_btagDDBvL==1){
fat_btagDDBvL_1->Fill(FatJet_btagDDBvL);}
        
        
        
        
        
        
        
         if (FatJet_btagDDCvB>0.1){fat_btagDDCvB_01->Fill(FatJet_btagDDCvB);}      
         if (FatJet_btagDDCvB>0.15){fat_btagDDCvB_015->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.2){fat_btagDDCvB_02->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.25){fat_btagDDCvB_025->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.3){fat_btagDDCvB_03->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.35){fat_btagDDCvB_035->Fill(FatJet_btagDDCvB);}  
         if (FatJet_btagDDCvB>0.4){fat_btagDDCvB_04->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.45){fat_btagDDCvB_045->Fill(FatJet_btagDDCvB);}       
         if (FatJet_btagDDCvB>0.5){fat_btagDDCvB_05->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.55){fat_btagDDCvB_055->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.6){fat_btagDDCvB_06->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.65){fat_btagDDCvB_065->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.7){fat_btagDDCvB_07->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.75){fat_btagDDCvB_075->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.8){fat_btagDDCvB_08->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.85){fat_btagDDCvB_085->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.9){fat_btagDDCvB_09->Fill(FatJet_btagDDCvB);}        
         if (FatJet_btagDDCvB>0.95){fat_btagDDCvB_095->Fill(FatJet_btagDDCvB);}  
         if (FatJet_btagDDCvB==1){fat_btagDDCvB_1->Fill(FatJet_btagDDCvB);}
        
       ////////////////// FatJet_ParticleNetMD_probQCDc////////////////////
        
        
     
 
        
        
        
        
 //////////////bbvsLight////////////7
 
 if (FatJet_deepTagMD_bbvsLight>0.1){
fatJet_deepTagMD_bbvsLight_01->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.15){
fatJet_deepTagMD_bbvsLight_015->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.2){
fatJet_deepTagMD_bbvsLight_02->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.25){
fatJet_deepTagMD_bbvsLight_025->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.3){
fatJet_deepTagMD_bbvsLight_03->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.35){
fatJet_deepTagMD_bbvsLight_035->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.4){
fatJet_deepTagMD_bbvsLight_04->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight>0.45){
fatJet_deepTagMD_bbvsLight_045->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.5){
fatJet_deepTagMD_bbvsLight_05->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.55){
fatJet_deepTagMD_bbvsLight_055->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.6){
fatJet_deepTagMD_bbvsLight_06->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.65){
fatJet_deepTagMD_bbvsLight_065->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.7){
fatJet_deepTagMD_bbvsLight_07->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.75){
fatJet_deepTagMD_bbvsLight_075->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.8){
fatJet_deepTagMD_bbvsLight_08->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.85){
fatJet_deepTagMD_bbvsLight_085->Fill(FatJet_deepTagMD_bbvsLight);}
         if (FatJet_deepTagMD_bbvsLight>0.9){
fatJet_deepTagMD_bbvsLight_09->Fill(FatJet_deepTagMD_bbvsLight);}
   if (FatJet_deepTagMD_bbvsLight>0.95){
fatJet_deepTagMD_bbvsLight_095->Fill(FatJet_deepTagMD_bbvsLight);}

         if (FatJet_deepTagMD_bbvsLight==1){
fatJet_deepTagMD_bbvsLight_1->Fill(FatJet_deepTagMD_bbvsLight);}

       
        
 ///////ccvsLight/////////
 
if (FatJet_deepTagMD_ccvsLight>0.1){
fatJet_deepTagMD_ccvsLight_01->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.15){
fatJet_deepTagMD_ccvsLight_015->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.2){
fatJet_deepTagMD_ccvsLight_02->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.25){
fatJet_deepTagMD_ccvsLight_025->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.3){
fatJet_deepTagMD_ccvsLight_03->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.35){
fatJet_deepTagMD_ccvsLight_035->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.4){
fatJet_deepTagMD_ccvsLight_04->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight>0.45){
fatJet_deepTagMD_ccvsLight_045->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.5){
fatJet_deepTagMD_ccvsLight_05->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.55){
fatJet_deepTagMD_ccvsLight_055->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.6){
fatJet_deepTagMD_ccvsLight_06->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.65){
fatJet_deepTagMD_ccvsLight_065->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.7){
fatJet_deepTagMD_ccvsLight_07->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.75){
fatJet_deepTagMD_ccvsLight_075->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.8){
fatJet_deepTagMD_ccvsLight_08->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.85){
fatJet_deepTagMD_ccvsLight_085->Fill(FatJet_deepTagMD_ccvsLight);}
         if (FatJet_deepTagMD_ccvsLight>0.9){
fatJet_deepTagMD_ccvsLight_09->Fill(FatJet_deepTagMD_ccvsLight);}
   if (FatJet_deepTagMD_ccvsLight>0.95){
fatJet_deepTagMD_ccvsLight_095->Fill(FatJet_deepTagMD_ccvsLight);}

         if (FatJet_deepTagMD_ccvsLight==1){
fatJet_deepTagMD_ccvsLight_1->Fill(FatJet_deepTagMD_ccvsLight);}       
        
}




 //ParticleNetMD_probQCDb=fatJet_ParticleNetMD_probQCDb->Integral();
//cout<<"ParticleNetMD_probQCDb ="<<ParticleNetMD_probQCDb<<endl;

//cut_0_1=fatJet_ParticleNetMD_probQCDb_01->Integral();
//eff_0_1_fatJet_ParticleNetMD_probQCDb=cut_0_1/ParticleNetMD_probQCDb;
//cout<<" cut_0_1 ="<< cut_0_1 <<endl;
//cout<<"efficiency ="<< eff_0_1_fatJet_ParticleNetMD_probQCDb<<endl;


    btagDDBvL.push_back((fat_btagDDBvL_01->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_015->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_02->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_025->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_03->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_035->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_04->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_045->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_05->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_055->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_06->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_065->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_07->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_075->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_08->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_085->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_09->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_095->Integral())/(fatJet_btagDDBvL->Integral()));
    btagDDBvL.push_back((fat_btagDDBvL_1->Integral())/(fatJet_btagDDBvL->Integral()));
    
    btagDDCvL.push_back((fat_btagDDCvL_01->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_015->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_02->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_025->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_03->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_035->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_04->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_045->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_05->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_055->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_06->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_065->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_07->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_075->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_08->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_085->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_09->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_095->Integral())/(fatJet_btagDDCvL->Integral()));
    btagDDCvL.push_back((fat_btagDDCvL_1->Integral())/(fatJet_btagDDCvL->Integral()));
    
    btagDDCvB.push_back((fat_btagDDCvB_01->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_015->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_02->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_025->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_03->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_035->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_04->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_045->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_05->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_055->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_06->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_065->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_07->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_075->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_08->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_085->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_09->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_095->Integral())/(fatJet_btagDDCvB->Integral()));
    btagDDCvB.push_back((fat_btagDDCvB_1->Integral())/(fatJet_btagDDCvB->Integral()));
    
    
    


    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_01->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_015->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_02->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_025->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_03->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_035->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_04->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_045->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_05->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_055->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_06->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_065->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_07->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_075->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_08->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_085->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_09->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_095->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));
    bbvsLight.push_back((fatJet_deepTagMD_bbvsLight_1->Integral())/(fatJet_deepTagMD_bbvsLight->Integral()));

    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_01->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_015->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_02->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_025->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_03->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_035->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_04->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_045->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_05->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_055->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_06->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_065->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
   ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_07->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_075->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_08->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_085->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_09->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
   ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_095->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));
    ccvsLight.push_back((fatJet_deepTagMD_ccvsLight_1->Integral())/(fatJet_deepTagMD_ccvsLight->Integral()));



    TCanvas *canvas2 = new TCanvas("canvas2", "", 1900, 1100);
                canvas2->cd();
    
    
TMultiGraph  *mg  = new TMultiGraph();
	
	
	TGraph g1(TGraph(btagDDBvL.size(), &btagDDBvL[0], &btagDDCvB[0]));
	g1.SetLineColor(kRed);
	mg->Add(&g1);
    mg->GetXaxis()->SetTitle("score btagDDCvB");
    mg->GetYaxis()->SetTitle("score btagDDCvL");
    
    mg->GetXaxis()->SetRangeUser(0.,1.);
    //mg->SetLogy();
    mg->GetYaxis()->SetRangeUser(0.,1.);
	mg->Draw("ALP");
     canvas2->SaveAs("fatJet_btagDDBvL_roc.pdf");
                delete canvas2;
  
	
 TCanvas *canvas3 = new TCanvas("canvas3", "", 1900, 1100);
                canvas3->cd();
    
    
TMultiGraph  *mg1  = new TMultiGraph();
	
	
	TGraph g2(TGraph(ccvsLight.size(), &ccvsLight[0], &bbvsLight[0]));
	g2.SetLineColor(kRed);
	mg1->Add(&g2);
    mg1->GetXaxis()->SetTitle("ccvsLight");
  mg1->GetYaxis()->SetTitle("bbvsLight");
	mg1->Draw("ALP");
     canvas3->SaveAs("fatJet_deepTagMD_ccvsLight_roc.pdf");
                delete canvas3;

    
    
    
    
/*
 TCanvas *canvas1 = new TCanvas("canvas1", "", 1900, 1100);
                canvas1->cd();

               fatJet_ParticleNetMD_probQCDb->GetYaxis()->SetTitle("NEvents");
               fatJet_ParticleNetMD_probQCDb->GetYaxis()->SetTitleSize(0.05);
               fatJet_ParticleNetMD_probQCDb->GetYaxis()->SetTitleOffset(0.80);
               fatJet_ParticleNetMD_probQCDb->GetXaxis()->SetTitle("score");
               fatJet_ParticleNetMD_probQCDb->GetXaxis()->SetTitleSize(0.05);
               fatJet_ParticleNetMD_probQCDb->GetXaxis()->SetTitleOffset(0.90);
               fatJet_ParticleNetMD_probQCDb->SetLineWidth(2);
               fatJet_ParticleNetMD_probQCDb->SetLineColor(kOrange+7);
 gPad->SetGrid();
    
               fatJet_ParticleNetMD_probQCDb->SetTitle("fatJet_ParticleNetMD");
    
               fatJet_ParticleNetMD_probQCDb->SetStats(0);
	



	fatJet_ParticleNetMD_probQCDb->Draw("hist");
        TLegend* leg1 = new TLegend(0.50,0.70,0.890,0.89);
        leg1->AddEntry( fatJet_ParticleNetMD_probQCDb, "fatJet_ParticleNetMD_probQCDbb");

        leg1->SetTextSize(0.03);
        leg1->SetTextFont(42);
        leg1->SetBorderSize(0);
        leg1->Draw("same");
        canvas1->SaveAs("fatJet_ParticleNetMD_QCDb_roc.pdf");
                delete canvas1;
*/


      
        



fatJet_ParticleNetMD_probQCDb->Write();       
fatJet_ParticleNetMD_probQCDbb->Write();         
fatJet_ParticleNetMD_probQCDc->Write();          
fatJet_ParticleNetMD_probQCDcc->Write();         
fatJet_ParticleNetMD_probQCDothers->Write();     
fatJet_ParticleNetMD_probXbb->Write();           
fatJet_ParticleNetMD_probXcc->Write();           
fatJet_ParticleNetMD_probXqq->Write();           
fatJet_btagCMVA->Write();                        
fatJet_btagCSVV2->Write();                       
fatJet_btagDDBvL->Write();                       
fatJet_btagDDBvL_noMD->Write();                 
fatJet_btagDDCvB->Write();                  
fatJet_btagDDCvB_noMD->Write();                  
fatJet_btagDDCvL->Write();                 
fatJet_btagDDCvL_noMD->Write();                  
fatJet_btagDeepB->Write();                       
fatJet_btagHbb->Write();                         
fatJet_deepTagMD_H4qvsQCD->Write();              
fatJet_deepTagMD_HbbvsQCD->Write();              
fatJet_deepTagMD_TvsQCD->Write();                
fatJet_deepTagMD_WvsQCD->Write();                
fatJet_deepTagMD_ZHbbvsQCD->Write();             
fatJet_deepTagMD_ZHccvsQCD->Write();             
fatJet_deepTagMD_ZbbvsQCD->Write();              
fatJet_deepTagMD_ZvsQCD->Write();                
fatJet_deepTagMD_bbvsLight->Write();             
fatJet_deepTagMD_ccvsLight->Write();             
fatJet_deepTag_H->Write();                      
fatJet_deepTag_QCD->Write();                   
fatJet_deepTag_QCDothers->Write();             
fatJet_deepTag_TvsQCD->Write();                   
fatJet_deepTag_WvsQCD->Write();                 
fatJet_deepTag_ZvsQCD->Write();                
aK15Puppi_ParticleNetMD_probQCD->Write();         
aK15Puppi_ParticleNetMD_probXbb->Write();         
aK15Puppi_ParticleNetMD_probXcc->Write();         
aK15Puppi_ParticleNetMD_probXqq->Write();         
aK15Puppi_btagCSVV2->Write();                     
aK15Puppi_btagDeepB->Write();                     
aK15Puppi_btagJP->Write();                        


    
    
    
    
    



fnew->Close();

}
