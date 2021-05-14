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

 

void draw1d3hist(TH1F* h1, string nc=""){

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



void roccurves_2L_iteration() {


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
    
    
    TH1F *fatJet_ParticleNetMD_probXbbvsQCD= new TH1F("tagger6vsQCD","",50,0,1);  
    
     TH1F *fatJet_ParticleNetMD_probXccvsQCD= new TH1F("tagger7vsQCD","",50,0,1);  
    
     TH1F *fatJet_ParticleNetMD_probUDSvsQCD= new TH1F("tagger8vsQCD","",50,0,1);  
    
    TH1F *fatJet_ParticleNetMD_probWvsQCD= new TH1F("taggerWvsQCD","",50,0,1);  
    
    TH1F * fatJet_ParticleNetMD_probQCD= new TH1F("QCD","",50,0,1); 
    
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



TH1F *fatJet_ParticleNetMD_probQCDc_01= new TH1F("tagger_FatJet_ParticleNetMD_probQCDc_01","",50,0,1);
    



TH1F *fatJet_deepTagMD_bbvsLight_01= new TH1F("tagger_FatJet_deepTagMD_bbvsLight_01","",50,0,1);



TH1F *fatJet_deepTagMD_ccvsLight_01= new TH1F("tagger_FatJet_deepTagMD_ccvsLight_01","",50,0,1);


    
    
TH1F *fat_btagDDBvL_01= new TH1F("tagger_btagDDBvL_01","",50,0,1);


        
TH1F *fat_btagDDCvL_01= new TH1F("tagger_btagDDCvL_01","",50,0,1);

 
        
TH1F *fat_btagDDCvB_01= new TH1F("tagger_btagDDCvB_01","",50,0,1);
   
    
    
 
    
    
    
    
    
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
        
        
      /*  
        Float_t FatJet_ParticleNetMD_probQCDb,       
FatJet_ParticleNetMD_probQCDbb,      
FatJet_ParticleNetMD_probQCDc,       
FatJet_ParticleNetMD_probQCDcc,      
FatJet_ParticleNetMD_probQCDothers,
        */
       
         fatJet_ParticleNetMD_probQCD->Fill(FatJet_ParticleNetMD_probQCDb+FatJet_ParticleNetMD_probQCDbb+FatJet_ParticleNetMD_probQCDc+FatJet_ParticleNetMD_probQCDcc+FatJet_ParticleNetMD_probQCDothers);
        
     fatJet_ParticleNetMD_probXbbvsQCD->Fill(FatJet_ParticleNetMD_probXbb/(FatJet_ParticleNetMD_probXbb+FatJet_ParticleNetMD_probQCDb+FatJet_ParticleNetMD_probQCDbb+FatJet_ParticleNetMD_probQCDc+FatJet_ParticleNetMD_probQCDcc+FatJet_ParticleNetMD_probQCDothers));
    
     fatJet_ParticleNetMD_probXccvsQCD->Fill(FatJet_ParticleNetMD_probXcc/(FatJet_ParticleNetMD_probXcc+FatJet_ParticleNetMD_probQCDb+FatJet_ParticleNetMD_probQCDbb+FatJet_ParticleNetMD_probQCDc+FatJet_ParticleNetMD_probQCDcc+FatJet_ParticleNetMD_probQCDothers));
    
     fatJet_ParticleNetMD_probUDSvsQCD->Fill(FatJet_ParticleNetMD_probXqq/(FatJet_ParticleNetMD_probXqq+FatJet_ParticleNetMD_probQCDb+FatJet_ParticleNetMD_probQCDbb+FatJet_ParticleNetMD_probQCDc+FatJet_ParticleNetMD_probQCDcc+FatJet_ParticleNetMD_probQCDothers));
    
     fatJet_ParticleNetMD_probWvsQCD->Fill((FatJet_ParticleNetMD_probXbb+FatJet_ParticleNetMD_probXcc)/(FatJet_ParticleNetMD_probXcc+FatJet_ParticleNetMD_probXbb+FatJet_ParticleNetMD_probQCDb+FatJet_ParticleNetMD_probQCDbb+FatJet_ParticleNetMD_probQCDc+FatJet_ParticleNetMD_probQCDcc+FatJet_ParticleNetMD_probQCDothers));
        
        
        
       
        
        
    }  
        
        // if sigHist is a TH1F* and bkgHist is a TH1F* (with the same number of bins!)

int nbins = fatJet_btagDDCvB->GetNbinsX();

// get the total integrals for each histogram
float fatJet_btagDDCvB_integral = fatJet_btagDDCvB->Integral(1,nbins);
float fatJet_btagDDBvL_integral = fatJet_btagDDBvL->Integral(1,nbins);

float fatJet_btagDDBvL_noMD_integral = fatJet_btagDDBvL_noMD->Integral(1,nbins);    
float fatJet_btagDDCvB_noMD_integral = fatJet_btagDDCvB_noMD->Integral(1,nbins);
float fatJet_btagDDCvL_noMD_integral = fatJet_btagDDCvL_noMD->Integral(1,nbins); 
float fatJet_btagDDCvL_integral = fatJet_btagDDCvL->Integral(1,nbins);   

// create containers sig = x points, bkg = y points
std::vector<float> fatJet_btagDDCvBPoints(nbins);
std::vector<float> fatJet_btagDDBvLPoints(nbins);
 
std::vector<float> fatJet_btagDDBvL_noMDPoints(nbins);
std::vector<float> fatJet_btagDDCvB_noMDPoints(nbins);
std::vector<float> fatJet_btagDDCvL_noMDPoints(nbins);
std::vector<float> fatJet_btagDDCvLPoints(nbins);    

for ( int i = 0; i < nbins; ++i ) {
 
  float fatJet_btagDDCvB_slice_integral = fatJet_btagDDCvB->Integral(nbins-i,nbins);
  float fatJet_btagDDBvL_slice_integral = fatJet_btagDDBvL->Integral(nbins-i,nbins);
    
    float fatJet_btagDDBvL_noMD_slice_integral = fatJet_btagDDBvL_noMD->Integral(nbins-i,nbins);
    float fatJet_btagDDCvB_noMD_slice_integral = fatJet_btagDDCvB_noMD->Integral(nbins-i,nbins);
    float fatJet_btagDDCvL_noMD_slice_integral = fatJet_btagDDCvL_noMD->Integral(nbins-i,nbins);
    float fatJet_btagDDCvL_slice_integral = fatJet_btagDDCvL->Integral(nbins-i,nbins);
    
    
  fatJet_btagDDCvBPoints.push_back(fatJet_btagDDCvB_slice_integral/fatJet_btagDDCvB_integral);
  fatJet_btagDDBvLPoints.push_back(fatJet_btagDDBvL_slice_integral/fatJet_btagDDBvL_integral);
    
    fatJet_btagDDBvL_noMDPoints.push_back(fatJet_btagDDBvL_noMD_slice_integral/fatJet_btagDDBvL_noMD_integral);
    fatJet_btagDDCvB_noMDPoints.push_back(fatJet_btagDDCvB_noMD_slice_integral/fatJet_btagDDCvB_noMD_integral);
    fatJet_btagDDCvL_noMDPoints.push_back(fatJet_btagDDCvL_noMD_slice_integral/fatJet_btagDDCvL_noMD_integral);
    fatJet_btagDDCvLPoints.push_back(fatJet_btagDDCvL_slice_integral/fatJet_btagDDCvL_integral);
    
}
   
    
   
    
    int Nbins = fatJet_ParticleNetMD_probXbbvsQCD->GetNbinsX();
    
    
float fatJet_ParticleNetMD_probXbbvsQCD_integral = fatJet_ParticleNetMD_probXbbvsQCD->Integral(1,Nbins);
float fatJet_ParticleNetMD_probXccvsQCD_integral = fatJet_ParticleNetMD_probXccvsQCD->Integral(1,Nbins);

float fatJet_ParticleNetMD_probUDSvsQCD_integral = fatJet_ParticleNetMD_probUDSvsQCD->Integral(1,Nbins);    
float fatJet_ParticleNetMD_probWvsQCD_integral = fatJet_ParticleNetMD_probWvsQCD->Integral(1,Nbins);
float fatJet_ParticleNetMD_probQCD_integral = fatJet_ParticleNetMD_probQCD->Integral(1,Nbins); 
     
 // create containers sig = x points, bkg = y points

std::vector<float> fatJet_ParticleNetMD_probXbbvsQCDPoints(Nbins);
std::vector<float> fatJet_ParticleNetMD_probXccvsQCDPoints(Nbins);
std::vector<float> fatJet_ParticleNetMD_probUDSvsQCDPoints(Nbins);
std::vector<float> fatJet_ParticleNetMD_probWvsQCDPoints(Nbins);
std::vector<float> fatJet_ParticleNetMD_probQCDPoints(Nbins);    
    
    
    for ( int i = 0; i < Nbins; ++i ) {
 

    
    float fatJet_ParticleNetMD_probXbbvsQCD_slice_integral = fatJet_ParticleNetMD_probXbbvsQCD->Integral(Nbins-i,Nbins);
    float fatJet_ParticleNetMD_probXccvsQCD_slice_integral = fatJet_ParticleNetMD_probXccvsQCD->Integral(Nbins-i,Nbins);
    float fatJet_ParticleNetMD_probUDSvsQCD_slice_integral = fatJet_ParticleNetMD_probUDSvsQCD->Integral(Nbins-i,Nbins);
    float fatJet_ParticleNetMD_probWvsQCD_slice_integral = fatJet_ParticleNetMD_probWvsQCD->Integral(Nbins-i,Nbins);
     float fatJet_ParticleNetMD_probQCD_slice_integral = fatJet_ParticleNetMD_probQCD->Integral(Nbins-i,Nbins);
        
    fatJet_ParticleNetMD_probXbbvsQCDPoints.push_back(fatJet_ParticleNetMD_probXbbvsQCD_slice_integral/fatJet_ParticleNetMD_probXbbvsQCD_integral);
    fatJet_ParticleNetMD_probXccvsQCDPoints.push_back(fatJet_ParticleNetMD_probXccvsQCD_slice_integral/fatJet_ParticleNetMD_probXccvsQCD_integral);
    fatJet_ParticleNetMD_probUDSvsQCDPoints.push_back(fatJet_ParticleNetMD_probUDSvsQCD_slice_integral/fatJet_ParticleNetMD_probUDSvsQCD_integral);
    fatJet_ParticleNetMD_probWvsQCDPoints.push_back(fatJet_ParticleNetMD_probWvsQCD_slice_integral/fatJet_ParticleNetMD_probWvsQCD_integral);
     fatJet_ParticleNetMD_probQCDPoints.push_back(fatJet_ParticleNetMD_probQCD_slice_integral/fatJet_ParticleNetMD_probQCD_integral);
    
    
    
    
}

     /*
     TH1F *fatJet_ParticleNetMD_probXbbvsQCD= new TH1F("tagger6vsQCD","",50,0,1);  
    
     TH1F *fatJet_ParticleNetMD_probXccvsQCD= new TH1F("tagger7vsQCD","",50,0,1);  
    
     TH1F *fatJet_ParticleNetMD_probUDSvsQCD= new TH1F("tagger8vsQCD","",50,0,1);  
    
    TH1F *fatJet_ParticleNetMD_probWvsQCD= new TH1F("taggerWvsQCD","",50,0,1);  
    
    TH1F * fatJet_ParticleNetMD_probQCD= new TH1F("QCD","",50,0,1); 
    */
     TCanvas *canvas1 = new TCanvas("canvas1", "", 1900, 1100);
    canvas1->cd(); 
     TGraph *d = new TGraph(fatJet_ParticleNetMD_probXbbvsQCDPoints.size(),&fatJet_ParticleNetMD_probXbbvsQCDPoints[0],&fatJet_ParticleNetMD_probXccvsQCDPoints[0]);
     TGraph *d1 = new TGraph(fatJet_ParticleNetMD_probXccvsQCDPoints.size(),&fatJet_ParticleNetMD_probUDSvsQCDPoints[0],&fatJet_ParticleNetMD_probXccvsQCDPoints[0]);
     TGraph *d2 = new TGraph(fatJet_ParticleNetMD_probXccvsQCDPoints.size(),&fatJet_ParticleNetMD_probWvsQCDPoints[0],&fatJet_ParticleNetMD_probXccvsQCDPoints[0]);
     
     d->SetName("d");
     d1->SetName("d1");
     d2->SetName("d2");    
    d->SetTitle("ROC");
    d->SetLineColor(kBlack);
    d->GetYaxis()->SetTitle("signal efficiency");
    d->GetXaxis()->SetTitle("background efficiency");
    d->Draw("ALP");        
    d1->SetLineColor(kRed);   
	d1->Draw("same"); 
    d2->SetLineColor(kGreen);
    d2->Draw("same");
    gPad->SetLogy();
   auto legend1 = new TLegend( 0.2,0.2,0.4,0.4);

   legend1->SetHeader("ParticleNet Taggers","C"); // option "C" allows to center the header
  legend1->AddEntry("d","XccvQCD vs. XbbvQCD ");
   legend1->AddEntry("d1","XccvQCD vs. udsvQCD");
   legend1->AddEntry("d2","XccvQCD vs. WvQCD");
  
   legend1->Draw();
               
 canvas1->SaveAs("fatJet_ParticleNet_roc.pdf");
                delete canvas1;
    
    
   /* 
    
    TCanvas *canvas3 = new TCanvas("canvas3", "", 1900, 1100);
    canvas3->cd();    

     TGraph *g = new TGraph(fatJet_btagDDCvBPoints.size(),&fatJet_btagDDBvLPoints[0],&fatJet_btagDDCvBPoints[0]);
     
    
    
    g->SetTitle("ROC");
    g->SetLineColor(kRed);
    g->GetXaxis()->SetTitle("btagDDBvL");
    g->GetYaxis()->SetTitle("btagDDCvB");
	g->Draw("ALP"); 
       
               
 canvas3->SaveAs("fatJet_BvLvsCvB_roc.pdf");
                delete canvas3;
    
    */
    
    
    
    
    TCanvas *canvas = new TCanvas("canvas", "", 1900, 1100);
    canvas->cd(); 
     TGraph *g = new TGraph(fatJet_btagDDCvBPoints.size(),&fatJet_btagDDBvLPoints[0],&fatJet_btagDDCvBPoints[0]);
     TGraph *g1 = new TGraph(fatJet_btagDDBvLPoints.size(),&fatJet_btagDDBvLPoints[0],&fatJet_btagDDCvLPoints[0]);
     TGraph *g2 = new TGraph(fatJet_btagDDBvL_noMDPoints.size(),&fatJet_btagDDBvL_noMDPoints[0],&fatJet_btagDDCvL_noMDPoints[0]);
     TGraph *g3 = new TGraph(fatJet_btagDDCvB_noMDPoints.size(),&fatJet_btagDDCvL_noMDPoints[0],&fatJet_btagDDCvB_noMDPoints[0]);
     TGraph *g4 = new TGraph(fatJet_btagDDCvBPoints.size(),&fatJet_btagDDCvLPoints[0],&fatJet_btagDDCvBPoints[0]);
     g->SetName("g");
     g1->SetName("g1");
     g2->SetName("g2");
     g3->SetName("g3");
     g4->SetName("g4");
    
    
    g->SetTitle("ROC");
    g->SetLineColor(kBlack);
    g->GetYaxis()->SetTitle("signal efficiency");
    g->GetXaxis()->SetTitle("background efficiency");
    g->Draw("ALP");
        
    g1->SetLineColor(kRed);
   
	g1->Draw("same"); 
     
    g2->SetLineColor(kGreen);
    g2->Draw("same"); 
    
    g3->SetLineColor(kBlue);
    g3->Draw("same"); 
    
    g4->SetLineColor(kOrange);
    g4->Draw("same"); 
    gPad->SetLogy();
    // auto legend = new TLegend(0.1,0.90,0.27,0.6);
   auto legend = new TLegend( 0.2,0.2,0.4,0.4);
   legend->SetHeader("Taggers","C"); // option "C" allows to center the header
  legend->AddEntry("g","CvB vs. BvL");
   legend->AddEntry("g1","CvL vs. BvL");
   legend->AddEntry("g2","CvL_noMD vs. BvL_noMD");
   legend->AddEntry("g3","CvB_noMD vs. CvL_noMD");
   legend->AddEntry("g4","CvB vs. CvL");
   legend->Draw();
               
 canvas->SaveAs("fatJet_CvB_roc.pdf");
                delete canvas;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

 /* for(int l=0;l<numybins;l++){


                                std::vector<TF2*> f2;

                                NameString.str("");
                               NameString<<"fat_btagDDBvL"<<ybins[l]<<"_"<<ybins[l+1];
                                draw1d3hist(hist[NameString.str().c_str()]);
        
        
        
         }
        
        
        
        
        */
        
        
fnew->Close();

}