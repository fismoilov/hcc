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
#include <TDatabasePDG.h>
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

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////analysis for AK15 AK8 AK4 Jets//////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


void akjet_efficiency_2L() {


Float_t  LeadAK8_pt[50],LeadAK8_eta[50],LeadAK8_phi[50],LeadAK8_mass[50],LeadAK8_msoftdrop[50], HJ1_eta[50],HJ1_phi[50],HJ1_pt[50],HJ2_eta[50],
HJ2_phi[50], HJ2_pt[50],HJ1_HJ2_dR_noFSR[50],Jet_eta[50],Jet_mass[50],Jet_phi[50],Jet_pt[50],H_pt[50],H_eta[50],H_phi[50];

Float_t GenPart_eta[150],GenPart_mass[150],GenPart_phi[150],GenPart_pt[150];
 int  GenPart_genPartIdxMother[150], GenPart_pdgId[150], GenPart_status[150], GenPart_statusFlags[150];
UInt_t nGenPart;
Float_t LeadAK15_pt[50],LeadAK15_eta[50],LeadAK15_phi[50],LeadAK15_mass[50],LeadAK15_msoftdrop[50];
Float_t  nLeadAK15;

Float_t Electron_pt[50],Electron_phi[50],Electron_eta[50],Electron_mass[50], Muon_pt[50], Muon_eta[50], Muon_phi[50], Muon_mass[50]; 




Float_t JetPt_1[50],JetPt_2[50];
Float_t DeepJet_CvsL_1[50], DeepJet_CvsB_1[50],DeepJet_CvsL_2[50],DeepJet_CvsB_2[50],H_mass[50],V_mass[50],V_pt[50],HVdPhi[50];

Float_t  nLeadAK8;
UInt_t nJet;

int isZmm,isZee;
int controlSample;
int hJetInd1,hJetInd2;
bool twoResolvedJets;
UInt_t  nElectron, nMuon;



	TChain *t = new TChain("Events");
 // TChain *chain = new TChain("tree");
   for (int fileNum=1;fileNum < 12;fileNum++) {
      t->AddFile(Form("output_ZH125ToCC_ZLL_powheg_%d.root", fileNum));
//std::cout<<"File number = "<<fileNum<<std::endl; 
 }
  t->SetBranchAddress("twoResolvedJets",&twoResolvedJets);
        t->SetBranchAddress("isZmm",&isZmm);
        t->SetBranchAddress("isZee",&isZee);
        t->SetBranchAddress("JetPt_1",JetPt_1);
        t->SetBranchAddress("JetPt_2",JetPt_2);
        t->SetBranchAddress("DeepJet_CvsL_1",DeepJet_CvsL_1);
        t->SetBranchAddress("DeepJet_CvsL_2",DeepJet_CvsL_2);
        t->SetBranchAddress("DeepJet_CvsB_1",DeepJet_CvsB_1);
        t->SetBranchAddress("DeepJet_CvsB_2",DeepJet_CvsB_2);
        t->SetBranchAddress("hJetInd1",&hJetInd1);
        t->SetBranchAddress("hJetInd2",&hJetInd2);

        t->SetBranchAddress( "GenPart_pt",GenPart_pt);
        t->SetBranchAddress( "GenPart_phi",GenPart_phi);
        t->SetBranchAddress( "GenPart_eta",GenPart_eta);
        t->SetBranchAddress( "nGenPart",&nGenPart);
        t->SetBranchAddress( "GenPart_mass",GenPart_mass);

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
        t->SetBranchAddress("H_mass",H_mass);
        t->SetBranchAddress("V_pt",V_pt);
        t->SetBranchAddress("HVdPhi",HVdPhi);
        t->SetBranchAddress("V_mass",V_mass);

        t->SetBranchAddress("HJ1_HJ2_dR_noFSR",HJ1_HJ2_dR_noFSR);

	t->SetBranchAddress("LeadAK8_pt",LeadAK8_pt);
	t->SetBranchAddress("LeadAK8_eta",LeadAK8_eta);
	t->SetBranchAddress("LeadAK8_phi",LeadAK8_phi);
	t->SetBranchAddress("LeadAK8_mass",LeadAK8_mass);
	t->SetBranchAddress("nLeadAK8",&nLeadAK8);
        t->SetBranchAddress("LeadAK8_msoftdrop",&LeadAK8_msoftdrop);

        t->SetBranchAddress("LeadAK15_pt",LeadAK15_pt);
        t->SetBranchAddress("LeadAK15_eta",LeadAK15_eta);
        t->SetBranchAddress("LeadAK15_phi",LeadAK15_phi);
        t->SetBranchAddress("LeadAK15_mass",LeadAK15_mass);
	t->SetBranchAddress("LeadAK15_msoftdrop",LeadAK15_msoftdrop);
        t->SetBranchAddress("nLeadAK15",&nLeadAK15);

        t->SetBranchAddress("nJet",&nJet);
        t->SetBranchAddress("Jet_eta",Jet_eta);
        t->SetBranchAddress("Jet_phi",Jet_phi);
        t->SetBranchAddress("Jet_pt",Jet_pt);
        t->SetBranchAddress("Jet_mass",Jet_mass);


        t->SetBranchAddress("H_eta",H_eta);
        t->SetBranchAddress("H_phi",H_phi);
        t->SetBranchAddress("H_pt",H_pt);
        t->SetBranchAddress("H_mass",H_mass);

        t->SetBranchAddress("HJ1_eta",HJ1_eta);
        t->SetBranchAddress("HJ1_phi",HJ1_phi);
        t->SetBranchAddress("HJ1_pt",HJ1_pt);

        t->SetBranchAddress("HJ2_eta",HJ2_eta);
        t->SetBranchAddress("HJ2_phi",HJ2_phi);
        t->SetBranchAddress("HJ2_pt",HJ2_pt);



//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////To plot histogramms read  .root file below////////////////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////



	TFile *fnew = new TFile("histogrammi2.root", "RECREATE");


	TH1F *ak_eta = new TH1F("LeadAK8_eta","",100,-3.5,3.5);
	TH1F *ak_phi = new TH1F("LeadAK8_phi","",100,-3.5,3.5);
	TH1F *ak_pt = new TH1F("LeadAK8_pt","",20,0,1000);
	TH1F *ak_mass = new TH1F("LeadAK8_mass","",720,20,300);
	TH1F *ak_msoftdrop= new TH1F("LeadAK8_msoftdrop","",50,0,300);
        TH1F *ak_num = new TH1F("nLeadAK8","",100,0,50);

        TH1F *ak15_eta = new TH1F("LeadAK15_eta","",100,-3.5,3.5);
        TH1F *ak15_phi = new TH1F("LeadAK15_phi","",100,-3.5,3.5);
        TH1F *ak15_pt = new TH1F("LeadAK15_pt","",20,0,1000);
     
 TH1F *genhiggspt=new TH1F("genhiggspt","",16,200,1000);
 TH1F *h_pt=new TH1F("h_pt","",100,0,1000);


 TH1F *higgs_mass_50_100_pt=new TH1F("higgs_mass_50_100_pt","",100,50,200);
 TH1F *higgs_mass_100_200_pt=new TH1F("higgs_mass_100_200_pt","",100,50,200);
 TH1F *higgs_mass_200_300_pt=new TH1F("higgs_mass_200_300_pt","",100,50,200);
 TH1F *higgs_mass_300_400_pt=new TH1F("higgs_mass_300_400_pt","",100,50,200);
 TH1F *higgs_mass_400_500_pt=new TH1F("higgs_mass_400_500_pt","",100,50,200);
 TH1F *higgs_mass_500_600_pt=new TH1F("higgs_mass_500_600_pt","",100,50,200);
 TH1F *higgs_mass_600_pt=new TH1F("higgs_mass_600_pt","",100,50,200);








	TH1F *ak15_mass = new TH1F("LeadAK15_mass","",50,20,300);
        TH1F *ak15_msoftdrop = new TH1F("LeadAK15_msoftdrop","",50,0,300);
	
	TH2F *ak8_etaphi = new TH2F("ak8_etaphi","",50,-2.4,2.4,50,-3.5,3.5);
        TH2F *ak15_etaphi = new TH2F("ak15_etaphi","",50,-2.4,2.4,50,-3.5,3.5);



        TH1F *jet_eta = new TH1F("Jet_eta","",100,-3.5,3.5);
        TH1F *jet_phi = new TH1F("Jet_phi","",100,-3.5,3.5);
        TH1F *jet_pt = new TH1F("Jet_pt","",100,0,600);
        TH1F *jet_mass = new TH1F("Jet_mass","",100,20,600);

        TH1F *zmumu_mass = new TH1F("zmumu_mass","",100,20,300);
        TH1F *hja_eta = new TH1F("HJ1_eta","",100,-3.5,3.5);
        TH1F *hja_phi = new TH1F("HJ1_phi","",100,-3.5,3.5);
        TH1F *hja_pt = new TH1F("HJ1_pt","",100,0,600);
       
        TH1F *hjb_eta = new TH1F("HJ2_eta","",100,-3.5,3.5);
        TH1F *hjb_phi = new TH1F("HJ2_phi","",100,-3.5,3.5);
        TH1F *hjb_pt = new TH1F("HJ2_pt","",100,0,600);
        TH1F *hj_mass = new TH1F("hj_mass","",100,0,600);
        TH2F *ak_jet = new TH2F("ak8_jet","",100,20,300,100,20,300);
        TH2F *hja = new TH2F("hja","",50,-2.4,2.4,50,-3.5,3.5);
        TH2F *hjb = new TH2F("hjb","",50,-2.4,2.4,50,-3.5,3.5);

  TH2F *higgs_ak8 = new TH2F("higgs_ak8","",100,350,700,100,350,700);

 TH2F *true_higgs_ak15 = new TH2F("true_higgs_ak15","",1000,200,1000,1000,200,1000);
 TH2F *true_higgs_ak4 = new TH2F("true_higgs_ak4","",1000,200,1000,1000,200,1000);




////For Efficiency/////////////


TH1F *genhiggsall= new TH1F("genhiggsall","",20,0,1000);
TH1F *genhiggsall_jet= new TH1F("genhiggsall_jet","",20,0,1000);

TH1F *genhiggs_with_ak15= new TH1F("genhiggs_with_ak15","",20,0,1000);
TH1F *genhiggs_with_ak4= new TH1F("genhiggs_with_ak4","",20,0,1000);
TH1F *genhiggs_with_ak8= new TH1F("genhiggs_with_ak8","",20,0,1000);



///////refining the binning///////



   const Int_t NBINS = 14;
   float  edges[NBINS + 1] = {1,30,50,100,150,200,250,300,350,400,450,500,600,700,1000};



TH1F *genhiggsall_varbin= new TH1F("genhiggsall_varbin","",NBINS,edges);
TH1F *genhiggsall_jet_varbin= new TH1F("genhiggsall_jet_varbin","",NBINS,edges);
TH1F *genhiggs_with_ak15_varbin= new TH1F("genhiggs_with_ak15_varbin","",NBINS,edges);
TH1F *genhiggs_with_ak4_varbin= new TH1F("genhiggs_with_ak4_varbin","",NBINS,edges);
TH1F *genhiggs_with_ak8_varbin= new TH1F("genhiggs_with_ak8_varbin","",NBINS,edges);

TH1F *genhiggs_with_ak15_deltaR_0_3_varbin= new TH1F("genhiggs_with_ak15_deltaR_0_3_varbin","",NBINS,edges);
TH1F *genhiggs_with_ak4_deltaR_0_3_varbin= new TH1F("genhiggs_with_ak4_deltaR_0_3_varbin","",NBINS,edges);
TH1F *genhiggs_with_ak8_deltaR_0_3_varbin= new TH1F("genhiggs_with_ak8_deltaR_0_3_varbin","",NBINS,edges);


TH1F *genhiggs_eta_all= new TH1F("genhiggs_eta_all","",20,-2.4,2.4);
TH1F *genhiggs_eta_all_jet= new TH1F("genhiggs_eta_all_jet","",20,-2.4,2.4);
TH1F *genhiggs_eta_with_ak15_deltaR_0_3= new TH1F("genhiggs_eta_with_ak15_deltaR_0_3","",20,-2.4,2.4);
TH1F *genhiggs_eta_with_ak4_deltaR_0_3= new TH1F("genhiggs_eta_with_ak4_deltaR_0_3","",20,-2.4,2.4);
TH1F *genhiggs_eta_with_ak8_deltaR_0_3= new TH1F("genhiggs_eta_with_ak8_deltaR_0_3","",20,-2.4,2.4);


TH1F * genhiggs_with_ak4_no_Vpt= new TH1F("genhiggs_with_ak15_no_Vpt","",20,0,1000);
TH1F * genhiggs_with_ak4_no_sel= new TH1F("genhiggs_with_ak15_no_sel","",20,0,1000);


TH1F *genhiggs_with_ak15_deltaR_0_2=new TH1F("genhiggs_with_ak15_deltaR_0_2","",20,0,1000);
TH1F *genhiggs_with_ak15_deltaR_0_3=new TH1F("genhiggs_with_ak15_deltaR_0_3","",20,0,1000);
TH1F *genhiggs_with_ak15_deltaR_0_5=new TH1F("genhiggs_with_ak15_deltaR_0_5","",20,0,1000);

TH1F *genhiggs_with_ak8_deltaR_0_2=new TH1F("genhiggs_with_ak8_deltaR_0_2","",20,0,1000);
TH1F *genhiggs_with_ak8_deltaR_0_3=new TH1F("genhiggs_with_ak8_deltaR_0_3","",20,0,1000);
TH1F *genhiggs_with_ak8_deltaR_0_5=new TH1F("genhiggs_with_ak8_deltaR_0_5","",20,0,1000);


 TH1F *genhiggs_with_ak4_deltaR_0_2=new TH1F("genhiggs_with_ak4_deltaR_0_2","",20,0,1000);
 TH1F *genhiggs_with_ak4_deltaR_0_3=new TH1F("genhiggs_with_ak4_deltaR_0_3","",20,0,1000);
 TH1F *genhiggs_with_ak4_deltaR_0_5=new TH1F("genhiggs_with_ak4_deltaR_0_5","",20,0,1000);


 TH1F *genhiggs_with_ak4_no_Vpt_deltaR_0_2=new TH1F("genhiggs_with_ak4_no_Vpt_deltaR_0_2","",20,0,1000);
 TH1F *genhiggs_with_ak4_no_Vpt_deltaR_0_3=new TH1F("genhiggs_with_ak4_no_Vpt_deltaR_0_3","",20,0,1000);
 TH1F *genhiggs_with_ak4_no_Vpt_deltaR_0_5=new TH1F("genhiggs_with_ak4_no_Vpt_deltaR_0_5","",20,0,1000);




// const Int_t nbins = 14;
//   float  EDGES[nbins + 1] = {0,30,50,100,150,200,250,300,350,400,450,500,600,700,1000};

TH1F *frac_ak4_ak15=new TH1F("frac_ak4_ak15","",100,0,1000);
TH1F *frac_ak4_ak8=new TH1F("frac_ak4_ak8","",20,0,1000);
TH1F *frac_ak4_ak15_ak4=new TH1F("frac_ak4_ak15_ak4","",20,0,1000);
TH1F *frac_ak4_ak8_ak4=new TH1F("frac_ak4_ak8_ak4","",20,0,1000);
TH1F *ak15_nosel_pt=new TH1F("ak15_nosel_pt","",20,0,1000);
TH1F *ak8_nosel_pt=new TH1F("ak8_nosel_pt","",20,0,1000);







// TH2F *dR_ak8_ak4 = new TH2F("dR_ak8_ak4","",50,0,5,50,0,5);       

// TH2F *dR_ak15_ak4 = new TH2F("dR_ak15_ak4","",50,0,5,50,0,5);

 TH1F *dR_ak8_ak4 = new TH1F("dR_ak8_ak4","",100,0,5); 
 TH1F *dR_ak15_ak4 = new TH1F("dR_ak15_ak4","",100,0,5);

 TH1F *dR_Zee_ak8 = new TH1F("dR_Zee_ak8","",100,0,5);
 TH1F *dR_Zmumu_ak8 = new TH1F("dR_Zmumu_ak8","",100,0,5);
 TH1F *dR_Zee_ak15 = new TH1F("dR_Zee_ak15","",100,0,5);
 TH1F *dR_Zmumu_ak15 = new TH1F("dR_Zmumu_ak15","",100,0,5);

TH1F *dR_ak8_electron = new TH1F("dR_ak8_electron","",100,0,5);
TH1F *dR_ak8_muon = new TH1F("dR_ak8_muon","",100,0,5);
TH1F *dR_ak15_electron = new TH1F("dR_ak15_electron","",100,0,5);
TH1F *dR_ak15_muon = new TH1F("dR_ak15_muon","",100,0,5);

TH1F *dR_ak8_jet = new TH1F("dR_ak8_jet","",100,0,5);


TH1F *ak8_pt_dR_above0_5=new TH1F("ak8_pt_dR_above0_5","",100,200,1000);
TH1F *ak8_pt_dR_less0_5=new TH1F("ak8_pt_dR_less0_5","",100,200,1000);
TH1F *ak8_msoftdrop_dR_above0_5=new TH1F("ak8_msoftdrop_dR_above0_5","",100,50,300);
TH1F *ak8_msoftdrop_dR_less0_5=new TH1F("ak8_msoftdrop_dR_less0_5","",100,50,300);

 TH1F *dR_ak8_ak4_above0_5 = new TH1F("dR_ak8_ak4_above0_5","",100,0,5);
 TH1F *dR_ak15_ak4_above0_5 = new TH1F("dR_ak15_ak4_above0_5","",100,0,5);
 TH1F *dR_ak8_ak4_less0_5 = new TH1F("dR_ak8_ak4_less0_5","",100,0,5);
 TH1F *dR_ak15_ak4_less0_5 = new TH1F("dR_ak15_ak4_less0_5","",100,0,5);




 TH1F * dR_GenHiggs_ak8=new TH1F("dR_GenHiggs_ak8","",100,0,5);
 TH1F * dR_GenHiggs_ak15=new TH1F("dR_GenHiggs_ak15","",100,0,5);
 TH1F * dR_GenHiggs_ak4=new TH1F("dR_GenHiggs_ak4","",100,0,5);





TH1F *ak15_pt_dR_above0_5=new TH1F("ak15_pt_dR_above0_5","",100,200,1000);
TH1F *ak15_pt_dR_less0_5=new TH1F("ak15_pt_dR_less0_5","",100,200,1000);
TH1F *ak15_msoftdrop_dR_above0_5=new TH1F("ak15_msoftdrop_dR_above0_5","",100,50,300);
TH1F *ak15_msoftdrop_dR_less0_5=new TH1F("ak15_msoftdrop_dR_less0_5","",100,50,300);


TH1F *ak_msoftdrop_pt250= new TH1F("ak_msoftdrop_pt250","",50,0,300);
TH1F *ak_msoftdrop_pt300= new TH1F("ak_msoftdrop_pt300","",50,0,300);

TH1F *ak15_msoftdrop_pt250= new TH1F("ak15_msoftdrop_pt250","",50,0,300);
TH1F *ak15_msoftdrop_pt300= new TH1F("ak15_msoftdrop_pt300","",50,0,300);


TH1F *ak_msoftdrop_200_300_pt=new TH1F("ak_msoftdrop_200_300_pt","",50,0,300);
TH1F *ak_msoftdrop_300_400_pt=new TH1F("ak_msoftdrop_300_400_pt","",50,0,300);
TH1F *ak_msoftdrop_400_500_pt=new TH1F("ak_msoftdrop_400_500_pt","",50,0,300);
TH1F *ak_msoftdrop_500_pt=new TH1F("ak_msoftdrop_500_pt","",50,0,300);


TH1F *ak15_msoftdrop_200_300_pt=new TH1F("ak15_msoftdrop_200_300_pt","",50,0,300);
TH1F *ak15_msoftdrop_300_400_pt=new TH1F("ak15_msoftdrop_300_400_pt","",50,0,300);
TH1F *ak15_msoftdrop_400_500_pt=new TH1F("ak15_msoftdrop_400_500_pt","",50,0,300);
TH1F *ak15_msoftdrop_500_pt=new TH1F("ak15_msoftdrop_500_pt","",50,0,300);

/////Fitting function













int ak8jets=0;
int ak4jets=0;
int ak15jets=0;



//int ak4_selected_true=0;
//int ak4_selected_boosted_wrt_ak8=0;
//int ak4_selected_boosted_wrt_ak15=0;
//int ak8_selected_true=0;
//int ak15_selected_true=0;


/*
const int npar = 1;
TF1 *ftf = new TF1("ftf","fitf", 20, 300,npar);
double f2params[npar]={100};

ftf->SetParameters(f2params);
std::cout<<ftf->GetParameter(0)<<std::endl;
*/




	int NumberOfEvents =68060;
	std::cout<< t->GetEntries()<<std::endl;

	for (int i=0; i<NumberOfEvents;i++) {

		if ((i % 10000) == 0)
		{
			std::cout<<"Event = "<<i<<std::endl;
		}
		t->GetEntry(i);
	


        TLorentzVector LeadAK15;
	TLorentzVector LeadAK8;
        
        TLorentzVector HJ1;
        TLorentzVector HJ2;		
        
        TLorentzVector Hcc;
        TLorentzVector AK8Higgs,AK15Higgs;

        TLorentzVector Zboson;
        TLorentzVector Zee,Zmumu,Jet,Electron1, Electron2, Muon1, Muon2; 

        TLorentzVector GenHiggs;


		LeadAK8.SetPtEtaPhiM(LeadAK8_pt[0],LeadAK8_eta[0],LeadAK8_phi[0],LeadAK8_mass[0]);
		LeadAK15.SetPtEtaPhiM(LeadAK15_pt[0],LeadAK15_eta[0],LeadAK15_phi[0],LeadAK15_mass[0]);

                AK8Higgs.SetPtEtaPhiM(LeadAK8_pt[0],LeadAK8_eta[0],LeadAK8_phi[0],LeadAK8_mass[0]);
                AK15Higgs.SetPtEtaPhiM(LeadAK15_pt[0],LeadAK15_eta[0],LeadAK15_phi[0],LeadAK15_mass[0]);

                Hcc.SetPtEtaPhiM(H_pt[0],H_eta[0],H_phi[0],H_mass[0]);
	
        	HJ1.SetPtEtaPhiM(HJ1_pt[0],HJ1_eta[0],HJ1_phi[0],Jet_mass[0]);
                HJ2.SetPtEtaPhiM(HJ2_pt[0],HJ2_eta[0],HJ2_phi[0],Jet_mass[0]);

  
                Jet.SetPtEtaPhiM(Jet_pt[0],Jet_eta[0],Jet_phi[0],Jet_mass[0]);
                Electron1.SetPtEtaPhiM(Electron_pt[0],Electron_eta[0],Electron_phi[0],Electron_mass[0]);
                Electron2.SetPtEtaPhiM(Electron_pt[1],Electron_eta[1],Electron_phi[1],Electron_mass[1]);
                  
                Muon1.SetPtEtaPhiM(Muon_pt[0],Muon_eta[0],Muon_phi[0],Muon_mass[0]);
                Muon2.SetPtEtaPhiM(Muon_pt[1],Muon_eta[1],Muon_phi[1],Muon_mass[1]);


Zee=Electron1 + Electron2;
Zmumu=Muon1+Muon2;

zmumu_mass->Fill(Zmumu.M());
/*
int events_have_tworesolvedjets=0;
int events_pass_ak4_higgs_selection=0;
int events_ak15_pt_above_200=0;
int events_pass_both_selection=0;


int events_ak15_pt_above_250=0;
int events_ak15_pt_above_300=0;
*/


//////////////////Binning to plot "efficiency"/////////////////////////





///True Higgs////////// ///// ///// ///// ///// ///// ///// ///// 
/////// ///// ///// ///// ///// ///// ///// ///// ///// ///// ///  



 			for (int k=0;k<nGenPart;k++){

 if(GenPart_pdgId[k]==25&&GenPart_genPartIdxMother[k]==0){
	GenHiggs.SetPtEtaPhiM(GenPart_pt[k],GenPart_eta[k],GenPart_phi[k],GenPart_mass[k]);

		genhiggspt->Fill(GenHiggs.Pt());

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////Denominator for resolved topology efficiency////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(GenHiggs.Pt()>0&&(isZmm || isZee)&&JetPt_1[0]>20 &&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&JetPt_2[0]>20 &&V_mass[0]>75 && V_mass[0]<105){
		genhiggsall_jet->Fill(GenHiggs.Pt());
		genhiggs_eta_all_jet->Fill(GenHiggs.Eta());
                genhiggsall_jet_varbin->Fill(GenHiggs.Pt());

                           }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Denominator for boosted topology efficiency/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                 if(GenHiggs.Pt()>0&&(isZmm || isZee)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&V_mass[0]>75 && V_mass[0]<105){                

                 genhiggsall->Fill(GenHiggs.Pt());
                 genhiggsall_varbin->Fill(GenHiggs.Pt());
                 genhiggs_eta_all->Fill(GenHiggs.Eta());
		
          			     }

			   }


true_higgs_ak15->Fill(GenHiggs.Pt(),LeadAK15.Pt());
true_higgs_ak4->Fill(GenHiggs.Pt(),Hcc.Pt());



							// }
						} // end of for loop for GenPart


	if(LeadAK15.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK15.DeltaR(GenHiggs)<0.1){
		genhiggs_with_ak15->Fill(GenHiggs.Pt());
		}

	if(LeadAK8.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK8.DeltaR(GenHiggs)<0.1){
		genhiggs_with_ak8->Fill(GenHiggs.Pt());
		}

	if(LeadAK15.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK15.DeltaR(GenHiggs)<0.2){
		genhiggs_with_ak15_deltaR_0_2->Fill(GenHiggs.Pt());
		}

	if(LeadAK8.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK8.DeltaR(GenHiggs)<0.2){
		genhiggs_with_ak8_deltaR_0_2->Fill(GenHiggs.Pt());
		}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Numerator for boosted topology efficiency/////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	if(LeadAK15.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4  &&V_mass[0]>75&&V_mass[0]<105&&LeadAK15.DeltaR(GenHiggs)<0.3){
		genhiggs_with_ak15_deltaR_0_3->Fill(GenHiggs.Pt());
		genhiggs_with_ak15_deltaR_0_3_varbin->Fill(GenHiggs.Pt());
		genhiggs_eta_with_ak15_deltaR_0_3->Fill(GenHiggs.Eta());
                
 		 }

	if(LeadAK8.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4&&V_mass[0]>75&&V_mass[0]<105&&LeadAK8.DeltaR(GenHiggs)<0.3){
		genhiggs_with_ak8_deltaR_0_3->Fill(GenHiggs.Pt());
		genhiggs_with_ak8_deltaR_0_3_varbin->Fill(GenHiggs.Pt());
 		genhiggs_eta_with_ak8_deltaR_0_3->Fill(GenHiggs.Eta());
                 }
//////////////////////////////////////////////////////////


        if(LeadAK15.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK15.DeltaR(GenHiggs)<0.5){
                genhiggs_with_ak15_deltaR_0_5->Fill(GenHiggs.Pt());
                }

        if(LeadAK8.Pt()>0&&GenHiggs.Pt()>0&&(isZmm || isZee)&&V_mass[0]>75 && V_mass[0]<105&&V_pt[0]>=60&&LeadAK8.DeltaR(GenHiggs)<0.5){
                genhiggs_with_ak8_deltaR_0_5->Fill(GenHiggs.Pt());
                }




//if (Hcc.Pt()>0&&GenHiggs.Pt()>0){
//genhiggs_with_ak4_no_sel->Fill(GenHiggs.Pt());
//}


ak15_nosel_pt->Fill(LeadAK15.Pt());
ak8_nosel_pt->Fill(LeadAK8.Pt());

//if (Hcc.Pt()>200){
//h_pt->Fill(Hcc.Pt());}











		if (JetPt_1[0]>20 && JetPt_2[0]>20&&GenHiggs.Pt()>0){
			genhiggs_with_ak4_no_sel->Fill(GenHiggs.Pt());
			}


//////////////////////////////////////////////////
//Playing with AK4 to get true efficiency///////
//////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Numerator for resolved topology efficiency////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

            if(Hcc.Pt()>0&&
 (isZmm || isZee)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&JetPt_1[0]>20 && H_mass[0]>50 && H_mass[0]<200 && JetPt_2[0]>20  && V_mass[0]>75 && V_mass[0]<105 
 &&HVdPhi[0]>2.5&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.3){
                        genhiggs_with_ak4_deltaR_0_3->Fill(GenHiggs.Pt());
                        genhiggs_with_ak4_deltaR_0_3_varbin->Fill(GenHiggs.Pt());
			genhiggs_eta_with_ak4_deltaR_0_3->Fill(GenHiggs.Eta());      
			  }




		if( (twoResolvedJets=true)){

	if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
		DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0
		&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)                                  {

		
                        h_pt->Fill(Hcc.Pt());

                      


		if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.1){
			genhiggs_with_ak4->Fill(GenHiggs.Pt());
				}
                 if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.2){
                        genhiggs_with_ak4_deltaR_0_2->Fill(GenHiggs.Pt());
                                }
              //    if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.3){
                //       genhiggs_with_ak4_deltaR_0_3->Fill(GenHiggs.Pt());
                  //              }
  
                if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.5){
                        genhiggs_with_ak4_deltaR_0_5->Fill(GenHiggs.Pt());
                                }

  
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Numerator for Lost of signal search////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
                  
  ////// check fractions of selected AK4s and AK15 and AK8//////
  if (Hcc.Pt()>200&&LeadAK15.Pt()>300){
     frac_ak4_ak15->Fill(LeadAK15.Pt());
 frac_ak4_ak15_ak4->Fill(LeadAK15.Pt()); 
 }

   if (Hcc.Pt()>200&&LeadAK8.Pt()>200){
     frac_ak4_ak8->Fill(LeadAK8.Pt());
frac_ak4_ak8_ak4->Fill(LeadAK8.Pt()); 
 }


                        } // end of cut 

} // end of if resolved 


///////with no Vpt

if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0
&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 &&V_mass[0]>75 && V_mass[0]<105  && HVdPhi[0]>2.5)                                  {


	if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.1){
		genhiggs_with_ak4_no_Vpt->Fill(GenHiggs.Pt());
					}
		
         if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.2){
                        genhiggs_with_ak4_no_Vpt_deltaR_0_2->Fill(GenHiggs.Pt());
                                }
                  if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.3){
                        genhiggs_with_ak4_no_Vpt_deltaR_0_3->Fill(GenHiggs.Pt());
                                }

                  if(Hcc.Pt()>0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.5){
                        genhiggs_with_ak4_no_Vpt_deltaR_0_5->Fill(GenHiggs.Pt());
                                }


       


			}// end of cut no Vpt 
		









////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK8Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


		for(int k=0;k<nLeadAK8;k++)
			{
//std::cout<<"k = "<<k<<std::endl;
//std::cout<<"JetPt_1 = "<<k<<JetPt_1[0]<<" JetPt_2 = "<<JetPt_2[0]<<" DeepJet_CvsL_1 = "<<DeepJet_CvsL_1[0]<<" DeepJet_CvsB_1 = "<<DeepJet_CvsB_1[0]<<" DeepJet_CvsL_2 = "<<DeepJet_CvsL_2[0]<<" DeepJet_CvsB_2"<<DeepJet_CvsB_2[0]<<" H_mass = "<<H_mass[0]<<"V_mass = "<<V_mass[0]<<"V_pt = "<<V_pt[0]<<"HVdPhi = "<<HVdPhi[0]<<std::endl;
if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&

DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0 

&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)                                  {     

		   if (LeadAK8_msoftdrop[k]>50)              {
			 ak_num->Fill(nLeadAK8);	
			 ak_mass->Fill((LeadAK8).M());
                         ak_pt->Fill((LeadAK8).Pt());
                         ak_eta->Fill(LeadAK8.Eta());
                         ak_phi->Fill(LeadAK8.Phi());
			 ak_msoftdrop->Fill(LeadAK8_msoftdrop[k]);		 
      		         ak8_etaphi->Fill(LeadAK8.Eta(),LeadAK8.Phi());                       
                         higgs_ak8->Fill(Hcc.Pt(),LeadAK8.Pt());                      
                                              

               		                                	         }





							}// end of cut 



                                        } //end of loop    


////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK15Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(int k=0;k<nLeadAK15;k++)
                        {
 
if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&

 DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0

&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5){
 
       if (LeadAK15_msoftdrop[k]>50){
          ak15_pt->Fill((LeadAK15).Pt());
                        ak15_mass->Fill(LeadAK15.M());
                        ak15_eta->Fill(LeadAK15.Eta());
                        ak15_phi->Fill(LeadAK15.Phi());
 ak15_msoftdrop->Fill(LeadAK15_msoftdrop[k]);//}
                ak15_etaphi->Fill(LeadAK15.Eta(),LeadAK15.Phi());                                                                                                  														}



	}// end of cut

}// end of for loop
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK4Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
}








////////////////////Plotting graphs//////////////////////////////////////////////
 // Layout
/*
h->GetYaxis()->SetLabelOffset(0.01);
        h->GetYaxis()->SetTitleOffset(1.15);
        h->GetXaxis()->SetTitleOffset(1.15);
gPad->SetTopMargin(0.07);
        gPad->SetRightMargin(0.1);
        gPad->SetLeftMargin(0.15);
        gPad->SetBottomMargin(0.1);

gStyle->SetPadBorderMode(0);
 h2->SetMarkerStyle(37);
        h2->SetMarkerColor(kOrange);
        h2->SetMarkerSize(2.0);
 h1->SetStats(0);
h1->GetXaxis()->SetRangeUser(1,600);


gStyle->SetOptFit(111111);
*/







/*
 
 TCanvas *canvas20 = new TCanvas("canvas732567", "", 1900, 1100);
                canvas20->cd();
               true_higgs_ak15->GetYaxis()->SetTitle("GenHiggs P_{T}, GeV");
               true_higgs_ak15 ->GetYaxis()->SetTitleSize(0.05);
               true_higgs_ak15 ->GetYaxis()->SetTitleOffset(0.80);
               true_higgs_ak15 ->GetXaxis()->SetTitle("AK15 P_{T}, GeV");
               true_higgs_ak15 ->GetXaxis()->SetTitleSize(0.05);
               true_higgs_ak15 ->GetXaxis()->SetTitleOffset(0.90);
               true_higgs_ak15 ->SetLineWidth(2);
               true_higgs_ak15 ->SetLineColor(kOrange+7);
               true_higgs_ak15 ->SetFillColor(kOrange);
                gPad->SetGrid();
               true_higgs_ak15 ->SetTitle("True Higgs vs. AK15");
               true_higgs_ak15 ->Draw("colz");
            
         canvas20->SaveAs("true_higgs_ak15.pdf");
                delete canvas20;

TCanvas *canvas21 = new TCanvas("canvas732567", "", 1900, 1100);
                canvas21->cd();
               true_higgs_ak4->GetYaxis()->SetTitle("GenHiggs P_{T}, GeV");
               true_higgs_ak4 ->GetYaxis()->SetTitleSize(0.05);
               true_higgs_ak4 ->GetYaxis()->SetTitleOffset(0.80);
               true_higgs_ak4 ->GetXaxis()->SetTitle("H P_{T}, GeV");
               true_higgs_ak4 ->GetXaxis()->SetTitleSize(0.05);
               true_higgs_ak4 ->GetXaxis()->SetTitleOffset(0.90);
               true_higgs_ak4 ->SetLineWidth(2);
               true_higgs_ak4 ->SetLineColor(kOrange+7);
               true_higgs_ak4 ->SetFillColor(kOrange);
                gPad->SetGrid();
               true_higgs_ak4 ->SetTitle("True Higgs vs. H_pt");
               true_higgs_ak4 ->Draw("colz");

         canvas20->SaveAs("true_higgs_ak4.pdf");
                delete canvas20;

 TCanvas *canvas14 = new TCanvas("canvas732567", "", 1900, 1100);
                canvas14->cd();
                ak15_pt->GetYaxis()->SetTitle("NEvents");
                ak15_pt->GetYaxis()->SetTitleSize(0.05);
                ak15_pt->GetYaxis()->SetTitleOffset(0.80);
                ak15_pt->GetXaxis()->SetTitle("P_{t},GeV");
                ak15_pt->GetXaxis()->SetTitleSize(0.05);
                ak15_pt->GetXaxis()->SetTitleOffset(0.90);
                ak15_pt->SetLineWidth(2);
                ak15_pt->SetLineColor(kOrange+7);
                gPad->SetGrid();
                ak15_pt->SetTitle("Transverse Momentum");
ak15_pt->SetStats(0);
ak_pt->SetStats(0);
h_pt->SetStats(0);
ak15_pt->Sumw2();
ak15_pt->Divide(genhiggspt);
      ak15_pt->Draw("E");
      h_pt->Draw("same");
      ak_pt->Draw("same");

TLegend* leg66 = new TLegend(0.1,0.88,0.35,0.85);
        leg66->AddEntry( ak15_pt,"AK15Jet");
 leg66->AddEntry( ak_pt,"AK8Jet");
 leg66->AddEntry( h_pt,"resolved");
        leg66->SetTextSize(0.03);
        leg66->SetTextFont(42);
        leg66->SetBorderSize(0);
        leg66->Draw("same");
   canvas14->SaveAs("ak15_dividegenhiggs_pt.pdf");

                delete canvas14;





 TCanvas *c2 = new TCanvas("c2", "", 1900, 1100);
  TGraphAsymmErrors* h_eff_ak15 = new TGraphAsymmErrors(genhiggs_with_ak15_deltaR_0_2,genhiggsall);
  TGraphAsymmErrors* h_eff_ak8 = new TGraphAsymmErrors(genhiggs_with_ak8_deltaR_0_2,genhiggsall);
  TGraphAsymmErrors* h_eff_ak4 = new TGraphAsymmErrors(genhiggs_with_ak4_deltaR_0_2,genhiggsall);
 TGraphAsymmErrors* h_eff_ak4_no_Vpt_deltaR_0_2=new TGraphAsymmErrors(genhiggs_with_ak4_no_Vpt_deltaR_0_2,genhiggsall);
  c2->SetFillStyle(1001);
  c2->SetFillColor(kWhite);
  h_eff_ak15->SetLineColor(kBlack);
  h_eff_ak15->SetLineStyle(1); 
  h_eff_ak15->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_ak15->GetYaxis()->SetTitle("#varepsilon");
  h_eff_ak8->SetLineColor(kRed);
  h_eff_ak8->SetLineStyle(1);
  h_eff_ak4->SetLineColor(kBlue);
  h_eff_ak4->SetLineStyle(1);

 h_eff_ak4_no_Vpt_deltaR_0_2 ->SetLineColor(kGreen);
 h_eff_ak4_no_Vpt_deltaR_0_2->SetLineStyle(1);

  gStyle->SetOptStat(0); 
  h_eff_ak15->UseCurrentStyle();
  h_eff_ak15->Draw("AP");
  h_eff_ak8->Draw("same ][");
  h_eff_ak4->Draw("same ][");
h_eff_ak4_no_Vpt_deltaR_0_2->Draw("same ][");
  TLegend* leg2 = new TLegend(0.1,0.88,0.35,0.75);
        leg2->AddEntry(h_eff_ak15,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.2");
        leg2->AddEntry(h_eff_ak8,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.2");
        leg2->AddEntry(h_eff_ak4,"resolved,  #Delta R(H_{true},H_{reco})<0.2");      
       leg2->AddEntry(h_eff_ak4_no_Vpt_deltaR_0_2,"resolved no Vpt,  #Delta R(H_{true},H_{reco})<0.2");      
        leg2->SetTextSize(0.03);
        leg2->SetTextFont(42);
        leg2->SetBorderSize(0);
        leg2->Draw("same"); 
  c2->SaveAs("tgraphefficiency_ak15_deltaR_0_2.pdf");
  delete c2;
*/

TCanvas *c4 = new TCanvas("c4", "", 1900, 1100);
//TGraphAsymmErrors* ak4_ak15 = new TGraphAsymmErrors(frac_ak4_ak15,ak15_pt);
  c4->SetFillStyle(1001);
  c4->SetFillColor(kWhite);
  frac_ak4_ak15->GetXaxis()->SetTitle("p_{T} [GeV]");
  frac_ak4_ak15->GetYaxis()->SetTitle("fraction");
  frac_ak4_ak15->SetLineColor(kRed);
  frac_ak4_ak15->SetLineStyle(1);
  frac_ak4_ak15->SetMarkerStyle(26);
  frac_ak4_ak15->SetFillStyle(kRed);
  frac_ak4_ak15->SetMarkerSize(1.5);
  gStyle->SetOptStat(0);
/*
int num_int_frac_ak4_ak15=frac_ak4_ak15->Integral();
cout<<"num_int_frac_ak4_ak15  ="<<num_int_frac_ak4_ak15<<endl;

int num_int_3_h_pt=h_pt->Integral();
cout<<"num_int_3_h_pt  ="<<num_int_3_h_pt<<endl;

double  frac_ak4_ak15_h_pt;

frac_ak4_ak15_h_pt=num_int_frac_ak4_ak15/num_int_3_h_pt;
cout<<"Lost of signal in frac_ak4_ak15_h_pt = "<<frac_ak4_ak15_h_pt<<endl;
*/
frac_ak4_ak15->Sumw2();
h_pt->Sumw2();
frac_ak4_ak15->Divide(h_pt);
frac_ak4_ak15->GetYaxis()->SetRangeUser(0., 1.);
  frac_ak4_ak15->Draw("E1");

  TLegend* leg4 = new TLegend(0.1,0.91,0.3,0.85);
   leg4->AddEntry(frac_ak4_ak15,"AK4 & AK15");
        leg4->SetTextSize(0.03);
        leg4->SetTextFont(42);
        leg4->SetBorderSize(0);
        leg4->Draw("same");
  c4->SaveAs("tgraphefficiency_ak4_ak15.pdf");
  delete c4;
   
TCanvas *c61 = new TCanvas("c61","",1900,1100);
//TGraphAsymmErrors* ak4_ak8 = new TGraphAsymmErrors(frac_ak4_ak8,ak_pt);
  c61->SetFillStyle(1001);
  c61->SetFillColor(kWhite);
  frac_ak4_ak8->GetXaxis()->SetTitle("p_{T} [GeV]");
  frac_ak4_ak8->GetYaxis()->SetTitle("Fraction of Events");
  frac_ak4_ak8->SetLineColor(kRed);
  frac_ak4_ak8->SetLineStyle(1);
 frac_ak4_ak8->SetMarkerStyle(26);
  frac_ak4_ak8->SetFillStyle(kRed);
  frac_ak4_ak8->SetMarkerSize(1.5);
  gStyle->SetOptStat(0);
frac_ak4_ak8->Sumw2();
//h_pt->Sumw2();
frac_ak4_ak8->Divide(h_pt);
frac_ak4_ak8->GetYaxis()->SetRangeUser(0., 1.);
  frac_ak4_ak8->Draw("E1");
int num_int_2=frac_ak4_ak8->Integral();
cout<<"ak4_ak4_ak8 ="<<num_int_2<<endl;
  TLegend* leg61 = new TLegend(0.1,0.91,0.3,0.85);
   leg61->AddEntry(frac_ak4_ak8,"AK4 & AK8");
        leg61->SetTextSize(0.03);
        leg61->SetTextFont(42);
        leg61->SetBorderSize(0);
        leg61->Draw("same");
  c61->SaveAs("tgraphefficiency_ak4_ak8.pdf");
  delete c61;


TCanvas *c65 = new TCanvas("c65", "", 1900, 1100);
//TGraphAsymmErrors* ak4_ak15_nosel = new TGraphAsymmErrors(frac_ak4_ak15,ak15_nosel_pt);
  c65->SetFillStyle(1001);
  c65->SetFillColor(kWhite);
//TH1F*ak4_ak15_nosel=new TH1F("ak4_ak15_nosel","",20,0,1000);

  frac_ak4_ak15_ak4->GetXaxis()->SetTitle("p_{T} [GeV]");
  frac_ak4_ak15_ak4->GetYaxis()->SetTitle("Fraction of Events");
  frac_ak4_ak15_ak4->SetLineColor(kRed);
  frac_ak4_ak15_ak4->SetLineStyle(1);
  frac_ak4_ak15_ak4->SetMarkerStyle(26);
  frac_ak4_ak15_ak4->SetFillStyle(kRed);
  frac_ak4_ak15_ak4->SetMarkerSize(1.5);
frac_ak4_ak15_ak4->Sumw2();
ak15_nosel_pt->Sumw2();
frac_ak4_ak15_ak4->Divide(ak15_nosel_pt);
frac_ak4_ak15_ak4->GetYaxis()->SetRangeUser(0., 1.); 
 gStyle->SetOptStat(0);

  frac_ak4_ak15_ak4->Draw("E1");
int num_int=frac_ak4_ak15_ak4->Integral();
cout<<"ak4_ak4_ak15 ="<<num_int<<endl;

  TLegend* leg65 = new TLegend(0.1,0.91,0.3,0.85);
   leg65->AddEntry(frac_ak4_ak15_ak4,"AK4 & AK15");
        leg65->SetTextSize(0.03);
        leg65->SetTextFont(42);
        leg65->SetBorderSize(0);
        leg65->Draw("same");
  c65->SaveAs("tgraphefficiency_ak4_ak15_vs_nosel_ak15_pt.pdf");
  delete c65;

  TCanvas *c66 = new TCanvas("c66","",1900,1100);
//TGraphAsymmErrors* ak4_ak8_nosel = new TGraphAsymmErrors(frac_ak4_ak8,ak8_nosel_pt);
  c66->SetFillStyle(1001);
  c66->SetFillColor(kWhite);
  frac_ak4_ak8_ak4->GetXaxis()->SetTitle("p_{T} [GeV]");
  frac_ak4_ak8_ak4->GetYaxis()->SetTitle("Fraction of Events");
  frac_ak4_ak8_ak4->SetLineColor(kRed);
  frac_ak4_ak8_ak4->SetLineStyle(1);
  frac_ak4_ak8_ak4->SetMarkerStyle(26);
  frac_ak4_ak8_ak4->SetFillStyle(kRed);
  frac_ak4_ak8_ak4->SetMarkerSize(1.5);
frac_ak4_ak8_ak4->GetYaxis()->SetRangeUser(0., 1.);
frac_ak4_ak8_ak4->Sumw2();
ak8_nosel_pt->Sumw2();
frac_ak4_ak8_ak4->Divide(ak8_nosel_pt);
  gStyle->SetOptStat(0);
  frac_ak4_ak8_ak4->Draw("E1");
int num_int_1=frac_ak4_ak8_ak4->Integral();
cout<<"ak4_ak4_ak8_nosel_pt ="<<num_int_1<<endl;
TLegend* leg66 = new TLegend(0.1,0.91,0.3,0.85);
   leg66->AddEntry(frac_ak4_ak8_ak4,"AK4 & AK8");
        leg66->SetTextSize(0.03);
        leg66->SetTextFont(42);
        leg66->SetBorderSize(0);
        leg66->Draw("same");
  c66->SaveAs("tgraphefficiency_ak4_ak8_vs_nosel_ak8_pt.pdf");
  delete c66;


 TCanvas *c6 = new TCanvas("c6", "", 1900, 1100);
  TGraphAsymmErrors* h_eff_ak15_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_with_ak15_deltaR_0_3,genhiggsall);
  TGraphAsymmErrors* h_eff_ak8_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_with_ak8_deltaR_0_3,genhiggsall);
  TGraphAsymmErrors* h_eff_ak4_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_with_ak4_deltaR_0_3,genhiggsall_jet);


  c6->SetFillStyle(1001);
  c6->SetFillColor(kWhite);
  h_eff_ak15_deltaR_0_3->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_ak15_deltaR_0_3->GetYaxis()->SetTitle("fraction");
  h_eff_ak15_deltaR_0_3->SetLineColor(6);
  h_eff_ak15_deltaR_0_3->SetLineStyle(1);
  h_eff_ak8_deltaR_0_3->SetLineColor(kRed);
  h_eff_ak8_deltaR_0_3->SetLineStyle(1);

  h_eff_ak4_deltaR_0_3->SetLineColor(kBlue);
  h_eff_ak4_deltaR_0_3->SetLineStyle(1);


  gStyle->SetOptStat(0);
  h_eff_ak15_deltaR_0_3->Draw("AP");
  h_eff_ak8_deltaR_0_3->Draw("same ][");
  h_eff_ak4_deltaR_0_3->Draw("same ][");
   TLegend* leg6 = new TLegend(0.1,0.88,0.35,0.75);
        leg6->AddEntry(h_eff_ak15_deltaR_0_3,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.3");
        leg6->AddEntry(h_eff_ak8_deltaR_0_3,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.3");
        leg6->AddEntry(h_eff_ak4_deltaR_0_3,"resolved,  #Delta R(H_{true},H_{reco})<0.3");


        leg6->SetTextSize(0.03);
        leg6->SetTextFont(42);
        leg6->SetBorderSize(0);
        leg6->Draw("same");
  c6->SaveAs("tgraphefficiency_ak15_deltaR_0_3.pdf");
  delete c6;

/*
 TCanvas *c63 = new TCanvas("c6", "", 1900, 1100);
  TGraphAsymmErrors* h_eff_ak4_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_with_ak4_deltaR_0_3,genhiggsall_jet);
  c63->SetFillStyle(1001);
  c63->SetFillColor(kWhite);
  h_eff_ak4_deltaR_0_3->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_ak4_deltaR_0_3->GetYaxis()->SetTitle("fraction");
  h_eff_ak4_deltaR_0_3->SetLineColor(kBlue);
  h_eff_ak4_deltaR_0_3->SetLineStyle(1);

  gStyle->SetOptStat(0);
  h_eff_ak4_deltaR_0_3->Draw("AP");
  TLegend* leg63 = new TLegend(0.1,0.88,0.35,0.75);
        leg63->AddEntry(h_eff_ak4_deltaR_0_3,"resolved,  #Delta R(H_{true},H_{reco})<0.3");
        leg63->SetTextSize(0.03);
        leg63->SetTextFont(42);
        leg63->SetBorderSize(0);
        leg63->Draw("same");
  c63->SaveAs("tgraphefficiency_ak4_deltaR_0_3.pdf");
  delete c63;
*/





 TCanvas *c60 = new TCanvas("c60", "", 1900, 1100);
 
 TMultiGraph  *mg  = new TMultiGraph();
   
   
  TGraphAsymmErrors* h_eff_eta_ak15_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_eta_with_ak15_deltaR_0_3,genhiggs_eta_all);
  TGraphAsymmErrors* h_eff_eta_ak8_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_eta_with_ak8_deltaR_0_3,genhiggs_eta_all);
  TGraphAsymmErrors* h_eff_eta_ak4_deltaR_0_3 = new TGraphAsymmErrors(genhiggs_eta_with_ak4_deltaR_0_3,genhiggs_eta_all_jet);
  c60->SetFillStyle(1001);
  c60->SetFillColor(kWhite);
  h_eff_eta_ak15_deltaR_0_3 ->GetXaxis()->SetTitle("#eta");
  h_eff_eta_ak15_deltaR_0_3  ->GetYaxis()->SetTitle("fraction");
  h_eff_eta_ak15_deltaR_0_3->SetLineColor(6);
  h_eff_eta_ak15_deltaR_0_3->SetLineStyle(1);
  h_eff_eta_ak8_deltaR_0_3->SetLineColor(kRed);
  h_eff_eta_ak8_deltaR_0_3->SetLineStyle(1);
  h_eff_eta_ak4_deltaR_0_3->SetLineColor(kBlue);
  h_eff_eta_ak4_deltaR_0_3->SetLineStyle(1);
  gStyle->SetOptStat(0);
  mg->Add(h_eff_eta_ak15_deltaR_0_3);
  mg->Add(h_eff_eta_ak8_deltaR_0_3); 
  mg->Add(h_eff_eta_ak4_deltaR_0_3);
  mg->Draw("AP");

 // h_eff_eta_ak15_deltaR_0_3->Draw("AP");
 //  h_eff_eta_ak4_deltaR_0_3->Draw("same ][");
 // h_eff_eta_ak8_deltaR_0_3->Draw("same ][");
 
  TLegend* leg60 = new TLegend(0.1,0.61,0.30,0.50);
        leg60->AddEntry(h_eff_eta_ak15_deltaR_0_3,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.3");
        leg60->AddEntry(h_eff_eta_ak8_deltaR_0_3,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.3");
        leg60->AddEntry(h_eff_eta_ak4_deltaR_0_3,"resolved,  #Delta R(H_{true},H_{reco})<0.3");
        leg60->SetTextSize(0.03);
        leg60->SetTextFont(42);
        leg60->SetBorderSize(0);
        leg60->Draw("same");
  c60->SaveAs("tgraphefficiency_eta_ak15_deltaR_0_3.pdf");
  delete c60;

 TCanvas *c64 = new TCanvas("c64", "", 1900, 1100);

 TMultiGraph  *mg1  = new TMultiGraph();
  TGraphAsymmErrors* h_eff_ak15_deltaR_0_3_varbin = new TGraphAsymmErrors(genhiggs_with_ak15_deltaR_0_3_varbin,genhiggsall_varbin);
  TGraphAsymmErrors* h_eff_ak8_deltaR_0_3_varbin = new TGraphAsymmErrors(genhiggs_with_ak8_deltaR_0_3_varbin,genhiggsall_varbin);
  TGraphAsymmErrors* h_eff_ak4_deltaR_0_3_varbin = new TGraphAsymmErrors(genhiggs_with_ak4_deltaR_0_3_varbin,genhiggsall_jet_varbin);
  mg1->Add(h_eff_ak15_deltaR_0_3_varbin);

  mg1->Add(h_eff_ak8_deltaR_0_3_varbin);
  mg1->Add(h_eff_ak4_deltaR_0_3_varbin);
  mg1->GetXaxis()->SetTitle("p_{T} [GeV]");
  mg1->GetYaxis()->SetTitle("Efficiency");
mg1->GetYaxis()->SetRangeUser(0., 1.);  
c64->SetFillStyle(1001);
  c64->SetFillColor(kWhite);
//  h_eff_ak15_deltaR_0_3_varbin->GetXaxis()->SetTitle("p_{T} [GeV]");
//  h_eff_ak15_deltaR_0_3_varbin->GetYaxis()->SetTitle("fraction");
  h_eff_ak15_deltaR_0_3_varbin->SetLineColor(kGreen);
  h_eff_ak15_deltaR_0_3_varbin->SetLineStyle(1);
  h_eff_ak8_deltaR_0_3_varbin->SetLineColor(kRed);
  h_eff_ak8_deltaR_0_3_varbin->SetLineStyle(1);

  h_eff_ak4_deltaR_0_3_varbin->SetLineColor(kBlue);
  h_eff_ak4_deltaR_0_3_varbin->SetLineStyle(1);

h_eff_ak15_deltaR_0_3_varbin->SetMarkerStyle(25);
h_eff_ak8_deltaR_0_3_varbin->SetMarkerStyle(24);
h_eff_ak4_deltaR_0_3_varbin->SetMarkerStyle(26);

h_eff_ak15_deltaR_0_3_varbin->SetFillStyle(kGreen);
h_eff_ak8_deltaR_0_3_varbin->SetFillStyle(kRed);
h_eff_ak4_deltaR_0_3_varbin->SetFillStyle(kBlue);

 h_eff_ak15_deltaR_0_3_varbin->SetMarkerSize(1.5);
h_eff_ak8_deltaR_0_3_varbin->SetMarkerSize(1.5);
h_eff_ak4_deltaR_0_3_varbin->SetMarkerSize(1.5);
  gStyle->SetOptStat(0);
h_eff_ak4_deltaR_0_3_varbin->SetFillColor(16);
h_eff_ak4_deltaR_0_3_varbin->SetFillStyle(1001);
 
mg1->Draw("a5");
mg1->Draw("p");
   // h_eff_ak15_deltaR_0_3_varbin->Draw("a2");
   // h_eff_ak15_deltaR_0_3_varbin->Draw("p");

//  h_eff_ak15_deltaR_0_3_varbin->Draw("APS; ; 5 S=0.5");
//  h_eff_ak8_deltaR_0_3_varbin->Draw("same ][");
//  h_eff_ak4_deltaR_0_3_varbin->Draw("same ][");
//  h_eff_ak8_deltaR_0_3_varbin->Draw("same a2");
 // h_eff_ak4_deltaR_0_3_varbin->Draw("same a2");
   TLegend* leg64 = new TLegend(0.1,0.88,0.35,0.75);
        leg64->AddEntry(h_eff_ak15_deltaR_0_3_varbin,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.3");
        leg64->AddEntry(h_eff_ak8_deltaR_0_3_varbin,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.3");
        leg64->AddEntry(h_eff_ak4_deltaR_0_3_varbin,"resolved,  #Delta R(H_{true},H_{reco})<0.3");


        leg64->SetTextSize(0.03);
        leg64->SetTextFont(42);
        leg64->SetBorderSize(0);
        leg64->Draw("same");
  c64->SaveAs("tgraphefficiency_ak15_deltaR_0_3_varbin.pdf");
  delete c64;

/*
 TCanvas *c7 = new TCanvas("c7", "", 1900, 1100);
  TGraphAsymmErrors* h_eff_ak15_deltaR_0_1 = new TGraphAsymmErrors(genhiggs_with_ak15,genhiggsall);
  TGraphAsymmErrors* h_eff_ak8_deltaR_0_1 = new TGraphAsymmErrors(genhiggs_with_ak8,genhiggsall);
  TGraphAsymmErrors* h_eff_ak4_deltaR_0_1 = new TGraphAsymmErrors(genhiggs_with_ak4,genhiggsall);
 TGraphAsymmErrors* h_eff_ak4_no_Vpt_deltaR_0_1=new TGraphAsymmErrors(genhiggs_with_ak4_no_Vpt,genhiggsall); 
  c7->SetFillStyle(1001);
  c7->SetFillColor(kWhite);
  h_eff_ak15_deltaR_0_1->SetLineColor(kBlack);
  h_eff_ak15_deltaR_0_1->SetLineStyle(1);
  h_eff_ak15_deltaR_0_1->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_ak15_deltaR_0_1->GetYaxis()->SetTitle("#varepsilon");
  h_eff_ak8_deltaR_0_1->SetLineColor(kRed);
  h_eff_ak8_deltaR_0_1->SetLineStyle(1);
  h_eff_ak4_deltaR_0_1->SetLineColor(kBlue);
  h_eff_ak4_deltaR_0_1->SetLineStyle(1);
h_eff_ak4_no_Vpt_deltaR_0_1->SetLineColor(kGreen);
h_eff_ak4_no_Vpt_deltaR_0_1->SetLineStyle(1);
  gStyle->SetOptStat(0);
  h_eff_ak15_deltaR_0_1->UseCurrentStyle();
  h_eff_ak15_deltaR_0_1->Draw("AP");
  h_eff_ak8_deltaR_0_1->Draw("same ][");
  h_eff_ak4_deltaR_0_1->Draw("same ][");
h_eff_ak4_no_Vpt_deltaR_0_1->Draw("same ][");
  TLegend* leg7 = new TLegend(0.1,0.88,0.35,0.75);
        leg7->AddEntry(h_eff_ak15_deltaR_0_1,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.1");
        leg7->AddEntry(h_eff_ak8_deltaR_0_1,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.1");
        leg7->AddEntry(h_eff_ak4_deltaR_0_1,"resolved,  #Delta R(H_{true},H_{reco})<0.1");
       leg7->AddEntry(h_eff_ak4_no_Vpt_deltaR_0_1,"resolved no Vpt,  #Delta R(H_{true},H_{reco})<0.1");
        leg7->SetTextSize(0.03);
        leg7->SetTextFont(42);
        leg7->SetBorderSize(0);
        leg7->Draw("same");
  c7->SaveAs("tgraphefficiency_ak15_deltaR_0_1.pdf");
  delete c7;


 TCanvas *c8 = new TCanvas("c8", "", 1900, 1100);
  TGraphAsymmErrors* h_eff_ak15_deltaR_0_5 = new TGraphAsymmErrors(genhiggs_with_ak15_deltaR_0_5,genhiggsall);
  TGraphAsymmErrors* h_eff_ak8_deltaR_0_5 = new TGraphAsymmErrors(genhiggs_with_ak8_deltaR_0_5,genhiggsall);
  TGraphAsymmErrors* h_eff_ak4_deltaR_0_5 = new TGraphAsymmErrors(genhiggs_with_ak4_deltaR_0_5,genhiggsall);
 TGraphAsymmErrors* h_eff_ak4_no_Vpt_deltaR_0_5=new TGraphAsymmErrors(genhiggs_with_ak4_no_Vpt_deltaR_0_5,genhiggsall);
 c8->SetFillStyle(1001);
  c8->SetFillColor(kWhite);
  h_eff_ak15_deltaR_0_5->SetLineColor(kBlack);
  h_eff_ak15_deltaR_0_5->SetLineStyle(1);
  h_eff_ak15_deltaR_0_5->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_eff_ak15_deltaR_0_5->GetYaxis()->SetTitle("#varepsilon");
  h_eff_ak8_deltaR_0_5->SetLineColor(kRed);
  h_eff_ak8_deltaR_0_5->SetLineStyle(1);
  h_eff_ak4_deltaR_0_5->SetLineColor(kBlue);
  h_eff_ak4_deltaR_0_5->SetLineStyle(1);
 h_eff_ak4_no_Vpt_deltaR_0_5->SetLineColor(kGreen);
  h_eff_ak4_no_Vpt_deltaR_0_5->SetLineStyle(1);

  gStyle->SetOptStat(0);
  h_eff_ak15_deltaR_0_5->UseCurrentStyle();
  h_eff_ak15_deltaR_0_5->Draw("AP");
  h_eff_ak8_deltaR_0_5->Draw("same ][");
  h_eff_ak4_deltaR_0_5->Draw("same ][");
h_eff_ak4_no_Vpt_deltaR_0_5->Draw("same ][");
  TLegend* leg8 = new TLegend(0.1,0.88,0.35,0.75);
        leg8->AddEntry(h_eff_ak15_deltaR_0_5,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.5");
        leg8->AddEntry(h_eff_ak8_deltaR_0_5,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.5");
        leg8->AddEntry(h_eff_ak4_deltaR_0_5,"resolved,  #Delta R(H_{true},H_{reco})<0.5");
        leg8->AddEntry(h_eff_ak4_no_Vpt_deltaR_0_5,"resolved no Vpt,  #Delta R(H_{true},H_{reco})<0.5");
        leg8->SetTextSize(0.03);
        leg8->SetTextFont(42);
        leg8->SetBorderSize(0);
        leg8->Draw("same");
  c8->SaveAs("tgraphefficiency_ak15_deltaR_0_5.pdf");
  delete c8;
*/

/*
 TCanvas *canvas50 = new TCanvas("canvas50732567", "", 1900, 1100);
                canvas50->cd();
 canvas50->SetFillStyle(1001);
   canvas50->SetFillColor(kWhite);
                genhiggs_with_ak15->GetYaxis()->SetTitle("efficiency");
                genhiggs_with_ak15->GetYaxis()->SetTitleSize(0.05);
                genhiggs_with_ak15->GetYaxis()->SetTitleOffset(0.80);
                genhiggs_with_ak15->GetXaxis()->SetTitle("p_{T}, GeV");
                genhiggs_with_ak15->GetXaxis()->SetTitleSize(0.05);
                genhiggs_with_ak15->GetXaxis()->SetTitleOffset(0.90);
                genhiggs_with_ak15->SetLineWidth(2);
genhiggs_with_ak15->GetYaxis()->SetRangeUser(0., 1.);
               genhiggs_with_ak15->SetLineColor(kBlack);
                genhiggs_with_ak8->SetLineColor(kRed);
		genhiggs_with_ak4->SetLineColor(kBlue);
genhiggs_with_ak4_no_Vpt->SetLineColor(kOrange);
genhiggs_with_ak15->SetStats(0);
genhiggs_with_ak4->SetStats(0);
genhiggs_with_ak4_no_Vpt->SetStats(0);
genhiggsall->SetStats(0);
genhiggs_with_ak15->Sumw2();
genhiggs_with_ak4->Sumw2();
genhiggs_with_ak4_no_Vpt->Sumw2();
genhiggsall->Sumw2();
genhiggs_with_ak8->Sumw2();
genhiggs_with_ak8->Divide(genhiggsall);
genhiggs_with_ak15->Divide(genhiggsall);
genhiggs_with_ak4->Divide(genhiggsall);
genhiggs_with_ak4_no_Vpt->Divide(genhiggsall);
      genhiggs_with_ak15->Draw("E");
genhiggs_with_ak8->Draw("same");
      genhiggs_with_ak4->Draw("same");
genhiggs_with_ak4_no_Vpt->Draw("same");
TLegend* leg666 = new TLegend(0.1,0.78,0.42,0.95);
        leg666->AddEntry(genhiggs_with_ak15,"with AK15Jet,  #deltaR(H_{True},H_{AK15Jets})<0.1");
 leg666->AddEntry(genhiggs_with_ak8,"with AK8Jet,  #deltaR(H_{True},H_{AK8Jets})<0.1");
 leg666->AddEntry(genhiggs_with_ak4,"with resolved,  #deltaR(H_{True},H_{AK4Jets})<0.1");
 leg666->AddEntry(genhiggs_with_ak4_no_Vpt,"with resolved (no Vpt),  #deltaR(H_{True},H_{AK4Jets})<0.1");       
 leg666->SetTextSize(0.03);
        leg666->SetTextFont(42);
        leg666->SetBorderSize(0);
        leg666->Draw("same");
   canvas50->SaveAs("efficiency_genhiggs.pdf");
                delete canvas50;





 
 TCanvas *canvas67 = new TCanvas("canvas50732567", "", 1900, 1100);
                canvas67->cd();
 canvas67->SetFillStyle(1001);
   canvas67->SetFillColor(kWhite);
                genhiggs_with_ak4_no_sel->GetYaxis()->SetTitle("efficiency");
                genhiggs_with_ak4_no_sel->GetYaxis()->SetTitleSize(0.05);
                genhiggs_with_ak4_no_sel->GetYaxis()->SetTitleOffset(0.80);
                genhiggs_with_ak4_no_sel->GetXaxis()->SetTitle("p_{T}, GeV");
                genhiggs_with_ak4_no_sel->GetXaxis()->SetTitleSize(0.05);
                genhiggs_with_ak4_no_sel->GetXaxis()->SetTitleOffset(0.90);
                genhiggs_with_ak4_no_sel->SetLineWidth(2);
genhiggs_with_ak4_no_sel->GetYaxis()->SetRangeUser(0., 1.);
                genhiggs_with_ak4_no_sel->SetLineColor(kRed);
genhiggs_with_ak4_no_sel->SetStats(0);
genhiggs_with_ak4_no_sel->Sumw2();
genhiggs_with_ak4_no_sel->Divide(genhiggsall);
      genhiggs_with_ak4_no_sel->Draw("E");
TLegend* leg67 = new TLegend(0.1,0.78,0.42,0.65);
 leg67->AddEntry(genhiggs_with_ak4_no_sel,"with resolved (no Higgs selection)");
 leg67->SetTextSize(0.03);
        leg67->SetTextFont(42);
        leg67->SetBorderSize(0);
        leg67->Draw("same");
   canvas67->SaveAs("efficiency_genhiggs_ak4_no_sel.pdf");
                delete canvas67;

*/

//ak_msoftdrop->Fit(func,"R");
//ak15_msoftdrop->Fit(func_ak15,"R");

frac_ak4_ak15->Write();
 genhiggs_eta_with_ak4_deltaR_0_3->Write();
	ak15_mass->Write();
        ak15_pt->Write();
        ak15_eta->Write();
        ak15_phi->Write();
        ak15_msoftdrop->Write();


	ak_eta->Write();
	ak_phi->Write();
	ak_pt->Write();
        ak_mass->Write();
        ak_num->Write(); 
        ak_msoftdrop->Write();
 
genhiggsall->Write();
genhiggs_with_ak15->Write();
genhiggs_with_ak4->Write(); 
genhiggs_with_ak8->Write();


 genhiggs_with_ak4_no_Vpt->Write();
 genhiggs_with_ak4_no_sel->Write();


genhiggs_with_ak15_deltaR_0_2->Write();
genhiggs_with_ak15_deltaR_0_3->Write();
genhiggs_with_ak8_deltaR_0_2->Write();
genhiggs_with_ak8_deltaR_0_3->Write();

 genhiggs_with_ak4_deltaR_0_2->Write();
 genhiggs_with_ak4_deltaR_0_3->Write();

 genhiggs_with_ak4_no_Vpt_deltaR_0_2->Write();
 genhiggs_with_ak4_no_Vpt_deltaR_0_3->Write();


fnew->Close();

}






