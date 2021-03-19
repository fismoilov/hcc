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


void akjet_efficiency_ZH125ToCC_ZNuNu_0L() {


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


Int_t isZnn;
Int_t nJetsCloseToMET;
Float_t dPhi_MET_TkMET[50]; 
Float_t  MET_Pt[50];


	TChain *t = new TChain("Events");
 // TChain *chain = new TChain("tree");
   for (int fileNum=1;fileNum < 39;fileNum++) {
      t->AddFile(Form(" output_ZH125ToCC_ZNuNu_powheg_%d.root", fileNum));
//std::cout<<"File number = "<<fileNum<<std::endl; 
 }
  t->SetBranchAddress("twoResolvedJets",&twoResolvedJets);
  
  t->SetBranchAddress("nJetsCloseToMET",&nJetsCloseToMET);
        t->SetBranchAddress("isZnn",&isZnn);
        t->SetBranchAddress("dPhi_MET_TkMET",dPhi_MET_TkMET);
        t->SetBranchAddress("MET_Pt",MET_Pt);

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





TH1F *genhiggsall= new TH1F("genhiggsall","",20,0,1000);
TH1F *genhiggsall_jet= new TH1F("genhiggsall_jet","",20,0,1000);

TH1F *genhiggs_with_ak15= new TH1F("genhiggs_with_ak15","",20,0,1000);
TH1F *genhiggs_with_ak4= new TH1F("genhiggs_with_ak4","",20,0,1000);
TH1F *genhiggs_with_ak8= new TH1F("genhiggs_with_ak8","",20,0,1000);
TH1F *genhiggs_with_ak4_noVpt= new TH1F("genhiggs_with_ak4_noVpt","",20,0,1000);
TH1F *genhiggs_with_ak4_no_sel= new TH1F("genhiggs_with_ak4_no_sel","",20,0,1000);

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

TH1F *frac_ak4_ak15=new TH1F("frac_ak4_ak15","",20,0,1000);
TH1F *frac_ak4_ak8=new TH1F("frac_ak4_ak8","",20,0,1000);
TH1F *frac_ak4_ak15_ak4=new TH1F("frac_ak4_ak15_ak4","",20,0,1000);
TH1F *frac_ak4_ak8_ak4=new TH1F("frac_ak4_ak8_ak4","",20,0,1000);
TH1F *ak15_nosel_pt=new TH1F("ak15_nosel_pt","",20,0,1000);
TH1F *ak8_nosel_pt=new TH1F("ak8_nosel_pt","",20,0,1000);
 TH1F *h_pt=new TH1F("h_pt","",20,0,1000);

TH1F *genhiggs_with_ak15_deltaR_0_2=new TH1F("genhiggs_with_ak15_deltaR_0_2","",20,0,1000);
TH1F *genhiggs_with_ak15_deltaR_0_3=new TH1F("genhiggs_with_ak15_deltaR_0_3","",20,0,1000);
TH1F *genhiggs_with_ak15_deltaR_0_5=new TH1F("genhiggs_with_ak15_deltaR_0_5","",20,0,1000);

TH1F *genhiggs_with_ak8_deltaR_0_2=new TH1F("genhiggs_with_ak8_deltaR_0_2","",20,0,1000);
TH1F *genhiggs_with_ak8_deltaR_0_3=new TH1F("genhiggs_with_ak8_deltaR_0_3","",20,0,1000);
TH1F *genhiggs_with_ak8_deltaR_0_5=new TH1F("genhiggs_with_ak8_deltaR_0_5","",20,0,1000);


 TH1F *genhiggs_with_ak4_deltaR_0_2=new TH1F("genhiggs_with_ak4_deltaR_0_2","",20,0,1000);
 TH1F *genhiggs_with_ak4_deltaR_0_3=new TH1F("genhiggs_with_ak4_deltaR_0_3","",20,0,1000);
 TH1F *genhiggs_with_ak4_deltaR_0_5=new TH1F("genhiggs_with_ak4_deltaR_0_5","",20,0,1000);






       
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

	int NumberOfEvents =298868;
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



///////////////////////////////////////////////////////////////////
///True Higgs////////// ///// ///// ///// ///// ///// ///// /////// 
///////// ///// ///// ///// ///// ///// ///// ///// ///// ///// /// 

for (int k=0;k<nGenPart;k++){

 if(GenPart_pdgId[k]==25&&GenPart_genPartIdxMother[k]==0){
        GenHiggs.SetPtEtaPhiM(GenPart_pt[k],GenPart_eta[k],GenPart_phi[k],GenPart_mass[k]);

                if(GenHiggs.Pt()>0){
                        genhiggsall->Fill(GenHiggs.Pt());
                }


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Denominator for resolved topology efficiency////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

                if(GenHiggs.Pt()>0&&(isZnn)&&JetPt_1[0]>20 &&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&JetPt_2[0]>20){
                genhiggsall_jet->Fill(GenHiggs.Pt());
                genhiggsall_jet_varbin->Fill(GenHiggs.Pt());

                                  }

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Denominator for boosted  topology efficiency////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 if(GenHiggs.Pt()>0&&(isZnn)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4){

                 genhiggsall->Fill(GenHiggs.Pt());
                 genhiggsall_varbin->Fill(GenHiggs.Pt());

                                     }

                }
        }  // end of loop GenPart

       

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Numerator  for boosted  topology efficiency////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

         if(LeadAK15.Pt()>0&&GenHiggs.Pt()>0&&(isZnn)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&LeadAK15.DeltaR(GenHiggs)<0.3){
                genhiggs_with_ak15_deltaR_0_3->Fill(GenHiggs.Pt());
                genhiggs_with_ak15_deltaR_0_3_varbin->Fill(GenHiggs.Pt());

                 }

        if(LeadAK8.Pt()>0&&GenHiggs.Pt()>0&&(isZnn)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&LeadAK8.DeltaR(GenHiggs)<0.3){
                genhiggs_with_ak8_deltaR_0_3->Fill(GenHiggs.Pt());
                genhiggs_with_ak8_deltaR_0_3_varbin->Fill(GenHiggs.Pt());
                 }




ak15_nosel_pt->Fill(LeadAK15.Pt());
ak8_nosel_pt->Fill(LeadAK8.Pt());
if (Hcc.Pt()>200){
h_pt->Fill(Hcc.Pt());}

 if (JetPt_1[0]>20 && JetPt_2[0]>20&&GenHiggs.Pt()>0){
                        genhiggs_with_ak4_no_sel->Fill(GenHiggs.Pt());
                        }




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////Numerator  for resolved topology efficiency////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


            if(Hcc.Pt()>0&&
 (isZnn)&&GenHiggs.Eta()<2.4&&GenHiggs.Eta()>-2.4 &&JetPt_1[0]>20 && JetPt_2[0]>20 &&H_mass[0]>50 && H_mass[0]<200&&HVdPhi[0]>2.0&&GenHiggs.Pt()>0&&Hcc.DeltaR(GenHiggs)<0.3){
                        genhiggs_with_ak4_deltaR_0_3->Fill(GenHiggs.Pt());
                        genhiggs_with_ak4_deltaR_0_3_varbin->Fill(GenHiggs.Pt());
                          }





		if( (twoResolvedJets=true)){

if(( isZnn) && nJetsCloseToMET==0 && dPhi_MET_TkMET[0] < 0.5 && max(JetPt_1[0],JetPt_2[0]) > 60 && min(JetPt_1[0],JetPt_2[0]) > 35 && MET_Pt[0]>170 && H_pt[0]>120 && controlSample==0 &&
 H_mass[0]>50 && H_mass[0]<200 && DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.00 && DeepJet_CvsB_2[0]>0.00 && HVdPhi[0]>2.0)                                  {
events_pass_ak4_higgs_selection++;



if (Hcc.Pt()>200&&LeadAK15.Pt()>200){
     frac_ak4_ak15->Fill(LeadAK15.Pt());
 frac_ak4_ak15_ak4->Fill(LeadAK15.Pt());
 }

   if (Hcc.Pt()>200&&LeadAK8.Pt()>200){
     frac_ak4_ak8->Fill(LeadAK8.Pt());
frac_ak4_ak8_ak4->Fill(LeadAK8.Pt());
 }






		}
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK8Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////


		for(int k=0;k<nLeadAK8;k++)
			{
//std::cout<<"k = "<<k<<std::endl;
//std::cout<<"JetPt_1 = "<<k<<JetPt_1[0]<<" JetPt_2 = "<<JetPt_2[0]<<" DeepJet_CvsL_1 = "<<DeepJet_CvsL_1[0]<<" DeepJet_CvsB_1 = "<<DeepJet_CvsB_1[0]<<" DeepJet_CvsL_2 = "<<DeepJet_CvsL_2[0]<<" DeepJet_CvsB_2"<<DeepJet_CvsB_2[0]<<" H_mass = "<<H_mass[0]<<"V_mass = "<<V_mass[0]<<"V_pt = "<<V_pt[0]<<"HVdPhi = "<<HVdPhi[0]<<std::endl;
if(( isZnn) && nJetsCloseToMET==0 && dPhi_MET_TkMET[0] < 0.5 && max(JetPt_1[0],JetPt_2[0]) > 60 && min(JetPt_1[0],JetPt_2[0]) > 35 && MET_Pt[0]>170 && H_pt[0]>120 && controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.00 && DeepJet_CvsB_2[0]>0.00 && HVdPhi[0]>2.0)


                                  {     

		   if (LeadAK8_msoftdrop[k]>50&&LeadAK8.Pt()>200)              {
			 ak_num->Fill(nLeadAK8);	
			 ak_mass->Fill((LeadAK8).M());
                         ak_pt->Fill((LeadAK8).Pt());
                         ak_eta->Fill(LeadAK8.Eta());
                         ak_phi->Fill(LeadAK8.Phi());
			 ak_msoftdrop->Fill(LeadAK8_msoftdrop[k]);		 
      		         ak8_etaphi->Fill(LeadAK8.Eta(),LeadAK8.Phi());                       
                         higgs_ak8->Fill(Hcc.Pt(),LeadAK8.Pt());                      
                                              

               		                                         }
    		 if (LeadAK8_msoftdrop[k]>50&&LeadAK8.Pt()>250)             {
              	ak_msoftdrop_pt250->Fill(LeadAK8_msoftdrop[k]);}

		 if (LeadAK8_msoftdrop[k]>50&&LeadAK8.Pt()>300)             {
                ak_msoftdrop_pt300->Fill(LeadAK8_msoftdrop[k]);}
}



                                        } //end of loop    
////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK15Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(int k=0;k<nLeadAK15;k++)
                        {
 if(( isZnn) && nJetsCloseToMET==0 && dPhi_MET_TkMET[0] < 0.5 && max(JetPt_1[0],JetPt_2[0]) > 60 && min(JetPt_1[0],JetPt_2[0]) > 35 && MET_Pt[0]>170 && H_pt[0]>120 && controlSample==0 &&
 H_mass[0]>50 && H_mass[0]<200 && DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.00 && DeepJet_CvsB_2[0]>0.00 && HVdPhi[0]>2.0){
  if (LeadAK15_msoftdrop[k]>50){
                        ak15_mass->Fill(LeadAK15.M());
                        ak15_pt->Fill((LeadAK15).Pt());
                        ak15_eta->Fill(LeadAK15.Eta());
                        ak15_phi->Fill(LeadAK15.Phi());
 ak15_msoftdrop->Fill(LeadAK15_msoftdrop[k]);//}
                ak15_etaphi->Fill(LeadAK15.Eta(),LeadAK15.Phi());                                                                                                  														}
 	        if (LeadAK15_msoftdrop[k]>50&&LeadAK15.Pt()>250)             {
                ak15_msoftdrop_pt250->Fill(LeadAK15_msoftdrop[k]);}

                 if (LeadAK15_msoftdrop[k]>50&&LeadAK15.Pt()>300)             {
                ak15_msoftdrop_pt300->Fill(LeadAK15_msoftdrop[k]);}





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


TCanvas *c4 = new TCanvas("c4", "", 1900, 1100);
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
frac_ak4_ak15->Sumw2();
ak15_pt->Sumw2();
frac_ak4_ak15->Divide(ak15_pt);
frac_ak4_ak15->GetYaxis()->SetRangeUser(0., 1.);
  frac_ak4_ak15->Draw("E1");
int num_int_3=frac_ak4_ak15->Integral();
cout<<"ak4_ak4_ak15 ="<<num_int_3<<endl;

  TLegend* leg4 = new TLegend(0.1,0.91,0.3,0.85);
   leg4->AddEntry(frac_ak4_ak15,"AK4 & AK15");
        leg4->SetTextSize(0.03);
        leg4->SetTextFont(42);
        leg4->SetBorderSize(0);
        leg4->Draw("same");
  c4->SaveAs("tgraphefficiency_ak4_ak15_0L_channel.pdf");
  delete c4;

TCanvas *c61 = new TCanvas("c61","",1900,1100);
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
ak_pt->Sumw2();
frac_ak4_ak8->Divide(ak_pt);
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
  c61->SaveAs("tgraphefficiency_ak4_ak8_0L_channel.pdf");
  delete c61;
TCanvas *c65 = new TCanvas("c65", "", 1900, 1100);
 c65->SetFillStyle(1001);
  c65->SetFillColor(kWhite);
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
  c65->SaveAs("tgraphefficiency_ak4_ak15_vs_nosel_ak15_pt_0L_channel.pdf");
  delete c65;

 TCanvas *c66 = new TCanvas("c66","",1900,1100);
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
  c66->SaveAs("tgraphefficiency_ak4_ak8_vs_nosel_ak8_pt_0L_channel.pdf");
  delete c66;

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
 TLegend* leg64 = new TLegend(0.1,0.99,0.25,0.80);
        leg64->AddEntry(h_eff_ak15_deltaR_0_3_varbin,"AK15Jets,  #Delta R(H_{true},H_{AK15Jets})<0.3");
        leg64->AddEntry(h_eff_ak8_deltaR_0_3_varbin,"AK8Jets,  #Delta R(H_{true},H_{AK8Jets})<0.3");
        leg64->AddEntry(h_eff_ak4_deltaR_0_3_varbin,"resolved,  #Delta R(H_{true},H_{reco})<0.3");


        leg64->SetTextSize(0.03);
        leg64->SetTextFont(42);
        leg64->SetBorderSize(0);
        leg64->Draw("same");
  c64->SaveAs("tgraphefficiency_ak15_deltaR_0_3_varbin_0L_channel.pdf");
  delete c64;













//ak_msoftdrop->Fit(func,"R");
//ak15_msoftdrop->Fit(func_ak15,"R");

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
 

fnew->Close();

}






