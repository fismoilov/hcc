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

using namespace std;

 std::ostringstream NameString;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////analysis for AK15 AK8 AK4 Jets//////////////////////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////



void akjet() {


Float_t  LeadAK8_pt[50],LeadAK8_eta[50],LeadAK8_phi[50],LeadAK8_mass[50],LeadAK8_msoftdrop, HJ1_eta[50],HJ1_phi[50],HJ1_pt[50],HJ2_eta[50],
HJ2_phi[50], HJ2_pt[50],HJ1_HJ2_dR_noFSR[50],Jet_eta[50],Jet_mass[50],Jet_phi[50],Jet_pt[50],H_pt[50],H_eta[50],H_phi[50];

Float_t GenPart_eta[150],GenPart_mass[150],GenPart_phi[150],GenPart_pt[150];
 int  GenPart_genPartIdxMother[150], GenPart_pdgId[150], GenPart_status[150], GenPart_statusFlags[150];
UInt_t nGenPart;
Float_t LeadAK15_pt[50],LeadAK15_eta[50],LeadAK15_phi[50],LeadAK15_mass[50],LeadAK15_msoftdrop[50];
Float_t  nLeadAK15;


Float_t JetPt_1[50],JetPt_2[50];
Float_t DeepJet_CvsL_1[50], DeepJet_CvsB_1[50],DeepJet_CvsL_2[50],DeepJet_CvsB_2[50],H_mass[50],V_mass[50],V_pt[50],HVdPhi[50];

Float_t  nLeadAK8;
UInt_t nJet;

int isZmm,isZee;
int controlSample;
int hJetInd1,hJetInd2;


	TChain *t = new TChain("Events");
 // TChain *chain = new TChain("tree");
   for (int fileNum=1;fileNum < 12;fileNum++) {
      t->AddFile(Form("output_ZH125ToCC_ZLL_powheg_%d.root", fileNum));
  }

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
	TH1F *ak_pt = new TH1F("LeadAK8_pt","",100,400,1000);
	TH1F *ak_mass = new TH1F("LeadAK8_mass","",720,20,300);
	TH1F *ak_msoftdrop= new TH1F("LeadAK8_msoftdrop","",100,20,300);
        TH1F *ak_num = new TH1F("nLeadAK8","",100,0,50);

        TH1F *ak15_eta = new TH1F("LeadAK15_eta","",100,-3.5,3.5);
        TH1F *ak15_phi = new TH1F("LeadAK15_phi","",100,-3.5,3.5);
        TH1F *ak15_pt = new TH1F("LeadAK15_pt","",100,400,1000);
        TH1F *ak15_mass = new TH1F("LeadAK15_mass","",100,20,300);
        TH1F *ak15_msoftdrop = new TH1F("LeadAK15_msoftdrop","",100,20,300);
	
	TH2F *ak8_etaphi = new TH2F("ak8_etaphi","",50,-2.4,2.4,50,-3.5,3.5);
        TH2F *ak15_etaphi = new TH2F("ak15_etaphi","",50,-2.4,2.4,50,-3.5,3.5);



        TH1F *jet_eta = new TH1F("Jet_eta","",100,-3.5,3.5);
        TH1F *jet_phi = new TH1F("Jet_phi","",100,-3.5,3.5);
        TH1F *jet_pt = new TH1F("Jet_pt","",100,0,600);
        TH1F *jet_mass = new TH1F("Jet_mass","",100,20,600);


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
 TH1F *dR_Zboson_ak8 = new TH1F("dR_Zboson_ak8","",50,0,5);

       
int ak8jets=0;
int ak4jets=0;
int ak15jets=0;



//int ak4_selected_true=0;
//int ak4_selected_boosted_wrt_ak8=0;
//int ak4_selected_boosted_wrt_ak15=0;
//int ak8_selected_true=0;
//int ak15_selected_true=0;


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

		LeadAK8.SetPtEtaPhiM(LeadAK8_pt[0],LeadAK8_eta[0],LeadAK8_phi[0],LeadAK8_mass[0]);
		LeadAK15.SetPtEtaPhiM(LeadAK15_pt[0],LeadAK15_eta[0],LeadAK15_phi[0],LeadAK15_mass[0]);

                AK8Higgs.SetPtEtaPhiM(LeadAK8_pt[0],LeadAK8_eta[0],LeadAK8_phi[0],LeadAK8_mass[0]);
                AK15Higgs.SetPtEtaPhiM(LeadAK15_pt[0],LeadAK15_eta[0],LeadAK15_phi[0],LeadAK15_mass[0]);

                Hcc.SetPtEtaPhiM(H_pt[0],H_eta[0],H_phi[0],H_mass[0]);
		HJ1.SetPtEtaPhiM(HJ1_pt[0],HJ1_eta[0],HJ1_phi[0],Jet_mass[0]);
                HJ2.SetPtEtaPhiM(HJ2_pt[0],HJ2_eta[0],HJ2_phi[0],Jet_mass[0]);









/*
if((isZmm==1 || isZee==1) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
 DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0
&&controlSample==0 && //H_mass[0]>50 && H_mass[0]<200 &&
 V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)  {
	if(){
		ak4_selected_true++;
				    }
	if (LeadAK8.Pt()>200){ 
//std::cout<<"ak4jets ="<<ak4jets<<std::endl;
		ak8_selected_true++;
			     }
	if (LeadAK15.Pt()>200){
	        ak15_selected_true++;
                              }
}
*/

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////GenPart ZBoson////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
 for (int k=0; k<nGenPart;k++) {
//std::cout<<"GenPart_pdgId ="<<GenPart_pdgId[k]<<" GenPart_genPartIdxMother="<<GenPart_genPartIdxMother[k]<<
//"GenPart_pt ="<<GenPart_pt[k]<<"GenPart_mass ="<<GenPart_mass[k]<<std::endl;
 if(GenPart_pdgId[k]==23&&GenPart_genPartIdxMother[k]==0){
Zboson.SetPtEtaPhiM(GenPart_pt[k],GenPart_eta[k],GenPart_phi[k],GenPart_mass[k]);
//std::cout<<"nGenPart = "<< nGenPart <<"GenPart_pdgId ="<<GenPart_pdgId[k]<<" GenPart_genPartIdxMother="<<GenPart_genPartIdxMother[k]<<
//"GenPart_pt ="<<GenPart_pt[k]<<"GenPart_mass ="<<GenPart_mass[k]<<std::endl;
float R_Zboson_ak8;
if (LeadAK15.Pt()>400)  {
R_Zboson_ak8=Zboson.DeltaR(AK15Higgs);
dR_Zboson_ak8->Fill(R_Zboson_ak8);
			}
							}// end of cut
                               }// end of for loop
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////AK8Jets////////////////////////////// 
////////////////////////////////////////////////////////////////////////////////////////////////////////////
		 for(int k=0;k<nLeadAK8;k++)
			{
//std::cout<<"k = "<<k<<std::endl;
//std::cout<<"JetPt_1 = "<<k<<JetPt_1[0]<<" JetPt_2 = "<<JetPt_2[0]<<" DeepJet_CvsL_1 = "<<DeepJet_CvsL_1[0]<<" DeepJet_CvsB_1 = "<<DeepJet_CvsB_1[0]<<" DeepJet_CvsL_2 = "<<DeepJet_CvsL_2[0]<<" DeepJet_CvsB_2"<<DeepJet_CvsB_2[0]<<" H_mass = "<<H_mass[0]<<"V_mass = "<<V_mass[0]<<"V_pt = "<<V_pt[0]<<"HVdPhi = "<<HVdPhi[0]<<std::endl;
if((isZmm==1 || isZee==1) && JetPt_1[0]>20 && JetPt_2[0]>20 &&

 DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0 

&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)                                  {     

float R_ak8_ak4;
//std::cout<<"R_ak8_ak4 = "<<R_ak8_ak4<<std::endl;
R_ak8_ak4=AK8Higgs.DeltaR(Hcc);
dR_ak8_ak4->Fill(R_ak8_ak4); 


                   if (LeadAK8.Pt()>400)	        {
		         ak_num->Fill(nLeadAK8);	
			 ak_mass->Fill((LeadAK8).M());
                         ak_pt->Fill((LeadAK8).Pt());
                         ak_eta->Fill(LeadAK8.Eta());
                         ak_phi->Fill(LeadAK8.Phi());
			 ak_msoftdrop->Fill(LeadAK8_msoftdrop);		 
      			ak8_etaphi->Fill(LeadAK8.Eta(),LeadAK8.Phi());                       
                        higgs_ak8->Fill(Hcc.Pt(),LeadAK8.Pt());                      


     		                                         }
                                        } //end of cut     
                           }// end of for loop



////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK15Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
   for(int k=0;k<nLeadAK15;k++)
                        {
 
if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&

 DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0

&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5){

//std::cout<<"dR_ak15="<<dR_ak15<<std::endl;

float R_ak15_ak4;
R_ak15_ak4=AK15Higgs.DeltaR(Hcc);
dR_ak15_ak4->Fill(R_ak15_ak4);

 if (LeadAK15.Pt()>400){




//ak15jets++;
                        ak15_mass->Fill(LeadAK15.M());
                        ak15_pt->Fill((LeadAK15).Pt());
                        ak15_eta->Fill(LeadAK15.Eta());
                        ak15_phi->Fill(LeadAK15.Phi());
 // if (LeadAK15.M()<200&&LeadAK15.M()>50){
 ak15_msoftdrop->Fill(LeadAK15_msoftdrop[k]);//}
                ak15_etaphi->Fill(LeadAK15.Eta(),LeadAK15.Phi());                                                                                                  
														}
 				            		}// end of cut



}// end of for loop

////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////AK4Jets////////////////////////////// 
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

// not supported 

                        for(int k=0;k<nJet;k++)
               {





if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
 
 DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.40 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0

&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5){
 if (HJ1.Pt()+HJ2.Pt()>60)     {
 


if( ((HJ1+HJ2).M())>60&&((HJ1+HJ2).M())<200){

                        hj_mass->Fill((HJ1+HJ2).M());
}
 //                  ak4jets++;         
		         hja_pt->Fill((HJ1).Pt());
                         hja_eta->Fill(HJ1.Eta());
                         hja_phi->Fill(HJ1.Phi());
                         // change 17.12.2020 to hj1+hj2
                         hjb_pt->Fill((HJ1+HJ2).Pt());
                         hjb_eta->Fill(HJ2.Eta());
                         hjb_phi->Fill(HJ2.Phi());
                         hja->Fill(HJ1.Eta(),HJ1.Phi());
                         hjb->Fill(HJ2.Eta(),HJ2.Phi());

			if(LeadAK8.Pt()>350){ ak_jet->Fill(LeadAK8.M(),(HJ1+HJ2).M());
                                            }
		        	}
			}
                 }// end of for loop


        

}

// cout << "Number of AK8Jets = " << ak8jets;
// cout << "Number of AK15Jets = " << ak15jets;
//cout << "Number of AK4Jets = " << ak4jets;

/*
cout <<"ak4_selected_true =" <<ak4_selected_true;
cout <<"ak4_selected_boosted_wrt_ak8 =" << ak4_selected_boosted_wrt_ak8;
cout <<"ak4_selected_boosted_wrt_ak15 =" <<ak4_selected_boosted_wrt_ak15;
cout <<"ak8_selected_true =" <<ak8_selected_true;
cout <<"ak15_selected_true ="<<ak15_selected_true;

*/








TCanvas *canvas11 = new TCanvas("canvas11", "", 1900, 1100);
                canvas11->cd();

           higgs_ak8->GetYaxis()->SetTitle("NEvents");

               higgs_ak8 ->GetYaxis()->SetTitleSize(0.05);
             higgs_ak8 ->GetYaxis()->SetTitle("higgs_pt, GeV");
              higgs_ak8  ->GetYaxis()->SetTitleOffset(0.80);
               higgs_ak8 ->GetXaxis()->SetTitle("ak8_pt, GeV");
               higgs_ak8 ->GetXaxis()->SetTitleSize(0.05);
               higgs_ak8 ->GetXaxis()->SetTitleOffset(0.90);
 //              dR_ak8_ak4 ->SetLineWidth(2);
   //            dR_ak8_ak4 ->SetLineColor(5);
         gPad->SetGrid();
              higgs_ak8 ->SetTitle("higgs_ak8_pt");
               higgs_ak8 ->Draw("Colz");
           canvas11->SaveAs("higgs_ak8_pt.pdf");
                  delete canvas11;

TCanvas *canvas111 = new TCanvas("canvas11", "", 1900, 1100);
                canvas111->cd();

           dR_ak15_ak4->GetYaxis()->SetTitle("NEvents");

                dR_ak15_ak4->GetYaxis()->SetTitleSize(0.05);
              dR_ak15_ak4 ->GetYaxis()->SetTitle("ak15");
               dR_ak15_ak4 ->GetYaxis()->SetTitleOffset(0.80);
               dR_ak15_ak4 ->GetXaxis()->SetTitle("ak4");
               dR_ak15_ak4 ->GetXaxis()->SetTitleSize(0.05);
               dR_ak15_ak4 ->GetXaxis()->SetTitleOffset(0.90);
     //          dR_ak15_ak4 ->SetLineWidth(2);
       //        dR_ak15_ak4 ->SetLineColor(5);
         gPad->SetGrid();
               dR_ak15_ak4 ->SetTitle("dR_ak15_ak4");
                dR_ak15_ak4->Draw("Colz");
           canvas111->SaveAs("dR_ak15_ak4.pdf");
                  delete canvas111;


  TCanvas *canvas12 = new TCanvas("canvas1232567", "", 1900, 1100);
                canvas12->cd();
                ak_msoftdrop->GetYaxis()->SetTitle("NEvents");
//Double_t binContent = ak_msoftdrop->GetBinContent(bin);

//ak_msoftdrop->Smooth();
                ak_msoftdrop->GetYaxis()->SetTitleSize(0.05);
                ak_msoftdrop->GetYaxis()->SetTitleOffset(0.80);
                ak_msoftdrop->GetXaxis()->SetTitle("GeV");
                ak_msoftdrop->GetXaxis()->SetTitleSize(0.05);
                ak_msoftdrop->GetXaxis()->SetTitleOffset(0.90);
                ak_msoftdrop->SetLineWidth(2);
                ak_msoftdrop->SetLineColor(kOrange+7);
                ak_msoftdrop->SetFillColor(kOrange);
          gPad->SetGrid();
                ak_msoftdrop->SetTitle("msoftdrop");
         ak_msoftdrop->Draw("hist");
       /*  TLegend *leg52 = new TLegend(0.70,0.70,0.90,0.9,"msoftdrop");
                leg52->SetTextSize(0.05);
                leg52->SetBorderSize(0);
                leg52->SetFillColor(0);
        leg52->Draw();
         NameString.str("");

                NameString<<"AK8Jet";
                leg52->AddEntry(ak_msoftdrop,NameString.str().c_str(),"l");
*/
         canvas12->SaveAs("ak8_msoftdrop.pdf");
                delete canvas12;


 TCanvas *canvas13 = new TCanvas("canvas1232567", "", 1900, 1100);
                canvas13->cd();
                ak15_msoftdrop->GetYaxis()->SetTitle("NEvents");
                ak15_msoftdrop->GetYaxis()->SetTitleSize(0.05);
                ak15_msoftdrop->GetYaxis()->SetTitleOffset(0.80);
                ak15_msoftdrop->GetXaxis()->SetTitle("GeV");
                ak15_msoftdrop->GetXaxis()->SetTitleSize(0.05);
                ak15_msoftdrop->GetXaxis()->SetTitleOffset(0.90);
                ak15_msoftdrop->SetLineWidth(2);
                ak15_msoftdrop->SetLineColor(kOrange+7);
                ak15_msoftdrop->SetFillColor(kOrange);
          gPad->SetGrid();
                ak15_msoftdrop->SetTitle("msoftdrop");
         ak15_msoftdrop->Draw("hist");
         canvas13->SaveAs("ak15_msoftdrop.pdf");
                delete canvas13;

 TCanvas *canvas7 = new TCanvas("canvas732567", "", 1900, 1100);
                canvas7->cd();
                ak_pt->GetYaxis()->SetTitle("NEvents");

                ak_pt->GetYaxis()->SetTitleSize(0.05);
                ak_pt->GetYaxis()->SetTitleOffset(0.80);
                ak_pt->GetXaxis()->SetTitle("P_{t},GeV");
                ak_pt->GetXaxis()->SetTitleSize(0.05);
                ak_pt->GetXaxis()->SetTitleOffset(0.90);
                ak_pt->SetLineWidth(2);
                ak_pt->SetLineColor(kOrange+7);
                ak_pt-> SetFillColor(kOrange);

 gPad->SetGrid();

    ak_pt->SetTitle("Transverse Momentum");
                ak_pt->Draw("hist");
 canvas7->SaveAs("ak8_pt.pdf");
                delete canvas7;

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
		ak15_pt->SetFillColor(kOrange);
                gPad->SetGrid();
                ak15_pt->SetTitle("Transverse Momentum");
                ak15_pt->Draw("hist");



         canvas14->SaveAs("ak15_pt.pdf");
                delete canvas14;



  TCanvas *canvas1 = new TCanvas("canvas1", "", 1900, 1100);
                canvas1->cd();
                hj_mass->GetYaxis()->SetTitle("NEvents");
                hj_mass->GetYaxis()->SetTitleSize(0.05);
                hj_mass->GetYaxis()->SetTitleOffset(0.80);
                hj_mass->GetXaxis()->SetTitle("M, GeV");
                hj_mass->GetXaxis()->SetTitleSize(0.05);
                hj_mass->GetXaxis()->SetTitleOffset(0.90);
                hj_mass->SetLineWidth(2);
                hj_mass->SetLineColor(kOrange+7);
		hj_mass->SetFillColor(kOrange);

                 gPad->SetGrid();
                hj_mass->SetTitle("Resolved-Jet Topology Mass");
                hj_mass->Draw("hist");
canvas1->SaveAs("hjresolved_mass.pdf");
              delete canvas1;

 TCanvas *canvas9 = new TCanvas("canvas9", "", 1900, 1100);
                canvas9->cd();
                hjb_pt->GetYaxis()->SetTitle("NEvents");
                hjb_pt->GetYaxis()->SetTitleSize(0.05);
                hjb_pt->GetYaxis()->SetTitleOffset(0.80);
                hjb_pt->GetXaxis()->SetTitle("P_{t},GeV");
                hjb_pt->GetXaxis()->SetTitleSize(0.05);
                hjb_pt->GetXaxis()->SetTitleOffset(0.90);
                hjb_pt->SetLineWidth(2);
                hjb_pt->SetLineColor(kOrange+7);
                hjb_pt-> SetFillColor(kOrange);
                gPad->SetGrid();
                hjb_pt->SetTitle("Transverse Momentum");
                hjb_pt->Draw("hist");
         canvas9->SaveAs("hj1_2_pt.pdf");
                delete canvas9;


 TCanvas *canvas141 = new TCanvas("canvas1232567", "", 1900, 1100);
                canvas141->cd();
               dR_ak8_ak4 ->GetYaxis()->SetTitle("NEvents");
              dR_ak8_ak4  ->GetYaxis()->SetTitleSize(0.05);
               dR_ak8_ak4 ->GetYaxis()->SetTitleOffset(0.80);
              dR_ak8_ak4  ->GetXaxis()->SetTitle("dR");
             dR_ak8_ak4   ->GetXaxis()->SetTitleSize(0.05);
               dR_ak8_ak4 ->GetXaxis()->SetTitleOffset(0.90);
              dR_ak8_ak4  ->SetLineWidth(2);
              dR_ak8_ak4  ->SetLineColor(kOrange+7);
               dR_ak8_ak4 ->SetFillColor(kOrange);
          gPad->SetGrid();
              dR_ak8_ak4  ->SetTitle("dR");
        dR_ak8_ak4 ->Draw("hist");
         canvas141->SaveAs("deltar_ak8_ak4.pdf");
                delete canvas141;
 TCanvas *canvas1411 = new TCanvas("canvas1232567", "", 1900, 1100);
                canvas1411->cd();
               dR_ak15_ak4 ->GetYaxis()->SetTitle("NEvents");
              dR_ak15_ak4  ->GetYaxis()->SetTitleSize(0.05);
               dR_ak15_ak4 ->GetYaxis()->SetTitleOffset(0.80);
              dR_ak15_ak4  ->GetXaxis()->SetTitle("dR");
             dR_ak15_ak4   ->GetXaxis()->SetTitleSize(0.05);
               dR_ak15_ak4 ->GetXaxis()->SetTitleOffset(0.90);
              dR_ak15_ak4  ->SetLineWidth(2);
              dR_ak15_ak4  ->SetLineColor(kOrange+7);
               dR_ak15_ak4 ->SetFillColor(kOrange);
          gPad->SetGrid();
              dR_ak15_ak4  ->SetTitle("dR");
        dR_ak15_ak4 ->Draw("hist");
         canvas1411->SaveAs("deltar_ak15_ak4.pdf");
                delete canvas1411;

TCanvas *canvas20 = new TCanvas("canvas1232567", "", 1900, 1100);
                canvas20->cd();
               dR_Zboson_ak8 ->GetYaxis()->SetTitle("NEvents");
              dR_Zboson_ak8 ->GetYaxis()->SetTitleSize(0.05);
               dR_Zboson_ak8 ->GetYaxis()->SetTitleOffset(0.80);
              dR_Zboson_ak8 ->GetXaxis()->SetTitle("dR");
             dR_Zboson_ak8   ->GetXaxis()->SetTitleSize(0.05);
               dR_Zboson_ak8 ->GetXaxis()->SetTitleOffset(0.90);
              dR_Zboson_ak8  ->SetLineWidth(2);
              dR_Zboson_ak8 ->SetLineColor(kOrange+7);
               dR_Zboson_ak8 ->SetFillColor(kOrange);
          gPad->SetGrid();
              dR_Zboson_ak8 ->SetTitle("dR");
        dR_Zboson_ak8 ->Draw("hist");
         canvas20->SaveAs("deltar_Zboson_ak8.pdf");
                delete canvas20;





dR_ak8_ak4->Write();
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
 

        hj_mass->Write();
        hja_pt->Write();
        hja_eta->Write();
        hja_phi->Write();
        
        hjb_pt->Write();
        hjb_eta->Write();
        hjb_phi->Write();	
        ak_jet->Write();
	hja->Write();
	hjb->Write();
fnew->Close();

}






