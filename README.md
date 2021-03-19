# hcc
VHccAnalysis
akjet_efficiency_2L.cxx contains code to calculate efficiency in 2L channel 
akjet_efficiency_Wplus_minus_H125ToCC_WLNu_powheg_13_1L.cxx in 1L channel 
akjet_efficiency_ZH125ToCC_ZNuNu_0L.cxx in 0L channel






**Numerators and denominators to calculate integrals for Efficiency are commented as** 

- Denominator for resolved jet topology efficiency
- Denominator for boosted jet topology efficiency
- Numerator for resolved jet topology efficiency
- Numerator for boosted jet topology efficiency




**Calculation of signal lost is done in 2L channel** 
    
   
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////Numerator for Lost of signal search////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
       **Numerator**            
    ////// check fractions of selected AK4s and AK15 and AK8//////
    

	if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
		DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0
		&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)                                  {


     if (Hcc.Pt()>200&&LeadAK15.Pt()>300){
     frac_ak4_ak15->Fill(LeadAK15.Pt());
     frac_ak4_ak15_ak4->Fill(LeadAK15.Pt()); 
             }
     } //end of cut 
       ** Denominator**
       if((isZmm || isZee) && JetPt_1[0]>20 && JetPt_2[0]>20 &&
		DeepJet_CvsL_1[0]>0.225 && DeepJet_CvsB_1[0]>0.4 && DeepJet_CvsL_2[0]>0.0 && DeepJet_CvsB_2[0]>0.0
		&&controlSample==0 && H_mass[0]>50 && H_mass[0]<200 && V_mass[0]>75 && V_mass[0]<105 && V_pt[0]>=60 && HVdPhi[0]>2.5)                                  {

		
                        h_pt->Fill(Hcc.Pt());
                        }

**Plots and integral calculation**

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

	int num_int_frac_ak4_ak15=frac_ak4_ak15->Integral();
	cout<<"num_int_frac_ak4_ak15  ="<<num_int_frac_ak4_ak15<<endl;

	int num_int_3_h_pt=h_pt->Integral();
	cout<<"num_int_3_h_pt  ="<<num_int_3_h_pt<<endl;

	double  frac_ak4_ak15_h_pt;
	
	frac_ak4_ak15_h_pt=num_int_frac_ak4_ak15/num_int_3_h_pt;
	cout<<"Lost of signal in frac_ak4_ak15_h_pt = "<<frac_ak4_ak15_h_pt<<endl;

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


**For Efficiency**
**Plots are done with piece of code beginning from** 

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
