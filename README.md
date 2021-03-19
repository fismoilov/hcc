# hcc
VHccAnalysis
akjet_efficiency_2L.cxx contains code to calculate efficiency in 2L channel 
akjet_efficiency_Wplus_minus_H125ToCC_WLNu_powheg_13_1L.cxx for 1L channel 

Plots are done with piece of code beginning from 

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
