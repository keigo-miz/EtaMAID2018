void show() {
  TFile *file0 = new TFile("rt/out.root");
  TGraph *tg0 = (TGraph *)file0->Get("tg0");
  TGraph *tg1 = (TGraph *)file0->Get("tg1");

  TGraph *ref0 = new TGraph("e0p_re.dat");
  TGraph *ref1 = new TGraph("e0p_im.dat");
  ref0->SetMarkerStyle(20);
  ref1->SetMarkerStyle(20);
  ref0->SetMarkerColor(2);
  ref1->SetMarkerColor(2);

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 600);
  c1->Divide(2, 1);
  TH1 *frame0 = c1->cd(1)->DrawFrame(1480.01, -15.0, 1610.0, 0.0);
  frame0->SetTitle("Re E_{0#plus}  [N(1535)]");
  frame0->GetXaxis()->SetTitle("W (MeV)");
  tg0->Draw("l");
  ref0->Draw("p");
  TH1 *frame1 = c1->cd(2)->DrawFrame(1480.01, -0.0, 1610.0, 25.0);
  frame1->SetTitle("Im E_{0#plus}  [N(1535)]");
  frame1->GetXaxis()->SetTitle("W (MeV)");
  tg1->Draw("l");
  ref1->Draw("p");
}
