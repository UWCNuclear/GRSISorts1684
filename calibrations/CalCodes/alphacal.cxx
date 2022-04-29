// g++ alphacal.cxx -std=c++0x -I$GRSISYS/include -L$GRSISYS/libraries `grsi-config --cflags --all-libs` -lTDetector -lTGRSIDetector -lTSuppressed `root-config --cflags --libs` -lTreePlayer -lSpectrum -ltbb -o ACal

// Used to calibrate a silicon detector with a triple alpha source. Uses a peak search to find the peaks, fits them at the same time and extracts the centroids. It plots this in a TGraph and performs a calibration with the parameters output into the terminal 
#include <iostream>
#include <iomanip>
#include "TCutG.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TSpectrum.h"
#include "TChannel.h"
#include "TFragment.h"
#include "TMath.h"
#include <cmath>
using namespace std;

void peak_search(char const* infile, char const* calfile, char const* outfile) {

  Double_t si_nBins = 8192;
  Double_t si_min = 0;
  Double_t si_max = 16384;

  TList * list = new TList;

  TH1D * Si_singles[5000];
  char siname[20];
  for (int iii = 0; iii < 2000; iii++) {
    sprintf(siname, "Channel%d", iii);
    Si_singles[iii] = new TH1D(siname, Form("Si Singles d%1.1i", iii), si_nBins, si_min, si_max);
  }

  TFile * inputfile = new TFile(infile, "READ");
  if (!inputfile->IsOpen()) {
    printf("Opening file failed, aborting\n");
    return;
  }

  TChain * FragmentTree = (TChain * ) inputfile->Get("FragmentTree");
  printf("%i tree files, details:\n", FragmentTree->GetNtrees());
  FragmentTree->ls();
  TTree * tree = (TTree * ) FragmentTree->GetTree();
  Int_t nentries = FragmentTree->GetEntries();
  TFragment * frag = 0;
  FragmentTree->SetBranchAddress("TFragment", & frag);
  TChannel::ReadCalFile(calfile);
  int npeaks = 20;

  printf("Begin sort\n");
  int one;
  for (int jentry = 0; jentry < (nentries); jentry++) {
    tree->GetEntry(jentry);
    if (frag->GetChannelNumber() > 0) Si_singles[frag->GetChannelNumber()]->Fill(frag->GetEnergy());
    if (jentry % 10000 == 0) cout << setiosflags(ios::fixed) << "Si Entry " << jentry << " of " << nentries << ", " << 100 * jentry / nentries << "% complete" << "\r" << flush;
  }

  for (int iii = 0; iii < 2000; iii++) {
    for(int jjj=0; jjj<200; jjj++)Si_singles[iii]->SetBinContent(jjj, 0);
    if (Si_singles[iii]->Integral(0, si_nBins) > 10) list->Add(Si_singles[iii]);
  }

  cout << "Si Entry " << nentries << " of " << nentries << ", 100% Complete!\n";
  cout << "Histograms written, sorting complete" << endl;

  double si_calib_en[7] = {5156.59,5485.56,5804.82};
  double si_calib_err[7] = {0.14,0.12,0.05};

  double si_Array[7];
  double si_eArray[7];
  double sigma;
  double gc1;
  double cent =0;
  TF1 * si[2000];
  TF1 * si_ccurve[2000];
  TGraphErrors * si_gr[2000];
  for (int qq = 480; qq < 1024; qq++) {
    if (Si_singles[qq]->Integral(0, si_nBins) < 10) cout << qq << " " << "-1" << "\t" << "0" << "\t" << endl;
    if (Si_singles[qq]->Integral(0, si_nBins) < 10) continue;
    TSpectrum * s = new TSpectrum(2 * npeaks);
    Int_t nfound = s->Search(Si_singles[qq], 2, "", 0.3);
    Double_t * xpeaks = s->GetPositionX();
    if (nfound == 3) {
    //  si[qq] = new TF1("gp", "[3]*(exp(-((x-[4])^2/(2*[2]^2))))*(1+ROOT::Math::erf([1]*((x-[4]))/([2]*pow(2,0.5))))+[5]*(exp(-((x-[6])^2/(2*[2]^2))))*(1+ROOT::Math::erf([1]*((x-[6]))/([2]*pow(2,0.5))))+[7]*(exp(-((x-[8])^2/(2*[2]^2))))*(1+ROOT::Math::erf([1]*((x-[8]))/([2]*pow(2,0.5))))+[0]", 1000, 2000);
si[qq] = new TF1("gp", "[0]+[4]*((((exp(-((x-[5])^2/(2*[3]^2)))))*[2])+(exp((x-[5])/([1]))*(ROOT::Math::erfc(((x-[5])/([3]*2^(0.5)))+([3]/([1]*2^(0.5)))))*(1-[2])))+[6]*((((exp(-((x-[7])^2/(2*[3]^2)))))*[2])+(exp((x-[7])/([1]))*(ROOT::Math::erfc(((x-[7])/([3]*2^(0.5)))+([3]/([1]*2^(0.5)))))*(1-[2])))+[8]*((((exp(-((x-[9])^2/(2*[3]^2)))))*[2])+(exp((x-[9])/([1]))*(ROOT::Math::erfc(((x-[9])/([3]*2^(0.5)))+([3]/([1]*2^(0.5)))))*(1-[2])))",100,2000);

    for (int p = 0; p < nfound; p++) {
        Double_t xp = xpeaks[p];
        si[qq]->SetParLimits(4 + 2 * p, 10, 1e6);
        si[qq]->SetParameter(4 + 2 * p, 1000);
        si[qq]->SetParLimits(5 + 2 * p, xp - 10, xp + 10);
        si[qq]->SetParameter(5 + 2 * p, xp);
      }
      si[qq]->SetParLimits(1, 0.1, 10);
      si[qq]->SetParLimits(2, 0.001, 1);
      si[qq]->SetParLimits(3, 1, 15);
      si[qq]->SetParameter(1, 5);
      si[qq]->SetParameter(2, 0.5);
      si[qq]->SetParameter(3, 5);
      Si_singles[qq]->Fit(si[qq], "QR+");
      si_Array[0] = si[qq]->GetMaximumX(si[qq]->GetParameter(5)-10,si[qq]->GetParameter(5)+10);
      si_eArray[0] = si[qq]->GetParError(5);
      si_Array[1] = si[qq]->GetMaximumX(si[qq]->GetParameter(7)-10,si[qq]->GetParameter(7)+10);
      si_eArray[1] = si[qq]->GetParError(7);
      si_Array[2] = si[qq]->GetMaximumX(si[qq]->GetParameter(9)-10,si[qq]->GetParameter(9)+10);
      si_eArray[2] = si[qq]->GetParError(9);

    }
    if (nfound != 3) continue;
    sort(si_Array, si_Array + nfound);
    si_gr[qq] = new TGraphErrors(nfound, si_Array, si_calib_en, si_eArray, si_calib_err);
    si_gr[qq]->SetTitle(Form("Si Calibration_Curve_%d", qq));
    si_gr[qq]->Draw("A*");
    si_ccurve[qq] = new TF1("curve", "pol1");
    si_gr[qq]->Fit(si_ccurve[qq], "Q");
    if (nfound == 3) cout << qq << " " << si_ccurve[qq]->GetParameter(0) << "\t" << si_ccurve[qq]->GetParameter(1) << "\t" << sigma << endl;
    if (nfound != 3) cout << qq << " " << "-1" << "\t" << "0" << "\t" << "Peaks found " << nfound << endl;
    list->Add(si_gr[qq]);
   // while (cin.get() != '\n');
  }

  cout << "Writing histograms to " << outfile << endl;
  TFile * myfile = new TFile(outfile, "RECREATE");
  myfile->cd();
  list->Write();
  myfile->Close();

}

int main(int argc, char **argv){
	
	char const *ffile;
	char const *outfile;
	char const *calfile;
	printf("Starting sortcode\n");

	if(argc == 1)
	{
		cout << "Insufficient arguments, provide argument tree files" << endl;
		return 0;
	}

	else if(argc == 2)
	{
		ffile = argv[1];
		calfile = "/data2/sgillespie/TigressCodes/CalibrationFile.cal";
		outfile = "SiCal.root";
	}
	else if(argc == 3)
	{
		ffile   = argv[1];
		calfile = argv[2];
		outfile = "SiCal.root";
	}
	else if(argc == 4)
	{
		ffile   = argv[1];
		calfile = argv[2];
		outfile = argv[3];
	}
	else if(argc > 4)
	{
		printf("Too many arguments\n");
		return 0;
	}

	printf("Input file:%s\nCalibration file: %s\nOutput file: %s\n",ffile,calfile,outfile);

	peak_search(ffile,calfile, outfile);
	
	return 0;
}
