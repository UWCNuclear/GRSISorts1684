const int num_cores = 64; //number of cores in TIGRESS detector
Double_t num_bins = 8192; // 4096; //number of bins to read in from histograms
Double_t min_bin = 0; 
Double_t max_bin = 32768; //16384;

Int_t num_known_sources = 4; //number of sources that can be used for calibration

//names of usable sources
string source_names[] = {
"eu152", "ba133", "co56", "co60"
};

//number of peaks used for each source's calibration
Int_t num_peaks[] = {
10, 4, 5, 2
};

//sorts an array into ascending order
//preserves correspondance between the main array and the 'carry' array
void insertion_sort(Double_t array[], Double_t carry[], int len) {

	int i;
	for (i = 1; i < len; i++) {
		Double_t val = array[i];
		Double_t carry_val = carry[i];

		for (int j = i - 1; j >= 0; j--) {
			if (array[j+1] < array[j]) {
				array[j+1] = array[j];
				carry[j+1] = carry[j];

				array[j] = val;
				carry[j] = carry_val;				
			} else {
				break;
			}//else
		}//for
	}//for
}//insertion_sort

//prints an array of Doubles
void print_array(Double_t array[], Int_t len) {
	cout << "{";
	for (int i = 0; i < len; i++) {
		cout << array[i] << ", ";
	}//for
	cout << "}" << endl;
}//print_array

//loads the histograms from analysis root file
void load_histograms(char analysis_filepath[], char cal_filepath[], TH1F *hist[], Int_t num_cores) {
	
	//opening analysis root file
	TFile *input_file = new TFile(analysis_filepath, "READ");
	if (!input_file->IsOpen()) {
		cout << "Cannot open input file!" << endl;
		return 0;
	}//if
	TChain *analysis = (TChain *) input_file->Get("AnalysisTree");
	TTree *tree = (TTree *) analysis->GetTree();

	//opening calibration file that organizes data into channels
	TChannel::ReadCalFile(cal_filepath);

	Int_t num_entries = analysis->GetEntries()/10;	
	TTigress *tigr = 0;
	analysis->SetBranchAddress("TTigress", & tigr);

	//initializing histograms
	for (int i = 0; i < num_cores; i++) {
		hist[i] = new TH1F(Form("hist%d", i), Form("Gamma Singles Crystal %1.1i", i), num_bins, min_bin, max_bin);
	}//for

	cout << "Histograms created." << endl;

	//filling the histograms with analysis file data
	for (int i = 0; i < num_entries - 1; i++) {
		tree->GetEntry(i);
		for (int j = 0; j < tigr->GetMultiplicity(); j++) {
			TTigressHit *tigr_hit = tigr->GetTigressHit(j);
			hist[tigr_hit->GetArrayNumber()]->Fill(tigr_hit->GetEnergy());
		}//for
		
		if (i % 10000 == 0) {
			cout << setiosflags(ios::fixed) << "TIGRESS Entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
		}//if
	}//for

	cout << setiosflags(ios::fixed) << "TIGRESS Entry " << num_entries << " of " << num_entries << ", 100% complete" << "\r" << flush;		
}//load_histograms

//using gamma ray spectra and actual energy peaks, calculates a linear equation to calibrate the detector's outputs
void linear_calibration(TList *list, TH1F *hist[], Int_t num_peaks_used, Double_t energy[], Double_t energy_er[], Double_t gains[num_cores], Double_t offsets[num_cores]) {
	
	//arrays to store the centroids of the peaks
	Double_t centroids[num_cores][num_peaks_used];
	Double_t centroids_er[num_cores][num_peaks_used];

	ofstream fwhm("sigma_76.txt");

	for (int i = 0; i < num_cores; i++) {
		hist[i]->SetAxisRange(100, max_bin-100,"X");
		hist[i]->SetBinContent(1, 0);		
		Double_t intr = hist[i]->Integral(min_bin, max_bin, "width");		

		cout << intr << endl;

		if (intr < 1000) {
			cout << i << " FAILED!" << endl;
			for (int j = 0; j < num_peaks_used; j++) {
				centroids[i][j] = -1;
			}//for			
			continue;
		}//if

		//roughly locates the peaks in the spectrum		
		TSpectrum *spec = new TSpectrum(2 * num_peaks_used);
		Int_t num_found = spec->Search(hist[i], 2, "", 0.5);

		cout << "Found " << num_found << " peaks in histogram." << endl;
		
		//if too many or too few peaks have been found,
		//something has gone wrong and we move on to the next core
		if (num_found != num_peaks_used) {
			for (int j = 0; j < num_peaks_used; j++) {
				centroids[i][j] = -1;
			}//for
			continue;
		}//if
		
		//we fit a gaussian distribution around each peak
		//using the TSpectrum's data as a starting point		
		TF1 *fit[num_peaks_used];
		Double_t* x_pos_array = spec->GetPositionX();
		
		for (int j = 0; j < num_found; j++) {
			
			Double_t x_pos = x_pos_array[j];
			Int_t bin = hist[i]->GetXaxis()->FindBin(x_pos);
			Double_t y_pos = hist[i]->GetBinContent(bin);

			fit[j] = new TF1(Form("Fit %i-%i", i, j), 
							"[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))",
							x_pos - 50, x_pos + 50);

			fit[j]->SetParameters(y_pos, x_pos, 1, 15, 1, -1);
			fit[j]->SetParLimits(0, 10, 1e6); //area
			fit[j]->SetParLimits(1, x_pos - 15, x_pos + 10); //centroid
			fit[j]->SetParLimits(2, 0.2, 15); //sigma
			fit[j]->SetParLimits(4, 0.1, 100); //magnitude of step in background noise
			fit[j]->SetParLimits(5, -10, -0.1); //background noise constant

			//fitting the equation and storing the calculated centroids			
			hist[i]->Fit(fit[j], "RQ+");
			centroids[i][j] = fit[j]->GetParameter(1);
			centroids_er[i][j] = fit[j]->GetParError(1);
			
			cout << fit[j]->GetParameter(2) << endl;
			fwhm << i << " " << fit[j]->GetParameter(2) << " ";

			cout << "CENTROID PROCESSED: Graph " << i << " Guess: " << x_pos << " Actual: " << centroids[i][j] << endl;
		}//for

		fwhm << endl;

		//to make sure the centroids match up to the correct energies
		//we sort the centroids in ascending order
		//making sure the centroids_er keep the correspondance
		insertion_sort(centroids[i], centroids_er[i], num_peaks_used);
	}//for
	
	//we graph the centroids vs their corresponding energies
	//and fit a linear equation on the points
	for (int i = 0; i < num_cores; i++) {
		
		//if the histogram is empty, skip this core
		if (centroids[i][0] == -1) {
			gains[i] = -1;
			offsets[i] = -1;
			continue;
		}//if

		TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy, centroids_er[i], energy_er);
		gr->Draw("AP");
		list->Add(gr); //adding the graph to be saved to root file later

		TF1 *coeffFit = new TF1("coeffFit", "[0] + [1]*x");
		gr->Fit(coeffFit, "Q+");
		
		//storing linear equation parameters
		gains[i] = coeffFit->GetParameter(1);
		offsets[i] = coeffFit->GetParameter(0);

	}//for

	fwhm.close();

}//linear_calibration

//Main method to be executed by GRSISort
void linear_energy() {	

	/*******Initialization of Source Energy Peak Data*******

	Double_t** source_energy = new Double_t*[num_known_sources];
	source_energy[0] = new Double_t[9]{121.7818, 244.6976, 344.2789, 411.1171, 778.9066, 867.383, 964.082, 1112.080, 1299.148, 1408.013}; //eu152
	source_energy[1] = new Double_t[4]{276.4, 302.851, 356.013, 383.848}; //ba133
	source_energy[2] = new Double_t[5]{1771.3567, 2034.7097, 2598.500, 3253.5030, 3451.232}; //co56
	source_energy[3] = new Double_t[2]{1173.240, 1332.508}; //co60
	
	Double_t** source_energy_er = new Double_t*[num_known_sources];
	source_energy_er[0] = new Double_t[9]{0.0003, 0.0008, 0.0012, 0.0012, 0.0024, 0.003, 0.018, 0.003, 0.008, 0.003}; //eu152
	source_energy_er[1] = new Double_t[4]{0.0021, 0.0016, 0.0017, 0.0012}; //ba133
	source_energy_er[2] = new Double_t[5]{0.0039, 0.0047, 0.004, 0.0044, 0.004}; //co56
	source_energy_er[3] = new Double_t[2]{0.003, 0.004}; //co60

	************************************************************/

	Double_t co60_ener[2] = {1173.240, 1332.508};
	Double_t co60_ener_e[2] = {0.003, 0.004};

	TList *list = new TList;
	
	TH1F *lin_hist[num_cores]; //histograms for each core

	//linear calibration equation parameters for each core
	Double_t lin_gains[num_cores];
	Double_t lin_offsets[num_cores];

	//getting calibration source
	string source = "co60";
	//cout << "\nEnter the gamma source used in calibration. Use the source's first two letters and mass number (e.g. eu152, ba133).\n";
	//getline(cin, source);
	
	//getting analysis root file	
	string rootfile;
	cout << "\nEnter the file path for the analysis tree.\n";
	getline(cin, rootfile);

	char empty_gains_cal[70] = "EmptyGains.cal";

	//loading in histograms from analysis root file
	load_histograms(&rootfile[0], empty_gains_cal, lin_hist, num_cores);

	//adding histograms to be saved later
	for (int i = 0; i < num_cores; i++) {
		if (lin_hist[i]->GetMaximum() != 0) {		
			list->Add(lin_hist[i]);
		}//if
	}//for	
	
	//running through all possible sources and, if it's the source being used, use its energies for calibration	
	for (int i = 0; i < num_known_sources; i++) {
		if (source == source_names[i]) {
			linear_calibration(list, lin_hist, num_peaks[i], co60_ener, co60_ener_e, lin_gains, lin_offsets);
			break;		
		}//if
	}//for

	//write out linear equations parameters
	ofstream coeff_fw("lin_energy_coeff.txt");
	
	for (int i = 0; i < num_cores; i++) {
		coeff_fw << lin_gains[i] << endl;
		coeff_fw << lin_offsets[i] << endl;
	}//for

	coeff_fw.close();
		
	//write out histograms and linear fit graphs for review
	TFile * outfile = new TFile("lin_cal.root", "RECREATE");
	outfile->cd();
	list->Write();
	outfile->Close();

}//linear_energy_RF
