const int num_cores = 64; //number of cores in TIGRESS detector
Double_t num_bins = 8192; // 4096; //number of bins to read in from histograms
Double_t min_bin = 0; 
Double_t max_bin = 32768; //16384;

Int_t num_known_sources = 4; //number of sources that can be used for calibration

//names of usable sources
string source_names[] = {
"eu152", "ba133", "co56", "co60"
};

//{457.827, 244.698}; {646.354, 344.279}; {835.381, 411.117}; {1468.05, 778.907}; {1634.71, 867.383}; {1817.92, 964.082}; {2097.63, 1112.08}; {2445.63, 1299.15}; {2650.42, 1408.01}; {3338.08, 1771.36}; {3834.16, 2034.71}; {4890, 2598.5}; {6032.89, 3253.5}; {6350.24, 3451.23}; 

//IDs for each source
//makes it simple to know which source is being used
Int_t source_ids[] = {
2, 3, 5, 7, 11, 13, 17, 23, 29
};

//number of peaks used for each source's calibration
Int_t num_peaks[] = {
8, 4, 4, 2
};

//prints an array of Doubles
void print_array(Double_t array[], Int_t len) {
	cout << "{";
	for (int i = 0; i < len; i++) {
		cout << array[i] << ", ";
	}//for
	cout << "}" << endl;
}//print_array

//adds second array to the end of the first one
//note: first array must have enough space for second
void move_into_array(Double_t to[], Double_t from[], int start, int len) {
	
	for (int i = start; i < start + len; i++) {
		to[i] = from[i - start];
	}//for

}//move_into_array

//searches array of strings for a certain string
Int_t search_array(string array[], string search, int len) {
	for (int i = 0; i < len; i++) {
		if (search == array[i]) {
			return i;
		}//if
	}//for
	
	return -1;
}//search_array

//loads the histograms from analysis root file
void load_histograms(char analysis_filepath[], char cal_filepath[], TH1F *hist[], Int_t num_cores, Int_t source_count) {

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
		
	Int_t num_entries = analysis->GetEntries();
	TTigress *tigr = 0;
	analysis->SetBranchAddress("TTigress", & tigr);

	//initializes histograms
	for (int i = 0; i < num_cores; i++) {
		hist[i] = new TH1F(Form("Histogram - %d", i), Form("Gamma Singles Crystal %1.1i", i), num_bins, min_bin, max_bin);
	}//for

	cout << "Histograms created." << endl;
	
	//filling histograms with data from analysis root file
	for (int i = 0; i < (num_entries - 1); i++) {
		tree->GetEntry(i);
		for (int j = 0; j < tigr->GetMultiplicity(); j++) {
			TTigressHit *tigr_hit = tigr->GetTigressHit(j);
			hist[tigr_hit->GetArrayNumber()]->Fill(tigr_hit->GetEnergy());
		}//for
		
		if (i % 50000 == 0) {
			cout << setiosflags(ios::fixed) << "TIGRESS Entry " << i << " of " << num_entries << ", " << 100 * (i) / num_entries << "% complete" << "\r" << flush;
		}//if
	}//for

	//for (int i = 0; i < 64; i++)hist[i]->RebinX(2);
	cout << setiosflags(ios::fixed) << "TIGRESS Entry " << num_entries << " of " << num_entries << ", 100% complete" << "\r" << flush;	
	for(int j = 0; j < 64; j++) {
	for(int i = 0; i < num_bins; i++) {
		double bc = hist[j]->GetBinContent(i);
		if(bc==0) {
		hist[j]->SetBinError(i,1);
		} 
		else {
		hist[j]->SetBinError(i, pow(bc,0.5));
		}
	}	
	}
}//load_histograms

//fits an equation on the centroids vs energies for each core of the TIGRESS detector
void fit_equation(TList *list, Double_t energy[], Double_t energy_er[], Double_t** centroids, Double_t** centroids_er, Int_t num_peaks_used,
					Double_t thirds[num_cores], Double_t gains[num_cores], Double_t offsets[num_cores], char use_quad) {
	
	string equation;

	if (use_quad == 'Y' || use_quad == 'y') {
		equation = "[0] + [1]*x + [2]*x*x";
	} else {
		equation = "[0] + [1]*x";
	}//else

	ofstream quad_fit("quad_fit_points.txt");
	for (int i = 0; i < num_peaks_used; i++) {
		quad_fit << energy[i] << " " << energy_er[i] << endl;
	}//for

	for (int i = 0; i < num_cores; i++) {
			
		for (int j = 0; j < num_peaks_used; j++) {
			quad_fit << centroids[i][j] << " " << centroids_er[i][j] << endl;
		}//for

		//if this centroid is null, then skip this core	
		if (centroids[i][0] == -1 || centroids_er[i][0] == -1) {
			thirds[i] = 0;
			gains[i] = 0;
			offsets[i] = -1;
			continue;
		}//if
		//creating graph of centroids vs energy
		TGraphErrors *gr = new TGraphErrors(num_peaks_used, centroids[i], energy, centroids_er[i], energy_er);
		gr->SetName(Form("Equation Fit - %i", i));
//		gr->Draw("AP");

		list->Add(gr); //saving graph for future 

		//fitting quadratic on data points
		TF1 *coeffFit = new TF1("coeffFit", &equation[0]);
		gr->Fit(coeffFit, "Q");
		
		//saving curve parameters
		if (use_quad == 'Y' || use_quad == 'y') {		
			thirds[i] = coeffFit->GetParameter(2);
		} else {
			thirds[i] = 0;
		}//else
		gains[i] = coeffFit->GetParameter(1);
		offsets[i] = coeffFit->GetParameter(0);
	
		double x[num_peaks_used];
		double resid[num_peaks_used];
		
		for (int j = 0; j < num_peaks_used; j++) {
	      x[j] = j;
		  resid[j] = energy[j] - coeffFit->Eval(centroids[i][j]);  
		}//for
		
		TGraph *resid_g = new TGraph(num_peaks_used, x, resid);
		resid_g->SetName(Form("Fit Residuals - %i", i));
	  
	  	list->Add(resid_g);	

	}//for

}//fit_quadratic

//using linear calibration, estimates the location of the centroids
//and then refines that estimate
void find_centroids(TH1F *hist[], Double_t energy[], Double_t energy_er[], Double_t lin_gains[], Double_t lin_offsets[],
					Int_t start_pos, Int_t num_peaks_used, Double_t** centroids, Double_t** centroids_er) {

	for (int i = 0; i < num_cores; i++) {
		
		//if this linear equation is null, skip this core
		if (lin_gains[i] == -1 && lin_offsets[i] == -1) {
			for (int j = 0; j < num_peaks_used; j++) {
				centroids[i][j] = -1;
				centroids_er[i][j] = -1;
			}//for
		
			continue;
		}//if

		for (int j = 0; j < num_peaks_used; j++) {

			//using linear calibration
			//guess where this peak would be in this core's gamma ray spectrum
			Double_t x_guess = (energy[j] - lin_offsets[i]) / lin_gains[i];

			Int_t bin_guess = hist[i]->GetXaxis()->FindBin(x_guess);
			//Double_t y_guess = hist[i]->GetBinContent(bin_guess);
			hist[i]->SetAxisRange(x_guess-30, x_guess+30,"X");
			//to help with the peak fitting, move initial centroid guess to bin with most counts
            Int_t bin = hist[i]->GetMaximumBin();
            x_guess = hist[i]->GetXaxis()->GetBinCenter(bin);
			Double_t y_guess = hist[i]->GetBinContent(bin);			
//			x_guess = max_bin * hist[i]->GetMaximumBin() / num_bins;
			hist[i]->SetAxisRange(min_bin, max_bin,"X");

			TF1 *fit;
			fit = new TF1(Form("fit %i-%i", i, j),
//								   "[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))",
//                     				"[0]*(exp(-((x-[1])^2/(2*[2]^2))))*(1+ROOT::Math::erf([5]*((x-[1]))/([2]*pow(2,0.5)))) + [3] + [4]*x",
						"[0]*(exp(-((x-[1])^2/(2*[2]^2)))) + [3] + [4]*(0.5*(1-(ROOT::Math::erf(((x-[1])/([2]*2^(0.5)))))))",
 								   x_guess - 50, x_guess + 50);
			fit->SetParameters(y_guess, x_guess, 2, 14, 1);
			fit->SetParLimits(0, y_guess*0.8, y_guess*1.2); //area
			fit->SetParLimits(1, x_guess - 15, x_guess + 15); //centroid
			fit->SetParLimits(2, 0.2, 5); //sigma of gaussian distribution
	//		fit->SetParLimits(4, 0.1, 10000); //magnitude of step in background noise
	//		fit->SetParLimits(5, -10, -0.1); //background noise constant
			
			//fitting equation and saving centroids
			hist[i]->Fit(fit, "QR+");
			centroids[i][start_pos + j] = fit->GetParameter(1);
			centroids_er[i][start_pos + j] = fit->GetParError(1);

			cout << "Graph " << i << " Guess: " << x_guess << " Actual: " << centroids[i][start_pos + j] << " " << fit->GetChisquare()/fit->GetNDF() << endl;
			//cout << y_guess << "\t" << fit->GetParameter(0) << endl; 
		}//for
		
	}//for

}//find_centroids

void quad_energy() {	

	/*******Initialization of Source Energy Peak Data********/

	Double_t** source_energy = new Double_t*[num_known_sources];
	source_energy[0] = new Double_t[num_peaks[0]]{121.77,244.6976, 344.2785, 443.96, 778.9045, 964.082, 1112.080, 1408.013}; //eu152
	source_energy[1] = new Double_t[num_peaks[1]]{276.4, 302.851, 356.013, 383.848}; //ba133
	source_energy[2] = new Double_t[num_peaks[2]]{1771.3567, 2034.7097, 2598.500, 3253.5030}; //co56
	source_energy[3] = new Double_t[num_peaks[3]]{1173.240, 1332.508}; //co60
	
	Double_t** source_energy_er = new Double_t*[num_known_sources];
	source_energy_er[0] = new Double_t[num_peaks[0]]{0.0008,0.0008, 0.0012, 0.0012, 0.0024, 0.018, 0.003, 0.003}; //eu152
	source_energy_er[1] = new Double_t[num_peaks[1]]{0.0021, 0.0016, 0.0017, 0.0012}; //ba133
	source_energy_er[2] = new Double_t[num_peaks[2]]{0.0039, 0.0047, 0.004, 0.0044}; //co56
	source_energy_er[3] = new Double_t[num_peaks[3]]{0.003, 0.004}; //co60

	/*************************************************************/

	//number to store all the source IDs
	Int_t sources_used = 1;

	Int_t num_sources = 0; //number of sources used in this calibration
	Int_t num_peaks_used = 0; //number of peaks to process in this calibration
	
	TList *list = new TList;

	char use_quad;
	cout << "Use quadratic fit equation? (Y/N) ";
	cin >> use_quad;

	cin.clear();
    cin.ignore(INT_MAX, '\n');

	string source = "START";
	cout << "\nEnter the gamma sources used in calibration. Use the source's first two letters and mass number. (e.g. eu152, ba133).\nType 'END' when you've finished. Press <Enter> after each input.\n";
	while (source != "END" && source != "end") {

		Int_t sourceIdx = search_array(source_names, source, num_known_sources);

		if (sourceIdx != -1) {		
			if (sources_used % source_ids[sourceIdx] != 0) {
				sources_used *= source_ids[sourceIdx];
				num_peaks_used += num_peaks[sourceIdx];
				num_sources++;
			} else if (sources_used % source_ids[sourceIdx] == 0) {
				cout << "You already entered that source!" << endl;
			}//elseif
		} else if (source != "START") {
			cout << "\"" << source << "\"; That source isn't recognized!" << endl;
		}//else

		getline(cin, source);

	}//while

	if (num_sources == 0) {
		cout << "No sources inputted. Exiting program..." << endl;
		return;	
	}//if

	Double_t energy[num_peaks_used];
	Double_t energy_er[num_peaks_used];

	Double_t** centroids = new Double_t*[num_cores];
	Double_t** centroids_er = new Double_t*[num_cores];
	for (int i = 0; i < num_cores; i++) {
		centroids[i] = new Double_t[num_peaks_used];
		centroids_er[i] = new Double_t[num_peaks_used];
	}//for		

	TH1F *quad_hist[num_sources][num_cores];
	string sources_used_files[num_sources];
	Int_t count = 0;

	for (int i = 0; i < num_known_sources; i++) {
		
		if (sources_used % source_ids[i] == 0) {
			cout << "\nPlease enter the file path for the " << source_names[i] << " analysis tree.\n";
			getline(cin, sources_used_files[count++]);
	
			TFile * analysisfile = new TFile(&sources_used_files[count - 1][0], "READ");

			if (!analysisfile->IsOpen()) {
				i--;
				count--;
				continue;
			}//if
		}//if

	}//for

	Double_t lin_gains[num_cores];
	Double_t lin_offsets[num_cores];
					
	ifstream coeff_fr("lin_energy_coeff.txt");

	for (int i = 0; i < num_cores; i++) {
		coeff_fr >> lin_gains[i];
		coeff_fr >> lin_offsets[i];
	}//for	

	coeff_fr.close();	
			
	Int_t num_proc_peaks = 0;
	Int_t num_proc_sources = 0;
	for (int i = 0; i < num_known_sources; i++) {
		
		if (sources_used % source_ids[i] == 0) {
			string rootfile = sources_used_files[num_proc_sources];			

			move_into_array(energy, source_energy[i], num_proc_peaks, num_peaks[i]);
			move_into_array(energy_er, source_energy_er[i], num_proc_peaks, num_peaks[i]);

			char empty_gains_cal[70] = "EmptyGains.cal";

			load_histograms(&rootfile[0], empty_gains_cal, quad_hist[num_proc_sources], num_cores, i);

			for (int j = 0; j < num_cores; j++) {
				list->Add(quad_hist[num_proc_sources][j]);
			}//for

			find_centroids(quad_hist[num_proc_sources], source_energy[i], source_energy_er[i],
							lin_gains, lin_offsets, num_proc_peaks, num_peaks[i],
 							centroids, centroids_er);

			num_proc_peaks += num_peaks[i];
			num_proc_sources++;
		}//if

	}//for
	
	Double_t quad_thirds[num_cores];
	Double_t quad_gains[num_cores];
	Double_t quad_offsets[num_cores];

	fit_equation(list, energy, energy_er, centroids, centroids_er, num_peaks_used, quad_thirds, quad_gains, quad_offsets, use_quad);

	cout << "Writing fitted quadratic coefficients to 'quad_energy_coeff.txt'..." << endl;
	ofstream quad_out("quad_energy_coeff.txt");

	quad_out << "float non_lin[" << num_cores << "] = {" << quad_thirds[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_thirds[i];
	}//for

	quad_out << "};" << endl;
	
	quad_out << "float gain[" << num_cores << "] = {" << quad_gains[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_gains[i];
	}//for

	quad_out << "};" << endl;

	quad_out << "float offset[" << num_cores << "] = {" << quad_offsets[0];

	for (int i = 1; i < num_cores; i++) {
		quad_out << ", " << quad_offsets[i];
	}//for

	quad_out << "};" << endl;
	
	quad_out.close();

	TFile * outfile = new TFile("quad_cal.root", "RECREATE");
	outfile->cd();
	list->Write();
	outfile->Close();
}
