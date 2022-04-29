//g++ Sort.C -std=c++0x -I$GRSISYS/include -L$GRSISYS/lib `grsi-config --cflags --all-libs` -lTSuppressed `root-config --cflags --libs` -lTreePlayer -lMathMore -lSpectrum -lMinuit -lPyROOT -ltbb -o SortData

#define Sortcode_cxx

#include "Sort.h"

using namespace std;

void Sortcode::SortData(char const* infile, char const* calfile, char const* outfile, char const *beam, char const *target){

	Initialise();

	TCutG *kin_cutg;
	TCutG *target_cutg;
	//Krypton on Lead 

	if(strcmp(target,"196Pt")==0) {

		kin_cutg = new TCutG("Kr_S3_E_thetaCut",9);
		kin_cutg->SetPoint(0,17.8533,302311);
		kin_cutg->SetPoint(1,17.7174,219440);
		kin_cutg->SetPoint(2,30.7609,196608);
		kin_cutg->SetPoint(3,49.2391,155172);
		kin_cutg->SetPoint(4,49.2391,219440);
		kin_cutg->SetPoint(5,40.9511,278634);
		kin_cutg->SetPoint(6,23.8315,302311);
		kin_cutg->SetPoint(7,19.212,305694);
		kin_cutg->SetPoint(8,17.8533,302311);
		
		target_cutg = new TCutG("Target_E_theta",11);
		target_cutg->SetPoint(0,18.125,168702);
		target_cutg->SetPoint(1,17.3098,123884);
		target_cutg->SetPoint(2,27.7717,83294.1);
		target_cutg->SetPoint(3,39.4565,41858.5);
		target_cutg->SetPoint(4,45.0272,30019.7);
		target_cutg->SetPoint(5,49.6467,25791.6);
		target_cutg->SetPoint(6,49.1033,77374.8);
		target_cutg->SetPoint(7,40,121347);
		target_cutg->SetPoint(8,26.0054,168702);
		target_cutg->SetPoint(9,18.8043,174622);
		target_cutg->SetPoint(10,18.125,168702);
	
	} else if(strcmp(target,"208Pb")==0) {

		kin_cutg = new TCutG("Kr_S3_E_thetaCut",9);
		kin_cutg->SetPoint(0,17.8533,302311);
		kin_cutg->SetPoint(1,17.7174,219440);
		kin_cutg->SetPoint(2,30.7609,196608);
		kin_cutg->SetPoint(3,49.2391,155172);
		kin_cutg->SetPoint(4,49.2391,219440);
		kin_cutg->SetPoint(5,40.9511,278634);
		kin_cutg->SetPoint(6,23.8315,302311);
		kin_cutg->SetPoint(7,19.212,305694);
		kin_cutg->SetPoint(8,17.8533,302311);

		target_cutg = new TCutG("Target_E_theta",11);
		target_cutg->SetPoint(0,18.125,168702);
		target_cutg->SetPoint(1,17.3098,123884);
		target_cutg->SetPoint(2,27.7717,83294.1);
		target_cutg->SetPoint(3,39.4565,41858.5);
		target_cutg->SetPoint(4,45.0272,30019.7);
		target_cutg->SetPoint(5,49.6467,25791.6);
		target_cutg->SetPoint(6,49.1033,77374.8);
		target_cutg->SetPoint(7,40,121347);
		target_cutg->SetPoint(8,26.0054,168702);
		target_cutg->SetPoint(9,18.8043,174622);
		target_cutg->SetPoint(10,18.125,168702);

	}

	TCutG *S3Tigtime = new TCutG("TIG_S3_time",13);
	S3Tigtime->SetPoint(0,-57.8886,69.3677);
	S3Tigtime->SetPoint(1,57.4722,69.3677);
	S3Tigtime->SetPoint(2,94.4438,82.5806);
	S3Tigtime->SetPoint(3,100.996,122.219);
	S3Tigtime->SetPoint(4,94.2098,399.69);
	S3Tigtime->SetPoint(5,69.8741,485.574);
	S3Tigtime->SetPoint(6,51.6223,597.884);
	S3Tigtime->SetPoint(7,44.6024,809.29);
	S3Tigtime->SetPoint(8,29.6265,3986.99);
	S3Tigtime->SetPoint(9,-27.0009,4013.42);
	S3Tigtime->SetPoint(10,-48.9967,3808.62);
	S3Tigtime->SetPoint(11,-57.4206,1833.29);
	S3Tigtime->SetPoint(12,-57.8886,69.3677);

	double tgatemin = -80; // Note that this has a crude timegate, not a TCutG, so check these numbers - JH
	double tgatemax = 400;
//	double betarough = 0.11;
	double betarough = 0.05;
	double s3_x_offset = 0;
	double s3_y_offset = 0;
	double s3_phi_offset = -25.7*TMath::Pi()/180;

	TFile* inputfile = new TFile(infile,"READ");
	if(!inputfile->IsOpen()){
	  printf("Opening file failed, aborting\n");
	  return;
	}
	TChain* AnalysisTree = (TChain*)inputfile->Get("AnalysisTree");
	TTree* tree = (TTree*)AnalysisTree->GetTree();
	printf("Reading calibration file: %s\n",calfile);
	TChannel::ReadCalFile(calfile);
	TGRSIRunInfo *runinfo = new TGRSIRunInfo;
//	runinfo->SetHPGeArrayPosition(145.); // We're running at 145. - JH

	TTigress::SetArrayBackPos();//We're running at 145. - JS

//	TS3 *s3 = 0;
	TRF *rf = 0;
	TTigress *tigress = 0;	

	AnalysisTree->SetBranchAddress("TTigress",&tigress);
//	AnalysisTree->SetBranchAddress("TS3",&s3);
//	AnalysisTree->SetBranchAddress("TRF",&rf);

	long nentries = AnalysisTree->GetEntries();

	TTigressHit *tigress_hit, *tigress_hit2;
//	TS3Hit *s3hit, *ringhit, *sectorhit;
	TTigressHit *addback_hit, *addback_hit2;
	TGRSIDetectorHit *segment_hit;

	double mBeam;
	if(strcmp(beam,"80Kr")!=0)
		mBeam = 79.916379;
	else if(strcmp(beam,"82Sr")!=0)
		mBeam = 81.918402;
	else if(strcmp(beam,"80Sr")!=0)
		mBeam = 79.924521;
	else{
		printf("Beam not recognized, assuming 82Sr\n");
		mBeam = 81.918402;
	}
	mBeam = 81.918402;
	printf("Beam mass: %f\n",mBeam);	

	double EBeam;// = 4.17 * mBeam;	
	if(strcmp(beam,"80Kr")!=0 || strcmp(beam,"82Sr")!=0)
		EBeam = 4.17 * mBeam;
	else if(strcmp(beam,"80Sr")!=0)
		EBeam = 4.27 * mBeam;
	else{
		printf("Beam not recognized, assuming 82Sr energy\n");
		EBeam = 4.17 * mBeam;
	}
	printf("Beam energy: %f MeV\t%fMeV/u\n",EBeam,EBeam/mBeam);
	TReaction *Sr82 = new TReaction(beam,target,beam,target,EBeam,0,true);            	//Added by Stephen
	TReaction *Sr82_1stex = new TReaction(beam,target,beam,target,EBeam,0.616,true);  	//Added by Stephen
	TReaction *Sr82_2ndex = new TReaction(beam,target,beam,target,EBeam,1.256,true);  	//Added by Stephen
	TReaction *Sr82_3rdex = new TReaction(beam,target,beam,target,EBeam,1.436,true);  	//Added by Stephen
	TReaction *Pt196 = new TReaction(beam,target,beam,target,EBeam,0.355);            	//Added by Stephen


	int tig_ring[16][4] = {
		{ 0, 1, 1, 0 }, { 0, 1, 1, 0 }, { 0, 1, 1, 0 }, { 0, 1, 1, 0 },
		{ 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 },
		{ 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 }, { 2, 3, 3, 2 },
		{ 4, 5, 5, 4 }, { 4, 5, 5, 4 }, { 4, 5, 5, 4 }, { 4, 5, 5, 4 }
	};

	TVector3 reac_vec;
	TVector3 recoil_vec;
	TVector3 ejectile_vec;
	TVector3 *corrected_s3_vec;
	TVector3 s3pos;


	double M_82Sr = 81.918402 * 931494.0954; // Not used
	double M_196Pt = 195.9649515 * 931494.0954; // Not used
	double M_208Pb = 207.9766521 * 931494.0954; // Not used
	TGraph *kinline = Sr82->KinVsTheta(0,180,2);
	TGraph *kinline_Sr82_1 = Sr82_1stex->KinVsTheta(0, 180, 1, 2);  //Added by Stephen
	TGraph *kinline_Sr82_2 = Sr82_2ndex->KinVsTheta(0, 180, 1, 2);  //Added by Stephen
	TGraph *kinline_Sr82_3 = Sr82_3rdex->KinVsTheta(0, 180, 1, 2);  //Added by Stephen

	double E_kin[180], theta_kin[180], beta[180], beta2[180], beta3[180];
	double tempEkin, tempExc, tempBeta, tempThetalab, tempThetaCM, tempBeta_Tar;

	double innerradius=11.;
	double outerradius=35.;
	double targetdist =31.;
//	double targetdist =14.5;
	double strippitch = (outerradius-innerradius)/24.;
	double ringsize[24], ringnum[24];
	for(int i=0;i<24;i++){
		double innertheta = TMath::RadToDeg()*TMath::ATan((i*strippitch + innerradius)/targetdist);
		double outertheta = TMath::RadToDeg()*TMath::ATan(((i+1)*strippitch + innerradius)/targetdist); 
		ringsize[i] = outertheta - innertheta;
		ringnum[i] = i;
	}
	TGraph *ringsizespline = new TGraph(24,ringnum,ringsize);
	ringsizespline->SetName("RingSize");

	Long_t ts_start = 0;

	printf("Begin sort\n");

	TRandom3* rand = new TRandom3(0);

	int DetSet[24] = {0,0,0,0,0,0,1,1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3};


	for(long jentry=0;jentry<nentries;jentry++){
		
//		double rftime; // I don't think we're taking RF (there's no need to) so this can be removed - JH

		AnalysisTree->GetEntry(jentry);

//		s3->SetMultiHit();

//		rf_time->Fill(rf->Time()); // This'll confuse it if there's no RF - JH
//		rftime=rf->Time()/10; // Ditto - JH

//		tig_mult->Fill(tigress->GetMultiplicity());
		tigress->SetArrayBackPos();
		tigress->ResetAddback(); // Don't know if the TIGRESS addback scheme has changed, leave for now, could cause problems - JH
		//tigress->SetSegmentHits(false);

		// Raw histograms
		for(int i=0;i<tigress->GetMultiplicity();i++){

			tigress_hit = tigress->GetTigressHit(i);
			if(tigress_hit->GetEnergy()<10)
				continue;
//			double tigtime = (tigress_hit->GetTimeStamp()<<4 & 0x07ffffff) - tigress_hit->GetCfd(); //This might no longer work (not hugely important) - JH
//			tig_time->Fill(((tigress_hit->GetTimeStamp()<<4 & 0x07ffffff) - tigress_hit->GetCfd())); // Ditto - JH
			tig_E->Fill(tigress_hit->GetEnergy());
//			tig_E_dopp->Fill(tigress_hit->GetDoppler(betarough)); 
//			for(int j=0;j<s3->GetPixelMultiplicity();j++){ 
//				s3hit = s3->GetS3Hit(j); 
//				double s3time = (s3hit->GetTimeStamp()<<4 & 0x07ffffff) - s3hit->GetCfd(); // Again - JH
//				tig_s3_time->Fill(tigtime,s3time);
//				s3_time->Fill(s3time);
//				tig_s3_dCfd->Fill(tigress_hit->GetCfd()-s3hit->GetCfd());
//				tigE_dCfd->Fill(tigress_hit->GetCfd()-s3hit->GetCfd(),tigress_hit->GetEnergy());
//			}
		}
		for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
			addback_hit = tigress->GetAddbackHit(i);
			if(addback_hit->GetEnergy()<10)
				continue;
			addback_E->Fill(addback_hit->GetEnergy());
			addback_E_channel->Fill(addback_hit->GetEnergy(),addback_hit->GetArrayNumber()); // Nikita
//			addback_E_dopp->Fill(tigress_hit->GetDoppler(betarough));
		}
/*niki		for(int i=0;i<s3->GetRingMultiplicity();i++){
			ringhit = s3->GetRingHit(i); 
			if(ringhit->GetIsDownstream())
				downstream_rings->Fill(ringhit->GetRing(),ringhit->GetEnergy());
			else
				upstream_rings->Fill(ringhit->GetRing(),ringhit->GetEnergy());
			for(int j=0;j<s3->GetSectorMultiplicity();j++){
				sectorhit = s3->GetSectorHit(j);
				s3_dCfd->Fill(ringhit->GetCfd()-sectorhit->GetCfd()); 
				front_back_energy->Fill(ringhit->GetEnergy(),sectorhit->GetEnergy());
			}
		}
		for(int i=0;i<s3->GetSectorMultiplicity();i++){
			sectorhit = s3->GetSectorHit(i);
			if(sectorhit->GetIsDownstream())
				downstream_sectors->Fill(sectorhit->GetSector(),sectorhit->GetEnergy());
			else
				upstream_sectors->Fill(sectorhit->GetSector(),sectorhit->GetEnergy());
			
		}
		for(int i=0;i<s3->GetPixelMultiplicity();i++){
			s3hit = s3->GetS3Hit(i);
			s3_E->Fill(s3hit->GetEnergy());
			s3_E_theta->Fill(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy());
			if(s3hit->GetIsDownstream() && kin_cutg->IsInside(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy())){
				s3pos = s3hit->GetPosition(s3_phi_offset,true);			
				TVector3 iter, orig;
				orig = s3pos;
				/*for(int j=0;j<21;j++){
					for(int k=0;k<21;k++){
						iter = s3pos;
						iter.SetX(orig.X()+j*0.5-5);		
						iter.SetY(orig.Y()+k*0.5-5);
						thetaphi_offset_iter[j][k]->Fill(iter.Theta()*TMath::RadToDeg(),iter.Phi()*TMath::RadToDeg());
					}
				}*/
/*niki				hitmap_down->Fill(s3pos.X(),s3pos.Y());
				theta_phi_hitmap->Fill(s3pos.Theta()*TMath::RadToDeg(),s3pos.Phi()*TMath::RadToDeg());
				s3pos.SetX(s3pos.X()+s3_x_offset);
				s3pos.SetY(s3pos.Y()+s3_y_offset);
				theta_phi_hitmap_corr->Fill(s3pos.Theta()*TMath::RadToDeg(),s3pos.Phi()*TMath::RadToDeg());
				theta_phi_hitmap_corr_ring[s3hit->GetRing()]->Fill(s3pos.Theta()*TMath::RadToDeg(),s3pos.Phi()*TMath::RadToDeg());
				theta_phi_hitmap_corr_ring_group[DetSet[s3hit->GetRing()]]->Fill(s3pos.Theta()*TMath::RadToDeg(),s3pos.Phi()*TMath::RadToDeg());

			}
			else{
				s3pos = s3hit->GetPosition(s3_phi_offset,true);					
				hitmap_up->Fill(s3pos.X(),s3pos.Y());

			}
		}

		// Time gated spectra
		bool timegated = false;
		for(int j=0;j<s3->GetPixelMultiplicity();j++){
			s3hit = s3->GetS3Hit(j);
			if(s3hit->GetEnergy()<500)
				continue;
			double s3time = (s3hit->GetTimeStamp()<<4 & 0x07ffffff) - s3hit->GetCfd(); // Again - JH
			for(int i=0;i<tigress->GetAddbackMultiplicity();i++){
				addback_hit = tigress->GetAddbackHit(i);
				if(addback_hit->GetEnergy()<10)
					continue;

//				if(addback_hit->GetCfd()-s3hit->GetCfd()>tgatemin && addback_hit->GetCfd()-s3hit->GetCfd()<tgatemax){ // Time gated	

			double addtime = (addback_hit->GetTimeStamp()<<4 & 0x07ffffff) - addback_hit->GetCfd(); //This might no longer
				if(S3Tigtime->IsInside(addback_hit->GetCfd()-s3hit->GetCfd(),addback_hit->GetEnergy())) { //Time Gated 
					timegated = true;
					s3_E_time_gated->Fill(s3hit->GetEnergy());
					s3_E_theta_time_gated->Fill(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy());	
					s3_E_tig_E->Fill(s3hit->GetEnergy(),addback_hit->GetEnergy());			
					s3_E_dopp_E->Fill(s3hit->GetEnergy(),addback_hit->GetDoppler(betarough));	

					// You can detected either beam-like or target-like
					// and need to Doppler correct for detected beam-like
					// to target like and vice-versa.
					// This code (should) do that - JH

					reac_vec = (TVector3)s3hit->GetPosition(s3_phi_offset,true);                 //Added By Stephen
					recoil_vec.SetMagThetaPhi(1.,Sr82->ConvertThetaLab(reac_vec.Theta(),2,3),reac_vec.Phi()+TMath::Pi()); // Convert detected Sr to Pb/Pt - JH
					ejectile_vec.SetMagThetaPhi(1.,Sr82->ConvertThetaLab(reac_vec.Theta(),3,2),reac_vec.Phi()+TMath::Pi()); // Convert detector Pb/Pt to Sr - JH
					tig_dopp_s3_corr->Fill(addback_hit->GetDoppler(betarough,&reac_vec));   //Added By Stephen
                   			tempThetalab = reac_vec.Theta();   					//Added By Stephen - This shouldn't be necessary, you've already defined theta from your S3
                  			tempEkin = s3hit->GetEnergy();						//Added By Stephen

					if(addback_hit->GetDoppler(tempBeta,&reac_vec) > 545 && addback_hit->GetDoppler(tempBeta,&reac_vec) < 605) 
						E_Theta_Gamma_Gate->Fill(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy());
					s3_theta->Fill(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg());
					if(s3hit->GetIsDownstream() && kin_cutg->IsInside(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy()) ) { //Beam Kinematic Gate on Downstream  -- SG Added TCUTG
						tempBeta = Sr82->AnalysisBetaFromThetaLab(reac_vec.Theta(),2); // Beta for Sr, gated on Sr
						tempBeta_Tar = Sr82->AnalysisBetaFromThetaLab(recoil_vec.Theta(),3); // Beta for Target, gated on Sr
						tig_dopp_beam_beamcut->Fill(addback_hit->GetDoppler(tempBeta,&reac_vec));
						tig_dopp_target_beamcut->Fill(addback_hit->GetDoppler(tempBeta_Tar,&recoil_vec));
						tigE_angle_beamcut->Fill(addback_hit->GetEnergy(),reac_vec.Angle(addback_hit->GetPosition())*TMath::RadToDeg());
						tigdopp_angle_beamcut->Fill(addback_hit->GetDoppler(tempBeta,&reac_vec),reac_vec.Angle(addback_hit->GetPosition())*TMath::RadToDeg());
						tig_beam_beam_DetSet[DetSet[s3hit->GetRing()]]->Fill(addback_hit->GetDoppler(tempBeta,&reac_vec));
						tig_beam_target_DetSet[DetSet[s3hit->GetRing()]]->Fill(addback_hit->GetDoppler(tempBeta_Tar,&recoil_vec));
						for(int k=i+1;k<tigress->GetAddbackMultiplicity();k++){
							addback_hit2 = tigress->GetAddbackHit(k);
							double addtime2 = (addback_hit2->GetTimeStamp()<<4 & 0x07ffffff) - addback_hit2->GetCfd();
							if(S3Tigtime->IsInside(addback_hit->GetCfd()-s3hit->GetCfd(),addback_hit->GetEnergy())) {
								tig_dopp_gg_beam_beamcut->Fill(addback_hit->GetDoppler(tempBeta,&reac_vec),addback_hit2->GetDoppler(tempBeta,&reac_vec));
								tig_dopp_gg_beam_beamcut->Fill(addback_hit2->GetDoppler(tempBeta,&reac_vec),addback_hit->GetDoppler(tempBeta,&reac_vec));
								tig_dopp_gg_target_beamcut->Fill(addback_hit->GetDoppler(tempBeta_Tar,&recoil_vec),addback_hit2->GetDoppler(tempBeta_Tar,&recoil_vec));
								tig_dopp_gg_target_beamcut->Fill(addback_hit2->GetDoppler(tempBeta_Tar,&recoil_vec),addback_hit->GetDoppler(tempBeta_Tar,&recoil_vec));
							}
						}
					}
					else if(s3hit->GetIsDownstream() && target_cutg->IsInside(s3hit->GetPosition(s3_phi_offset,true).Theta()*TMath::RadToDeg(),s3hit->GetEnergy()) ) { //Target Kinematic Gate on Downstream -- SG Added TCUTG
						tempBeta = Sr82->AnalysisBetaFromThetaLab(ejectile_vec.Theta(),2); // Beta for Sr, gated on Target
						tempBeta_Tar = Sr82->AnalysisBetaFromThetaLab(reac_vec.Theta(),3); // Beta for target, gated on Target
						tig_dopp_beam_targetcut->Fill(addback_hit->GetDoppler(tempBeta,&ejectile_vec));
						tig_dopp_target_targetcut->Fill(addback_hit->GetDoppler(tempBeta_Tar,&reac_vec));
						tigE_angle_targetcut->Fill(addback_hit->GetEnergy(),ejectile_vec.Angle(addback_hit->GetPosition())*TMath::RadToDeg());
						tigdopp_angle_targetcut->Fill(addback_hit->GetDoppler(tempBeta,&ejectile_vec),ejectile_vec.Angle(addback_hit->GetPosition())*TMath::RadToDeg());
						tig_target_beam_DetSet[DetSet[s3hit->GetRing()]]->Fill(addback_hit->GetDoppler(tempBeta,&ejectile_vec));
						tig_target_target_DetSet[DetSet[s3hit->GetRing()]]->Fill(addback_hit->GetDoppler(tempBeta_Tar,&reac_vec));
						for(int k=i+1;k<tigress->GetAddbackMultiplicity();k++){
							addback_hit2 = tigress->GetAddbackHit(k);
							double addtime2 = (addback_hit2->GetTimeStamp()<<4 & 0x07ffffff) - addback_hit2->GetCfd();
							if(S3Tigtime->IsInside(addback_hit->GetCfd()-s3hit->GetCfd(),addback_hit->GetEnergy())) {
								tig_dopp_gg_beam_targetcut->Fill(addback_hit->GetDoppler(tempBeta,&ejectile_vec),addback_hit2->GetDoppler(tempBeta,&ejectile_vec));
								tig_dopp_gg_beam_targetcut->Fill(addback_hit2->GetDoppler(tempBeta,&ejectile_vec),addback_hit->GetDoppler(tempBeta,&ejectile_vec));
								tig_dopp_gg_target_targetcut->Fill(addback_hit->GetDoppler(tempBeta_Tar,&reac_vec),addback_hit2->GetDoppler(tempBeta_Tar,&reac_vec));
								tig_dopp_gg_target_targetcut->Fill(addback_hit2->GetDoppler(tempBeta_Tar,&reac_vec),addback_hit->GetDoppler(tempBeta_Tar,&reac_vec));
							}
						}

					}
					if(!s3hit->GetIsDownstream() && s3hit->GetEnergy()>20000) { //Beam Kinematic Gated Upstream S3 -- SG Added TCUTG
						tempBeta = Sr82->AnalysisBetaFromThetaLab(reac_vec.Theta(),2); // Beta for Sr, gated on Sr
						tempBeta_Tar = Sr82->AnalysisBetaFromThetaLab(recoil_vec.Theta(),3); // Beta for Target, gated on Sr
						tig_dopp_beam_upstream->Fill(addback_hit->GetDoppler(tempBeta,&reac_vec));
						tig_dopp_target_upstream->Fill(addback_hit->GetDoppler(tempBeta_Tar,&recoil_vec));
					}
				}
			}
			if(timegated){
				tig_E_time_gated->Fill(addback_hit->GetEnergy());
				doppler_E_time_gated->Fill(addback_hit->GetDoppler(betarough));
			}
		}	
		for(int i=0;i<tigress->GetMultiplicity();i++){
			tigress_hit = tigress->GetTigressHit(i);
			if(tigress_hit->GetEnergy()<10)
				continue;
			for(int j=0;j<s3->GetPixelMultiplicity();j++){
				s3hit = s3->GetS3Hit(j);
				if(tigress_hit->GetCfd()-s3hit->GetCfd()>-80 && tigress_hit->GetCfd()-s3hit->GetCfd()<100){ // Time gated	
					tig_dopp_s3_corr_noab->Fill(tigress_hit->GetDoppler(betarough,&reac_vec));
					tig_dopp_kin_corr_noab->Fill(tigress_hit->GetDoppler(tempBeta,&reac_vec));
				}
			}
		}	
niki */	
	  	if(jentry%10000 == 0) 
		    	cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry/nentries << "% complete" << "\r" << flush;    

	}

	cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";	

	cout << "Writing histograms to " << outfile << endl;

	TFile *myfile = new TFile(outfile,"RECREATE");
	TDirectory *rawdir = myfile->mkdir("RawHistograms");
	rawdir->cd();
	rawlist->Write();
	myfile->cd();
/*	TDirectory* timedir = myfile->mkdir("TimeGatedHistograms");
	timedir->cd();
	timelist->Write();
	myfile->cd();
	TDirectory* gateddir = myfile->mkdir("GatedHistograms");
	gateddir->cd();
	gatedlist->Write();
	myfile->cd();
	TDirectory* kindir = myfile->mkdir("Kinematics");
	kindir->cd();
	kinline->SetName("KinLine");
	kin_cutg->Write();
	S3Tigtime->Write();
	target_cutg->Write(); // Write cuts to file, useful later - JH
	kinline->Write();
	myfile->cd();
	TDirectory *s3_dir = myfile->mkdir("S3Spectra");
	s3_dir->cd();
	s3list->Write();
	//TDirectory *thetaphi_dir = s3_dir->mkdir("ThetaPhi_Ring");
	//thetaphi_dir->cd();
	//thetaphilist->Write();
	//s3_dir->cd();
	myfile->cd();
	TDirectory *combo_dir = myfile->mkdir("Combined_Rings");
	combo_dir->cd();	
	combolist->Write();
	myfile->cd();
	TDirectory *iter_dir = myfile->mkdir("ThetaPhi_Iterations");
	iter_dir->cd();
	thetaphi_iter->Write();
*/	myfile->cd();
	myfile->Close();

	cout << "Histograms written, sorting complete" << endl;

}

void Sortcode::SortFrag(char const* infile, char const* calfile, char const* outfile){

	cout << "Sorting fragments..." << endl;

	InitialiseFragments();

	cout << "Histograms created, opening " << infile << endl;

	TFile* inputfile = new TFile(infile,"READ");
	if(!inputfile->IsOpen()){
	  printf("Opening file failed, aborting\n");
	  return;
	}
	TChain* FragmentTree = (TChain*)inputfile->Get("FragmentTree");
	TTree* tree = (TTree*)FragmentTree->GetTree();
	printf("Reading calibration file: %s\n",calfile);
	TChannel::ReadCalFile(calfile);

	long nentries = FragmentTree->GetEntries();

	TFragment *frag = 0;	
	FragmentTree->SetBranchAddress("TFragment",&frag);

	cout << nentries << " entries in FragmentTree" << endl;

	for(long jentry=0;jentry<nentries;jentry++){

		FragmentTree->GetEntry(jentry);

		charge->Fill(frag->GetCharge()/20.);
	
		channel_charge->Fill(frag->GetChannel()->GetNumber(),frag->GetCharge()/20.);

	  	if(jentry%10000 == 0) 
		    	cout << setiosflags(ios::fixed) << "Entry " << jentry << " of " << nentries << ", " << 100 * jentry/nentries << "% complete" << "\r" << flush;    

	}

	cout << "Entry " << nentries << " of " << nentries << ", 100% Complete!\n";	


	cout << "Writing histograms to " << outfile << endl;

	TFile *outputfile = new TFile(outfile,"RECREATE");
	fraglist->Write();
	outputfile->Close();

}

int main(int argc, char** argv){

	Sortcode *mysort = new Sortcode();

	char const *infile;
	char const *outfile;
	char const *calfile;
	char const *fragfile;
	char const *beam;
	char const *target;
	printf("Starting sortcode\n");

	bool SortFragments = false;

	for(int i=0;i<argc;i++)
		printf("%s\t",argv[i]);
	printf("\n");

	// Input-chain-file, output-histogram-file
	if(argc == 1)
	{
		infile = "Data.root";
		calfile = "CalibrationFile.cal";
		outfile = "Histograms.root";
		beam = "82Sr";	
		target = "208Pb";
	}
	else if(argc == 2)
	{
		infile = argv[1];
		calfile = "CalibrationFile.cal";
		outfile = "Histograms.root";
		beam = "82Sr";	
		target = "208Pb";
	}
	else if(argc == 3)
	{
		infile = argv[1];
		calfile = argv[2];
		outfile = "Histograms.root";
		beam = "82Sr";	
		target = "208Pb";
	}
	else if(argc == 4)
	{
		infile = argv[1];
		calfile = argv[2];
		outfile = argv[3];
		beam = "82Sr";	
		target = "208Pb";
	}
	else if(argc == 5)
	{
		infile = argv[1];
		calfile = argv[2];
		outfile = argv[3];
		beam = argv[4];	
		target = "208Pb";
	}
	else if(argc == 6)
	{
		infile = argv[1];
		calfile = argv[2];
		outfile = argv[3];
		beam = argv[4];	
		target = argv[5];
	}
	else if(argc == 7)
	{
		infile = argv[1];
		calfile = argv[2];
		outfile = argv[3];
		beam = argv[4];	
		target = argv[5];
		fragfile = argv[6];
		SortFragments = true;
	}
	else if(argc > 7)
	{
		printf("Doh! Too many arguments\n");
		return 0;
	}

	printf("Input file: %s\nCalibration file: %s\nOutput file: %s\n",infile,calfile,outfile);
	printf("Beam: %s\nTarget: %s\n",beam,target);

	if(SortFragments){
		char tmpfraghist[128];
		sprintf(tmpfraghist,"Fragments_%s",outfile);
		const char *fraghist = tmpfraghist;
		mysort->SortFrag(fragfile,calfile,fraghist);	
	}	

	mysort->SortData(infile,calfile,outfile,beam,target);
	
	return 0;

}
