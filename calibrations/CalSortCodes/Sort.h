// ADD NEW HISTOGRAMS IN Declarations.h !!!!!!!!!!!!
#ifndef Sortcode_h
#define Sortcode_h

#include <iostream>
#include "TCutG.h"
#include <iomanip>
#include "TH1.h"
#include "TTree.h"
#include "TChain.h"
#include "TH2.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TList.h"
#include "TTip.h"
#include "TTigress.h"
#include "TS3.h"
#include "TS3Hit.h"
#include "TRF.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "Declarations.h"
#include "TGRSIRunInfo.h"
#include "TKinematics.h"
#include "TReaction.h"
#include "TNucleus.h"
#include "TPeak.h"
#include "TRandom3.h"
#include "TFragment.h"

using namespace std;

class Sortcode {

	public :

		Sortcode(){;} 
		void SortData(const char*, const char*, const char*, const char*, const char*);
		void Initialise();
		void SortFrag(const char*, const char*, const char*);	
		void InitialiseFragments();

};
#endif

#ifdef Sortcode_cxx

void Sortcode::Initialise(){

	printf("Start initializations\n");
  	printf("Creating list\n");
	rawlist = new TList; timelist = new TList; gatedlist = new TList; ringlist = new TList; pbgated = new TList; ringlist2 = new TList; pb_ringlist = new TList; pb_ringlist2 = new TList; laserlist = new TList; s3list = new TList; combolist = new TList; thetaphilist = new TList;
	printf("Creating histograms\n");

	// Raw histos:
	rf_time = new TH1F("RF_Time","RF_Time",128,0,128); rawlist->Add(rf_time);	
	tig_mult = new TH1F("TIG_multiplicity","TIG_multiplicity",10,0,10); rawlist->Add(tig_mult);
	tig_E = new TH1F("TIG_E","TIG_E",8192,0,8192); rawlist->Add(tig_E);
	addback_E = new TH1F("Addback_E","Addback_E",8192,0,8192); rawlist->Add(addback_E);
	addback_E_channel = new TH2F("Addback_E_channel","Addback_E_channel",8192,0,8192,64,0,64); rawlist->Add(addback_E_channel);
/*niki	tig_time = new TH1F("tig_time","tig_time",1024,0,8192); rawlist->Add(tig_time);
	s3_time = new TH1F("s3_time","s3_time",1024,0,8192); rawlist->Add(s3_time);
	s3_dCfd = new TH1F("s3_dCfd","s3_dCfd",1024,-1024,1024); rawlist->Add(s3_dCfd);
	front_back_energy = new TH2F("Front_Back_Energy","Front_Back_Energy",2048,0,524288,2048,0,524288); rawlist->Add(front_back_energy);
	s3_E = new TH1F("S3_E","S3_E",8192,0,65536); rawlist->Add(s3_E);
	tig_s3_time = new TH2F("Tig_S3_time","Tig_S3_time",1024,0,8192,1024,0,8192); rawlist->Add(tig_s3_time);
	tig_s3_dCfd = new TH1F("Tig_S3_dCfd","tig_s3_dCfd",1024,-1024,1024); rawlist->Add(tig_s3_dCfd);
	tigE_dCfd = new TH2F("Tig_E_dCfd","Tig_E_dCfd",1024,-1024,1024,2048,0,4096); rawlist->Add(tigE_dCfd);
	s3_E_theta = new TH2F("S3_E_Theta","S3_E_Theta",180,0,180,4096,0,524288); rawlist->Add(s3_E_theta);
	tig_E_dopp = new TH1F("Tigress_E_Doppler_Corrected","Tigress_E_Doppler_Corrected",8192,0,8192); rawlist->Add(tig_E_dopp);
	addback_E_dopp = new TH1F("Addback_E_Doppler_Corrected","Addback_E_Doppler_Corrected",8192,0,8192); rawlist->Add(addback_E_dopp);
	s3_ring_charge = new TH2F("S3_Ring_vs_Charge","S3_Ring_vs_Charge",24,0,24,4096,0,524288); rawlist->Add(s3_ring_charge);

	// Time gated histos:
	s3_E_time_gated = new TH1F("S3_E_time","S3_E_time",8192,0,524288); timelist->Add(s3_E_time_gated);
	s3_E_theta_time_gated = new TH2F("S3_E_Theta_time","S3_E_Theta_time",180,0,180,4096,0,524288); timelist->Add(s3_E_theta_time_gated);
	s3_E_tig_E = new TH2F("S3_E_Tigress_E","S3_E_Tigress_E",2048,0,524288,2048,0,2048); timelist->Add(s3_E_tig_E);
	s3_E_dopp_E = new TH2F("S3_E_Dopp_E","S3_E_Dopp_E",2048,0,524288,2048,0,2048); timelist->Add(s3_E_dopp_E);
	tig_E_time_gated = new TH1F("tig_E_time_gated","tig_E_time_gated",8192,0,8192); timelist->Add(tig_E_time_gated);
	doppler_E_time_gated = new TH1F("doppler_E_time_gated","doppler_E_time_gated",8192,0,8192); timelist->Add(doppler_E_time_gated);
	tig_dopp_s3_corr = new TH1F("TIGRESS_Doppler_S3_Corrected","TIGRESS_Doppler_S3_Corrected",8192,0,8192); timelist->Add(tig_dopp_s3_corr);	
	tig_dopp_kin_corr = new TH1F("TIGRESS_Doppler_kinematics_corrected","TIGRESS_Doppler_kinematics_corrected",8192,0,8192); timelist->Add(tig_dopp_kin_corr);
	s3_theta = new TH1F("s3_theta","s3_theta",180,0,180); timelist->Add(s3_theta);
	tig_dopp_s3_corr_noab = new TH1F("TIGRESS_Doppler_S3_Corrected_noab","TIGRESS_Doppler_S3_Corrected_noab",8192,0,8192); timelist->Add(tig_dopp_s3_corr_noab);	
	tig_dopp_kin_corr_noab = new TH1F("TIGRESS_Doppler_kinematics_corrected_noab","TIGRESS_Doppler_kinematics_corrected_noab",8192,0,8192); timelist->Add(tig_dopp_kin_corr_noab);

	tigE_angle_beamcut = new TH2F("TIGRESS_E_vs_angle_BeamCut","TIGRESS_E_vs_angle_BeamCut",2048,0,2048,180,0,180); gatedlist->Add(tigE_angle_beamcut);
	tigE_angle_targetcut = new TH2F("TIGRESS_E_vs_angle_TargetCut","TIGRESS_E_vs_angle_TargetCut",2048,0,2048,180,0,180); gatedlist->Add(tigE_angle_targetcut);
	tigdopp_angle_beamcut = new TH2F("TIGRESS_Doppler_vs_angle_BeamCut","TIGRESS_Doppler_vs_angle_BeamCut",2048,0,2048,180,0,180); gatedlist->Add(tigdopp_angle_beamcut);
	tigdopp_angle_targetcut = new TH2F("TIGRESS_Doppler_vs_angle_TargetCut","TIGRESS_Doppler_vs_angle_TargetCut",2048,0,2048,180,0,180); gatedlist->Add(tigdopp_angle_targetcut);

	tig_dopp_beam_beamcut = new TH1F("TIGRESS_E_Beam_Cut_Beam_Doppler","TIGRESS_E_Beam_Cut_Beam_Doppler",2048,0,4096); gatedlist->Add(tig_dopp_beam_beamcut);
	tig_dopp_target_beamcut = new TH1F("TIGRESS_E_Beam_Cut_Target_Doppler","TIGRESS_E_Beam_Cut_Target_Doppler",2048,0,4096); gatedlist->Add(tig_dopp_target_beamcut);
	tig_dopp_beam_upstream = new TH1F("Upstream_Beam_Doppler_Corrected","Upstream_Beam_Doppler_Corrected",8192,0,8192); gatedlist->Add(tig_dopp_beam_upstream);
	tig_dopp_target_upstream = new TH1F("Upstream_Target_Doppler_Corrected","Upstream_Target_Doppler_Corrected",8192,0,8192); gatedlist->Add(tig_dopp_target_upstream);
	tig_dopp_beam_targetcut = new TH1F("TIGRESS_E_Target_Cut_Beam_Doppler","TIGRESS_E_Target_Cut_Beam_Doppler",2048,0,4096); gatedlist->Add(tig_dopp_beam_targetcut);
	tig_dopp_target_targetcut = new TH1F("TIGRESS_E_Target_Cut_Target_Doppler","TIGRESS_E_Target_Cut_Target_Doppler",2048,0,4096); gatedlist->Add(tig_dopp_target_targetcut);
	
	tig_dopp_gg_beam_beamcut = new TH2F("TIGRESS_GG_Beam_Cut_Beam_Doppler","TIGRESS_GG_Beam_Cut_Beam_Doppler",2048,0,4096,2048,0,4096); gatedlist->Add(tig_dopp_gg_beam_beamcut);
	tig_dopp_gg_target_beamcut = new TH2F("TIGRESS_GG_Beam_Cut_Target_Doppler","TIGRESS_GG_Beam_Cut_Target_Doppler",2048,0,4096,2048,0,4096); gatedlist->Add(tig_dopp_gg_target_beamcut);
	tig_dopp_gg_beam_targetcut = new TH2F("TIGRESS_GG_Target_Cut_Beam_Doppler","TIGRESS_GG_Target_Cut_Beam_Doppler",2048,0,4096,2048,0,4096); gatedlist->Add(tig_dopp_gg_beam_targetcut);
	tig_dopp_gg_target_targetcut = new TH2F("TIGRESS_GG_Target_Cut_Target_Doppler","TIGRESS_GG_Target_Cut_Target_Doppler",2048,0,4096,2048,0,4096); gatedlist->Add(tig_dopp_gg_target_targetcut);

	s3_theta_fitted = new TH1F("S3_Theta_Fitted","S3_Theta_Fitted",180,0,180); gatedlist->Add(s3_theta_fitted);

	E_Theta_Gamma_Gate = new TH2F("E_Theta_Gamma_Gate","E_Theta_Gamma_Gate",180,0,180,2048,0,65536); timelist->Add(E_Theta_Gamma_Gate);

	hitmap_down = new TH2F("hitmap_down","Hitmap Downstream",200,-50,50,200,-50,50); s3list->Add(hitmap_down);
	hitmap_up = new TH2F("hitmap_up","Hitmap Upstream",200,-50,50,200,-50,50); s3list->Add(hitmap_up);

	downstream_rings = new TH2F("downstream_rings","downstream_rings",24,0,24,2048,0,524288); s3list->Add(downstream_rings);
	upstream_rings = new TH2F("upstream_rings","upstream_rings",24,0,24,2048,0,524288); s3list->Add(upstream_rings);
	downstream_sectors = new TH2F("downstream_sectors","downstream_sectors",32,0,32,2048,0,524288); s3list->Add(downstream_sectors);
	upstream_sectors = new TH2F("upstream_sectors","upstream_sectors",32,0,32,2048,0,524288); s3list->Add(upstream_sectors);


	for(int i=0;i<4;i++){
		char hname[64];
		sprintf(hname,"TIGRESS_Beam_Gated_Beam_Doppler_%i",i+1);
		tig_beam_beam_DetSet[i] = new TH1F(hname,hname,2048,0,2048); combolist->Add(tig_beam_beam_DetSet[i]);
	}
	for(int i=0;i<4;i++){
		char hname[64];
		sprintf(hname,"TIGRESS_Beam_Gated_Target_Doppler_%i",i+1);
		tig_beam_target_DetSet[i] = new TH1F(hname,hname,2048,0,2048); combolist->Add(tig_beam_target_DetSet[i]);
	}
	for(int i=0;i<4;i++){
		char hname[64];
		sprintf(hname,"TIGRESS_Target_Gated_Beam_Doppler_%i",i+1);
		tig_target_beam_DetSet[i] = new TH1F(hname,hname,2048,0,2048); combolist->Add(tig_target_beam_DetSet[i]);
	}
	for(int i=0;i<4;i++){
		char hname[64];
		sprintf(hname,"TIGRESS_Target_Gated_Target_Doppler_%i",i+1);
		tig_target_target_DetSet[i] = new TH1F(hname,hname,2048,0,2048); combolist->Add(tig_target_target_DetSet[i]);
	}
	
	theta_phi_hitmap = new TH2F("theta_phi_hitmap","theta_phi_hitmap",1000,10,60,800,-200,200); s3list->Add(theta_phi_hitmap);
	theta_phi_hitmap_corr = new TH2F("theta_phi_hitmap_corr","theta_phi_hitmap_corr",1000,10,60,800,-200,200); s3list->Add(theta_phi_hitmap_corr);

	for(int i=0;i<24;i++){
		char hname[128];
		sprintf(hname,"theta_phi_hitmap_corr_%i",i+1);
		theta_phi_hitmap_corr_ring[i] = new TH2F(hname,hname,1000,10,60,800,-200,200); thetaphilist->Add(theta_phi_hitmap_corr_ring[i]);
	}
	for(int i=0;i<6;i++){
		char hname[128];
		sprintf(hname,"Hitmap_PD_Rings_%i_to_%i",i*4 + 1,(i+1)*4 + 1);
		theta_phi_hitmap_corr_ring_group[i] = new TH2F(hname,hname,1000,10,60,800,-200,200); thetaphilist->Add(theta_phi_hitmap_corr_ring_group[i]);
	}


	thetaphi_iter = new TList;
	for(int i=0;i<21;i++){
		for(int j=0;j<21;j++){
			char hname[256];
			float x = i*0.5 - 5;
			float y = j*0.5 - 5;		
			sprintf(hname,"Hitmap_iter_%f_%f",x,y);
			thetaphi_offset_iter[i][j] = new TH2F(hname,hname,1000,10,60,800,-200,200); thetaphi_iter->Add(thetaphi_offset_iter[i][j]);

		}
	}
niki*/
}

void Sortcode::InitialiseFragments(){

	cout << "Creating fragment list" << endl;
	fraglist = new TList;
	
	cout << "Creating fragment histograms" << endl;
	charge = new TH1F("Charge","Charge",2000,0,5000000); fraglist->Add(charge);
	channel_charge = new TH2F("Channel_Charge","Channel_Charge",2000,0,2000,2000,0,5000000); fraglist->Add(channel_charge);
	//cout << channel_charge->GetNbinsX() << endl;

}

#endif
