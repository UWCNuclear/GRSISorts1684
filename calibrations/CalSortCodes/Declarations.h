TList *rawlist, *timelist, *gatedlist, *ringlist, *pbgated, *ringlist2, *pb_ringlist, *pb_ringlist2, *laserlist, *s3list, *combolist;
TList *thetaphilist;

TH1F *tig_mult, *addback_mult, *tig_E, *addback_E, *s3_E, *tig_E_dopp, *addback_E_dopp, *laser_amp;
TH1F *rf_time, *tig_time, *tig_s3_dCfd, *s3_time;

TH2F *front_back_energy, *tig_s3_time, *s3_E_theta, *laser_amp_addback_E, *laser_wf_overlay, *laser_amp_ts, *s3_ring_charge,*addback_E_channel;

TH2F *tigE_angle_beamcut, *tigdopp_angle_beamcut;
TH2F *tigE_angle_targetcut, *tigdopp_angle_targetcut;
TH2F *tigE_dCfd;

TH1F *s3_dCfd; 
TH2F *s3_E_theta_time_gated, *s3_E_tig_E, *s3_E_dopp_E;
TH1F *s3_E_time_gated, *tig_E_time_gated, *doppler_E_time_gated, *tig_dopp_s3_corr, *tig_dopp_kin_corr;
TH1F *s3_theta, *tig_dopp_s3_corr_noab, *tig_dopp_kin_corr_noab;

TH1F *tig_dopp_kin_kincut, *tig_dopp_kin_kincut_lasers, *tig_dopp_s3_corr_kincut;
TH1F *tig_dopp_kin_upstream;
TH2F *tig_E_S3_E_upstream, *tig_E_doppler_theta;

TH1F *tig_beam_beam_DetSet[4];
TH1F *tig_beam_target_DetSet[4];
TH1F *tig_target_beam_DetSet[4];
TH1F *tig_target_target_DetSet[4];

TH1F *tig_dopp_beam_beamcut;
TH1F *tig_dopp_target_beamcut;
TH1F *tig_dopp_beam_upstream; 
TH1F *tig_dopp_target_upstream; 
TH1F *tig_dopp_beam_targetcut;
TH1F *tig_dopp_target_targetcut;

TH2F *tig_dopp_gg_beam_beamcut; 
TH2F *tig_dopp_gg_target_beamcut;
TH2F *tig_dopp_gg_beam_targetcut;
TH2F *tig_dopp_gg_target_targetcut;

TH1F *Na_tig_E_downstream_ring[24];
TH1F *Na_tig_E_downstream_ring_ti_gated[24];
TH1F *Ti_tig_E_downstream_ring[24];
TH1F *Ti_tig_E_downstream_ring_ti_gated[24];

TH1F *tig_E_22Ne_Kin_Ti_gated, *tig_E_110Pd_Kin_Ti_gated;

TH1F *laser_waveform[1000];

TH2F *tig_dopp_kin_kincut_vs_laseramp;

TH1F *s3_theta_fitted;

TH2F *E_Theta_Gamma_Gate;

TH2F *hitmap_down, *hitmap_up;

TList *thetaphi_iter;
TH2F *thetaphi_offset_iter[21][21];

TH2F *downstream_rings, *upstream_rings;
TH2F *downstream_sectors, *upstream_sectors;

TH1F *four_ring_combo[6];
TH1F *four_ring_comboPD[6];

TH2F *theta_phi_hitmap, *theta_phi_hitmap_corr, *theta_phi_hitmap_corr_ring[24], *theta_phi_hitmap_corr_ring_group[6];

TList *fraglist;

TH2F *channel_charge;
TH1F *charge;
