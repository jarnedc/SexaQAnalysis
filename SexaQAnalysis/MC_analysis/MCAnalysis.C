
#include "./headers/HQClass.h"
#include "./headers/MCAnalysis.h"

///CONVENTIONS:
/// PV = primary vertex position
/// DV = decay vertex position (= daughterCV)
/// CV = creation vertex position (= momDV)

//Splitting up the plots in 3 categories:
//c1 = prompt production
//c2 = material interaction
//c3 = decays
//c123 = everything put together


void bookHistos(std::map<string, TH1F*>& histos_1d, std::map<string, TH2F*>& histos_2d, std::map<string, TEfficiency*>& histos_eff, TString c){

	bool c1 = (c=="c1");
	bool c2 = (c=="c2");
	bool c3 = (c=="c3");
	bool c123 = (c=="c123");

	c="";
	//c=c+"_";

	///XI PLUS MINUS
	/*
	histos_1d["inv_mass_Xi_plus_delta_z_1mm"] = new TH1F("inv_mass_Xi_plus_delta_z_1mm", "inv_mass_Xi_plus_delta_z_1mm; inv mass (GeV)", 100, 0.0, 10);
	histos_1d["inv_mass_Xi_min_delta_z_1mm"] = new TH1F("inv_mass_Xi_min_delta_z_1mm", "inv_mass_Xi_min_delta_z_1mm; inv mass (GeV)", 100, 0.0, 10);
	histos_1d["inv_mass_Xi_plus_delta_z_5mm"] = new TH1F("inv_mass_Xi_plus_delta_z_5mm", "inv_mass_Xi_plus_delta_z_5mm; inv mass (GeV)", 100, 0.0, 10);
	histos_1d["inv_mass_Xi_min_delta_z_5mm"] = new TH1F("inv_mass_Xi_min_delta_z_5mm", "inv_mass_Xi_min_delta_z_5mm; inv mass (GeV)", 100, 0.0, 10);
	
	histos_1d["inv_mass_Xi_plus"] = new TH1F("inv_mass_Xi_plus", "inv_mass_Xi_plus; inv mass (GeV)", 200, 0.0, 10);
	histos_1d["inv_mass_Xi_min"] = new TH1F("inv_mass_Xi_min", "inv_mass_Xi_min; inv mass (GeV)", 200, 0.0, 10);
	
	histos_1d["inv_mass_Xi_plus_zoom"] = new TH1F("inv_mass_Xi_plus_zoom", "inv_mass_Xi_plus_zoom; inv mass (GeV)", 100, 1, 2);
	histos_1d["inv_mass_Xi_min_zoom"] = new TH1F("inv_mass_Xi_min_zoom", "inv_mass_Xi_min_zoom; inv mass (GeV)", 100, 1, 2);
	*/
	
	///INV MASS LAMBDA, KSHORT, S
	histos_1d["mass_S_reco"] = new TH1F("mass_S_reco", "mass_S_reco; mass (GeV)", 400, -20, 20);
	histos_1d["mass_S_gen"] = new TH1F("mass_S_gen", c+"mass_S_gen; mass (GeV)", 400, -20, 20);
	
	histos_1d["inv_mass_Ks_L0_reco"] = new TH1F("inv_mass_Ks_L0_reco", "inv_mass_Ks_L0_reco; mass (GeV)", 200, 0, 20);
	histos_1d["inv_mass_Ks_L0_gen"] = new TH1F("inv_mass_Ks_L0_gen", c+"inv_mass_Ks_L0_gen; mass (GeV)", 200, 0, 20);
	
	
	///LAMBDA KSHORT POSITION
	/*
	histos_1d["Lambda_z"] = new TH1F("Lambda_z", "Lambda_z; z-postion (cm)", 100, -20, 20);
	histos_1d["Lambda_impact_param"] = new TH1F("Lambda_impact_param", "Lambda_impact_param; impact param (cm)", 100, 0, 0.5);
	
	histos_1d["Lambda_impact_param_closest_appr_vtx"] = new TH1F("Lambda_impact_param_closest_appr_vtx", "Lambda_impact_param_closest_appr_vtx; impact param (cm)", 100, 0, 0.5);
	histos_1d["Kshort_impact_param_closest_appr_vtx"] = new TH1F("Kshort_impact_param_closest_appr_vtx", "Kshort_impact_param_closest_appr_vtx; impact param (cm)", 100, 0, 0.5);
	
	histos_1d["Lambda_z_closest_approach_beam"] = new TH1F("Lambda_z_closest_approach_beam", "Lambda_z_closest_approach_beam; z-position (cm)", 100, -20, 20);
	*/

	///LAMBDA KSHORT EXTRAPOLATION CLOSEST APPROACH TO EACHOTHER
	/*
	histos_1d["lambda_kshort_dist_at_closest_appr"] = new TH1F("lambda_kshort_dist_at_closest_appr", "lambda_kshort_dist_at_closest_appr; 3D-distance (cm)", 300, 0, 5);
	histos_1d["lambda_kshort_delta_z_at_closest_appr"] = new TH1F("lambda_kshort_delta_z_at_closest_appr", "lambda_kshort_delta_z_at_closest_appr; z-distance (cm)", 300, 0, 5);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_delta_z"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_delta_z", "lambda_kshort_avrg_pos_to_prim_vtx_delta_z; z-distance (cm)", 300, 0, 12);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_dist"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_dist", "lambda_kshort_avrg_pos_to_prim_vtx_dist; 3D-distance (cm)", 300, 0, 12);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_impact_param"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_impact_param", "lambda_kshort_avrg_pos_to_prim_vtx_impact_param; (x,y)-distance (cm)", 300, 0, 5);
	
	histos_1d["lambda_kshort_dist_at_closest_appr_zoom"] = new TH1F("lambda_kshort_dist_at_closest_appr_zoom", "lambda_kshort_dist_at_closest_appr_zoom; 3D-distance (cm)", 300, 0, 1);
	histos_1d["lambda_kshort_delta_z_at_closest_appr_zoom"] = new TH1F("lambda_kshort_delta_z_at_closest_appr_zoom", "lambda_kshort_delta_z_at_closest_appr_zoom; z-distance (cm)", 300, 0, 1);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_delta_z_zoom"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_delta_z_zoom", "lambda_kshort_avrg_pos_to_prim_vtx_delta_z_zoom; z-distance (cm)", 300, 0, 1);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_dist_zoom"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_dist_zoom", "lambda_kshort_avrg_pos_to_prim_vtx_dist_zoom; 3D-distance (cm)", 300, 0, 1);
	histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_impact_param_zoom"] = new TH1F("lambda_kshort_avrg_pos_to_prim_vtx_impact_param_zoom", "lambda_kshort_avrg_pos_to_prim_vtx_impact_param_zoom; (x,y)-distance (cm)", 300, 0, 1);
	*/
	
	///Ks AND L0 STATUS 8 & STATUS 1 MOTHER INVESTIGATION
	
	histos_1d["status_8_Ks_L0_(grand)mother_pdgid"] = new TH1F("status_8_Ks_L0_(grand)mother_pdgid", c+"status_8_Ks_L0_(grand)mother_pdgid; pdgid of (grand)mothers of status 8 gen L0 and Ks", 3500, -0.5, 3499.5);
	histos_1d["status_8_Ks_(grand)mother_pdgid"] = new TH1F("status_8_Ks_(grand)mother_pdgid", c+"status_8_Ks_(grand)mother_pdgid; pdgid of (grand)mothers of status 8 gen Ks", 3500, -0.5, 3499.5);
	histos_1d["status_8_L0_(grand)mother_pdgid"] = new TH1F("status_8_L0_(grand)mother_pdgid", c+"status_8_L0_(grand)mother_pdgid; pdgid of (grand)mothers of status 8 gen L0", 3500, -0.5, 3499.5);

	//histos_1d["status_1_Ks_(grand)mother_pdgid"] = new TH1F("status_1_Ks_(grand)mother_pdgid", c+"status_1_Ks_(grand)mother_pdgid; pdgid of (grand)mothers of status 1 gen Ks", 10000500, -0.5, 10000499.5);
	//histos_1d["status_1_L0_(grand)mother_pdgid"] = new TH1F("status_1_L0_(grand)mother_pdgid", c+"status_1_L0_(grand)mother_pdgid; pdgid of (grand)mothers of status 1 gen L0", 10000500, -0.5, 10000499.5);


	///CATEGORY 1, 2 or 3 INVESTIGATION

	//histos_1d["status_8_pdgid"] = new TH1F("status_8_pdgid", "status_8_pdgid; pdgid of status 8 gen particles", 10000, -5000, 5000);



	///OTHER USEFUL DISTRIBUTIONS

	
	histos_1d["reco_Ks_L0_overlap"] = new TH1F("reco_Ks_L0_overlap", "reco_Ks_L0_overlap; Same (x,y,z) and (px,py,pz) or not?", 3, -0.25, 1.25);

	histos_1d["status_8_pdgid"] = new TH1F("status_8_pdgid", "status_8_pdgid; pdgid of status 8 gen particles", 10000, -5000, 5000);
	histos_1d["L0_with_correct_daughters_status"] = new TH1F("L0_with_correct_daughters_status", "L0_with_correct_daughters_status; status", 10, 0, 10);
	histos_1d["Ks_with_correct_daughters_status"] = new TH1F("Ks_with_correct_daughters_status", "Ks_with_correct_daughters_status; status", 10, 0, 10);

	///XI(1530) AND RESONANCE CANDIDATE INVESTIGATION

	histos_1d["Xi(1530)_0_gen_Eta"] = new TH1F("Xi(1530)_0_gen_Eta","Xi(1530)_0_gen_Eta;Eta;",100,-10,10);
	histos_1d["Xi(1530)_minus_gen_Eta"] = new TH1F("Xi(1530)_minus_gen_Eta","Xi(1530)_minus_gen_Eta;Eta;",100,-10,10);




	///RECO LAMBDA & KSHORT Pt, Eta, ... DISTRIBUTIONS
	
	histos_1d["L0_reco_amount_per_event"] = new TH1F("L0_reco_amount_per_event", "L0_reco_amount_per_event; Number of", 5, 0, 5);
	histos_1d["Ks_reco_amount_per_event"] = new TH1F("Ks_reco_amount_per_event", "Ks_reco_amount_per_event; Number of", 5, 0, 5);
	histos_1d["Ks_reco_amount_per_event_nodoublecount"] = new TH1F("Ks_reco_amount_per_event_nodoublecount", "Ks_reco_amount_per_event_nodoublecount; Number of", 5, 0, 5);

	histos_1d["L0_reco_Pt"] = new TH1F("L0_reco_Pt", "L0_reco_Pt; Pt (GeV)", 100, 0, (c1?10:c2?10:c3?10:c123?10:0));
	histos_1d["L0_reco_Eta"] = new TH1F("L0_reco_Eta", "L0_reco_Eta; Eta", 100, -5, 5);
	histos_1d["L0_reco_rapidity"] = new TH1F("L0_reco_rapidity", "L0_reco_rapidity; rapidity", 100, -5, 5);
	
	
	histos_1d["Ks_reco_Pt"] = new TH1F("Ks_reco_Pt", "Ks_reco_Pt; Pt (GeV)", 100, 0, (c1?10:c2?10:c3?10:c123?10:0));
	histos_1d["Ks_reco_Eta"] = new TH1F("Ks_reco_Eta", "Ks_reco_Eta; Eta", 100, -5, 5);
	histos_1d["Ks_reco_rapidity"] = new TH1F("Ks_reco_rapidity", "Ks_reco_rapidity; rapidity", 100, -5, 5);
	
	
	///LAMBDA DECAY POSITION AND RECO EFFICIENCY
	
	histos_1d["L0_gen_dxyz(DV_PV)"] = new TH1F("L0_gen_dxyz(DV_PV)", c+"L0_gen_dxyz(DV_PV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(DV_PV)"] = new TH1F("L0_gen_dxy(DV_PV)", c+"L0_gen_dxy(DV_PV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["L0_gen_dz(DV_PV)"] = new TH1F("L0_gen_dz(DV_PV)", c+"L0_gen_dz(DV_PV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["L0_gen_dxyz(DV_CV)"] = new TH1F("L0_gen_dxyz(DV_CV)", c+"L0_gen_dxyz(DV_CV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(DV_CV)"] = new TH1F("L0_gen_dxy(DV_CV)", c+"L0_gen_dxy(DV_CV); xy-distance (cm)", 300, 0, (c1?60:c2?60:c3?60:c123?60:0));
	histos_1d["L0_gen_dz(DV_CV)"] = new TH1F("L0_gen_dz(DV_CV)", c+"L0_gen_dz(DV_CV); z-distance (cm)", 500, -(c1?250:c2?250:c3?250:c123?250:0), (c1?250:c2?250:c3?250:c123?250:0));
	
	histos_1d["L0_gen_dxyz(CV_PV)"] = new TH1F("L0_gen_dxyz(CV_PV)",c+"L0_gen_dxyz(CV_PV);dxyz(CV_PV) (cm);",300,0,(c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(CV_PV)"] = new TH1F("L0_gen_dxy(CV_PV)",c+"L0_gen_dxy(CV_PV);dxy(CV_PV) (cm);",300,0,(c1?60:c2?60:c3?60:c123?60:0));
	histos_1d["L0_gen_dz(CV_PV)"] = new TH1F("L0_gen_dz(CV_PV)",c+"L0_gen_dz(CV_PV);dz(CV_PV) (cm);",500,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_1d["L0_gen_dxyz(DV_momCV)"] = new TH1F("L0_gen_dxyz(DV_momCV)", c+"L0_gen_dxyz(DV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(DV_momCV)"] = new TH1F("L0_gen_dxy(DV_momCV)", c+"L0_gen_dxy(DV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["L0_gen_dz(DV_momCV)"] = new TH1F("L0_gen_dz(DV_momCV)", c+"L0_gen_dz(DV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["L0_gen_dxyz(PV_momCV)"] = new TH1F("L0_gen_dxyz(PV_momCV)", c+"L0_gen_dxyz(PV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(PV_momCV)"] = new TH1F("L0_gen_dxy(PV_momCV)", c+"L0_gen_dxy(PV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["L0_gen_dz(PV_momCV)"] = new TH1F("L0_gen_dz(PV_momCV)", c+"L0_gen_dz(PV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["L0_gen_dxyz(CV_momCV)"] = new TH1F("L0_gen_dxyz(CV_momCV)", c+"L0_gen_dxyz(CV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["L0_gen_dxy(CV_momCV)"] = new TH1F("L0_gen_dxy(CV_momCV)", c+"L0_gen_dxy(CV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["L0_gen_dz(CV_momCV)"] = new TH1F("L0_gen_dz(CV_momCV)", c+"L0_gen_dz(CV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	
	histos_1d["L0_gen_reco_deltaR"] = new TH1F("L0_gen_reco_deltaR", c+"L0_gen_reco_deltaR; delta_R", 200, 0, 8);
	histos_1d["L0_gen_reco_deltaR_zoom"] = new TH1F("L0_gen_reco_deltaR_zoom", c+"L0_gen_reco_deltaR_zoom; delta_R", 200, 0, 0.5);
	
	
	histos_1d["L0_gen_Pt"] = new TH1F("L0_gen_Pt",c+"L0_gen_Pt;Pt (GeV);",100,0,(c1?5:c2?5:c3?5:c123?5:0));
	histos_1d["L0_gen_Eta"] = new TH1F("L0_gen_Eta",c+"L0_gen_Eta;Eta;",100,-10,10);
	histos_1d["L0_gen_rapidity"] = new TH1F("L0_gen_rapidity",c+"L0_gen_rapidity;rapidity;",100,-10,10);
	histos_1d["L0_gen_dxy(CV_origin)"] = new TH1F("L0_gen_dxy(CV_origin)",c+"L0_gen_dxy(CV_origin);dxy(CV_origin) (cm);",500,0,(c1?50:c2?50:c3?50:c123?50:0));
	
	
	histos_1d["L0_gen_dxy(daughtersCV_origin)"] = new TH1F("L0_gen_dxy(daughtersCV_origin)",c+"L0_gen_dxy(daughtersCV_origin);dxy(daughtersCV_origin) (cm);",500,0,(c1?50:c2?50:c3?50:c123?50:0));

	histos_eff["L0_recoeff_gen_dxyz(CV_PV)"] = new TEfficiency("L0_recoeff_gen_dxyz(CV_PV)",c+"L0_recoeff_gen_dxyz(CV_PV);dxyz(PV, CV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(CV_PV)"] = new TEfficiency("L0_recoeff_gen_dxy(CV_PV)",c+"L0_recoeff_gen_dxy(CV_PV);dxy(PV, CV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(CV_PV)"] = new TEfficiency("L0_recoeff_gen_dz(CV_PV)",c+"L0_recoeff_gen_dz(CV_PV);dz(PV, CV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_eff["L0_recoeff_gen_dxyz(DV_PV)"] = new TEfficiency("L0_recoeff_gen_dxyz(DV_PV)",c+"L0_recoeff_gen_dxyz(DV_PV);dxyz(DV, PV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(DV_PV)"] = new TEfficiency("L0_recoeff_gen_dxy(DV_PV)",c+"L0_recoeff_gen_dxy(DV_PV);dxy(DV, PV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(DV_PV)"] = new TEfficiency("L0_recoeff_gen_dz(DV_PV)",c+"L0_recoeff_gen_dz(DV_PV);dz(DV, PV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));
	
	histos_eff["L0_recoeff_gen_dxyz(DV_CV)"] = new TEfficiency("L0_recoeff_gen_dxyz(DV_CV)",c+"L0_recoeff_gen_dxyz(DV_CV);dxyz(DV, CV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(DV_CV)"] = new TEfficiency("L0_recoeff_gen_dxy(DV_CV)",c+"L0_recoeff_gen_dxy(DV_CV);dxy(DV, CV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(DV_CV)"] = new TEfficiency("L0_recoeff_gen_dz(DV_CV)",c+"L0_recoeff_gen_dz(DV_CV);dz(DV, CV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_eff["L0_recoeff_gen_dxyz(DV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxyz(DV_momCV)",c+"L0_recoeff_gen_dxyz(DV_momCV);dxyz(DV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(DV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxy(DV_momCV)",c+"L0_recoeff_gen_dxy(DV_momCV);dxy(DV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(DV_momCV)"] = new TEfficiency("L0_recoeff_gen_dz(DV_momCV)",c+"L0_recoeff_gen_dz(DV_momCV);dz(DV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));
	
	histos_eff["L0_recoeff_gen_dxyz(PV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxyz(PV_momCV)",c+"L0_recoeff_gen_dxyz(PV_momCV);dxyz(PV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(PV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxy(PV_momCV)",c+"L0_recoeff_gen_dxy(PV_momCV);dxy(PV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(PV_momCV)"] = new TEfficiency("L0_recoeff_gen_dz(PV_momCV)",c+"L0_recoeff_gen_dz(PV_momCV);dz(PV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_eff["L0_recoeff_gen_dxyz(CV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxyz(CV_momCV)",c+"L0_recoeff_gen_dxyz(CV_momCV);dxyz(CV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["L0_recoeff_gen_dxy(CV_momCV)"] = new TEfficiency("L0_recoeff_gen_dxy(CV_momCV)",c+"L0_recoeff_gen_dxy(CV_momCV);dxy(CV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["L0_recoeff_gen_dz(CV_momCV)"] = new TEfficiency("L0_recoeff_gen_dz(CV_momCV)",c+"L0_recoeff_gen_dz(CV_momCV);dz(CV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	
	histos_eff["L0_recoeff_gen_Pt"] = new TEfficiency("L0_recoeff_gen_Pt",c+"L0_recoeff_gen_Pt;Pt (GeV);reconstruction efficiency",50,0,10);
	histos_eff["L0_recoeff_gen_Eta"] = new TEfficiency("L0_recoeff_gen_Eta",c+"L0_recoeff_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency",50,-2.5,2.5);
	histos_eff["L0_recoeff_gen_rapidity"] = new TEfficiency("L0_recoeff_gen_rapidity",c+"L0_recoeff_gen_rapidity;rapidity;reconstruction efficiency",50,-2.5,2.5);


	histos_eff["L0_recoeff_min_gen_Pt_daughters"] = new TEfficiency("L0_recoeff_min_gen_Pt_daughters",c+"L0_recoeff_min_gen_Pt_daughters;Pt (GeV);reconstruction efficiency",50,0,10);
	

	histos_1d["L0_recoeff_passed_gen_Pt"] = new TH1F("L0_recoeff_passed_gen_Pt",c+"L0_recoeff_passed_gen_Pt;Pt (GeV);reconstruction efficiency passed",50,0,8);
	histos_1d["L0_recoeff_passed_gen_Eta"] = new TH1F("L0_recoeff_passed_gen_Eta",c+"L0_recoeff_passed_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency passed",50,-10,10);
	histos_1d["L0_recoeff_passed_gen_rapidity"] = new TH1F("L0_recoeff_passed_gen_rapidity",c+"L0_recoeff_passed_gen_rapidity;rapidity;reconstruction efficiency passed",50,-10,10);
	histos_1d["L0_recoeff_total_gen_Pt"] = new TH1F("L0_recoeff_total_gen_Pt",c+"L0_recoeff_total_gen_Pt;Pt (GeV);reconstruction efficiency total",50,0,8);
	histos_1d["L0_recoeff_total_gen_Eta"] = new TH1F("L0_recoeff_total_gen_Eta",c+"L0_recoeff_total_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency total",50,-10,10);
	histos_1d["L0_recoeff_total_gen_rapidity"] = new TH1F("L0_recoeff_total_gen_rapidity",c+"L0_recoeff_total_gen_rapidity;rapidity;reconstruction efficiency total",50,-10,10);

	
	histos_eff["L0_recoeff_gen_dxy(CV_origin)"] = new TEfficiency("L0_recoeff_gen_dxy(CV_origin)",c+"L0_recoeff_gen_dxy(CV_origin);dxy(CV_origin) (cm) ;reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));

	///KSHORT DECAY POSITION AND RECO EFFICIENCY
	
	histos_1d["Ks_gen_dxyz(DV_PV)"] = new TH1F("Ks_gen_dxyz(DV_PV)", c+"Ks_gen_dxyz(DV_PV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(DV_PV)"] = new TH1F("Ks_gen_dxy(DV_PV)", c+"Ks_gen_dxy(DV_PV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["Ks_gen_dz(DV_PV)"] = new TH1F("Ks_gen_dz(DV_PV)", c+"Ks_gen_dz(DV_PV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["Ks_gen_dxyz(DV_CV)"] = new TH1F("Ks_gen_dxyz(DV_CV)", c+"Ks_gen_dxyz(DV_CV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(DV_CV)"] = new TH1F("Ks_gen_dxy(DV_CV)", c+"Ks_gen_dxy(DV_CV); xy-distance (cm)", 300, 0, (c1?60:c2?60:c3?60:c123?60:0));
	histos_1d["Ks_gen_dz(DV_CV)"] = new TH1F("Ks_gen_dz(DV_CV)", c+"Ks_gen_dz(DV_CV); z-distance (cm)", 500, -(c1?250:c2?250:c3?250:c123?250:0), (c1?250:c2?250:c3?250:c123?250:0));
	
	histos_1d["Ks_gen_dxyz(CV_PV)"] = new TH1F("Ks_gen_dxyz(CV_PV)",c+"Ks_gen_dxyz(CV_PV);dxyz(CV_PV) (cm);",300,0,(c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(CV_PV)"] = new TH1F("Ks_gen_dxy(CV_PV)",c+"Ks_gen_dxy(CV_PV);dxy(CV_PV) (cm);",300,0,(c1?60:c2?60:c3?60:c123?60:0));
	histos_1d["Ks_gen_dz(CV_PV)"] = new TH1F("Ks_gen_dz(CV_PV)",c+"Ks_gen_dz(CV_PV);dz(CV_PV) (cm);",500,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_1d["Ks_gen_dxyz(DV_momCV)"] = new TH1F("Ks_gen_dxyz(DV_momCV)", c+"Ks_gen_dxyz(DV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(DV_momCV)"] = new TH1F("Ks_gen_dxy(DV_momCV)", c+"Ks_gen_dxy(DV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["Ks_gen_dz(DV_momCV)"] = new TH1F("Ks_gen_dz(DV_momCV)", c+"Ks_gen_dz(DV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["Ks_gen_dxyz(PV_momCV)"] = new TH1F("Ks_gen_dxyz(PV_momCV)", c+"Ks_gen_dxyz(PV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(PV_momCV)"] = new TH1F("Ks_gen_dxy(PV_momCV)", c+"Ks_gen_dxy(PV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["Ks_gen_dz(PV_momCV)"] = new TH1F("Ks_gen_dz(PV_momCV)", c+"Ks_gen_dz(PV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	histos_1d["Ks_gen_dxyz(CV_momCV)"] = new TH1F("Ks_gen_dxyz(CV_momCV)", c+"Ks_gen_dxyz(CV_momCV); xyz-distance (cm)", 300, 0, (c1?400:c2?400:c3?400:c123?400:0));
	histos_1d["Ks_gen_dxy(CV_momCV)"] = new TH1F("Ks_gen_dxy(CV_momCV)", c+"Ks_gen_dxy(CV_momCV); xy-distance (cm)", 300, 0, (c1?120:c2?120:c3?120:c123?120:0));
	histos_1d["Ks_gen_dz(CV_momCV)"] = new TH1F("Ks_gen_dz(CV_momCV)", c+"Ks_gen_dz(CV_momCV); z-distance (cm)", 500, -(c1?500:c2?500:c3?500:c123?500:0), (c1?500:c2?500:c3?500:c123?500:0));
	
	
	histos_1d["Ks_gen_reco_deltaR"] = new TH1F("Ks_gen_reco_deltaR", c+"Ks_gen_reco_deltaR; delta_R", 200, 0, 8);
	histos_1d["Ks_gen_reco_deltaR_zoom"] = new TH1F("Ks_gen_reco_deltaR_zoom", c+"Ks_gen_reco_deltaR_zoom; delta_R", 200, 0, 0.5);

	
	histos_1d["Ks_gen_Pt"] = new TH1F("Ks_gen_Pt",c+"Ks_gen_Pt;Pt (GeV);",100,0,5);
	histos_1d["Ks_gen_Eta"] = new TH1F("Ks_gen_Eta",c+"Ks_gen_Eta;Eta;",100,-10,10);
	histos_1d["Ks_gen_rapidity"] = new TH1F("Ks_gen_rapidity",c+"Ks_gen_rapidity;rapidity;",100,-10,10);

	histos_1d["Ks_gen_dxy(CV_origin)"] = new TH1F("Ks_gen_dxy(CV_origin)",c+"Ks_gen_dxy(CV_origin);dxy(CV_origin) (cm);",500,0,(c1?50:c2?50:c3?50:c123?50:0));
	
	histos_1d["Ks_gen_dxy(daughtersCV_origin)"] = new TH1F("Ks_gen_dxy(daughtersCV_origin)",c+"Ks_gen_dxy(daughtersCV_origin);dxy(daughtersCV_origin) (cm);",500,0,(c1?50:c2?50:c3?50:c123?50:0));


	histos_eff["Ks_recoeff_gen_dxyz(CV_PV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(CV_PV)",c+"Ks_recoeff_gen_dxyz(CV_PV);dxyz(PV, CV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(CV_PV)"] = new TEfficiency("Ks_recoeff_gen_dxy(CV_PV)",c+"Ks_recoeff_gen_dxy(CV_PV);dxy(PV, CV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(CV_PV)"] = new TEfficiency("Ks_recoeff_gen_dz(CV_PV)",c+"Ks_recoeff_gen_dz(CV_PV);dz(PV, CV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_eff["Ks_recoeff_gen_dxyz(DV_PV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(DV_PV)",c+"Ks_recoeff_gen_dxyz(DV_PV);dxyz(DV, PV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(DV_PV)"] = new TEfficiency("Ks_recoeff_gen_dxy(DV_PV)",c+"Ks_recoeff_gen_dxy(DV_PV);dxy(DV, PV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(DV_PV)"] = new TEfficiency("Ks_recoeff_gen_dz(DV_PV)",c+"Ks_recoeff_gen_dz(DV_PV);dz(DV, PV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));
	
	histos_eff["Ks_recoeff_gen_dxyz(DV_CV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(DV_CV)",c+"Ks_recoeff_gen_dxyz(DV_CV);dxyz(DV, CV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(DV_CV)"] = new TEfficiency("Ks_recoeff_gen_dxy(DV_CV)",c+"Ks_recoeff_gen_dxy(DV_CV);dxy(DV, CV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(DV_CV)"] = new TEfficiency("Ks_recoeff_gen_dz(DV_CV)",c+"Ks_recoeff_gen_dz(DV_CV);dz(DV, CV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	histos_eff["Ks_recoeff_gen_dxyz(DV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(DV_momCV)",c+"Ks_recoeff_gen_dxyz(DV_momCV);dxyz(DV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(DV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxy(DV_momCV)",c+"Ks_recoeff_gen_dxy(DV_momCV);dxy(DV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(DV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dz(DV_momCV)",c+"Ks_recoeff_gen_dz(DV_momCV);dz(DV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));
	
	histos_eff["Ks_recoeff_gen_dxyz(PV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(PV_momCV)",c+"Ks_recoeff_gen_dxyz(PV_momCV);dxyz(PV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(PV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxy(PV_momCV)",c+"Ks_recoeff_gen_dxy(PV_momCV);dxy(PV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(PV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dz(PV_momCV)",c+"Ks_recoeff_gen_dz(PV_momCV);dz(PV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));
	
	histos_eff["Ks_recoeff_gen_dxyz(CV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxyz(CV_momCV)",c+"Ks_recoeff_gen_dxyz(CV_momCV);dxyz(CV, momCV);reconstruction efficiency",50,0,(c1?500:c2?500:c3?500:c123?500:0)); 
	histos_eff["Ks_recoeff_gen_dxy(CV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dxy(CV_momCV)",c+"Ks_recoeff_gen_dxy(CV_momCV);dxy(CV, momCV);reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));
	histos_eff["Ks_recoeff_gen_dz(CV_momCV)"] = new TEfficiency("Ks_recoeff_gen_dz(CV_momCV)",c+"Ks_recoeff_gen_dz(CV_momCV);dz(CV, momCV);reconstruction efficiency",50,-(c1?250:c2?250:c3?250:c123?250:0),(c1?250:c2?250:c3?250:c123?250:0));

	
	histos_eff["Ks_recoeff_gen_Pt"] = new TEfficiency("Ks_recoeff_gen_Pt",c+"Ks_recoeff_gen_Pt;Pt (GeV);reconstruction efficiency",50,0,10);
	histos_eff["Ks_recoeff_gen_Eta"] = new TEfficiency("Ks_recoeff_gen_Eta",c+"Ks_recoeff_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency",50,-2.5,2.5);
	histos_eff["Ks_recoeff_gen_rapidity"] = new TEfficiency("Ks_recoeff_gen_rapidity",c+"Ks_recoeff_gen_rapidity;rapidity;reconstruction efficiency",50,-2.5,2.5);

	histos_eff["Ks_recoeff_min_gen_Pt_daughters"] = new TEfficiency("Ks_recoeff_min_gen_Pt_daughters",c+"Ks_recoeff_min_gen_Pt_daughters;Pt (GeV);reconstruction efficiency",50,0,10);

	
	histos_1d["Ks_recoeff_passed_gen_Pt"] = new TH1F("Ks_recoeff_passed_gen_Pt",c+"Ks_recoeff_passed_gen_Pt;Pt (GeV);reconstruction efficiency passed",50,0,8);
	histos_1d["Ks_recoeff_passed_gen_Eta"] = new TH1F("Ks_recoeff_passed_gen_Eta",c+"Ks_recoeff_passed_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency passed",50,-10,10);
	histos_1d["Ks_recoeff_passed_gen_rapidity"] = new TH1F("Ks_recoeff_passed_gen_rapidity",c+"Ks_recoeff_passed_gen_rapidity;rapidity;reconstruction efficiency passed",50,-10,10);
	histos_1d["Ks_recoeff_total_gen_Pt"] = new TH1F("Ks_recoeff_total_gen_Pt",c+"Ks_recoeff_total_gen_Pt;Pt (GeV);reconstruction efficiency total",50,0,8);
	histos_1d["Ks_recoeff_total_gen_Eta"] = new TH1F("Ks_recoeff_total_gen_Eta",c+"Ks_recoeff_total_gen_Eta;Pseudorapidity (Eta);reconstruction efficiency total",50,-10,10);
	histos_1d["Ks_recoeff_total_gen_rapidity"] = new TH1F("Ks_recoeff_total_gen_rapidity",c+"Ks_recoeff_total_gen_rapidity;rapidity;reconstruction efficiency total",50,-10,10);

	
	histos_eff["Ks_recoeff_gen_dxy(CV_origin)"] = new TEfficiency("Ks_recoeff_gen_dxy(CV_origin)",c+"Ks_recoeff_gen_dxy(CV_origin);dxy(CV_origin) (cm) ;reconstruction efficiency",50,0,(c1?150:c2?150:c3?150:c123?150:0));

	///CREATION VERTEX SCATTER PLOTS
	
	histos_2d["Ks_gen_dx(CV_origin)_dy(CV_origin)"] = new TH2F("Ks_gen_dx(CV_origin)_dy(CV_origin)", c+"Ks_gen_dx(CV_origin)_dy(CV_origin); dx(CV_origin) (cm); dy(CV_origin) (cm)", 200, -80, 80, 200, -80, 80);
	histos_2d["L0_gen_dx(CV_origin)_dy(CV_origin)"] = new TH2F("L0_gen_dx(CV_origin)_dy(CV_origin)", c+"L0_gen_dx(CV_origin)_dy(CV_origin); dx(CV_origin) (cm); dy(CV_origin) (cm)", 200, -80, 80, 200, -80, 80);
	
	histos_2d["Ks_gen_dr(CV_origin)_dz(CV_origin)"] = new TH2F("Ks_gen_dr(CV_origin)_dz(CV_origin)", c+"Ks_gen_dr(CV_origin)_dz(CV_origin); dz(CV_origin) (cm); dr(CV_origin) (cm)", 300, -250, 250, 300, 0, 80);
	histos_2d["L0_gen_dr(CV_origin)_dz(CV_origin)"] = new TH2F("L0_gen_dr(CV_origin)_dz(CV_origin)", c+"L0_gen_dr(CV_origin)_dz(CV_origin); dz(CV_origin) (cm); dr(CV_origin) (cm)", 300, -250, 250, 300, 0, 80);


	histos_2d["Ks_gen_dx(CV_PV)_dy(CV_PV)"] = new TH2F("Ks_gen_dx(CV_PV)_dy(CV_PV)", c+"Ks_gen_dx(CV_PV)_dy(CV_PV); dx(CV_PV) (cm); dy(CV_PV) (cm)", 200, -80, 80, 200, -80, 80);
	histos_2d["L0_gen_dx(CV_PV)_dy(CV_PV)"] = new TH2F("L0_gen_dx(CV_PV)_dy(CV_PV)", c+"L0_gen_dx(CV_PV)_dy(CV_PV); dx(CV_PV) (cm); dy(CV_PV) (cm)", 200, -80, 80, 200, -80, 80);
	
	histos_2d["Ks_gen_dr(CV_PV)_dz(CV_PV)"] = new TH2F("Ks_gen_dr(CV_PV)_dz(CV_PV)", c+"Ks_gen_dr(CV_PV)_dz(CV_PV); dz(CV_PV) (cm); dr(CV_PV) (cm)", 300, -250, 250, 300, 0, 80);
	histos_2d["L0_gen_dr(CV_PV)_dz(CV_PV)"] = new TH2F("L0_gen_dr(CV_PV)_dz(CV_PV)", c+"L0_gen_dr(CV_PV)_dz(CV_PV); dz(CV_PV) (cm); dr(CV_PV) (cm)", 300, -250, 250, 300, 0, 80);


	histos_2d["Ks_gen_dx(CV_momCV)_dy(CV_momCV)"] = new TH2F("Ks_gen_dx(CV_momCV)_dy(CV_momCV)", c+"Ks_gen_dx(CV_momCV)_dy(CV_momCV); dx(CV_momCV) (cm); dy(CV_momCV) (cm)", 200, -80, 80, 200, -80, 80);
	histos_2d["L0_gen_dx(CV_momCV)_dy(CV_momCV)"] = new TH2F("L0_gen_dx(CV_momCV)_dy(CV_momCV)", c+"L0_gen_dx(CV_momCV)_dy(CV_momCV); dx(CV_momCV) (cm); dy(CV_momCV) (cm)", 200, -80, 80, 200, -80, 80);
	
	histos_2d["Ks_gen_dr(CV_momCV)_dz(CV_momCV)"] = new TH2F("Ks_gen_dr(CV_momCV)_dz(CV_momCV)", c+"Ks_gen_dr(CV_momCV)_dz(CV_momCV); dz(CV_momCV) (cm); dr(CV_momCV) (cm)", 300, -250, 250, 300, 0, 80);
	histos_2d["L0_gen_dr(CV_momCV)_dz(CV_momCV)"] = new TH2F("L0_gen_dr(CV_momCV)_dz(CV_momCV)", c+"L0_gen_dr(CV_momCV)_dz(CV_momCV); dz(CV_momCV) (cm); dr(CV_momCV) (cm)", 300, -250, 250, 300, 0, 80);


	histos_1d["L0_dr_dz_range_mothers"] = new TH1F("L0_dr_dz_range_mothers",c+"L0_dr_dz_range_mothers;mothers pdgid",7600,-3800,3800);
	histos_1d["Ks_dr_dz_range_mothers"] = new TH1F("Ks_dr_dz_range_mothers",c+"Ks_dr_dz_range_mothers;mothers pdgid",7600,-3800,3800);
	histos_1d["L0_dr_dz_range_(grand)mothers"] = new TH1F("L0_dr_dz_range_(grand)mothers",c+"L0_dr_dz_range_(grand)mothers;(grand)mothers pdgid",7600,-3800,3800);
	histos_1d["Ks_dr_dz_range_(grand)mothers"] = new TH1F("Ks_dr_dz_range_(grand)mothers",c+"Ks_dr_dz_range_(grand)mothers;(grand)mothers pdgid",7600,-3800,3800);

	///DELTA PHI, THETA, ETA DISTRIBUTION INVESTIGATION FOR COMPARISON TO DATA sCAND ONES

	histos_1d["L0_Ks_gen_delta_Eta"] = new TH1F("L0_Ks_gen_delta_Eta",c+"L0_Ks_gen_delta_Eta;delta Eta;",100,-5,5);
	histos_1d["L0_Ks_gen_delta_Phi"] = new TH1F("L0_Ks_gen_delta_Phi",c+"L0_Ks_gen_delta_Phi;delta Phi;",100,-5,5);
	histos_1d["L0_Ks_gen_delta_R"] = new TH1F("L0_Ks_gen_delta_R",c+"L0_Ks_gen_delta_R;delta R;",100,0,10);

	histos_1d["L0_Ks_gen_delta_Eta_dCV_within_1mm"] = new TH1F("L0_Ks_gen_delta_Eta_dCV_within_1mm",c+"L0_Ks_gen_delta_Eta_dCV_within_1mm;delta Eta;",100,-5,5);
	histos_1d["L0_Ks_gen_delta_Phi_dCV_within_1mm"] = new TH1F("L0_Ks_gen_delta_Phi_dCV_within_1mm",c+"L0_Ks_gen_delta_Phi_dCV_within_1mm;delta Phi;",100,-5,5);
	histos_1d["L0_Ks_gen_delta_R_dCV_within_1mm"] = new TH1F("L0_Ks_gen_delta_R_dCV_within_1mm",c+"L0_Ks_gen_delta_R_dCV_within_1mm;delta R;",100,0,10);

}

void MCAnalysis(TString category){//"c1", "c2", "c3" or "c123"
	
	///SETTING UP TREE AND HISTOGRAMS

	//TString category = "c123"; //c1, c2, c3, c123

	Int_t Int_category = (category=="c1")? 1 :(category=="c2")? 2 : (category=="c3")? 3 :(category=="c123")? 123 : 0;

	cout <<"category: "<<category<<endl;

    std::map<string, TH1F*> histos_1d;
    std::map<string, TH2F*> histos_2d;
    std::map<string, TEfficiency*> histos_eff;

    bookHistos(histos_1d, histos_2d, histos_eff, category);

    TChain *tree = new TChain("tree/HexaQAnalysis");    
    
	tree->Add("../trees/tree_fixed_GEANT.root");
	//tree->Add("../test/largeskimmedminbias/tree_fixed_GEANT.root"); 
	 
	cout<<"full 10M events minbias + skimmer"<<endl;	 

    HQClass hqhand;
    hqhand.Init(tree);

    Int_t Nentries = hqhand.fChain->GetEntries();
    std::cout<<"Processing "<<Nentries<<" entries"<<std::endl;
    boost::progress_display show_progress( Nentries );
    
    //tree->Print();

    for(Int_t entry = 0; entry < Nentries; ++entry){ //main loop over the events
			
		

		++show_progress;
		//if(entry>20) continue;
		

		hqhand.GetEntry(entry); //hqhand holds all stuff of a particular event
		
		///DECLARATIONS AND DEFINITIONS OF VARIABLES THAT CHANGE EVERY EVENT
		
		
		
		Int_t nGen = hqhand.gen_charge->size(); //number of generated MC particles
		Int_t nLambda=hqhand.nLambda; //number of reconstructed Lambda0s
		Int_t nKshort=hqhand.nKshort; //number of reconstructed Kshorts
		Int_t nTrack=hqhand.nTrack; //number of tracks
		
		nr_reco_L0+=nLambda;
		nr_reco_Ks+=nKshort;
		
		PV.SetXYZ(hqhand.vtx_x->at(0),hqhand.vtx_y->at(0),hqhand.vtx_z->at(0)); //primary vertex position
		
		
		///OVERLAP/DOUBLE COUNTING OF RECO KSHORTS AND LAMBDAS?
		/*
		for(Int_t i=0; i<nLambda; i++){ //loop over reco Lambdas
			for(Int_t j=0; j<nKshort; j++){ //loop over Kshorts
				
				histos_1d["reco_Ks_L0_overlap"]->Fill( //Histogram filled with 0s or 1s
				hqhand.lambda_x->at(i) == hqhand.kshort_x->at(j) &&
				hqhand.lambda_y->at(i) == hqhand.kshort_y->at(j) &&
				hqhand.lambda_z->at(i) == hqhand.kshort_z->at(j) &&
				hqhand.lambda_px->at(i) == hqhand.kshort_px->at(j) &&
				hqhand.lambda_py->at(i) == hqhand.kshort_py->at(j) &&
				hqhand.lambda_pz->at(i) == hqhand.kshort_pz->at(j)
				);
				
			}
		}*/

		///XI(1530) AND RESONANCE CANDIDATE INVESTIGATION

		for(Int_t gen=0; gen<nGen; gen++){

			if(abs(hqhand.gen_pdgid->at(gen))==Xi1530_pdgid) histos_1d["Xi(1530)_0_gen_Eta"]->Fill(hqhand.gen_eta->at(gen));
			if(abs(hqhand.gen_pdgid->at(gen))==Xi1530_pdgid2) histos_1d["Xi(1530)_minus_gen_Eta"]->Fill(hqhand.gen_eta->at(gen));


		}

		///DELTA PHI, THETA, ETA DISTRIBUTION INVESTIGATION FOR COMPARISON TO DATA sCAND ONES

		for(Int_t gen=0; gen<nGen; gen++){

			if(abs(hqhand.gen_pdgid->at(gen))!=Lambda_pdgid) continue;

			Float_t phi_L0 = hqhand.gen_phi->at(gen);
			Float_t eta_L0 = hqhand.gen_phi->at(gen);
			TVector3 CV_L0(hqhand.gen_x->at(gen), hqhand.gen_y->at(gen), hqhand.gen_z->at(gen));

			for(Int_t gen2=0; gen2<nGen; gen2++){

				if(abs(hqhand.gen_pdgid->at(gen2))!=Kshort_pdgid) continue;

				Float_t phi_Ks = hqhand.gen_phi->at(gen2);
				Float_t eta_Ks = hqhand.gen_phi->at(gen2);
				TVector3 CV_Ks(hqhand.gen_x->at(gen2), hqhand.gen_y->at(gen2), hqhand.gen_z->at(gen2));

				histos_1d["L0_Ks_gen_delta_Eta"]->Fill(eta_L0 - eta_Ks);
				histos_1d["L0_Ks_gen_delta_Phi"]->Fill(phi_L0 - phi_Ks);
				histos_1d["L0_Ks_gen_delta_R"]  ->Fill(deltaR(eta_L0, eta_Ks, phi_L0, phi_Ks ));


				TVector3 delta_CV = CV_L0 - CV_Ks;
				if(delta_CV.Mag() > 0.1) continue; //only look at L0 and Ks that have a CV within 1 cm

				histos_1d["L0_Ks_gen_delta_Eta_dCV_within_1mm"] ->Fill(eta_L0 - eta_Ks);
				histos_1d["L0_Ks_gen_delta_Phi_dCV_within_1mm"] ->Fill(phi_L0 - phi_Ks);
				histos_1d["L0_Ks_gen_delta_R_dCV_within_1mm"]   ->Fill(deltaR(eta_L0, eta_Ks, phi_L0, phi_Ks ));




			}

		}
		
		

		
		for(Int_t gen=0; gen<nGen; gen++){
			
			
			///STATUS 8 VS STATUS 1 DOUBLE COUNTING INVESTIGATION
			//Conclusion: there seems to be no double counting
			/*
			if(hqhand.gen_status->at(gen)==8){
				for(Int_t gen2=0; gen2<nGen; gen2++){
					if(hqhand.gen_status->at(gen2)!=1) continue;

					//cout<<hqhand.gen_status->at(gen)<<" "<<hqhand.gen_status->at(gen2)<<endl;

					if(	//hqhand.gen_pdgid->at(gen)==hqhand.gen_pdgid->at(gen2) && 
						hqhand.gen_x->at(gen)==hqhand.gen_x->at(gen2) //&&
						//hqhand.gen_y->at(gen)==hqhand.gen_y->at(gen2) &&
						//hqhand.gen_z->at(gen)==hqhand.gen_z->at(gen2) &&
						//hqhand.gen_px->at(gen)==hqhand.gen_px->at(gen2) &&
						//hqhand.gen_py->at(gen)==hqhand.gen_py->at(gen2) &&
						//hqhand.gen_pz->at(gen)==hqhand.gen_pz->at(gen2)
						) cout<<"double count"<<endl;


				}
			}
			*/



			///Ks AND L0 STATUS 8 & STATUS 1 MOTHER INVESTIGATION
			Int_t mom_index = hqhand.gen_mom->at(gen);
			Int_t mom_pdgid = hqhand.gen_pdgid->at(mom_index);	
			Int_t grandmom_index = hqhand.gen_mom->at(mom_index);
			Int_t grandmom_pdgid = hqhand.gen_pdgid->at(grandmom_index);
			Int_t pdgid = hqhand.gen_pdgid->at(gen);
			Int_t grandmom_mom_pdgid = (abs(mom_pdgid)==abs(pdgid))? grandmom_pdgid:mom_pdgid; //if the mother is the particle itself again, then use the grandmother pdgid

			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			delta_CV_PV = CV - PV;
			Float_t dr = delta_CV_PV.Mag();
			if(!check_category(grandmom_mom_pdgid, dr, Int_category)) continue; 

			
			if(abs(hqhand.gen_pdgid->at(gen))==Lambda_pdgid){

				if(hqhand.gen_status->at(gen) == 8) histos_1d["status_8_L0_(grand)mother_pdgid"]->Fill(abs(grandmom_mom_pdgid));
				//if(hqhand.gen_status->at(gen) == 1) histos_1d["status_1_L0_(grand)mother_pdgid"]->Fill(abs(grandmom_mom_pdgid));
			
				
			} else if(abs(hqhand.gen_pdgid->at(gen))==Kshort_pdgid){
	
				if(hqhand.gen_status->at(gen) == 8) histos_1d["status_8_Ks_(grand)mother_pdgid"]->Fill(abs(grandmom_mom_pdgid));
				//if(hqhand.gen_status->at(gen) == 1) histos_1d["status_1_Ks_(grand)mother_pdgid"]->Fill(abs(grandmom_mom_pdgid));
			
			}
		}
			
			
		
		
		///RECO LAMBDA & KSHORT Pt, Eta, ... DISTRIBUTIONS
		
		histos_1d["L0_reco_amount_per_event"]->Fill(nLambda);
		for(Int_t reco=0; reco<nLambda; reco++){
			//Calculate eta, phi and Pt from Px, Py, Pz and mass
			lorentz.SetPxPyPzE( hqhand.lambda_px->at(reco), hqhand.lambda_py->at(reco), hqhand.lambda_pz->at(reco), sqrt(pow(hqhand.lambda_m->at(reco),2)+pow(hqhand.lambda_px->at(reco),2)+pow(hqhand.lambda_py->at(reco),2)+pow(hqhand.lambda_pz->at(reco),2) )); 

			histos_1d["L0_reco_Pt"]	->Fill(lorentz.Pt());  
			histos_1d["L0_reco_Eta"]->Fill(lorentz.Eta());
			histos_1d["L0_reco_rapidity"]->Fill(lorentz.Rapidity());
			
			
		}
		histos_1d["Ks_reco_amount_per_event"]->Fill(nKshort);
		Int_t nKshort_nodoublecount = 0;
		for(Int_t reco=0; reco<nKshort; reco++){
			//Think about double counting of recostructed L0 and Ks --> Veto Ks, not L0
			Int_t double_count=0;
			for(Int_t reco2=0; reco2<nLambda; reco2++){
				if(hqhand.lambda_x->at(reco2)==hqhand.kshort_x->at(reco)) double_count=1;
			}
			if(double_count==1) continue; //veto Kshorts that have already been reconstructed as Lambda0s
			nKshort_nodoublecount++;
			lorentz.SetPxPyPzE( hqhand.kshort_px->at(reco), hqhand.kshort_py->at(reco), hqhand.kshort_pz->at(reco), sqrt(pow(hqhand.kshort_m->at(reco),2)+pow(hqhand.kshort_px->at(reco),2)+pow(hqhand.kshort_py->at(reco),2)+pow(hqhand.kshort_pz->at(reco),2) )); 

			histos_1d["Ks_reco_Pt"]	->Fill(lorentz.Pt()); 
			histos_1d["Ks_reco_Eta"]->Fill(lorentz.Eta());
			histos_1d["Ks_reco_rapidity"]->Fill(lorentz.Rapidity());

		}
		histos_1d["Ks_reco_amount_per_event_nodoublecount"]->Fill(nKshort_nodoublecount);


		

		
		
		///LAMBDA DECAY POSITION AND RECO EFFICIENCY
		
		
		for(Int_t gen=0; gen<nGen; gen++){ //loop over gen particles
				
			//mother, grandmother and daughters' indices and pdgids
			Int_t d1_index = hqhand.gen_d1->at(gen); //d1 & d2 are daughters of a given particle
			Int_t d2_index = hqhand.gen_d2->at(gen); //index = index of a particle in the vectors of gen particle properties
			Int_t mom_index = hqhand.gen_mom->at(gen);
			Int_t d1_pdgid = hqhand.gen_pdgid->at(d1_index);
			Int_t d2_pdgid = hqhand.gen_pdgid->at(d2_index);
			Int_t mom_pdgid = hqhand.gen_pdgid->at(mom_index);		
			Int_t grandmom_index = hqhand.gen_mom->at(mom_index); //grandmother is the mother of the mother of a given particle
			Int_t grandmom_pdgid = hqhand.gen_pdgid->at(grandmom_index);
			if(abs(mom_pdgid)==Lambda_pdgid) mom_index=grandmom_index; //if the mother is the particle itself again, then look at the grandmother


			Int_t pdgid = hqhand.gen_pdgid->at(gen);
			Int_t grandmom_mom_pdgid = (abs(mom_pdgid)==abs(pdgid))? grandmom_pdgid:mom_pdgid; //if the mother is the particle itself again, then use the grandmother pdgid




			//CHECKING SOME DISTIBUTIONS BEFORE CUTS
			
			if(hqhand.gen_status->at(gen)==8) histos_1d["status_8_pdgid"]->Fill(hqhand.gen_pdgid->at(gen));
			
			
			
			//APPLYING CUTS
			 
			//Only look at gen (anti)Lambda0s
			if(abs(hqhand.gen_pdgid->at(gen))!=Lambda_pdgid) continue; 
			nr_gen_L0++;

			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			delta_CV_PV = CV - PV;
			Float_t dr = delta_CV_PV.Mag();
			if(check_category(grandmom_mom_pdgid, dr, Int_category)) nr_gen_L0_category++;
			
			//must decay to proton+- and Pion+-	as daughters
			if(not(abs(d1_pdgid)==Piplus_pdgid && abs(d2_pdgid)==Proton_pdgid) && not(abs(d2_pdgid)==Piplus_pdgid && abs(d1_pdgid)==Proton_pdgid) ) continue;
			nr_gen_L0_decay++;
			histos_1d["L0_with_correct_daughters_status"]->Fill(hqhand.gen_status->at(gen));

			//must be part of the given category (prompt production, decays of mesons/baryons, material interactions)
			if(!check_category(grandmom_mom_pdgid, dr, Int_category)) continue; 
			nr_gen_L0_decay_category++;

			
			//VERTEX POSITION RELATED STUFF
			//CV = creation vertex position; DV = decay vertex position; PV = primary vertex position; mom=mother; d12 = daughters
			
			momCV.SetXYZ(hqhand.gen_x->at(mom_index),hqhand.gen_y->at(mom_index),hqhand.gen_z->at(mom_index));
			momDV.SetXYZ(hqhand.gen_x->at(hqhand.gen_d1->at(mom_index)),hqhand.gen_y->at(hqhand.gen_d1->at(mom_index)),hqhand.gen_z->at(hqhand.gen_d1->at(mom_index)));
			
			DV.SetXYZ(hqhand.gen_x->at(d1_index),hqhand.gen_y->at(d1_index),hqhand.gen_z->at(d1_index));
			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			d1CV.SetXYZ(hqhand.gen_x->at(d1_index),hqhand.gen_y->at(d1_index),hqhand.gen_z->at(d1_index));
			d2CV.SetXYZ(hqhand.gen_x->at(d2_index),hqhand.gen_y->at(d2_index),hqhand.gen_z->at(d2_index));
			
			delta_DV_PV = DV - PV;
			delta_DV_momCV = DV - momCV;
			delta_CV_momCV = CV - momCV;
			delta_PV_momCV = PV - momCV;
			delta_DV_momDV = DV - momDV;	
			delta_CV_PV = CV - PV;
			
			//PLOTS
			//L0 = Lambda0; gen = it is about gen particles and not reco particles or tracks; 
			
			lorentz2.SetPxPyPzE(hqhand.gen_px->at(gen),hqhand.gen_py->at(gen),hqhand.gen_pz->at(gen),hqhand.gen_energy->at(gen));
						
			histos_1d["L0_gen_dxyz(DV_PV)"]		->Fill(delta_DV_PV.Mag());
			histos_1d["L0_gen_dxy(DV_PV)"]		->Fill(xy_norm(&delta_DV_PV));
			histos_1d["L0_gen_dz(DV_PV)"]		->Fill(delta_DV_PV.Z());
			histos_1d["L0_gen_dxyz(DV_CV)"]	    ->Fill(delta_DV_momDV.Mag());
			histos_1d["L0_gen_dxy(DV_CV)"]	    ->Fill(xy_norm(&delta_DV_momDV));
			histos_1d["L0_gen_dz(DV_CV)"]	    ->Fill(delta_DV_momDV.Z());
			histos_1d["L0_gen_dxyz(CV_PV)"]		->Fill(delta_CV_PV.Mag());
			histos_1d["L0_gen_dxy(CV_PV)"]		->Fill(xy_norm(&delta_CV_PV));
			histos_1d["L0_gen_dz(CV_PV)"]		->Fill(delta_CV_PV.Z());

			histos_1d["L0_gen_dxyz(DV_momCV)"]	->Fill(delta_DV_momCV.Mag());
			histos_1d["L0_gen_dxy(DV_momCV)"]	->Fill(xy_norm(&delta_DV_momCV));
			histos_1d["L0_gen_dz(DV_momCV)"]	->Fill(delta_DV_momCV.Z());
			histos_1d["L0_gen_dxyz(CV_momCV)"]	->Fill(delta_CV_momCV.Mag());
			histos_1d["L0_gen_dxy(CV_momCV)"]	->Fill(xy_norm(&delta_CV_momCV));
			histos_1d["L0_gen_dz(CV_momCV)"]	->Fill(delta_CV_momCV.Z());
			histos_1d["L0_gen_dxyz(PV_momCV)"]	->Fill(delta_PV_momCV.Mag());
			histos_1d["L0_gen_dxy(PV_momCV)"]	->Fill(xy_norm(&delta_PV_momCV));
			histos_1d["L0_gen_dz(PV_momCV)"]	->Fill(delta_PV_momCV.Z());

			histos_1d["L0_gen_Pt"]				->Fill(hqhand.gen_pt->at(gen));
			histos_1d["L0_gen_Eta"]				->Fill(hqhand.gen_eta->at(gen));
			histos_1d["L0_gen_rapidity"]		->Fill(lorentz2.Rapidity());
			histos_1d["L0_gen_dxy(CV_origin)"]	->Fill(xy_norm(&CV));
			

			
			if(d1_index !=0) histos_1d["L0_gen_dxy(daughtersCV_origin)"]->Fill(xy_norm(&d1CV));
			if(d2_index !=0) histos_1d["L0_gen_dxy(daughtersCV_origin)"]->Fill(xy_norm(&d2CV));
			
			//PLOTTING OF AND INVESTIGATION IN SCATTER PLOTS
			
			histos_2d["L0_gen_dx(CV_origin)_dy(CV_origin)"] ->Fill(CV.X(), CV.Y()); 
			histos_2d["L0_gen_dr(CV_origin)_dz(CV_origin)"] ->Fill(CV.Z(), xy_norm(&CV)); //dr and dz switched
			
			
			histos_2d["L0_gen_dx(CV_PV)_dy(CV_PV)"] ->Fill(delta_CV_PV.X(), delta_CV_PV.Y()); 
			histos_2d["L0_gen_dr(CV_PV)_dz(CV_PV)"] ->Fill(delta_CV_PV.Z(), xy_norm(&delta_CV_PV)); //dr and dz switched

			histos_2d["L0_gen_dx(CV_momCV)_dy(CV_momCV)"] ->Fill(delta_CV_momCV.X(), delta_CV_momCV.Y()); 
			histos_2d["L0_gen_dr(CV_momCV)_dz(CV_momCV)"] ->Fill(delta_CV_momCV.Z(), xy_norm(&delta_CV_momCV)); //dr and dz switched
			
			if(abs(delta_CV_PV.Z())<20 && abs(xy_norm(&delta_CV_PV))<20) histos_1d["L0_dr_dz_range_mothers"] ->Fill(hqhand.gen_pdgid->at(hqhand.gen_mom->at(gen)));
			if(abs(delta_CV_PV.Z())<20 && abs(xy_norm(&delta_CV_PV))<20) histos_1d["L0_dr_dz_range_(grand)mothers"] ->Fill(hqhand.gen_pdgid->at(mom_index));


			//ADDITIONAL CUTS ON GEN PARTICLES BEFORE MAKING RECO EFFICIENCY PLOTS
			
			//Cut on eta of the particle or its daughters
			if(abs(hqhand.gen_eta->at(gen))>max_eta || abs(hqhand.gen_eta->at(d1_index))>max_eta || abs(hqhand.gen_eta->at(d2_index))>max_eta) continue;
					
			for(Int_t reco=0; reco<nLambda; reco++){ //loop over reco Lambda0s
				
				
				//Calculate eta, phi and Pt from Px, Py, Pz and mass
				lorentz.SetPxPyPzE( hqhand.lambda_px->at(reco), hqhand.lambda_py->at(reco), hqhand.lambda_pz->at(reco), sqrt(pow(hqhand.lambda_m->at(reco),2)+pow(hqhand.lambda_px->at(reco),2)+pow(hqhand.lambda_py->at(reco),2)+pow(hqhand.lambda_pz->at(reco),2) )); 
				
				reco_eta=lorentz.Eta();
				reco_phi=lorentz.Phi();
				reco_pt=lorentz.Pt();
				
				//RECONSTRUCTION EFFICIENCY PLOTS
				//A recostruction is successful if a reco partcile and gen particle lie withing a cone defined by deltaR < dRmax
				
				dR=deltaR(hqhand.gen_eta->at(gen), reco_eta ,hqhand.gen_phi->at(gen), reco_phi ); //calculate deltaR
				
				histos_eff["L0_recoeff_gen_Pt"]				 ->Fill(dR<dRmax, hqhand.gen_pt->at(gen));
				histos_eff["L0_recoeff_min_gen_Pt_daughters"]->Fill(dR<dRmax, min(hqhand.gen_pt->at(d1_index), hqhand.gen_pt->at(d2_index)));
				
				if(hqhand.gen_pt->at(gen) < 1.2) continue; //Pt cut applies to all recoeff related plots, except for the ones in function of Pt

				histos_1d["L0_gen_reco_deltaR"]->Fill(dR);
				histos_1d["L0_gen_reco_deltaR_zoom"]->Fill(dR);
				
				histos_eff["L0_recoeff_gen_dxyz(DV_PV)"]	->Fill(dR<dRmax, delta_DV_PV.Mag());
				histos_eff["L0_recoeff_gen_dxy(DV_PV)"]		->Fill(dR<dRmax, xy_norm(&delta_DV_PV));
				histos_eff["L0_recoeff_gen_dz(DV_PV)"]		->Fill(dR<dRmax, delta_DV_PV.Z());
				histos_eff["L0_recoeff_gen_dxyz(CV_PV)"]	->Fill(dR<dRmax, delta_CV_PV.Mag());
				histos_eff["L0_recoeff_gen_dxy(CV_PV)"]		->Fill(dR<dRmax, xy_norm(&delta_CV_PV));
				histos_eff["L0_recoeff_gen_dz(CV_PV)"]		->Fill(dR<dRmax, delta_CV_PV.Z());
				histos_eff["L0_recoeff_gen_dxyz(DV_CV)"]	->Fill(dR<dRmax, delta_DV_momDV.Mag());
				histos_eff["L0_recoeff_gen_dxy(DV_CV)"]	    ->Fill(dR<dRmax, xy_norm(&delta_DV_momDV));
				histos_eff["L0_recoeff_gen_dz(DV_CV)"]	    ->Fill(dR<dRmax, delta_DV_momDV.Z());

				histos_eff["L0_recoeff_gen_dxyz(DV_momCV)"]	->Fill(dR<dRmax, delta_DV_momCV.Mag());
				histos_eff["L0_recoeff_gen_dxy(DV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_DV_momCV));
				histos_eff["L0_recoeff_gen_dz(DV_momCV)"]	->Fill(dR<dRmax, delta_DV_momCV.Z());
				histos_eff["L0_recoeff_gen_dxyz(PV_momCV)"]	->Fill(dR<dRmax, delta_PV_momCV.Mag());
				histos_eff["L0_recoeff_gen_dxy(PV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_PV_momCV));
				histos_eff["L0_recoeff_gen_dz(PV_momCV)"]	->Fill(dR<dRmax, delta_PV_momCV.Z());
				histos_eff["L0_recoeff_gen_dxyz(CV_momCV)"]	->Fill(dR<dRmax, delta_CV_momCV.Mag());
				histos_eff["L0_recoeff_gen_dxy(CV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_CV_momCV));
				histos_eff["L0_recoeff_gen_dz(CV_momCV)"]	->Fill(dR<dRmax, delta_CV_momCV.Z());

				histos_eff["L0_recoeff_gen_Eta"]			->Fill(dR<dRmax, hqhand.gen_eta->at(gen));
				histos_eff["L0_recoeff_gen_rapidity"]		->Fill(dR<dRmax, lorentz2.Rapidity());
				histos_eff["L0_recoeff_gen_dxy(CV_origin)"]	->Fill(dR<dRmax, xy_norm(&CV));
				
			 
				
			}
					
		}
		
		if(entry == Nentries - 1){
			//extract the 'numerator' and 'denominator' from the Pt and Eta efficiency plots
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_Pt"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["L0_recoeff_passed_gen_Pt"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_passed_gen_Pt"]->GetTitle());
			temp->SetYTitle("");
			histos_1d["L0_recoeff_passed_gen_Pt"] = temp;
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_Eta"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["L0_recoeff_passed_gen_Eta"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_passed_gen_Eta"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["L0_recoeff_passed_gen_Eta"] = temp;
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_rapidity"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["L0_recoeff_passed_gen_rapidity"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_passed_gen_rapidity"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["L0_recoeff_passed_gen_rapidity"] = temp;
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_Pt"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["L0_recoeff_total_gen_Pt"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_total_gen_Pt"]->GetTitle());
			temp->SetYTitle("");
			histos_1d["L0_recoeff_total_gen_Pt"] = temp;
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_Eta"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["L0_recoeff_total_gen_Eta"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_total_gen_Eta"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["L0_recoeff_total_gen_Eta"] = temp;
			
			temp = (TH1F *) histos_eff["L0_recoeff_gen_rapidity"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["L0_recoeff_total_gen_rapidity"]->GetName());
			temp->SetTitle(histos_1d["L0_recoeff_total_gen_rapidity"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["L0_recoeff_total_gen_rapidity"] = temp;
		}



		
		///KSHORT DECAY POSITION AND RECO EFFICIENCY
		
		for(Int_t gen=0; gen<nGen; gen++){ //loop over gen particles
			
			//mother, grandmother and daughters' indices and pdgids
			Int_t d1_index = hqhand.gen_d1->at(gen);
			Int_t d2_index = hqhand.gen_d2->at(gen);
			Int_t mom_index = hqhand.gen_mom->at(gen);
			Int_t d1_pdgid = hqhand.gen_pdgid->at(d1_index);
			Int_t d2_pdgid = hqhand.gen_pdgid->at(d2_index);
			Int_t mom_pdgid = hqhand.gen_pdgid->at(mom_index);
			Int_t grandmom_index = hqhand.gen_mom->at(mom_index);
			Int_t grandmom_pdgid = hqhand.gen_pdgid->at(grandmom_index);
			if(abs(mom_pdgid)==Kshort_pdgid) mom_index=grandmom_index; //if the mother is the particle itself, then look at the grandmother

			Int_t pdgid = hqhand.gen_pdgid->at(gen);
			Int_t grandmom_mom_pdgid = (abs(mom_pdgid)==abs(pdgid))? grandmom_pdgid:mom_pdgid; //if the mother is the particle itself again, then use the grandmother pdgid

			//if((abs(hqhand.gen_pdgid->at(gen))==Kshort_pdgid2) && hqhand.gen_status->at(gen)==8) cout<<" d1_index: "<<d1_index<<" d2_index: "<<d2_index<<" d1: "<<d1_pdgid<<" d2: "<<d2_pdgid<<endl;
			//if(hqhand.gen_status->at(gen)==8 && abs(d1_pdgid)==Piplus_pdgid && abs(d2_pdgid)==Piplus_pdgid && abs(hqhand.gen_pdgid->at(gen))!=310) cout_pdgid_name(hqhand.gen_pdgid->at(gen));
			
			
			//APPLYING CUTS

			//Only look at gen Kshorts of pdgid 310 (and possibly also 311 or 130)
			if((abs(hqhand.gen_pdgid->at(gen))!=Kshort_pdgid)) continue;// && (abs(hqhand.gen_pdgid->at(gen))!=Kshort_pdgid2) && (abs(hqhand.gen_pdgid->at(gen))!=Kshort_pdgid3)) continue; 
			nr_gen_Ks++;	

			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			delta_CV_PV = CV - PV;
			Float_t dr = delta_CV_PV.Mag();
			if(check_category(grandmom_mom_pdgid, dr, Int_category)) nr_gen_Ks_category++;


			//must decay to pi+- pi+- (not pi0 pi0)
			if(not(abs(d1_pdgid)==Piplus_pdgid && abs(d2_pdgid)==Piplus_pdgid) ) continue;
			nr_gen_Ks_decay++;

			histos_1d["Ks_with_correct_daughters_status"]->Fill(hqhand.gen_status->at(gen));

			//must be part of the given category (prompt production, decays of mesons/baryons, material interactions)
			if(!check_category(grandmom_mom_pdgid, dr, Int_category)) continue; 

			nr_gen_Ks_decay_category++;
			
			//VERTEX POSITION RELATED STUFF
			//CV = creation vertex position; DV = decay vertex position; PV = primary vertex position; mom=mother; d12 = daughters
			
			momCV.SetXYZ(hqhand.gen_x->at(mom_index),hqhand.gen_y->at(mom_index),hqhand.gen_z->at(mom_index));
			momDV.SetXYZ(hqhand.gen_x->at(hqhand.gen_d1->at(mom_index)),hqhand.gen_y->at(hqhand.gen_d1->at(mom_index)),hqhand.gen_z->at(hqhand.gen_d1->at(mom_index)));
			
			DV.SetXYZ(hqhand.gen_x->at(d1_index),hqhand.gen_y->at(d1_index),hqhand.gen_z->at(d1_index));
			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			d1CV.SetXYZ(hqhand.gen_x->at(d1_index),hqhand.gen_y->at(d1_index),hqhand.gen_z->at(d1_index));
			d2CV.SetXYZ(hqhand.gen_x->at(d2_index),hqhand.gen_y->at(d2_index),hqhand.gen_z->at(d2_index));
			
			delta_DV_PV = DV - PV;
			delta_DV_momCV = DV - momCV;
			delta_CV_momCV = CV - momCV;
			delta_PV_momCV = PV - momCV;
			delta_DV_momDV = DV - momDV;	
			delta_CV_PV = CV - PV;
			
			//PLOTS
			//Ks = Kshort; gen = it is about gen particles and not reco particles or tracks;
			
			lorentz2.SetPxPyPzE(hqhand.gen_px->at(gen),hqhand.gen_py->at(gen),hqhand.gen_pz->at(gen),hqhand.gen_energy->at(gen));
						
			histos_1d["Ks_gen_dxyz(DV_PV)"]		->Fill(delta_DV_PV.Mag());
			histos_1d["Ks_gen_dxy(DV_PV)"]		->Fill(xy_norm(&delta_DV_PV));
			histos_1d["Ks_gen_dz(DV_PV)"]		->Fill(delta_DV_PV.Z());
			histos_1d["Ks_gen_dxyz(DV_CV)"]	    ->Fill(delta_DV_momDV.Mag());
			histos_1d["Ks_gen_dxy(DV_CV)"]	    ->Fill(xy_norm(&delta_DV_momDV));
			histos_1d["Ks_gen_dz(DV_CV)"]	    ->Fill(delta_DV_momDV.Z());
			histos_1d["Ks_gen_dxyz(CV_PV)"]		->Fill(delta_CV_PV.Mag());
			histos_1d["Ks_gen_dxy(CV_PV)"]		->Fill(xy_norm(&delta_CV_PV));
			histos_1d["Ks_gen_dz(CV_PV)"]		->Fill(delta_CV_PV.Z());

			histos_1d["Ks_gen_dxyz(DV_momCV)"]	->Fill(delta_DV_momCV.Mag());
			histos_1d["Ks_gen_dxy(DV_momCV)"]	->Fill(xy_norm(&delta_DV_momCV));
			histos_1d["Ks_gen_dz(DV_momCV)"]	->Fill(delta_DV_momCV.Z());
			histos_1d["Ks_gen_dxyz(CV_momCV)"]	->Fill(delta_CV_momCV.Mag());
			histos_1d["Ks_gen_dxy(CV_momCV)"]	->Fill(xy_norm(&delta_CV_momCV));
			histos_1d["Ks_gen_dz(CV_momCV)"]	->Fill(delta_CV_momCV.Z());
			histos_1d["Ks_gen_dxyz(PV_momCV)"]	->Fill(delta_PV_momCV.Mag());
			histos_1d["Ks_gen_dxy(PV_momCV)"]	->Fill(xy_norm(&delta_PV_momCV));
			histos_1d["Ks_gen_dz(PV_momCV)"]	->Fill(delta_PV_momCV.Z());

			histos_1d["Ks_gen_Pt"]				->Fill(hqhand.gen_pt->at(gen));
			histos_1d["Ks_gen_Eta"]				->Fill(hqhand.gen_eta->at(gen));
			histos_1d["Ks_gen_rapidity"]		->Fill(lorentz2.Rapidity());
			histos_1d["Ks_gen_dxy(CV_origin)"]	->Fill(xy_norm(&CV));
			
		
						
			if(d1_index !=0) histos_1d["Ks_gen_dxy(daughtersCV_origin)"]->Fill(xy_norm(&d1CV));
			if(d2_index !=0) histos_1d["Ks_gen_dxy(daughtersCV_origin)"]->Fill(xy_norm(&d2CV));
			
			//PLOTTING OF AND INVESTIGATION IN SCATTER PLOTS
			
			histos_2d["Ks_gen_dx(CV_origin)_dy(CV_origin)"] -> Fill(CV.X(), CV.Y()); 
			histos_2d["Ks_gen_dr(CV_origin)_dz(CV_origin)"] -> Fill(CV.Z(), xy_norm(&CV)); //dr and dz switched
			
			histos_2d["Ks_gen_dx(CV_PV)_dy(CV_PV)"] -> Fill(delta_CV_PV.X(), delta_CV_PV.Y()); 
			histos_2d["Ks_gen_dr(CV_PV)_dz(CV_PV)"] -> Fill(delta_CV_PV.Z(), xy_norm(&delta_CV_PV)); //dr and dz switched

			histos_2d["Ks_gen_dx(CV_momCV)_dy(CV_momCV)"] -> Fill(delta_CV_momCV.X(), delta_CV_momCV.Y()); 
			histos_2d["Ks_gen_dr(CV_momCV)_dz(CV_momCV)"] -> Fill(delta_CV_momCV.Z(), xy_norm(&delta_CV_momCV)); //dr and dz switched
			
			if(abs(delta_CV_PV.Z())<20 && abs(xy_norm(&delta_CV_PV))<20) histos_1d["Ks_dr_dz_range_mothers"] ->Fill(hqhand.gen_pdgid->at(hqhand.gen_mom->at(gen)));
			if(abs(delta_CV_PV.Z())<20 && abs(xy_norm(&delta_CV_PV))<20) histos_1d["Ks_dr_dz_range_(grand)mothers"] ->Fill(hqhand.gen_pdgid->at(mom_index));
			
			//ADDITIONAL CUTS ON GEN PARTICLES BEFORE MAKING RECO EFFICIENCY PLOTS
			
			//Cut on eta of the particle or its daughters
			if(abs(hqhand.gen_eta->at(gen))>max_eta || abs(hqhand.gen_eta->at(d1_index))>max_eta || abs(hqhand.gen_eta->at(d2_index))>max_eta) continue;
		
			for(Int_t reco=0; reco<nKshort; reco++){ //loop over reco Kshorts
				
				//Think about double counting of recostructed L0 and Ks --> Veto Ks, not L0
				Int_t double_count=0;
				for(Int_t reco2=0; reco2<nLambda; reco2++){ //loop over reco Lambda0s
					if(hqhand.lambda_x->at(reco2) == hqhand.kshort_x->at(reco)) double_count = 1;
				}
				if(double_count ==1) continue; //veto Kshorts that have already been reconstructed as Lambda0s
				
				
				//Calculate eta, phi and Pt from Px, Py, Pz and mass
				lorentz.SetPxPyPzE( hqhand.kshort_px->at(reco) , hqhand.kshort_py->at(reco), hqhand.kshort_pz->at(reco), sqrt(pow(hqhand.kshort_m->at(reco),2)+pow(hqhand.kshort_px->at(reco),2)+pow(hqhand.kshort_py->at(reco),2)+pow(hqhand.kshort_pz->at(reco),2) )); 
				reco_eta=lorentz.Eta();
				reco_phi=lorentz.Phi();
				reco_pt=lorentz.Pt();
								
				//RECONSTRUCTION EFFICIENCY PLOTS
				//A recostruction is successful if a reco partcile and gen particle lie withing a cone defined by deltaR < dRmax
				
				dR=deltaR(hqhand.gen_eta->at(gen), reco_eta ,hqhand.gen_phi->at(gen), reco_phi ); //calculating deltaR
				
				histos_eff["Ks_recoeff_gen_Pt"]				 ->Fill(dR<dRmax, hqhand.gen_pt->at(gen));
				histos_eff["Ks_recoeff_min_gen_Pt_daughters"]->Fill(dR<dRmax, min(hqhand.gen_pt->at(d1_index), hqhand.gen_pt->at(d2_index)));
				
				if(hqhand.gen_pt->at(gen) < 0.6) continue; //Pt cut applies to all recoeff related plots, except for the ones in function of Pt
				
				histos_1d["Ks_gen_reco_deltaR"]->Fill(dR);
				histos_1d["Ks_gen_reco_deltaR_zoom"]->Fill(dR);
				
				histos_eff["Ks_recoeff_gen_dxyz(DV_PV)"]	->Fill(dR<dRmax, delta_DV_PV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(DV_PV)"]		->Fill(dR<dRmax, xy_norm(&delta_DV_PV));
				histos_eff["Ks_recoeff_gen_dz(DV_PV)"]		->Fill(dR<dRmax, delta_DV_PV.Z());
				histos_eff["Ks_recoeff_gen_dxyz(CV_PV)"]	->Fill(dR<dRmax, delta_CV_PV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(CV_PV)"]		->Fill(dR<dRmax, xy_norm(&delta_CV_PV));
				histos_eff["Ks_recoeff_gen_dz(CV_PV)"]		->Fill(dR<dRmax, delta_CV_PV.Z());
				histos_eff["Ks_recoeff_gen_dxyz(DV_CV)"]	->Fill(dR<dRmax, delta_DV_momDV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(DV_CV)"]	    ->Fill(dR<dRmax, xy_norm(&delta_DV_momDV));
				histos_eff["Ks_recoeff_gen_dz(DV_CV)"]	    ->Fill(dR<dRmax, delta_DV_momDV.Z());

				histos_eff["Ks_recoeff_gen_dxyz(DV_momCV)"]	->Fill(dR<dRmax, delta_DV_momCV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(DV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_DV_momCV));
				histos_eff["Ks_recoeff_gen_dz(DV_momCV)"]	->Fill(dR<dRmax, delta_DV_momCV.Z());
				histos_eff["Ks_recoeff_gen_dxyz(PV_momCV)"]	->Fill(dR<dRmax, delta_PV_momCV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(PV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_PV_momCV));
				histos_eff["Ks_recoeff_gen_dz(PV_momCV)"]	->Fill(dR<dRmax, delta_PV_momCV.Z());
				histos_eff["Ks_recoeff_gen_dxyz(CV_momCV)"]	->Fill(dR<dRmax, delta_CV_momCV.Mag());
				histos_eff["Ks_recoeff_gen_dxy(CV_momCV)"]	->Fill(dR<dRmax, xy_norm(&delta_CV_momCV));
				histos_eff["Ks_recoeff_gen_dz(CV_momCV)"]	->Fill(dR<dRmax, delta_CV_momCV.Z());

				histos_eff["Ks_recoeff_gen_Eta"]			->Fill(dR<dRmax, hqhand.gen_eta->at(gen));
				histos_eff["Ks_recoeff_gen_rapidity"]		->Fill(dR<dRmax, lorentz2.Rapidity());
				histos_eff["Ks_recoeff_gen_dxy(CV_origin)"]	->Fill(dR<dRmax, xy_norm(&CV));
						 
				
			}
					
		}
		
		if(entry == Nentries - 1){
			//extract the 'numerator' and 'denominator' from the Pt and Eta efficiency plots
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_Pt"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["Ks_recoeff_passed_gen_Pt"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_passed_gen_Pt"]->GetTitle());
			temp->SetYTitle("");
			histos_1d["Ks_recoeff_passed_gen_Pt"] = temp;
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_Eta"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["Ks_recoeff_passed_gen_Eta"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_passed_gen_Eta"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["Ks_recoeff_passed_gen_Eta"] = temp;
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_rapidity"]->GetCopyPassedHisto();
			temp->SetName(histos_1d["Ks_recoeff_passed_gen_rapidity"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_passed_gen_rapidity"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["Ks_recoeff_passed_gen_rapidity"] = temp;
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_Pt"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["Ks_recoeff_total_gen_Pt"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_total_gen_Pt"]->GetTitle());
			temp->SetYTitle("");
			histos_1d["Ks_recoeff_total_gen_Pt"] = temp;
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_Eta"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["Ks_recoeff_total_gen_Eta"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_total_gen_Eta"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["Ks_recoeff_total_gen_Eta"] = temp;
			
			temp = (TH1F *) histos_eff["Ks_recoeff_gen_rapidity"]->GetCopyTotalHisto();
			temp->SetName(histos_1d["Ks_recoeff_total_gen_rapidity"]->GetName());
			temp->SetTitle(histos_1d["Ks_recoeff_total_gen_rapidity"]->GetTitle());
			temp->SetYTitle("");
			temp->SetAxisRange(-3.0, 3.0, "X");
			histos_1d["Ks_recoeff_total_gen_rapidity"] = temp;
		}
				
		
		///LAMBDA KSHORT INV MASS & S MASS
		
		TLorentzVector proton;
		proton.SetPxPyPzE(0.0,0.0,0.0,Mproton); //proton at rest
		
		//Mass for reco Ks and L0
		for(Int_t l=0; l<hqhand.nLambda; l++){ 
			
			TLorentzVector lambda;
			lambda.SetPxPyPzE( hqhand.lambda_px->at(l), hqhand.lambda_py->at(l), hqhand.lambda_pz->at(l), sqrt(pow(hqhand.lambda_m->at(l),2)+pow(hqhand.lambda_px->at(l),2)+pow(hqhand.lambda_py->at(l),2)+pow(hqhand.lambda_pz->at(l),2) )); 

			for(Int_t k=0; k<hqhand.nKshort; k++){ 
				
				if(hqhand.lambda_px->at(l) == hqhand.kshort_px->at(k)) continue; //veto Ks that overlap with L0
				
				TLorentzVector kshort;
				kshort.SetPxPyPzE( hqhand.kshort_px->at(k) , hqhand.kshort_py->at(k), hqhand.kshort_pz->at(k), sqrt(pow(hqhand.kshort_m->at(k),2)+pow(hqhand.kshort_px->at(k),2)+pow(hqhand.kshort_py->at(k),2)+pow(hqhand.kshort_pz->at(k),2) )); 
				
				TLorentzVector sum=kshort+lambda;
				TLorentzVector sum2=kshort+lambda-proton;//sum2.M() = inv mass van S
					
				histos_1d["mass_S_reco"]->Fill(sum2.M());
				histos_1d["inv_mass_Ks_L0_reco"]->Fill(sum.M());
				
			}
		}

		//Mass for gen Ks and L0
		for(Int_t l=0; l<nGen; l++){ 
			
			if(abs(hqhand.gen_pdgid->at(l))!=Lambda_pdgid) continue;
		
			Int_t mom_index = hqhand.gen_mom->at(l);
			Int_t mom_pdgid = hqhand.gen_pdgid->at(mom_index);		
			Int_t grandmom_index = hqhand.gen_mom->at(mom_index); //grandmother is the mother of the mother of a given particle
			Int_t grandmom_pdgid = hqhand.gen_pdgid->at(grandmom_index);
			Int_t pdgid = hqhand.gen_pdgid->at(l);
			Int_t grandmom_mom_pdgid = (abs(mom_pdgid)==abs(pdgid))? grandmom_pdgid:mom_pdgid; //if the mother is the particle itself again, then use the grandmother pdgid

			Int_t gen=l;
			CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
			delta_CV_PV = CV - PV;
			Float_t dr = delta_CV_PV.Mag();

			Int_t L0isC3 = 0; //In case of category 3: only one of the two (Ks, L0) has to be of that category, otherwise no statistics, since chance of both belonging to c3 is tiny
			if(!check_category(grandmom_mom_pdgid, dr, Int_category)){
				if(Int_category!=3) continue; 
			}
			if(get_category(grandmom_mom_pdgid, dr)==3) L0isC3=1;


			
			TLorentzVector lambda;
			lambda.SetPxPyPzE( hqhand.gen_px->at(l) , hqhand.gen_py->at(l), hqhand.gen_pz->at(l), hqhand.gen_energy->at(l)); 

			for(Int_t k=0; k<nGen; k++){ 
				
				if(abs(hqhand.gen_pdgid->at(k))!=Kshort_pdgid) continue;


				mom_index = hqhand.gen_mom->at(k);
				mom_pdgid = hqhand.gen_pdgid->at(mom_index);		
				grandmom_index = hqhand.gen_mom->at(mom_index); //grandmother is the mother of the mother of a given particle
				grandmom_pdgid = hqhand.gen_pdgid->at(grandmom_index);
				pdgid = hqhand.gen_pdgid->at(k);
				grandmom_mom_pdgid = (abs(mom_pdgid)==abs(pdgid))? grandmom_pdgid:mom_pdgid; //if the mother is the particle itself again, then use the grandmother pdgid

				Int_t gen=k;
				CV.SetXYZ(hqhand.gen_x->at(gen),hqhand.gen_y->at(gen),hqhand.gen_z->at(gen));
				delta_CV_PV = CV - PV;
				Float_t dr = delta_CV_PV.Mag();

				Int_t KsisC3 = 0;
				if(!check_category(grandmom_mom_pdgid, dr, Int_category)){
					if(Int_category!=3) continue; 
				}
				if(get_category(grandmom_mom_pdgid, dr)==3) KsisC3=1;

				if(Int_category==3 && !KsisC3 && !L0isC3) continue;

				TLorentzVector kshort;
				kshort.SetPxPyPzE( hqhand.gen_px->at(k) , hqhand.gen_py->at(k), hqhand.gen_pz->at(k), hqhand.gen_energy->at(k)); 
				
				TLorentzVector sum=kshort+lambda;
				TLorentzVector sum2=kshort+lambda-proton;//sum2.M() = inv mass van S
					
				histos_1d["mass_S_gen"]->Fill(sum2.M());
				histos_1d["inv_mass_Ks_L0_gen"]->Fill(sum.M());
				
			}
		}
		
		




		///LAMBDA KSHORT EXTRAPOLATION CLOSEST APPROACH TO EACHOTHER
		/*
		//Extrapolate Lambda and Kshort in direction of their unit momentum vector to their point of closest appraoch to eachother
		//First answer in: https://math.stackexchange.com/questions/538958/finding-coordinates-of-closest-approach    
		
		TMatrixD X(2,1); //A*X=B --> X=A^-1 * B
		TMatrixD A(2,2);
		TMatrixD B(2,1); //(rows,cols)
	

		

		if(hqhand.nLambda!=0 && hqhand.nKshort!=0){
			
			TVector3 vtx_position(hqhand.vtx_x->at(0),hqhand.vtx_y->at(0),hqhand.vtx_z->at(0));	
			Float_t dist_ptcl1_ptcl2;  //this gets minimized
			Float_t dist_ptcl1_ptcl2_prev=100000000000.0;
			Float_t min_dist_ptcl1_ptcl2=1000000000000.0;

			TVector3 delta_pos;
			TVector3 delta_pos_initial;
			TVector3 delta_pos_at_closest_appr;
			TVector3 average_pos_at_closest_appr; //average position of lambda and kshort at closest approach
			TVector3 delta_avrg_pos_and_prim_vtx;
			


			for(Int_t l=0; l<hqhand.nLambda; l++){ //loop over the lambdas
				
				
				TVector3 ptcl1_position(hqhand.lambda_x->at(l),hqhand.lambda_y->at(l),hqhand.lambda_z->at(l));
				TVector3 ptcl1_momentum(hqhand.lambda_px->at(l),hqhand.lambda_py->at(l),hqhand.lambda_pz->at(l));
				TVector3 ptcl1_momentum_unit=ptcl1_momentum.Unit();
				TVector3 ptcl1_pos_extrap;
				TVector3 ptcl1_pos_at_closest_appr;

				for(Int_t k=0; k<hqhand.nKshort; k++){ //loop over the kshorts
					
					
					
					TVector3 ptcl2_position(hqhand.kshort_x->at(k),hqhand.kshort_y->at(k),hqhand.kshort_z->at(k));
					TVector3 ptcl2_momentum(hqhand.kshort_px->at(k),hqhand.kshort_py->at(k),hqhand.kshort_pz->at(k));
					TVector3 ptcl2_momentum_unit=ptcl2_momentum.Unit();
					TVector3 ptcl2_pos_extrap;
					TVector3 ptcl2_pos_at_closest_appr;
					
					
					if(ptcl1_position==ptcl2_position) continue; //Ignore double counting of Lambdas and Kshorts
					
					
					delta_pos_initial=ptcl1_position-ptcl2_position;

					A(0,0)=ptcl1_momentum_unit.Mag2();
					A(1,1)=-ptcl2_momentum_unit.Mag2();
					A(0,1)=-ptcl1_momentum_unit.Dot(ptcl2_momentum_unit);
					A(1,0)=ptcl1_momentum_unit.Dot(ptcl2_momentum_unit);
					B(0,0)=-delta_pos_initial.Dot(ptcl1_momentum_unit);
					B(1,0)=-delta_pos_initial.Dot(ptcl2_momentum_unit);
					
					//A.Print();
					
					A.Invert();
					
					X=A*B;
					
					ptcl1_pos_at_closest_appr=ptcl1_position+X(0,0)*ptcl1_momentum_unit;
					ptcl2_pos_at_closest_appr=ptcl2_position+X(1,0)*ptcl2_momentum_unit;
					delta_pos_at_closest_appr=ptcl1_pos_at_closest_appr-ptcl2_pos_at_closest_appr;
					average_pos_at_closest_appr=(ptcl1_pos_at_closest_appr+ptcl2_pos_at_closest_appr)*0.5;
					
					

					
							
					//~ for(Int_t alfa=-10000; alfa<0; alfa++){ //extrapoleer ptcl1 positie volgens momentum unitvector
						//~ 
						//~ ptcl1_pos_extrap=ptcl1_position + (float(alfa)/50.0) * ptcl1_momentum_unit;
						//~ 
						//~ for(Int_t beta=-10000; beta<0; beta++){ //extrapoleer ptcl2 positie volgens momentum unitvector
						//~ 
							//~ ptcl2_pos_extrap=ptcl2_position + (float(beta)/50.0) * ptcl2_momentum_unit;
							//~ 
							//~ 
							//~ delta_pos = ptcl2_pos_extrap - ptcl1_pos_extrap;
							//~ 
							//~ dist_ptcl1_ptcl2=delta_pos.Mag(); 
							//~ 
							//~ 
							//~ 
							//~ if(dist_ptcl1_ptcl2<min_dist_ptcl1_ptcl2){ //at the closest approach up until now
							//~ 
								//~ 
								//~ 
								//~ average_pos_at_closest_appr=(ptcl1_pos_extrap+ptcl2_pos_extrap)*0.5;
								//~ 
								//~ delta_pos_at_closest_appr=delta_pos;
								//~ 
								//~ min_dist_ptcl1_ptcl2 = dist_ptcl1_ptcl2;
								//~ //cout<<"  "<<min_dist_ptcl1_ptcl2;
							//~ }
							//~ 
							//~ dist_ptcl1_ptcl2_prev = dist_ptcl1_ptcl2;
							//~ 
							//~ 
						//~ }									
					//~ }
					
					//cout <<endl<<"delta_pos_at_closest_appr.Mag(): " <<delta_pos_at_closest_appr.Mag()<<endl;	

					delta_avrg_pos_and_prim_vtx = average_pos_at_closest_appr - vtx_position; //separation vector of average pos to prim vtx

					histos_1d["lambda_kshort_dist_at_closest_appr"]->Fill(delta_pos_at_closest_appr.Mag()); //distance between extrapolated kshort and lambda at their closest approach
					histos_1d["lambda_kshort_delta_z_at_closest_appr"]->Fill(delta_pos_at_closest_appr.Z()); //z-coord distance between extrapolated lambda and kshort at closest appr
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_delta_z"]->Fill(delta_avrg_pos_and_prim_vtx.Z()); //z-coord distance of average of extrap lambda and kshort position at closest appr to prim vtx
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_dist"]->Fill(delta_avrg_pos_and_prim_vtx.Mag()); //distance of average of extrap lambda and kshort position at closest appr to prim vtx
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_impact_param"]->Fill(sqrt(pow(delta_avrg_pos_and_prim_vtx.X(),2)+pow(delta_avrg_pos_and_prim_vtx.Y(),2))); //(X,Y)-coord distance of average of extrap lambda and kshort position at closest appr to prim vtx

					histos_1d["lambda_kshort_dist_at_closest_appr_zoom"]->Fill(delta_pos_at_closest_appr.Mag()); //distance between extrapolated kshort and lambda at their closest approach
					histos_1d["lambda_kshort_delta_z_at_closest_appr_zoom"]->Fill(delta_pos_at_closest_appr.Z()); //z-coord distance between extrapolated lambda and kshort at closest appr
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_delta_z_zoom"]->Fill(delta_avrg_pos_and_prim_vtx.Z()); //z-coord distance of average of extrap lambda and kshort position at closest appr to prim vtx
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_dist_zoom"]->Fill(delta_avrg_pos_and_prim_vtx.Mag()); //distance of average of extrap lambda and kshort position at closest appr to prim vtx
					histos_1d["lambda_kshort_avrg_pos_to_prim_vtx_impact_param_zoom"]->Fill(sqrt(pow(delta_avrg_pos_and_prim_vtx.X(),2)+pow(delta_avrg_pos_and_prim_vtx.Y(),2))); //(X,Y)-coord distance of average of extrap lambda and kshort position at closest appr to prim vtx

				

				}								
			}
		}

		*/

		
		///LAMBDA POSITION EXTRAPOLATION TO VTX
		
		/*
		if(hqhand.nLambda!=0 && hqhand.nKshort!=0){
		for(Int_t l=0; l<hqhand.nLambda; l++){ //loop over the lambdas
			

			histos_1d["inv_mass_Lambda"]->Fill(hqhand.lambda_m->at(l));
			
			
			TVector3 ptcl_position(hqhand.lambda_x->at(l),hqhand.lambda_y->at(l),hqhand.lambda_z->at(l));
			TVector3 vtx_position(hqhand.vtx_x->at(0),hqhand.vtx_y->at(0),hqhand.vtx_z->at(0));	
			TVector3 delta_pos;
			TVector3 ptcl_momentum(hqhand.lambda_px->at(l),hqhand.lambda_py->at(l),hqhand.lambda_pz->at(l));
			TVector3 ptcl_momentum_unit=ptcl_momentum.Unit();
			TVector3 ptcl_pos_extrap;
			
			Float_t dist_to_vtx; //this gets minimized
			Float_t dist_to_vtx_prev=100000000000.0;
			Float_t min_dist_to_vtx=0;//its minimum value
			Float_t ptcl_impact_param=-1;

						
			for(Int_t alfa=-4000; alfa<0; alfa++){ //extrapoleer positie volgens momentum unitvector
				
				ptcl_pos_extrap=ptcl_position + (float(alfa)/20.0) * ptcl_momentum_unit;
				
				delta_pos = vtx_position - ptcl_pos_extrap;

				dist_to_vtx=delta_pos.Mag(); 
				//if(abs(alfa)%2!=0 || abs(alfa)%3!=0 || abs(alfa)%5!=0) cout << " "<<dist_to_vtx ;
				
				if(dist_to_vtx<dist_to_vtx_prev){ 
					ptcl_impact_param = sqrt(pow(ptcl_pos_extrap.X(),2)+pow(ptcl_pos_extrap.Y(),2)); 
					min_dist_to_vtx = dist_to_vtx;
				}
				
				dist_to_vtx_prev=dist_to_vtx;
	
			
			}	
			

			//cout <<endl<<"min_dist_to_vtx: " << min_dist_to_vtx <<endl;
			

			histos_1d["Lambda_impact_param_closest_appr_vtx"]->Fill(ptcl_impact_param);	
			
		}
		}

		
		///KSHORT POSITION EXTRAPOLATION TO VTX
		
		if(hqhand.nLambda!=0 && hqhand.nKshort!=0){
		for(Int_t l=0; l<hqhand.nKshort; l++){ //loop over the kshorts
			
			histos_1d["inv_mass_Kshort"]->Fill(hqhand.kshort_m->at(l));
			
			
			TVector3 ptcl_position(hqhand.kshort_x->at(l),hqhand.kshort_y->at(l),hqhand.kshort_z->at(l));
			TVector3 vtx_position(hqhand.vtx_x->at(0),hqhand.vtx_y->at(0),hqhand.vtx_z->at(0));	
			TVector3 delta_pos;
			TVector3 ptcl_momentum(hqhand.kshort_px->at(l),hqhand.kshort_py->at(l),hqhand.kshort_pz->at(l));
			TVector3 ptcl_momentum_unit=ptcl_momentum.Unit();
			TVector3 ptcl_pos_extrap;
			
			Float_t dist_to_vtx; //this gets minimized
			Float_t dist_to_vtx_prev=100000000000.0;
			Float_t min_dist_to_vtx=0;//its minimum value
			Float_t ptcl_impact_param=-1;

						
			for(Int_t alfa=-4000; alfa<0; alfa++){ //extrapoleer positie volgens momentum unitvector
				
				ptcl_pos_extrap=ptcl_position + (float(alfa)/20.0) * ptcl_momentum_unit;
				
				delta_pos = vtx_position - ptcl_pos_extrap;

				dist_to_vtx=delta_pos.Mag(); 
				//if(abs(alfa)%2!=0 || abs(alfa)%3!=0 || abs(alfa)%5!=0) cout << " "<<dist_to_vtx ;
				
				if(dist_to_vtx<dist_to_vtx_prev){ 
					ptcl_impact_param = sqrt(pow(ptcl_pos_extrap.X(),2)+pow(ptcl_pos_extrap.Y(),2)); 
					min_dist_to_vtx = dist_to_vtx;
				}
				
				dist_to_vtx_prev=dist_to_vtx;
	
			
			}	
			

			//cout <<endl<<"min_dist_to_vtx: " << min_dist_to_vtx <<endl;
			

			histos_1d["Kshort_impact_param_closest_appr_vtx"]->Fill(ptcl_impact_param);	
			
		}
		}
		*/		
				
		///LAMBDA POSITION EXTRAP TO ORIGIN & INV MASS Xi+-
		
		/*
		Int_t Xi_present=0; //First loop through all gen particles to see if Xi is present

		if(hqhand.nLambda!=0 && Xi_present==1){ //if there is at least one lambda and gen Xi+-
			for(Int_t l=0; l<hqhand.nLambda; l++){ //loop over the lambdas
				
				

				histos_1d["Lambda_z"]->Fill(hqhand.lambda_z->at(l));

				//vind lambda_z waarvoor het lambda het dichtste bij de beamline is, m.a.w. waarvoor |(lambda_x, lambda_y)| minimaal is.

				TVector3 lambda_position(hqhand.lambda_x->at(l),hqhand.lambda_y->at(l),hqhand.lambda_z->at(l));	
				TVector3 lambda_momentum(hqhand.lambda_px->at(l),hqhand.lambda_py->at(l),hqhand.lambda_pz->at(l));
				TVector3 lambda_momentum_unit=lambda_momentum.Unit();
				TVector3 lambda_pos_extrap;
				Float_t norm_x_y; //this gets minimized
				Float_t norm_x_y_prev=100000000000.0;
				Float_t min_z=1.0; //resulting z for which it is minimized
				Float_t min_norm_x_y=0;//impact parameter = closest approach distance to beamline

							
				for(Int_t alfa=-3000; alfa<0; alfa++){ //extrapoleer positie van lambda volgens momentum unitvector
					
					lambda_pos_extrap=lambda_position + (float(alfa)/20.0) * lambda_momentum_unit;

					norm_x_y=sqrt(pow(lambda_pos_extrap.X(),2)+pow(lambda_pos_extrap.Y(),2));
					
					if(norm_x_y<norm_x_y_prev){ min_z=lambda_pos_extrap.Z(); min_norm_x_y = norm_x_y;}
					
					norm_x_y_prev=norm_x_y;
		
				
				}	
				

				//cout <<"min_norm_x_y: " << min_norm_x_y <<endl;
				//cout << "min_z: " << min_z<<endl;

				histos_1d["Lambda_z_closest_approach_beam"]->Fill(min_z);
				histos_1d["Lambda_impact_param"]->Fill(min_norm_x_y);		

				
				//put the lambda into a TLorentzVector
				lambda.SetPxPyPzE( hqhand.lambda_px->at(l) , hqhand.lambda_py->at(l), hqhand.lambda_pz->at(l), sqrt(pow(hqhand.lambda_m->at(l),2)+pow(hqhand.lambda_px->at(l),2)+pow(hqhand.lambda_py->at(l),2)+pow(hqhand.lambda_pz->at(l),2) )); 

				
				for(Int_t trk=0; trk<hqhand.nTrack; trk++){ //Loop over the tracks
					
					//put the track into a TLorentzVector (0.13957 is the pi+- mass in GeV)
					pion.SetPtEtaPhiM( hqhand.track_pt->at(trk) , hqhand.track_eta->at(trk), hqhand.track_phi->at(trk), 0.13957 ); 
					sum=pion+lambda; //sum the TLorentzVectors
					
					
					Float_t delta_z = abs(min_z-hqhand.track_z->at(trk)); //in cm
					histos_1d["delta_z"]->Fill(delta_z);

					if(hqhand.track_charge->at(trk)==1){ //if track charge=1, assume it is a pi+


						
						if(delta_z<=0.1) histos_1d["inv_mass_Xi_plus_delta_z_1mm"]->Fill(sum.M()); //fill the histo with the invariant mass
						if(delta_z<=0.5) histos_1d["inv_mass_Xi_plus_delta_z_5mm"]->Fill(sum.M());
						
						histos_1d["inv_mass_Xi_plus"]->Fill(sum.M()); //fill the histo with the invariant mass
						histos_1d["inv_mass_Xi_plus_zoom"]->Fill(sum.M());
						histos_1d["Positive_track_z_Lambda_Xi_present"]->Fill(hqhand.track_z->at(trk));
						
					} else if(hqhand.track_charge->at(trk)==-1){ //if track charge=-1, assume it is a pi-

						
						if(delta_z<=0.1) histos_1d["inv_mass_Xi_min_delta_z_1mm"]->Fill(sum.M());
						if(delta_z<=0.5) histos_1d["inv_mass_Xi_min_delta_z_5mm"]->Fill(sum.M());
						
						histos_1d["inv_mass_Xi_min"]->Fill(sum.M());
						histos_1d["inv_mass_Xi_min_zoom"]->Fill(sum.M());
						histos_1d["Negative_track_z_Lambda_Xi_present"]->Fill(hqhand.track_z->at(trk));
						
				
					} 
				}
			}	
		}
		*/
		


    }//for loop over entries/events

    /*
	pdgid_amount_histo(histos_1d["status_8_Ks_(grand)mother_pdgid"], 0);
	pdgid_amount_histo(histos_1d["status_8_L0_(grand)mother_pdgid"], 0);
	pdgid_amount_histo(histos_1d["status_1_Ks_(grand)mother_pdgid"], 0);
	pdgid_amount_histo(histos_1d["status_1_L0_(grand)mother_pdgid"], 0);
	*/
	
	//pdgid_amount_histo(histos_1d["status_1_L0_(grand)mother_pdgid"], 1);
	//pdgid_amount_histo(histos_1d["status_1_L0_(grand)mother_pdgid"], 2);
	
	
    
	cout<<endl;
	cout<<"Nentries: "<<Nentries<<endl;
	cout<<"avrg_nr_gen_L0: "<<nr_gen_L0/Nentries<<endl;
	cout<<"avrg_nr_gen_Ks: "<<nr_gen_Ks/Nentries<<endl;
	cout<<"avrg_nr_gen_L0_decay: "<<nr_gen_L0_decay/Nentries<<endl;
	cout<<"avrg_nr_gen_Ks_decay: "<<nr_gen_Ks_decay/Nentries<<endl;
	cout<<"avrg_nr_gen_L0_decay_category: "<<nr_gen_L0_decay_category/Nentries<<endl;
	cout<<"avrg_nr_gen_Ks_decay_category: "<<nr_gen_Ks_decay_category/Nentries<<endl;
	cout<<"avrg_nr_gen_L0_category: "<<nr_gen_L0_category/Nentries<<endl;
	cout<<"avrg_nr_gen_Ks_category: "<<nr_gen_Ks_category/Nentries<<endl;
	cout<<"avrg_nr_reco_L0: "<<nr_reco_L0/Nentries<<endl;
	cout<<"avrg_nr_reco_Ks: "<<nr_reco_Ks/Nentries<<endl;
	cout<<endl;
	cout<<"nr_gen_L0: "<<nr_gen_L0<<endl;
	cout<<"nr_gen_Ks: "<<nr_gen_Ks<<endl;
	cout<<"nr_gen_L0_decay: "<<nr_gen_L0_decay<<endl;
	cout<<"nr_gen_Ks_decay: "<<nr_gen_Ks_decay<<endl;
	cout<<"nr_gen_L0_decay_category: "<<nr_gen_L0_decay_category<<endl;
	cout<<"nr_gen_Ks_decay_category: "<<nr_gen_Ks_decay_category<<endl;
	cout<<"nr_gen_L0_category: "<<nr_gen_L0_category<<endl;
	cout<<"nr_gen_Ks_category: "<<nr_gen_Ks_category<<endl;
	cout<<"nr_reco_L0: "<<nr_reco_L0<<endl;
	cout<<"nr_reco_Ks: "<<nr_reco_Ks<<endl;
    
	
	///WRITE PLOTS TO OUTPUT FILE
	
  
	TString outfile = TString("../rootfiles/HQ_plots_minbias_10M_events_skimmed_fixedGEANT")+TString("_")+category+TString(".root");
	TFile * fout = new TFile(outfile, "RECREATE");
	writeHistos(fout, histos_1d, histos_2d, histos_eff);
	fout->Close();
   
    
} 
