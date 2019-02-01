// change done 7 July 2018 -- add >= 1 bjets requirement
#include "ETauAnalysis.h"
using namespace std;

ETauAnalysis::ETauAnalysis(const edm::ParameterSet& iConfig)

{

	isoTau = iConfig.getParameter<string>("isolationTau");
	//	DataHistosForPU = iConfig.getParameter<string>("DataHistos");
	//	MCHistosForPU = iConfig.getParameter<string>("MCHistos");
	//	BTagInputFile = iConfig.getParameter<string>("BTagInputFile");
	InputFile = iConfig.getParameter<string>("InputFile");
	DYSample = iConfig.getParameter<string>("DYSample");                                                                                   

	MuonIDIsoSF = iConfig.getParameter<string>("MuonIDIsoSF");
	MuonTriggerSF = iConfig.getParameter<string>("MuonTriggerSF");
	FakeRateSSfile = iConfig.getParameter<string>("FakeRateSSfile");
	FakeRateDYJetsfile = iConfig.getParameter<string>("FakeRateDYJetsfile");
	size_t npos = -1;
	if(DYSample.find("SE") != npos || DYSample.find("SingleElectron") != npos || DYSample.find("singlephoton") != npos ) { isdata = true; } else { isdata = false;}
	std::cout<<"file:"<<DYSample<<"\tis data \t"<<isdata<<std::endl;
	muscale=1.00;

	TFile *fMuSF = TFile::Open(Form(MuonIDIsoSF.c_str()));
	fhDMuMediumSF = (TH2F*)(fMuSF->Get("hEff_Ele27OR115OR175"));

	delete fMuSF;

	TFile* fake_file = TFile::Open(Form(FakeRateSSfile.c_str()),"R");
	feta_barrel_taupt_300_jetpt_150 = (TH2F*)(fake_file->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_0_150"));
	feta_endcap_taupt_300_jetpt_150 = (TH2F*)(fake_file->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_0_150"));

	feta_barrel_taupt_300_jetpt_300 = (TH2F*)(fake_file->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_150_300"));
	feta_endcap_taupt_300_jetpt_300 = (TH2F*)(fake_file->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_150_300"));

	feta_barrel_taupt_300_jetpt_1000 = (TH2F*)(fake_file->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_300_1000"));
	feta_endcap_taupt_300_jetpt_1000 = (TH2F*)(fake_file->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_300_1000"));

	feta_barrel_taupt_1000_jetpt_300 = (TH2F*)(fake_file->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_300_1000_jetpt_0_300"));
	feta_endcap_taupt_1000_jetpt_300 = (TH2F*)(fake_file->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_300_1000_jetpt_0_300"));

	feta_barrel_taupt_1000_jetpt_1000 = (TH2F*)(fake_file->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_300_1000_jetpt_300_1000"));
	feta_endcap_taupt_1000_jetpt_1000 = (TH2F*)(fake_file->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_300_1000_jetpt_300_1000"));


	delete fake_file;

	TFile* fake_file_DY = TFile::Open(Form(FakeRateDYJetsfile.c_str()),"R");

	DY_feta_barrel_taupt_300_jetpt_150 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_0_150"));
	DY_feta_endcap_taupt_300_jetpt_150 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_0_150"));

	DY_feta_barrel_taupt_300_jetpt_300 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_150_300"));
	DY_feta_endcap_taupt_300_jetpt_300 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_150_300"));

	DY_feta_barrel_taupt_300_jetpt_1000 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_0_300_jetpt_300_1000"));
	DY_feta_endcap_taupt_300_jetpt_1000 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_0_300_jetpt_300_1000"));

	DY_feta_barrel_taupt_1000_jetpt_300 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_300_1000_jetpt_0_300"));
	DY_feta_endcap_taupt_1000_jetpt_300 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_300_1000_jetpt_0_300"));

	DY_feta_barrel_taupt_1000_jetpt_1000 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_barrel_taupt_jetpt_pass_taupt_300_1000_jetpt_300_1000"));
	DY_feta_endcap_taupt_1000_jetpt_1000 = (TH2F*)(fake_file_DY->Get("hratio_data_eta_endcap_taupt_jetpt_pass_taupt_300_1000_jetpt_300_1000"));


	isOS = iConfig.getParameter<bool>("isOS");
	isSS = iConfig.getParameter<bool>("isSS");

	std::cout<<"reached end:"<<std::endl;

	delete fake_file_DY;
}


ETauAnalysis::~ETauAnalysis()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
ETauAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;
	nEventsRaw = nEventsStored = mc_nEventsWeighted = nEventsiihe = 0;
	//	while(!file_db1.eof()) {
	//
	//		file_db1>>datafile;
	//		if(file_db1.eof()) break;
	std::cout<<"datafile::"<<InputFile<<std::endl;
	TFile *file_in=TFile::Open(InputFile.c_str(),"READ");

	TTree* treePtr = (TTree*) file_in->Get("IIHEAnalysis"); 
	//          evCounter = (TH1F*) file_in->Get("HTauTauTree/Counters");
	tree = new IIHEAnalysis (treePtr);   

	std::cout << "entries in tree : "<< tree->GetEntries() << std::endl;
	for (int iEntry = 0; iEntry < tree->GetEntries(); iEntry++)
	{     
		//		cout<<"Entries"<<iEntry<<"\t";                               
		tree->GetEntry(iEntry);
		bool trigger_ele(false), trigger_photon(false);
		h_Events_Before_Skim->Fill(1.);
		// gen filters
		reject_event=false;
		size_t npos = -1;
		Mass=0.;

		gen_tau_had_pt.clear();
		gen_tau_had_eta.clear();
		gen_tau_had_phi.clear();
		gen_tau_had_energy.clear();
		gen_tau_had_charge.clear();

		if(!isdata) {

			for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

				TLorentzVector gen_part, gen_part2, gen_part3;
				gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
				bool isMuon = abs(tree->mc_pdgId->at(iGen))==15  ? true : false ;
				unsigned int  moth_ind = tree->mc_mother_index->at(iGen).at(0);
				bool ishadronicdecay(false);
				int neutrino = 0;
				bool leptau= false;
				if( isMuon) {

					for (unsigned int iGen2 = 0; iGen2 < tree->mc_px->size(); iGen2++){

						gen_part2.SetPxPyPzE(tree->mc_px->at(iGen2),tree->mc_py->at(iGen2),tree->mc_pz->at(iGen2),tree->mc_energy->at(iGen2));
						if(fabs(tree->mc_pdgId->at(iGen2))==11 || fabs(tree->mc_pdgId->at(iGen2))==13) {
							if((tree->mc_mother_index->at(iGen2).at(0)) > 0) {
								if(fabs(tree->mc_pdgId->at(tree->mc_mother_index->at(iGen2).at(0))) == 15 && int(tree->mc_mother_index->at(iGen2).at(0)) == int(iGen) ) {
									leptau=true;
									break;
								}
							}
						}
					}
				}

				if(!(leptau)) {
					for (unsigned int iGen2 = 0; iGen2 < tree->mc_px->size(); iGen2++){
						gen_part3.SetPxPyPzE(tree->mc_px->at(iGen2),tree->mc_py->at(iGen2),tree->mc_pz->at(iGen2),tree->mc_energy->at(iGen2));
						if(fabs(tree->mc_pdgId->at(iGen2))== 16 || fabs(tree->mc_pdgId->at(iGen2))== 12 || fabs(tree->mc_pdgId->at(iGen2))== 14 ) {
							if((tree->mc_mother_index->at(iGen2).at(0)) > 0) {
								if(fabs(tree->mc_pdgId->at(tree->mc_mother_index->at(iGen2).at(0))) == 15. && int(tree->mc_mother_index->at(iGen2).at(0)) == int(iGen) ) {


									neutrino++ ;
								}
							}
						}
					}
				}
				if(neutrino == 1) {
					gen_part = gen_part - gen_part3; 
					gen_tau_had_pt.push_back(gen_part.Pt());
					gen_tau_had_eta.push_back(gen_part.Eta());
					gen_tau_had_phi.push_back(gen_part.Phi());
					gen_tau_had_energy.push_back(gen_part.E());
					gen_tau_had_charge.push_back(tree->mc_charge->at(iGen));
				} else {

					gen_tau_had_pt.push_back(-9);
					gen_tau_had_eta.push_back(-9);
					gen_tau_had_phi.push_back(-9);
					gen_tau_had_energy.push_back(-9);
					gen_tau_had_charge.push_back(-9);



				}
			}



		}

		if(!isdata &&  (DYSample.find("DYinc") != npos || DYSample.find("DY") != npos ||  DYSample.find("TT") != npos || DYSample.find("WW")  != npos || DYSample.find("TTinc") != npos ||  DYSample.find("WWinc") != npos  )) {


			vector<TLorentzVector> lep;
			lep.clear();
			reject_event=false;

			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {
				if ( (fabs(tree->LHE_pdgid->at(iLHE)) == 11 || fabs(tree->LHE_pdgid->at(iLHE)) == 13 || fabs(tree->LHE_pdgid->at(iLHE)) == 15)) {
					TLorentzVector l1_p4 ; 
					l1_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					lep.push_back(l1_p4);
				}
			}
			cout<<"size of lep:"<< lep.size()<< std::endl;
			if( lep.size() ==  2 ) { //continue;
				Mass =  (lep.at(0)+ lep.at(1)).M();
				if( DYSample.find("DYinc") != npos || DYSample.find("DYJets_inc")  != npos ){  if (  (lep.at(0)+ lep.at(1)).M() > 400.) { reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M();}}

				if( DYSample.find("WWinc") != npos || DYSample.find("WW_inc") != npos ) { if(  (lep.at(0)+ lep.at(1)).M() > 200.) {reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M(); }}

				if( DYSample.find("TTinc") != npos || DYSample.find("TTJets_inc") != npos || DYSample.find("TTTo2L2Nu_TuneCUETP8M2_ttHtranche3") != npos) { if(  (lep.at(0)+ lep.at(1)).M() > 500.) {reject_event = true; Mass=  (lep.at(0)+ lep.at(1)).M(); }}
			} else { reject_event = false;}
		} else if(!isdata &&  (DYSample.find("WJets") != npos || DYSample.find("WJetsinc") != npos )) {

			reject_event=false;

			TLorentzVector l_p4, nu_p4, lnu_p4;
			l_p4.SetPxPyPzE(0, 0, 0, 0);
			nu_p4.SetPxPyPzE(0, 0, 0, 0);
			lnu_p4.SetPxPyPzE(0, 0, 0, 0);
			int l_pdgid = 0, nu_pdgid = 0;
			Mass = 0;
			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {

				if (  (fabs(tree->LHE_pdgid->at(iLHE)) == 11 || fabs(tree->LHE_pdgid->at(iLHE)) == 13 || fabs(tree->LHE_pdgid->at(iLHE)) == 15))  {
					l_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					l_pdgid = tree->LHE_pdgid->at(iLHE);
				}
				else if (abs(tree->LHE_pdgid->at(iLHE)) == 12 || abs(tree->LHE_pdgid->at(iLHE)) == 14 || abs(tree->LHE_pdgid->at(iLHE)) == 16) {
					nu_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
					nu_pdgid = tree->LHE_pdgid->at(iLHE);
				}
			}
			if(fabs(l_pdgid)+1 == fabs(nu_pdgid)) {
				lnu_p4 = l_p4 + nu_p4;
				Mass=  lnu_p4.Pt();
				if( Mass < 600.) std::cout<<"pt:::"<< Mass << "::SIgn evt weight::"<< tree->mc_w_sign << std::endl;
				if(DYSample.find("WJetsinc") != npos  || DYSample.find("WJetsToLNu_TuneCUETP8M1") != npos) {if (lnu_p4.Pt() > 100) { reject_event = true; Mass=  lnu_p4.Pt();}}

			} else { reject_event = false;}
		} else {reject_event = false; } 

		if(!isdata) h_Fill_Mass_Gen_toChk_before->Fill(Mass);
		h_Events_Before_Skim->Fill(1.);
		if(reject_event) continue;

		if(!isdata) {h_Events_After_GenFilter->Fill(1.);}
		//			std::cout<<"=========================="<<std::endl;
		if(!isdata) h_Fill_Mass_Gen_toChk->Fill(Mass,tree->mc_w_sign);

		if( !isdata) {  h_Count_Taus->Fill(GenTaus(), tree->mc_w_sign); }
		if(DYSample.find("singlephoton") != npos && isdata ) {

			if(!(!tree->trig_HLT_Ele27_WPTight_Gsf_accept && (!tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept) &&  (tree->trig_HLT_Photon175_accept))) continue;
		}
		else if(DYSample.find("SE") != npos && isdata ){
			if( !( tree->trig_HLT_Ele27_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept )) continue;
		} else {
			if( !( tree->trig_HLT_Ele27_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept || tree->trig_HLT_Photon175_accept)) continue;      
		}


		if(isdata) {h_Events_After_GenFilter->Fill(1.);}

		if(!tree->trig_Flag_goodVertices_accept) continue;
		if(!tree->trig_Flag_globalTightHalo2016Filter_accept) continue;
		if(!tree->trig_Flag_HBHENoiseFilter_accept) continue;
		if(!tree->trig_Flag_HBHENoiseIsoFilter_accept) continue;
		if(!tree->trig_Flag_EcalDeadCellTriggerPrimitiveFilter_accept) continue;
		if(!tree->trig_Flag_BadPFMuonFilter_accept) continue;
		if(!tree->trig_Flag_BadChargedCandidateFilter_accept) continue;
		if(isdata) if(!tree->trig_Flag_eeBadScFilter_accept) continue;


		// loop for OS and SS events
		// events filters
		GoodElecs.clear();
		GoodElecsIndex.clear();
		GoodTaus.clear();
		GoodTausIndex.clear();
		vector<pair<TLorentzVector, unsigned int>> ElePairs;
		ElePairs.clear();
		vector<pair<TLorentzVector, unsigned int>> TauPairs;
		TauPairs.clear();

		sumweight =0;
		weighthis=1.;

		if(!isdata) {pu_weight =PU_reReco_Morind17::MC_pileup_weight(tree->mc_trueNumInteractions, 0, "all"); } else {pu_weight = 1;}
		sumweight = sumweight+tree->mc_weight;
		if(!isdata) weighthis = tree->mc_w_sign*pu_weight; 
		if(isdata) { weighthis = 1.;}

		if(!isdata &&  ( DYSample.find("TT") != npos || DYSample.find("TTinc") != npos )) {

			double weight_top1 = 1.0;
			//double weight_top2 = 1.0;
			int tops=0;
			for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {

				if ( fabs(tree->LHE_pdgid->at(iLHE)) == 6) {

					weight_top1 =  weight_top1 * (exp(-1.08872 - 0.0119998* tree->LHE_Pt->at(iLHE)) + 0.895139);
					tops++;
				}

			}
			weighthis = weighthis* sqrt(weight_top1);
//			std::cout<<"how many gen tops:"<< tops<< ": weight:" << weighthis << std::endl;
		}

		double masspair = 0;
		TLorentzVector FirstObj(0.,0.,0.,0), SecondObj(0.,0.,0.,0);
		unsigned int first_index = -1;
		unsigned int second_index = -1;

		TLorentzVector FirstObj_FF(0.,0.,0.,0), SecondObj_FF(0.,0.,0.,0);
		first_index_FF= -1;
		second_index_FF = -1;
		double masspair_FF= 0;
		jet_index_FF =  -1;

		for (unsigned int dau1index = 0; dau1index < tree->gsf_caloEnergy->size(); dau1index++){
			float ET1 = tree->gsf_caloEnergy->at(dau1index)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(dau1index)))) ;


			TLorentzVector DauEle ;
			DauEle.SetPtEtaPhiM(ET1, tree->gsf_eta->at(dau1index), tree->gsf_phi->at(dau1index),m_el);
			if(!( ET1 > 50.  && fabs(tree->gsf_sc_eta->at(dau1index)) < 2.5) ) continue;
			if(!(tree->gsf_VIDHEEP7->at(dau1index))) continue;

			for (unsigned int dau2index = 0; dau2index < tree->tau_pt->size(); dau2index++){

				if(!(tree->tau_pt->at(dau2index) > 30.)) continue;
				if(!(fabs(tree->tau_eta->at(dau2index)) < 2.3)) continue;


				if(!(tree->tau_decayModeFinding->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_byVLooseIsolationMVArun2v1DBoldDMwLT->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstMuonLoose3->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstElectronTightMVA6->at(dau2index) > 0.5)) continue;

				TLorentzVector DauTau ;
				DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));

				if(!(DauEle.DeltaR(DauTau) > 0.5 ))  continue;
				TVector2 METcorr;
				METcorr.SetMagPhi(tree->MET_T1Txy_Pt,tree->MET_T1Txy_phi); 
				double mt_cut = mTCalculation(METcorr.Px(), METcorr.Py(), DauEle.Px(), DauEle.Py(), DauEle.Pt());


				if( mt_cut > 120. ) continue;

				//					if(tree->gsf_charge->at(dau1index) * tree->tau_charge->at(dau2index) != 1 ) continue;
				double applySF= false;
				int genindex = -1;
				unsigned int matchgen_obj_type;
				if(!isdata) { if( (MatchingToGenTaus(DauTau, genindex) && genindex!=-1) || matchedToGenObjetcs(DauTau, matchgen_obj_type)) {applySF = true;} }
				if(isdata) {applySF= true;}
				//				applySF= true;
				if(applySF) {


					if(tree->tau_byTightIsolationMVArun2v1DBoldDMwLT->at(dau2index) > 0.5) {

						if( (DauEle + DauTau).M() > masspair) {
							masspair = (DauEle +DauTau).M();
							FirstObj = DauEle;
							SecondObj = DauTau;
							first_index = dau1index;
							second_index =dau2index;
						}
					} 
					if ( !(tree->tau_byTightIsolationMVArun2v1DBoldDMwLT->at(dau2index) > 0.5)) {
						unsigned int jet_indx = 999;
						if(TauMatchedToJet(DauTau, jet_indx) && jet_indx != 999) {
							if( (DauEle + DauTau).M() > masspair_FF) {
								masspair_FF = (DauEle +DauTau).M();
								FirstObj_FF = DauEle;
								SecondObj_FF = DauTau;
								first_index_FF = dau1index;
								second_index_FF =   dau2index;
								jet_index_FF = jet_indx;
							}
						}
					}

				} // applySF

			}
		}


		/////////// wo matching...... ///
		//



		double weight = weighthis;
		double weight_pass = weighthis;
		TVector2 metV2;
		// use of MET final collection
		// patPFMetT1Txy
		//	metV.SetPxPyPzE(tree->MET_T1Txy_Px, tree->MET_T1Txy_Py, 0. , tree->MET_T1Txy_Pt);
		metV2.SetMagPhi(tree->MET_T1Txy_Pt,tree->MET_T1Txy_phi);;
		TLorentzVector metV;
		metV.SetPxPyPzE(metV2.Px(), metV2.Py(), 0. , tree->MET_T1Txy_Pt);

		if(int(first_index_FF) != -1 && int(second_index_FF) != -1 ) {
			if( (!isdata) )  { weight = weight * GetEfficiency(fabs(FirstObj_FF.Eta()), FirstObj_FF.Pt(), fhDMuMediumSF); }



			if(!isdata) { double fake_wt =  MutoTauFR(SecondObj_FF);
				weight = weight * fake_wt;
			}

			//			double fr_weight_ss_DY = 0;
			//			if(TauMatchedToJet(SecondObj_FF, jet_indx) && jet_indx != 999) {

			double fr_weight_ss_DY = TauPtBinSF_CR5( SecondObj_FF.Pt(), SecondObj_FF.Eta(), SecondObj_FF.Pt()/ tree->jet_pt->at(jet_index_FF)); 


			//			}
			double unc_weight;
			unc_weight=1.;
			double tauSF= 0.95;

			for (unsigned int t=0; t <5; t++) {
				if(t==0) { unc_weight =tauSF;}


				if( t == 1) { unc_weight = (1+ 0.05)*tauSF; }
				if( t == 2) { unc_weight = (1- 0.05)*tauSF; }
				if( t == 3) { unc_weight = (1+ ((0.05*SecondObj.Pt())/1000.))*tauSF; }
				if( t == 4) { unc_weight = (1- ((0.35*SecondObj.Pt())/1000.))*tauSF; }

				if(isdata) unc_weight = 1.;
				//	weight = unc_weight* weight;

				h_Fill_NV_TauFake[t]->Fill(tree->pv_n,weight*unc_weight*fr_weight_ss_DY);
				h_EleTauCharge_TauFake[t]->Fill(tree->gsf_charge->at(first_index_FF) * tree->tau_charge->at(second_index_FF),weight*unc_weight*fr_weight_ss_DY);
				h_FillmT_Ele_TauFake[t]->Fill(mTCalculation(metV.Px(), metV.Py(), FirstObj_FF.Px(), FirstObj_FF.Py(), FirstObj_FF.Pt()),weight*unc_weight*fr_weight_ss_DY); 
				h_FillmT_Tau_TauFake[t]->Fill(mTCalculation(metV.Px(), metV.Py(), SecondObj_FF.Px(), SecondObj_FF.Py(), SecondObj_FF.Pt()),weight*unc_weight*fr_weight_ss_DY); 
				h_Fill_DPhi_Ele_Met_TauFake[t]->Fill(deltaPhi(FirstObj_FF.Phi(), metV.Phi()),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_DPhi_Tau_Met_TauFake[t]->Fill(deltaPhi(SecondObj_FF.Phi(), metV.Phi()),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_DPhi_Ele_Tau_TauFake[t]->Fill(deltaPhi(FirstObj_FF.Phi(), SecondObj_FF.Phi()),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_EleTauMass_TauFake[t]->Fill( (FirstObj_FF + SecondObj_FF).M(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_TauMETMass_TauFake[t]->Fill( ( SecondObj_FF+  TLorentzVector(metV)).M(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_TotalMass_TauFake[t]->Fill( ( FirstObj_FF + SecondObj_FF + metV ).M(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_CollMass_TauFake[t]->Fill(GetCollinearMass(FirstObj_FF,SecondObj_FF,metV), weight*unc_weight*fr_weight_ss_DY);
				h_Fill_Met_TauFake[t]->Fill(metV.Pt(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_taupt_TauFake[t]->Fill(SecondObj_FF.Pt(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_taueta_TauFake[t]->Fill(SecondObj_FF.Eta(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_mupt_TauFake[t]->Fill(tree->gsf_et->at(first_index_FF),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_mupt_SC_TauFake[t]->Fill(FirstObj_FF.Pt(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_mueta_TauFake[t]->Fill(FirstObj_FF.Eta(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_mueta_SC_TauFake[t]->Fill(tree->gsf_sc_eta->at(first_index_FF),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_metphi_TauFake[t]->Fill(metV.Phi(),weight*unc_weight*fr_weight_ss_DY);

				h_Fill_tauphi_TauFake[t]->Fill(SecondObj_FF.Phi(),weight*unc_weight*fr_weight_ss_DY);

				h_Fill_muphi_SC_TauFake[t]->Fill(tree->gsf_sc_phi->at(first_index_FF),weight*unc_weight*fr_weight_ss_DY);

				h_Fill_tauphi_muphi_TauFake[t]->Fill(SecondObj_FF.Phi(),tree->gsf_sc_phi->at(first_index_FF),weight*unc_weight*fr_weight_ss_DY);

				h_Fill_tauphi_metphi_TauFake[t]->Fill(SecondObj_FF.Phi(),metV.Phi(),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_muphi_metphi_TauFake[t]->Fill(tree->gsf_sc_phi->at(first_index_FF),metV.Phi(),weight*unc_weight*fr_weight_ss_DY);


			}
		} // fake rate plots DONE ====


		if(int(first_index) != -1 && int(second_index) != -1 ) {
			if( (!isdata) )  { weight_pass = weight_pass * GetEfficiency(fabs(FirstObj.Eta()), FirstObj.Pt(), fhDMuMediumSF); }

			if(!isdata)  {double fake_wt =  MutoTauFR(SecondObj);
				weight_pass = weight_pass * fake_wt;
			}


			double unc_weight;
			unc_weight=1.;
			double tauSF= 0.95;
			for (unsigned int t=0; t <5; t++) {
				if(t==0) { unc_weight =tauSF;}


				if( t == 1) { unc_weight = (1 + 0.05)*tauSF; }
				if( t == 2) { unc_weight = (1 - 0.05)*tauSF; }

				if( t == 3) { unc_weight = (1 + ((0.05*SecondObj.Pt())/1000.))*tauSF; }
				if( t == 4) { unc_weight = (1 - ((0.35*SecondObj.Pt())/1000.))*tauSF; }

				if(isdata) unc_weight = 1.;

				h_Fill_NV_TauPass[t]->Fill(tree->pv_n,weight_pass*unc_weight);

				h_EleTauCharge_TauPass[t]->Fill(tree->gsf_charge->at(first_index) * tree->tau_charge->at(second_index),weight_pass*unc_weight);
				h_FillmT_Ele_TauPass[t]->Fill(mTCalculation(metV.Px(), metV.Py(), FirstObj.Px(), FirstObj.Py(), FirstObj.Pt()),weight_pass*unc_weight);
				h_FillmT_Tau_TauPass[t]->Fill(mTCalculation(metV.Px(), metV.Py(), SecondObj.Px(), SecondObj.Py(), SecondObj.Pt()),weight_pass*unc_weight);
				h_Fill_DPhi_Ele_Met_TauPass[t]->Fill(deltaPhi(FirstObj.Phi(), metV.Phi()),weight_pass*unc_weight);
				h_Fill_DPhi_Tau_Met_TauPass[t]->Fill(deltaPhi(SecondObj.Phi(), metV.Phi()),weight_pass*unc_weight);
				h_Fill_DPhi_Ele_Tau_TauPass[t]->Fill(deltaPhi(FirstObj.Phi(), SecondObj.Phi()),weight_pass*unc_weight);
				h_Fill_EleTauMass_TauPass[t]->Fill( (FirstObj + SecondObj).M(),weight_pass*unc_weight);
				h_Fill_TauMETMass_TauPass[t]->Fill( ( SecondObj +  metV).M(),weight_pass*unc_weight);
				h_Fill_TotalMass_TauPass[t]->Fill( (FirstObj + SecondObj + metV ).M(),weight_pass*unc_weight);
				h_Fill_CollMass_TauPass[t]->Fill(GetCollinearMass(FirstObj,SecondObj ,metV), weight_pass*unc_weight);
				h_Fill_Met_TauPass[t]->Fill(metV.Pt(),weight_pass*unc_weight);
				h_Fill_taupt_TauPass[t]->Fill(SecondObj.Pt(),weight_pass*unc_weight);
				h_Fill_taueta_TauPass[t]->Fill(SecondObj.Eta(),weight_pass*unc_weight);
				h_Fill_mupt_TauPass[t]->Fill(tree->gsf_et->at(first_index),weight_pass*unc_weight);
				h_Fill_mupt_SC_TauPass[t]->Fill(FirstObj.Pt(),weight_pass*unc_weight);
				h_Fill_mueta_TauPass[t]->Fill(FirstObj.Eta(),weight_pass*unc_weight);
				h_Fill_mueta_SC_TauPass[t]->Fill(tree->gsf_sc_eta->at(first_index),weight_pass*unc_weight);
				h_Fill_metphi_TauPass[t]->Fill(metV.Phi(),weight_pass*unc_weight);

			}
		}





	}






	file_in->Close();
	//	}    
	//	file_db1.close();     



#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif

}
// ------------ method called once each job just before starting event loop  ------------


	void 
ETauAnalysis::beginJob()
{
	//	fIn = TFile::Open(InputFile.c_str());
	std::cout<<"InputFile::"<<InputFile<< std::endl;
	//	file_db1.open(InputFile);                                                                                        


	h_Count_Taus = fs->make<TH1D>("h_Count_Taus","h_Count_Taus",10,0,10); 
	h_Count_Taus->Sumw2();
	h_Events_Before_Skim = fs->make<TH1D>("h_Events_Before_Skim","h_Events_Before_Skim",2,0,2);
	h_Events_After_Skim = fs->make<TH1D>("h_Events_After_Skim","h_Events_After_Skim",2,0,2);
	h_Events_After_GenFilter = fs->make<TH1D>("h_Events_After_GenFilter","h_Events_After_GenFilter",2,0,2);
	h_Fill_Mass_Gen_toChk = fs->make<TH1D>("h_Fill_Mass_Gen_toChk","h_Fill_Mass_Gen_toChk",3000,0,3000);
	h_Fill_Mass_Gen_toChk_before = fs->make<TH1D>("h_Fill_Mass_Gen_toChk_before","h_Fill_Mass_Gen_toChk_before",3000,0,3000);

	h_Fill_Mass_Gen_toChk->Sumw2();
	h_Fill_Mass_Gen_toChk_before->Sumw2();
	h_Events_Before_Skim->Sumw2();
	h_Events_After_Skim->Sumw2();
	h_Events_After_GenFilter->Sumw2();


	//	h_Fill_bjet_pt_ETauFake = fs->make<TH1D>("h_Fill_bjet_pt_ETauFake","h_Fill_bjet_pt_ETauFake", 1000,0,1000);
	//	h_Fill_bjet_eta_ETauFake = fs->make<TH1D>("h_Fill_bjet_eta_ETauFake","h_Fill_bjet_eta_ETauFake",100,-5,5);
	//	h_Fill_bjet_disc_ETauFake = fs->make<TH1D>("h_Fill_bjet_disc_ETauFake","h_Fill_bjet_disc_ETauFake",100,0,1);

	//	h_Fill_bjet_pt_ETauFake->Sumw2();
	//	h_Fill_bjet_eta_ETauFake->Sumw2();
	//	h_Fill_bjet_disc_ETauFake->Sumw2();

	string systName[] = {"Normal","Normal_TauUp","Normal_TauDown","Normal_TauUp_high", "Normal_TauDown_high"};
	int numH = sizeof(systName)/sizeof(string);

	char histname[256], histname2[256];
	for(int j = 0; j <numH; j++){  

		sprintf(histname,"h_EleTauCharge_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_ETauFake_%s",systName[j].c_str());
		h_EleTauCharge_ETauFake[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_ETauFake[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_ETauFake_%s",systName[j].c_str());
		h_FillmT_Ele_ETauFake[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_ETauFake[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_ETauFake_%s",systName[j].c_str());
		h_FillmT_Tau_ETauFake[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_ETauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_ETauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_ETauFake[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_ETauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_ETauFake[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_ETauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_ETauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_ETauFake[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_ETauFake_%s",systName[j].c_str());
		h_Fill_EleTauMass_ETauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_ETauFake_%s",systName[j].c_str());
		h_Fill_TauMETMass_ETauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_ETauFake_%s",systName[j].c_str());
		h_Fill_TotalMass_ETauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_ETauFake_%s",systName[j].c_str());
		h_Fill_Met_ETauFake[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_ETauFake_%s",systName[j].c_str());
		h_Fill_CollMass_ETauFake[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_ETauFake_%s",systName[j].c_str());
		h_Fill_taupt_ETauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_ETauFake_%s",systName[j].c_str());
		h_Fill_taueta_ETauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_ETauFake_%s",systName[j].c_str());
		h_Fill_mupt_ETauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_ETauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_ETauFake_%s",systName[j].c_str());

		h_Fill_mupt_SC_ETauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_ETauFake_%s",systName[j].c_str());
		h_Fill_mueta_ETauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_ETauFake_%s",systName[j].c_str());
		h_Fill_mueta_SC_ETauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_ETauFake_%s",systName[j].c_str());	 
		h_Fill_metphi_ETauFake[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_ETauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_ETauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_ETauFake_%s",systName[j].c_str());	
		h_Fill_NV_ETauFake[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_ETauFake[j]->Sumw2();

		///

		sprintf(histname,"h_EleTauCharge_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauFakeM_%s",systName[j].c_str());
		h_EleTauCharge_TauFakeM[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauFakeM[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauFakeM_%s",systName[j].c_str());
		h_FillmT_Ele_TauFakeM[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauFakeM_%s",systName[j].c_str());
		h_FillmT_Tau_TauFakeM[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauFakeM[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauFakeM_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauFakeM[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauFakeM_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauFakeM[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauFakeM[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauFakeM_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauFakeM[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauFakeM_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauFakeM[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauFakeM_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauFakeM[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauFakeM_%s",systName[j].c_str());
		h_Fill_TotalMass_TauFakeM[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauFakeM_%s",systName[j].c_str());
		h_Fill_Met_TauFakeM[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauFakeM_%s",systName[j].c_str());
		h_Fill_CollMass_TauFakeM[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauFakeM_%s",systName[j].c_str());
		h_Fill_taupt_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauFakeM_%s",systName[j].c_str());
		h_Fill_taueta_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauFakeM_%s",systName[j].c_str());
		h_Fill_mupt_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauFakeM[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauFakeM_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauFakeM_%s",systName[j].c_str());
		h_Fill_mueta_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauFakeM_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauFakeM_%s",systName[j].c_str());	 
		h_Fill_metphi_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauFakeM[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauFakeM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauFakeM_%s",systName[j].c_str());	
		h_Fill_NV_TauFakeM[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauFakeM[j]->Sumw2();




		///

		/// _DY


		sprintf(histname,"h_EleTauCharge_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauFakeM_DY_%s",systName[j].c_str());
		h_EleTauCharge_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauFakeM_DY[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauFakeM_DY_%s",systName[j].c_str());
		h_FillmT_Ele_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauFakeM_DY_%s",systName[j].c_str());
		h_FillmT_Tau_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauFakeM_DY[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauFakeM_DY[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauFakeM_DY[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_TotalMass_TauFakeM_DY[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_Met_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_CollMass_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_taupt_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_taueta_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_mupt_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauFakeM_DY[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauFakeM_DY_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_mueta_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauFakeM_DY_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauFakeM_DY_%s",systName[j].c_str());	 
		h_Fill_metphi_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauFakeM_DY[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauFakeM_DY_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauFakeM_DY_%s",systName[j].c_str());	
		h_Fill_NV_TauFakeM_DY[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauFakeM_DY[j]->Sumw2();




		///


		///

		sprintf(histname,"h_EleTauCharge_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauFake_%s",systName[j].c_str());
		h_EleTauCharge_TauFake[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauFake[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauFake_%s",systName[j].c_str());
		h_FillmT_Ele_TauFake[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauFake[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauFake_%s",systName[j].c_str());
		h_FillmT_Tau_TauFake[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauFake[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauFake[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauFake_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauFake[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauFake_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauFake_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauFake_%s",systName[j].c_str());
		h_Fill_TotalMass_TauFake[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauFake_%s",systName[j].c_str());
		h_Fill_Met_TauFake[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauFake_%s",systName[j].c_str());
		h_Fill_CollMass_TauFake[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauFake_%s",systName[j].c_str());
		h_Fill_taupt_TauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauFake_%s",systName[j].c_str());
		h_Fill_taueta_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauFake_%s",systName[j].c_str());
		h_Fill_mupt_TauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauFake_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauFake[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauFake_%s",systName[j].c_str());
		h_Fill_mueta_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauFake_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauFake_%s",systName[j].c_str());	 
		h_Fill_metphi_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauFake_%s",systName[j].c_str());	
		h_Fill_NV_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_tauphi_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_tauphi_TauFake_%s",systName[j].c_str());

		h_Fill_tauphi_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_tauphi_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_muphi_SC_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_muphi_SC_TauFake_%s",systName[j].c_str());

		h_Fill_muphi_SC_TauFake[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_muphi_SC_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_tauphi_muphi_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_tauphi_muphi_TauFake_%s",systName[j].c_str());
		h_Fill_tauphi_muphi_TauFake[j] =  fs->make<TH2D>(histname,histname2,100,-3.4,3.4,100,-3.4,3.4);
		h_Fill_tauphi_muphi_TauFake[j]->Sumw2();

		sprintf(histname,"h_Fill_tauphi_metphi_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_tauphi_metphi_TauFake_%s",systName[j].c_str());

		h_Fill_tauphi_metphi_TauFake[j] =  fs->make<TH2D>(histname,histname2,100,-3.4,3.4,100,-3.4,3.4);
		h_Fill_tauphi_metphi_TauFake[j]->Sumw2();


		sprintf(histname,"h_Fill_muphi_metphi_TauFake_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_muphi_metphi_TauFake_%s",systName[j].c_str());

		h_Fill_muphi_metphi_TauFake[j] =  fs->make<TH2D>(histname,histname2,100,-3.4,3.4,100,-3.4,3.4);
		h_Fill_muphi_metphi_TauFake[j]->Sumw2();




		////


		sprintf(histname,"h_EleTauCharge_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauPass_%s",systName[j].c_str());
		h_EleTauCharge_TauPass[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauPass[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauPass_%s",systName[j].c_str());
		h_FillmT_Ele_TauPass[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauPass[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauPass_%s",systName[j].c_str());
		h_FillmT_Tau_TauPass[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauPass[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauPass_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauPass[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauPass_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauPass[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauPass[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauPass_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauPass[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauPass_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauPass[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauPass_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauPass[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauPass_%s",systName[j].c_str());
		h_Fill_TotalMass_TauPass[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauPass_%s",systName[j].c_str());
		h_Fill_Met_TauPass[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauPass_%s",systName[j].c_str());
		h_Fill_CollMass_TauPass[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauPass_%s",systName[j].c_str());
		h_Fill_taupt_TauPass[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauPass_%s",systName[j].c_str());
		h_Fill_taueta_TauPass[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauPass_%s",systName[j].c_str());
		h_Fill_mupt_TauPass[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauPass[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauPass_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauPass[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauPass_%s",systName[j].c_str());
		h_Fill_mueta_TauPass[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauPass_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauPass[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauPass_%s",systName[j].c_str());	 
		h_Fill_metphi_TauPass[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauPass[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauPass_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauPass_%s",systName[j].c_str());	
		h_Fill_NV_TauPass[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauPass[j]->Sumw2();


		/////


		sprintf(histname,"h_EleTauCharge_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauPassM_%s",systName[j].c_str());
		h_EleTauCharge_TauPassM[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauPassM[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauPassM_%s",systName[j].c_str());
		h_FillmT_Ele_TauPassM[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauPassM[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauPassM_%s",systName[j].c_str());
		h_FillmT_Tau_TauPassM[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauPassM[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauPassM_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauPassM[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauPassM_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauPassM[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauPassM[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauPassM_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauPassM[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauPassM_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauPassM[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauPassM_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauPassM[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauPassM_%s",systName[j].c_str());
		h_Fill_TotalMass_TauPassM[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauPassM_%s",systName[j].c_str());
		h_Fill_Met_TauPassM[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauPassM_%s",systName[j].c_str());
		h_Fill_CollMass_TauPassM[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauPassM_%s",systName[j].c_str());
		h_Fill_taupt_TauPassM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauPassM_%s",systName[j].c_str());
		h_Fill_taueta_TauPassM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauPassM_%s",systName[j].c_str());
		h_Fill_mupt_TauPassM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauPassM[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauPassM_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauPassM[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauPassM_%s",systName[j].c_str());
		h_Fill_mueta_TauPassM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauPassM_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauPassM[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauPassM_%s",systName[j].c_str());	 
		h_Fill_metphi_TauPassM[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauPassM[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauPassM_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauPassM_%s",systName[j].c_str());	
		h_Fill_NV_TauPassM[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauPassM[j]->Sumw2();



		///


		sprintf(histname,"h_EleTauCharge_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauPassM_SS_%s",systName[j].c_str());
		h_EleTauCharge_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauPassM_SS[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauPassM_SS_%s",systName[j].c_str());
		h_FillmT_Ele_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauPassM_SS_%s",systName[j].c_str());
		h_FillmT_Tau_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauPassM_SS[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauPassM_SS[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauPassM_SS[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_TotalMass_TauPassM_SS[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_Met_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_CollMass_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_taupt_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_taueta_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_mupt_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauPassM_SS[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauPassM_SS_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_mueta_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauPassM_SS_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauPassM_SS_%s",systName[j].c_str());	 
		h_Fill_metphi_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauPassM_SS[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauPassM_SS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauPassM_SS_%s",systName[j].c_str());	
		h_Fill_NV_TauPassM_SS[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauPassM_SS[j]->Sumw2();


		///

		sprintf(histname,"h_EleTauCharge_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_EleTauCharge_TauPassM_OS_%s",systName[j].c_str());
		h_EleTauCharge_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2,3,-1,2);
		h_EleTauCharge_TauPassM_OS[j]->Sumw2();


		sprintf(histname,"h_FillmT_Ele_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Ele_TauPassM_OS_%s",systName[j].c_str());
		h_FillmT_Ele_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Ele_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_FillmT_Tau_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_FillmT_Tau_TauPassM_OS_%s",systName[j].c_str());
		h_FillmT_Tau_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2, 1000,0,1000);
		h_FillmT_Tau_TauPassM_OS[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Met_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Met_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Met_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2,100, 0, 4);
		h_Fill_DPhi_Ele_Met_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_DPhi_Tau_Met_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Tau_Met_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_DPhi_Tau_Met_TauPassM_OS[j]  = fs->make<TH1D>(histname,histname2, 100, 0, 4);
		h_Fill_DPhi_Tau_Met_TauPassM_OS[j]->Sumw2();


		sprintf(histname,"h_Fill_DPhi_Ele_Tau_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_DPhi_Ele_Tau_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_DPhi_Ele_Tau_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2,  100, 0, 5);
		h_Fill_DPhi_Ele_Tau_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_EleTauMass_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_EleTauMass_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_EleTauMass_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_EleTauMass_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_TauMETMass_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TauMETMass_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_TauMETMass_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2, 2000, 0,2000);
		h_Fill_TauMETMass_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_TotalMass_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_TotalMass_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_TotalMass_TauPassM_OS[j] = fs->make<TH1D>(histname,histname2, 2000, 0 , 2000);
		h_Fill_TotalMass_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_Met_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_Met_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_Met_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_Met_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_CollMass_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_CollMass_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_CollMass_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2, 1000, 0 , 1000);
		h_Fill_CollMass_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_taupt_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taupt_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_taupt_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_taupt_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_taueta_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_taueta_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_taueta_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_taueta_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_mupt_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_mupt_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_TauPassM_OS[j]->Sumw2();


		sprintf(histname,"h_Fill_mupt_SC_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mupt_SC_TauPassM_OS_%s",systName[j].c_str());

		h_Fill_mupt_SC_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,500,0,500);
		h_Fill_mupt_SC_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_mueta_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_mueta_SC_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_mueta_SC_TauPassM_OS_%s",systName[j].c_str());
		h_Fill_mueta_SC_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,100,-2.5,2.5);
		h_Fill_mueta_SC_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_metphi_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_metphi_TauPassM_OS_%s",systName[j].c_str());	 
		h_Fill_metphi_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,100,-3.4,3.4);
		h_Fill_metphi_TauPassM_OS[j]->Sumw2();

		sprintf(histname,"h_Fill_NV_TauPassM_OS_%s",systName[j].c_str());
		sprintf(histname2,"h_Fill_NV_TauPassM_OS_%s",systName[j].c_str());	
		h_Fill_NV_TauPassM_OS[j] =  fs->make<TH1D>(histname,histname2,100,0,100);
		h_Fill_NV_TauPassM_OS[j]->Sumw2();

		//

		//

	}

	//////

}


// ------------ method called once each job just after ending the event loop  ------------
	void 
ETauAnalysis::endJob() 
{

	std::cout<<"Total No of nEventsRaw"<<":\t"<<nEventsRaw<<"\t: Stored events : \t"<<nEventsStored<<"\t weighted events :\t"<<mc_nEventsWeighted<<"\t sum of weights : \t"<<nEventsiihe<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ETauAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


///

bool ETauAnalysis::PassEleSelections(int eleindex, TLorentzVector ele, int CRele){

	//	if(!(tree->gsf_et->at(eleindex) > 50.)) return false;
	//	if(!(fabs(ele.Eta()) < 2.5)) return false;
	// HEEP ID

	if( CRele == 0 ) if(!tree->gsf_VIDHEEP7->at(eleindex)) return false;
	//////////////////	if( CRele == 1 ) if(!(PassPreSelections(eleindex, ele) && (tree->gsf_VIDHEEP7->at(eleindex)== false))) return false;

	return true;

}


bool ETauAnalysis::PassPreSelections(int iele, TLorentzVector ele) {
	bool pass_presel=false;
	bool preselectedEle=false;
	bool ismatchedTrig=false;

	if(fabs(tree->gsf_sc_eta->at(iele)) < 1.4442) {

		if(tree->gsf_sigmaIetaIeta->at(iele) < 0.013 && tree->gsf_hadronicOverEm->at(iele) < 0.15 && tree->gsf_nLostInnerHits->at(iele) <= 1 && fabs(tree->gsf_dxy_firstPVtx->at(iele)) < 0.02) { preselectedEle = true;}

	}       else if(fabs(tree->gsf_sc_eta->at(iele)) > 1.566 && fabs(tree->gsf_sc_eta->at(iele)) < 2.5) {
		if(tree->gsf_sigmaIetaIeta->at(iele) < 0.034 && tree->gsf_hadronicOverEm->at(iele) < 0.10 && tree->gsf_nLostInnerHits->at(iele) <= 1 && fabs(tree->gsf_dxy_firstPVtx->at(iele)) < 0.05) { preselectedEle = true;}

	}


	double DRtmp = 0.1;
	goodHLTElecIndex1 = 9999 ;
	vector<float> trig_double_ele_eta, trig_double_ele_phi;
	trig_double_ele_eta.clear();
	trig_double_ele_phi.clear();

	if(tree->ev_run >= 276453 && tree->ev_run <= 278822 ) {
		for(unsigned int trig1 = 0; trig1 < tree->trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta->size() ;  trig1++) {
			trig_double_ele_eta.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_eta->at(trig1));
			trig_double_ele_phi.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_hltEle33CaloIdLGsfTrkIdVLDPhiFilter_phi->at(trig1));
		}
	} else {


		for(unsigned int trig1 = 0; trig1 < tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLMWPMS2Filter_eta->size() ;  trig1++) {
			trig_double_ele_eta.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLMWPMS2Filter_eta->at(trig1));
			trig_double_ele_phi.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEG33CaloIdLMWPMS2Filter_phi->at(trig1));

		}
	}

	for(unsigned int trig1 = 0; trig1 < trig_double_ele_eta.size();  trig1++) {
		if(dR(tree->gsf_sc_eta->at(iele), tree->gsf_sc_phi->at(iele), trig_double_ele_eta.at(trig1), trig_double_ele_phi.at(trig1)  < DRtmp )) { DRtmp = dR(tree->gsf_sc_eta->at(iele), tree->gsf_sc_phi->at(iele),  trig_double_ele_eta.at(trig1), trig_double_ele_phi.at(trig1) ) ;

			ismatchedTrig = true; goodHLTElecIndex1=trig1 ; 
		}
	}


	if(preselectedEle && ismatchedTrig) pass_presel=true; 

	return pass_presel;
}

/*
   bool ETauAnalysis::PassHEEP(int eleindex) {
   bool isHeep(false);

   double ET = tree->gsf_caloEnergy->at(eleindex)*sin(2.*atan(exp(-1.*tree->gsf_eta->at(eleindex))));
   if ( (ET > 35.)  && fabs(tree->gsf_sc_eta->at(eleindex)) < 1.4442  &&
   fabs(tree->gsf_deltaEtaSeedClusterTrackAtVtx->at(eleindex)) < 0.004     &&
   fabs(tree->gsf_deltaPhiSuperClusterTrackAtVtx->at(eleindex)) < 0.06     &&
   tree->gsf_hadronicOverEm->at(eleindex) < (0.05 + 1./ tree->gsf_sc_energy->at(eleindex)) &&
   (tree->gsf_full5x5_e1x5->at(eleindex)/tree->gsf_full5x5_e5x5->at(eleindex) > 0.83 || tree->gsf_full5x5_e2x5Max->at(eleindex)/tree->gsf_full5x5_e5x5->at(eleindex) > 0.94) &&
   tree->gsf_nLostInnerHits->at(eleindex) < 2  &&
   fabs(tree->gsf_dxy_firstPVtx->at(eleindex)) < 0.02 &&
   tree->gsf_dr03EcalRecHitSumEt->at(eleindex) + tree->gsf_dr03HcalDepth1TowerSumEt->at(eleindex) < 2 + 0.03 * ET + 0.28 * tree->ev_fixedGridRhoFastjetAll   && tree->gsf_dr03TkSumPtHEEP7->at(eleindex) < 5) isHeep = true;

// endcap

if ( ET > 35  && (fabs(tree->gsf_sc_eta->at(eleindex)) > 1.566  && (fabs(tree->gsf_sc_eta->at(eleindex)) < 2.5) )&&
fabs(tree->gsf_deltaEtaSeedClusterTrackAtVtx->at(eleindex)) < 0.006                    &&
fabs(tree->gsf_deltaPhiSuperClusterTrackAtVtx->at(eleindex)) < 0.06                     &&
tree->gsf_hadronicOverEm->at(eleindex) < (0.05 + 5./ tree->gsf_sc_energy->at(eleindex))           &&
tree->gsf_full5x5_sigmaIetaIeta->at(eleindex) < 0.03                                         &&
tree->gsf_nLostInnerHits->at(eleindex) < 2                                               &&
fabs(tree->gsf_dxy_firstPVtx->at(eleindex)) < 0.05                       &&
(( ET < 50 && tree->gsf_dr03EcalRecHitSumEt->at(eleindex) + tree->gsf_dr03HcalDepth1TowerSumEt->at(eleindex) < 2.5 + 0.28 * tree->ev_fixedGridRhoFastjetAll) || 
( ET > 50 && tree->gsf_dr03EcalRecHitSumEt->at(eleindex) + tree->gsf_dr03HcalDepth1TowerSumEt->at(eleindex) < 2.5 + 0.03 * (ET-50) + 0.28 * tree->ev_fixedGridRhoFastjetAll)) &&
tree->gsf_dr03TkSumPtHEEP7->at(eleindex) < 5) isHeep = true;

return isHeep;


}
*/


bool ETauAnalysis::PassTauSelections( int tauindex, TLorentzVector tau){

	if(!(tau.Pt() > 30.)) return false;
	if(!(fabs(tau.Eta()) < 2.3)) return false;

	// discrimination against e/mu

	if(!(tree->tau_decayModeFinding->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_byVLooseIsolationMVArun2v1DBoldDMwLT->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_againstMuonLoose3->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_againstElectronTightMVA6->at(tauindex) > 0.5)) return false;

	return true;
}



bool ETauAnalysis::OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize) {
	if(dR(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi()) <= conesize) return true;
	else return false;
}  

float ETauAnalysis::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

} 

float ETauAnalysis::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}




float ETauAnalysis::mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
	float mt = -1;
	float pX = mupx+metx;
	float pY = mupy+mety;
	float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}      

/*
   float ETauAnalysis::PZetaVis( int muindex, int tauindex){
   float pzetavis;
   pzetavis = 999;
   TLorentzVector tau, mu;  
   tau.SetPtEtaPhiE(tree->tau_pt->at(tauindex),tree->tau_eta->at(tauindex),tree->tau_phi->at(tauindex),tree->tau_energy->at(tauindex));
   mu.SetPtEtaPhiM(tree->->at(muindex), tree->mu_gt_eta->at(muindex), tree->mu_gt_phi->at(muindex),mu_mass);
   float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
   float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
   float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2));
   zetax = zetax/zetaR;
   zetay = zetay/zetaR;

   float visPx = mu.Px() + tau.Px() ;
   float visPy = mu.Py() + tau.Py() ;

   pzetavis = visPx*zetax + visPy*zetay;
   return pzetavis;

   }

   float ETauAnalysis::PZeta(int muindex, int tauindex , float metpx, float metpy){
   float pzeta;
   pzeta = 999;
   TLorentzVector tau, mu;                                         
   tau.SetPxPyPzE(tree->tau_px->at(tauindex),tree->tau_py->at(tauindex),tree->tau_pz->at(tauindex),tree->tau_energy->at(tauindex));
   mu.SetPtEtaPhiM(tree->mu_gt_pt->at(muindex), tree->mu_gt_eta->at(muindex), tree->mu_gt_phi->at(muindex),mu_mass);

   float zetax = TMath::Cos(mu.Phi()) + TMath::Cos(tau.Phi()) ;
   float zetay = TMath::Sin(mu.Phi()) + TMath::Sin(tau.Phi()) ;
   float zetaR = TMath::Sqrt(pow(zetax,2) + pow(zetay,2)); 
   zetax = zetax/zetaR;                                    
   zetay = zetay/zetaR;         

   float vPx = mu.Px() + tau.Px()+metpx ;
   float vPy = mu.Py() + tau.Py()+metpy ;

   pzeta = vPx*zetax + vPy*zetay;
   return pzeta;


   }    
   */

void ETauAnalysis::HistoFiller(TH1F *histo, double value, double weight){                                    

	histo->Fill(value, weight);

}  



void ETauAnalysis::initializePileupInfo() {

	// Filenames must be c_strings below. Here is the conversion from strings to c_strings
	//         //   // As you can see above cstr1 corresponds to MC and cstr2 corresponds to data.
	//
	char * cstr1;
	char * cstr2;

	//Filenames must be c_strings below. Here is the conversion from strings
	cstr1 = new char [MCHistosForPU.size()+1];
	strcpy (cstr1, MCHistosForPU.c_str());
	cstr2 = new char [DataHistosForPU.size()+1];
	strcpy (cstr2, DataHistosForPU.c_str());

	TFile *file1 = TFile::Open(cstr1);
	TH1* histmc = dynamic_cast<TH1*>(file1->Get("pileup"));
	if(!histmc) {throw std::runtime_error("failed to extract histogram");}
	for(int bin=0; bin<=(histmc->GetXaxis()->GetNbins() + 1); bin++) {
		hPUmc->SetBinContent(bin,histmc->GetBinContent(bin));
	}      
	file1->Close();

	TFile *file2 = TFile::Open(cstr2);
	//file:/uscms_data/d3/aman30/RUN2_50ns/CMSSW_7_4_15/src/Plots/Code/Data50ns.root");
	TH1* histdata = dynamic_cast<TH1*>(file2->Get("pileup"));
	if(!histdata) {throw std::runtime_error("failed to extract histogram");}
	for(int bin=0; bin<=(histdata->GetXaxis()->GetNbins() + 1); bin++) {
		hPUdata->SetBinContent(bin,histdata->GetBinContent(bin));
	}
	file2->Close();
}        

double ETauAnalysis::getPileupWeight(float ntruePUInt) {
	int bin;
	double MCintegral;
	double MCvalue;
	double Dataintegral;
	double Datavalue;

	bin = hPUmc->GetBin(ntruePUInt+1);
	MCvalue = hPUmc->GetBinContent(bin);
	MCintegral = hPUmc->Integral();
	Datavalue = hPUdata->GetBinContent(bin);
	Dataintegral = hPUdata->Integral();
	if((MCvalue * Dataintegral) != 0) {pu_weight = (Datavalue * MCintegral) / (MCvalue * Dataintegral);}
	else {pu_weight = 1.0;}
	return pu_weight;
}       


bool ETauAnalysis::MatchingToGenMuons(TLorentzVector tau, int &genindex){
	double DRact = 0.5;
	genindex=-1;  
	bool ismatched=false;
	for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

		TLorentzVector gen_part;
		gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));
		bool isMuon = abs(tree->mc_pdgId->at(iGen))==13  ? true : false ;    //change names to Muon
		// muon from Z == not checked for now
		if(isMuon) {
			if(tau.DeltaR(gen_part) < DRact) {

				DRact=tau.DeltaR(gen_part);
				genindex=iGen;
				ismatched=true;

			}

		}     
	}             
	return ismatched;
}


int ETauAnalysis::GenTaus(){

	int HadronicTau = 0;
	for (unsigned int iGen = 0; iGen < gen_tau_had_pt.size(); iGen++){
		TLorentzVector gen_part;
		gen_part.SetPtEtaPhiE(gen_tau_had_pt.at(iGen), gen_tau_had_eta.at(iGen), gen_tau_had_phi.at(iGen), gen_tau_had_energy.at(iGen));
		if( gen_part.Pt() > 15.) { HadronicTau++; }


	}
	return HadronicTau;
}

bool ETauAnalysis::MatchingToGenTaus(TLorentzVector tau, int &genindex){
	double DRact = 0.2;                    
	genindex=-1;     
	bool ismatched=false;           
	for (unsigned int iGen = 0; iGen < gen_tau_had_pt.size(); iGen++){
		TLorentzVector gen_part3;       
		gen_part3.SetPtEtaPhiE(gen_tau_had_pt.at(iGen), gen_tau_had_eta.at(iGen), gen_tau_had_phi.at(iGen), gen_tau_had_energy.at(iGen));
		if( gen_part3.Pt() > 15.) {
			if(tau.DeltaR(gen_part3) < DRact) {                       

				DRact=tau.DeltaR(gen_part3);                      
				genindex=iGen;                                   
				ismatched=true;                                  

			}                                                        
		}
	}                                                                
	return ismatched;                                                        
} 

bool ETauAnalysis::FillChain(TChain *chain, const char* inputFileList) {

	ifstream infile(inputFileList);
#if 0
	std::string buffer;
#endif
	if (!infile.is_open()) {
		std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
		return kFALSE;
	}  

	std::cout << "TreeUtilities : FillChain " << std::endl;
	char buffer[255];
	while (infile) {
		infile.getline(buffer, 255);                // delim defaults to '\n'
		if (!infile.good()) break;
		std::cout << "Adding " << buffer << " to chain" << std::endl;
		chain->Add(buffer);
	}  
#if 0
	while (1) {
		infile >> buffer;
		if (!infile.good()) break;
		chain->Add(buffer.c_str());
	}  
#endif
	std::cout << "No. of Entries in this tree : " << chain->GetEntries() << std::endl;
	infile.close();
	return kTRUE;
} 

float ETauAnalysis::GetEfficiency(float eta, float pt, TH2F *hist) 
{

	double sf = 1.;
	int nbins = hist->GetXaxis()->GetNbins();

	if (eta > (hist->GetXaxis()->GetBinLowEdge(nbins) +  hist->GetXaxis()->GetBinWidth(nbins)))
		eta =  hist->GetXaxis()->GetBinLowEdge(nbins) + (hist->GetXaxis()->GetBinWidth(nbins)/2.0);

	nbins = hist->GetYaxis()->GetNbins();
	if (pt > (hist->GetYaxis()->GetBinLowEdge(nbins) +  hist->GetYaxis()->GetBinWidth(nbins)))
		pt =  hist->GetYaxis()->GetBinLowEdge(nbins) + (hist->GetYaxis()->GetBinWidth(nbins)/2.0);


	Int_t binX = hist->GetXaxis()->FindFixBin(eta);
	Int_t binY = hist->GetYaxis()->FindFixBin(pt);
	if (pt >= 0 && pt <=1000){
		sf = hist->GetBinContent(binX, binY);}

	return sf;
}

double ETauAnalysis::GetCollinearMass(const TLorentzVector &tau, const TLorentzVector &mu,  const TLorentzVector MET) {

	double METproj=fabs((MET.Px()*tau.Px()+MET.Py()*tau.Py())/tau.Pt());
	double xth=1;
	if((tau.Pt()+METproj)!=0) xth=tau.Pt()/(tau.Pt()+METproj);
	//now calculate the visibsle mass
	double mass_vis=(tau+mu).M();
	return mass_vis/sqrt(xth);

}

double ETauAnalysis::TauPtBinSF_CR5(double taupt, double taueta, double tauratio){ 

	double SF=0;
	double fake_rate=1.;
	double fake_rate_corr=1.;

	/*
	   if(fabs(taueta) < 1.46 ) {
	   if(taupt >= 30.  && taupt < 40.) { fake_rate= 0.213875;}
	   else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.226474;}
	   else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.227316;}
	   else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.229475;}
	   else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.234918;}
	   else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.23698;}
	   else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.240892;}
	   else if(taupt >= 120.  && taupt < 150.) { fake_rate= 0.236748;}
	   else if(taupt >= 150.  && taupt < 300.) { fake_rate= 0.267417;}
	   else if(taupt >= 300.  && taupt < 1000.) { fake_rate= 0.254386;}

	   if(tauratio >= 0 && tauratio < 0.3) { fake_rate_corr= 0.312808;}
	   else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.51794 ;}
	   else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.631975 ;}
	   else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.718813 ;}
	   else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.851849 ;}
	   else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.02398 ;}
	   else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=1.38971 ;}
	   else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=1.94469 ;}
	   else if(tauratio >= 1.0 && tauratio < 3.0) { fake_rate_corr=1.25749 ;}
	   else if(tauratio >= 3 ) { fake_rate_corr=1.25749 ;}

	   } else if (fabs(taueta) > 1.56 ) {

	   if(taupt >= 30.  && taupt < 40.) { fake_rate= 0.229114;} 
	   else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.242096;}
	   else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.246945;}
	   else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.228754;}
	   else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.236739;}
	   else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.215014;}
	   else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.231663;}
	   else if(taupt >= 120.  && taupt < 150.) { fake_rate= 0.223865;}
	   else if(taupt >= 150.  && taupt < 300.) { fake_rate= 0.2009;}
	   else if(taupt >= 300.  && taupt < 1000.) { fake_rate= 0.147059;}



	   if(tauratio >= 0 && tauratio < 0.3) { fake_rate_corr= 0.312311;}
	   else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.634086 ;} 
	   else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.716417 ;} 
	   else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.792882 ;} 
	   else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.973854 ;} 
	   else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.3048;} 
	   else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=1.85971;} 
	   else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=1.27056 ;} 
	   else if(tauratio >= 1.0 && tauratio < 3.0) { fake_rate_corr=1.27591 ;} 
	   else if(tauratio >= 3 ) { fake_rate_corr=1.27591 ;} 

	   }
	   */
	if(fabs(taueta) < 1.46 ) {
		if(taupt >= 30.  && taupt < 40.) { fake_rate= 0.26154 ;}
		else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.276199;}
		else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.269829;}
		else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.263042;}
		else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.267563;}
		else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.257829;}
		else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.250019;}
		else if(taupt >= 120.  && taupt < 150.) { fake_rate= 0.232584;}
		else if(taupt >= 150.  && taupt < 300.) { fake_rate=0.253978;}
		else if(taupt >= 300.  && taupt < 1000.) { fake_rate=0.236631;}

		if(tauratio >= 0 && tauratio < 0.3) { fake_rate_corr= 0.245353;}
		else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.427461 ;}
		else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.528457 ;}
		else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.609129 ;}
		else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.738999;}
		else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.909256 ;}
		else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=1.30165 ;}
		else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=2.02029 ;}
		else if(tauratio >= 1.0 && tauratio < 3.0) { fake_rate_corr=1.17963 ;}
		else if(tauratio >= 3 ) { fake_rate_corr=1.17963 ;}

	} else {// if (fabs(taueta) > 1.56 ) {
		if(taupt >= 30.  && taupt < 40.) { fake_rate=0.292741;}
		else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.311107;}
		else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.314325;}
		else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.284192;}
		else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.290438;}
		else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.248139;}
		else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.275715;}
		else if(taupt >= 120.  && taupt < 150.) { fake_rate= 0.272504;}
		else if(taupt >= 150.  && taupt < 300.) { fake_rate= 0.209448;}
		else if(taupt >= 300.  && taupt < 1000.) { fake_rate= 0.164226;}

		if(tauratio >= 0. && tauratio < 0.3) { fake_rate_corr=0.23544 ;}
		else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.52569 ;}
		else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.603695 ;}
		else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.679337 ;}
		else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.865115;}
		else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.24529;}
		else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=2.0318 ;}
		else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=1.18442 ;}
		else if(tauratio >= 1.0 && tauratio < 3.0) {fake_rate_corr=1.21508 ;}
		else if(tauratio >= 3 ) { fake_rate_corr=1.21508 ;}
	}

	//		std::cout<<"fake_rate:"<<fake_rate<<"fake_rate_corr:"<<fake_rate_corr<<"\t";
	SF = fake_rate*fake_rate_corr;
	//	SF = fake_rate;
	return SF;

	}

	double ETauAnalysis::ElePtBinSF_CR5(double elept, double eleSCeta) {
		double SF=1;

		if(fabs(eleSCeta) < 1.4442) { if (elept >= 35. && elept < 131.6) { SF= 0.140 - (0.0029*elept) + (2.56*0.00001* pow(elept,2)) - (8.48 * 0.00000001 * pow(elept,3));} else if (elept >= 131.6 && elept < 359.3) { SF=0.020 - (0.00013*elept) + (3.50*0.0000001* pow(elept,2)) - (2.90 * 0.0000000001 * pow(elept,3));} else { SF=0.00514 + (4.73*0.0000001*elept);}} 

		if(fabs(eleSCeta) > 1.566 && fabs(eleSCeta) < 2.0) { if (elept >= 35. && elept < 125) {SF= 0.1012 - (0.00094*elept) + (3.37*0.000001* pow(elept,2));} else if (elept >= 125.0 && elept < 226.3) {SF= 0.0488 - (11.37*0.00001*elept);} else { SF =0.0241 + (1.24*0.000001*elept);}} 

		if(fabs(eleSCeta) > 2.0 && fabs(eleSCeta) < 2.5) {  if (elept >= 35. && elept < 152) { SF= 0.0622 - (0.00012*elept);} else {SF= 0.0387;}}

		return SF;
	}

	/*
	   bool ETauAnalysis::CheckBit (int number, int bitpos)
	   {
	   bool res = number & (1 << bitpos);
	   return res;
	   }

*/
	double ETauAnalysis::MutoTauFR(TLorentzVector tau) {
		int Otype=0;
		double DR=0.2;
		double SF_fake =1.;
		double SF_fake_1 =1.;

		// muon faking tau
		for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

			TLorentzVector gen_part;
			gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));


			bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 ) && (tree->mc_status_flags->at(iGen) ==1 )) ? true : false ;

			bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) ==1)) ? true : false ;

			bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) ==6))   ? true : false ;


			bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) && ( tree->mc_status_flags->at(iGen) == 6))   ? true : false ;
			if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
				if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {Otype=1;}
					else if(isElectronPrompt) {Otype=2;} else if( isMuonTau) {Otype=3;} else if(isElectronTau) {Otype=4;} 
				}
			}
		} 

		if(Otype == 1 || Otype == 3) {

			if( fabs(tau.Eta()) > 0 && fabs(tau.Eta()) < 0.4) { SF_fake = 1.47; }
			if( fabs(tau.Eta()) > 0.4 && fabs(tau.Eta()) < 0.8) { SF_fake = 1.55; }
			if( fabs(tau.Eta()) > 0.8 && fabs(tau.Eta()) < 1.2) { SF_fake = 1.33; }
			if( fabs(tau.Eta()) > 1.2 && fabs(tau.Eta()) < 1.7) { SF_fake = 1.72; }
			if( fabs(tau.Eta()) > 1.7 && fabs(tau.Eta()) < 2.3) { SF_fake = 2.50; }


		}  else if ( Otype == 2 || Otype == 4) {

			if( fabs(tau.Eta())  < 1.460) { SF_fake = 1.21; }
			if( fabs(tau.Eta())  > 1.558) { SF_fake = 1.38; }


		} else { SF_fake = 1.;}

		return SF_fake;
	}
	bool ETauAnalysis::matchedToGenObjetcs(TLorentzVector tau, unsigned int &typey) {
		bool ismatched= false;
		double DR=0.2;
		typey = 0;

		for (unsigned int iGen = 0; iGen < tree->mc_px->size(); iGen++){

			TLorentzVector gen_part;
			gen_part.SetPxPyPzE(tree->mc_px->at(iGen),tree->mc_py->at(iGen),tree->mc_pz->at(iGen),tree->mc_energy->at(iGen));


			bool isMuonPrompt = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==13 ) && (tree->mc_status_flags->at(iGen) ==1 )) ? true : false ;

			bool isElectronPrompt = ( gen_part.Pt() >  8. && ( abs(tree->mc_pdgId->at(iGen))==11 ) && ( tree->mc_status_flags->at(iGen) ==1)) ? true : false ;

			bool isElectronTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))==11 ) && (tree->mc_status_flags->at(iGen) ==6))   ? true : false ;


			bool isMuonTau = ( gen_part.Pt() >  8. && (abs(tree->mc_pdgId->at(iGen))== 13   ) && ( tree->mc_status_flags->at(iGen) == 6))   ? true : false ;

			if(isMuonPrompt || isElectronTau || isElectronPrompt || isMuonTau ) {
				if( tau.DeltaR(gen_part) < DR ) { if(isMuonPrompt) {typey=1;} 
					else if(isElectronPrompt) {typey=2;} else if( isMuonTau) {typey=3;} else if(isElectronTau) {typey=4;} ismatched = true;  break;}


			}
		}
		return ismatched;

	}


	double ETauAnalysis::FakeRate_SSMtLow(double taupt, double jetpt, TH2F *hist) {
		double SF=0.2;
		double reweight = 0;
		int iBin = hist->FindBin(taupt, jetpt);
		SF = hist->GetBinContent(iBin);
		if (SF != 1) reweight = SF/(1-SF);
		return reweight;

	}




	bool ETauAnalysis::TauMatchedToJet(TLorentzVector tau_p4, unsigned int &matched_jet_indx) {

		bool matched_to_reco_jet=false;
		TLorentzVector jet_p4(0.,0.,0.,0.);
		for (unsigned int ijet = 0; ijet < tree->jet_pt->size(); ijet++){
			if(!(tree->jet_pt->at(ijet) > 0.)) continue;
			if(!(fabs(tree->jet_eta->at(ijet)) < 2.3)) continue;
			if(!(tree->jet_isJetIDLoose->at(ijet))) continue;
			TLorentzVector jet_p4_tmp;
			jet_p4_tmp.SetPtEtaPhiE(tree->jet_pt->at(ijet), tree->jet_eta->at(ijet), tree->jet_phi->at(ijet), tree->jet_energy->at(ijet));
			if(!(tau_p4.DeltaR(jet_p4_tmp) < 0.2)) continue;
			matched_to_reco_jet=true;
			jet_p4=jet_p4_tmp;
			matched_jet_indx = ijet;
			break;

		}
		return matched_to_reco_jet;
	}



	bool ETauAnalysis::LHEinfoPt(double &Pt_Obj) { 
		bool ismatched_obj=false;
		TLorentzVector l_p4, nu_p4, lnu_p4;
		l_p4.SetPxPyPzE(0, 0, 0, 0);
		nu_p4.SetPxPyPzE(0, 0, 0, 0);
		lnu_p4.SetPxPyPzE(0, 0, 0, 0);
		int l_pdgid = 0, nu_pdgid = 0;

		for (unsigned int iLHE = 0; iLHE < tree->LHE_Pt->size(); ++iLHE) {

			if (  (fabs(tree->LHE_pdgid->at(iLHE)) == 11 || fabs(tree->LHE_pdgid->at(iLHE)) == 13 || fabs(tree->LHE_pdgid->at(iLHE)) == 15))  {
				l_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));
				l_pdgid = tree->LHE_pdgid->at(iLHE);
			}       
			if (abs(tree->LHE_pdgid->at(iLHE)) == 12 || abs(tree->LHE_pdgid->at(iLHE)) == 14 || abs(tree->LHE_pdgid->at(iLHE)) == 16) {
				nu_p4.SetPtEtaPhiE(tree->LHE_Pt->at(iLHE),tree->LHE_Eta->at(iLHE),tree->LHE_Phi->at(iLHE),tree->LHE_E->at(iLHE));   
				nu_pdgid = tree->LHE_pdgid->at(iLHE);
			}       


		}
		if(fabs(l_pdgid)+1 == fabs(nu_pdgid)) {
			lnu_p4 = l_p4 + nu_p4;
			Pt_Obj = lnu_p4.Pt();
			ismatched_obj=true;
		}

		return ismatched_obj;
	}


	double ETauAnalysis::FakeRate(double taupt, double jetpt) {
		double SF=0.2;

		TFile* fake_file = new TFile("Reweighting/fakerate.root","R");

		TString taupt_string = "_tau_pt_";

		if( taupt >= 30. && taupt < 50.)        taupt_string += "30_50";
		else if( taupt >= 50. && taupt < 100.)  taupt_string += "50_100";
		else if( taupt >= 100. )                taupt_string += "100";

		double reweight = 1;

		if (jetpt >= 1000) jetpt = 999;

		TString hname = "hpass_jetpt_pass_data_total_tauindiff" + taupt_string;
		TH1F* h_fake = (TH1F*) fake_file->Get(hname);
		int iBin = h_fake->FindBin(jetpt);
		SF = h_fake->GetBinContent(iBin);
		if (SF != 1) reweight = SF/(1-SF);

		return reweight;

	}
	DEFINE_FWK_MODULE(ETauAnalysis);
