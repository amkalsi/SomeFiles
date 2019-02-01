// change done 7 July 2018 -- add >= 1 bjets requirement
#include "ETauAnalysis_2017.h"
using namespace std;

ETauAnalysis_2017::ETauAnalysis_2017(const edm::ParameterSet& iConfig)

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



	isOS = iConfig.getParameter<bool>("isOS");
	isSS = iConfig.getParameter<bool>("isSS");

	//	std::cout<<"reached end:"<<std::endl;

}


ETauAnalysis_2017::~ETauAnalysis_2017()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
ETauAnalysis_2017::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	using namespace edm;
	nEventsRaw = nEventsStored = mc_nEventsWeighted = nEventsiihe = 0;


	////	std::cout<<"datafile::"<<datafile<<std::endl;
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

		if(!isdata) { weighthis = tree->mc_w_sign;}  else { weighthis = 1.0;}

		if(!isdata) {h_Events_After_GenFilter->Fill(1.);}
		//			std::cout<<"=========================="<<std::endl;
		if(!isdata) h_Fill_Mass_Gen_toChk->Fill(Mass,tree->mc_w_sign);

		if( !isdata) {  h_Count_Taus->Fill(GenTaus(), tree->mc_w_sign); }

		if(tree->trig_HLT_Ele35_WPTight_Gsf_accept) {h_fill_ele35->Fill(1);}
		if(tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept) {h_fill_ele115->Fill(1);}
		if(tree->trig_HLT_Photon200_accept) {h_fill_photon->Fill(1);}

		if(DYSample.find("singlephoton") != npos && isdata ) {

			if(!(!tree->trig_HLT_Ele35_WPTight_Gsf_accept && (!tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept) &&  (tree->trig_HLT_Photon200_accept))) continue;
		}

		if(DYSample.find("SE") != npos && isdata ){
			if( !( tree->trig_HLT_Ele35_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept )) continue;
		} 
		if(!isdata){
			if( !( tree->trig_HLT_Ele35_WPTight_Gsf_accept || tree->trig_HLT_Ele115_CaloIdVT_GsfTrkIdT_accept || tree->trig_HLT_Photon200_accept)) continue;
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

		//if(!isdata) {pu_weight =PU_reReco_Morind17::MC_pileup_weight(tree->mc_trueNumInteractions, 0, "all"); } else {pu_weight = 1;}


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
				if(!(tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstMuonLoose3->at(dau2index) > 0.5)) continue;
				if(!(tree->tau_againstElectronTightMVA6->at(dau2index) > 0.5)) continue;

				TLorentzVector DauTau ;
				DauTau.SetPxPyPzE(tree->tau_px->at(dau2index), tree->tau_py->at(dau2index), tree->tau_pz->at(dau2index),tree->tau_energy->at(dau2index));

				if(!(DauEle.DeltaR(DauTau) > 0.5 ))  continue;
				TVector2 METcorr;
				METcorr.SetMagPhi(tree->MET_eefix_Pt,tree->MET_eefix_phi); 
				double mt_cut = mTCalculation(METcorr.Px(), METcorr.Py(), DauEle.Px(), DauEle.Py(), DauEle.Pt());



				if( mt_cut < 120. ) continue;

				bool applySF;
				applySF= false;
				int genindex=-1;
				unsigned int matchgen_obj_type;
				if(!isdata) { if( (MatchingToGenTaus(DauTau, genindex) && genindex!=-1) || matchedToGenObjetcs(DauTau, matchgen_obj_type)) {applySF = true;} }
				if(isdata) {applySF= true;}
				if(!applySF) continue;

				if(tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5) {

					if( (DauEle + DauTau).M() > masspair) {
						masspair = (DauEle +DauTau).M();
						FirstObj = DauEle;
						SecondObj = DauTau;
						first_index = dau1index;
						second_index =dau2index;
					}
				}

				if(!(tree->tau_byTightIsolationMVArun2017v2DBnewDMwLT2017->at(dau2index) > 0.5)) {
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
				} // fail
			}
		} // loop ends

		double weight = weighthis;
		double weight_pass = weighthis;
		TVector2 metV2;


		metV2.SetMagPhi(tree->MET_eefix_Pt,tree->MET_eefix_phi);;
		TLorentzVector metV;
		metV.SetPxPyPzE(metV2.Px(), metV2.Py(), 0. , tree->MET_eefix_Pt);

		if(int(first_index_FF) != -1 && int(second_index_FF) != -1 ) {
//			if( (!isdata) )  { weight = weight * GetEfficiency(fabs(FirstObj_FF.Eta()), FirstObj_FF.Pt(), fhDMuMediumSF); }



			if(!isdata) { double fake_wt =  MutoTauFR(SecondObj_FF);
				weight = weight * fake_wt;
			}

			double fr_weight_ss_DY = TauPtBinSF_CR5( SecondObj_FF.Pt(), SecondObj_FF.Eta(), SecondObj_FF.Pt()/ tree->jet_pt->at(jet_index_FF));
			double unc_weight;
			unc_weight=1.;
			double tauSF= 1;

			for (unsigned int t=0; t <5; t++) {
				if(t==0) { unc_weight =tauSF;}


				if( t == 1) { unc_weight = (1+ 0.05)*tauSF; }
				if( t == 2) { unc_weight = (1- 0.05)*tauSF; }
				if( t == 3) { unc_weight = (1+ ((0.05*SecondObj.Pt())/1000.))*tauSF; }
				if( t == 4) { unc_weight = (1- ((0.35*SecondObj.Pt())/1000.))*tauSF; }

				if(isdata) unc_weight = 1.;
				h_Fill_NV_TauFake[t]->Fill(tree->pv_n,weight*unc_weight*fr_weight_ss_DY);
				h_EleTauCharge_TauFake[t]->Fill(tree->gsf_charge->at(first_index_FF) * tree->tau_charge->at(second_index_FF),weight*unc_weight*fr_weight_ss_DY);
				h_FillmT_Ele_TauFake[t]->Fill(mTCalculation(metV.Px(), metV.Py(), FirstObj_FF.Px(), FirstObj_FF.Py(), FirstObj_FF.Pt()),weight*unc_weight*fr_weight_ss_DY);
				h_FillmT_Tau_TauFake[t]->Fill(mTCalculation(metV.Px(), metV.Py(), SecondObj_FF.Px(), SecondObj_FF.Py(), SecondObj_FF.Pt()),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_DPhi_Ele_Met_TauFake[t]->Fill(deltaPhi(FirstObj_FF.Phi(), metV.Phi()),weight*unc_weight*fr_weight_ss_DY);
				h_Fill_DPhi_Tau_Met_TauFake[t]->Fill(deltaPhi(SecondObj_FF.Phi(), metV.Phi()),weight*unc_weight*fr_weight_ss_DY);
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
//			if( (!isdata) )  { weight_pass = weight_pass * GetEfficiency(fabs(FirstObj.Eta()), FirstObj.Pt(), fhDMuMediumSF); }

			if(!isdata)  {double fake_wt =  MutoTauFR(SecondObj);
				weight_pass = weight_pass * fake_wt;
			}


			double unc_weight;
			unc_weight=1.;
			double tauSF= 1;
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
	delete tree;





	file_in->Close();



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
ETauAnalysis_2017::beginJob()
{
	//	fIn = TFile::Open(InputFile.c_str());
	std::cout<<"InputFile::"<<InputFile<< std::endl;


	h_Count_Taus = fs->make<TH1D>("h_Count_Taus","h_Count_Taus",10,0,10); 
	h_Count_Taus->Sumw2();
	h_Events_Before_Skim = fs->make<TH1D>("h_Events_Before_Skim","h_Events_Before_Skim",2,0,2);
	h_Events_After_Skim = fs->make<TH1D>("h_Events_After_Skim","h_Events_After_Skim",2,0,2);
	h_Events_After_GenFilter = fs->make<TH1D>("h_Events_After_GenFilter","h_Events_After_GenFilter",2,0,2);
	h_Fill_Mass_Gen_toChk = fs->make<TH1D>("h_Fill_Mass_Gen_toChk","h_Fill_Mass_Gen_toChk",1000,0,1000);
	h_Fill_Mass_Gen_toChk->Sumw2();
	h_Events_Before_Skim->Sumw2();
	h_Events_After_Skim->Sumw2();
	h_Events_After_GenFilter->Sumw2();
	h_fill_ele35 = fs->make<TH1D>("h_fill_ele35","h_fill_ele35",2,0,2);
	h_fill_ele115 = fs->make<TH1D>("h_fill_ele115","h_fill_ele115",2,0,2);
	h_fill_photon = fs->make<TH1D>("h_fill_photon","h_fill_photon",2,0,2);

	h_fill_ele35->Sumw2();
	h_fill_ele115->Sumw2();
	h_fill_photon->Sumw2();
	string h_names[] = {"taupt_jetpt_pass_OS","taupt_jetpt_fail_OS","taupt_jetpt_pass_SS","taupt_jetpt_fail_SS"};
	unsigned int numH_name = sizeof(h_names)/sizeof(string);

	string dms[] = {"DM0", "DM1", "DM10"};
	unsigned int numH_dms = sizeof(dms)/sizeof(string);

	string etaarr[] = {"barrel","endcap"};
	unsigned int num_etaarr = sizeof(etaarr)/sizeof(string);

	for (unsigned int i = 0; i< numH_name; ++i) {
		for (unsigned int k = 0; k< numH_dms ; ++k) {
			for (unsigned int l = 0; l< num_etaarr; ++l) {
				TString nname = h_names[i]+"_"+dms[k]+"_"+etaarr[l];

				TString nname1 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"tau_pt";
				TString nname2 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"jet_pt";
				TString nname3 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio";
				TString nname4 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_1";
				TString nname5 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_2";
				TString nname6 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_3";
				TString nname7 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_4";
				TString nname8 = h_names[i]+"_"+dms[k]+"_"+etaarr[l]+"_ratio_5";


				hh[i][k][l] =  fs->make<TH2F>(nname, nname, 1000, 0, 1000, 1000, 0, 1000) ;
				h1[i][k][l] = fs->make<TH1F>(nname1, nname1, 1000, 0, 1000) ;
				h2[i][k][l] = fs->make<TH1F>(nname2, nname2, 1000, 0, 1000) ;
				h3[i][k][l] = fs->make<TH1F>(nname3, nname3, 1000, 0, 10) ;
				h4[i][k][l] = fs->make<TH1F>(nname4, nname4, 1000, 0, 10) ;
				h5[i][k][l] = fs->make<TH1F>(nname5, nname5, 1000, 0, 10) ;
				h6[i][k][l] = fs->make<TH1F>(nname6, nname6, 1000, 0, 10) ;
				h7[i][k][l] = fs->make<TH1F>(nname7, nname7, 1000, 0, 10) ;
				h8[i][k][l] = fs->make<TH1F>(nname8, nname8, 1000, 0, 10) ;

				hh[i][k][l]->Sumw2();
				h1[i][k][l]->Sumw2();
				h2[i][k][l]->Sumw2();
				h3[i][k][l]->Sumw2();
				h4[i][k][l]->Sumw2();
				h5[i][k][l]->Sumw2();
				h6[i][k][l]->Sumw2();
				h7[i][k][l]->Sumw2();
				h8[i][k][l]->Sumw2();

			}
		}
	}


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
ETauAnalysis_2017::endJob() 
{

	std::cout<<"Total No of nEventsRaw"<<":\t"<<nEventsRaw<<"\t: Stored events : \t"<<nEventsStored<<"\t weighted events :\t"<<mc_nEventsWeighted<<"\t sum of weights : \t"<<nEventsiihe<<std::endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ETauAnalysis_2017::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}


///

bool ETauAnalysis_2017::PassEleSelections(int eleindex, TLorentzVector ele, int CRele){

	//	if(!(tree->gsf_et->at(eleindex) > 50.)) return false;
	//	if(!(fabs(ele.Eta()) < 2.5)) return false;
	// HEEP ID

	if( CRele == 0 ) if(!tree->gsf_VIDHEEP7->at(eleindex)) return false;
	//////////////////	if( CRele == 1 ) if(!(PassPreSelections(eleindex, ele) && (tree->gsf_VIDHEEP7->at(eleindex)== false))) return false;

	return true;

}


bool ETauAnalysis_2017::PassPreSelections(int iele, TLorentzVector ele) {
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
		for(unsigned int trig1 = 0; trig1 < tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta->size() ;  trig1++) {
			trig_double_ele_eta.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta->at(trig1));
			trig_double_ele_phi.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi->at(trig1));
		}
	} else {


		for(unsigned int trig1 = 0; trig1 < tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta->size() ;  trig1++) {
			trig_double_ele_eta.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_eta->at(trig1));
			trig_double_ele_phi.push_back(tree->trig_HLT_DoubleEle33_CaloIdL_MW_hltEle33CaloIdLMWPMS2Filter_phi->at(trig1));

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
   bool ETauAnalysis_2017::PassHEEP(int eleindex) {
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


bool ETauAnalysis_2017::PassTauSelections( int tauindex, TLorentzVector tau){

	if(!(tau.Pt() > 30.)) return false;
	if(!(fabs(tau.Eta()) < 2.3)) return false;

	// discrimination against e/mu

	if(!(tree->tau_decayModeFinding->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_byVLooseIsolationMVArun2017v2DBnewDMwLT2017->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_againstMuonLoose3->at(tauindex) > 0.5)) return false;
	if(!(tree->tau_againstElectronTightMVA6->at(tauindex) > 0.5)) return false;

	return true;
}



bool ETauAnalysis_2017::OverLap05(TLorentzVector l1 , TLorentzVector l2, float conesize) {
	if(dR(l1.Eta(), l1.Phi(), l2.Eta(), l2.Phi()) <= conesize) return true;
	else return false;
}  

float ETauAnalysis_2017::deltaPhi( float a, float b) {
	float result = a-b;
	while (result > M_PI) result -= 2* M_PI;
	while (result <= -M_PI) result += 2* M_PI;
	return (fabs(result));

} 

float ETauAnalysis_2017::dR(float l1eta, float l1phi, float l2eta, float l2phi ) {
	float deta = l1eta - l2eta;
	float dphi = deltaPhi(l1phi,l2phi);
	return sqrt(deta*deta + dphi*dphi);
}




float ETauAnalysis_2017::mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
	float mt = -1;
	float pX = mupx+metx;
	float pY = mupy+mety;
	float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
	mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
	return mt;

}      

/*
   float ETauAnalysis_2017::PZetaVis( int muindex, int tauindex){
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

   float ETauAnalysis_2017::PZeta(int muindex, int tauindex , float metpx, float metpy){
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

void ETauAnalysis_2017::HistoFiller(TH1F *histo, double value, double weight){                                    

	histo->Fill(value, weight);

}  



void ETauAnalysis_2017::initializePileupInfo() {

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

double ETauAnalysis_2017::getPileupWeight(float ntruePUInt) {
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


bool ETauAnalysis_2017::MatchingToGenMuons(TLorentzVector tau, int &genindex){
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


int ETauAnalysis_2017::GenTaus(){

	int HadronicTau = 0;
	for (unsigned int iGen = 0; iGen < tree->mc_tau_had_pt->size(); iGen++){
		TLorentzVector gen_part;
		gen_part.SetPtEtaPhiE(tree->mc_tau_had_pt->at(iGen), tree->mc_tau_had_eta->at(iGen), tree->mc_tau_had_phi->at(iGen), tree->mc_tau_had_energy->at(iGen));
		if( gen_part.Pt() > 15.) { HadronicTau++; }


	}
	return HadronicTau;
}

bool ETauAnalysis_2017::MatchingToGenTaus(TLorentzVector tau, int &genindex){
	double DRact = 0.2;                    
	genindex=-1;     
	bool ismatched=false;           
	for (unsigned int iGen = 0; iGen < tree->mc_tau_had_pt->size(); iGen++){
		TLorentzVector gen_part3;       
		gen_part3.SetPtEtaPhiE(tree->mc_tau_had_pt->at(iGen), tree->mc_tau_had_eta->at(iGen), tree->mc_tau_had_phi->at(iGen), tree->mc_tau_had_energy->at(iGen));
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

bool ETauAnalysis_2017::FillChain(TChain *chain, const char* inputFileList) {

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

float ETauAnalysis_2017::GetEfficiency(float eta, float pt, TH2F *hist) 
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

double ETauAnalysis_2017::GetCollinearMass(const TLorentzVector &tau, const TLorentzVector &mu,  const TLorentzVector MET) {

	double METproj=fabs((MET.Px()*tau.Px()+MET.Py()*tau.Py())/tau.Pt());
	double xth=1;
	if((tau.Pt()+METproj)!=0) xth=tau.Pt()/(tau.Pt()+METproj);
	//now calculate the visibsle mass
	double mass_vis=(tau+mu).M();
	return mass_vis/sqrt(xth);

}

double ETauAnalysis_2017::TauPtBinSF_CR5(double taupt, double taueta, double tauratio){ 
 double SF=0;
        double fake_rate=1.;
        double fake_rate_corr=1.;

 if(fabs(taueta) < 1.46 ) {
                if(taupt >= 30.  && taupt < 40.) { fake_rate= 0.216369;}
                else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.207415;}
                else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.197045;}
                else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.190562;}
                else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.184203;}
                else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.180642;}
                else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.195444;}
                else if(taupt >= 120.  && taupt < 150.) { fake_rate=0.180619;}
                else if(taupt >= 150.  && taupt < 300.) { fake_rate=0.199716;}
                else if(taupt >= 300.  && taupt < 1000.) { fake_rate=0.471206;}

                if(tauratio >= 0 && tauratio < 0.3) { fake_rate_corr= 0.181674;}
                else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.397122 ;}
                else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.480164 ;}
                else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.601786 ;}
                else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.773796;}
                else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.10783 ;}
                else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=1.93134 ;}
                else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=1.73882 ;}
                else if(tauratio >= 1.0 && tauratio < 3.0) { fake_rate_corr=1.16325 ;}
                else if(tauratio >= 3 ) { fake_rate_corr=1.16325 ;}
} else {
// check biining could be different
  if(taupt >= 30.  && taupt < 40.) { fake_rate=0.230194;}
                else if(taupt >= 40.  && taupt < 50.) { fake_rate= 0.226093;}
                else if(taupt >= 50.  && taupt < 60.) { fake_rate= 0.217827;}
                else if(taupt >= 60.  && taupt < 70.) { fake_rate= 0.216928;}
                else if(taupt >= 70.  && taupt < 80.) { fake_rate= 0.226077;}
                else if(taupt >= 80.  && taupt < 100.) { fake_rate= 0.197912;}
                else if(taupt >= 100.  && taupt < 120.) { fake_rate= 0.208331;}
                else if(taupt >= 120.  && taupt < 150.) { fake_rate= 0.237798;}
              
  else if(taupt >= 150.  && taupt < 200.) { fake_rate= 0.194021;}
// change here
                else if(taupt >= 200.  && taupt < 1000.) { fake_rate= 0.196989;}

                if(tauratio >= 0. && tauratio < 0.3) { fake_rate_corr=0.16431 ;}
                else if(tauratio >= 0.3 && tauratio < 0.4) { fake_rate_corr=0.509698 ;}
                else if(tauratio >= 0.4 && tauratio < 0.5) { fake_rate_corr=0.60256 ;}
                else if(tauratio >= 0.5 && tauratio < 0.6) { fake_rate_corr=0.727782 ;}
                else if(tauratio >= 0.6 && tauratio < 0.7) { fake_rate_corr=0.966557;}
                else if(tauratio >= 0.7 && tauratio < 0.8) { fake_rate_corr=1.5371;}
                else if(tauratio >= 0.8 && tauratio < 0.9) { fake_rate_corr=1.61498 ;}
                else if(tauratio >= 0.9 && tauratio < 1.0) { fake_rate_corr=1.03456 ;}
                else if(tauratio >= 1.0 && tauratio < 3.0) {fake_rate_corr=1.07144 ;}
                else if(tauratio >= 3 ) { fake_rate_corr=1.07144 ;}


}

        SF = fake_rate*fake_rate_corr;

return SF;

}

double ETauAnalysis_2017::ElePtBinSF_CR5(double elept, double eleSCeta) {
	double SF=1;

	if(fabs(eleSCeta) < 1.4442) { if (elept >= 35. && elept < 131.6) { SF= 0.140 - (0.0029*elept) + (2.56*0.00001* pow(elept,2)) - (8.48 * 0.00000001 * pow(elept,3));} else if (elept >= 131.6 && elept < 359.3) { SF=0.020 - (0.00013*elept) + (3.50*0.0000001* pow(elept,2)) - (2.90 * 0.0000000001 * pow(elept,3));} else { SF=0.00514 + (4.73*0.0000001*elept);}} 

	if(fabs(eleSCeta) > 1.566 && fabs(eleSCeta) < 2.0) { if (elept >= 35. && elept < 125) {SF= 0.1012 - (0.00094*elept) + (3.37*0.000001* pow(elept,2));} else if (elept >= 125.0 && elept < 226.3) {SF= 0.0488 - (11.37*0.00001*elept);} else { SF =0.0241 + (1.24*0.000001*elept);}} 

	if(fabs(eleSCeta) > 2.0 && fabs(eleSCeta) < 2.5) {  if (elept >= 35. && elept < 152) { SF= 0.0622 - (0.00012*elept);} else {SF= 0.0387;}}

	return SF;
}

/*
   bool ETauAnalysis_2017::CheckBit (int number, int bitpos)
   {
   bool res = number & (1 << bitpos);
   return res;
   }

*/
double ETauAnalysis_2017::MutoTauFR(TLorentzVector tau) {
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

		if( fabs(tau.Eta()) > 0 && fabs(tau.Eta()) < 0.4) { SF_fake = 1.06; }
		if( fabs(tau.Eta()) > 0.4 && fabs(tau.Eta()) < 0.8) { SF_fake = 1.02; }
		if( fabs(tau.Eta()) > 0.8 && fabs(tau.Eta()) < 1.2) { SF_fake = 1.10; }
		if( fabs(tau.Eta()) > 1.2 && fabs(tau.Eta()) < 1.7) { SF_fake = 1.03; }
		if( fabs(tau.Eta()) > 1.7 && fabs(tau.Eta()) < 2.3) { SF_fake = 1.94; }


	}  else if ( Otype == 2 || Otype == 4) {

		if( fabs(tau.Eta())  < 1.460) { SF_fake = 1.96; }
		if( fabs(tau.Eta())  > 1.558) { SF_fake = 1.66; }


	} else { SF_fake = 1.;}

	return SF_fake;
}
bool ETauAnalysis_2017::matchedToGenObjetcs(TLorentzVector tau, unsigned int &typey) {
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


double ETauAnalysis_2017::FakeRate_SSMtLow(double taupt, double jetpt, TH2F *hist) {
	double SF=0.2;
	double reweight = 0;
	int iBin = hist->FindBin(taupt, jetpt);
	SF = hist->GetBinContent(iBin);
	if (SF != 1) reweight = SF/(1-SF);
	return reweight;

}




bool ETauAnalysis_2017::TauMatchedToJet(TLorentzVector tau_p4, unsigned int &matched_jet_indx) {

	bool matched_to_reco_jet=false;
	TLorentzVector jet_p4(0.,0.,0.,0.);
	for (unsigned int ijet = 0; ijet < tree->jet_pt->size(); ijet++){
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


double ETauAnalysis_2017::FakeRate(double taupt, double jetpt) {
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
DEFINE_FWK_MODULE(ETauAnalysis_2017);
