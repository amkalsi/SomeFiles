// An analyzer to produce the b-tagged jet / b-tagged + mu PT ratio.
// USAGE:
// gSystem->Load("<path_to_delphes>/Delphes/libDelphes.so");
// .x DelphesAnalyzer.C("input_file.root", "saved_histogram_name.root", "plot_name", nbins);
// Set env using folllowing lines before
/*
   export LD_LIBRARY_PATH=/Users/amkalsi/Downloads/MG5_aMC_v2_9_2/Delphes:$LD_LIBRARY_PATH
   export ROOT_INCLUDE_PATH=/Users/amkalsi/Downloads/MG5_aMC_v2_9_2//Delphes/external/
   export ROOT_INCLUDE_PATH=/Users/amkalsi/Downloads/MG5_aMC_v2_9_2/Delphes/classes/:$ROOT_INCLUDE_PATH
   export PYTHONPATH=$PYTHONPATH:$ROOT_INCLUDE_PATH
 */

// * run using following command //

//root -l DelphesAnalyzer.C'("tag_3_delphes_events.root", "saved_histogram_name.root", "plot_name", 1000)'

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "math.h"
#include "TLorentzVector.h"
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
    class ExRootTreeReader;
    class ExRootResult;
#endif
float mTCalculation(float metx, float mety, float mupx, float mupy, float mupt);

float mTCalculation(float metx, float mety, float mupx, float mupy, float mupt){
        float mt = -1;
        float pX = mupx+metx;
        float pY = mupy+mety;
        float et = mupt + TMath::Sqrt(metx*metx + mety*mety);
        mt = TMath::Sqrt(et*et-(pX*pX + pY*pY));
        return mt;

}


    void DelphesAnalyzer(const char *file_name, const char *save_name, const char *plot_name, int nbins) {
        gSystem->Load("libDelphes");

        // Create a chain of root trees.
        TChain chain("Delphes");
        chain.Add(file_name);

        // Create object of class ExRootTreeReader.
        ExRootTreeReader *tree_reader = new ExRootTreeReader(&chain);
        Long64_t number_of_entries = tree_reader->GetEntries();

        // Get pointers to branches used in this analysis.
        //  TClonesArray *branch_mu = tree_reader->UseBranch("Muon");
        TClonesArray *branch_jet = tree_reader->UseBranch("Jet");
        TClonesArray *branch_met = tree_reader->UseBranch("MissingET");

        // Book histograms.
        TH1F  *h_jet_pt = new TH1F("h_jet_pt","h_jet_pt",100,0,1000);
        TH1F  *h_bjet_pt = new TH1F("h_bjet_pt","h_bjet_pt",100,0,1000);
        TH1F  *h_jet_eta = new TH1F("h_jet_eta","h_jet_eta",100,-5,5);
        TH1F  *h_bjet_eta = new TH1F("h_bjet_eta","h_bjet_eta",100,-5,5);

        TH1F  *h_jet_phi = new TH1F("h_jet_phi","h_jet_phi",100,-5,5);
        TH1F  *h_bjet_phi = new TH1F("h_bjet_phi","h_bjet_phi",100,-5,5);


        TH1F  *h_MET = new TH1F("h_MET","h_MET",150,0,1500);
        TH1F  *h_MT = new TH1F("h_MT","h_MT",200,0,2000);

/*

https://github.com/sethzenz/Delphes/blob/master/examples/GeneralExample.C

  // Status code 3 (+high pt leptons, b+t quarks, etc)  particle collection
    if (verbose && branchGenParticle) {
      for (int i = 0 ; i < branchGenParticle->GetEntries() ; i++ ) {
    GenParticle *part = (GenParticle*) branchGenParticle->At(i);
    cout << "     Status code" << part->Status << " generator particle PID Pt Eta Phi Z T (at origin) "  << part->PID << " "
         << part->PT << " " << part->Eta << " " << part->Phi << " " << part->Z << " " << part->T << endl;
      }
    }
    */


        //TH1F *h_pt_ratio = new TH1F(plot_name, "", nbins, 0., 1.);
        //h_pt_ratio->GetXaxis()->SetTitle("Pt(b)/(Pt(b) + Pt(mu))");
        //h_pt_ratio->GetYaxis()->SetTitle("Nb events");

        printf("next event with size %lld \n", number_of_entries);

        // Event loop.
        for(Int_t entry = 0; entry < number_of_entries; ++entry) {

            tree_reader->ReadEntry(entry);

            // Declare physics objects.
            // Muon *muon;
            Jet *jet;
            MissingET *met;

            // Declare running indices.
            bool btag_found = false;
            bool non_b_jet_found = false;
            //Double_t muon_pt = 0.0;
            //Double_t muon_phi, muon_eta;
            Double_t btag_pt, btag_phi, btag_eta, nonb_jet_pt, nonb_jet_eta, nonb_jet_phi;
            double btag_px, btag_py, nonb_jet_px, nonb_jet_py;
            //Int_t muon_n = branch_mu->GetEntries();
            Int_t jet_n = branch_jet->GetEntries();
            Int_t met_n = branch_met->GetEntries();

            if (jet_n > 0) {

                // Loop over jets and find the leading b-jet.
                for (Int_t i = 0; i < jet_n; i++) {
                    jet = (Jet*) branch_jet->At(i);
                    if (jet->BTag == 1 && !btag_found) {
                        // add b-jet pt and eta cut
                        if( jet->PT > 20 && fabs(jet->Eta) < 2.5 ) {
                            btag_found = true;
                            btag_pt = jet->PT;
                            btag_phi = jet->Phi;
                            btag_eta = jet->Eta;
                            btag_px = jet->PT*cos(jet->Phi);
                            btag_py = jet->PT*sin(jet->Phi);

                        }
                    }

                    // look for non-bjets with pt > 100

                    if(jet->BTag == 0 && jet->PT > 100 && !non_b_jet_found) {
                        non_b_jet_found= true;
                        nonb_jet_pt =  jet->PT;
                        nonb_jet_eta = jet->Eta;
                        nonb_jet_phi = jet->Phi;
                        nonb_jet_px = jet->PT*cos(jet->Phi);
                        nonb_jet_py = jet->PT*sin(jet->Phi);

                    }
                }

                // Loop over muons and find the highest pT muon.
                /*
                   for (Int_t i = 0; i < muon_n; i++) {
                   muon = (Muon*) branch_mu->At(i);
                   if (muon->PT > muon_pt) {  // Get the maximum PT muon.
                   muon_pt = muon->PT;
                   printf("new muon pt = %f \n", muon_pt);
                   muon_phi = muon->Phi;
                   muon_eta = muon->Eta;
                   }
                   }
                 */

                float sum_met, metx, mety;
                for (Int_t i = 0; i < met_n; i++) {
                 //   std::cout<<"i.."<<i<<":";
                    met = (MissingET*) branch_met->At(i);
                    sum_met = met->MET;
                    metx = met->MET*cos(met->Phi);
                    mety = met->MET*sin(met->Phi);

                }
                //      printf("MET = %f \n", sum_met);

                // Make cuts.

                // transverse mass of jet and MET i.e. MT variable


                if (btag_found && sum_met > 100 && non_b_jet_found) {

                    // Fill Histograms

float mt_value = mTCalculation(metx, mety,nonb_jet_px , nonb_jet_py, nonb_jet_pt);
                    h_jet_pt->Fill(nonb_jet_pt);
                    h_jet_eta->Fill(nonb_jet_eta);
                    h_jet_phi->Fill(nonb_jet_phi);

                    h_bjet_pt->Fill(btag_pt);
                    h_bjet_eta->Fill(btag_eta);
                    h_bjet_phi->Fill(btag_phi);

                    h_MET->Fill(sum_met);
                    h_MT->Fill(mt_value);


                }  // End if (btag_found)
            }  // End if (muon_n > 0 && jet_n > 0)

        }

        // Draw hist and save.
        // TCanvas *c = new TCanvas();

        h_jet_pt->Scale(1./h_jet_pt->Integral());
        h_jet_eta->Scale(1./h_jet_eta->Integral());
        h_jet_phi->Scale(1./h_jet_phi->Integral());
        h_bjet_pt->Scale(1./h_bjet_pt->Integral());
        h_bjet_eta->Scale(1./h_bjet_eta->Integral());
        h_bjet_phi->Scale(1./h_bjet_phi->Integral());
        h_MET->Scale(1./h_MET->Integral());
        h_MT->Scale(1./h_MT->Integral());


        TFile *f = new TFile(save_name,"RECREATE");
        h_jet_pt->Write("",TObject::kOverwrite);
        h_jet_eta->Write("",TObject::kOverwrite);
        h_jet_phi->Write("",TObject::kOverwrite);
        h_bjet_pt->Write("",TObject::kOverwrite);
        h_bjet_eta->Write("",TObject::kOverwrite);
        h_bjet_phi->Write("",TObject::kOverwrite);
        h_MET->Write("",TObject::kOverwrite);
        h_MT->Write("",TObject::kOverwrite);


    }
