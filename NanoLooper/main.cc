#include "main.h"
#include "Nano.h" // Contains the definition of the "nt" object that reads the NanoAOD
#include "MCTools.h" // Contains the definition of the dumpGenParticleInfos();
#include "ElectronSelections.h" // Contains the definitions of ttH::electronID
#include "MuonSelections.h" // Contains the definitions of ttH::muonID
#include "Config.h" // Contains the definitions of gconf
#include "rooutil.h" // Contains the definitions of various help functions that start with "RooUtil"
#include "cxxopts.h" // Contains how to parse argc and argv 


//=================================================================================================
// Objects
//=================================================================================================
namespace Obj
{

    //_______________________________________________________
    // Electron data structure
    struct Elec
    {
        LV p4;
        int isTight;
        int jetIdx; // Overlapping jet index
        int pdgid;
    };

    //_______________________________________________________
    // Muon data structure
    struct Muon
    {
        LV p4;
        int isTight;
        int jetIdx; // Overlapping jet index
        int pdgid;
    };


    //_______________________________________________________
    // Jet data structure
    struct Jet
    {
        LV p4;
        int JetIdx; // Overlapping jet index
        float hbbScore;
        float isBtagScore;
        int isBtagTight;
        int isBtagLoose;
        double mass;
        float tau1;
        float tau2;
        float tau3;
        float tau4;
        float n2b1;
        float n3b1;
        float wQCDScore;
        float zQCDScore;
    };
    
}

//=================================================================================================
// Key analysis variables
// NOTE: variables have "_" ending
//=================================================================================================
namespace Analysis
{

    //_______________________________________________________
    // Per event weight (normalized to 1/fb of data)
    Double_t scale1fb_;

    //_______________________________________________________
    // Integrated Luminosity
    Double_t lumi_;

    //_______________________________________________________
    // Per event weight (normalized * lumi)
    Double_t wgt_;

    //_______________________________________________________
    // Year of the data
    Int_t year_;

    //_______________________________________________________
    // Flags used in the analysis
    bool is_vbs;
    bool is_hbb;
    bool is_wjj;

    //_______________________________________________________
    // Signal's generator level kinematics
    LV gen_ijet_;
    LV gen_jjet_;
    LV gen_jet0_;
    LV gen_jet1_;
    std::vector<LV> gen_W_;
    std::vector<LV> gen_Z_;
    LV gen_H_;
    LV gen_b0_;
    LV gen_b1_;
    LV gen_Wb0_;
    LV gen_Wb1_;

    //_______________________________________________________
    // Reconstructed objects
    std::vector<Obj::Elec> elecs_;
    std::vector<Obj::Muon> muons_;
    std::vector<Obj::Jet> jets_;
    std::vector<LV> leptons_;
    std::vector<LV> VBFjets_;
    std::vector<LV> bjets_;
    LV leadingJet_;
    LV subleadingJet_;

    //________________________________________________________
    // Fat Jets reconstruction
    std::vector<Obj::Jet> fatJets_;
    std::vector<float> tau_(4);
    std::vector<float> mets_(3);
    LV hbbFatJet_;
    LV wjjFatJet_;
    float n2b1;
    float n3b1;
    float wvsQCD;
    float zvsQCD;
 

    //_______________________________________________________
    // Set per event weight (normalized to 1/fb of data)
    void setScale1fb(Double_t scale1fb) { scale1fb_ = scale1fb; }

    //_______________________________________________________
    // Set integrated luminosity
    void setLumi(Double_t lumi) { lumi_ = lumi; }

    //_______________________________________________________
    // Set the year of the data
    void setYear(Int_t year) { year_ = year; }

    //_______________________________________________________
    // Set up the NanoCORE's common configuration service tool
    void setConfig()
    {
        gconf.GetConfigs(nt.year());
    }

    //_______________________________________________________
    // Clear all analysis variables
    void clearAnalysisVariables()
    {
        // Flags used in the anlaysis
        is_vbs = false;
        is_hbb = false;
        is_wjj = false;

        // Signal's generator level kinematics
        gen_ijet_ = LV();
        gen_jjet_ = LV();
        gen_jet0_ = LV();
        gen_jet1_ = LV();
        gen_H_ = LV();
        gen_b0_ = LV();
        gen_b1_ = LV();
        gen_Wb0_ = LV();
        gen_Wb1_ = LV();
        

        // Reconstructed objects
        gen_W_.clear();
        gen_Z_.clear();
        elecs_.clear();
        muons_.clear();
        jets_.clear();
        fatJets_.clear();
        tau_.clear();
        mets_.clear();
        leptons_.clear();
        VBFjets_.clear();
        bjets_.clear();
        leadingJet_ = LV();
        subleadingJet_ = LV();
        hbbFatJet_ = LV(); 
        wjjFatJet_ = LV();
        n2b1 = 0;
        n3b1 = 0.5;
        wvsQCD = 0;
        zvsQCD = 0;

    }

    //_______________________________________________________
    // Compute the event weight for this event
    void computeEventWeight()
    {
        // There are some cases where the event weight is negative
        Double_t sign_of_the_weight = ((nt.Generator_weight() > 0) - (nt.Generator_weight() < 0));
        Double_t per_event_weight = sign_of_the_weight * Analysis::scale1fb_ * Analysis::lumi_;
        wgt_ = per_event_weight;
    }

    //_______________________________________________________
    // Compute signal generator level kinematics
    void computeSignalGeneratorLevelKinematics()
    {
        is_vbs = nt.GenPart_status()[2] == 23;
        if (is_vbs)
        {
            // LV ordered [NObjects];
            gen_ijet_ = nt.GenPart_p4()[2];
            gen_jjet_ = nt.GenPart_p4()[3];
            gen_jet0_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_ijet_ : gen_jjet_;
            gen_jet1_ = gen_ijet_.pt() > gen_jjet_.pt() ? gen_jjet_ : gen_ijet_;
            for (unsigned int i = 0; i < 7; ++i) {
                if (abs(nt.GenPart_pdgId()[i] == 24)) { gen_W_.push_back(nt.GenPart_p4()[i]);}
                if (abs(nt.GenPart_pdgId()[i] == 23)) { gen_Z_.push_back(nt.GenPart_p4()[i]);}
                if (abs(nt.GenPart_pdgId()[i] == 25)) { gen_H_ = nt.GenPart_p4()[i];} 
            }
            // Select b quarks and hbb
            std::vector<LV> h_decay_p4s;
            std::vector<int> h_decay_pdgIds;
            // std::vector<int> h_decay_statuses;
            for (unsigned int igen = 0; igen < nt.GenPart_p4().size(); ++igen)
            {
                int imom = nt.GenPart_genPartIdxMother()[igen];
                if (abs(nt.GenPart_pdgId()[imom]) == 25 and nt.GenPart_status()[imom] == 62)
                {
                    // int pdgid = nt.GenPart_pdgId()[igen];
                    // int status = nt.GenPart_status()[igen];
                    h_decay_p4s.push_back(nt.GenPart_p4()[igen]);
                    h_decay_pdgIds.push_back(nt.GenPart_pdgId()[igen]);
                    // h_decay_statuses.push_back(nt.GenPart_status()[igen]);
                }
            }
            if (h_decay_pdgIds.size() != 0) {
                if (abs(h_decay_pdgIds[0]) == 5) {
                    is_hbb = true;
                }
            }
            if (is_hbb) 
            {
                gen_b0_ = h_decay_p4s[0].pt() > h_decay_p4s[1].pt() ? h_decay_p4s[0] : h_decay_p4s[1];
                gen_b1_ = h_decay_p4s[0].pt() > h_decay_p4s[1].pt() ? h_decay_p4s[1] : h_decay_p4s[0];
            }

            // Select jet quarks and wjj
            std::vector<LV> w_decay_p4s;
            std::vector<int> w_decay_pdgIds;
            // std::vector<int> w_decay_statuses;
            for (unsigned int igen = 0; igen < nt.GenPart_p4().size(); ++igen)
            {
                int imom = nt.GenPart_genPartIdxMother()[igen];
                if (abs(nt.GenPart_pdgId()[imom]) == 24 and nt.GenPart_status()[imom] == 62)
                {
                    // int pdgid = nt.GenPart_pdgId()[igen];
                    // int status = nt.GenPart_status()[igen];
                    w_decay_p4s.push_back(nt.GenPart_p4()[igen]);
                    w_decay_pdgIds.push_back(nt.GenPart_pdgId()[igen]);
                    // w_decay_statuses.push_back(nt.GenPart_status()[igen]);
                }
            }
            if (w_decay_pdgIds.size() != 0) {
                if (abs(w_decay_pdgIds[0]) >= 1 and abs(w_decay_pdgIds[0]) <= 5 and (w_decay_pdgIds[0] * w_decay_pdgIds[1]) < 0){
                    is_wjj = true;
                }
            }
            if (is_wjj) 
            {
                if (w_decay_p4s.size() == 2) {
                    gen_Wb0_ = w_decay_p4s[0].pt() > w_decay_p4s[1].pt() ? w_decay_p4s[0] : w_decay_p4s[1];
                    gen_Wb1_ = w_decay_p4s[0].pt() > w_decay_p4s[1].pt() ? w_decay_p4s[1] : w_decay_p4s[0];
                }
            }
        }
    }

    //_______________________________________________________
    // Select MET
    void readMET()
    {
        mets_[0] = nt.MET_pt();
        mets_[1] = nt.MET_phi();
        mets_[2] = nt.MET_sumEt();
    }

   

    //_______________________________________________________
    // Select electrons
    void selectElectrons()
    {
        for (unsigned int iel = 0; iel < nt.Electron_pt().size(); ++iel)
        {
            if (ttH::electronID(iel, ttH::IDveto, nt.year()))
            {
                Obj::Elec this_elec;
                this_elec.p4 = nt.Electron_p4()[iel];
                this_elec.isTight = ttH::electronID(iel, ttH::IDtight, nt.year());
                this_elec.jetIdx = nt.Electron_jetIdx()[iel];
                this_elec.pdgid = nt.Electron_pdgId()[iel];
                elecs_.push_back(this_elec);
                leptons_.push_back(this_elec.p4);
            }
        }
    }

    //_______________________________________________________
    // Select muons
    void selectMuons()
    {
        // Select muons
        for (unsigned int imu = 0; imu < nt.Muon_pt().size(); ++imu)
        {
            if (ttH::muonID(imu, ttH::IDveto, nt.year()))
            {
                Obj::Muon this_muon;
                this_muon.p4 = nt.Muon_p4()[imu];
                this_muon.isTight = ttH::muonID(imu, ttH::IDtight, nt.year());
                this_muon.jetIdx = nt.Muon_jetIdx()[imu];
                this_muon.pdgid = nt.Muon_pdgId()[imu];
                muons_.push_back(this_muon);
                leptons_.push_back(this_muon.p4);
            }
        }
    }

    // _______________________________________________________
    // Select Fatjets
    void selectFatJets()
    {
        for (unsigned int ifatjet = 0; ifatjet < nt.FatJet_pt().size(); ++ifatjet)
        {
            Obj::Jet this_fatJet;
            // Overlap check against good electrons
            bool isOverlap = false;
            for (unsigned int ilep = 0; ilep < elecs_.size(); ++ilep)
            {
                if (RooUtil::Calc::DeltaR(elecs_[ilep].p4, nt.FatJet_p4()[ifatjet]) < 0.8)
                {
                    isOverlap = true;
                    break;
                }
            }

             // Overlap check against good muons
            for (unsigned int ilep = 0; ilep < muons_.size(); ++ilep)
            {
                if (RooUtil::Calc::DeltaR(muons_[ilep].p4, nt.FatJet_p4()[ifatjet]) < 0.8)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap 
            // Then skip
            if (isOverlap)
                continue;
           
            // Keep fat jets softdropmass above 40 GeV
            float boostedMass = nt.FatJet_msoftdrop()[ifatjet]; 
            this_fatJet.mass = boostedMass;
            if (not (boostedMass > 40.))
                continue;

            // p4
            TLorentzVector p4;
            p4.SetPtEtaPhiM(nt.FatJet_pt()[ifatjet], nt.FatJet_eta()[ifatjet], nt.FatJet_phi()[ifatjet], 
                            nt.FatJet_msoftdrop()[ifatjet]);
            this_fatJet.p4 = RooUtil::Calc::getLV(p4);
            this_fatJet.hbbScore = nt.FatJet_particleNet_HbbvsQCD()[ifatjet];
            this_fatJet.JetIdx = -999;
            this_fatJet.tau1 = nt.FatJet_tau1()[ifatjet];
            this_fatJet.tau2 = nt.FatJet_tau2()[ifatjet];
            this_fatJet.tau3 = nt.FatJet_tau3()[ifatjet];
            this_fatJet.tau4 = nt.FatJet_tau4()[ifatjet];
            this_fatJet.n2b1 = nt.FatJet_n2b1()[ifatjet];
            this_fatJet.n3b1 = nt.FatJet_n3b1()[ifatjet];
            this_fatJet.wQCDScore = nt.FatJet_particleNet_WvsQCD()[ifatjet];
            this_fatJet.zQCDScore = nt.FatJet_particleNet_ZvsQCD()[ifatjet];
            fatJets_.push_back(this_fatJet);
        }
        
        // FatJet Selection
        if (fatJets_.size() > 0) 
        { 
            float maxhbbscore = -999;
            int maxHbbNo;
        // Prioritize Hbb selection
            for (unsigned int ifatjet = 0; ifatjet < fatJets_.size(); ifatjet++) {
                if (fatJets_[ifatjet].hbbScore >= maxhbbscore) {
                    maxHbbNo = ifatjet;
                    maxhbbscore = fatJets_[ifatjet].hbbScore;
                }
            }
            hbbFatJet_ = fatJets_[maxHbbNo].p4;
            tau_[0] = fatJets_[maxHbbNo].tau1;
            tau_[1] = fatJets_[maxHbbNo].tau2;
            tau_[2] = fatJets_[maxHbbNo].tau3;
            tau_[3] = fatJets_[maxHbbNo].tau4;
            n2b1 = fatJets_[maxHbbNo].n2b1;
            n3b1 = fatJets_[maxHbbNo].n3b1;
            wvsQCD = fatJets_[maxHbbNo].wQCDScore;
            zvsQCD = fatJets_[maxHbbNo].zQCDScore;
        }
    }


    //_______________________________________________________
    // Select jets
    void selectJets()
    {
        // Loop over the jets
        for (unsigned int ijet = 0; ijet < nt.Jet_pt().size(); ++ijet)
        {

            // Read jet p4
            const LV& jet_p4 = nt.Jet_p4()[ijet];

            if (nt.year() == 2016)
            {
                if (nt.Jet_jetId()[ijet] < 1) // For 2016 apparently it's >= 1
                    continue;
            }
            else
            {
                if (nt.Jet_jetId()[ijet] < 2) // "Tight" ID requirement while for others >= 2
                    continue;
            }

            // if (nt.year() == 2017)
            // {
            //     if (nt.Jet_puId()[ijet] < 7) // For 2017 "111" (pass tight)
            //         continue;
            // }

            // Overlap check against good electrons
            bool isOverlap = false;
            for (unsigned int ilep = 0; ilep < elecs_.size(); ++ilep)
            {
                if (elecs_[ilep].jetIdx == (int) ijet)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap check against good muons
            for (unsigned int ilep = 0; ilep < muons_.size(); ++ilep)
            {
                if (muons_[ilep].jetIdx == (int) ijet)
                {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap check against good fatjets
            for (unsigned int ilep = 0; ilep < fatJets_.size(); ++ilep) {
                if (RooUtil::Calc::DeltaR(fatJets_[ilep].p4, jet_p4) <= 0.8) {
                    isOverlap = true;
                    break;
                }
            }

            // Overlap 
            // Then skip
            if (isOverlap)
                continue;

            // We keep jets above 30 GeV only
            if (not (jet_p4.pt() > 30.) )
                continue;

            Obj::Jet this_jet;
            this_jet.p4 = nt.Jet_p4()[ijet];
            this_jet.isBtagScore = nt.Jet_btagDeepFlavB()[ijet];
            this_jet.isBtagTight = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_tight;
            this_jet.isBtagLoose = nt.Jet_btagDeepFlavB()[ijet] > gconf.WP_DeepFlav_loose;
            jets_.push_back(this_jet);
            

        }
    }

        
    //___________________________________________________________
    // Select jet quarks
    void classifyJet()
    {
        
        if (jets_.size() == 0)
            return;
        double largestEta1 = -999;
        double largestEta2 = 999;
        double largestPt = 0;
        double largerPt = 0;
        int leading = 0;
        int subleading = 0;
        int VBFIndex1 = 0;
        int VBFIndex2 = 0;
        for (unsigned int i = 0; i < jets_.size(); i++) {
            LV jeti = jets_[i].p4;
            if (jeti.Eta() > 0 && jeti.Eta() >= largestEta1) {
                VBFIndex1 = i;
                largestEta1 = jeti.Eta();
            }
            else if (jeti.Eta() < largestEta2 && jeti.Eta() < 0) {
                VBFIndex2 = i;
                largestEta2 = jeti.Eta();
            }
            if (jeti.Pt() > largestPt) {
                leading = i;
                largestPt = jeti.Pt();
            }
            else if (jeti.Pt() < largestPt && jeti.Pt() > largerPt) {
                subleading = i;
                largerPt = jeti.Pt();
            }
        } 
        leadingJet_ = jets_[leading].p4;
        subleadingJet_ = jets_[subleading].p4;

        VBFjets_.push_back(jets_[VBFIndex1].p4);
        VBFjets_.push_back(jets_[VBFIndex2].p4);
    }

    //_______________________________________________________
    // Run the analysis for this event
    void runAlgorithms()
    {
        // Beginning of every event we should clear the analysis variables
        Analysis::clearAnalysisVariables();

        // Compute per event weight
        Analysis::computeEventWeight();

        // Compute signal generator level kinematics
        Analysis::computeSignalGeneratorLevelKinematics();

        // Select electrons and muons
        Analysis::readMET();
        Analysis::selectElectrons();
        Analysis::selectMuons();
        Analysis::selectFatJets();
        Analysis::selectJets();
        Analysis::classifyJet();

    }

}

//=================================================================================================
// Cutflows
// Naming convention: "h_NAME_" (i.e. Starts with "h_" and ends with "_")
//=================================================================================================
namespace Cutflow
{
    // cutflow
    TH1F* h_cutflow_;
    TH1F* h_cutflow_unwgt_;

    //_______________________________________________________
    // Cutflow
    enum Cuts
    {
        kNoSelection = 0,
        kVbs,
        kHbb,
        kWjj,
        kTwoLightLeptons,
        kOneHbbFatJet,
        kN3B1,
        kHbbScore,
        kAtLeastTwoPt30Jets,
        kdRVBF,
        kNCuts,
    };



    //_______________________________________________________
    // Book the histograms
    void bookCutflow()
    {
        // Cutflow histogram for including gen_level
        h_cutflow_ = new TH1F("h_cutflow", "Cutflow of the selections", Cuts::kNCuts, 0, Cuts::kNCuts);

    }

    //_______________________________________________________
    // Fill cutflow
    void fillCutflow(Cuts cutstage)
    {
        h_cutflow_->Fill(cutstage, Analysis::wgt_);
    } 

    //_______________________________________________________
    // write cutflow
    void writeCutflow(TFile* ofile)
    {
        ofile->cd();
        h_cutflow_->Write();
    } 

}

//=================================================================================================
// Histograms
// Naming convention: "h_NAME_" (i.e. Starts with "h_" and ends with "_")
//=================================================================================================
namespace Hist
{

    // pt
    TH1F* h_gen_ptZ_;
    TH1F* h_gen_ptW_;
    TH1F* h_gen_ptH_;

    // VBF quarks
    TH1F* h_gen_deltaEta_;
    TH1F* h_gen_massVBF_;

    // b quarks
    TH1F* h_gen_massbQsystem_;
    TH1F* h_gen_ptb0_;
    TH1F* h_gen_ptb1_;
    
    // MET_
    TH1F* ptMET_;
    TH1F* phiMET_;
    TH1F* EtMET_;

    // leptons_
    TH1F* ptLep1;
    TH1F* ptLep2;
    TH1F* re_deltaEtaVBF;
    TH1F* etaLep1;
    TH1F* etaLep2;
    TH1F* phiLep1;
    TH1F* phiLep2;
    TH1F* massDiLep;
    TH1F* ptDiLep;
    TH1F* dRLep; 
    TH1F* deltaEtaLep;
    TH1F* deltaEtaLep1FatJet_;
    TH1F* deltaEtaLep2FatJet_;
    TH1F* nLep;

    // jets
    TH1F* nJets_;
    TH1F* ptDijet;
    TH1F* etaVBFjet1;
    TH1F* etaVBFjet2;
    TH1F* ptVBFjet1;
    TH1F* ptVBFjet2;
    TH1F* mpRatioVBFjet1_;
    TH1F* mpRatioVBFjet2_;
    TH1F* dRVBF;
    TH1F* VBFjetMass;
    TH1F* ptLeading_;
    TH1F* ptSubleading_;
    TH1F* etaLeading_;
    TH1F* dRjets_;
    TH1F* deltaEtajets_;
    TH1F* etaSubleading_;
    TH1F* phiLeading_;
    TH1F* phiSubleading_;
    TH1F* dRVBF1FatJet_;
    TH1F* dRVBF2FatJet_;
    TH1F* deltaEtaVBF1FatJet_;
    TH1F* deltaEtaVBF2FatJet_;
    
    // dR between reconstructed and gen-level
    TH1F* dRleadingVBF_;
    TH1F* dRsubleadingVBF;
    TH1F* dRsubleadingbq;
    TH1F* dRleadingbq;
    TH1F* dRhiggs;

    // FatJets
    TH1F* ptHbb_;
    TH1F* etaHbb_;
    TH1F* phiHbb_;
    TH1F* massHbb;
    TH1F* nFatjets_;
    TH1F* softdropmass_;

    // FatJets substructure observables
    TH1F* tau1_;
    TH1F* tau2_;
    TH1F* tau3_;
    TH1F* tau4_;
    TH1F* tau21_;
    TH1F* tau32_;
    TH1F* tau41_;
    TH1F* tau31_;
    TH1F* tau42_;
    TH1F* tau43_;
    TH1F* n2b1_;
    TH1F* n3b1_;

    // FatJetScore
    TH1F* hbbScore_;
    TH1F* hjetScore_;
    TH1F* zScore_;
    TH1F* wScore_;

    // s-hat variable
    TH1F* massZH_;
    TH1F* massZHzoom_;
    TH1F* ST_;
    TH1F* KT_;

    //_______________________________________________________
    // Book the histograms
    void bookHistograms()
    {
        // dR
        dRleadingVBF_ = new TH1F("dRleadingVBF", "delta R of leading VBF quarks", 1080, 0, 5);
        dRsubleadingbq = new TH1F("dRsubleadingbq", "delta R of subleading bottom quarks", 1080, 0, 5);
        dRsubleadingVBF = new TH1F("dRsubleadingVBF", "delta R of subleading VBF quarks", 1080, 0, 5);
        dRleadingbq = new TH1F("dRleadingbq", "delta R of leading bottom quarks", 1080, 0, 5);
        dRhiggs = new TH1F("dRhiggs", "delta R of gen higgs and tagged reco hbbfatjet", 1080, 0, 5);
        // Pt's of the bosons
        h_gen_ptZ_ = new TH1F("h_gen_ptZ", "Transverse momentum of Z bosons", 1080, 0, 1800);
        h_gen_ptW_ = new TH1F("h_gen_ptW", "Transverse momentum of W bosons", 1080, 0, 1800);
        h_gen_ptH_ = new TH1F("h_gen_ptH", "Transverse momentum of H bosons", 1080, 0, 1800);
        // VBF quarks kinematics
        h_gen_deltaEta_ = new TH1F("h_gen_deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
        h_gen_massVBF_ = new TH1F("h_gen_massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);

        // bQ kinematics
        h_gen_massbQsystem_ = new TH1F("h_gen_massbQsystem", "Mass of bottom quark system", 1080, 0, 200);
        h_gen_ptb0_ = new TH1F("h_gen_ptb0", "pt of leading bottom quarks", 1080, 0, 400);
        h_gen_deltaEta_ = new TH1F("h_gen_deltaEta", "Delta Eta of VBF quarks", 1080, 0, 10);
        h_gen_massVBF_ = new TH1F("h_gen_massVBF", "Invariant mass of VBF quark system", 1080, 0, 3500);

        //______________________________________________________
        // MET
        ptMET_ = new TH1F("ptMET", "missing transverse momentum", 1080, 0, 1000);
        phiMET_ = new TH1F("phiMET", "missing transverse momentum phi", 1080, -5, 5);
        EtMET_ = new TH1F("EtMET", "missing transverse momentum Et", 1080, 0, 4500);
        
        //______________________________________________________
        // Leptons kinematics
        ptLep1 = new TH1F("ptLep1", "pt of leading leptons_", 1080, 0, 600);
        ptLep2 = new TH1F("ptLep2", "pt of subleading leptons_", 1080, 0, 600);
        etaLep1 = new TH1F("etaLep1", "eta of leading leptons_", 1080, -5, 5);
        etaLep2 = new TH1F("etaLep2", "eta of subleading leptons_", 1080, -5, 5);
        phiLep1 = new TH1F("phiLep1", "phi of leading leptons_", 1080, -5, 5);
        phiLep2 = new TH1F("phiLep2", "phi of subleading leptons_", 1080, -5, 5);
        massDiLep = new TH1F("massDiLep", "mass of dileptons", 1080, 0, 400);
        ptDiLep = new TH1F("ptDiLep", "pt of the dilepton system", 1080, 0, 2500);
        dRLep = new TH1F("dRLep", "delta R between two leptons", 1080, 0, 5);
        deltaEtaLep = new TH1F("deltaEtaLep", "delta Eta between two leptons", 1080, 0, 5);
        deltaEtaLep1FatJet_ = new TH1F("deltaEtaLep1FatJet", "delta Eta between lepton1 and Fatjet", 1080, 0, 5);
        deltaEtaLep2FatJet_ = new TH1F("deltaEtaLep2FatJet", "delta Eta between lepton2 and Fatjet", 1080, 0, 5);
        nLep = new TH1F("nLep", "number of leptons", 5, 0, 5);

        //______________________________________________________
        // jets kinematics
        etaVBFjet1 = new TH1F("etaVBFjet1", "eta of leading VBF jet", 1080, -5, 5); 
        etaVBFjet2 = new TH1F("etaVBFjet2", "eta of subleading VBF jet", 1080, -5, 5); 
        re_deltaEtaVBF = new TH1F("re_deltaEtaVBF", "Delta Eta of VBF quarks (Reconstructed)", 1080, 0, 10);
        ptVBFjet1 = new TH1F("ptVBFjet1", "pt of leading VBF jets", 1080, 0, 500);
        ptVBFjet2 = new TH1F("ptVBFjet2", "pt of subleading VBF jets", 1080, 0, 500);
        mpRatioVBFjet1_ = new TH1F("mpRatioVBFjet1", "mass/pt of VBF jet1", 1080, 0, 1);
        mpRatioVBFjet2_ = new TH1F("mpRatioVBFjet2", "mass/pt of VBF jet2", 1080, 0, 1);
        ptDijet = new TH1F("ptDijet", "pt of dijet", 1080, 0, 1000);
        VBFjetMass = new TH1F("VBFjetMass", "Invariant mass of VBF quark system", 1080, 0, 3500);
        dRVBF = new TH1F("dRVBF", "delta R of VBF jets", 1080, 0, 10);
        dRVBF1FatJet_ = new TH1F("dRVBF1FatJet", "delta R between VBF jet1 and the Hbb fatjet", 100, 0, 10);
        dRVBF2FatJet_ = new TH1F("dRVBF2FatJet", "delta R between VBF jet2 and the Hbb fatjet", 100, 0, 10);
        deltaEtaVBF1FatJet_ = new TH1F("deltaEtaVBF1FatJet", "delta eta between VBF jet1 and the Hbb fatjet", 100, 0, 10);
        deltaEtaVBF2FatJet_ = new TH1F("deltaEtaVBF2FatJet", "delta eta between VBF jet2 and the Hbb fatjet", 100, 0, 10);
        ptLeading_ = new TH1F("ptLeading", "pt of leading jet quarks", 1080, 0, 500);
        ptSubleading_ = new TH1F("ptSubleading", "pt of subleading jet quarks", 1080, 0, 500);
        etaLeading_ = new TH1F("etaLeading", "eta of leading jet quarks", 1080, -5, 5);
        dRjets_ = new TH1F("dRjets", "delta R of two leading jets", 1080, 0, 5);
        deltaEtajets_ = new TH1F("deltaEtajets", "delta Eta of two leading jets", 1080, 0, 10);
        etaSubleading_ = new TH1F("etaSubleading", "eta of subleading jet", 1080, -5, 5);
        phiLeading_ = new TH1F("phiLeading", "phi of leading jet", 1080, -5, 5);
        phiSubleading_ = new TH1F("phiSubleading", "phi of subleading jet", 1080, -5, 5);
        nJets_ = new TH1F("nJets", "number of ak4 jets", 8, 0, 8);

        //______________________________________________________
        // fatJets Kinematics
        ptHbb_ = new TH1F("ptHbb", "pt of Hbb fatjet", 1080, 0, 1800);
        etaHbb_ = new TH1F("etaHbb", "eta of Hbb fatjet", 1080, -5, 5);
        phiHbb_ = new TH1F("phiHbb", "phi of Hbb fatjet", 1080, -5, 5);
        massHbb = new TH1F("massHbb", "mass of Hbb fatjet", 1080, 0, 200);
        nFatjets_ = new TH1F("nFatjets", "Number of fatjets", 5, 0, 5);
        softdropmass_ = new TH1F("softdropmass", "softdrop mass of hbb fatjets", 1080, 0, 200);

        //______________________________________________________
        // fatJets substructure observables
        tau1_ = new TH1F("tau1", "Nsubjettiness 1 axis", 1080, 0, 1); 
        tau2_ = new TH1F("tau2", "Nsubjettiness 2 axis", 1080, 0, 1); 
        tau3_ = new TH1F("tau3", "Nsubjettiness 3 axis", 1080, 0, 1); 
        tau4_ = new TH1F("tau4", "Nsubjettiness 4 axis", 1080, 0, 1); 
        tau21_ = new TH1F("tau21", "Nsubjettiness ratio", 1080, 0, 1); 
        tau32_ = new TH1F("tau32", "Nsubjettiness ratio 32", 1080, 0, 1); 
        tau41_ = new TH1F("tau41", "Nsubjettiness ratio 41", 1080, 0, 1); 
        tau31_ = new TH1F("tau31", "Nsubjettiness ratio 31", 1080, 0, 1); 
        tau42_ = new TH1F("tau42", "Nsubjettiness ratio 42", 1080, 0, 1); 
        tau43_ = new TH1F("tau43", "Nsubjettiness ratio 43", 1080, 0, 1); 
        n2b1_ = new TH1F("n2b1", "N2 with b1", 1080, -5, 5); 
        n3b1_ = new TH1F("n3b1", "N3 with b1", 1080, -5, 5); 

        //______________________________________________________
        // fatJets Score
        hbbScore_ = new TH1F("hbbScore", "hbb score of fatjets", 1080, 0, 1);
        hjetScore_ = new TH1F("hjetScore", "the hbb score of the selected hbb fatjet", 1080, 0, 1);
        wScore_ = new TH1F("wScore", "the particleNet w score of the fatjet with highest hbb", 1080, 0, 1);
        zScore_ = new TH1F("zScore", "the particleNet z score of the fatjet with highest hbb", 1080, 0, 1);

        //_______________________________________________________
        // s-hat variable
        massZH_ = new TH1F("massZH", "Invariant mass of the ZH system", 1080, 0, 3500);
        massZHzoom_ = new TH1F("massZHzoom", "Invariant mass of the ZH system", 1080, 0, 500);
        ST_ = new TH1F("ST", "Scalar sum of VVH system", 1080, 0, 3500);
        KT_ = new TH1F("KT", "Scalar sum of Hbb pt and VBF jets pt", 1080, 0, 4500);
    }

    //_______________________________________________________
    // Fill fatjet Histograms
    void fillHbbFatJetKinematicHistograms()
    {
        ptHbb_->Fill(Analysis::hbbFatJet_.Pt(), Analysis::wgt_);
        etaHbb_->Fill(Analysis::hbbFatJet_.Eta(), Analysis::wgt_);
        phiHbb_->Fill(Analysis::hbbFatJet_.Phi(), Analysis::wgt_);
        massHbb->Fill(Analysis::hbbFatJet_.M(), Analysis::wgt_);
        nFatjets_->Fill(Analysis::fatJets_.size(), Analysis::wgt_);
        for (unsigned int j = 0; j < Analysis::fatJets_.size(); j++) {
            hbbScore_->Fill(Analysis::fatJets_[j].hbbScore, Analysis::wgt_);
            softdropmass_->Fill(Analysis::fatJets_[j].mass, Analysis::wgt_);
        }
        tau1_->Fill(Analysis::tau_[0], Analysis::wgt_);
        tau2_->Fill(Analysis::tau_[1], Analysis::wgt_);
        tau3_->Fill(Analysis::tau_[2], Analysis::wgt_);
        tau4_->Fill(Analysis::tau_[3], Analysis::wgt_);
        tau21_->Fill(Analysis::tau_[1]/Analysis::tau_[0], Analysis::wgt_);
        tau32_->Fill(Analysis::tau_[2]/Analysis::tau_[1], Analysis::wgt_);
        tau41_->Fill(Analysis::tau_[3]/Analysis::tau_[0], Analysis::wgt_);
        tau31_->Fill(Analysis::tau_[2]/Analysis::tau_[0], Analysis::wgt_);
        tau42_->Fill(Analysis::tau_[3]/Analysis::tau_[1], Analysis::wgt_);
        tau43_->Fill(Analysis::tau_[3]/Analysis::tau_[2], Analysis::wgt_);
        n2b1_->Fill(Analysis::n2b1, Analysis::wgt_);
        n3b1_->Fill(Analysis::n3b1, Analysis::wgt_);

        float max_hbbScore = -999;
        for (unsigned int j = 0; j < Analysis::fatJets_.size(); j++) {
            if (max_hbbScore < Analysis::fatJets_[j].hbbScore)
            {
                max_hbbScore = Analysis::fatJets_[j].hbbScore;
            }
        }
        hjetScore_->Fill(max_hbbScore, Analysis::wgt_);
        wScore_->Fill(Analysis::wvsQCD, Analysis::wgt_);
        zScore_->Fill(Analysis::zvsQCD, Analysis::wgt_);
    }

    //_______________________________________________________
    // Fill MET histograms
    void fillMETHistograms() {
        ptMET_->Fill(Analysis::mets_[0], Analysis::wgt_);
        phiMET_->Fill(Analysis::mets_[1], Analysis::wgt_);
        EtMET_->Fill(Analysis::mets_[2], Analysis::wgt_);
    }

    //_______________________________________________________
    // Fill Lepton histograms
    void fillLeptonsKinematicHistograms()
    {
        LV lep1 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[0] : Analysis::leptons_[1];
        LV lep2 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[1] : Analysis::leptons_[0];
        ptLep1->Fill(lep1.Pt(), Analysis::wgt_);
        ptLep2->Fill(lep2.Pt(), Analysis::wgt_);
        etaLep1->Fill(lep1.Eta(), Analysis::wgt_);
        etaLep2->Fill(lep2.Eta(), Analysis::wgt_);
        phiLep1->Fill(lep1.Phi(), Analysis::wgt_);
        phiLep2->Fill(lep2.Phi(), Analysis::wgt_);
        massDiLep->Fill((lep1+lep2).M(), Analysis::wgt_);
        ptDiLep->Fill((lep1+lep2).Pt(), Analysis::wgt_);
        dRLep->Fill(RooUtil::Calc::DeltaR(lep1, lep2), Analysis::wgt_);
        deltaEtaLep->Fill(TMath::Abs(lep1.Eta()-lep2.Eta()), Analysis::wgt_);
        deltaEtaLep1FatJet_->Fill(TMath::Abs(lep1.Eta()-Analysis::hbbFatJet_.Eta()), Analysis::wgt_);
        deltaEtaLep2FatJet_->Fill(TMath::Abs(lep2.Eta()-Analysis::hbbFatJet_.Eta()), Analysis::wgt_);
        nLep->Fill(Analysis::leptons_.size(), Analysis::wgt_);
    }



    void fillJetsKinematicHistograms()
    {
        LV VBF1 = (Analysis::VBFjets_[0]).Pt() > (Analysis::VBFjets_[1]).Pt() ? Analysis::VBFjets_[0] : Analysis::VBFjets_[1];
        LV VBF2 = (Analysis::VBFjets_[0]).Pt() > (Analysis::VBFjets_[1]).Pt() ? Analysis::VBFjets_[1] : Analysis::VBFjets_[0];
        ptLeading_->Fill(Analysis::leadingJet_.Pt(), Analysis::wgt_);
        ptSubleading_->Fill(Analysis::subleadingJet_.Pt(), Analysis::wgt_);
        etaLeading_->Fill(Analysis::leadingJet_.Eta(), Analysis::wgt_);
        dRjets_->Fill(RooUtil::Calc::DeltaR(Analysis::leadingJet_, Analysis::subleadingJet_), Analysis::wgt_);
        deltaEtajets_->Fill(TMath::Abs(Analysis::leadingJet_.Eta()-Analysis::subleadingJet_.Eta()), Analysis::wgt_);
        etaSubleading_->Fill(Analysis::subleadingJet_.Eta(), Analysis::wgt_);
        phiLeading_->Fill(Analysis::leadingJet_.Phi(), Analysis::wgt_);
        phiSubleading_->Fill(Analysis::subleadingJet_.Phi(), Analysis::wgt_);

        etaVBFjet1->Fill(VBF1.Eta(), Analysis::wgt_);
        etaVBFjet2->Fill(VBF2.Eta(), Analysis::wgt_);
        ptVBFjet1->Fill(VBF1.Pt(), Analysis::wgt_);
        ptVBFjet2->Fill(VBF2.Pt(), Analysis::wgt_);
        mpRatioVBFjet1_->Fill(VBF1.M()/VBF1.Pt(), Analysis::wgt_);
        mpRatioVBFjet2_->Fill(VBF2.M()/VBF2.Pt(), Analysis::wgt_);


        nJets_->Fill(Analysis::jets_.size(), Analysis::wgt_);
        dRVBF->Fill(RooUtil::Calc::DeltaR(VBF1, VBF2), Analysis::wgt_);
        dRVBF1FatJet_->Fill(RooUtil::Calc::DeltaR(VBF1, Analysis::hbbFatJet_), Analysis::wgt_);
        dRVBF2FatJet_->Fill(RooUtil::Calc::DeltaR(VBF2, Analysis::hbbFatJet_), Analysis::wgt_);
        deltaEtaVBF1FatJet_->Fill(TMath::Abs(VBF1.Eta()-Analysis::hbbFatJet_.Eta()), Analysis::wgt_);
        deltaEtaVBF2FatJet_->Fill(TMath::Abs(VBF2.Eta()-Analysis::hbbFatJet_.Eta()), Analysis::wgt_);
        ptDijet->Fill((VBF1+VBF2).Pt(), Analysis::wgt_);
        re_deltaEtaVBF->Fill(TMath::Abs(Analysis::VBFjets_[0].Eta()-Analysis::VBFjets_[1].Eta()), Analysis::wgt_);
        VBFjetMass->Fill((Analysis::VBFjets_[0] + Analysis::VBFjets_[1]).M(), Analysis::wgt_);
        
        

        // Figuring out which combination of VBF jets is good
        float dR11 = RooUtil::Calc::DeltaR(VBF1, Analysis::gen_jet0_);
        float dR22 = RooUtil::Calc::DeltaR(VBF2, Analysis::gen_jet1_);
        float dR12 = RooUtil::Calc::DeltaR(VBF2, Analysis::gen_jet0_);
        float dR21 = RooUtil::Calc::DeltaR(VBF1, Analysis::gen_jet1_);

        LV VBF1_matched = dR11 + dR22 < dR12 + dR21 ? VBF1 : VBF2;
        LV VBF2_matched = dR11 + dR22 < dR12 + dR21 ? VBF2 : VBF1;

        // dR between reco- and gen- level
        dRleadingVBF_->Fill(RooUtil::Calc::DeltaR(VBF1_matched, Analysis::gen_jet0_), Analysis::wgt_);
        dRsubleadingVBF->Fill(RooUtil::Calc::DeltaR(VBF2_matched, Analysis::gen_jet1_), Analysis::wgt_);
        dRhiggs->Fill(RooUtil::Calc::DeltaR(Analysis::hbbFatJet_, Analysis::gen_H_), Analysis::wgt_);

    }
    

    // Fill gen-level histograms
    void fillGenLevelHistograms()
    {
        h_gen_massVBF_->Fill((Analysis::gen_ijet_ + Analysis::gen_jjet_).M(), Analysis::wgt_);
        h_gen_deltaEta_->Fill(TMath::Abs(Analysis::gen_ijet_.Eta() - Analysis::gen_jjet_.Eta()), Analysis::wgt_);
        for (unsigned int j = 0; j < Analysis::gen_Z_.size(); j++) { h_gen_ptZ_->Fill(Analysis::gen_Z_[j].Pt(), Analysis::wgt_);}
        for (unsigned int i = 0; i < Analysis::gen_W_.size(); i++) { h_gen_ptW_->Fill(Analysis::gen_W_[i].Pt(), Analysis::wgt_);}
        h_gen_ptH_->Fill(Analysis::gen_H_.Pt(), Analysis::wgt_);
        h_gen_massbQsystem_->Fill((Analysis::gen_b0_+Analysis::gen_b1_).M(), Analysis::wgt_); 
        h_gen_ptb0_->Fill(Analysis::gen_b0_.Pt(), Analysis::wgt_);
        h_gen_ptb1_->Fill(Analysis::gen_b1_.Pt(), Analysis::wgt_);
    }

    // Fill s-hat kinematic variables
    void fillSHatHistograms()
    {
        massZH_->Fill((Analysis::leptons_[0] + Analysis::leptons_[1] + Analysis::hbbFatJet_).M(), Analysis::wgt_);
        massZHzoom_->Fill((Analysis::leptons_[0] + Analysis::leptons_[1] + Analysis::hbbFatJet_).M(), Analysis::wgt_);
        ST_->Fill(Analysis::leptons_[0].pt() + Analysis::leptons_[1].pt() + Analysis::hbbFatJet_.pt() + Analysis::wjjFatJet_.pt(), Analysis::wgt_);
        KT_->Fill(Analysis::VBFjets_[0].pt() + Analysis::VBFjets_[1].pt() + Analysis::hbbFatJet_.pt(), Analysis::wgt_);
    }
    

 
    //_______________________________________________________
    // Write the histograms
    void writeHistograms(TFile* ofile, bool gen_level)
    {

        // Done with the analysis
        ofile->cd();

        if (gen_level)
        {
        // Gen-level
            h_gen_ptb0_->Write();
            h_gen_ptb1_->Write();
            h_gen_massbQsystem_->Write();
            h_gen_massVBF_->Write();
            h_gen_deltaEta_->Write();
            h_gen_ptH_->Write();
            h_gen_ptZ_->Write();
            h_gen_ptW_->Write();
        // difference between reco- and gen-
            dRleadingVBF_->Write();
            dRsubleadingVBF->Write();
            dRsubleadingbq->Write();
            dRleadingbq->Write();
            dRhiggs->Write();
        }
        softdropmass_->Write();
        tau1_->Write();
        tau2_->Write();
        tau3_->Write();
        tau4_->Write();
        tau21_->Write();
        tau32_->Write();
        tau41_->Write();
        tau43_->Write();
        tau42_->Write();
        tau31_->Write();
        n2b1_->Write();
        n3b1_->Write();
        ptMET_->Write();
        phiMET_->Write();
        EtMET_->Write();
        ptLep1->Write();
        ptLep2->Write();
        re_deltaEtaVBF->Write();
        etaLep1->Write();
        etaLep2->Write();
        phiLep1->Write();
        phiLep2->Write();
        massDiLep->Write();
        ptDiLep->Write();
        dRLep->Write();
        deltaEtaLep->Write();
        deltaEtaLep1FatJet_->Write();
        deltaEtaLep2FatJet_->Write();
        nLep->Write();
        nJets_->Write();
        etaVBFjet1->Write();
        etaVBFjet2->Write();
        ptVBFjet2->Write();
        ptVBFjet1->Write();
        mpRatioVBFjet1_->Write();
        mpRatioVBFjet2_->Write();
        VBFjetMass->Write();
        ptDijet->Write();
        dRVBF->Write();
        dRVBF1FatJet_->Write();
        dRVBF2FatJet_->Write();
        deltaEtaVBF1FatJet_->Write();
        deltaEtaVBF2FatJet_->Write();
        ptLeading_->Write();
        ptSubleading_->Write();
        etaLeading_->Write();
        dRjets_->Write();
        deltaEtajets_->Write();
        etaSubleading_->Write();
        phiLeading_->Write();
        phiSubleading_->Write();
        etaHbb_->Write();
        phiHbb_->Write();
        massHbb->Write();
        
        // fatJets Kinematics
        ptHbb_->Write();
        nFatjets_->Write();
        hbbScore_->Write();
        hjetScore_->Write();
        wScore_->Write();
        zScore_->Write();
        massZH_->Write();
        ST_->Write();
        KT_->Write();
        massZHzoom_->Write();
    }

}


        
        
//====================================================================
// Main function
//====================================================================
int main(int argc, char** argv) 
{
    // Setting 
    std::string input_path, output_path;
    float scale1fb;
    bool gen_level;
    int n_events;
     
    // Set the year 
    Analysis::setYear(2018); // TODO: Update properly in the future. For now it's a Placeholder!
    nt.SetYear(2018);
 
    // Grand option setting
    cxxopts::Options options("\n  $ doAnalysis",  "\n         **********************\n         *                    *\n         *       Looper       *\n         *                    *\n         **********************\n");

    // Read the options
    options.add_options()
        ("i,input"       , "Comma separated input file list OR if just a directory is provided it will glob all in the directory BUT must end with '/' for the path", cxxopts::value<std::string>())
        ("o,output"      , "Output file name"                                                                                    , cxxopts::value<std::string>())
        ("n,nevents"     , "N events to loop over"                                                                               , cxxopts::value<int>()->default_value("-1"))
        ("s,scale1fb"    , "pass scale1fb of the sample"                                                                         , cxxopts::value<float>())
    	("g,gen_level"   , "pass boolean for whether a gen-level analysis should be done"					 , cxxopts::value<bool>()->default_value("false"))
        ("h,help"        , "Print help")
        ;

    auto result = options.parse(argc, argv);

    //_______________________________________________________________________________
    // --help
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    } 

    //_______________________________________________________________________________
    // --input 
    if (result.count("input"))
    {
	input_path = result["input"].as<std::string>();
    }
    else 
    {
	std::cout << options.help() << std::endl;
	std::cout << "No input provided" << std::endl;
        exit(1);
    } 
   
    //_______________________________________________________________________________
    // --nevents
    n_events = result["nevents"].as<int>();

    //_______________________________________________________________________________
    // --gen_level
    gen_level = result["gen_level"].as<bool>();

    //_______________________________________________________________________________
    // --scale1fb
    if (result.count("scale1fb"))
    {
        scale1fb = result["scale1fb"].as<float>();
    }

    //_______________________________________________________________________________
    // --output
    if (result.count("output"))
    {
	output_path = result["output"].as<std::string>();
    }
    else 
    {
	std::cout << options.help() << std::endl;
	std::cout << "No output provided" << std::endl;
        exit(1);
    } 
    


    // Create your output root file
    TFile* output_file = new TFile(output_path.c_str(), "recreate");


    // Create Histograms
    Hist::bookHistograms();

    // Create Cutflow Histogram
    Cutflow::bookCutflow();

    // Set scale 1fb (the per event weight normalized for 1/fb)
    Analysis::setScale1fb(scale1fb); 

    // Set the luminosity
    Analysis::setLumi(137); // TODO: Update properly in the future. For now it's a Placeholder!

   // Set up year dependent configuration of the analysis that are POG specific from CMS
    Analysis::setConfig();

    // Input Path (RooUtil::Looper can accept comma separated list)
    // input_file_paths.push_back(input_path.c_str());
    // Use RooUtil::StringUtil to concatenate the input file paths
    // TString input = RooUtil::StringUtil::join(input_file_paths, ",");

    // Name of the TTree in the input root files.
    TString input_tree_name = "Events";

    // Create the TChain that holds the TTree's of the NanoAOD ntuples.
    TChain* events_tchain = RooUtil::FileUtil::createTChain(input_tree_name, input_path);

    // Create a looper that will handle event looping
    RooUtil::Looper<Nano> looper;

    // Initializer the looper
    looper.init(events_tchain, &nt, n_events);

    // Tracking cutflow number;
    int cut1_events = 0;
    int cut2_events = 0;
    int cut3_events = 0;
    int cut4_events = 0;
    int cut5_events = 0;
    int cut6_events = 0;
    int cut7_events = 0;
    int cut8_events = 0;
    
    // Loop through events
    while (looper.nextEvent())
    {
	    // dumpGenParticleInfos();
          
        // Run the analysis algorithms (selections, computing variables, etc.)
        Analysis::runAlgorithms();
        // Fill the "counter" histogram
        Cutflow::fillCutflow(Cutflow::Cuts::kNoSelection);
        
	
        if (gen_level) {
            // Cut#1: Check if vbs
            if (Analysis::is_vbs == false) { continue;}
            Cutflow::fillCutflow(Cutflow::Cuts::kVbs);
            cut1_events ++;
            // Cut#2: Check if hbb
            if (Analysis::is_hbb == false) { continue;}
            Cutflow::fillCutflow(Cutflow::Cuts::kHbb);
            cut2_events ++;
            // Cut#3: Check if wjj
            if (Analysis::is_wjj == false) { continue;}
            Cutflow::fillCutflow(Cutflow::Cuts::kWjj);
            cut3_events ++;
        }

        
        // Cut#3: Require that there are exactly two leptons
        if (not (Analysis::leptons_.size() ==2 )) {continue; }
        if (not ((Analysis::elecs_.size() == 2) || (Analysis::muons_.size() == 2))) { continue; }

        // Cut#3: Require that the two leptons have OS (opposite-sign)
        int is_os = false;
        if (Analysis::elecs_.size() == 2)
        {
            is_os = Analysis::elecs_[0].pdgid * Analysis::elecs_[1].pdgid < 0;
        }
        else if (Analysis::muons_.size() == 2)
        {
            is_os = Analysis::muons_[0].pdgid * Analysis::muons_[1].pdgid < 0;
        }
        else if (Analysis::muons_.size() == 1 && Analysis::elecs_.size() == 1)
        {
            is_os = Analysis::muons_[0].pdgid * Analysis::elecs_[0].pdgid < 0;
        }

        if (not (is_os)) { continue; }
        // Cut#3: Require that the leading leptons Pt >= 40, subleading >= 30
        LV leading = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[0] : Analysis::leptons_[1];
        LV subleading = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[1] : Analysis::leptons_[0];
            
        if (not (leading.Pt() >= 40 && subleading.Pt() >= 30)) {continue;}
        if (not (TMath::Abs(leading.Eta()<2.4))) {continue; } 

        Cutflow::fillCutflow(Cutflow::Cuts::kTwoLightLeptons);
        cut4_events ++;
        
        
        // Cut#4: Require at least one fatjet with softdropmass > 40 GeV
        if (not (Analysis::fatJets_.size() >= 1 ) ) { continue;}
        cut5_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kOneHbbFatJet);

        // Cut#5: Require n3b1 > 0.8
        if (not (Analysis::n3b1 > 0.8)) {continue;}
        Cutflow::fillCutflow(Cutflow::Cuts::kN3B1);

        // Cut#6: Require that there are at least 2 pt > 30 GeV jets
        if (not (Analysis::jets_.size() >= 2)) { continue; }
        Cutflow::fillCutflow(Cutflow::Cuts::kAtLeastTwoPt30Jets);
        cut7_events ++;



        /* 
        // Cut#5: Require the Hbb score > 0.8
        float maxhbbscore = Analysis::fatJets_[0].hbbScore;
        int maxHbbNo = -999;
        for (unsigned int ifatjet = 1; ifatjet < Analysis::fatJets_.size(); ifatjet++) {
            if (Analysis::fatJets_[ifatjet].hbbScore >= maxhbbscore) {
                maxHbbNo = ifatjet;
                maxhbbscore = Analysis::fatJets_[ifatjet].hbbScore;
            }
        }
        if (maxhbbscore < 0.8) { continue;}
        cut6_events ++;
        Cutflow::fillCutflow(Cutflow::Cuts::kHbbScore);
        

        // Cut#7: Require dRVBF > 3.5
        if (not (RooUtil::Calc::DeltaR(Analysis::VBFjets_[0], Analysis::VBFjets_[1]) > 3.5)) {continue; }
        Cutflow::fillCutflow(Cutflow::Cuts::kdRVBF);
        */

        Hist::fillSHatHistograms();
        Hist::fillMETHistograms();
        Hist::fillLeptonsKinematicHistograms();
        Hist::fillHbbFatJetKinematicHistograms();
        Hist::fillJetsKinematicHistograms();


        /*
        // Cut#7: Require mjj > 500, deltaEtajj > 3
        if (not ((Analysis::VBFjets_[0] + Analysis::VBFjets_[1]).M() > 500)) { continue;}
        if (not (TMath::Abs(Analysis::VBFjets_[0].Eta()-Analysis::VBFjets_[1].Eta()) > 3)) { continue;}
        Cutflow::fillCutflow(Cutflow::Cuts::kMjj);
        cut8_events ++;


        
        // Cut#9: Require dRLep<1.5
        LV lep1 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[0] : Analysis::leptons_[1];
        LV lep2 = (Analysis::leptons_[0]).Pt() > (Analysis::leptons_[1]).Pt() ? Analysis::leptons_[1] : Analysis::leptons_[0];
        if (not (RooUtil::Calc::DeltaR(lep1, lep2) <= 1.5)) { continue;}
        Cutflow::fillCutflow(Cutflow::Cuts::kDrlep);

        // Cut#10: Require massZH>250
        if (not ((Analysis::leptons_[0] + Analysis::leptons_[1] + Analysis::hbbFatJet_).M() >= 215)) { continue; } 
        Cutflow::fillCutflow(Cutflow::Cuts::kZhmass);


        // Cut#8: Require etaHbb<1.6
        /*
        if (not (Analysis::hbbFatJet_.Eta() <= 1.3)) { continue;}
        Cutflow::fillCutflow(Cutflow::Cuts::kEtahbb);

        */
        // Fill Lepton Histogram;

        // All cuts have passed
        // Now fill the histograms
        if (gen_level) { Hist::fillGenLevelHistograms(); }

 

    }

    // Write out the histograms
    Hist::writeHistograms(output_file, gen_level);

    // Write out the cutflow histogram
    Cutflow::writeCutflow(output_file);

    // Produce cutflow events number
    cout << "The first cut (vbs) produces " << cut1_events << " events." << endl;
    cout << "The second cut (hbb) produces " << cut2_events << " events." << endl;
    cout << "The third cut (wjj) produces " << cut3_events << " events." << endl;
    cout << "The fourth cut (2 OS leptons: ee or uu) produces " << cut4_events << "events." << endl;
    cout << "The fifth cut (at least one fatjet) produces " << cut5_events << " events." << endl;
    cout << "The sixth cut (hbb score) produces " << cut6_events << " events." << endl;
    cout << "The sixth cut (2 jets) produces " << cut7_events << " events." << endl;
    cout << "The sixth cut (mjj > 500, deltaEtajj) produces " << cut7_events << " events." << endl;

    // Close the file
    output_file->Close();


    return 0;
}
