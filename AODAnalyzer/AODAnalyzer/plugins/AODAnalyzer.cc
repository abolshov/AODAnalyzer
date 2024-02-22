// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "TLorentzVector.h"
static constexpr double PI = 3.14159;

class AODAnalyzer : public edm::one::EDAnalyzer<> 
{
	public:
		explicit AODAnalyzer(const edm::ParameterSet&);
		~AODAnalyzer() override;
		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		void beginJob() override;
		void analyze(const edm::Event&, const edm::EventSetup&) override;
		void endJob() override;

		edm::EDGetTokenT<reco::GenParticleCollection> m_genpart_token; 
		edm::EDGetTokenT<std::vector<reco::GenJet>> m_genjet_token;

		enum class GENPART { H, h_bb, h_WW, b, bbar, Whad, Wlep, q1, q2, l, nu, Count };
		enum class GENJET { b, bbar, q1, q2, Count }; // either make enum class or do something else
		std::vector<reco::Candidate*> m_signal;
		std::vector<reco::Candidate*> m_matched_jets;

		inline bool ValidSignal() const noexcept { return std::none_of(m_signal.begin(), m_signal.end(), [](reco::Candidate const* p) { return p == nullptr; }); }
		bool FindSignal(reco::GenParticleCollection const& genpart);
		int Match(GENPART target, std::vector<reco::GenJet> const& genjets);
};

AODAnalyzer::AODAnalyzer(const edm::ParameterSet& iConfig)
    : m_genpart_token(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"))),
	  m_genjet_token(consumes<std::vector<reco::GenJet>>(iConfig.getUntrackedParameter<edm::InputTag>("genjets"))),
	  m_signal(static_cast<int>(GENPART::Count), nullptr),
	  m_matched_jets(static_cast<int>(GENJET::Count), nullptr) {}


AODAnalyzer::~AODAnalyzer() 
{
	if (!m_signal.empty())
	{
		for (auto& p_cand: m_signal)
		{
			if (p_cand) delete p_cand;
		}
	}

	if (!m_matched_jets.empty())
	{
		for (auto& p_jet: m_matched_jets)
		{
			if (p_jet) delete p_jet;
		}
	}
}


void AODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
	std::cout << "analyze\n";
	edm::Handle<reco::GenParticleCollection> genpart_handle;
	iEvent.getByToken(m_genpart_token, genpart_handle);
  	edm::Handle<std::vector<reco::GenJet>> genjet_handle;
	iEvent.getByToken(m_genjet_token, genjet_handle);

	std::vector<reco::GenParticle> const& genpart = *genpart_handle;
	std::vector<reco::GenJet> const& genjet = *genjet_handle;

	// "dereference" handle to get const& to underlying type
	if (FindSignal(genpart))
	{
		std::cout << "\tFound signal in collection of size " << genpart.size() << "\n";
		for (auto const pcand: m_signal)
		{
			if (pcand)
			{
				std::cout << "\t" << pcand->pdgId() << "\n";
			}
			else
			{
				std::cout << "\tnullptr\n";
			}
		}

		int q1_match = Match(GENPART::q1, genjet);
		if (q1_match != -1)
		{
			// m_matched_jets[static_cast<int>(GENJET::q1)] = genjet[q1_match].clone();
			std::cout << "\tMatched jet with energy "<< genjet[q1_match].energy() << " to light quark 1 with energy " << m_signal[static_cast<int>(GENPART::q1)]->energy() << "\n";
		}
	}
}

void AODAnalyzer::beginJob() {}
void AODAnalyzer::endJob() {}
void AODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {}

bool AODAnalyzer::FindSignal(reco::GenParticleCollection const& genpart)
{
	auto heavy_higgs = std::find_if(genpart.rbegin(), genpart.rend(), [](reco::Candidate const& cand) { return cand.pdgId() == 35; });
	if (heavy_higgs == genpart.rend()) return false;
	m_signal[static_cast<int>(GENPART::H)] = heavy_higgs->clone();

	if (heavy_higgs->numberOfDaughters() != 2) return false;

	// heavy higgs decays immediately only to pair of standard higgses
	if (heavy_higgs->daughter(0)->pdgId() != 25 || heavy_higgs->daughter(1)->pdgId() != 25) return false;

	reco::Candidate const* higgs_1 = heavy_higgs->daughter(0); 
	reco::Candidate const* higgs_2 = heavy_higgs->daughter(1); 

	if (higgs_1->numberOfDaughters() != 2 || higgs_2->numberOfDaughters() != 2) return false;

	// deal with daughters of first higgs and first higgs
	reco::Candidate const* daughter_1 = higgs_1->daughter(0); 
	reco::Candidate const* daughter_2 = higgs_1->daughter(1);

	reco::Candidate const* w_plus = nullptr;
	reco::Candidate const* w_minus = nullptr;

	if (std::abs(daughter_1->pdgId()) == std::abs(daughter_2->pdgId()))
	{
		if (std::abs(daughter_1->pdgId()) == 5)
		{
			m_signal[static_cast<int>(GENPART::bbar)] = (daughter_1->pdgId() == -5) ? daughter_1->clone() : daughter_2->clone();
			m_signal[static_cast<int>(GENPART::b)] = (daughter_1->pdgId() == 5) ? daughter_1->clone() : daughter_2->clone();
			m_signal[static_cast<int>(GENPART::h_bb)] = higgs_1->clone();
		}
		else if (std::abs(daughter_1->pdgId()) == 24)
		{
			w_plus = (daughter_1->pdgId() == 24) ? daughter_1 : daughter_2;
			w_minus = (daughter_1->pdgId() == 24) ? daughter_2 : daughter_1;
			m_signal[static_cast<int>(GENPART::h_WW)] = higgs_1->clone();
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	} 

	// deal with daughters of second higgs and second higgs
	daughter_1 = higgs_2->daughter(0); 
	daughter_2 = higgs_2->daughter(1);

	if (std::abs(daughter_1->pdgId()) == std::abs(daughter_2->pdgId()))
	{
		if (std::abs(daughter_1->pdgId()) == 5)
		{
			m_signal[static_cast<int>(GENPART::bbar)] = (daughter_1->pdgId() == -5) ? daughter_1->clone() : daughter_2->clone();
			m_signal[static_cast<int>(GENPART::b)] = (daughter_1->pdgId() == 5) ? daughter_1->clone() : daughter_2->clone();
			m_signal[static_cast<int>(GENPART::h_bb)] = higgs_2->clone();
		}
		else if (std::abs(daughter_1->pdgId()) == 24)
		{
			w_plus = (daughter_1->pdgId() == 24) ? daughter_1 : daughter_2;
			w_minus = (daughter_1->pdgId() == 24) ? daughter_2 : daughter_1;
			m_signal[static_cast<int>(GENPART::h_WW)] = higgs_2->clone();
		}
		else
		{
			return false;
		}
	}
	else
	{
		return false;
	}

	// find light quarks and leptons
	if (w_plus->numberOfDaughters() != 2 || w_minus->numberOfDaughters() != 2) return false;

	//deal with daughters of W plus
	daughter_1 = w_plus->daughter(0);
	daughter_2 = w_plus->daughter(1);

	// W plus decays hadronically if pdgIds of daughters is less than b quark
	bool hadronic = (std::abs(daughter_1->pdgId()) < 5 && std::abs(daughter_2->pdgId()) < 5);
	bool leptonic = (std::abs(daughter_1->pdgId()) <= 14 && std::abs(daughter_1->pdgId()) >= 11) && (std::abs(daughter_2->pdgId()) <= 14 && std::abs(daughter_2->pdgId()) >= 11);
	if (hadronic)
	{
		m_signal[static_cast<int>(GENPART::Whad)] = w_plus->clone();
		m_signal[static_cast<int>(GENPART::q1)] = w_plus->daughter(0)->clone();
		m_signal[static_cast<int>(GENPART::q2)] = w_plus->daughter(1)->clone();
	}
	else if (leptonic)
	{
		m_signal[static_cast<int>(GENPART::Wlep)] = w_plus->clone();
		m_signal[static_cast<int>(GENPART::l)] = w_plus->daughter(0)->clone();
		m_signal[static_cast<int>(GENPART::nu)] = w_plus->daughter(1)->clone();
	}
	else
	{
		return false;
	}

	//deal with daughters of W minus
	daughter_1 = w_minus->daughter(0);
	daughter_2 = w_minus->daughter(1);

	// W minus decays hadronically if pdgIds of daughters is less than b quark
	hadronic = (std::abs(daughter_1->pdgId()) < 5 && std::abs(daughter_2->pdgId()) < 5);
	leptonic = (std::abs(daughter_1->pdgId()) <= 14 && std::abs(daughter_1->pdgId()) >= 11) && (std::abs(daughter_2->pdgId()) <= 14 && std::abs(daughter_2->pdgId()) >= 11);
	if (hadronic)
	{
		m_signal[static_cast<int>(GENPART::Whad)] = w_minus->clone();
		m_signal[static_cast<int>(GENPART::q1)] = w_minus->daughter(0)->clone();
		m_signal[static_cast<int>(GENPART::q2)] = w_minus->daughter(1)->clone();
	}
	else if (leptonic)
	{
		m_signal[static_cast<int>(GENPART::Wlep)] = w_minus->clone();
		m_signal[static_cast<int>(GENPART::l)] = w_minus->daughter(0)->clone();
		m_signal[static_cast<int>(GENPART::nu)] = w_minus->daughter(1)->clone();
	}
	else
	{
		return false;
	}
	
	if (!ValidSignal()) return false;

	return true;
}

int AODAnalyzer::Match(GENPART target, std::vector<reco::GenJet> const& genjets)
{
	int idx = static_cast<int>(target);
	double t_phi = m_signal[idx]->phi();
	double t_eta = m_signal[idx]->eta();
	
	TLorentzVector part;

	double best_dR = 0.5;
	int best_match = -1;
	for (int i = 0; i < static_cast<int>(genjets.size()); ++i)
	{
		double j_phi = genjets[i].phi();
		double j_eta = genjets[i].eta();
		double deta = t_eta - j_eta;
		double dphi = (std::abs(t_phi - j_phi) <= PI) ? std::abs(t_phi - j_phi) : std::abs(t_phi - j_phi) - PI;
		double dR = std::sqrt(deta*deta + dphi*dphi);
		if (dR < best_dR && dR < 0.4)
		{
			best_dR = dR;
			best_match = i;
		}
	}

	return best_match;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODAnalyzer);
