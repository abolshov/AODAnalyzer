// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iomanip>

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
#include "TH1.h"


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
		static std::ofstream m_file;
		std::vector<reco::Candidate*> m_signal;
		std::vector<reco::Candidate*> m_matched_jets;
		double m_dr_thresh = 0.4;
		static int m_event;

		// std::vector<TLorentzVector> m_stable_particles;

		inline bool ValidSignal() const noexcept { return std::none_of(m_signal.begin(), m_signal.end(), [](reco::Candidate const* p) { return p == nullptr; }); }
		bool FindSignal(reco::GenParticleCollection const& genpart);
		int Match(GENPART target, std::vector<reco::GenJet> const& genjets);
		std::vector<std::pair<int, TLorentzVector>> FindStableParticles(reco::GenParticleCollection const& genpart);
		std::vector<TLorentzVector> ClusterAK4(std::vector<std::pair<int, TLorentzVector>> const& stable);
};

// static variables initialization
std::ofstream AODAnalyzer::m_file("debug.csv");
int AODAnalyzer::m_event = 1;

AODAnalyzer::AODAnalyzer(const edm::ParameterSet& iConfig)
    : m_genpart_token(consumes<reco::GenParticleCollection>(iConfig.getUntrackedParameter<edm::InputTag>("genparticles"))),
	  m_genjet_token(consumes<std::vector<reco::GenJet>>(iConfig.getUntrackedParameter<edm::InputTag>("genjets"))),
	  m_signal(static_cast<int>(GENPART::Count), nullptr),
	  m_matched_jets(static_cast<int>(GENJET::Count), nullptr),
	  m_dr_thresh(iConfig.getUntrackedParameter<double>("dr_thresh")) {}


AODAnalyzer::~AODAnalyzer() 
{
	for (auto& p_cand: m_signal)
	{
		delete p_cand;
	}

	
	for (auto& p_jet: m_matched_jets)
	{
		delete p_jet;
	}
}


void AODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) 
{
	std::cout << "analyze event " << m_event << "\n";
	++m_event;
	std::cout << std::setprecision(3);
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
		std::cout << "\tHeavy higgs Mass = " << m_signal[static_cast<int>(GENPART::H)]->mass() << "\n";

		int q1_match = Match(GENPART::q1, genjet);
		int q2_match = Match(GENPART::q2, genjet);

		int q1_idx = static_cast<int>(GENPART::q1);
		reco::Candidate const* ptr_q1 = m_signal[q1_idx];
		TLorentzVector q1(ptr_q1->px(), ptr_q1->py(), ptr_q1->pz(), ptr_q1->energy());

		int q2_idx = static_cast<int>(GENPART::q2);
		reco::Candidate const* ptr_q2 = m_signal[q2_idx];
		TLorentzVector q2(ptr_q2->px(), ptr_q2->py(), ptr_q2->pz(), ptr_q2->energy());

		int w_had_idx = static_cast<int>(GENPART::Whad);
		reco::Candidate const* ptr_Whad = m_signal[w_had_idx];
		TLorentzVector w_had(ptr_Whad->px(), ptr_Whad->py(), ptr_Whad->pz(), ptr_Whad->energy());

		int w_lep_idx = static_cast<int>(GENPART::Wlep);
		reco::Candidate const* ptr_Wlep = m_signal[w_lep_idx];
		TLorentzVector w_lep(ptr_Wlep->px(), ptr_Wlep->py(), ptr_Wlep->pz(), ptr_Wlep->energy());

		int lep_idx = static_cast<int>(GENPART::l);
		reco::Candidate const* ptr_lep = m_signal[lep_idx];
		TLorentzVector lep(ptr_lep->px(), ptr_lep->py(), ptr_lep->pz(), ptr_lep->energy());

		int nu_idx = static_cast<int>(GENPART::nu);
		reco::Candidate const* ptr_nu = m_signal[nu_idx];
		TLorentzVector nu(ptr_nu->px(), ptr_nu->py(), ptr_nu->pz(), ptr_nu->energy());

		int b_idx = static_cast<int>(GENPART::b);
		reco::Candidate const* ptr_b = m_signal[b_idx];
		TLorentzVector b(ptr_b->px(), ptr_b->py(), ptr_b->pz(), ptr_b->energy());

		int bbar_idx = static_cast<int>(GENPART::bbar);
		reco::Candidate const* ptr_bbar = m_signal[bbar_idx];
		TLorentzVector bbar(ptr_bbar->px(), ptr_bbar->py(), ptr_bbar->pz(), ptr_bbar->energy());
		
		if (q1_match != -1 && q2_match != -1)
		{
			double r1 = genjet[q1_match].pt()/q1.Pt();
			double r2 = genjet[q2_match].pt()/q2.Pt();
			bool j1_mismatch = (r1 < 0.8 || r1 > 1.2);
			bool j2_mismatch = (r2 < 0.8 || r2 > 1.2);

			if (j1_mismatch || j2_mismatch)
			{
				std::cout << "\tq1_Pt = " << q1.Pt() << "\n";
				std::cout << "\tq2_Pt = " << q2.Pt() << "\n";
				std::cout << "\tbq_Pt = " << b.Pt() << "\n";
				std::cout << "\tbbar_Pt = " << bbar.Pt() << "\n";
				std::cout << "\tdR(q1, q2) = " << q1.DeltaR(q2) << "\n";
				std::cout << "\tdR(q1, b) = " << q1.DeltaR(b) << "\n";
				std::cout << "\tdR(q1, bbar) = " << q1.DeltaR(bbar) << "\n";
				std::cout << "\tdR(q2, b) = " << q2.DeltaR(b) << "\n";
				std::cout << "\tdR(q2, bbar) = " << q2.DeltaR(bbar) << "\n";
				std::cout << "\tr1 = " << r1 << ", r2 = " << r2 << "\n";

				auto stable = FindStableParticles(genpart);

				std::cout << "\tStable particle info:\n";
				std::cout << "\tid\tdR(q1, p)\tdR(q2, p)\tPt\n";
				// TLorentzVector tot_stable;
				for (auto const& [id, p4]: stable)
				{
					// tot_stable += p4;
					std::cout << "\t" << id << "\t" << q1.DeltaR(p4) << "\t\t" << q2.DeltaR(p4) << "\t\t" << p4.Pt() << "\n";
				}

				// std::cout << "\tConsistency check: E(W) = " << w_had.E() << ", E(stable) = " << tot_stable.E() << "\n";

				auto jets = ClusterAK4(stable);
				std::cout << "\tnumber of stable partcles = " << stable.size() << "\n";
				std::cout << "\tcustom_size = " << jets.size() << ", genjets_size = " << genjet.size() << "\n";

				std::cout << "\tcustom jet pt: ";
				for (auto& cj: jets)
				{
					std::cout << cj.Pt() << " ";
				}
				std::cout << "\n";

				if (j1_mismatch)
				{
					std::cout << "\tJet 1:\n";
					for (auto& jet: jets)
					{
						double dR = jet.DeltaR(q1);
						if (dR < m_dr_thresh)
						{
							std::cout << "\tcustom: r = " << jet.Pt()/q1.Pt() << "\n";
						}
					}
					std::cout << "\tbuilt-in: r = " << r1 << "\n";
				}

				if (j2_mismatch)
				{
					std::cout << "\tJet 2:\n";
					for (auto& jet: jets)
					{
						double dR = jet.DeltaR(q2);
						if (dR < m_dr_thresh)
						{
							std::cout << "\tcustom: r = " << jet.Pt()/q2.Pt() << "\n";
						}
					}
					std::cout << "\tbuilt-in: r = " << r2 << "\n";
				}
			}
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

	std::cout << "\tgenjet coll size: " << genjets.size() << "\n";
	
	TLorentzVector part;
	part.SetPtEtaPhiM(m_signal[idx]->pt(), m_signal[idx]->eta(), m_signal[idx]->phi(), m_signal[idx]->mass());

	double best_dR = 0.5;
	int best_match = -1;
	for (int i = 0; i < static_cast<int>(genjets.size()); ++i)
	{
		TLorentzVector jet;
		jet.SetPtEtaPhiM(genjets[i].pt(), genjets[i].eta(), genjets[i].phi(), genjets[i].mass());
		double dR = jet.DeltaR(part);
		// std::cout << "\t jet #" << i << ": E = " << jet.E() << ", dR = " << dR << "\n";
		if (dR < best_dR && dR < m_dr_thresh)
		{
			best_dR = dR;
			best_match = i;
		}
	}

	return best_match;
}

std::vector<std::pair<int, TLorentzVector>> AODAnalyzer::FindStableParticles(reco::GenParticleCollection const& genpart)
{
	int sz = genpart.size();
	std::vector<int> forbidden{12, 14, 16};
	std::vector<std::pair<int, TLorentzVector>> res;
	for (int i = 0; i < sz; ++i)
	{
		int abs_id = std::abs(genpart[i].pdgId()); 
		// bool from_w = false;
		// reco::Candidate const* mmother = genpart[i].mother();
		// while (mmother)
		// {
		// 	if (mmother->pdgId() == m_signal[static_cast<int>(GENPART::Whad)]->pdgId())
		// 	{
		// 		from_w = true;
		// 		break;
		// 	}
		// 	mmother = mmother->mother();
		// }

		// if (!from_w)
		// {
		// 	continue;
		// }

		auto it = std::find(forbidden.begin(), forbidden.end(), abs_id);
		if (it != forbidden.end())
		{
			continue;
		}
		// if no daughters or not neutrino
		if (genpart[i].numberOfDaughters() == 0 && genpart[i].status() == 1)
		{
			TLorentzVector p(genpart[i].px(), genpart[i].py(), genpart[i].pz(), genpart[i].energy());
			int pdg_id = genpart[i].pdgId();
			res.push_back({pdg_id, p});
		}
	}
	return res;
}

std::vector<TLorentzVector> AODAnalyzer::ClusterAK4(std::vector<std::pair<int, TLorentzVector>> const& stable)
{
	// std::cout << "ClusterAK4\n";
	std::vector<TLorentzVector> res;
	auto d_ij = [](TLorentzVector const& pi, TLorentzVector const& pj)
	{
		double dR = pi.DeltaR(pj)*pi.DeltaR(pj);
		double min_inv_pt = std::min(1/(pi.Pt()*pi.Pt()), 1/(pj.Pt()*pj.Pt()));
		return min_inv_pt*(dR*dR)/(0.4*0.4);
	};

	int sz = stable.size();
	std::vector<TLorentzVector> p4;
	std::vector<int> id;
	id.reserve(sz);
	p4.reserve(sz);
	for (int i = 0; i < sz; ++i)
	{
		auto& [pid, p] = stable[i];
		p4.push_back(p);
		id.push_back(pid);
	}

	while (!p4.empty())
	{
		// if only one object is left:
		// if it is a particle ignore it
		// if it's a result of merging add to result
		int cur_size = p4.size();
		if (cur_size == 1)
		{
			if (id[0] != 0 && p4[0].Pt() > 5)
			{
				res.push_back(p4[0]);
			}
			break;
		}

		// here at least two objects are present in the collection
		double min_d_iB = 1/(p4[0].Pt()*p4[0].Pt());
		int min_iB = 0;

		double min_d_ij = d_ij(p4[0], p4[1]);
		int min_i = 0;
		int min_j = 1;

		for (int i = 0; i < cur_size; ++i)
		{
			double d_iB = 1/(p4[i].Pt()*p4[i].Pt());
			if (d_iB < min_d_iB)
			{
				min_d_iB = d_iB;
				min_iB = i;
			}

			for (int j = i + 1; j < cur_size; ++j)
			{
				if (i == j)
				{
					continue;
				}
				double cur_d_ij = d_ij(p4[i], p4[j]);
				if (cur_d_ij < min_d_ij)
				{
					min_d_ij = cur_d_ij;
					min_i = i;
					min_j = j;
				}
			}
		}

		// std::cout << "min d_iB = " << min_d_iB 
		// 		  << ":\n\tat = " << min_iB << ", pdgId = " << id[min_iB] << ", pt = " << p4[min_iB].Pt() << "\n";

		// std::cout << "min d_ij = " << min_d_ij
		// 		  << ":\n\tat (" << min_i << ", " << min_j << "), pdgId = (" << id[min_i] << ", " << id[min_j] << "), "
		// 		  << "pt = (" << p4[min_i].Pt() << ", " << p4[min_j].Pt() << ")\n";

		// if (std::min(min_d_iB, min_d_ij) > 0.6) break;

		if (min_d_ij < min_d_iB)
		{
			if (min_i == min_j)
			{
				std::cout << "\tCOMBINING SELF!!!\n";
			}

			// std::cout << "\tcombine (pdgId, idx) = (" << id[min_i] << ", " << min_i << "), pt = " << p4[min_i].Pt() 
			// 		  << " with (pdgId, idx) = (" << id[min_j] << ", " << min_j << "), pt = " << p4[min_j].Pt() << "\n";
			
			int keep = std::min(min_i, min_j);
			int remove = std::max(min_i, min_j);
			p4[keep] += p4[remove];
			auto itr_p4 = p4.begin() + remove;
			p4.erase(itr_p4);
			id[keep] = 0;
			auto itr_id = id.begin() + remove;
			id.erase(itr_id);

			// std::cout << "\tresult pt = " << p4[keep].Pt() << "\n";

			// combine i and j to min(i, j)
			// remove max(i, j)
		}
		else
		{
			if (id[min_iB] != 0 || p4[min_iB].Pt() < 5)
			{
				auto itr_p4 = p4.begin() + min_iB;
				auto itr_id = id.begin() + min_iB;
				p4.erase(itr_p4);
				id.erase(itr_id);
				continue;
			}
			// std::cout << "\tobjetct (pdgID, idx) = (" << id[min_iB] << ", " << min_iB << ") is a jet with energy " << p4[min_iB].E() << "\n";
			// std::cout << "\tformed jet with pt = " << p4[min_iB].Pt() << "\n";
			res.push_back(p4[min_iB]);
			auto itr_p4 = p4.begin() + min_iB;
			p4.erase(itr_p4);
			auto itr_id = id.begin() + min_iB;
			id.erase(itr_id);
			// std::cout << "================================================================================\n";
		}
	}
	return res;
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODAnalyzer);
