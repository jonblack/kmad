#include "seq_data.h"

#include <boost/regex.hpp>
#include <algorithm>


seq_data::SequenceData seq_data::process_fasta_data(
    const fasta::FastaData& fasta_data, 
    const f_config::FeatureSettingsMap& f_set, bool gapped, 
    const std::map<std::string, double>& probabilities) {
  seq_data::SequenceData s;
  // probabilities are a concatenation of probs from fasta data and probs
  // from the config file (if there are key repeats - thereshoudl be none -
  // then the one from the config file is kept)
  s.probabilities = probabilities;
  s.probabilities.insert(fasta_data.probabilities.begin(),
                         fasta_data.probabilities.end());
  s.sequences = fasta_data.sequences;
  if (!gapped) {
    s.sequences = remove_gaps(s.sequences);
  }
  for (auto feat_it = f_set.begin(); feat_it != f_set.end(); ++feat_it) {
    assign_feature_by_pattern(s.sequences, feat_it->second.pattern,
        feat_it->first);
    for (auto& seq : feat_it->second.positions) {
      if ((signed)s.sequences.size() > seq.seq_no && seq.seq_no >= 0) {
        for (auto& pos : seq.positions) {
          if ((signed)s.sequences[seq.seq_no].residues.size() > pos && pos >= 0) {
            s.sequences[seq.seq_no].residues[pos].features.push_back(
                feat_it->first);
          }
          else {
            std::cout << "Warning: feature positions should be in range: 1 - "
                          << "sequence length, feature " << feat_it->first
                          << " cannot be annotated at position " << pos
                          << " in sequence " << seq.seq_no << std::endl;
          }
        }
      }
      else {
        std::cout << "Warning: sequence numbers should be in range: 1 - "
                      << "number of sequences (" << s.sequences.size()
                      << "), feature " << feat_it->first
                      << " cannot be annotated in sequence "
                      << seq.seq_no << std::endl;
      }
    }
  }
  s.feature_list = make_feature_list(s.sequences);
  return s;
}


FeatureNamesList seq_data::make_feature_list(
    const fasta::SequenceList& sequences) {
  FeatureNamesList feature_list = {"ptm_phosph0", "ptm_phosph1",
                                   "ptm_phosph2", "ptm_phosph3",
                                   "ptm_phosphP", "ptm_acet0",
                                   "ptm_acet1", "ptm_acet2",
                                   "ptm_acet3", "ptm_Nglyc0",
                                   "ptm_Nglyc1", "ptm_Nglyc2",
                                   "ptm_Nglyc3", "ptm_amid0",
                                   "ptm_amid1", "ptm_amid2",
                                   "ptm_amid3", "ptm_hydroxy0",
                                   "ptm_hydroxy1", "ptm_hydroxy2",
                                   "ptm_hydroxy3", "ptm_methyl0",
                                   "ptm_methyl1", "ptm_methyl2",
                                   "ptm_methyl3", "ptm_Oglyc0",
                                   "ptm_Oglyc1", "ptm_Oglyc2",
                                   "ptm_Oglyc3", "ptm_cys_bridge0",
                                   "strct_a_helix", "strct_turn",
                                   "strct_b_ladder", "strct_b_bridge",
                                   "strct_310_helix", "strct_pi_helix",
                                   "strct_b_ladder"}; 
  for (auto& seq : sequences) {
    for (auto& res : seq.residues) {
      for (auto& feat_name : res.features) {
        if (std::find(feature_list.begin(), feature_list.end(), feat_name) 
            == feature_list.end()) {
          feature_list.push_back(feat_name);
        }
      }
    }
  }
  return feature_list;
}


fasta::SequenceList seq_data::remove_gaps(
    const fasta::SequenceList& sequences) {
  fasta::SequenceList s = sequences;
  for (auto& seq : s) {
    seq.residues.clear();
  }
  for (size_t i = 0; i < sequences.size(); ++i) {
    for (size_t j = 0; j < sequences[i].residues.size(); ++j) {
      if (sequences[i].residues[j].codon[0] != '-') {
        s[i].residues.push_back(sequences[i].residues[j]);
      }
    }
  }
  return s;
}

void seq_data::assign_feature_by_pattern(fasta::SequenceList& sequences,
                                         const std::string& pattern,
                                         const std::string& feat_name)
{
  if (pattern.size() > 0) {
    boost::regex re(pattern);
    for (size_t i = 0; i < sequences.size(); ++i) {
      std::string seq = fasta::make_string(sequences[i]);
      std::string seq_nogaps = seq;
      seq_nogaps.erase(std::remove(seq_nogaps.begin(), seq_nogaps.end(), '-'),
          seq_nogaps.end());
      for(auto it = boost::sregex_iterator(seq_nogaps.begin(), seq_nogaps.end(),
            re);
              it != boost::sregex_iterator();
                   ++it)
      {
        int match_start = find_real_pos(seq, it->position());
        int match_end = find_real_pos(seq, match_start + it->str().size());
        for (int j = match_start; j < match_end; ++j) {
          if (sequences[i].residues[j].codon[0] != '-') {
            sequences[i].residues[j].features.push_back(feat_name);
          }
        }
      }
    }
 }
}


int seq_data::find_real_pos(const std::string& sequence, int position) {
  int pos = 0;
  size_t i = 0;
  while (i < sequence.size() && pos < position) {
    if (sequence[i] != '-') {
      ++pos;
    }
    ++i;
  }
  return i;
}


bool seq_data::compare_alignments(const std::vector<fasta::SequenceList>& al1,
    const std::vector<fasta::SequenceList>& al2) {
  bool result = true;
  if (al1.size() != al2.size() || al1[0].size() != al2[0].size()) {
    result = false;
  }
  std::vector<std::string> al1_strings;
  std::vector<std::string> al2_strings;
  for (auto& item : al1) {
    for (auto& seq : item) {
      al1_strings.push_back(fasta::make_string(seq));
    }
  }
  for (auto& item : al2) {
    for (auto& seq : item) {
      al2_strings.push_back(fasta::make_string(seq));
    }
  }
  if (al1_strings != al2_strings) {
    result = false;
  }
  return result;
}
