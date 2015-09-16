#include "profile.h"

#include <boost/range/numeric.hpp>

#include <algorithm>

namespace {
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
  static const profile::SimilarityScoresMap BLOSUM = {
    {'A', { 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0}},
    {'R', {-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3}},
    {'N', {-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3}},
    {'D', {-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3}},
    {'C', { 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1}},
    {'Q', {-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2}},
    {'E', {-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2}},
    {'G', { 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3}},
    {'H', {-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3}},
    {'I', {-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3}},
    {'L', {-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1}},
    {'K', {-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2}},
    {'M', {-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1}},
    {'F', {-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1}},
    {'P', {-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2}},
    {'S', { 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2}},
    {'T', { 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0}},
    {'W', {-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3}},
    {'Y', {-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1}},
    {'V', { 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4}}};
  static const profile::SimilarityScoresMap DISORDER = {
    {'A', { 9,  1,  2,  2,  1,  3,  2,  2,  1,  1,  1,  2,  0,  0,  3,  3,  4, -1,  0,  3}},
    {'R', { 1, 10,  1,  0,  1,  3,  1,  2,  3,  0,  0,  5, -1, -2,  0,  1,  1,  2, -1,  0}},
    {'N', { 2,  1, 11,  4,  1,  3,  2,  2,  3,  0, -1,  3, -1, -1,  1,  4,  3, -2,  1,  0}},
    {'D', { 2,  0,  4, 10, -2,  2,  5,  2,  2, -1, -2,  1, -2, -2,  1,  2,  2, -3, -1,  0}},
    {'C', { 1,  1,  1, -2, 17,  0, -3,  1,  0,  0,  0, -1, -2,  1, -1,  2,  1,  3,  2,  1}},
    {'Q', { 3,  3,  3,  2,  0, 11,  4,  0,  4,  0,  1,  3,  0, -1,  2,  2,  2,  1,  0,  1}},
    {'E', { 2,  1,  2,  5, -3,  4,  9,  1,  1, -1, -1,  3, -1, -2,  0,  1,  1, -2, -2,  0}},
    {'G', { 2,  2,  2,  2,  1,  0,  1, 10,  0, -2, -2,  0, -2, -2,  0,  2,  1,  0, -2,  0}},
    {'H', { 1,  3,  3,  2,  0,  4,  1,  0, 13,  0,  0,  1, -1,  2,  1,  1,  1,  1,  4,  0}},
    {'I', { 1,  0,  0, -1,  0,  0, -1, -2,  0, 12,  5,  0,  4,  4,  0,  0,  2,  1,  2,  7}},
    {'L', { 1,  0, -1, -2,  0,  1, -1, -2,  0,  5, 10, -1,  4,  5,  1,  0,  1,  2,  2,  4}},
    {'K', { 2,  5,  3,  1, -1,  3,  3,  0,  1,  0, -1, 10, -1, -2,  0,  1,  2, -2, -1,  0}},
    {'M', { 0, -1, -1, -2, -2,  0, -1, -2, -1,  4,  4, -1, 13,  2, -2, -1,  1,  0,  0,  3}},
    {'F', { 0, -2, -1, -2,  1, -1, -2, -2,  2,  4,  5, -2,  2, 13, -1,  0,  0,  6,  8,  3}},
    {'P', { 3,  0,  1,  1, -1,  2,  0,  0,  1,  0,  1,  0, -2, -1, 11,  2,  2, -1, -2,  1}},
    {'S', { 3,  1,  4,  2,  2,  2,  1,  2,  1,  0,  0,  1, -1,  0,  2,  9,  4, -1,  0,  1}},
    {'T', { 4,  1,  3,  2,  1,  2,  1,  1,  1,  2,  1,  2,  1,  0,  2,  4, 10, -2,  0,  3}},
    {'W', {-1,  2, -2, -3,  3,  1, -2,  0,  1,  1,  2, -2,  0,  6, -1, -1, -2, 18,  6,  0}},
    {'Y', { 0, -1,  1, -1,  2,  0, -2, -2,  4,  2,  2, -1,  0,  8, -2,  0,  0,  6, 14,  1}}, 
    {'V', { 3,  0,  0,  0,  1,  1,  0,  0,  0,  7,  4,  0,  3,  3,  1,  1,  3,  0,  1, 11}}};
}

profile::ProfileMap profile::create_score_profile(
    const fasta::SequenceList& sequences, const std::string& sbst_mat) {
  profile::ProfileMap p = create_profile(sequences);

  // Convert the profile occurrences to probabilities.
  for (auto& occ: p) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      ++i;
    }
  }
  const profile::SimilarityScoresMap* sim_scores;
  if (sbst_mat == "BLOSUM") {
    sim_scores = &BLOSUM;
  } else {
    sim_scores = &DISORDER;
  }
  for (unsigned i = 0; i < p['A'].size(); ++i) {
    std::vector<double> score_column(ALPHABET.size(), 0); 
    for (auto &prob: p) {
      std::vector<double> sbst_column = sim_scores->at(prob.first);
      for (size_t k = 0; k < sbst_column.size(); ++k) {
        score_column[k] += sbst_column[k] * prob.second[i];
      }
    }
    // replace the old column (with probs) with the new column (with scores)
    for (size_t j = 0; j < ALPHABET.size(); ++j) {
      p[ALPHABET[j]][i] = score_column[j];
    }
  }
  p['B'] = std::vector<double>(sequences[0].residues.size(), 0);
  p['Z'] = std::vector<double>(sequences[0].residues.size(), 0);
  p['X'] = std::vector<double>(sequences[0].residues.size(), 0);
  for (size_t i = 0; i < sequences[0].residues.size(); ++i) {
    for (auto& aa : ALPHABET ) {
      p['X'][i] += p[aa][i] / 20;
    }
    p['B'][i] = (p['D'][i] + p['N'][i]) / 2;
    p['Z'][i] = (p['Q'][i] + p['E'][i]) / 2;
  }
  return p;
}


profile::ProfileMap profile::create_profile(
    const fasta::SequenceList& sequences) {
  profile::ProfileMap p;

  // Initialise profile map with all letters  except for B, X, Z to 
  // a vector of the correct
  // size, with all values set to 0.
  for (size_t i = 0; i < ALPHABET.size(); ++i) {
    p[ALPHABET[i]] = std::vector<double>(sequences[0].residues.size(), 0);
  }

  for (size_t i = 0; i < sequences[0].residues.size(); ++i) {
    for (size_t j = 0; j < sequences.size(); ++j) {
      char amino_acid = sequences[j].residues[i].codon[0];

      // if (amino_acid == '-') {
      //   continue;
      // }

      if (amino_acid == 'B') {
        assert(p.find('D') != p.end() && p.find('N') != p.end());
        p['D'][i] += 0.5;
        p['N'][i] += 0.5;
      } else if (amino_acid == 'Z') {
        assert(p.find('E') != p.end() && p.find('Q') != p.end());
        p['E'][i] += 0.5;
        p['Q'][i] += 0.5;
      } else if (amino_acid == 'X' || p.find(amino_acid) == p.end()) {
        for (auto& kv: p) {
          kv.second[i] += 0.05;
        }
      } else if (amino_acid != '-'){
          p[amino_acid][i] += 1.0;
      }
    }
  }
  return p;
}


double profile::get_score(const profile::ProfileMap& p, int position,
                          char aa) {
  double result = 0;
  if (p.find(aa) != p.end()) {
    result = p.at(aa)[position];
  } else if ( aa == 'B') {
    result = p.at('D')[position] * 0.5 + p.at('N')[position] * 0.5;
  } else if ( aa == 'Z') {
    result = p.at('E')[position] * 0.5 + p.at('Q')[position] * 0.5;
  } else {
    for (auto& kv: p) {
      result += kv.second[position] * 0.05;
    }
  }
  return result;
}

double profile::profile_score(const profile::ProfileMap& profile1,
                              const profile::ProfileMap& profile2,
                              int pos1, int pos2, 
                              const profile::SimilarityScoresMap* sim_scores)
{
  double result = 0;
  int aa2_index = 0;
  double pos1_sum = profile::calc_sum(profile1, pos1);
  double pos2_sum = profile::calc_sum(profile2, pos2);
  for (auto i = profile1.begin(); i != profile1.end(); ++i) {
    char aa1 = i->first;
    for (auto j = profile2.begin(); j != profile2.end(); ++j) {
      result += (i->second[pos1] / pos1_sum) * (j->second[pos2] / pos2_sum) * sim_scores->at(aa1)[aa2_index];
      aa2_index++; 
    }
  }
  return result;
}

double profile::calc_sum(const profile::ProfileMap& profile, int pos) {
  double sum = 0;
  for (auto i = profile.begin(); i != profile.end(); ++i) {
    sum += i->second[pos];
  }
  return sum;
}

// based on the two dummy sequences representing profile
// alignment create a new profile
profile::ProfileMap profile::make_new_profile(
    const fasta::SequenceList& sequences,
    const profile::ProfileMap& profile1, const profile::ProfileMap& profile2)
{
  profile::ProfileMap new_profile;
  for (unsigned pos = 0; pos < sequences.size(); ++pos) {
    if (sequences[0].residues[pos].codon[0] != '-'
        && sequences[1].residues[pos].codon[0] != '-') {
      for (auto it = profile1.begin(); it != profile1.end(); ++it) {
        double new_item = (it->second[pos] + profile2.at(it->first)[pos]) / 2;
        new_profile.at(it->first).push_back(new_item);
      }
    } else if (sequences[0].residues[pos].codon[0] == '-') {
      for (auto it = profile1.begin(); it != profile1.end(); ++it) {
        double new_item = it->second[pos];
        new_profile.at(it->first).push_back(new_item);
      }
    } else {
      for (auto it = profile2.begin(); it != profile2.end(); ++it) {
        double new_item = it->second[pos];
        new_profile.at(it->first).push_back(new_item);
      }
    }
  }
  return new_profile;
}

profile::ProfileMap profile::make_scores(profile::ProfileMap& p,
    const fasta::SequenceList& sequences,
    const std::string& sbst_mat) 
{
  for (auto& occ: p) {
    size_t i = 0;
    for (auto& v: occ.second) {
      occ.second[i] = v / sequences.size();
      ++i;
    }
  }
  const profile::SimilarityScoresMap* sim_scores;
  if (sbst_mat == "BLOSUM") {
    sim_scores = &BLOSUM;
  } else {
    sim_scores = &DISORDER;
  }
  for (unsigned i = 0; i < p['A'].size(); ++i) {
    std::vector<double> score_column(ALPHABET.size(), 0); 
    for (auto &prob: p) {
      std::vector<double> sbst_column = sim_scores->at(prob.first);
      for (size_t k = 0; k < sbst_column.size(); ++k) {
        score_column[k] += sbst_column[k] * prob.second[i];
      }
    }
    // replace the old column (with probs) with the new column (with scores)
    for (size_t j = 0; j < ALPHABET.size(); ++j) {
      p[ALPHABET[j]][i] = score_column[j];
    }
  }
  p.at('B') = std::vector<double>(sequences[0].residues.size(), 0);
  p.at('Z') = std::vector<double>(sequences[0].residues.size(), 0);
  p.at('X') = std::vector<double>(sequences[0].residues.size(), 0);
  for (size_t i = 0; i < sequences[0].residues.size(); ++i) {
    for (auto& aa : ALPHABET ) {
      p.at('X')[i] += p.at(aa)[i] / 20;
    }
    p.at('B')[i] = (p.at('D')[i] + p.at('N')[i]) / 2;
    p.at('Z')[i] = (p.at('Q')[i] + p.at('E')[i]) / 2;
  }
  return p;
}
