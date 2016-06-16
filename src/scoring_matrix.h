#ifndef SCORINGMATRIX_H
#define SCORINGMATRIX_H


#include "fasta.h"
#include "types.h"

namespace t = types;

typedef std::vector<double> ScoringMatrixRow;
typedef std::vector<ScoringMatrixRow> SingleScoringMatrix;
typedef std::vector<int> ValueCoords;

class ScoringMatrix {
public:
  ///
  /// Constructor; creates scoring matrices for the two sequences
  /// (sequence and a profile)
  /// @param s1size length of the pseudo sequence (=profile)
  /// @param s2size lengtn of the sequence that will be aligned to the profile
  /// @param pen gap opening penalty
  /// @param extensionPenalty gap extension penalty
  /// @param endPenalty penalty for gaps at the beginning and the end
  ///
  ScoringMatrix(int profile_length, int sequence_length,
	  t::SettingsMap& aln_params);
  ///
  /// Fills in the scoring matrices m_matrix_v, m_matrix_g, m_matrix_h
  ///
  void calculate_scores(const fasta::Sequence& sequence,
                        const profile::ProfileMap& profile,
                        const FeatureScores& f_profile, int codon_length);
  ///
  /// traces back the alignment path in the scoring matrices
  ///
  fasta::SequenceList backtrace_alignment_path(
      const fasta::Sequence& sequence,
      const profile::ProfileMap& profile,
      const FeatureScores& f_profile,
      int codon_length);
  SingleScoringMatrix get_V_matrix();
private:
  ///
  /// finds the best score either in the last column or in the last row of the
  /// V matrix (takes the end gap penaltie into account)
  ///
  ValueCoords find_best_score();
  int m_i_length;
  int m_j_length;
  double m_gap_opening;
  double m_gap_extension;
  double m_end_pen;
  bool m_no_feat;
  SingleScoringMatrix m_matrix_v, m_matrix_g, m_matrix_h;
};

#endif /* SCORINGMATRIX_H */
