#include "d_matrix.h"
#include "msa.h"

d_matrix::DistanceMatrix d_matrix::make_distance_matrix(
    const f_config::FeatureSettingsMap& f_set,
    const seq_data::SequenceData& sequence_data, double gap_open_pen,
    double end_pen, double gap_ext_pen, double domain_modifier,
    double motif_modifier, double ptm_modifier, double strct_modifier,
    int codon_length, const bool no_feat, const std::string& sbst_mat)
{
  d_matrix::DistanceMatrix dist_matrix;
  // identity of the 1st one to itself
  //
  size_t seq_number = sequence_data.sequences.size();

  // assign a NxN matrxi of zeroes to dist_matrix.values
  std::vector<double> row(seq_number, 0);
  dist_matrix.values.assign(seq_number, row);

  make_profile_lists(dist_matrix.profiles, dist_matrix.feature_profiles,
                     domain_modifier, motif_modifier, ptm_modifier,
                     strct_modifier, sbst_mat, no_feat, f_set, sequence_data);

  fasta::Sequence aligned_seq_uppercase; 
  //pairwise alignment with lowercase characters where chars were removed
  fasta::Sequence aligned_seq_with_lower; 
  bool gapped = true;
  for (size_t i = 0; i < seq_number; ++i) {
    for (size_t j = 0; j < seq_number; ++j) {
    // aligned_sequence: vector
    // first element is a dummy polyA sequence to indicate where are the gaps
    // in the profile,
    // second one is with lowercase where the gaps
    // were cut out
    if (i !=j) {
      fasta::SequenceList aligned_sequence = msa::align_pairwise(
          sequence_data.sequences[i], dist_matrix.profiles[j],
          dist_matrix.feature_profiles[j], gap_open_pen, end_pen, 
          gap_ext_pen, codon_length, gapped, no_feat);


      double identity = msa::calc_identity(aligned_sequence[0],
                                           aligned_sequence[1],
                                           sequence_data.sequences[j]);
      dist_matrix.values[i][j] = identity;
      dist_matrix.values[j][i] = identity;
    }
    }
  }

  return dist_matrix;
}

void d_matrix::make_profile_lists(ProfileList& profiles,
    FeatureProfileList& f_profiles, double domain_modifier,
    double motif_modifier, double ptm_modifier, double strct_modifier,
    const std::string& sbst_mat, const bool no_feat,
    const f_config::FeatureSettingsMap& f_set,
    const seq_data::SequenceData& sequence_data)
{
  bool fade_out = false;
  for (auto& seq : sequence_data.sequences) {
      fasta::SequenceList tmp_seq_list = {seq};
      FeatureScores f_profile(sequence_data.feature_list, domain_modifier,
                              ptm_modifier, motif_modifier, strct_modifier,
                              sequence_data.probabilities);
      profile::ProfileMap profile = profile::create_score_profile(
          tmp_seq_list, sbst_mat);
      std::vector<double> fake_identities = {1};
      if (!no_feat) {
        f_profile.update_scores(tmp_seq_list, f_set, fake_identities,
            fade_out=false);
      }
      profiles.push_back(profile);
      f_profiles.push_back(f_profile);
  }
}
