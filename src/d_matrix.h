#ifndef DMATRIX_H
#define DMATRIX_H

#include "feature_scores.h"
#include "profile.h"
#include "seq_data.h"


typedef std::vector<profile::ProfileMap> ProfileList;
typedef std::vector<FeatureScores> FeatureProfileList;
namespace d_matrix {

  struct DistanceMatrix {
    ProfileList profiles;
    FeatureProfileList feature_profiles;
    std::vector<std::vector<double> > values;
  };

  DistanceMatrix make_distance_matrix(
      const f_config::FeatureSettingsMap& f_set,
      const seq_data::SequenceData& sequence_data, double gap_open_pen,
      double end_pen, double gap_ext_pen, double domain_modifier,
      double motif_modifier, double ptm_modifier, double strct_modifier,
      int codon_length, const bool no_feat, const std::string& sbst_mat);


  void join_nodes(DistanceMatrix& dist_matrix, int node_index1,
                  int node_index2);


  void make_profile_lists(
      ProfileList& profiles,
      FeatureProfileList& f_profiles, double domain_modifier,
      double motif_modifier, double ptm_modifier, double strct_modifier,
      const std::string& sbst_mat, const bool no_feat,
      const f_config::FeatureSettingsMap& f_set,
      const seq_data::SequenceData& sequence_data);

  DistanceMatrix update_matrix(const DistanceMatrix& dist_matrix,
      int node1, int node2, const profile::ProfileMap& new_profile);


  std::pair<int, int> find_closest_nodes(const DistanceMatrix& dist_matrix);
  DistanceMatrix update_matrix(DistanceMatrix& dist_matrix, int node1,
      int node2, const profile::ProfileMap& profile);
}

#endif /* DMATRIX_H */
