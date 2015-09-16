#define BOOST_TEST_DYN_LINK

#include "src/d_matrix.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>


BOOST_AUTO_TEST_SUITE(test_d_matrix)

BOOST_AUTO_TEST_CASE(test_make_distance_matrix)
{
  int codon_length = 1;
  fasta::Sequence s1 = fasta::make_sequence("d", "ASTSTAAAAA", codon_length);
  fasta::Sequence s2 = fasta::make_sequence("d", "ASSSTAAAAA", codon_length);
  fasta::Sequence s3 = fasta::make_sequence("d", "ATTSTAAAAA", codon_length);
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
  std::map<std::string, double> probabilities;
  f_config::FeatureSettingsMap f_set;
  fasta::SequenceList sequences = {s1, s2, s3};
  int domain_modifier = 4;
  int motif_modifier = 3;
  int ptm_modifier = 10;
  int strct_modifier = 0;
  double gap_open_pen = -5;
  double gap_ext_pen = -1;
  double end_pen = -1;
  seq_data::SequenceData sequence_data;
  sequence_data.sequences = sequences;
  sequence_data.feature_list = feature_list;
  std::vector<fasta::SequenceList> alignment;
  std::string sbst_mat = "BLOSUM";
  bool no_feat = false;
  d_matrix::DistanceMatrix d = d_matrix::make_distance_matrix(
      f_set, sequence_data, gap_open_pen, end_pen, gap_ext_pen, 
      domain_modifier, motif_modifier, ptm_modifier,
      strct_modifier, codon_length, no_feat, sbst_mat);

  for (auto& row : d.values) {
    for (auto& item : row) {
      std::cout << item << " ";
    }
    std::cout << std::endl;
  }
}


BOOST_AUTO_TEST_SUITE_END()
