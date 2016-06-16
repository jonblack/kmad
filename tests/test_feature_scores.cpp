#define BOOST_TEST_DYN_LINK

#include "src/fasta.h"
#include "src/feature_scores.h"
#include "src/f_config.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/auto_unit_test.hpp>

#include <cmath>


BOOST_AUTO_TEST_SUITE(test_feature_scores)

BOOST_AUTO_TEST_CASE(test_update_scores)
{
  FeatureNamesList feature_list = {"p_phosph0", "p_phosph1",
                                   "p_phosph2", "p_phosph3",
                                   "p_phosphP", "p_acet0",
                                   "p_acet1", "p_acet2",
                                   "p_acet3", "p_Nglyc0",
                                   "p_Nglyc1", "p_Nglyc2",
                                   "p_Nglyc3", "p_amid0",
                                   "p_amid1", "p_amid2",
                                   "p_amid3", "p_hydroxy0",
                                   "p_hydroxy1", "p_hydroxy2",
                                   "p_hydroxy3", "p_methyl0",
                                   "p_methyl1", "p_methyl2",
                                   "p_methyl3", "p_Oglyc0",
                                   "p_Oglyc1", "p_Oglyc2",
                                   "p_Oglyc3", "p_cys_bridge0",
                                   "s_a_helix", "s_turn",
                                   "s_b_ladder", "s_b_bridge",
                                   "s_310_helix", "s_pi_helix",
                                   "s_b_ladder", "m_aa",
                                   "m_ab", "m_ac",
                                   "d_aa", "d_ac",
                                   "USR_feature1", "USR_feature2"};

  // SEQUENCE S1
  fasta::Residue r1_1("AAAAdaa", std::vector<std::string>({"p_phosphP",
                                                          "m_aa"}));
  fasta::Residue r1_2("MAAAAAA", std::vector<std::string>({}));
  fasta::Residue r1_3("EAaadaa", std::vector<std::string>({"d_aa"}));
  fasta::Residue r1_4("LAAAAAA", std::vector<std::string>({}));
  // SEQUENCE S2
  fasta::Residue r2_1("AAAAdaa", std::vector<std::string>({"p_phosphP",
                                                          "m_aa"}));
  fasta::Residue r2_2("EAAAZAA", std::vector<std::string>({"p_phosph0"}));
  fasta::Residue r2_3("EAaaZAA", std::vector<std::string>({"p_phosph0",
                                                           "d_aa"}));
  fasta::Residue r2_4("KAAAZAA", std::vector<std::string>({"p_phosph0"}));
  // SEQUENCE S3
  fasta::Residue r3_1("AAAAZaa", std::vector<std::string>({"p_phosph0",
                                                           "m_aa"}));
  fasta::Residue r3_2("MAAAZAA", std::vector<std::string>({"p_phosph0"}));
  fasta::Residue r3_3("EAacZAA", std::vector<std::string>({"p_phosph0",
                                                           "d_ac"}));
  fasta::Residue r3_4("LAaaZaa", std::vector<std::string>({"p_phosph0",
                                                           "m_aa",
                                                           "d_aa"}));
  // SEQUENCE S4
  fasta::Residue r4_1("MAAAZab", std::vector<std::string>({"p_phosph0",
                                                           "m_ab"}));
  fasta::Residue r4_2("MAAAZAA", std::vector<std::string>({"p_phosph0"}));
  fasta::Residue r4_3("KAaaZAA", std::vector<std::string>({"p_phosph0",
                                                           "d_aa"}));
  fasta::Residue r4_4("LAAAAaa", std::vector<std::string>({"m_aa"}));

  // SEQUENCE S5
  fasta::Residue r5_1("AAAAAac", std::vector<std::string>({"m_ac"}));
  fasta::Residue r5_2("MAAAZAA", std::vector<std::string>({"p_phosph0"}));
  fasta::Residue r5_3("AAAAAAA", std::vector<std::string>({"USR_feature2"}));
  fasta::Residue r5_4("LAAAAAA", std::vector<std::string>({"USR_feature1"}));

  fasta::Sequence s1 = fasta::make_sequence({r1_1, r1_2, r1_3, r1_4});
  fasta::Sequence s2 = fasta::make_sequence({r2_1, r2_2, r2_3, r2_4});
  fasta::Sequence s3 = fasta::make_sequence({r3_1, r3_2, r3_3, r3_4});
  fasta::Sequence s4 = fasta::make_sequence({r4_1, r4_2, r4_3, r4_4});
  fasta::Sequence s5 = fasta::make_sequence({r5_1, r5_2, r5_3, r5_4});

  fasta::SequenceList sequences = {s1, s2, s3, s4, s5};

  std::unordered_map<std::string, double> probs = {{"m_aa", 1.0},
                                         {"m_ab", 0.5},
                                         {"m_ac", 0.8}};

  t::SettingsMap aln_params = {
          {"ptm", 10},
          {"domain", 4},
          {"motif", 3},
          {"strct", 0}
  };
  FeatureScores f_profile(feature_list, aln_params, probs);
  f_config::FeatureSettings settings1;
  settings1.add_score = 1;
  settings1.subtract_score = 0;
  settings1.add_features = {"USR_feature1"};
  f_config::FeatureSettings settings2;
  settings2.add_score = 1;
  settings2.subtract_score = 1;
  settings2.add_features = {"USR_feature2"};
  settings2.subtract_features = {"USR_feature1"};
  f_config::FeatureSettingsMap s_map;
  s_map = {{"USR_feature1", settings1},
           {"USR_feature2", settings2}};
  std::vector<double> identities(sequences.size(), 1.0);
  bool fade_out = false;
  f_profile.update_scores(sequences, s_map, identities, fade_out);
  std::unordered_map<std::string, Scores> p = f_profile.get_scores();
  std::unordered_map<std::string, Scores> expected_scores;
  expected_scores = {{"p_phosph0", {5.2, 8, 6, 4}},
                     {"p_phosph1", {4.68, 7.2, 5.4, 3.6}},
                     {"p_phosph2", {4.16, 6.4, 4.8, 3.2}},
                     {"p_phosph3", {3.64, 5.6, 4.2, 2.8}},
                     {"p_phosphP", {1.56, 2.4, 1.8, 1.2}},
                     {"p_acet0", {0, 0, 0, 0}},
                     {"p_acet1", {0, 0, 0, 0}},
                     {"p_acet2", {0, 0, 0, 0}},
                     {"p_acet3", {0, 0, 0, 0}},
                     {"p_Nglyc0", {0, 0, 0, 0}},
                     {"p_Nglyc1", {0, 0, 0, 0}},
                     {"p_Nglyc2", {0, 0, 0, 0}},
                     {"p_Nglyc3", {0, 0, 0, 0}},
                     {"p_amid0", {0, 0, 0, 0}},
                     {"p_amid1", {0, 0, 0, 0}},
                     {"p_amid2", {0, 0, 0, 0}},
                     {"p_amid3", {0, 0, 0, 0}},
                     {"p_hydroxy0", {0, 0, 0, 0}},
                     {"p_hydroxy1", {0, 0, 0, 0}},
                     {"p_hydroxy2", {0, 0, 0, 0}},
                     {"p_hydroxy3", {0, 0, 0, 0}},
                     {"p_methyl0", {0, 0, 0, 0}},
                     {"p_methyl1", {0, 0, 0, 0}},
                     {"p_methyl2", {0, 0, 0, 0}},
                     {"p_methyl3", {0, 0, 0, 0}},
                     {"p_Oglyc0", {0, 0, 0, 0}},
                     {"p_Oglyc1", {0, 0, 0, 0}},
                     {"p_Oglyc2", {0, 0, 0, 0}},
                     {"p_Oglyc3", {0, 0, 0, 0}},
                     {"p_cys_bridge0", {0, 0, 0, 0}},
                     {"s_a_helix", {0, 0, 0, 0}},
                     {"s_turn", {0, 0, 0, 0}},
                     {"s_b_ladder", {0, 0, 0, 0}},
                     {"s_b_bridge", {0, 0, 0, 0}},
                     {"s_310_helix", {0, 0, 0, 0}},
                     {"s_pi_helix", {0, 0, 0, 0}},
                     {"s_b_ladder", {0, 0, 0, 0}},
                     {"m_aa", {1.8, 0, 0, 1.2}},
                     {"m_ab", {0.3, 0, 0, 0}},
                     {"m_ac", {0.48, 0, 0, 0}},
                     {"d_aa", {0, 0, 1.6, 0.8}},
                     {"d_ac", {0, 0, -1.6, -0.8}},
                     {"USR_feature2", {0, 0,  0.2, -0.2}},
                     {"USR_feature1", {0, 0, 0, 0.2}}};
  for (auto& feat: feature_list){
    for (size_t i = 0; i < p[feat].size(); ++i) {
      BOOST_CHECK(std::abs(p[feat][i] - expected_scores[feat][i]) < 0.001);
    }
  }
  feature_list = {"p_phosph9", "m_aa", "d_aa"};
  aln_params = {
          {"ptm", 10},
          {"domain", 4},
          {"motif", 3},
          {"strct", 0}
  };
  f_profile = FeatureScores(feature_list, aln_params, probs);
  // SEQUENCE S1
  r1_1 = fasta::Residue ("AAAAdaa", std::vector<std::string>({"p_phosph9",
                                                              "m_aa"}));
  r1_2 = fasta::Residue ("MAAAAAA", std::vector<std::string>({}));
  r1_3 = fasta::Residue ("EAaadaa", std::vector<std::string>({"d_aa"}));
  r1_4 = fasta::Residue ("LAAAAAA", std::vector<std::string>({}));
  s1 = fasta::make_sequence({r1_1, r1_2, r1_3, r1_4});
  sequences = {s1};
  BOOST_CHECK_THROW(f_profile.update_scores(sequences, s_map, identities,
                                            fade_out),
                    std::invalid_argument);

}

BOOST_AUTO_TEST_SUITE_END()
