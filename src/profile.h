#ifndef PROFILE_H
#define PROFILE_H

#include "fasta.h"


namespace profile {
  typedef std::map<char, std::vector<double>> SimilarityScoresMap;
  typedef std::map<char, std::vector<double>> ProfileMap;

  ProfileMap create_profile(const fasta::SequenceList& sequences);
  ProfileMap create_score_profile(const fasta::SequenceList& sequences,
                                  const std::string& sbst_mat);
  double get_score(const ProfileMap& p, int position, char aa);
  double profile_score(const ProfileMap& profile1, const ProfileMap& profile2,
                       int pos1, int pos2, const SimilarityScoresMap*);
  double calc_sum(const ProfileMap& profile, int pos);

  ProfileMap make_new_profile(fasta::SequenceList sequences,
    const fasta::SequenceList& sequences,
    const ProfileMap& profile1, const ProfileMap& profile2);
  ProfileMap make_scores(ProfileMap& p, const std::string& sbst_mat);
}

#endif /* PROFILE_H */
