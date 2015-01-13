#ifndef SEQDATA_H
#define SEQDATA_H

#include "fasta.h"
#include "f_config.h"

#include <iostream>


namespace seq_data {
  struct SequenceData {
    fasta::SequenceList sequences;
    std::map<std::string, double> probabilities;  
    std::vector<std::string> feature_list;
  };

  SequenceData process_fasta_data(const fasta::FastaData& fasta_data,
                                  const f_config::FeatureSettingsMap& f_set);
  std::vector<std::string> make_feature_list(
      const fasta::SequenceList& sequences);
}

#endif /* SEQDATA_H */
