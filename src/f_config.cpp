#include "f_config.h"

#include <algorithm>
#include <iostream>

namespace lcg = libconfig;

f_config::FeatureSettingsMap f_config::ConfParser::parse_conf_file(
    const std::string& filename,
    std::map<std::string, double>& probabilities) {
  lcg::Config cnfg;
  try
  {
    cnfg.readFile(filename.c_str());
  }
  catch (const lcg::FileIOException &fioex) {
    std::cerr << "I/O error while reading file." << std::endl;
    throw;
  }
  catch (const lcg::ParseException &pex) {
    std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
              << " - " << pex.getError() << std::endl;
    throw;
  }
  f_config::RawFeatureSettingsMap raw_map = process_config(cnfg,
                                                           probabilities);
  return process_settings(raw_map);
}


f_config::RawFeatureSettingsMap f_config::ConfParser::process_config(
    const lcg::Config& cnfg, std::map<std::string, double>& probabilities) {
  f_config::RawFeatureSettingsMap feat_config;
  try
  {
    const lcg::Setting& root  = cnfg.getRoot();
    const lcg::Setting& features = root["feature_settings"]["usr_features"];
    int count = features.getLength();
    for (int i = 0; i < count; ++i) {
      RawFeatureSettings feat_set;
      const lcg::Setting& feature = features[i];
      std::string name;
      if (!feature.lookupValue("name", name))
        continue;
      if (name.substr(0, 5) != "motif"
          && name.substr(0, 6) != "domain"
          && name.substr(0, 3) != "ptm") {
        feature.lookupValue("tag", feat_set.tag);
        bool found_add_score = feature.lookupValue("add_score",
                                                   feat_set.add_score);
        bool found_sbtrct_score = feature.lookupValue("subtract_score",
                                                      feat_set.subtract_score);
        if (!(found_add_score || found_sbtrct_score))
          continue;
        feature.lookupValue("pattern", feat_set.pattern);

        lcg::Setting& add_features_set = feature["add_features"];
        for (int j = 0; j < add_features_set.getLength(); ++j) {
          feat_set.add_features.push_back(add_features_set[j]);
        }
        lcg::Setting& add_tags_set = feature["add_tags"];
        for (int j = 0; j < add_tags_set.getLength(); ++j) {
          feat_set.add_tags.push_back(add_tags_set[j]);
        }
        lcg::Setting& add_exceptions_set = feature["add_exceptions"];
        for (int j = 0; j < add_exceptions_set.getLength(); ++j) {
          feat_set.add_exceptions.push_back(add_exceptions_set[j]);
        }
        lcg::Setting& subtract_features_set = feature["subtract_features"];
        for (int j = 0; j < subtract_features_set.getLength(); ++j) {
          feat_set.subtract_features.push_back(subtract_features_set[j]);
        }
        lcg::Setting& subtract_tags_set = feature["subtract_tags"];
        for (int j = 0; j < subtract_tags_set.getLength(); ++j) {
          feat_set.subtract_tags.push_back(subtract_tags_set[j]);
        }
        lcg::Setting& subtract_exceptions_set = feature["subtract_exceptions"];
        for (int j = 0; j < subtract_exceptions_set.getLength(); ++j) {
          feat_set.subtract_exceptions.push_back(subtract_exceptions_set[j]);
        }
      }
      else if (name.substr(0, 5) == "motif") {
        bool found_add_score = feature.lookupValue("add_score",
                                                   probabilities[name]);
        if (!found_add_score) {
          std::cout << "Motif " << name << " has no probability assigned,"
            << " setting it to 0" << std::endl;
        }
      }
      lcg::Setting& positions_set = feature["positions"];
      for (int j = 0; j <  positions_set.getLength(); ++j) {
        FeaturePositions feat_pos;
        positions_set[j].lookupValue("seq", feat_pos.seq_no);
          --feat_pos.seq_no;
        lcg::Setting& single_pos_set = positions_set[j]["pos"];
        for (int k = 0; k < single_pos_set.getLength(); ++k) {
          feat_pos.positions.push_back(single_pos_set[k]);
        }
        for (auto& pos : feat_pos.positions) {
          --pos; 
        }
        feat_set.positions.push_back(feat_pos);
      }
      feat_config[name] = feat_set;
    }
  }
  catch(const lcg::SettingNotFoundException &nfex)
  {
    std::cerr << "Setting not found" << std::endl;
    throw;
  }
  return feat_config;
}


f_config::FeatureSettingsMap f_config::ConfParser::process_settings(
    const f_config::RawFeatureSettingsMap raw_map) {
  f_config::FeatureSettingsMap processed_map;
  for (auto feat_it = raw_map.begin(); feat_it != raw_map.end(); ++feat_it) {
     FeatureSettings processed_settings; 
     processed_settings.add_score = feat_it->second.add_score;
     processed_settings.subtract_score = feat_it->second.subtract_score;
     processed_settings.positions = feat_it->second.positions;
     processed_settings.pattern = feat_it->second.pattern;
     ///
     /// filter out the features from 'add_features' and 'subtract_features'
     /// that have no settings (are not in keys in the RawFetaureSettings map)
     ///
     if (feat_it->first.substr(0, 5) != "motif"
         && feat_it->first.substr(0, 6) != "domain"
         && feat_it->first.substr(0, 3) != "ptm") {

     for (auto& feat_name : feat_it->second.add_features) {
       if (raw_map.find(feat_name) != raw_map.end()) {
         processed_settings.add_features.push_back("USR_"+ feat_name);
       }
     }

     for (auto& feat_name : feat_it->second.subtract_features) {
       if (raw_map.find(feat_name) != raw_map.end()) {
         processed_settings.subtract_features.push_back("USR_"+ feat_name);
       }
     }
     ///
     /// add all features that contain the given tag, unless they are mentioned
     /// in the exceptions
     ///
     for (auto& feat_tag : feat_it->second.add_tags) {
       for (auto feat_it_j = raw_map.begin(); feat_it_j != raw_map.end();
            ++feat_it_j) {
         if (feat_it_j->second.tag == feat_tag
             && std::find(feat_it->second.add_exceptions.begin(),
                          feat_it->second.add_exceptions.end(), 
                          feat_it_j->first) 
                == feat_it->second.add_exceptions.end()
             && std::find(processed_settings.add_features.begin(), 
                          processed_settings.add_features.end(),
                          "USR_" + feat_it_j->first) 
                == processed_settings.add_features.end()) {
           processed_settings.add_features.push_back(
               "USR_" + feat_it_j->first);
         }
       }
     }

     for (auto& feat_tag : feat_it->second.subtract_tags) {
       for (auto feat_it_j = raw_map.begin(); feat_it_j != raw_map.end();
            ++feat_it_j) {
         if (feat_it_j->second.tag == feat_tag
             && std::find(feat_it->second.subtract_exceptions.begin(),
                          feat_it->second.subtract_exceptions.end(),
                          feat_it_j->first) 
                == feat_it->second.subtract_exceptions.end()
             && std::find(processed_settings.subtract_features.begin(), 
                          processed_settings.subtract_features.end(),
                          "USR_" + feat_it_j->first) 
                == processed_settings.subtract_features.end()) {
           processed_settings.subtract_features.push_back(
               "USR_" + feat_it_j->first);
         }
       }
     }
     ///
     /// add the feature to the new map 
     ///
     processed_map["USR_" + feat_it->first] = processed_settings;
     } else {
      processed_map[feat_it->first] = processed_settings;
     }
  }
  return processed_map;
}
