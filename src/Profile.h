#ifndef PROFILE_H
#define PROFILE_H

#include "types.h"
#include <iostream>
#include <string>
#include <vector>
class Residue;
class Profile{
public:
  ///
  /// constructor; creates a Profile object with profile matrix mat
  /// 
	Profile(profile_matrix mat);
  ///
  /// constructor; creates an empty profile
  ///
	Profile();
  ///
  /// builds a profile from the matrix of occurences and the substitution 
  /// matrix with appropriate weights
  ///
	void processProfile(sequenceList& alignment, 
                      const identitiesList& sequenceIdentityValues, 
                      bool weightsModeOn);
  ///
  /// returns the profile matrix
  ///
	profile_matrix getMatrix() const;
  ///
  /// returns a score for amino acid aAcid on a certain position
  ///
	double getElement(int position, char aAcid);
  ///
  /// returns a score for amino acid with index aAcidInt on a certain position
  ///
	double getElement(int aAcidInt, int position);
private:
	void createProfile(sequenceList& alignment,
                     const identitiesList& sequenceIdentityValues, 
                     bool weightsModeOn);
	profile_matrix m_prfMatrix;
};

#endif /* PROFILE_H */
