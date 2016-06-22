#ifndef SBST_MAT_H
#define SBST_MAT_H


// TODO: Look at Maarten's code to see what he did with this kind of stuff.
namespace substitution_matrix {
  typedef std::unordered_map<char, std::vector<double>> SimilarityScoresMap;
  static const std::vector<char> ALPHABET = {'A','R','N','D','C','Q','E','G',
                                             'H','I','L','K','M','F','P','S',
                                             'T','W','Y','V'};
  static const SimilarityScoresMap BLOSUM = {
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
  static const SimilarityScoresMap DISORDER = {
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
#endif /* SBTS_MAT_H */
