#include <iostream>
#include <string>
#include <vector>
class FeaturesProfile{
public:
	FeaturesProfile(int, int, std::vector<std::string>, std::vector<double>);
	//getters
	double getElement(int, std::string);
	double getElement(int, int);	
	double getScore(int,std::string);
	std::vector<std::vector<double> > getMatrix();
	std::vector<std::string> getMotifsIDs();
	std::vector<double> getMotifsProbs();
	double getDomainScore();
	double getPhosphScore();
	//setters
	void setMatrix(std::vector<std::vector<double> >);
	void printProfile();
	void createProfile(const std::vector< std::vector<std::string> >&,const std::vector<bool>&, const std::vector<double>&, bool, int);
	void expandListOfFeatures(const std::vector<std::string>&, int);
private:
	std::vector< std::vector<double> > prfMatrix;
	void countOccurences(const std::vector< std::vector<std::string> >&,const std::vector<bool>&, int);
	void countOccurences(const std::vector< std::vector<std::string> >&,const std::vector<double>&, int);
	int findFeaturesIndex(std::string);
	std::string name(std::string,int);
	int domainScore;
	int phosphScore;
	std::vector<std::string> motifs_ids;
	std::vector<double> motifs_probs;
};
