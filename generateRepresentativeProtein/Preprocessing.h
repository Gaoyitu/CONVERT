#pragma once
#include<string>
#include<fstream>
#include<vector>
using namespace std;
class Preprocessing
{
public:
	Preprocessing();
	~Preprocessing();
	int getNumberOfProteins(string theNameOfFile);
	vector<string> getSequencesInSameFamily(string familyId, string theNameOfInputFile);
	vector<string> getSequencesInSameFamilyUpdate(string familyId, string theNameOfInputFile);
	vector<string> getSequencesInSameSuperFamily(string superFamilyId, string theNameOfInputFile);
	void formatTheFile(string theNameOfInputFile, string theNameOfOutputFile);
	string getTheFamliyIdOfProteinSequenceInFasta(string proteinSequenceInSCOPInformation);
	vector<string> getTheAllFamilyIdInData(string theNameOfFile);
	vector<string> getAllSuperFamilyInData(string theNameOfFile);
	string getTheLengthestSequenceInSameFamily(vector<string> theSameFamily);
	string getTheShortestSequenceInSameFamily(vector<string> theSameFamily);
	string getTheSuitableAlignmentSequenceInSameFamily(vector<string> theSameFamily, int averageLength,
		int theLengthOfTheLengthestSequence, int theLengthOfTheShortestSequence);
	int getTheAverageLengthOfSameFamilySequence(vector<string> theSameFamily);
private:
};


