#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<fstream>
#include<iomanip>
#include<set>
#include<map>
#include<unordered_set>
#include"Preprocessing.h"
#include"Needleman_Wunsch.h"
#include<Windows.h>
using namespace std;
//The input file is a FASTA form file of protein sequences, and the output file contains representing protein sequences of all families.
void generateRepresentativeProtein(string theNameOfInputFile, string theNameOfOutputFile);

int main()
{
	generateRepresentativeProtein("SCOP40Format.txt","SCOP40representativeProtein_1.txt");
	generateRepresentativeProtein("SCOP95Format.txt", "SCOP95representativeProtein_1.txt");
	generateRepresentativeProtein("SCOPe40Format.txt", "SCOPe40representativeProtein_1.txt");
	generateRepresentativeProtein("SCOPe95Format.txt", "SCOPe95representativeProtein_1.txt");
	system("pause");
	return EXIT_SUCCESS;
}

void generateRepresentativeProtein(string theNameOfInputFile, string theNameOfOutputFile)
{
	Preprocessing p;

	Needleman_Wunsch nw;
	string thePathOfOutputFile = "./" + theNameOfOutputFile;
	ofstream ofs(thePathOfOutputFile, ios::out);
	vector<string> allProteinId;
	vector<string> allSequencesInSameFamily;
	allProteinId = p.getTheAllFamilyIdInData(theNameOfInputFile);//put all family id in a vector

	vector<string> alignmentResultInSameFamily;
	string theAlignmentSequence;
	int theLengthOfShortestSequence;
	int theLengthOfLongthestSequence;
	int averageLength;
	

	for (int i = 0; i < allProteinId.size(); ++i)
	{

		allSequencesInSameFamily = p.getSequencesInSameFamily(allProteinId[i], theNameOfInputFile);

		if (allSequencesInSameFamily.size() > 1)
		{
			theLengthOfShortestSequence = p.getTheShortestSequenceInSameFamily(allSequencesInSameFamily).size();
			theLengthOfLongthestSequence = p.getTheLengthestSequenceInSameFamily(allSequencesInSameFamily).size();
			if (theLengthOfLongthestSequence - theLengthOfShortestSequence < 20)
			{
				theAlignmentSequence = p.getTheShortestSequenceInSameFamily(allSequencesInSameFamily);
				for (int j = 0; j < allSequencesInSameFamily.size(); ++j)
				{
					alignmentResultInSameFamily.push_back(nw.semiAlignmentResult(theAlignmentSequence,
						allSequencesInSameFamily[j]));
				}
				ofs << allProteinId[i] << endl;
				ofs << nw.generatePseudoProtein(alignmentResultInSameFamily) << endl;
			}
			else
			{
				averageLength = p.getTheAverageLengthOfSameFamilySequence(allSequencesInSameFamily);
				theAlignmentSequence = p.getTheSuitableAlignmentSequenceInSameFamily(allSequencesInSameFamily, averageLength,
					theLengthOfLongthestSequence, theLengthOfShortestSequence);
				for (int j = 0; j < allSequencesInSameFamily.size(); ++j)
				{
					if (allSequencesInSameFamily[j].size() >= theAlignmentSequence.size())
					{
						alignmentResultInSameFamily.push_back(nw.semiAlignmentResult(theAlignmentSequence, allSequencesInSameFamily[j]));
					}
					else
					{
						alignmentResultInSameFamily.push_back(nw.alignmentResultAddGap(allSequencesInSameFamily[j], theAlignmentSequence));
					}

				}
				ofs << allProteinId[i] << endl;
				ofs << nw.generatePseudoProtein(alignmentResultInSameFamily) << endl;
			}
		}
		else
		{
			ofs << allProteinId[i] << endl;
			ofs << allSequencesInSameFamily[0] << endl;
		}
		alignmentResultInSameFamily.clear();
		allSequencesInSameFamily.clear();
		vector<string>().swap(alignmentResultInSameFamily);
		vector<string>().swap(allSequencesInSameFamily);
		cout << "the " << allProteinId[i] <<" family has been finished."<< endl;
	}
	ofs.close();
}


