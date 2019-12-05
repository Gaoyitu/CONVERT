#pragma once
#include<vector>
#include<iostream>
#include<string>
#include<iomanip>
#include"ElementInDynamicProgrammingArray.h"
#include"Preprocessing.h"
using namespace std;
class Needleman_Wunsch
{
public:
	Needleman_Wunsch();
	~Needleman_Wunsch();
	//The normal needleman-wunsch algorithm
	vector<string> alignmentResult(string theFirstProteinSequence, string theSecondProteinSequence);
	//Semi-global alignment needleman-wunsch algorithm. The penalty for initial vacancy of long sequence is set to 0,
	//and the traceback begins with the maximum score at the end of long sequence.The short sequence is the reference sequence for alignment
	//theFirstProteinSequence is the short sequence£¬theSecondProteinSequence is the long sequence
	string semiAlignmentResult(string theFirstProteinSequence, string theSecondProteinSequence);
	//Semi-global alignment needleman-wunsch algorithm.
	//The long sequence is the reference sequence for alignment.
	//theFirstProteinSequence is the short sequence£¬theSecondProteinSequence is the long sequence
	string alignmentResultAddGap(string theFirstProteinSequence, string theSecondProteinSequence);
	double semiAlignmentScore(string theFirstProteinSequence, string theSecondProteinSequence);
	double alignmentScore(string theFirstProteinSequence, string theSecondProteinSequence);
	vector<int> sequenceToArray(string sequence);
	string arrayToSequence(vector<int> sequenceArray);
	string generatePseudoProtein(vector<string> alignmentProteinResultArray);
};


