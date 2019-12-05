#include "Needleman_Wunsch.h"




Needleman_Wunsch::Needleman_Wunsch()
{
}


Needleman_Wunsch::~Needleman_Wunsch()
{
}

vector<string> Needleman_Wunsch::alignmentResult(string theFirstProteinSequence, string theSecondProteinSequence)
{
	//int *firstSequenceToArray = new int[theFirstProteinSequence.size()];

	vector<string> result;
	const double gapOpenPenalty = -10;
	const double gapExtendPenalty = -0.5;
	const string up = "upper";
	const string le = "left";
	const string upLe = "upperLeft";
	const int BLOSUM62[20][20] = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3 - 1,3},
	{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
	{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}
	};

	vector<int> theFirstSequenceArray = this->sequenceToArray(theFirstProteinSequence);
	vector<int> theSecondSequenceArray = this->sequenceToArray(theSecondProteinSequence);


	/*vector<char> theAlignmentResultOfFirstSequence;
	vector<char> theAlignmentResultOfSecondSequence;*/
	string theAlignmentResultOfFirstSequence;
	string theAlignmentResultOfSecondSequence;
	ElementInDynamicProgrammingArray **dynamicProgrammingArray =
		new ElementInDynamicProgrammingArray*[theFirstProteinSequence.size() + 1];
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i] = new ElementInDynamicProgrammingArray[theSecondProteinSequence.size() + 1];
	}



	dynamicProgrammingArray[0][0] = ElementInDynamicProgrammingArray(false, 0, "wu");
	dynamicProgrammingArray[1][0] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, up);
	dynamicProgrammingArray[0][1] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, le);

	for (int i = 2; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i][0] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[i - 1][0].getValue() + gapExtendPenalty, up);
	}
	for (int i = 2; i < (theSecondProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[0][i] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[0][i - 1].getValue() + gapExtendPenalty, le);
	}

	double upperLeft = 0;
	double left = 0;
	double upper = 0;

	for (int i = 0; i < theFirstSequenceArray.size(); i++)
	{
		for (int j = 0; j < theSecondSequenceArray.size(); j++)
		{

			upperLeft = dynamicProgrammingArray[i][j].getValue() + BLOSUM62[theFirstSequenceArray[i]][theSecondSequenceArray[j]];

			if (dynamicProgrammingArray[i][j + 1].getIsGap() == true)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i][j + 1].getIsGap() == false)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapOpenPenalty;
			}
			if (dynamicProgrammingArray[i + 1][j].getIsGap() == true)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i + 1][j].getIsGap() == false)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapOpenPenalty;
			}

			if (upperLeft >= left && upperLeft >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(false, upperLeft, upLe);
			}

			else if (left > upperLeft && left >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, left, le);
			}

			else if (upper > upperLeft && upper > left)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, upper, up);
			}
		}
	}
	//for (int i = 0; i<(theFirstProteinSequence.size()+1); i++)
	//{
	//	for (int j = 0; j<(theSecondProteinSequence.size()+1); j++)
	//	{
	//		//cout << dynamicProgrammingArray[i][j].getValue() << ","<<dynamicProgrammingArray[i][j].getIsGap()<<","<<dynamicProgrammingArray[i][j].getFromWhere()<<"  ";
	//		cout << dynamicProgrammingArray[i][j].getValue() << "   ";
	//	}
	//	cout << endl;
	//}
	int theLengthOfFirstSequence = theFirstProteinSequence.size();
	int theLengthOfSecondSequence = theSecondProteinSequence.size();

	while (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare("wu") != 0)
	{
		if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(upLe) == 0)
		{
			theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			theAlignmentResultOfSecondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(up) == 0)
		{
			theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			theAlignmentResultOfSecondSequence.push_back('-');
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(le) == 0)
		{
			theAlignmentResultOfFirstSequence.push_back('-');
			theAlignmentResultOfSecondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}

		/*for (int i = 0; i<theAlignmentResultOfFirstSequence.size(); ++i)
		{
			cout << theAlignmentResultOfFirstSequence[i];
		}
		cout << endl;
		for (int i =0;i< theAlignmentResultOfSecondSequence.size(); ++i)
		{
			cout << theAlignmentResultOfSecondSequence[i];
		}
		cout << endl;*/

	}

	string theAlignResutlOfPositiveDirection1;
	string theAlignResutlOfPositiveDirection2;
	for (int i = theAlignmentResultOfFirstSequence.size() - 1; i >= 0; --i)
	{
		theAlignResutlOfPositiveDirection1 = theAlignResutlOfPositiveDirection1 + theAlignmentResultOfFirstSequence.at(i);
	}
	for (int i = theAlignmentResultOfSecondSequence.size() - 1; i >= 0; --i)
	{
		theAlignResutlOfPositiveDirection2 = theAlignResutlOfPositiveDirection2 + theAlignmentResultOfSecondSequence.at(i);
	}
	result.push_back(theAlignResutlOfPositiveDirection1);
	result.push_back(theAlignResutlOfPositiveDirection2);

	/*for (int i = theAlignmentResultOfFirstSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfFirstSequence[i];
	}
	cout << endl;
	for (int i = theAlignmentResultOfSecondSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfSecondSequence[i];
	}
	cout << endl;*/

	//cout << dynamicProgrammingArray[theFirstProteinSequence.size()][theSecondProteinSequence.size()].getValue();
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		delete[] dynamicProgrammingArray[i];
		dynamicProgrammingArray[i] = NULL;
	}
	delete[] dynamicProgrammingArray;
	dynamicProgrammingArray = NULL;


	return result;
}

string Needleman_Wunsch::semiAlignmentResult(string theFirstProteinSequence, string theSecondProteinSequence)
{
	//vector<string> result;
	const double gapOpenPenalty = -10;
	const double gapExtendPenalty = -0.5;
	const string up = "upper";
	const string le = "left";
	const string upLe = "upperLeft";
	const int BLOSUM62[20][20] = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3 - 1,3},
	{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
	{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}
	};

	vector<int> theFirstSequenceArray = this->sequenceToArray(theFirstProteinSequence);
	vector<int> theSecondSequenceArray = this->sequenceToArray(theSecondProteinSequence);


	/*vector<char> theAlignmentResultOfFirstSequence;
	vector<char> theAlignmentResultOfSecondSequence;*/
	//string theAlignmentResultOfFirstSequence;
	string theAlignmentResultOfSecondSequence;
	ElementInDynamicProgrammingArray **dynamicProgrammingArray =
		new ElementInDynamicProgrammingArray*[theFirstProteinSequence.size() + 1];
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i] = new ElementInDynamicProgrammingArray[theSecondProteinSequence.size() + 1];
	}



	dynamicProgrammingArray[0][0] = ElementInDynamicProgrammingArray(false, 0, "wu");
	dynamicProgrammingArray[1][0] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, up);
	dynamicProgrammingArray[0][1] = ElementInDynamicProgrammingArray(true, 0, le);

	for (int i = 2; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i][0] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[i - 1][0].getValue() + gapExtendPenalty, up);
	}
	for (int i = 2; i < (theSecondProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[0][i] = ElementInDynamicProgrammingArray(true,
			0, le);
	}

	double upperLeft = 0;
	double left = 0;
	double upper = 0;

	for (int i = 0; i < theFirstSequenceArray.size(); i++)
	{
		for (int j = 0; j < theSecondSequenceArray.size(); j++)
		{

			upperLeft = dynamicProgrammingArray[i][j].getValue() + BLOSUM62[theFirstSequenceArray[i]][theSecondSequenceArray[j]];

			if (dynamicProgrammingArray[i][j + 1].getIsGap() == true)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i][j + 1].getIsGap() == false)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapOpenPenalty;
			}
			if (dynamicProgrammingArray[i + 1][j].getIsGap() == true)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i + 1][j].getIsGap() == false)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapOpenPenalty;
			}

			if (upperLeft >= left && upperLeft >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(false, upperLeft, upLe);
			}

			else if (left > upperLeft && left >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, left, le);
			}

			else if (upper > upperLeft && upper > left)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, upper, up);
			}
		}
	}

	int theLengthOfFirstSequence = theFirstProteinSequence.size();
	int theLengthOfSecondSequence = theSecondProteinSequence.size();
	int theMaxValueIndex = theLengthOfSecondSequence;
	double theMaxValue = dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getValue();
	for (int index = theLengthOfSecondSequence; index > 0; --index)
	{
		if (dynamicProgrammingArray[theLengthOfFirstSequence][index - 1].getValue() > theMaxValue)
		{
			theMaxValue = dynamicProgrammingArray[theLengthOfFirstSequence][index - 1].getValue();
			theMaxValueIndex = index - 1;
		}
	}
	for (int c = theLengthOfSecondSequence; c > theMaxValueIndex; --c)
	{
		dynamicProgrammingArray[theLengthOfFirstSequence][c].setFormWhere(le);
	}
	//for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	//{
	//	for (int j = 0; j < (theSecondProteinSequence.size() + 1); j++)
	//	{
	//		//cout << dynamicProgrammingArray[i][j].getValue() << ","<<dynamicProgrammingArray[i][j].getIsGap()<<","<<dynamicProgrammingArray[i][j].getFromWhere()<<"  ";
	//		cout << setw(8) << dynamicProgrammingArray[i][j].getValue() << setw(7) << dynamicProgrammingArray[i][j].getFromWhere() << "   ";
	//	}
	//	cout << endl;
	//}
	while (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare("wu") != 0)
	{

		if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(upLe) == 0)
		{
			//theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			theAlignmentResultOfSecondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(up) == 0)
		{
			//theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			theAlignmentResultOfSecondSequence.push_back('-');
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(le) == 0)
		{
			//theAlignmentResultOfFirstSequence.push_back('-');
			//theAlignmentResultOfSec0ondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}

		/*for (int i = 0; i<theAlignmentResultOfFirstSequence.size(); ++i)
		{
			cout << theAlignmentResultOfFirstSequence[i];
		}
		cout << endl;
		for (int i =0;i< theAlignmentResultOfSecondSequence.size(); ++i)
		{
			cout << theAlignmentResultOfSecondSequence[i];
		}
		cout << endl;*/

	}
	//string theAlignResutlOfPositiveDirection1;
	string theAlignResutlOfPositiveDirection2;
	/*for (int i=theAlignmentResultOfFirstSequence.size()-1;i>=0; --i)
	{
		theAlignResutlOfPositiveDirection1 = theAlignResutlOfPositiveDirection1 + theAlignmentResultOfFirstSequence.at(i);
	}*/
	for (int i = theAlignmentResultOfSecondSequence.size() - 1; i >= 0; --i)
	{
		theAlignResutlOfPositiveDirection2 = theAlignResutlOfPositiveDirection2 + theAlignmentResultOfSecondSequence.at(i);
	}
	//result.push_back(theAlignResutlOfPositiveDirection1);
	//result.push_back(theAlignResutlOfPositiveDirection2);

	/*for (int i = theAlignmentResultOfFirstSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfFirstSequence[i];
	}
	cout << endl;
	for (int i = theAlignmentResultOfSecondSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfSecondSequence[i];
	}
	cout << endl;*/

	//cout << dynamicProgrammingArray[theFirstProteinSequence.size()][theSecondProteinSequence.size()].getValue();
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		delete[] dynamicProgrammingArray[i];
		dynamicProgrammingArray[i] = NULL;
	}
	delete[] dynamicProgrammingArray;
	dynamicProgrammingArray = NULL;


	return theAlignResutlOfPositiveDirection2;
}

string Needleman_Wunsch::alignmentResultAddGap(string theFirstProteinSequence, string theSecondProteinSequence)
{
	const double gapOpenPenalty = -10;
	const double gapExtendPenalty = -0.5;
	const string up = "upper";
	const string le = "left";
	const string upLe = "upperLeft";
	const int BLOSUM62[20][20] = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3 - 1,3},
	{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
	{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}
	};

	vector<int> theFirstSequenceArray = this->sequenceToArray(theFirstProteinSequence);
	vector<int> theSecondSequenceArray = this->sequenceToArray(theSecondProteinSequence);


	/*vector<char> theAlignmentResultOfFirstSequence;
	vector<char> theAlignmentResultOfSecondSequence;*/
	string theAlignmentResultOfFirstSequence;
	//string theAlignmentResultOfSecondSequence;
	ElementInDynamicProgrammingArray **dynamicProgrammingArray =
		new ElementInDynamicProgrammingArray*[theFirstProteinSequence.size() + 1];
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i] = new ElementInDynamicProgrammingArray[theSecondProteinSequence.size() + 1];
	}



	dynamicProgrammingArray[0][0] = ElementInDynamicProgrammingArray(false, 0, "wu");
	dynamicProgrammingArray[1][0] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, up);
	dynamicProgrammingArray[0][1] = ElementInDynamicProgrammingArray(true, 0, le);

	for (int i = 2; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i][0] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[i - 1][0].getValue() + gapExtendPenalty, up);
	}
	for (int i = 2; i < (theSecondProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[0][i] = ElementInDynamicProgrammingArray(true,
			0, le);
	}

	double upperLeft = 0;
	double left = 0;
	double upper = 0;

	for (int i = 0; i < theFirstSequenceArray.size(); i++)
	{
		for (int j = 0; j < theSecondSequenceArray.size(); j++)
		{

			upperLeft = dynamicProgrammingArray[i][j].getValue() + BLOSUM62[theFirstSequenceArray[i]][theSecondSequenceArray[j]];

			if (dynamicProgrammingArray[i][j + 1].getIsGap() == true)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i][j + 1].getIsGap() == false)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapOpenPenalty;
			}
			if (dynamicProgrammingArray[i + 1][j].getIsGap() == true)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i + 1][j].getIsGap() == false)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapOpenPenalty;
			}

			if (upperLeft >= left && upperLeft >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(false, upperLeft, upLe);
			}

			else if (left > upperLeft && left >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, left, le);
			}

			else if (upper > upperLeft && upper > left)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, upper, up);
			}
		}
	}
	//for (int i = 0; i<(theFirstProteinSequence.size()+1); i++)
	//{
	//	for (int j = 0; j<(theSecondProteinSequence.size()+1); j++)
	//	{
	//		//cout << dynamicProgrammingArray[i][j].getValue() << ","<<dynamicProgrammingArray[i][j].getIsGap()<<","<<dynamicProgrammingArray[i][j].getFromWhere()<<"  ";
	//		cout << dynamicProgrammingArray[i][j].getValue() << "   ";
	//	}
	//	cout << endl;
	//}
	int theLengthOfFirstSequence = theFirstProteinSequence.size();
	int theLengthOfSecondSequence = theSecondProteinSequence.size();
	int theMaxValueIndex = theLengthOfSecondSequence;
	double theMaxValue = dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getValue();
	for (int index = theLengthOfSecondSequence; index > 0; --index)
	{
		if (dynamicProgrammingArray[theLengthOfFirstSequence][index - 1].getValue() > theMaxValue)
		{
			theMaxValue = dynamicProgrammingArray[theLengthOfFirstSequence][index - 1].getValue();
			theMaxValueIndex = index - 1;
		}
	}
	for (int c = theLengthOfSecondSequence; c > theMaxValueIndex; --c)
	{
		dynamicProgrammingArray[theLengthOfFirstSequence][c].setFormWhere(le);
	}
	//for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	//{
	//	for (int j = 0; j < (theSecondProteinSequence.size() + 1); j++)
	//	{
	//		//cout << dynamicProgrammingArray[i][j].getValue() << ","<<dynamicProgrammingArray[i][j].getIsGap()<<","<<dynamicProgrammingArray[i][j].getFromWhere()<<"  ";
	//		cout << setw(8) << dynamicProgrammingArray[i][j].getValue() << setw(7) << dynamicProgrammingArray[i][j].getFromWhere() << "   ";
	//	}
	//	cout << endl;
	//}
	while (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare("wu") != 0)
	{

		if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(upLe) == 0)
		{
			theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			//theAlignmentResultOfSecondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(up) == 0)
		{
			//theAlignmentResultOfFirstSequence.push_back(theFirstProteinSequence.at(theLengthOfFirstSequence - 1));
			//theAlignmentResultOfSecondSequence.push_back('-');
			theLengthOfFirstSequence = theLengthOfFirstSequence - 1;
		}
		else if (dynamicProgrammingArray[theLengthOfFirstSequence][theLengthOfSecondSequence].getFromWhere().compare(le) == 0)
		{
			theAlignmentResultOfFirstSequence.push_back('-');
			//theAlignmentResultOfSec0ondSequence.push_back(theSecondProteinSequence.at(theLengthOfSecondSequence - 1));
			theLengthOfSecondSequence = theLengthOfSecondSequence - 1;
		}

		/*for (int i = 0; i<theAlignmentResultOfFirstSequence.size(); ++i)
		{
			cout << theAlignmentResultOfFirstSequence[i];
		}
		cout << endl;
		for (int i =0;i< theAlignmentResultOfSecondSequence.size(); ++i)
		{
			cout << theAlignmentResultOfSecondSequence[i];
		}
		cout << endl;*/

	}
	string theAlignResutlOfPositiveDirection1;
	//string theAlignResutlOfPositiveDirection2;
	/*for (int i=theAlignmentResultOfFirstSequence.size()-1;i>=0; --i)
	{
		theAlignResutlOfPositiveDirection1 = theAlignResutlOfPositiveDirection1 + theAlignmentResultOfFirstSequence.at(i);
	}*/
	for (int i = theAlignmentResultOfFirstSequence.size() - 1; i >= 0; --i)
	{
		theAlignResutlOfPositiveDirection1 = theAlignResutlOfPositiveDirection1 + theAlignmentResultOfFirstSequence.at(i);
	}
	//result.push_back(theAlignResutlOfPositiveDirection1);
	//result.push_back(theAlignResutlOfPositiveDirection2);

	/*for (int i = theAlignmentResultOfFirstSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfFirstSequence[i];
	}
	cout << endl;
	for (int i = theAlignmentResultOfSecondSequence.size()-1; i >= 0; --i)
	{
		cout << theAlignmentResultOfSecondSequence[i];
	}
	cout << endl;*/

	//cout << dynamicProgrammingArray[theFirstProteinSequence.size()][theSecondProteinSequence.size()].getValue();
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		delete[] dynamicProgrammingArray[i];
		dynamicProgrammingArray[i] = NULL;
	}
	delete[] dynamicProgrammingArray;
	dynamicProgrammingArray = NULL;


	return theAlignResutlOfPositiveDirection1;
}

double Needleman_Wunsch::semiAlignmentScore(string theFirstProteinSequence, string theSecondProteinSequence)
{

	double score;
	const double gapOpenPenalty = -10;
	const double gapExtendPenalty = -0.5;
	const string up = "upper";
	const string le = "left";
	const string upLe = "upperLeft";
	const int BLOSUM62[20][20] = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3 - 1,3},
	{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
	{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}
	};

	vector<int> theFirstSequenceArray = this->sequenceToArray(theFirstProteinSequence);
	vector<int> theSecondSequenceArray = this->sequenceToArray(theSecondProteinSequence);

	string theAlignmentResultOfSecondSequence;
	ElementInDynamicProgrammingArray **dynamicProgrammingArray =
		new ElementInDynamicProgrammingArray*[theFirstProteinSequence.size() + 1];
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i] = new ElementInDynamicProgrammingArray[theSecondProteinSequence.size() + 1];
	}



	dynamicProgrammingArray[0][0] = ElementInDynamicProgrammingArray(false, 0, "wu");
	dynamicProgrammingArray[1][0] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, up);
	dynamicProgrammingArray[0][1] = ElementInDynamicProgrammingArray(true, 0, le);

	for (int i = 2; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i][0] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[i - 1][0].getValue() + gapExtendPenalty, up);
	}
	for (int i = 2; i < (theSecondProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[0][i] = ElementInDynamicProgrammingArray(true,
			0, le);
	}

	double upperLeft = 0;
	double left = 0;
	double upper = 0;

	for (int i = 0; i < theFirstSequenceArray.size(); i++)
	{
		for (int j = 0; j < theSecondSequenceArray.size(); j++)
		{

			upperLeft = dynamicProgrammingArray[i][j].getValue() + BLOSUM62[theFirstSequenceArray[i]][theSecondSequenceArray[j]];

			if (dynamicProgrammingArray[i][j + 1].getIsGap() == true)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i][j + 1].getIsGap() == false)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapOpenPenalty;
			}
			if (dynamicProgrammingArray[i + 1][j].getIsGap() == true)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i + 1][j].getIsGap() == false)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapOpenPenalty;
			}

			if (upperLeft >= left && upperLeft >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(false, upperLeft, upLe);
			}

			else if (left > upperLeft && left >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, left, le);
			}

			else if (upper > upperLeft && upper > left)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, upper, up);
			}
		}
	}
	score = dynamicProgrammingArray[theFirstProteinSequence.size()][theSecondProteinSequence.size()].getValue();
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		delete[] dynamicProgrammingArray[i];
		dynamicProgrammingArray[i] = NULL;
	}
	delete[] dynamicProgrammingArray;
	dynamicProgrammingArray = NULL;

	return score;
}

double Needleman_Wunsch::alignmentScore(string theFirstProteinSequence, string theSecondProteinSequence)
{
	double score;
	const double gapOpenPenalty = -10;
	const double gapExtendPenalty = -0.5;
	const string up = "upper";
	const string le = "left";
	const string upLe = "upperLeft";
	const int BLOSUM62[20][20] = { {4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0},
	{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3},
	{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3},
	{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3},
	{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1},
	{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2},
	{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2},
	{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3},
	{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3},
	{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3 - 1,3},
	{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1},
	{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2},
	{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1},
	{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1},
	{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2},
	{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2},
	{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0},
	{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3},
	{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1},
	{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4}
	};

	vector<int> theFirstSequenceArray = this->sequenceToArray(theFirstProteinSequence);
	vector<int> theSecondSequenceArray = this->sequenceToArray(theSecondProteinSequence);

	ElementInDynamicProgrammingArray **dynamicProgrammingArray =
		new ElementInDynamicProgrammingArray*[theFirstProteinSequence.size() + 1];
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i] = new ElementInDynamicProgrammingArray[theSecondProteinSequence.size() + 1];
	}



	dynamicProgrammingArray[0][0] = ElementInDynamicProgrammingArray(false, 0, "wu");
	dynamicProgrammingArray[1][0] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, up);
	dynamicProgrammingArray[0][1] = ElementInDynamicProgrammingArray(true, gapOpenPenalty, le);

	for (int i = 2; i < (theFirstProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[i][0] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[i - 1][0].getValue() + gapExtendPenalty, up);
	}
	for (int i = 2; i < (theSecondProteinSequence.size() + 1); i++)
	{
		dynamicProgrammingArray[0][i] = ElementInDynamicProgrammingArray(true,
			dynamicProgrammingArray[0][i - 1].getValue() + gapExtendPenalty, le);
	}

	double upperLeft = 0;
	double left = 0;
	double upper = 0;

	for (int i = 0; i < theFirstSequenceArray.size(); i++)
	{
		for (int j = 0; j < theSecondSequenceArray.size(); j++)
		{

			upperLeft = dynamicProgrammingArray[i][j].getValue() + BLOSUM62[theFirstSequenceArray[i]][theSecondSequenceArray[j]];

			if (dynamicProgrammingArray[i][j + 1].getIsGap() == true)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i][j + 1].getIsGap() == false)
			{
				upper = dynamicProgrammingArray[i][j + 1].getValue() + gapOpenPenalty;
			}
			if (dynamicProgrammingArray[i + 1][j].getIsGap() == true)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapExtendPenalty;
			}
			else if (dynamicProgrammingArray[i + 1][j].getIsGap() == false)
			{
				left = dynamicProgrammingArray[i + 1][j].getValue() + gapOpenPenalty;
			}

			if (upperLeft >= left && upperLeft >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(false, upperLeft, upLe);
			}

			else if (left > upperLeft && left >= upper)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, left, le);
			}

			else if (upper > upperLeft && upper > left)
			{
				dynamicProgrammingArray[i + 1][j + 1] = ElementInDynamicProgrammingArray(true, upper, up);
			}
		}
	}

	score = dynamicProgrammingArray[theFirstProteinSequence.size()][theSecondProteinSequence.size()].getValue();
	for (int i = 0; i < (theFirstProteinSequence.size() + 1); i++)
	{
		delete[] dynamicProgrammingArray[i];
		dynamicProgrammingArray[i] = NULL;
	}
	delete[] dynamicProgrammingArray;
	dynamicProgrammingArray = NULL;


	return score;
}

vector<int> Needleman_Wunsch::sequenceToArray(string sequence)
{
	vector<int> sequenceArray;
	int i = 0;
	for (i; i < sequence.size(); ++i)
	{
		switch (sequence.at(i))
		{
		case 'a':
			sequenceArray.push_back(0);
			break;
		case 'r':
			sequenceArray.push_back(1);
			break;
		case 'n':
			sequenceArray.push_back(2);
			break;
		case 'd':
			sequenceArray.push_back(3);
			break;
		case 'c':
			sequenceArray.push_back(4);
			break;
		case 'q':
			sequenceArray.push_back(5);
			break;
		case 'e':
			sequenceArray.push_back(6);
			break;
		case 'g':
			sequenceArray.push_back(7);
			break;
		case 'h':
			sequenceArray.push_back(8);
			break;
		case 'i':
			sequenceArray.push_back(9);
			break;
		case 'l':
			sequenceArray.push_back(10);
			break;
		case 'k':
			sequenceArray.push_back(11);
			break;
		case 'm':
			sequenceArray.push_back(12);
			break;
		case 'f':
			sequenceArray.push_back(13);
			break;
		case 'p':
			sequenceArray.push_back(14);
			break;
		case 's':
			sequenceArray.push_back(15);
			break;
		case 't':
			sequenceArray.push_back(16);
			break;
		case 'w':
			sequenceArray.push_back(17);
			break;
		case 'y':
			sequenceArray.push_back(18);
			break;
		case 'v':
			sequenceArray.push_back(19);
			break;
		default:
			sequenceArray.push_back(0);
		}
	}
	return sequenceArray;
}

string Needleman_Wunsch::arrayToSequence(vector<int> sequenceArray)
{
	string sequence;
	int i = 0;
	for (i; i < sequenceArray.size(); i++)
	{
		switch (sequenceArray[i])
		{
		case 0:
			sequence = sequence + 'a';
			break;
		case 1:
			sequence = sequence + 'r';
			break;
		case 2:
			sequence = sequence + 'n';
			break;
		case 3:
			sequence = sequence + 'd';
			break;
		case 4:
			sequence = sequence + 'c';
			break;
		case 5:
			sequence = sequence + 'q';
			break;
		case 6:
			sequence = sequence + 'e';
			break;
		case 7:
			sequence = sequence + 'g';
			break;
		case 8:
			sequence = sequence + 'h';
			break;
		case 9:
			sequence = sequence + 'i';
			break;
		case 10:
			sequence = sequence + 'l';
			break;
		case 11:
			sequence = sequence + 'k';
			break;
		case 12:
			sequence = sequence + 'm';
			break;
		case 13:
			sequence = sequence + 'f';
			break;
		case 14:
			sequence = sequence + 'p';
			break;
		case 15:
			sequence = sequence + 's';
			break;
		case 16:
			sequence = sequence + 't';
			break;
		case 17:
			sequence = sequence + 'w';
			break;
		case 18:
			sequence = sequence + 'y';
			break;
		case 19:
			sequence = sequence + 'v';
			break;
		}
	}
	return sequence;
}

string Needleman_Wunsch::generatePseudoProtein(vector<string> alignmentProteinResultArray)
{
	string pseudoProtein;
	Preprocessing pre;
	//string theLengthSequence = pre.getTheLengthestSequenceInSameFamily(alignmentProteinResultArray);

	int theLengthSequence = alignmentProteinResultArray[0].size();
	int i = 0;
	int j = 0;


	/*for (i; i < alignmentProteinResultArray.size(); i++)
	{
		for (j = alignmentProteinResultArray[i].size(); j < theLengthSequence.size(); j++)
		{
			alignmentProteinResultArray[i] = alignmentProteinResultArray[i] + '-';
		}
	}*/

	/*for (i = 0; i < alignmentProteinResultArray.size(); i++)
	{
		cout << alignmentProteinResultArray[i] << endl;
	}*/

	vector<int> countTheTimeOfAnAcid;
	countTheTimeOfAnAcid.resize(20);
	for (i = 0; i < 20; ++i)
	{
		countTheTimeOfAnAcid[i] = 0;
	}

	int count = 0;
	int index = 0;

	for (i = 0; i < theLengthSequence; ++i)
	{
		for (j = 0; j < alignmentProteinResultArray.size(); ++j)
		{
			switch (alignmentProteinResultArray[j].at(i))
			{
			case 'a':
				++countTheTimeOfAnAcid[0];
				break;
			case 'r':
				++countTheTimeOfAnAcid[1];
				break;
			case 'n':
				++countTheTimeOfAnAcid[2];
				break;
			case 'd':
				++countTheTimeOfAnAcid[3];
				break;
			case 'c':
				++countTheTimeOfAnAcid[4];
				break;
			case 'q':
				++countTheTimeOfAnAcid[5];
				break;
			case 'e':
				++countTheTimeOfAnAcid[6];
				break;
			case 'g':
				++countTheTimeOfAnAcid[7];
				break;
			case 'h':
				++countTheTimeOfAnAcid[8];
				break;
			case 'i':
				++countTheTimeOfAnAcid[9];
				break;
			case 'l':
				++countTheTimeOfAnAcid[10];
				break;
			case 'k':
				++countTheTimeOfAnAcid[11];
				break;
			case 'm':
				++countTheTimeOfAnAcid[12];
				break;
			case 'f':
				++countTheTimeOfAnAcid[13];
				break;
			case 'p':
				++countTheTimeOfAnAcid[14];
				break;
			case 's':
				++countTheTimeOfAnAcid[15];
				break;
			case 't':
				++countTheTimeOfAnAcid[16];
				break;
			case 'w':
				++countTheTimeOfAnAcid[17];
				break;
			case 'y':
				++countTheTimeOfAnAcid[18];
				break;
			case 'v':
				++countTheTimeOfAnAcid[19];
				break;
			}
		}
		for (int m = 0; m < 20; ++m)
		{
			if (count < countTheTimeOfAnAcid[m])
			{
				count = countTheTimeOfAnAcid[m];
				index = m;
			}
		}

		switch (index)
		{
		case 0:
			pseudoProtein = pseudoProtein + 'a';
			break;
		case 1:
			pseudoProtein = pseudoProtein + 'r';
			break;
		case 2:
			pseudoProtein = pseudoProtein + 'n';
			break;
		case 3:
			pseudoProtein = pseudoProtein + 'd';
			break;
		case 4:
			pseudoProtein = pseudoProtein + 'c';
			break;
		case 5:
			pseudoProtein = pseudoProtein + 'q';
			break;
		case 6:
			pseudoProtein = pseudoProtein + 'e';
			break;
		case 7:
			pseudoProtein = pseudoProtein + 'g';
			break;
		case 8:
			pseudoProtein = pseudoProtein + 'h';
			break;
		case 9:
			pseudoProtein = pseudoProtein + 'i';
			break;
		case 10:
			pseudoProtein = pseudoProtein + 'l';
			break;
		case 11:
			pseudoProtein = pseudoProtein + 'k';
			break;
		case 12:
			pseudoProtein = pseudoProtein + 'm';
			break;
		case 13:
			pseudoProtein = pseudoProtein + 'f';
			break;
		case 14:
			pseudoProtein = pseudoProtein + 'p';
			break;
		case 15:
			pseudoProtein = pseudoProtein + 's';
			break;
		case 16:
			pseudoProtein = pseudoProtein + 't';
			break;
		case 17:
			pseudoProtein = pseudoProtein + 'w';
			break;
		case 18:
			pseudoProtein = pseudoProtein + 'y';
			break;
		case 19:
			pseudoProtein = pseudoProtein + 'v';
			break;
		}



		for (int n = 0; n < 20; ++n)
		{
			countTheTimeOfAnAcid[n] = 0;
		}
		index = 0;
		count = 0;
	}
	return pseudoProtein;
}
