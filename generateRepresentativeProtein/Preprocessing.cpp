#include "Preprocessing.h"



Preprocessing::Preprocessing()
{
}


Preprocessing::~Preprocessing()
{
}

int Preprocessing::getNumberOfProteins(string theNameOfFile)
{
	string str = "./" + theNameOfFile;
	ifstream ifs(str, ios::in);

	int count = 0;
	char c;
	while ((c = ifs.get()) != EOF)
	{
		if (c == '>')
		{
			count++;
		}
	}
	ifs.close();
	return count;
}

vector<string> Preprocessing::getSequencesInSameFamily(string familyId, string theNameOfInputFile)
{
	string thePathOfInputFile = "./" + theNameOfInputFile;
	//string thePathOfOutputFile = "./" + theNameOfOutputFile;
	ifstream ifs(thePathOfInputFile, ios::in);
	//ofstream ofs(thePathOfOutputFile, ios::out);

	string ALine;
	string familyIdOfOneProtein;
	vector<string> theSequencesInSameFamily;

	int firstBlankIndex;
	int secondBlankIndex;
	int guard = 0;
	while (!ifs.eof())
	{
		getline(ifs, ALine);
		if (ALine.at(0) == '>')
		{
			firstBlankIndex = ALine.find(' ');
			secondBlankIndex = ALine.find(' ', firstBlankIndex + 1);
			familyIdOfOneProtein = ALine.substr(firstBlankIndex + 1, secondBlankIndex - firstBlankIndex - 1);
			if (familyIdOfOneProtein.compare(familyId) == 0)
			{
				guard = 1;
			}
			if (familyIdOfOneProtein.compare(familyId) != 0)
			{
				guard = 0;
			}
		}
		if (ALine.at(0) != '>'&& guard == 1)
		{
			theSequencesInSameFamily.push_back(ALine);
			//ofs << ALine;
			//ofs << "\n";
		}
	}

	ifs.close();
	//ofs.close();
	return theSequencesInSameFamily;
}

vector<string> Preprocessing::getSequencesInSameFamilyUpdate(string familyId, string theNameOfInputFile)
{
	string thePathOfInputFile = "./" + theNameOfInputFile;
	//string thePathOfOutputFile = "./" + theNameOfOutputFile;
	ifstream ifs(thePathOfInputFile, ios::in);
	//ofstream ofs(thePathOfOutputFile, ios::out);

	string ALine;
	string familyIdOfOneProtein;
	vector<string> theSequencesInSameFamily;

	int firstBlankIndex;
	int secondBlankIndex;
	int guard = 0;
	int stopGuard = 0;
	while (!ifs.eof())
	{
		getline(ifs, ALine);
		if (ALine.at(0) == '>')
		{
			firstBlankIndex = ALine.find(' ');
			secondBlankIndex = ALine.find(' ', firstBlankIndex + 1);
			familyIdOfOneProtein = ALine.substr(firstBlankIndex + 1, secondBlankIndex - firstBlankIndex - 1);
			if (familyIdOfOneProtein.compare(familyId) == 0)
			{
				guard = 1;
				stopGuard = 1;
			}
			if (familyIdOfOneProtein.compare(familyId) != 0)
			{
				guard = 0;
			}

		}
		if (ALine.at(0) != '>'&& guard == 1)
		{
			theSequencesInSameFamily.push_back(ALine);
			//ofs << ALine;
			//ofs << "\n";
		}
		if (guard == 0 && stopGuard == 1)
		{
			break;
		}
	}

	ifs.close();
	//ofs.close();
	return theSequencesInSameFamily;
}

vector<string> Preprocessing::getSequencesInSameSuperFamily(string superFamilyId, string theNameOfInputFile)
{
	string thePathOfInputFile = "./" + theNameOfInputFile;
	vector<string> theSequencesInSameSuperFamily;
	//string thePathOfOutputFile = "./" + theNameOfOutputFile;
	ifstream ifs(thePathOfInputFile, ios::in);
	//ofstream ofs(thePathOfOutputFile, ios::out);

	string ALine;
	string superFamilyIdOfOneProtein;

	int firstBlankIndex;
	int thirdPotIndex;
	int guard = 0;
	while (!ifs.eof())
	{
		getline(ifs, ALine);
		if (ALine.at(0) == '>')
		{
			firstBlankIndex = ALine.find(' ');
			thirdPotIndex = ALine.find('.', (ALine.find('.', (ALine.find('.') + 1)) + 1));
			superFamilyIdOfOneProtein = ALine.substr(firstBlankIndex + 1, thirdPotIndex - firstBlankIndex - 1);
			if (superFamilyIdOfOneProtein.compare(superFamilyId) == 0)
			{
				guard = 1;
			}
			else
			{
				guard = 0;
			}
		}
		if (ALine.at(0) != '>'&& guard == 1)
		{
			theSequencesInSameSuperFamily.push_back(ALine);
			//ofs << ALine;
			//ofs << "\n";
		}
	}

	ifs.close();

	return theSequencesInSameSuperFamily;
	//ofs.close();
}

void Preprocessing::formatTheFile(string theNameOfInputFile, string theNameOfOutputFile)
{
	string thePathOfInputFile = "./" + theNameOfInputFile;
	string thePathOfOutputFile = "./" + theNameOfOutputFile;
	ifstream ifs(thePathOfInputFile, ios::in);
	ofstream ofs(thePathOfOutputFile, ios::out);

	char c[1024];
	int guard = 0;
	while (!ifs.eof())
	{
		ifs.getline(c, 1024, '\n');
		if (c[0] == '>')
		{
			if (guard == 1)
			{
				ofs << '\n';
			}
			ofs << c;
			ofs << '\n';

		}
		if (c[0] != '>')
		{
			ofs << c;
		}
		guard = 1;

	}

	ifs.close();
	ofs.close();
}

string Preprocessing::getTheFamliyIdOfProteinSequenceInFasta(string proteinSequenceInSCOPInformation)
{
	int firstBlankIndex = proteinSequenceInSCOPInformation.find(' ');
	int secondBlankIndex = proteinSequenceInSCOPInformation.find(' ', firstBlankIndex + 1);
	string famliyId = proteinSequenceInSCOPInformation.substr(firstBlankIndex + 1,
		secondBlankIndex - firstBlankIndex - 1);
	return famliyId;
}

vector<string> Preprocessing::getTheAllFamilyIdInData(string theNameOfFile)
{
	string thePathOfFile = "./" + theNameOfFile;
	ifstream ifs(thePathOfFile, ios::in);
	string str;
	vector<string> familyId;
	string familyIdOfOneProtein;
	int iteration = 0;
	int guard = 0;
	while (!ifs.eof())
	{
		getline(ifs, str);

		if (str.at(0) == '>')
		{
			familyIdOfOneProtein = getTheFamliyIdOfProteinSequenceInFasta(str);
			for (iteration; iteration < familyId.size(); iteration++)
			{
				if (familyId[iteration].compare(familyIdOfOneProtein) == 0)
				{
					guard = 1;
					break;
				}
			}
			if (guard == 0)
			{
				familyId.push_back(familyIdOfOneProtein);
			}
		}
		guard = 0;
	}
	iteration = 0;
	ifs.close();
	return familyId;
}

vector<string> Preprocessing::getAllSuperFamilyInData(string theNameOfFile)
{
	vector<string> familyId = getTheAllFamilyIdInData(theNameOfFile);
	vector<string> superFamily;
	int iteration = 0;
	int i = 0;
	int guard = 0;
	string str;
	for (iteration; iteration < familyId.size(); iteration++)
	{
		str = familyId[iteration].substr(0, familyId[iteration].find('.',
			(familyId[iteration].find('.', (familyId[iteration].find('.') + 1)) + 1)));
		for (i; i < superFamily.size(); i++)
		{
			if (superFamily[i].compare(str) == 0)
			{
				guard = 1;
				break;
			}
		}
		if (guard == 0)
		{
			superFamily.push_back(str);
		}
		guard = 0;
	}
	i = 0;
	iteration = 0;

	return superFamily;
}

string Preprocessing::getTheLengthestSequenceInSameFamily(vector<string> theSameFamily)
{
	string theLengthestSequence = theSameFamily[0];
	int i = 1;

	for (i; i < theSameFamily.size(); i++)
	{
		if (theSameFamily[i].length() > theLengthestSequence.length())
		{
			theLengthestSequence = theSameFamily[i];
		}
	}
	return theLengthestSequence;
}

string Preprocessing::getTheShortestSequenceInSameFamily(vector<string> theSameFamily)
{
	string theShortestSequence = theSameFamily[0];
	int i = 1;

	for (i; i < theSameFamily.size(); i++)
	{
		if (theSameFamily[i].length() < theShortestSequence.length())
		{
			theShortestSequence = theSameFamily[i];
		}
	}
	return theShortestSequence;
}

string Preprocessing::getTheSuitableAlignmentSequenceInSameFamily(vector<string> theSameFamily, int averageLength, int theLengthOfTheLengthestSequence, int theLengthOfTheShortestSequence)
{
	string suitableSequence = "";
	int differenceValue = (theLengthOfTheLengthestSequence - averageLength) / 4;
	if (theSameFamily.size() > 2)
	{
		for (int i = 0; i < theSameFamily.size(); ++i)
		{
			if (theSameFamily[i].size() > averageLength + differenceValue &&
				theSameFamily[i].size() < theLengthOfTheLengthestSequence - differenceValue)
			{
				suitableSequence = theSameFamily[i];
				break;
			}
		}
		if (suitableSequence == "")
		{
			for (int i = 0; i < theSameFamily.size(); ++i)
			{
				if (theSameFamily[i].size() > theLengthOfTheShortestSequence &&
					theSameFamily[i].size() < theLengthOfTheLengthestSequence)
				{
					suitableSequence = theSameFamily[i];
					break;
				}
			}
		}
		if (suitableSequence == "")
		{
			suitableSequence = theSameFamily[1];
		}
	}
	else
	{
		suitableSequence = theSameFamily[0];
	}

	return suitableSequence;
}

int Preprocessing::getTheAverageLengthOfSameFamilySequence(vector<string> theSameFamily)
{
	int averageLength = 0;
	for (int i = 0; i < theSameFamily.size(); ++i)
	{
		averageLength += theSameFamily[i].size();
	}
	averageLength = averageLength / theSameFamily.size();
	return averageLength;
}

