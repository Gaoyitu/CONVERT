#pragma once
#include<string>
using std::string;
class ElementInDynamicProgrammingArray
{
public:
	ElementInDynamicProgrammingArray();
	ElementInDynamicProgrammingArray(bool isGap, double value, string fromWhere);
	bool getIsGap();
	double getValue();
	string getFromWhere();
	void setFormWhere(string direction);
	~ElementInDynamicProgrammingArray();
private:
	bool isGap;
	double value;
	string fromWhere;
};


