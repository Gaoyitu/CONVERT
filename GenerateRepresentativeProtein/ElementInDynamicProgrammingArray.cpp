#include "ElementInDynamicProgrammingArray.h"






ElementInDynamicProgrammingArray::ElementInDynamicProgrammingArray()
{
}

ElementInDynamicProgrammingArray::ElementInDynamicProgrammingArray(bool isGap, double value, string fromWhere)
{
	this->isGap = isGap;
	this->value = value;
	this->fromWhere = fromWhere;
}


bool ElementInDynamicProgrammingArray::getIsGap()
{
	return this->isGap;
}

double ElementInDynamicProgrammingArray::getValue()
{
	return this->value;
}

string ElementInDynamicProgrammingArray::getFromWhere()
{
	return this->fromWhere;
}

void ElementInDynamicProgrammingArray::setFormWhere(string direction)
{
	this->fromWhere = direction;
}

ElementInDynamicProgrammingArray::~ElementInDynamicProgrammingArray()
{
}
