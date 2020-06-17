#pragma once

template <class Type>
class nodeType
{
public:
	Type info;
	nodeType<Type> *link = nullptr;
	
};