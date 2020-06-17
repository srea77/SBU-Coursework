#pragma once
#include "nodeType.h"

using namespace std;

template<class Type>
class linkedListType
{
public:
	linkedListType();

	const linkedListType<Type>& operator=(const linkedListType<Type>&);
	void initializeList();
	bool isEmptyList() const;
	void print() const;
	int length() const;
	void destroyList();
	void insertAtEnd(const Type& value); // this also works, in essence, as addQ (for Queues)
										 // originally created to be called within the
										 // insertAtIndex() method when the last link is reached
	void insertAtIndex(const Type& value, int index);
	void deleteAtIndex(int index);
	void deleteValue(Type value);
	int search(const Type& newItem);
	void copyList(const linkedListType<Type>& otherList);

	void push(const Type& value);
	void pop();
	Type top();

	void addQ(const Type& value); // essentially an Alias for insertAtEnd() for clarity (Queue use)
	void deleteQ();
	Type front();

	nodeType<Type> *first;
	nodeType<Type> *last;
	mutable int count;

};



template <class Type>
linkedListType<Type>::linkedListType()
{
	first = nullptr;
	last = nullptr;
	count = 0;
}

template<class Type>
void linkedListType<Type>::initializeList()
{
	destroyList();
}


template<class Type>
bool linkedListType<Type>::isEmptyList() const
{
	if (first == nullptr)
	{
		return true;
	}
	else {
		return false;
	}
}


template<class Type>
void linkedListType<Type>::destroyList() {
	nodeType<Type> *temp;
	last->link = nullptr;

	while (first != nullptr)
	{
		temp = first;
		if (first->link != nullptr)
		{
			first = first->link;
			delete temp;
			temp = nullptr;
		}
		else
		{
			first = nullptr;
			delete first;
		}

	}
	last = nullptr;
	delete last;
	count = 0;
}


template<class Type>
void linkedListType<Type>::print() const
{
	if (this->isEmptyList()) {
		std::cout << "The list is empty." << std::endl;
	}

	nodeType<Type>* current = first;
	for(int i = 0; i <length(); i++){
		std::cout << current->info << " ";
		current = current->link;
	}
	std::cout << std::endl;

}

template<class Type>
int linkedListType<Type>::length() const
{
	if (first != nullptr)
	{
		nodeType<Type>* temp = first->link;
		count = 1;
		while (temp != first && temp != nullptr)
		{
			temp = temp->link;
			count++;
		}
	}

	else
	{
		count = 0;
	}
	
	return count;
}


template<class Type>
int linkedListType<Type>::search(const Type& newItem)
{
	nodeType<Type>* current = new nodeType<Type>;
	current = first;
	for (int i = 1; i <= this->length(); i++) {
		if (current->info == newItem) {
			return i;
		}
		else {
			current = current->link;
		}
	}
	return 0;
}


template<class Type>
void linkedListType<Type>::insertAtIndex(const Type& value, int index) {

	if (index == 1)
	{
		push(value);
		return;
	}

	nodeType<Type>* newNode = new nodeType<Type>;
	newNode->info = value;
	newNode->link = first;

	nodeType<Type>* current = new nodeType<Type>;
	nodeType<Type>* prevLink = new nodeType<Type>;

	if (first == nullptr) {
		first = newNode;
		last = newNode;
		count++;
	}

	current = first;


	if (index > 1 && index < length() + 1)
	{
		count++;
		for (int i = 1; i < index; i++)
		{
			if (current->link != nullptr)
			{
				prevLink = current;
				current = current->link;
			}
		}
		prevLink->link = newNode;
		newNode->link = current;
	}

	else if (index == length() + 1)
	{
		this->insertAtEnd(value); // ensures that once  the last node is created, it links back to the first node;

	}

	else
	{
		std::cout << "List index larger than list size!  Please expand the list or pick a smaller value." << std::endl;
	}
}


template <class Type>
void linkedListType<Type>::deleteAtIndex(int index)
{
	nodeType<Type>* current = first;

	if (index != 1 && length() > 1)
	{
		for (int i = 1; i < index - 1; i++)
		{
			current = current->link;
		}
		nodeType<Type>* temp = current->link;
		current->link = current->link->link;
		delete temp;
		std::cout << "Deleted." << std::endl;

	}

	else
	{
		if (length() > 1)
		{
			current = first;
			first = first->link;
			current = nullptr;
			delete current;
			last->link = first;
			count--;
		}

		else
		{
			first = nullptr;
			delete first;
		}

		std::cout << "Deleted." << std::endl;

	}
}


template <class Type>
void linkedListType<Type>::deleteValue(Type value)
{
	nodeType<Type>* current = first;
	nodeType<Type>* trailing = nullptr; // doubly linked-list implemented here, as trailing ptr is needed to link
										// the previous node to the next node (skipping over and then deleting the 
										// node whose correspond value is selected by the user
	bool found;	// boolean that denotes whether or not the value to be deleted has been found (used for traversing through list)


	if (first == nullptr) { // if the first element of the list has a memory location 
		std::cout << "List is empty; nothing to delete." << std::endl;
	}

	else
	{
		current = first; // current node begins at the beginning of the list
		found = false;  // item not yet found



		for (int i = 0; i < this->length(); i++) // until the "end" of the list is reached, i.e. the list is looped back to the beginning
		{
			if (first->info == value) {
				nodeType<Type>* temp = first;
				current = current->link;
				first = first->link;
				count--;

			}

			if (current->info == value) // if the current value equals the value to be deleted
			{
				found = true;
				break;
			}

			else
			{
				trailing = current;  // otherwise, the trailing node now becomes the current node
				current = current->link; // the current node becomes the next item in the list
			}

			if (current == nullptr) // if we reach a point where the list ends, found will still be false (as set by default),  
			{						// then the item must not be in the list
				std::cout << "Item is not in the list." << std::endl;
				break;
			}

			else // if current is not null
			{
				if (current->info == value) // if we've found our value
				{

					if (first == current)  // and if this value is found at the first node in the list
					{
						first = first->link; // the first element is now the second element
						if (first == nullptr) // if first is now null, 
						{
							last = nullptr;  // then it must be beyond the last element in the list, making last = first which equals nullptr
						}

						delete current; // delete the current node
					}

					else // otherwise, if we've found our value, but it is NOT the first element in the list,
					{
						trailing->link = current->link; //move trailing up by one
						if (current == last) { // if we've reached the last element in the list
							last = trailing;   // set last to the trailing node
							delete current;
						}
					}
					count--; // subtract from count, which denotes the total number of elements in the list
				}

			}
		}
	}

}


template<class Type>
void linkedListType<Type>::insertAtEnd(const Type& value) { // originally made as a special case for the insertAtIndex() function, but was also later used for stack method 'addQ()'
	nodeType<Type>* newNode = new nodeType<Type>;
	newNode->info = value;
	newNode->link = first;

	if (first == nullptr) 
	{
		first = newNode;
		last = newNode;
		count++;
	}

	else
	{
		last->link = newNode;
		last = newNode;
		count++;
	}

	last->link = first; //make the linked list circular
}



template<class Type>
void linkedListType<Type>::copyList(const linkedListType<Type>& otherList) {

	nodeType<Type> *newNode;
	nodeType<Type> *current;

	if (first != nullptr) {
		destroyList();
	}

	if (otherList.first == nullptr) {
		first = nullptr;
		last = nullptr;
		count = 0;
	}

	else
	{
		current = otherList.first;
		count = otherList.count;

		first = new nodeType<Type>;
		first->info = current->info;
		first = current;

		newNode = new nodeType<Type>;
		for (int i = 0; i < otherList.length(); i++)
		{
			newNode->info = current->info;
			newNode->link = current->link;
			current = current->link;


		}
	}
}


template<class Type>
const linkedListType<Type>& linkedListType<Type>::operator=(const linkedListType<Type>&otherList)
{
	if (this != &otherList)
	{
		copyList(otherList);  // simply use the previously created copyList to set the second list equal to the existing first list.
	}
	return *this;  // return this linked list
}


template<class Type>
void linkedListType<Type>::push(const Type& value)
{
	nodeType<Type>* newNode = new nodeType<Type>;

	if (first == nullptr) {
		newNode->info = value;
		newNode->link = first;
		first = newNode;
		count++;
	}

	else {
		newNode->info = value;
		newNode->link = first;

		nodeType<Type>* temp = new nodeType<Type>;
		temp = first;
		temp->info = first->info;

		first = newNode;
		first->info = newNode->info;

		newNode = temp;
		newNode->info = temp->info;
		count++;
		last->link = first;
	}

	if (last == nullptr) {
		last = newNode;
	}

}

template<class Type>
void linkedListType<Type>::pop()
{
	if (length() > 1)
	{
		nodeType<Type>* current = first;
		first = first->link;
		current = nullptr;
		delete current;
		last->link = first; // making the list circular so that instead of having to sift through every element in the list to reach the last element, 
							// which becomes increasingly timely as the list gets longer, we instead simply insert a pointer linking the first and last
							// element so that the time it takes to reach the last element is constant.
		count--;
	}

	else if (length() ==1)
	{
		destroyList();
		std::cout << "Destroyed." << std::endl;
	}

	else
	{
		std::cout << "List is empty." << std::endl;
	}
	
}

template<class Type>
Type linkedListType<Type>::top()
{
		return this->first->info;
}


template<class Type>
void linkedListType<Type>::addQ(const Type& value) // essentially just an alias for the same method to denote the use of the Queue
{ 
	this->insertAtEnd(value);
}

template<class Type>
void linkedListType<Type>::deleteQ()
{
	if (this->length() > 1)
	{
		deleteAtIndex(this->length());
	}

	else
	{
		deleteAtIndex(1);
	}
}

template<class Type>
Type linkedListType<Type>::front()
{
	if (this->isEmptyList())
	{
		std::cout << "The list is empty." << std::endl;
	}

	else
	{
		nodeType<Type>* current = first;
		for (int i = 0; i < length()-1; i++)
		{
			current = current->link;
		}
		return current->info;
	}


}
