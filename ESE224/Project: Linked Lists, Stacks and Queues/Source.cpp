#include <iostream>
#include "linkedListType.h"
using namespace std;


void someCodeDetails() {
	cout << "A few things to note: For the Search() function in the linkedListType .cpp, since the assignment \n" <<
		"does not specify which integer to return if the specified value does not exist in the list, I have \n" <<
		"decided to reserve '0' for the list element to be returned if the specified value does not exist, \n" <<
		"meaning that my list elements start at index 1, not 0.\n\n";
	cout << "With respect to double-link and circular implementation, comments are provided in the linkedListType \n" <<
		"class file to denote when these components are implemented.\n\n";

}

void menuOptions() {
	cout << "Enter 2 to determine if list 1 (L1) is empty." << endl;
	cout << "Enter 3 to destroy L1" << endl;
	cout << "Enter 4 to print L1" << endl;
	cout << "Enter 5 to find the length of L1" << endl;
	cout << "Enter 6 to search the list for a given item (value)" << endl;
	cout << "Enter 7 to insert an item with a given item/integer value at the "
		 << "specified position" << endl;
	cout << "Enter 8 to add an element to the end of the list" << endl;
	cout << "Enter 9 to delete an item at a specified position in the list" << endl;
	cout << "Enter 10 to delete the first occurence of an item with the entered value"
		 << endl;
	cout << "Enter 11 to copy L1 to L2" << endl;
	cout << "Enter 12 to add an element to the beginning of the list " << endl;
	cout << "Enter 13 to remove the first element from L1 (pop)" << endl;
	cout << "Enter 14 to view the first element from L1 (top)" << endl;
	cout << "Enter 15 to insert an item at the end of the list (addQ())" << endl;
	cout << "Enter 16 to delete the last element from the list (deleteQ())" << endl;
	cout << "Enter 17 to print the last element of the list (front())" << endl;
	cout << endl;
	cout << "Enter 0 to quit.";
	cout << endl << endl;
}

int main() {
	
	someCodeDetails();
	cout << "\n\n";
	menuOptions();

	linkedListType<int>* L1 = new linkedListType<int>;
	linkedListType<int>* L2 = new linkedListType<int>;
	int choice;
	bool empty;
	choice = 1;
	int menuCount = 0; // used to repeat the menuOptions() method after a certain number of inputs are put in, 
					   // so that user does not have to scroll up to read options

	do {
		if (menuCount > 5)
		{
			cout << "Re-listing options: " << endl;
			menuOptions();
			menuCount = 0;
		}

		menuCount++;
		cout << "User input: ";
		cin >> choice;
		switch (choice) {
			int index, value;

		case 0:
			return 0;

		case 2:
			empty = L1->isEmptyList();
			if (empty)
			{
				cout << "L1 is empty." << endl;
			}
			else
			{
				cout << "L1 is NOT empty." << endl;
			}
			break;

		case 3:
			if (L1->length() > 0)
			{
				L1->destroyList();
			}
			cout << "Destroyed.";
			break;

		case 4:
			cout << "L1: ";
			L1->print();
			break;

		case 5:
			cout << "List length = " << L1->length();
			break;

		case 6:
			cout << "Enter the value to be found in the list (if it exists): ";
			cin >> value;
			cout << endl;
			index = L1->search(value);
			if (index != 0)
			{
				cout << value << " found at " << index;
			}
			else
			{
				cout << value << " is not in the list.";
			}
			cout << endl;
			break;

		case 7:
			cout << "Enter the value to be inserted and its index (i.e, 1 for first element, 2 for second, etc.)";
			cin >> value >> index;
			L1->insertAtIndex(value, index);
			break;

		case 8:
			cout << "Enter the value: ";
			cin >> value;
			L1->insertAtEnd(value);
			break;

		case 9:
			cout << "Enter the position whose element will be deleted: ";
			cin >> index;
			if (index > 0 && index <= L1->length())
			{
				L1->deleteAtIndex(index);
			}
			else
			{
				cout << "Please enter a valid index (between 1 and the list length; you "
					 << "can enter '5' to display the length)" << endl;
			}
			
			break;

		case 10:
			cout << "Enter the value to delete: ";
			cin >> value;
			L1->deleteValue(value);
			break;

		case 11:
			L2 = L1;
			cout << "Copied successfuly." << endl;
			cout << "L2: ";
			L2->print();
			break;

		case 12:
			cout << "Enter the value: ";
			cin >> value;
			L1->push(value);
			break;

		case 13:
			L1->pop();
			break;

		case 14:
			if (L1->length() > 0)
			{
				cout << "First element of L1: " << L1->top();
			}

			else
			{
				cout << "Nothing to return (list is empty)." << endl;
			}
			break;

		case 15:
			cout << "Enter the value to be inserted: ";
			cin >> value;
			L1->addQ(value);
			break;

		case 16: 
			L1->deleteQ();
			break;

		case 17: 
			cout << "Last element = " << L1->front();
			break;

		default:
			cout << "Please pick an integer between 2 and 17";
			break;

		}
		cout << endl;

	} while (choice != 0);

}