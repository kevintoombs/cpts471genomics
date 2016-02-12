// Program1-Toombs.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;

struct config 
{
	int matchScore;
	int mismatchScore;
	int startGapScore;
	int continueGapScore;
};

struct DP_cell 
{
	int score;
};

bool parseFasta(char *argv[], string&, string&);
config getConfig(char *argv[]);

int main(int argc, char *argv[])
{
	string sequence1, sequence2;
	config c;

	if (!parseFasta(argv, sequence1, sequence2)) 
	{
		return 1;
	}

	c = getConfig(argv);
	cout << endl << sequence1 << endl;
	cout << endl << sequence2 << endl;

	
	
	DP_cell DPTable[11][11];

    return 0;
}

bool parseFasta(char *argv[], string &s1, string &s2) 
{
	ifstream fasta;
	if (argv[1]) 
	{
		fasta.open(argv[1], ios::in);
		if (!fasta.good()) 
		{
			cout << endl << argv[1] << " not found." << endl;
			return false;
		}
		cout << endl << argv[1] << " opened." << endl;
		
		char ch;
		int state = 0;

		//TODO:: should turn the state functionality into a recursive function that returns an array of string
		while (fasta.get(ch)) 
		{
			if (ch == '\n' && state == 0)
			{
				fasta.get(ch);
				state = 1;
				cout << endl <<  "String 1 start: " << endl;
			}
			if (state == 1) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					s1 += toupper(ch);
				}
				cout << ch;
				if (ch == '>') 
				{
					state = 2;
					cout << endl << "String 1 end. " << endl;
				}
			}
			if (ch == '\n' && state == 2) 
			{
				fasta.get(ch);
				state = 3;
				cout << endl << "String 2 start: " << endl;
			}
			if (state == 3) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					s2 += toupper(ch);
				}
				cout << ch;
				if (ch == '>') {
					state = 4;
				}
			}
		}
		if (state == 3) 
		{
			cout << endl << "String 2 end. " << endl;
		}
		else 
		{
			cout << endl << "Invalid fasta input." << endl;
		}
		fasta.close();
		return true;
	} 
	else 
	{
		cout << endl << "No input file specified. Correct usage is  $ <executable name> "
			<< "<input file containing both s1 and s2> <0: global, 1: local> <optional: path to parameters config file>"
			<< endl;
		return false;
	}
}

config getConfig(char *argv[]) 
{
	ifstream conFile;
	string filename;

	if (argv[3])
	{
		filename = argv[3];
	}
	else
	{
		filename = "parameters.config";
	}

	config c;
	return c;
}

