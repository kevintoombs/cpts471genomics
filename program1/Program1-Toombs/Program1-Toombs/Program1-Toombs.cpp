// Program1-Toombs.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

using namespace std;

struct config 
{
	int matchScore = 0;
	int mismatchScore = 0;
	int startGapScore = 0; //h
	int continueGapScore = 0; //g
};

struct DP_cell 
{
	int S = 0;// numeric_limits<int>::min();
	int D = 0;// numeric_limits<int>::min();
	int I = 0;// numeric_limits<int>::min();
};

struct DP_table
{
	string sequence1, sequence2;
	vector< vector<DP_cell> > t;
	config c;
	int alightmentType;
	tuple<int, int> maxPair;
};

bool parseFasta(char *argv[], string&, string&);
config getConfig(int argc, char *argv[]);
int getAlignmentType(int argc, char *argv[]);
void buildTable(DP_table &t);
void initTable(DP_table &t);
void calcTable(DP_table &t);
int maximum(int S, int D, int I, int alignmentType);
int subFunction(char a, char b, config c);
void printTable(DP_table &t);
void retrace(DP_table &t);
int direction(DP_cell);
int cellMax(DP_cell c, int alignmentType);

int main(int argc, char *argv[])
{
	DP_table t;

	string sequence1, sequence2;
	
	int aligntmentType = getAlignmentType(argc, argv);
	config c = getConfig(argc, argv); //cout << c.matchScore << c.mismatchScore << c.startGapScore << c.startGapScore;

	if (!parseFasta(argv, sequence1, sequence2)) 
	{
		return 1;
	} cout << endl << sequence1 << endl << "123423456" << endl << sequence2 << endl;
	
	t.sequence1 = sequence1;
	t.sequence2 = sequence2;
	t.c = c;
	t.alightmentType = aligntmentType;

	buildTable(t); // cout << t.t.size() << " " << t.t[3].size() << endl;
	calcTable(t);
	retrace(t);

    return 0;
}

void buildTable(DP_table &t)
{
	cout << endl << "Bulding table." << endl;
	DP_cell c;
	t.t.resize(t.sequence1.length()+1, vector<DP_cell>(t.sequence2.length()+1, c) );


	cout << "Done building." << endl;
	return;
}

void calcTable(DP_table &t)
{
	initTable(t);
	cout << endl << "Calculating Table." << endl;
	int i, j, maxValue = 0;
	tuple<int, int> maxPair;
	for (i = 1; i <= t.sequence1.length(); i++)
	{
		for (j = 1; j <= t.sequence2.length(); j++)
		{
			t.t[i][j].S = maximum
				(t.t[i - 1][j - 1].S + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.t[i - 1][j - 1].D + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.t[i - 1][j - 1].I + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.alightmentType);
			t.t[i][j].D = maximum
				(t.t[i - 1][j].S + t.c.startGapScore + t.c.continueGapScore,
					t.t[i - 1][j].D + t.c.continueGapScore,
					t.t[i - 1][j].I + t.c.startGapScore + t.c.continueGapScore,
					t.alightmentType);
			t.t[i][j].I = maximum
				(t.t[i][j - 1].S + t.c.startGapScore + t.c.continueGapScore,
					t.t[i][j - 1].D + t.c.continueGapScore + t.c.continueGapScore,
					t.t[i][j - 1].I + t.c.startGapScore,
					t.alightmentType);
			if (maximum(t.t[i][j].S, t.t[i][j].D, t.t[i][j].I, t.alightmentType) > maxValue)
			{
				maxValue = maximum(t.t[i][j].S, t.t[i][j].D, t.t[i][j].I, t.alightmentType);
				maxPair = make_tuple(i, j);
				t.maxPair = maxPair;
			}
		}
	} i--; j--;// cout << "Cell: " << i << "," << j << " max: " << maximum(t.t[i][j].S, t.t[i][j].D, t.t[i][j].I, t.alightmentType) << endl;
	if (t.alightmentType == 0)
	{
		cout << " maximum global allignment: " << maximum(t.t[i][j].S, t.t[i][j].D, t.t[i][j].I, t.alightmentType) << endl;
	}
	else
	{
		cout << " maximum local allignment: " << maxValue << endl;
	}
	//printTable(t);
	
}

void initTable(DP_table &t)
{
	cout << endl << "Initializing Table." << endl;
	int h = t.c.startGapScore;
	int g = t.c.continueGapScore;

	t.t[0][0].S = 0;

	t.t[1][0].D = h;
	t.t[1][0].S = -1147483648;
	t.t[1][0].I = -1147483648;

	t.t[0][1].I = h;
	t.t[0][1].S = -1147483648;
	t.t[0][1].D = -1147483648;

	for (int i = 2; i < t.sequence1.length(); i++) 
	{
		t.t[i][0].D = t.t[i-1][0].D + g;
		t.t[i][0].S = -1147483648;
		t.t[i][0].I = -1147483648;
	}
	for (int j = 2; j < t.sequence2.length(); j++) 
	{
		t.t[0][j].I = t.t[0][j - 1].I + g;
		t.t[0][j].S = -1147483648;
		t.t[0][j].D = -1147483648;
	} 

	cout << "Done initializing." << endl;
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
				cout << "Sequence 1 start. ";
			}
			if (state == 1) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					s1 += toupper(ch);
				}
				//cout << ch;
				if (ch == '>') 
				{
					state = 2;
					cout << "Sequence 1 end.";
				}
			}
			if (ch == '\n' && state == 2) 
			{
				fasta.get(ch);
				state = 3;
				cout << endl << "Sequence 2 start. ";
			}
			if (state == 3) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					s2 += toupper(ch);
				}
				//cout << ch;
				if (ch == '>') {
					state = 4;
				}
			}
		}
		if (state == 3) 
		{
			cout << "Sequence 2 end." << endl;
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

config getConfig(int argc, char *argv[]) 
{
	ifstream conFile;
	string filename;
	config c;

	if (argc == 4)
	{
		cout << argv[3];
		filename = argv[3];
		cout << endl << filename << " opened." << endl;
	}
	else
	{
		filename = "parameters.config";
		cout << "Default file: " << filename << " opened." << endl;
	}

	conFile.open(filename, ios::in);
	if (!conFile.good())
	{
		cout << endl << filename << " not found. All paremters default to 0." << endl;
		return c;
	}

	string line;
	while (getline(conFile, line))
	{
		stringstream s(line);
		string tmp1;
		string tmp2;
		while (!s.eof()) {
			s >> tmp1;
			s >> tmp2;
		} //cout << tmp1 << tmp2 << endl;
		if (tmp1 == "match")
		{
			c.matchScore = stoi(tmp2);
			cout << tmp1 << " = " << tmp2 << endl;
		}
		if (tmp1 == "mismatch")
		{
			c.mismatchScore = stoi(tmp2);
			cout << tmp1 << " = " << tmp2 << endl;
		}
		if (tmp1 == "h")
		{
			c.startGapScore = stoi(tmp2);
			cout << tmp1 << " = " << tmp2 << endl;
		}
		if (tmp1 == "g")
		{
			c.continueGapScore = stoi(tmp2);
			cout << tmp1 << " = " << tmp2 << endl;
		}
	}
	cout << "Parameters read." << endl;
	
	return c;
}

int getAlignmentType(int argc, char *argv[])
{
	if (argc >= 3)
	{
		int i = stoi(argv[2]);
		cout << "Alignment type (0:global, 1:local): " << i << endl << endl;
		return i;
	}
	else
	{
		cout << endl << "Correct usage is  $ <executable name> <input file containing both s1 and s2> "
			<< "<0: global, 1: local> <optional: path to parameters config file>"
			<< endl;
	}

	return 0;
}

int cellMax(DP_cell c, int alignmentType)
{
	int max = c.S;
	if (c.D > max)
	{
		max = c.D;
	}
	if (c.I > max)
	{
		max = c.I;
	}
	if (alignmentType == 1 && max < 0) max = 0;
	return max;
}

int maximum(int S, int D, int I, int alignmentType)
{
	int max = S;
	if (D > max) 
	{ 
		max = D;
	}
	if (I > max)
	{
		max = I;
	}
	if (alignmentType == 1 && max < 0) max = 0;
	return max;
}

int subFunction(char a, char b, config c)
{
	if (a == b)
	{
		//cout << "Returning: " << c.matchScore << endl;
		return c.matchScore;
	}
	else
	{
		//cout << "Returning: " << c.mismatchScore << endl;
		return c.mismatchScore;
	}
}

void printTable(DP_table &t)
{
	for (int i = 0; i < 100; i++)
	{
		for (int j = 0; j < 80; j++)
		{
			int m = maximum(t.t[i][j].S, t.t[i][j].D, t.t[i][j].I,t.alightmentType);
			if (m >= 0) printf("%3d", m);
			else cout << "   ";
		}
		cout << endl;
	}
}

void retrace(DP_table &t)
{
	cout << endl << "retracing." << endl;
	int matches = 0, mismatches = 0, gaps = 0, openingGaps = 0;
	int lastDir = 0;
	int lastValue = 1;
	int i = t.sequence1.length();
	int j = t.sequence2.length();

	if (t.alightmentType == 0)
	{
		while (i >= 0 || j >= 0)
		{
			cout << cellMax(t.t[i][j], t.alightmentType) << endl;
			int dir = direction(t.t[i][j]);
			if (dir == 1)
			{
				j--;
				if (lastDir == dir)
				{
					gaps++;
				}
				else
				{
					openingGaps++;
					gaps++;
				}
			}
			if (dir == 3)
			{
				i--;
				if (lastDir == dir)
				{
					gaps++;
				}
				else
				{
					openingGaps++;
					gaps++;
				}
			}

			if (dir == 2)
			{
				j--;
				i--;
				if (i >= 0 && j >= 0)
				{
					if (subFunction(t.sequence1[i], t.sequence2[j], t.c) > 0)
					{
						matches++;
					}
					else
					{
						mismatches++;
					}
				}
			}
			lastDir = dir;
		}
	}

	if (t.alightmentType == 1)
	{
		while (lastValue != 0)
		{
			lastValue = cellMax(t.t[i][j], t.alightmentType);
			cout << lastValue << endl;
			
			int dir = direction(t.t[i][j]);
			if (dir == 1)
			{
				j--;
				if (lastDir == dir)
				{
					gaps++;
				}
				else
				{
					openingGaps++;
					gaps++;
				}
			}
			if (dir == 3)
			{
				i--;
				if (lastDir == dir)
				{
					gaps++;
				}
				else
				{
					openingGaps++;
					gaps++;
				}
			}

			if (dir == 2)
			{
				j--;
				i--;
				if (i >= 0 && j >= 0)
				{
					if (subFunction(t.sequence1[i], t.sequence2[j], t.c) > 0)
					{
						matches++;
					}
					else
					{
						mismatches++;
					}
				}
			}
			lastDir = dir;
		}
	}


	cout << "Matches: " << matches << endl;
	cout << "Mismatches: " << mismatches << endl;
	cout << "Opening Gaps: " << openingGaps<< endl;
	cout << "Continuing gaps: " << gaps << endl;
	cout << endl << "done retracing." << endl;
	int calcedScore = matches*t.c.matchScore + mismatches*t.c.mismatchScore + openingGaps*t.c.startGapScore + gaps*t.c.continueGapScore;
	printf("%i*%i+%i*%i+%i*%i+%i*%i=%i", matches, t.c.matchScore, mismatches, t.c.mismatchScore, openingGaps, t.c.startGapScore, gaps, t.c.continueGapScore, calcedScore);
}

int direction(DP_cell c) 
{
	int dir = 2;
	int max = c.S;

	if (c.D > max)
	{
		max = c.D;
		dir = 3;
	}
	if (c.I > max)
	{
		max = c.I;
		dir = 1;
	}
	return dir;
}