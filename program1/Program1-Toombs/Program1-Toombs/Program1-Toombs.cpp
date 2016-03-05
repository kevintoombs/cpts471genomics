// Program1-Toombs.cpp : Defines the entry point for the console application.
// 6695 global/local 8475

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
	int S = -1147483648;// numeric_limits<int>::min();
	int D = -1147483648;// numeric_limits<int>::min();
	int I = -1147483648;// numeric_limits<int>::min();
};

struct DP_table
{
	string sequence1, sequence2;
	string id1, id2;
	vector< vector<DP_cell> > t;
	config c;
	int alightmentType;
	tuple<int, int> maxPair;
};

bool parseFasta(char *argv[], DP_table &t);
config getConfig(int argc, char *argv[]);
int getAlignmentType(int argc, char *argv[]);
int main(int argc, char * argv[]);
void buildTable(DP_table &t);
void initTable(DP_table &t);
void calcTable(DP_table &t);
int maximum(int S, int D, int I, int alignmentType);
int subFunction(char a, char b, config c);
void printTable(DP_table &t);
void retrace(DP_table &t);
void recursivelyPrintChildren(DP_table &t, int i, int j);
int direction(DP_table &t, int i, int j);
int maximum2(int S, int D, int I, int alignmentType, DP_cell &c);
int cellMax(DP_cell c);
void testDirection(int lastValue, DP_cell to, config c, int dir, int i, int j, DP_table &t);

int main(int argc, char *argv[])
{

	DP_table t;
	t.alightmentType = getAlignmentType(argc, argv);
	t.c = getConfig(argc, argv); 
	if (!parseFasta(argv, t) )
	{
		cin.ignore();
		return 1;
	} 

	printf("Scores: Match = %i, Mismatch = %i, h = %i, g = %i\n", t.c.matchScore, t.c.mismatchScore, t.c.startGapScore, t.c.continueGapScore);
	cout << endl;

	cout << "Sequence 1 = \"" << t.id1;
	printf("\", length = %i characters\n", (int)t.sequence1.length());
	cout << "Sequence 2 = \"" << t.id2;
	printf("\", length = %i characters\n", (int)t.sequence2.length());
	cout << endl;
	
	buildTable(t);
	calcTable(t);
	
	retrace(t);
	
	//printTable(t);

	cin.ignore();

    return 0;
}

void buildTable(DP_table &t)
{
	printf("Building Table.");
	DP_cell c;
	t.t.resize(t.sequence1.length()+1, vector<DP_cell>(t.sequence2.length()+1, c) );
	cout << endl;
	initTable(t);

	return;
}

void calcTable(DP_table &t)
{
	
	int maxValue = 0;
	tuple<int, int> maxPair;

	printf("Calculating Table[rows remaining]:");
	for (size_t i = 1; i <= t.sequence1.length(); i++)
	{
		if (i % ((t.sequence1.length() / 100) + 1) == 0) printf("[%i]", (int)t.sequence1.length() - (int)i);
		for (size_t j = 1; j <=  t.sequence2.length(); j++)
		{
			int sSub = subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c);
			DP_cell subCell = t.t[i - 1][j - 1];
			DP_cell deleteCell = t.t[i - 1][j];
			DP_cell insertCell = t.t[i][j - 1];
			

			t.t[i][j].S = maximum (subCell.S + sSub, subCell.D + sSub, subCell.I + sSub, t.alightmentType);
			t.t[i][j].D = maximum
				(deleteCell.S + t.c.startGapScore + t.c.continueGapScore,
					deleteCell.D + t.c.continueGapScore,
					deleteCell.I + t.c.startGapScore + t.c.continueGapScore,
					t.alightmentType);
			t.t[i][j].I = maximum
				(insertCell.S + t.c.startGapScore + t.c.continueGapScore,
					insertCell.D + t.c.continueGapScore + t.c.continueGapScore,
					insertCell.I + t.c.continueGapScore,
					t.alightmentType);

			if (t.alightmentType == 1)
			{
				int thisMax = cellMax(t.t[i][j]);
				if ( thisMax > maxValue)
				{
					maxValue = thisMax;
					maxPair = make_tuple(i, j);
					t.maxPair = maxPair;
				}
			}
		}
	} 
	cout << endl << endl;
	
}

void initTable(DP_table &t)
{
	int h = t.c.startGapScore;
	int g = t.c.continueGapScore;

	t.t[0][0].S = 0;
	t.t[1][0].D = -1147483648;
	t.t[1][0].I = -1147483648;

	for (size_t i = 1; i <= t.sequence1.length(); i++) 
	{
		t.t[i][0].S = -1147483648;
		t.t[i][0].D = h + (int)i * g;
		t.t[i][0].I = -1147483648;
	}
	for (size_t j = 1; j <= t.sequence2.length(); j++) 
	{
		t.t[0][j].S = -1147483648;
		t.t[0][j].D = -1147483648;
		t.t[0][j].I = h + (int)j * g;
	} 
}

bool parseFasta(char *argv[], DP_table &t)
{
	ifstream fasta;
	if (argv[1]) 
	{
		fasta.open(argv[1], ios::in);
		if (!fasta.good()) 
		{
			return false;
		}
		
		char ch;
		int state = 0;
		int getId = 0;

		//probably overly complicated state based method to parse the fasta file
		while (fasta.get(ch)) 
		{

			//should always evaluate true at first 
			if (state == 0 && getId == 1)
			{
				if (isspace(ch))
				{
					//cout << "not alnum" << endl;
					getId = 0;
				}
				else if (!isspace(ch))
				{
					t.id1 += ch;
					//cout << t.id1 << endl;
				}
			}

			if (ch == '>' && state == 0 && getId == 0 )
			{
				getId = 1;
			}

			

			//after the first line
			if (ch == '\n' && state == 0)
			{
				state = 1;
			}

			if (state == 2 && getId == 1)
			{
				if (isspace(ch))
				{
					//cout << "not alnum" << endl;
					getId = 0;
				}
				else if (!isspace(ch))
				{
					t.id2 += ch;
					//cout << t.id2 << endl;
				}
			}
			
			if (state == 1) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					t.sequence1 += toupper(ch);
					//cout << t.sequence1 << endl;
					//cout << t.sequence1.length() << endl;
				}
				if (ch == '>') 
				{
					state = 2;
					getId = 1;
				}
			}

			


			if (ch == '\n' && state == 2) 
			{
				state = 3;
			}

			if (state == 3) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					t.sequence2 += toupper(ch);
					//cout << t.sequence2 << endl;
					//cout << t.sequence2.length() << endl;
				}
				//cout << ch;
				if (ch == '>') {
					state = 4;
				}
			}

			
		}

		//cout << t.id1 << ',' << t.id2 << endl;
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
		}
		if (tmp1 == "mismatch")
		{
			c.mismatchScore = stoi(tmp2);
		}
		if (tmp1 == "h")
		{
			c.startGapScore = stoi(tmp2);
		}
		if (tmp1 == "g")
		{
			c.continueGapScore = stoi(tmp2);
		}
	}
	
	return c;
}

int getAlignmentType(int argc, char *argv[])
{
	if (argc >= 3)
	{
		int i = stoi(argv[2]);
		//cout << "Alignment type (0:global, 1:local): " << i << endl << endl;
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

int cellMax(DP_cell c)
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

/*
int maximum2(int S, int D, int I, int alignmentType, DP_cell &c)
{
	int max = S;
	c.dir = 2;
	if (D > max)
	{
		max = D;
		c.dir = 3;
	}
	if (I > max)
	{
		max = I;
		c.dir = 1;
	}
	if (alignmentType == 1 && max < 0) max = 0;
	return max;
}
*/

int subFunction(char a, char b, config c)
{
	if (a == b)
	{
		return c.matchScore;
	}
	else
	{
		return c.mismatchScore;
	}
}

void printTable(DP_table &t)
{
	int m = 0;
	int a = t.alightmentType;
	int s1length = (int)t.sequence1.length();
	int s2length = (int)t.sequence2.length();

	for (int i = 0; i <= 20/*s1length*/; i++) //do min(100, slength)
	{
		for (int j = 0; j <= 20/*s2length*/; j++) //do min(100, slength)
		{
			DP_cell c = t.t[i][j];
			m = cellMax(c);
			if (/*m >= 0*/true) printf("%3d", m);
			else cout << "   ";
		}
		cout << endl;
	}
}

void retrace(DP_table &t)
{
	//cout << "retrace\n";
	int matches = 0, mismatches = 0, gaps = 0, openingGaps = 0;
	int lastDir = -1;
	int dir = -1;
	stack<char> s1, s2, r;
	int lastValue = 1;
	int i = (int)t.sequence1.length();
	int j = (int)t.sequence2.length();

	if (t.alightmentType == 1)
	{
		i = get<0>(t.maxPair);
		j = get<1>(t.maxPair);
	}



	while (i > 0 || j > 0)
	{
		if (t.alightmentType == 1 && cellMax(t.t[i][j]) == 0)
		{
			DP_cell debug = t.t[i][j];
			DP_cell sub = t.t[i - 1][j - 1];
			DP_cell del = t.t[i - 1][j];
			DP_cell in = t.t[i][j - 1];

			break;
		}
		lastDir = dir;
		dir = direction(t, i, j);
		if (dir == 1)
		{
			j--;
			s1.push('-');
			s2.push(t.sequence2[j]);
			r.push(' ');
			if (lastDir == dir)
			{
				cout << "in" << endl;
				gaps++;
			}
			else
			{
				openingGaps++;
				gaps++;
				//cout << "start gap" << endl;
			}
		}
		if (dir == 3)
		{
			i--;
			s2.push('-');
			r.push(' ');
			s1.push(t.sequence1[i]);
			if (lastDir == dir)
			{
				gaps++;
				cout << "del" << endl;
			}
			else
			{
				openingGaps++;
				gaps++;
				//cout << "start gap" << endl;
			}
		}

		if (dir == 2)
		{
			j--;
			i--;
				
			if (i >= 0 && j >= 0)
			{
				s1.push(t.sequence1[i]);
				s2.push(t.sequence2[j]);
				if (subFunction(t.sequence1[i], t.sequence2[j], t.c) > 0)
				{
					matches++;
					r.push('|');
					//cout << "match" << endl;
				}
				else
				{
					mismatches++;
					r.push(' ');
					//cout << "mismatch" << endl;
				}
			}
		}	
	}

	cout << endl << endl;
	
	while (!s1.empty())
	{
		int to60 = 0;

		// printf("s1 %5d ", );
		while (to60 != 60)
		{
			if (!s1.empty()) cout << s1.top();
			if (!s1.empty()) s1.pop();
			to60++;
		} 
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!r.empty()) cout <<  r.top();
			if (!r.empty()) r.pop();
			to60++;
		} 
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!s2.empty()) cout << s2.top();
			if (!s2.empty()) s2.pop();
			to60++;
		}
		to60 = 0;
		cout << endl << endl;
	}
	cout << endl;


	cout << "Report:\n\n";
	if (t.alightmentType == 0)/**/ printf("Global optimal score = %i.\n", cellMax(t.t[t.sequence1.size()][t.sequence2.size()]));
	if (t.alightmentType == 1)/**/ printf("Local optimal score = %i found at %i,%i.\n", cellMax(t.t[get<0>(t.maxPair)][get<1>(t.maxPair)]), get<0>(t.maxPair), get<1>(t.maxPair));
	printf("Number of:  matches = %i, mismatches = %i, gaps = %i, opening gaps = %i\n", matches, mismatches, gaps, openingGaps);
	float total = (float)gaps + (float)matches + (float)mismatches;
	float identities = (float)matches + (float)mismatches;
	float p1 = 100 * identities / total;
	float p2 = 100 * gaps / total;
	printf("Identities = %.0f/%.0f (%2.1f percent), Gaps = %i/%.0f (%2.1f percent)\n", identities, total, p1, gaps, total, p2);
	printf("Sanity Check: %i = ", matches * t.c.matchScore + mismatches * t.c.mismatchScore + gaps * t.c.continueGapScore + openingGaps * t.c.startGapScore);
	printf("%i * %i + %i * %i + %i * %i + %i * %i\n\n", matches , t.c.matchScore , mismatches , t.c.mismatchScore , gaps ,t.c.continueGapScore , openingGaps , t.c.startGapScore);
}

int direction(DP_table &t, int i, int j)
{
	//2 = S, 3 = D, 1 = I 
	if (j == 0) return 3;
	if (i == 0) return 1;

	DP_cell c = t.t[i][j];
	DP_cell fromCell; 
	int max = cellMax(c);

	int dir = 0;

	if (i != 0 && j != 0)
	{
		int sSub = subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c);
		DP_cell subCell = t.t[i - 1][j - 1];
		int testValue1 = cellMax(subCell) + sSub;

		if (max == testValue1)
		{
			fromCell = subCell;
			dir = 2;
		}
	}
	if (i != 0)
	{
		DP_cell deleteCell = t.t[i - 1][j];
		int testValue1 = deleteCell.S + t.c.startGapScore + t.c.continueGapScore;
		int testValue2 = deleteCell.D + t.c.continueGapScore;
		int testValue3 = deleteCell.I + t.c.startGapScore + t.c.continueGapScore;
		
		if (max == deleteCell.S + t.c.startGapScore + t.c.continueGapScore ||
			max == deleteCell.D + t.c.continueGapScore ||
			max == deleteCell.I + t.c.startGapScore + t.c.continueGapScore)
		/*if (max == testValue1 || testValue2 || testValue3)*/
		{
			fromCell = deleteCell;
			dir = 3;
		}
	}
	if (j != 0)
	{
		DP_cell insertCell = t.t[i][j - 1];
		int testValue1 = insertCell.S + t.c.startGapScore + t.c.continueGapScore;
		int testValue2 = insertCell.D + t.c.startGapScore + t.c.continueGapScore;
		int testValue3 = insertCell.I + t.c.continueGapScore;

		if (max == testValue1 || max == testValue2 || max == testValue3)
		{
			fromCell = insertCell;
			dir = 1;
		}
	}
	
	if (dir == 0) 
	{
		exit(2);
	}

	testDirection(max, c, t.c, dir, i, j, t);
	return dir;
}

//a + b = c
void testDirection(int lastValue, DP_cell to, config c, int dir, int i, int j, DP_table &t)
{
	bool y = true;
	char a = t.sequence1[i];
	char b = t.sequence2[j];

	if (dir == 1) 
	{
		int sSub = subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c);
		y = cellMax(to) == lastValue;
	}
	else if (dir == 2)
	{
		y = true;

	}
	else if (dir == 3)
	{
		y = true;
	}
	if (!y)
	{
		cout << "[" << lastValue << " sub(" << 0 << ")" << cellMax(to) << "]";
		printf("\n@@@@X@@@@\n");
	}
}
