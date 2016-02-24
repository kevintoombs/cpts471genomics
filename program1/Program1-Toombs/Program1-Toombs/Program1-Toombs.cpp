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
	int dir = 0;
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
void buildTable(DP_table &t);
void initTable(DP_table &t);
void calcTable(DP_table &t);
int maximum(int S, int D, int I, int alignmentType);
int subFunction(char a, char b, config c);
void printTable(DP_table &t);
void retrace(DP_table &t);
int direction(DP_cell &c);
int maximum2(int S, int D, int I, int alignmentType, DP_cell &c);
int cellMax(DP_cell c, int alignmentType);

int main(int argc, char *argv[])
{
	DP_table t;
	t.alightmentType = getAlignmentType(argc, argv);
	t.c = getConfig(argc, argv); 
	if (!parseFasta(argv, t) )
	{
		return 1;
	} 

	printf("Scores: Match = %i, Mismatch = %i, h = %i, g = %i\n", t.c.matchScore, t.c.matchScore, t.c.startGapScore, t.c.continueGapScore);
	cout << endl;

	printf("Sequence 1 = %s, length = %i characters\n", "s1", t.sequence1.length());
	printf("Sequence 2 = %s, length = %i characters\n", "s2", t.sequence2.length());
	cout << endl;
	
	buildTable(t);
	calcTable(t);
	retrace(t);

    return 0;
}

void buildTable(DP_table &t)
{
	printf("Building Table.");
	DP_cell c;
	t.t.resize(t.sequence1.length()+1, vector<DP_cell>(t.sequence2.length()+1, c) );
	cout << endl;

	return;
}

void calcTable(DP_table &t)
{
	if (t.alightmentType == 0) initTable(t);
	int i, j, maxValue = 0;
	tuple<int, int> maxPair;

	printf("Calculating Table[characters remaining]:");
	for (i = 1; i <= t.sequence1.length(); i++)
	{
		if (i % 15 == 0) printf("[%i]", t.sequence1.length() - i);
		for (j = 1; j <=  t.sequence2.length(); j++)
		{
			t.t[i][j].S = maximum2
				(t.t[i - 1][j - 1].S + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.t[i - 1][j - 1].D + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.t[i - 1][j - 1].I + subFunction(t.sequence1[i - 1], t.sequence2[j - 1], t.c),
					t.alightmentType, t.t[i - 1][j - 1]);
			t.t[i][j].D = maximum2
				(t.t[i - 1][j].S + t.c.startGapScore + t.c.continueGapScore,
					t.t[i - 1][j].D + t.c.continueGapScore,
					t.t[i - 1][j].I + t.c.startGapScore + t.c.continueGapScore,
					t.alightmentType, t.t[i - 1][j]);
			t.t[i][j].I = maximum2
				(t.t[i][j - 1].S + t.c.startGapScore + t.c.continueGapScore,
					t.t[i][j - 1].D + t.c.continueGapScore + t.c.continueGapScore,
					t.t[i][j - 1].I + t.c.startGapScore,
					t.alightmentType, t.t[i][j - 1]);
			if (t.alightmentType == 1)
			{
				if (cellMax(t.t[i][j], t.alightmentType) > maxValue)
				{
					maxValue = cellMax(t.t[i][j], t.alightmentType);
					maxPair = make_tuple(i, j);
					t.maxPair = maxPair;
				}
			}
		}
	} 
	cout << endl << endl;
	t.t[i-1][j-1].dir = direction(t.t[i-1][j-1]);

	
}

void initTable(DP_table &t)
{
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

		while (fasta.get(ch)) 
		{

			if (ch == '>' && state == 0 && getId == 0 )
			{
				getId = 1;
			}
			if (ch == '>' && getId == 1)
			{
				if (ch == ' ') 
				{ 
					getId = 2;
					break;
				}
				t.id1 += ch;
			}
			if (ch == '\n' && state == 0)
			{
				fasta.get(ch);
				state = 1;
			}
			if (state == 1) 
			{
				if (ch == 'A' || ch == 'C' || ch == 'G' || ch == 'T' || ch == 'a' || ch == 'c' || ch == 'g' || ch == 't')
				{
					t.sequence1 += toupper(ch);
				}
				if (ch == '>') 
				{
					state = 2;
				}
			}
			if (getId == 2 && state == 2)
			{
				if (ch == ' ')
				{
					getId = 3;
					break;
				}
				t.id2 += ch;
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
				}
				//cout << ch;
				if (ch == '>') {
					state = 4;
				}
			}
		}
		if (state == 3) 
		{
			
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
	int matches = 0, mismatches = 0, gaps = 0, openingGaps = 0;
	int lastDir = -1;
	int dir = -1;
	stack<char> s1, s2, r;
	int lastValue = 1;
	int i = t.sequence1.length();
	int j = t.sequence2.length();

	if (t.alightmentType == 1)
	{
		int i = get<0>(t.maxPair);
		int j = get<1>(t.maxPair);
	}


	while (i >= 0 && j >= 0)
	{
		if (t.alightmentType == 1 && cellMax(t.t[i][j], t.alightmentType) == 0)
		{
			break;
		}
		lastDir = dir;
		dir = t.t[i][j].dir;
		if (dir == 1)
		{
			j--;
			s1.push('-');
			s2.push(t.sequence2[j]);
			r.push(' ');
			if (lastDir == dir)
			{
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
				//cout << "gap" << endl;
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
	
	while (!s1.empty())
	{
		int to60 = 0;

		// printf("s1 %5d ", );
		while (to60 != 60)
		{
			if (!s1.empty()) cout << ' ' << s1.top();
			if (!s1.empty()) s1.pop();
			to60++;
		} 
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!r.empty()) cout << ' ' << r.top();
			if (!r.empty()) r.pop();
			to60++;
		} 
		to60 = 0;
		cout << endl;
		while (to60 != 60)
		{
			if (!s2.empty()) cout << ' ' << s2.top();
			if (!s2.empty()) s2.pop();
			to60++;
		}
		to60 = 0;
		cout << endl;
	}
	cout << endl;


	cout << "Report:\n\n";
	if (t.alightmentType == 0) printf("Global optimal score = %i\n\n", t.t[t.sequence1.size()][t.sequence2.size()]);
	if (t.alightmentType == 1) printf("Local optimal score = %i\n\n", t.t[get<0>(t.maxPair)][get<1>(t.maxPair)]);
	printf("Number of:  matches = %i, mismatches = %i, gaps = %i, opening gaps = %i\n", matches, mismatches, gaps, openingGaps);
	float total = gaps + matches + mismatches;
	float identities = matches + mismatches;
	float p1 = 100 * identities / total;
	float p2 = 100 * gaps / total;
	printf("Identities = %.0f/%.0f (%2.1f percent), Gaps = %i/%.0f (%2.1f percent)\n", identities, total, p1, gaps, total, p2);
}

int direction(DP_cell &c) 
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