#ifndef __MARKOV_H__
#define __MARKOV_H__

#include <string>
#include <vector>
#include <map>

using namespace std;

// a class for markov model
class markov_model {
    public:
        int order; // 0,1,2
        string alphabet; // ACGT or others
        map<char,double> f1; // mononucleotide frequency
        map<string,double> f2; // dinucleotide frequency
        map<string,double> f3; // trinucleotide frequency
        map<char,map<char,double> > first_order; // f2/f1
        map<string, map<char,double> > second_order; // f3/f2
        markov_model();
        markov_model(int,string,vector<string>);
		markov_model(string,string);
        map<string,double> probs(vector<string>);
        void print(string);
};

#endif
