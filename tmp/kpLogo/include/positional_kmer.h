#ifndef __POSITIONALKMER_H__
#define __POSITIONALKMER_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array

#include <stdio.h>
#include <stdlib.h>

using namespace std;

class positional_kmer
{
public:
    string seq; // sequence, IUPAC code
    int pos;
    int size;
    double weight;
    int group; // if part of another, will have the same group number
    positional_kmer();
    positional_kmer(const positional_kmer &a);
    positional_kmer(string seq, int pos, int size, double weight=0, int group=0);
    bool equals(positional_kmer a);
    bool is_part_of(positional_kmer a);
    bool is_present_in(string seq);
    const positional_kmer &operator=(const positional_kmer &a);
    string as_string(string del="_");
    vector<positional_kmer> nondegenerate(map<char,string> degenerate_map);
};

vector<positional_kmer> positional_kmer_vector_from_string(string str, map<char,string> degenerate_map);

#endif
