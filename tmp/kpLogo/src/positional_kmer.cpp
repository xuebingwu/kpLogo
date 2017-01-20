#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array
#include <algorithm>

#include <stdio.h>
#include <stdlib.h>

#include "utility.h"
#include "positional_kmer.h"
#include "sequence.h"

#include <boost/algorithm/string.hpp> // erase_all


// build a vector from seq:pos:shift
// position: 1 based
// allow degenerate
vector<positional_kmer> positional_kmer_vector_from_string(string str, map<char,string> degenerate_map)
{
	vector<positional_kmer> pkmers;
	// remove space
	boost::erase_all(str," ");	
	vector<string> all_pkmers = string_split(str,",");
	for(int i=0;i<all_pkmers.size();i++)
	{
		vector<string> flds = string_split(all_pkmers[i],":");// seq, pos, shift
		if (flds.size() == 1){
			cerr << "ERROR: incorrect positional kmer string: "+all_pkmers[i] << endl;
			exit(1);
		} else if (flds.size() == 2){
			flds.push_back("0");
		}
		vector<string> seqs = expand_degenerate_kmer(flds[0], degenerate_map);
		int pos = stoi(flds[1])-1;
		int size = stoi(flds[2])+flds[0].size();
		for(int j=0;j<seqs.size();j++) pkmers.push_back(positional_kmer(seqs[j],pos,size));
	}
	return pkmers;
}

positional_kmer::positional_kmer(string seq1, int pos1, int size1, double weight1, int group1){
    seq = seq1;
    pos = pos1;
    size = size1;
    weight = weight1;
    group = group1;
}

positional_kmer::positional_kmer()
{
  // allocate variables
    seq = "";
    pos = -1;
    size = -1;
    weight = -1.0;
    group = -1;
}

positional_kmer::positional_kmer(const positional_kmer &a)
{
  // allocate variables
  positional_kmer();
  // copy values
  operator = (a);
}

bool positional_kmer::equals(positional_kmer a)
{
    return seq == a.seq && pos == a.pos && size == a.size;
}

bool positional_kmer::is_part_of(positional_kmer a)
{
    // only if a is part of this at the same position
    if (pos >= a.pos && ( (pos + seq.size()) <= (a.pos + a.seq.size()) ) && seq.size() != a.seq.size())
        if (a.seq.substr(pos - a.pos, seq.size()) == seq)
            return true;
    return false;
}

// is found in a sequence
bool positional_kmer::is_present_in(string sequence)
{
    return sequence.substr(pos,size).find(seq) != std::string::npos;
}

const positional_kmer &positional_kmer::operator = (const positional_kmer &a)
{
    seq = a.seq;
    pos = a.pos;
    size = a.size;
    weight = a.weight;
    group = a.group;
    return *this;
}

string positional_kmer::as_string(string del/*="_"*/)
{
    return seq+del+to_string(pos)+del+to_string(size)+del+to_string(weight)+del+to_string(group);
}

vector<positional_kmer> positional_kmer::nondegenerate(map<char,string> degenerate_map)
{
    vector<positional_kmer> res;
    vector<string> seqs = expand_degenerate_kmer(seq, degenerate_map);
    for(int i=0;i<seqs.size();i++) res.push_back(positional_kmer(seqs[i], pos,size, weight, group));
    return res;
}


