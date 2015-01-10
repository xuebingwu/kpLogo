#ifndef __STRUCTURE_H__
#define __STRUCTURE_H__

#include <fstream>

#include "sequence.h"

using namespace std;


void RNALfold_to_RNAfold(string infile, string outfile, int min_hairpin_length=50, int ext=20);
void RNALfold_to_RNAfold_filter_overlap(string infile, string outfile, int min_hairpin_length=50, int ext=20,int max_overlap_allowed=20);

void convert_RNALfold_to_RNAfold_output(string filename);

void hairpin_scoring(string infile, string outfile, vector<positional_kmer> model_str,vector<positional_kmer> model_top,vector<positional_kmer> model_bot,vector<positional_kmer> model_loop, int model_length);

int remove_short_stem(string &structure, int max_paired_length=4);

int remove_short_stem_from_file(string inputfile, string outputfile, int max_paired_length=4,int min_paired_bases_left=50);

int keep_longest_stem(string struc);

#endif