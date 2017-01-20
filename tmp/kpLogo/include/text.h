
#ifndef __TEXT_H__
#define __TEXT_H__

#include <string>
#include <vector>
#include "utility.h"
#include <map>


using namespace std;

int load_scores(string filename, vector<string> &names, vector<double> &scores, int c1=1, int c2=2);

void insert_header(string filename, string header);

int count_lines(string filename);

int split_file_for_cross_validation(string input, string output, int nfold);

int select_multi_lines(string inputfile, string idfile, string outputfile, int nline=2, int id_col=1, string id_prefix=">");

void intersectTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool subtract=false);

void mergeTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool header=false, string fill="none");

void remove_duplicates(string input, string output, int col, int max, string sort_opts="", bool print_rank=false);

int find_unique_lines(string input, string output, int col);

// load weight data, column c1 is id, c2 is the weight
map<string,double> load_weight(string filename, int c1, int c2);

#endif

