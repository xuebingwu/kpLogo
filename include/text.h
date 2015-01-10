
#ifndef __TEXT_H__
#define __TEXT_H__

#include <string>
#include <vector>

using namespace std;

int count_lines(string filename);

int split_file_for_cross_validation(string input, string output, int nfold);

void select_multi_lines(string inputfile, string idfile, string outputfile, int nline=2, int id_col=1, string id_prefix=">");

void intersectTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool subtract=false);

void mergeTab(string file1, string file2, string outputfile, unsigned col1=0, unsigned col2=0, bool header=false, string fill="none");

void remove_duplicates(string input, string output, int col, int max, string sort_opts);

int find_unique_lines(string input, string output, int col);

#endif

