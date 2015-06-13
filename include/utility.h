
#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <string>
#include <vector>
#include <set>

#include <sstream>
#include <iomanip> // setprecision

#include <boost/algorithm/string.hpp>

using namespace std;

string combine_spaces(string str,char ch=' ');

// generate a random string of letters and numbers of certain length
string random_string( int length );

string to_upper(string str);

string to_string(vector<string> str, string del="\t");
string to_string(vector<int> str, string del="\t");
string to_string(set<string> str, string del="\t");

template <typename T>
string to_string_with_precision(const T a_value, const int n = 3)
{
    ostringstream out;
    out << setprecision(n) << a_value;
    return out.str();
}


// split string
vector<string>  string_split(string str, string separator="\t,| ");

vector<string> set_overlap(set<string> first, set<string> second);

set<string> set_subtract(set<string> s1, set<string> s2);

// get current time, in the format: Thu Mar 15 21:06:57 2012
string current_time();

// write message to standard error starting with time
void message(string text, bool stdout=false);

// run system command
void system_run(string cmd);

// run R script, can be multiple lines, such as
/*

    string script = 
    "pdf('"+outputfile+"',width=10,height=5) \n"
    "con  <- file('"+inputfile+"', open = 'r') \n"
    "while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) { \n"
    "   myVector <- (strsplit(oneLine, '\\t')) \n"
    "   x=myVector[[1]] \n"
    "   name=x[1] \n"
    "   pos=x[2] \n"
    "   counts=as.numeric(x[3:length(x)])/"+to_string(nSeq)+"\n"
    "   plot(counts,type='h',xlab='position',ylab='frequency',main=paste(name,pos,sep=',')) \n"
    "} \n"
    "close(con) \n"
    "dev.off() \n";

*/
void R_run(string script,bool clean=true, string Rcmd="R CMD BATCH");

void sort_file_and_add_header(string filename, string header, string sort_options);


#endif
