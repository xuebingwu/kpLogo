
#ifndef __UTILITY_H__
#define __UTILITY_H__

#include <string>
#include <vector>
#include <set>


using namespace std;



// generate a random string of letters and numbers of certain length
string random_string( size_t length );

string to_upper(string str);

string to_string(vector<string> str, string del="\t");
string to_string(set<string> str, string del="\t");


// split string
vector<string>  string_split(string str, string separator);

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
