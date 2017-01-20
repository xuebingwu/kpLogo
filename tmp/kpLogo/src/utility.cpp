#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

#include "utility.h"
#include <boost/regex.hpp>


using namespace std;


//

// combine multiple spaces into one
bool BothAreSpaces(char lhs, char rhs) { return (lhs == rhs) && (lhs == ' '); }
bool BothAreTabs(char lhs, char rhs) { return (lhs == rhs) && (lhs == '\t'); }

string combine_spaces(string str, char ch)
{
	std::string::iterator new_end;
	if(ch == ' ') new_end = std::unique(str.begin(), str.end(), BothAreSpaces);
	else new_end = std::unique(str.begin(), str.end(), BothAreTabs);
	str.erase(new_end, str.end());   
	return str;
}

// generate a random string of letters and numbers of certain length
/*
string random_string( size_t length )
{
    srand(time(NULL));
    auto randchar = []() -> char
    {
        const char charset[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = (sizeof(charset) - 1);
        return charset[ rand() % max_index ];
    };
    string str(length,0);
    generate_n( str.begin(), length, randchar );
    return str;
}
*/

/**/
string random_string( int len) {
	string s;
    static const char alphanum[] =
        "0123456789"
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz";

    for (int i = 0; i < len; ++i) {
        s += alphanum[rand() % (sizeof(alphanum) - 1)];
    }
	return s;
}
/**/

string to_string(vector<int> str, string del/*="\t"*/)
{
	string res = to_string(str[0]);
	for(unsigned i=1;i<str.size();i++)
		res += del + to_string(str[i]);
	return res;
}


string to_string(vector<string> str, string del/*="\t"*/)
{
	string res = str[0];
	for(unsigned i=1;i<str.size();i++)
		res += del + str[i];
	return res;
}

string to_string(set<string> str, string del/*="\t"*/)
{
	string res;
	set<string>::iterator it;
	for(it=str.begin();it!=str.end();it++)
		res += *it + del;
	return res;
}

vector<string> set_overlap(set<string> first, set<string> second)
{
	vector<string> v(first.size());
	vector<string>::iterator it = set_intersection (first.begin(),first.end(), second.begin(), second.end(), v.begin());
	v.resize(it-v.begin());
	return v;      
}

set<string> set_subtract(set<string> s1, set<string> s2)
{
	set<string> res;
	set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),std::inserter(res, res.end()));
	return res;
}

/*
// split string
vector<string>  string_split(string str, string separator){
    size_t found;
    vector<string> results;
    found = str.find_first_of(separator);
    while(found != string::npos){
        if(found > 0){
            results.push_back(str.substr(0,found));
        }
        str = str.substr(found+1);
        found = str.find_first_of(separator);
    }
    if(str.length() > 0){
        results.push_back(str);
    }
    return results;
}
*/

bool ends_with(string const & value, string const & ending)
{
    if (ending.size() > value.size()) return false;
    return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

vector<string>  string_split(string str, string separator/*="\t,| "*/){
	vector<string> parts;
	boost::split(parts, str, boost::is_any_of(separator));
	return parts;
}

string to_upper(string strToConvert)
{
    std::transform(strToConvert.begin(), strToConvert.end(), strToConvert.begin(), ::toupper);
    return strToConvert;
}

string current_time()
{// get current time, in the format:
 // Thu Mar 15 21:06:57 2012
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  string str = asctime (timeinfo);
  return str.substr(0,str.size()-1);
}

// write message to standard error
void message(string text, bool stdout/*=false*/, string pre/*=""*/)
{
	if (stdout) cout << pre + "["<<current_time()<<"] " + text << endl;
	else cerr << pre + "["<<current_time()<<"] " + text << endl;
}

// run system command
void system_run(string cmd)
{
    system(cmd.c_str());
}

// sort a file with no header, add header
void sort_file_and_add_header(string filename, string header, string sort_options)
{
	// tmp file name, a random string
	string tmp = random_string(10);
	// mv data to tmp
	system_run("mv "+filename+" "+tmp);
	// write header to output file
	ofstream fout;
	fout.open(filename.c_str());
	fout << header << endl;
	fout.close();
	// sort and append data to header
	system_run("sort "+tmp+" "+sort_options+" >> "+filename);
	// remove tmp file
	system_run("rm "+tmp);
}

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
void R_run(string script, bool clean/*=true*/, string Rcmd/*="R CMD BATCH"*/)
{
    // create tmp R script file
	srand(time(NULL));
    string tmp = random_string(11)+".r";
    ofstream out;
    out.open(tmp.c_str());
    out << script;
    out.close();

    // run the script
    string cmd = Rcmd + " "+tmp;
    system(cmd.c_str());

    // remove the script and .Rout file
	if(clean)
	{
    	cmd = "rm "+tmp+" "+tmp+".Rout";
    	system(cmd.c_str());    
	}	
}

