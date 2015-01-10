#include <fstream>

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <algorithm>

#include "utility.h"
#include "text.h"
#include "time.h"
#include <map>
#include <iostream>



// split files for crossvalidation
// if each fold has less than 2 samples,do leave-one-out
int split_file_for_cross_validation(string input, string output, int nfold)
{
	// determine the number of lines
	ifstream fin(input.c_str());
	string line;
	int total = 0;
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		total ++;
	}
	fin.close();
	
	int step = floor(0.5 + total / nfold);
	
	// leave one out cv:LOOCV
	if (step < 2) 
	{
		nfold = total;
		step = 1;
	}
	
	vector<int> bins;
	for(int i=0;i < nfold-1;i++)
		for(int j=0;j<step;j++)
			bins.push_back(i);
	for(int i= bins.size(); i < total; i++) bins.push_back(nfold-1);
	
	//cout << "total=" << total << ",step=" << step << ",bins=" << bins.size() << endl; 
	srand(time(NULL));
	random_shuffle(bins.begin(),bins.end());
	
    // split files
	for(int i=0;i<nfold;i++)
	{
		string trainfile = output+".train."+to_string(i);
		string testfile = output+".test."+to_string(i);
		
		ofstream train(trainfile.c_str());
		ofstream test(testfile.c_str());
		
		int n = 0;
		fin.open(input.c_str());
		while(fin)	
		{
			getline(fin,line);
			if(line.length() == 0) continue;
			if (i == bins[n]) test << line << endl;
			else train << line << endl;
			n++;
		}
		fin.close();
		train.close();
		test.close();
	}
	return nfold;
}

// input file contine multi-line records such as fasta/fastq or others
// col is 1 based
void select_multi_lines(string inputfile, string idfile, string outputfile, int nline, int id_col, string id_prefix)
{
	id_col = id_col - 1;
	
	set<string> ids;
	
	ifstream fin;
	fin.open(idfile.c_str());

	string line;
	vector<string> flds;
	
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		flds = string_split(line,"\t");
		ids.insert(id_prefix+flds[id_col]);
	}
	fin.close();
	
	ofstream out;
	out.open(outputfile.c_str());
	
	fin.open(inputfile.c_str());
	while(fin)
	{
		getline(fin,line);
		if(line.length()==0) continue;
		if(line[0] == '>')
		{
			if (ids.find(line) != ids.end())
			{
				out << line << endl; // id
				for (int i =0;i< nline-1;i++)
				{
					getline(fin,line);
					out << line << endl; //
				}
			} 
			else
			{
				for (int i =0;i< nline-1;i++)
				{
					getline(fin,line);
				}
			}
		}
	}
	fin.close();
	out.close();
}


// subtract file 2 from file 1, based on shared key column
void intersectTab(string file1, string file2, string outputfile, unsigned col1/*=0*/, unsigned col2/*=0*/, bool subtract/*=false*/)
{
    // read file2 first 
    set<string> data2;
    vector<string> flds;
	string line;
	
	ifstream f2;
	f2.open(file2.c_str());
    
	while(f2)
	{
		getline(f2,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		data2.insert(flds[col2]);
	}
	f2.close();
	
	ifstream f1;
	f1.open(file1.c_str());
	
	ofstream out;
	out.open(outputfile.c_str());
    
	unsigned total = 0;
	unsigned comm = 0;
	unsigned diff = 0;
	
	while(f1)
	{
		getline(f1,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		total ++;
		//cout << flds[col1] << endl;
		if (data2.find(flds[col1]) != data2.end()) // shared
		{
		    if(subtract == false) 
			{
				out << line << endl; // output common lines
				comm ++;
			}
			
		}
		else if (subtract) // not shared and output difference
		{
			out << line << endl;
			diff ++;
		}
	}          
	
	if (subtract) cout << diff << " of " << total << " lines unique to " << file1 << endl;
	else cout << comm << " of " << total << " lines shared" << endl;
	
    f1.close();  
	out.close();         
}

void mergeTab(string file1, string file2, string outputfile, unsigned col1/*=0*/, unsigned col2/*=0*/, bool header/*=false*/, string fill/*="None"*/)
{
	/*
    // add file2 to file1 side by side, using column col1 in file1 and column col2 in file2 to match rows
    // col1 and col2 are 0-based, but the program may be 1-based
    // combine header directly if header = true
    // for lines in file1 but not file2, fill with 'fill' if fill != "none"
	*/
	
    // read file2 first 
    map<string,string> data2;
    vector<string> flds;
	string line;
	
	ifstream f2;
	f2.open(file2.c_str());

	string header2;
	if(header) getline(f2,header2);
    
	while(f2)
	{
		getline(f2,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		data2[flds[col2]] = line;
	}
	f2.close();
	
    unsigned n2 = data2.size();
    unsigned nc = flds.size();
	
	//cout << n2 << "," << nc << endl;
	
    // match and add to file1
	ifstream f1;
	f1.open(file1.c_str());
	
	ofstream out;
	out.open(outputfile.c_str());
	
	string header1;
	if(header) 
	{
		getline(f1,header1);
		header1 = header1 + "\t" +header2;
	}
    //cout << header1 << endl;
	
    unsigned n1 = 0; // lines in file 1
    unsigned n3 = 0; // common lines
    
    string fillline = fill;
	if (fill != "none") 
	{
		for (int i=1;i<nc;i++) fillline += "\t" + fill;
	}

	while(f1)
	{
		getline(f1,line);
		if(line.length() == 0) continue;
		flds = string_split(line,"\t");
		n1 ++;
		//cout << flds[col1] << endl;
		if (data2.find(flds[col1]) != data2.end())
		{
			n3 ++;
			out << line << "\t" << data2[flds[col1]] << endl;
		}
		else if (fill != "none")
		{
			out << line+ "\t" + fillline << endl;
		}
	}          
    f1.close();      
	out.close();     
}
		
void remove_duplicates(string input, string output, int col, int max, string sort_opts="")
{
    if (sort_opts.size()>0)
    {
        string cmd = "sort "+sort_opts+" "+input+" > xxx.tmp";
        system(cmd.c_str());
        input="xxx.tmp";
    }


    col = col - 1;
    ifstream fin;
    ofstream fout;
    fin.open(input.c_str());
    fout.open(output.c_str());
    string line;
    vector<string> flds;

    string curr = "";
    int k = max;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
        flds = string_split(line,"\t");
        if (curr != flds[col])
        {
            fout << line << endl;;
            k = max - 1;
            curr = flds[col];
        }
        else if (k>0)
        {
            fout << line << endl;
            k--;
        }
    }

    if(sort_opts.size()>0) system("rm xxx.tmp");

}

int count_lines(string filename)
{
	ifstream fin(filename.c_str());
	string line;
	int n = -1;
	while(fin) 
	{
        getline(fin,line);
		n++;
	}
	fin.close();
	return n;
}
// find lines the key column is unique
int find_unique_lines(string input, string output, int col)
{
    string cmd = "sort -k"+to_string(col)+","+to_string(col)+ " "+input+" > xxx.tmp";
	//message("sort input file: "+ cmd);
    system(cmd.c_str());
    input="xxx.tmp";

    col = col - 1;
    ifstream fin;
    ofstream fout;
    fin.open(input.c_str());
    fout.open(output.c_str());
    string line;
    vector<string> flds;

    getline(fin,line);
    flds = string_split(line,"\t");
    string pre = flds[col];
	string pre_line = line;
	
	bool duplicated = false;

	int n = 0;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
        flds = string_split(line,"\t");
        if (pre != flds[col]) // a new line
        {
			if(duplicated == false)  
			{
				fout << pre_line << endl;
				n ++; 
			}
			pre = flds[col];
			pre_line = line;
			duplicated = false;
		}
		else duplicated = true;
    }

	
   	system("rm xxx.tmp");

	return n;
}


