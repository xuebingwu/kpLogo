#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array

//#include <boost/regex.hpp>

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/lexical_cast.hpp>

	
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>       /* log10 */
#include <algorithm>
#include "stat.h"
#include "utility.h"
#include "sequence.h"
#include "container.h"
#include "markov.h"
	
extern "C"{
#include "ushuffle.h"
}

using namespace std;

// number of identical bases
int sequence_similarity(string a, string b)
{
	int s = 0;
	if(a.size() != b.size()) return -1;
	for(int i=0;i<a.size();i++)
	{
		if(a[i] == b[i]) s++;
	}
	return s;
}

// calculate and write pairwise similarity matrix of
void pairwise_sequence_similarity_matrix(vector<string> seqs, string filename)
{
	ofstream out(filename.c_str());
	map<string,int> data;
	for(int i =0;i<seqs.size();i++)
	{
		for(int j=0;j<seqs.size();j++)
		{
			if (j<i) out << data[to_string(j)+","+to_string(i)] << "\t";
			else if (i == j) out << seqs[i].size() << "\t";
			else
			{
				data[to_string(i)+","+to_string(j)] = sequence_similarity(seqs[i],seqs[j]);
				out << data[to_string(i)+","+to_string(j)] << "\t";
			}
		}
		out << endl;
	}
	out.close();
}

// load weighted sequence file into two vectors: seqs and weights
// the first three columns should be: id, seq, weight (tab-delimited)
// if cWeight < 0, load unweighted sequences
void load_sequences_from_tabular(string filename, vector<string> &seqs, vector<double> &weights, int cSeq/*=1*/, int cWeight/*=2*/) {
	// minimum number of columns in input
	int nCol = max(cSeq,cWeight)+1;
	
	ifstream fin;
	fin.open(filename.c_str());

    string line;
    vector<string> flds;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		line.erase(line.find_last_not_of(" \n\r\t")+1);
		
		if(line[0] == '#') continue;

        flds = string_split(to_upper(line));
		if (flds.size() < nCol) message("too few columns in line: "+line);
		else
		{
			seqs.push_back(flds[cSeq]);
			if (cWeight>=0) weights.push_back(stof(flds[cWeight]));
		}
	}
	fin.close();
}

map<string,string> vector2map(vector<string> seqs)
{
	map<string,string> res;
	for (int i=0;i<seqs.size();i++) res[to_string(i)] = seqs[i];
	return res;
}

// load weighted sequence file into two vectors: seqs and weights
// the first three columns should be: id, seq, weight (tab-delimited)
void load_weighted_sequences_to_vectors(string filename, vector<string> &seqs, vector<double> &weights, int cSeq/*=1*/, int cWeight/*=2*/) {
	// minimum number of columns in input
	int nCol = max(cSeq,cWeight)+1;
	
	ifstream fin;
	fin.open(filename.c_str());

    string line;
    vector<string> flds;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

		if(line[0] == '#') continue;

        flds = string_split(line);
		if (flds.size() < nCol) message("too few columns in line: "+line);
		else
		{
			seqs.push_back(flds[cSeq]);
			weights.push_back(stof(flds[cWeight]));
		}
	}
	fin.close();
}

// load ranked sequences to vectors
// default: no header, first column is id, second column is sequence
// c: column for sequence, 0-based, i.e. first column is 0
vector<string> load_ranked_sequences_to_vectors(string filename, int cSeq/*=1*/){
	vector<string> seqs;
	
	ifstream fin;
	fin.open(filename.c_str());

    string line;
    vector<string> flds;
	
	if (cSeq<0) // input is fasta
	{
	    while(fin)
	    {
	        getline(fin,line);
	        if (line.length() == 0)
	            continue;
			if(line[0] == '#') continue;
			
			if(line[0] != '>') seqs.push_back(line);
		}
		fin.close();
		return seqs;	
	}
	

	// read sequences into a vector
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		if(line[0] == '#') continue;

        flds = string_split(line);
		if (flds.size() < cSeq+1) message("skip lines with not enough columns: "+line);
		else seqs.push_back(flds[cSeq]);
	}
	fin.close();
	
	return seqs;
}

// remove sequences with size differ from lSeq
// if lSeq = 0 (default), use the length of the first 
vector<int> filter_sequences_by_size(vector<string> &seqs, int lSeq/*=0*/)
{
	vector<int> removed;
	
	if (lSeq == 0) lSeq = seqs[0].size();
	
	for(int i=seqs.size()-1;i>=0;i--)
	{
		if(seqs[i].size() != lSeq)
		{
			seqs.erase(seqs.begin()+i);
			removed.push_back(i);
		}
	}
	return removed;
}


// find significant positional kmer from ranked sequences
// use non-parametric U test
// input: a vector of ranked sequence with the same length, only ACGT 
int find_significant_kmer_from_ranked_sequences(
	vector<string> seqs, 
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	bool Bonferroni/*=true*/,
	int min_shift/*=0*/, 
	int max_shift/*=0*/,
	int startPos/*=0*/,
	int minCount/*=5*/ )
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	int nSig = 0;
		
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<kmers.size();i++) // for each kmer
	{
		// expand a degenerate kmer to all possible element/exact kmers
		//cout << kmers[i] << endl;
        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		int k = kmers[i].size();
		for (int shift = min_shift; shift <= max_shift; shift ++)
		{
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				// ranks of sequences with this kmer at this position
				vector<int> ranks;
				// for each sequence, 
				for( int j=0;j<nSeq;j++) 
				{
					//find if any of the expanded kmer is present at position pos
					for( int n=0;n<exp_kmers.size();n++)
					{
						/* speed not affected by shift */
						size_t found = seqs[j].substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
						if (found!=std::string::npos)
						{
							// add this sequence's rank to sample 1
							ranks.push_back(j);
							// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
							break;
						} /**/
						// another implementation, slower, linear with shift
						/*
						bool found = false;
						for(int s=0;s<=shift;s++)
						{
							if(seqs[j].substr(pos+s,k) == exp_kmers[n]) 
							{
								ranks.push_back(j);
								//cout << exp_kmers[n] <<","<< seqs[j].substr(pos+s,k) << endl;
								found = true;
								break;
							}
						}
						if (found) break;
						*/
					}
				}
				// now ranks includes ranks of all sequences containing kmer i at position pos
				// if less than 3 sequence contain this kmer, or less than 3 sequence don't have this kmer,
				// just go on to the next position
				if (ranks.size() < minCount || ranks.size() > nSeq-minCount ) 
				{
					continue;
				}
				//cout << kmers[i] << "\t" << pos << "\t" << ranks.size() << endl;
				array<double,2> utest = Mann_Whitney_U_test(ranks, nSeq);
				if (utest[1] < pCutoff)
				{
					double pB = min(1.0,utest[1]*nTest);
					if (Bonferroni && pB > pCutoff) continue;
					nSig ++;
					int position = pos-startPos+2;
					if (position < 1) position -= 1;
		            outstream << kmers[i] << "\t" << position  << "\t" << shift << "\t" << utest[0] << "\t" << -log10(utest[1]) << "\t" << -log10(pB) << endl;
				}
			}
		}
	}

    return nSig;
} // end of function

// find significant positional kmer from weighted sequences
// use t-test
// input: a vector of weighted sequence with the same length, only ACGT 
// 
int find_significant_kmer_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	bool Bonferroni/*=true*/,
	int min_shift/*=0*/, 
	int max_shift/*=2*/,
	int startPos/*=0*/,
	int minCount/*=5*/) 
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	int nSig = 0;
		
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<kmers.size();i++) // for each kmer
	{
		// expand a degenerate kmer to all possible element/exact kmers		
        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		int k = kmers[i].size();
		for (int shift = min_shift; shift <= max_shift; shift ++)
		{
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				//debug cout << kmers[i] << "@" << pos << endl;
			
				// whether each sequence is positive
				vector<bool> positive;
					// for each sequence, 
				for( int j=0;j<nSeq;j++) 
				{
					//debug cout << "seq " << j << endl;
				
					// initialize
					bool present = false;
					//find if any of the expanded kmer is present at position pos
					for( int n=0;n<exp_kmers.size();n++)
					{
						//debug cout << "exact kmer: " << exp_kmers[n] << endl;
						/* speed not affected by shift */
						size_t found = seqs[j].substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
						if (found!=std::string::npos)
						{
							// add this sequence's rank to sample 1
							//weights1.push_back(weights[j]);
							// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
							present = true; 
							break;
						} /**/
					}
				
					//debug cout << "present = " << present << endl;
				
					positive.push_back(present);
				}
				vector<double> weights1,weights2;
				for( int j=0;j<nSeq;j++)
				{
					if(positive[j]) weights1.push_back(weights[j]);
					else weights2.push_back(weights[j]);
				}
				// now weights includes weights of all sequences containing kmer i at position pos
				// if less than 3 sequences with/without this kmer, just go on to the next position
				// in such cases the kmer is unlikely to be significant
				// also the t.test will not work well. it requires at least 2 data points in each sample
				if (weights1.size() < minCount || weights2.size() < minCount ) 
				{
					continue;
				}
				array<double,6> ttest = t_test(weights1,weights2,false);
				if (ttest[1] < pCutoff)
				{
					// Bonferoni correction
					double pB = min(1.0,ttest[1]*nTest);
					if (Bonferroni && pB > pCutoff) continue;
					nSig ++;
					int position = pos-startPos+2;
					if (position < 1) position -= 1;
		            outstream << kmers[i] << "\t" << position << "\t" << shift << "\t" << ttest[0] << "\t" << -log10(ttest[1]) << "\t" << -log10(pB) << "\t" << weights1.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << weights2.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
				}
			}
		}
	}

    return nSig;
} // end of function


// only calculate local enrichment, no need for backgrounintd
int find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
	vector<string> seqs,
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	bool Bonferroni/*=true*/,
	int min_shift/*=0*/, 
	int max_shift/*=2*/,
	int startPos/*=0*/,
	int minCount/*=5*/) 
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	int nSig = 0;
		
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<kmers.size();i++) // for each kmer
	{
		// expand a degenerate kmer to all possible element/exact kmers		
        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		int k = kmers[i].size();
		for (int shift = min_shift; shift <= max_shift; shift ++)
		{
			vector<int> counts = {};
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				int n1 = 0;
					// for each sequence, 
				for( int j=0;j<nSeq;j++) 
				{
					//debug cout << "seq " << j << endl;
					//find if any of the expanded kmer is present at position pos
					for( int n=0;n<exp_kmers.size();n++)
					{
						//debug cout << "exact kmer: " << exp_kmers[n] << endl;
						/* speed not affected by shift */
						size_t found = seqs[j].substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
						if (found!=std::string::npos)
						{
							// add this sequence's rank to sample 1
							//weights1.push_back(weights[j]);
							// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
							n1 ++ ; 
							break;
						} /**/
					}
				}
				counts.push_back(n1);
			}
			
	        int total_counts = sum(counts);
			
			double expected = double(total_counts)/counts.size();
			
			// calculate background f2
			double f2 = expected / nSeq;
			
			//double sigma = sd(counts) / nSeq;
			for(int x=0;x<counts.size();x++)
			{
				if (counts[x] < minCount || counts[x] > nSeq - minCount) continue;
				//double f1 = float(counts[x])/nSeq;
				//double z = (f1 - f2)/sigma;
				//double p = p_value_for_z_score(z);
				double p = binom_test(nSeq,counts[x],f2);

	            if (p < pCutoff)
	            {					
					if(p == 0) p = 1e-16;

	                double corrected_p = min(1.0,p*nTest);	
					if (Bonferroni && corrected_p > pCutoff) continue;
					
	                nSig++;
					double f1 = float(counts[x])/nSeq;
	                double z = (counts[x] - expected) / sqrt(expected*(1-f2));
	                double local_r = f1/f2;
					int position = x-startPos+2;
					if (position < 1) position -= 1;
	                outstream << kmers[i] << "\t" << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << local_r << "\t"  << local_r << endl;
	            }
			}
		}
	}

    return nSig;
} // end of function

// generate paired kmer,
vector<paired_kmer> generate_paired_kmers (
	string alphabet,
int seq1_len,
int seq2_len,
int max_dist,
int min_dist,
int max_shift,
int min_shift
){
	vector<paired_kmer> paired_kmers;
	vector<string> kmers1 = generate_kmers(seq1_len,alphabet);	
	vector<string> kmers2 = generate_kmers(seq2_len,alphabet);	
	for (int i=0;i<kmers1.size();i++)
	{
		for (int j=0;j<kmers2.size();j++)
		{
			for (int d = min_dist; d <= max_dist; d++)
			{
				for (int s= min_shift; s <= max_shift; s++)
				{
					paired_kmer a(kmers1[i],kmers2[j],d, s, 0 , 0);
					paired_kmers.push_back(a);
				}
			}
		}
	}
	return paired_kmers;
}
	

// pairwise
int find_significant_pairs_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<paired_kmer> paired_kmers, 
	string outfile, 
	int nTest, 
	double pCutoff/*=0.05*/, 
	bool Bonferroni/*=true*/,
	int startPos/*=0*/,
	int minCount/*=5*/) 
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// the number of significant positional kmers found
	int nSig = 0;
		
	// will output found significant positional kmers
	ofstream outstream;
	outstream.open(outfile.c_str());

	// start of kmer counting and test
	for( int i=0; i<paired_kmers.size();i++) // for each kmer
	{
		int k = paired_kmers[i].len();
		for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
		{
			vector<bool> positive;
			// whether each sequence is positive
			// for each sequence, 
			for( int j=0;j<nSeq;j++) 
			{	
				bool present = false;	
				for (int offset = 0; offset <= paired_kmers[i].shift; offset ++)
				{
					if (pos + k + offset > lSeq - 1) break;
					if (seqs[j].substr(pos+offset,paired_kmers[i].seq1.size()) == paired_kmers[i].seq1 && seqs[j].substr(pos+offset+paired_kmers[i].dist,paired_kmers[i].seq2.size()) == paired_kmers[i].seq2) 
					{
						present = true;
						break;
					}
				}
				positive.push_back(present);
			}
			vector<double> weights1,weights2;
			for( int j=0;j<nSeq;j++)
			{
				if(positive[j]) weights1.push_back(weights[j]);
				else weights2.push_back(weights[j]);
			}
			// now weights includes weights of all sequences containing kmer i at position pos
			// if less than 3 sequences with/without this kmer, just go on to the next position
			// in such cases the kmer is unlikely to be significant
			// also the t.test will not work well. it requires at least 2 data points in each sample
			if (weights1.size() < minCount || weights2.size() < minCount) 
			{
				continue;
			}
			array<double,6> ttest = t_test(weights1,weights2,false);
			if (ttest[1] < pCutoff)
			{
				// Bonferoni correction
				double pB = min(1.0,ttest[1]*nTest);
				if (Bonferroni && pB > pCutoff) continue;
				nSig++;
				int position = pos-startPos+2;
				if (position < 1) position -= 1;
	            outstream << paired_kmers[i].as_string("_") << "\t" << position << "\t" << paired_kmers[i].shift+k << "\t" << ttest[0] << "\t" << -log10(ttest[1]) << "\t" << -log10(pB) << "\t" << weights1.size() << "\t" << ttest[2] << "\t" << ttest[3] << "\t" << weights2.size() << "\t" << ttest[4] << "\t"  << ttest[5] << endl;
			}
		}
	}
    return nSig;
} // end of function

// nucleotide plot from PKA2 weighted output
// no shift, k=1
// plot column 
void plot_nucleotide_profile(string infile, string outfile, int lSeq, int column, int startPos){
	
	string script = "# R script \n"
	"lSeq="+to_string(lSeq)+" \n"
	"startPos="+to_string(startPos)+"\n"
	"data = numeric(4*lSeq)\n"
	"dict = new.env()\n"
	"dict[[\"A\"]] = 1\n"
	"dict[[\"C\"]] = 2\n"
	"dict[[\"G\"]] = 3\n"
	"dict[[\"T\"]] = 4\n"

	"x=read.table('"+infile+"',header=F)\n"
	"k=nchar(as.character(x[,1]))\n"
	"x[x[,2]<0,2] = x[x[,2]<0,2] + 1\n"
	"x[,2] = x[,2] + startPos - 2\n"
	"x = x[x[,3]==0 & k==1,]\n"
	"pv = x[,"+to_string(column)+"] * ((x[,4]>0)*2-1)\n"

	"pos=numeric(nrow(x))\n"
	"for (i in 1:nrow(x)){\n"
	"	pos[i] = x[i,2]*4+dict[[as.character(x[i,1])]]\n"
	"}\n"
	"data[pos] = pv\n"

	"maxy = max(abs(data))\n"
	"color=c('red','green','blue','yellow') \n"
	"label=c('A','C','G','T')\n"
	"pdf('"+outfile+"',width=10,height=5)\n"
	"bp = barplot(data,col=rep(color,4),ylab='disfavored <- log10(p) -> favored',ylim=c(-maxy,maxy)*1.3,border=NA)\n"
	
	"sub = data != 0 \n"
	"text(bp[sub],data[sub],labels=rep(label,lSeq)[sub],pos= (data[sub]>0)*2+1,offset=0.2,col='gray',cex=0.75) \n"

	"for(i in 1:lSeq){\n"
	"	abline(v=bp[i*4+1]/2+bp[i*4]/2,col='gray',cex=0.1,lty=3)\n"
	"	text(bp[i*4-2]/2+bp[i*4-1]/2,-maxy*1.2,labels=(i-startPos+1 - (i-startPos+1 < 1)),srt=90,cex=0.75)\n"
	"}\n"

	"legend('topright',legend=label,col=color,pch=15,bty='n')\n"

	"dev.off() \n";
	
	R_run(script);
}

void plot_most_significant_kmers(string infile, string outfile, int lSeq, int column, int startPos)
{
	
	string script = "# R script \n"
		
		"dnaColor = new.env()\n"
		"dnaColor[[\"A\"]] = \"red\"\n"
		"dnaColor[[\"C\"]] = \"green\"\n"
		"dnaColor[[\"G\"]] = \"blue\"\n"
		"dnaColor[[\"T\"]] = \"magenta\"\n"

		"color_dna <-function(x,y,txt,col,cex) {  \n"
		"    letters=unlist(strsplit(txt,'')) \n"
		"    thisy<-y \n"
		"    for(txtstr in 1:length(letters)) { \n"
		"        text(x,thisy,letters[txtstr],col=col[[letters[txtstr]]],adj=0,srt=90,cex=cex) \n"
		"        thisy<-thisy+strheight(letters[txtstr],cex=cex)*1.1 \n"
		"    } \n"
		"} \n"
			
	"lSeq="+to_string(lSeq)+" \n"
	"startPos="+to_string(startPos)+"\n"		
    "column = "+to_string(column)+"\n"
	"x=read.table('"+infile+"',header=F)\n"
	"x[x[,2]<0,2] = x[x[,2]<0,2] + 1\n"
	"x[,2] = x[,2] + startPos - 2\n"		
	"pv = x[,column] * ((x[,4]>0)*2-1)\n"
	"pos = x[,2]+1\n"
	"y = numeric(lSeq)\n"
	"labels = character(lSeq)\n"
	"y[pos] = pv\n"
	"labels[pos] = as.character(x[,1])\n"
	"\n"
	"maxy = max(y)\n"	
	"miny = min(y)\n"
	"maxabs = maxy-miny\n"
	"pdf('"+outfile+"',width=10,height=5)\n"
	"bp = barplot(y,col='gray',ylab='disfavored <- log10(p) -> favored',ylim=c(miny-0.5*maxabs,maxy+0.5*maxabs),border=NA)\n"
	
	"#sub = (y != 0) \n"
	"#text(bp[sub],y[sub]+ ((y[sub]>0)*2-1)*(nchar(labels[sub])^0.65)*maxy/20,labels=labels[sub],col='blue',cex=0.75,srt=90) \n"
						
	"for(i in 1:lSeq){\n"
	"   if(y[i] > 0 ){\n"
	"    	color_dna(bp[i],y[i]+0.5*strheight(\"A\",cex=0.75),labels[i],dnaColor,0.75)\n"
	"   }else if(y[i] < 0){\n"
	"    color_dna(bp[i],y[i]-strheight(\"A\",cex=0.75)*(0.5+nchar(labels[i])),labels[i],dnaColor,0.75)\n"	
	"   }\n"
	"	text(bp[i],miny-0.4*maxabs,labels= (i-startPos+1 - (i-startPos+1 < 1)),srt=90,cex=0.75)\n"
	"}\n"

	"#legend('topright',legend=label,col=color,pch=15,bty='n')\n"

	"dev.off() \n";
	
	R_run(script);
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


paired_kmer::paired_kmer(string seqa, string seqb, int distance, int offset, int position, double score)
{
	seq1 = seqa;
	seq2 = seqb;
	dist = distance;
	pos = position;
	shift=offset;
	weight=score;
}

paired_kmer::paired_kmer() 
{
  // allocate variables
	seq1 = "A";
	seq2 = "A";
	dist = 1;
	shift=0;
	pos=0;
	weight=0.0;
}

paired_kmer::paired_kmer(const paired_kmer &a) 
{
  // allocate variables
  paired_kmer();
  // copy values
  operator = (a);
}

// everything equal except weight
bool paired_kmer::equals(paired_kmer a, bool identical_pos, bool identical_shift)
{
	if (identical_pos && pos != a.pos) return false;
	if (identical_shift && shift != a.shift) return false;
	return seq1 == a.seq1 && dist == a.dist && seq2 == a.seq2;
}


const paired_kmer &paired_kmer::operator = (const paired_kmer &a)
{
	seq1 = a.seq1;
	seq2 = a.seq2;
	dist = a.dist;
	pos = a.pos;
	shift = a.shift;
	weight = a. weight;
	return *this;
}

string paired_kmer::as_string(string del/*="_"*/, bool add_shift/*=false*/,bool add_pos/*=false*/, bool add_weight/*=false*/)
{
	string res = seq1+del+to_string(dist)+del+seq2;
	if (add_shift) res += del + to_string(shift);
	if (add_pos) res += del + to_string(pos);
	if (add_weight) res += del + to_string(weight);
	return res;
}

int paired_kmer::len()
{
	return dist+int(seq2.size());
}

int paired_kmer::gap()
{
	return dist - int(seq1.size());
}

// note input should be sorted by weight
vector<positional_kmer> build_model_from_PKA_output(string filename, int startPos)
{
	int cKmer = 0;
	int cStart = 1;
	int cShift = 2;
	int cStat = 3;
	int cWeight = 4; // p, -log10 
		
	vector<string> kmers;
	vector<int> starts, sizes;
	vector<double> weights;
	string kmer;
		
	map<char,string> iupac = define_IUPAC();
	
	set <string> ranked_kmer_ids;
	vector<positional_kmer> ranked_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;

	int group = 0;
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

		if(line[0] == '#') continue;
		
        flds = string_split(line);
		
		double weight = stof(flds[cWeight]);
				
		if(stof(flds[cStat])<0) weight = -weight;
		
		int pos = stoi(flds[cStart]);
		if(pos < 0) pos = pos + 1;
		pos = pos - 2 + startPos;
		
		vector<string> tmp = expand_degenerate_kmer(flds[cKmer],iupac);
		for(int i=0;i<tmp.size();i++)
		{
			string id = tmp[i]+","+flds[cStart]+","+flds[cShift];
			//debug 
			//cout << id << endl;
			// if new kmer instead of the same kmer with lower weight (can happen when degenerate base is used)
			if(ranked_kmer_ids.find(id) == ranked_kmer_ids.end()) 
			{
				ranked_kmer_ids.insert(id);
				positional_kmer x(tmp[i],pos,stoi(flds[cShift])+tmp[i].size(),weight,0);
				/*
				// tried to remove longer but weaker motif, or ignore shorter when stronger present, worse in training
				// try test with crossvalidation
				bool found = false;
				for(int j=0;j<ranked_kmers.size();j++)
				{
					if (ranked_kmers[j].is_part_of(x)) 
					{
						found = true;
						// debug cout << "ignore " << x.as_string() << " given " << ranked_kmers[j].as_string() << endl;
						break;
					}
				}
				if (found) continue;
				// if shorter but weaker, assign the same group number
				found = false;
				for(int j=0;j<ranked_kmers.size();j++)
				{					
					if (x.is_part_of(ranked_kmers[j])) 
					{
						found = true;
						x.group = ranked_kmers[j].group;
						// debug cout << "add " << x.as_string() << " to group " << to_string(x.group) << endl;
						break;
					}
				}
				if (!found) 
				{
					group++;
					x.group = group;
				}
				*/
				ranked_kmers.push_back(x);
			}
		}
	}
	return ranked_kmers;
}

// no degenerate nucleotides
vector<paired_kmer> build_paired_kmer_model(string filename)
{
	int cKmer = 0; // first column is kmer
	int cStart = 1; // second column is start 
	int cStat = 3; // statistics, tells you sign of signficance
	int cWeight = 4; // -log10(p)

	vector<paired_kmer> paired_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

		if(line[0] == '#') continue;
        flds = string_split(line);
		
		string kmer = flds[cKmer];
		int pos = stoi(flds[cStart]);
		double score = stof(flds[cWeight]);
		if(stof(flds[cStat]) < 0 ) score = -score;
		
		vector<string> tmp = string_split(kmer,"_");
		string seq1 = tmp[0];
		string seq2 = tmp[2];
		int dist = stoi(tmp[1]);
		int shift = stoi(tmp[3]);
		
		paired_kmer x(seq1, seq2, dist, shift, pos, score);
		paired_kmers.push_back(x);
	}
	fin.close();
	return paired_kmers;
}

void save_model_to_file(vector<positional_kmer> ranked_kmers, string filename)
{
	ofstream out;
	out.open(filename.c_str());
	
	for( int i=0;i<ranked_kmers.size();i++)
		out << ranked_kmers[i].as_string("\t") << endl;
	
	out.close();
}

vector<positional_kmer> load_model_from_file(string filename)
{
	vector<positional_kmer> ranked_kmers;
	
	ifstream fin;
	fin.open(filename.c_str());
	
	string line;
	vector<string> flds;
	
    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        flds = string_split(line);
		positional_kmer x(flds[0],stoi(flds[1]),stoi(flds[2]),stof(flds[3]),stoi(flds[4]));
		ranked_kmers.push_back(x);
	}
	return ranked_kmers;
}


// ignore group information
double score_sequence_using_PKA_model(vector<positional_kmer> ranked_kmers,  string seq)
{
	double score = 0;
	for( int i=0;i<ranked_kmers.size();i++)
	{
		//cout << i << "\t"<< ranked_kmers[i].as_string() ;
		//cout << "\t" << seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size) ;
		size_t found = seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(ranked_kmers[i].seq);
		if (found != std::string::npos) 
		{
			score += ranked_kmers[i].weight;
			//cout << "\t" << score;
		}
		//cout << endl;
	}
	return score;
}


double score_sequence_using_paired_kmer_model(vector<paired_kmer> model,  string seq)
{
	double score = 0;
	int lSeq = seq.size();
	for( int i=0;i<model.size();i++)
	{
		bool present = false;
		int l = model[i].len();	
		
		for (int offset = 0; offset <= model[i].shift; offset ++)
		{
			if (model[i].pos + l + offset > lSeq - 1) break;
			if (seq.substr(model[i].pos+offset,model[i].seq1.size()) == model[i].seq1 && seq.substr(model[i].pos+offset+model[i].dist,model[i].seq2.size()) == model[i].seq2) 
			{
				present = true;
				break;
			}
		}
		if (present) 
		{
			score += model[i].weight;
			//cout << "\t" << score;
		}
		//cout << endl;
	}
	return score;
}

// use group information, worse
double score_sequence_using_PKA_model_use_group(vector<positional_kmer> ranked_kmers,  string seq)
{
	double score = 0;
	set<int> scored;
	for( int i=0;i<ranked_kmers.size();i++)
	{
		//cout << ranked_kmers[i].as_string() ;
		//cout << "\t" << seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size) ;
		size_t found = seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(ranked_kmers[i].seq);
		if (found!=std::string::npos) 
		{
			if(scored.find(ranked_kmers[i].group) == scored.end())  // if the group was not used before
			{
				score += ranked_kmers[i].weight;
				scored.insert(ranked_kmers[i].group);
				//cout << "\t" << score;
			}
		}
		//cout << endl;
	}
	return score;
}

// for both PKA and PKA2
void score_fasta_using_PKA_model(string seqfile, string outputfile, vector<positional_kmer> ranked_kmers)
{
    ifstream fin(seqfile.c_str());
	ofstream fout(outputfile.c_str());
  
    string name,seq;
    while(fin.good())
    {
		ReadOneSeqFromFasta(fin,name,seq);
     	double score = score_sequence_using_PKA_model(ranked_kmers, seq);
  		fout << name << "\t" << seq << "\t" << score << endl;		
    }
    fin.close();
	fout.close();
}

// sequence in column col, 1 based
void score_tabular_using_PKA_model(string tabfile, int col, string outputfile, vector<positional_kmer> ranked_kmers)
{
	col = col - 1;
	
    ifstream fin(tabfile.c_str());
	ofstream fout(outputfile.c_str());
  
    string line;
	vector<string> flds;
	
    while(fin.good())
    {
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
     	double score = score_sequence_using_PKA_model(ranked_kmers, flds[col]);
  		fout << line << "\t" << score << endl;		
    }
    fin.close();
	fout.close();
}


// count and write feature matrix, format
// seqid, label (i.e. 1/0)
void save_feature_matrix(map<string,string> seqs, vector<string> kmers, string outfile, string label, bool append/*=false*/, int shift/*=0*/)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs.begin()->second.size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	if (append) 
	{
		outstream.open(outfile.c_str(),ios::app);
	} 
	else
	{		
		outstream.open(outfile.c_str());
		// output header
		outstream << "SeqID\t" << label;
		for( int i =0; i< kmers.size();i++)
		{
			for(int j=0;j< lSeq-kmers[i].size()+1; j++)
			{
				outstream << "\t" << kmers[i] << "_" << j << "_" << shift ;
			}
		}
		outstream << endl;
	}


	// start of kmer counting and test
	for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++) // each sequence is a line
	{ 
		outstream << it->first << "\t" << label; 
		for( int i=0; i<kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
			int k = kmers[i].size();
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				//debug cout << kmers[i] << "@" << pos << endl;
				// initialize
				int present = 0;
				//find if any of the expanded kmer is present at position pos
				for( int n=0;n<exp_kmers.size();n++)
				{
					//debug cout << "exact kmer: " << exp_kmers[n] << endl;
					/* speed not affected by shift */
					size_t found = it->second.substr(pos,k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
					if (found!=std::string::npos)
					{
						// add this sequence's rank to sample 1
						//weights1.push_back(weights[j]);
						// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
						present = 1; 
						break;
					} /**/
				}
				outstream << "\t" << present;
			}
		}
		outstream << endl;
	}
} // end of function

// read PKA output into two vectors
void read_significant_positional_kmer_from_file(string inputfile, vector<string> &kmers, vector<int> &positions)
{
	ifstream fin;
	fin.open(inputfile.c_str());
	
	string line;
	vector<string> flds;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		if(line[0] == '#') continue;

        flds = string_split(line);
		
		kmers.push_back(flds[0]);
		positions.push_back(stoi(flds[2]));
	}
	fin.close();
}

// read PKA2 output into two vectors
void read_significant_positional_kmer_from_PKA2_output(string inputfile, vector<string> &kmers, vector<int> &positions)
{
	ifstream fin;
	fin.open(inputfile.c_str());
	
	string line;
	vector<string> flds;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		if(line[0]=='#') continue;

        flds = string_split(line);
		
		kmers.push_back(flds[0]);
		positions.push_back(stoi(flds[1]));
	}
	fin.close();
}


void significant_feature_matrix(map<string,string> seqs, vector<string> kmers, vector<int> positions, string outfile, string label, bool append/*=false*/, int shift/*=0*/)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs.begin()->second.size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	if (append) 
	{
		outstream.open(outfile.c_str(),ios::app);
	} 
	else
	{		
		outstream.open(outfile.c_str());
		// output header
		outstream << "SeqID\tLabel";
		for( int i =0; i< kmers.size();i++)
		{
			outstream << "\t" << kmers[i] << "_" << positions[i] << "_" << shift ;
		}
		outstream << endl;
	}


	// start of kmer counting and test
	for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++) // each sequence is a line
	{ 
		outstream << it->first << "\t" << label; 
		for( int i=0; i<kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
			int k = kmers[i].size();
			
			int present = 0;
			//find if any of the expanded kmer is present at position pos
			for( int n=0;n<exp_kmers.size();n++)
			{
				//debug cout << "exact kmer: " << exp_kmers[n] << endl;
				/* speed not affected by shift */
				size_t found = it->second.substr(positions[i],k+shift).find(exp_kmers[n]); // if found any kmer allowing shift 
				if (found!=std::string::npos)
				{
					// add this sequence's rank to sample 1
					//weights1.push_back(weights[j]);
					// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
					present = 1; 
					break;
				} /**/
			}
			outstream << "\t" << present;
		}
		outstream << endl;
	}
	outstream.close();
} // end of function


// should work for ranked kmer with degenerate nucleotides
void significant_feature_matrix_PKA2(vector<string> seqs, vector<double> weights, vector<positional_kmer> ranked_kmers, string outfile)
{
    int nSeq = seqs.size();		// total number of sequences
	int lSeq = seqs[0].size(); // length of the first sequence, assume all have the same length
	
	//message(to_string(nSeq)+" sequences");
	//message(to_string(lSeq)+" length");
	
	// a map defining IUPAC degenerate nucletodes
	map<char,string> define_iupac = define_IUPAC();

	ofstream outstream;
	outstream.open(outfile.c_str());
	// output header
	outstream << "SeqID\tScore";
	for( int i =0; i< ranked_kmers.size();i++)
	{
		outstream << "\t" << ranked_kmers[i].as_string() ;
	}
	outstream << endl;


	// start of kmer counting and test
	for(int j=0;j<seqs.size();j++) // each sequence is a line
	{ 
		outstream << "Seq-"<< j << "\t" << weights[j]; 
		for( int i=0; i<ranked_kmers.size();i++) // for each kmer
		{
			// expand a degenerate kmer to all possible element/exact kmers		
	        vector<string> exp_kmers = expand_degenerate_kmer(ranked_kmers[i].seq,define_iupac);
			int k = ranked_kmers[i].seq.size();
			
			int present = 0;
			//find if any of the expanded kmer is present at position pos
			for( int n=0;n<exp_kmers.size();n++)
			{
				//debug cout << "exact kmer: " << exp_kmers[n] << endl;
				/* speed not affected by shift */
				size_t found = seqs[j].substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(exp_kmers[n]); // if found any kmer allowing shift 
				if (found!=std::string::npos)
				{
					// add this sequence's rank to sample 1
					//weights1.push_back(weights[j]);
					// stop looking for the next expanded kmer in the same sequence, continue to the next sequence
					present = 1; 
					break;
				} /**/
			}
			outstream << "\t" << present;
		}
		outstream << endl;
	}
	outstream.close();
} // end of function

// PKA: remove overlapping motifs
// input: all significant motifs
int non_overlapping_sig_motifs(string inputfile, string outputfile)
{
    // load file, store pos, kmer, z-score, line number
    // from top
    vector<int> position;
    vector<int> length;

    vector<bool> removed; // indicator whether removed or not

    ifstream fin;
    fin.open(inputfile.c_str());

    string line;
    vector<string> flds;

    while(fin)
    {  
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		if(line[0]=='#') continue;

        flds = string_split(line);
        length.push_back(flds[0].size());
        position.push_back(stoi(flds[2]));
        removed.push_back(false);
    }
    fin.close();
    
    //
    int current_best = 0; // position of currently best kmer
    int next_best = 1;
    while(next_best>0) // if next best can be found
    {
        // keep the current best motif
        removed[current_best] = false;
        next_best = -1;
        // remove from remaining set thsoe overlapping with the top motif
        int start = position[current_best];
        int end = start + length[current_best] - 1;
        for(int i = current_best+1;i<removed.size();i++)
        {
            if (removed[i] == false)
            {
                int start2 = position[i];
                int end2 = start2+ length[i] - 1;
                if ( ((start <= start2) && (end >= start2)) || ((start <= end2) && (end >= end2))) // overlap
                {
                    removed[i] = true;        
                }
                else if (next_best < 0) next_best = i; 
            }
        }
        current_best = next_best;
    }
    
    // output selected lines
    ofstream fout;
    fout.open(outputfile.c_str());
    fin.open(inputfile.c_str());
    getline(fin,line);
    int n = -1;
    int total_selected = 0;
    while(fin)
    {
        getline(fin,line);
        if(line.length() == 0) continue;
        n++;
        if (removed[n] == false) 
        {
            fout << line << endl;
            total_selected++;
        }
    }
    fin.close();
    fout.close();
    return total_selected;
}


// count substring in a map of sequences, i.e. from fasta
// used in generating markov model
int countSubstringInSeqs(vector<string>seqs, string sub)
{
    int count = 0;
    for(int i=0;i<seqs.size();i++)
    {
        count += findall(seqs[i],sub).size();
    } 
    return count;
}



set<int> findall(string seq, string motif)
{
    // find all match positions of motif in seq
    // allowing overlap 
    set<int> allpos; // to store all found positions
    size_t pos = 0; //start at pos
    pos = seq.find(motif); // the first match
    while(pos != string::npos)
    {
        allpos.insert(pos);
        pos = seq.find(motif,pos+1); // start search for the next match from pos + 1
    }
    return allpos;
}


// dict for IUPAC degenerate nucleotides
map<char,string> define_IUPAC()
{
    map<char,string> m;

    m['A']="A";
    m['C']="C";
    m['G']="G";
    m['T']="T";
    m['U']="T";
    m['R']="AG";
    m['Y']="CT";
    m['M']="AC";
    m['K']="GT";
    m['W']="AT";
    m['S']="CG";
    m['B']="CGT";
    m['D']="AGT";
    m['H']="ACT";
    m['V']="ACG";
    m['N']="ACGT";

    return m;
}

/*
// for output only
map<char,string> interpret_IUPAC()
{
    map<char,string> m;

    m['A']="A";
    m['C']="C";
    m['G']="G";
    m['T']="T";
    m['U']="T";
    m['R']="AG";
    m['Y']="CT";
    m['M']="AC";
    m['K']="GT";
    m['W']="AT";
    m['S']="CG";
    m['B']="a";
    m['D']="c";
    m['H']="g";
    m['V']="t";
    m['N']="N";

    return m;
}
*/

// generate all possible combinations of letters in alphabet of fixed length k, i.e. kmer
vector<string> generate_kmers(int k, string alphabet)
{
    vector<string> kmers;
    
    // first position, each letter in alphabet
    for(int i = 0; i< alphabet.size(); i++)
    {
        string s(1,alphabet[i]);
        kmers.push_back(s);
    }

    // for the rest k-1 positions, fill all possible letters 
    for(int i = 0; i< k-1;i++)
    {
        int n = kmers.size();
        // for each current kmer, append all possible letters
        for(int j=0; j< n; j++)
        {
            string tmp = kmers[j];
            kmers[j] = kmers[j]+alphabet[0];
            for( int m =0;m< alphabet.size()-1;m++)
            {
                string s(1,alphabet[m+1]);
                kmers.push_back(tmp+s);
            }
        }
    }
    return kmers;
}

// degenerate kmers based on IUPAC, DNA only
// remove those with terminal N, which is not a kmer, but (k-i)mer
vector<string> degenerate_kmer(int k, string alphabet/*ACGTRYMKWSBDHVN*/)
{
    vector<string> tmp = generate_kmers(k,alphabet);
    // remove those with terminal N
    vector<string> dkmers;
    for( int i = 0; i < tmp.size(); i++)
    {  
        if (tmp[i][0] != 'N' && tmp[i][k-1] != 'N') dkmers.push_back(tmp[i]);
    }
    return dkmers;
}

// convert one degenerate kmer to regular expression, e.g. ARG -> A[AG]G
string degenerate_kmer_to_regex(string kmer,map<char,string> iupac)
{
    string str;
    for( int i  = 0; i< kmer.size();i++)
    {
        if (iupac[kmer[i]].size() > 1) str = str + "["+iupac[kmer[i]]+"]";
        else str = str + iupac[kmer[i]];
    }    
    return str;
}

// expand a degenerate kmer to all possible exact kmers
vector<string> expand_degenerate_kmer(string seq, map<char,string> iupac)
{
    vector<string> res (1,"");
    vector<string> tmp;

    for( int i=0;i<seq.size();i++)
    {  
        string exp = iupac[seq[i]];
        int L = res.size();
        tmp.clear();
        for( int j=0;j<L;j++)
        {  
            for( int k=0;k<exp.size();k++)
            {  
                tmp.push_back(res[j]+exp[k]);
            }
        }
        res = tmp;
    }
    return res;
}

// implant a motif to a set of sequences
// motif can contain degenerate nucleotides
void implant_motif(map<string,string> &seqs, int position, string motif, double fraction)
{
    int k = motif.size();
    int nSeq = seqs.size() * fraction; // number of sequences to be implanted
    // expand degenerate nucleotides
    vector<string> kmers = expand_degenerate_kmer(motif, define_IUPAC());

    // implant into the first nSeq
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {
        it->second.replace(position,k,kmers[rand()%kmers.size()]);
        nSeq--;
        if (nSeq==0) return ;     
    }
}


// build position weight matrix for a positional kmer, can include flanking sequence
// haven't conder startPos, assume 1 start
boost::numeric::ublas::matrix<double> position_weight_matrix_from_PKA_output(string filename, int seqLen, int startPos, int cScore)
{
	int cSeq = 0;
	int cPos = 1; 
	int cStat = 3;
 	cScore -= 1;
	
	// initialize an empty matrix
    boost::numeric::ublas::matrix<double> pwm (4,seqLen);
    for (int i = 0; i < pwm.size1(); ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<string,int> letter2pos;
    letter2pos["A"] = 0;
    letter2pos["C"] = 1;
    letter2pos["G"] = 2;
    letter2pos["T"] = 3;
 
 	ifstream fin(filename.c_str());
	string line;
	vector<string> flds;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		if(flds[cSeq] == "A" || flds[cSeq] == "C" || flds[cSeq] == "G" || flds[cSeq] == "T")
		{
			double score = stof(flds[cScore]);
			if (stof(flds[cStat]) < 0) 
			{
				if (cScore != cStat) score = -score;
			}
			int pos = stoi(flds[cPos]);
			if (pos < 0) pos += 1;
			pos = pos + startPos - 2;
			//cout << pos << "," << line << endl;
			pwm(letter2pos[flds[cSeq]],pos) = score;
		}
	}
	fin.close();

    return pwm;
}



// one text all the same color
string postscript_text(string text, double x, double y, double width, double height, string color, double rotate)
{
	string res = "\n"
    "gsave \n "
    "  "+to_string(x)+" "+to_string(y)+" moveto \n"
    "  "+to_string(width) +" "+to_string(height) + " scale \n "
	"  "+to_string(rotate)+ " rotate\n"
    "  "+color + " setrgbcolor ("+text+") dup stringwidth pop 2 div neg 0 rmoveto show  \n"
    "grestore \n"
	"\n";	
	
	return res;
}

string postscript_line(double x1, double y1, double x2, double y2)
{
	string res = "\n"
	"gsave \n"
	"     newpath \n"
	"     " + to_string(x1) + " " + to_string(y1) +" moveto \n"
	"     " + to_string(x2) + " " + to_string(y2) +" lineto \n"
	"stroke \n\n";
	
	return res;
}

void generate_ps_logo_from_pwm(boost::numeric::ublas::matrix<double> pwm, string filename, string alphabet, map<char,string> colors, double score_cutoff, int startPos, int fontsize, string y_label, bool information_content, double max_scale)
{
	// fontsize
	//int fontsize = 20;
	
	// the actual size of each character is smaller than fontsize
	double scalex = 0.75;
	double scaley = 0.75;
	
	// max score will be used to normalize
	// this will be the max value in pwm after normalization
	// double max_scale = 6;
	
	// distance between pos and neg, for plotting coordinates
	
	double xstep = fontsize * scalex;
	double ystep = fontsize * scaley;
	double coord_size = fontsize * scaley * 2;
	
	int L = pwm.size2();
	
	if (information_content) 
	{
		vector<double> info_content = matrix_column_information_content(pwm);
		for (int i=0;i<pwm.size1();i++)
		{
			for (int j=0;j<pwm.size2();j++)
			{
				pwm(i,j) *= info_content[j];
			}
		}
	}

	double maxv = matrix_max(pwm);
	double minv = matrix_min(pwm);
	vector<double> col_max = matrix_column_max(pwm);
	vector<double> col_min = matrix_column_min(pwm);
		
	vector<bool> significant;
	for (int i=0;i<L;i++)
	{
		if (col_max[i] > score_cutoff || col_min[i] < -score_cutoff) significant.push_back(true);
		else significant.push_back(false);
	}
	

	
	double absmax = max(maxv,-minv);
		
	pwm = pwm / absmax * max_scale;
	
	vector<double> colum_sum_pos = matrix_column_sum(pwm, 1);
	
	vector<double> colum_sum_neg = matrix_column_sum(pwm, -1);
	//for (int i=0;i<colum_sum_neg.size();i++) cout << colum_sum_neg[i] << ",";
	//cout << y_label << endl;

	double height_pos = max(colum_sum_pos)  ;
	double height_neg = max(0.0,-min(colum_sum_neg)) ;	
	
    //print_matrix(pwm);
	
	//cout << height_pos << endl;
	//cout << height_neg << endl;
	
	
	
	// origin in the plot: left, center
	int x0 = 5 * xstep;
	int y0 = (height_neg + 2) * ystep;
	
	// ps file header, define plot window
	string header = ""
		"%!PS-Adobe-3.0 EPSF-3.0            \n"
		"%%Title: Sequence Logo : Logo      \n"
		"%%Creator: PKA 					\n"
		"%%CreationDate: "+current_time()+" \n"
		"%%BoundingBox:   0  0  "+to_string(L * xstep + x0 * 1.5)+" "+to_string( (height_pos+height_neg + 6) * ystep )+" \n"
		"%%Pages: 0                         \n"
		"%%DocumentFonts:                   \n"
		"%%EndComments                      \n"
		"\n"
		"/Helvetica-Bold findfont "+to_string(fontsize)+" scalefont setfont\n";
	
	ofstream out(filename.c_str());
	out << header << endl;
	
    // draw y axis, +
    out << postscript_line(x0-xstep*1.5, y0+coord_size,  x0-xstep*1.5, y0 + coord_size + height_pos * ystep);
    //  out << postscript_text("0", x0-xstep, y0+ystep/2, 0.6,0.6, "0 0 0", 0);
    double y_max = height_pos * absmax / max_scale;
    out << postscript_text(to_string_with_precision(y_max), x0 - xstep*1.5 , y0+coord_size + (height_pos+0.4) * ystep, 0.6,0.6, "0 0 0", 0);
    // draw y axis, -
    if(minv < 0)
    {
        out << postscript_line(x0-xstep*1.5, y0,  x0-xstep*1.5, y0 - height_neg * ystep);
        //out << postscript_text("0", x0-xstep*1.5, y0-ystep/2, 1, 1, "0 0 0", 0);
        double ymax = height_neg* absmax / max_scale;
        out << postscript_text(to_string_with_precision(ymax), x0-xstep*1.5, y0-(height_neg+1) * ystep , 0.6,0.6, "0 0 0", 0);
        out << postscript_text(y_label, x0-xstep*3, y0+coord_size/2+(height_pos - height_neg)*ystep/2, 0.6, 0.6, "0 0 0", 90);
    }
    // yaxis label
	else out << postscript_text(y_label, x0-xstep*3, y0+coord_size+(height_pos - height_neg)*ystep/2, 0.6, 0.6, "0 0 0", 90);

	
	// x axis
	for(int i = 0 ; i < L; i++) 
	{
		int pos = i - startPos + 2;
		if (pos < 1) pos -= 1;
		string color_sig = "0.5 0.5 0.5";
		if (significant[i]) color_sig = "1 0 0";
		out << postscript_text(to_string(pos),x0 + xstep * (i+0.3) ,y0+ystep,0.6,0.6,color_sig,90);
	}
	
	// plot logo at each position
	for (int i=0;i<L; i++)
	{
		double x = x0 + xstep * i;
		double y_pos = y0+coord_size;
		double y_neg = y0;
		double y;

		// sort and return index, ascending order
		vector<double> w;
		for(int j=0;j<alphabet.size();j++) w.push_back(pwm(j,i));
		vector<size_t> idx = sort_indexes(w);
		//for(int j=0;j<4;j++) cout << i << "\t" << letters[idx[j]] << "\t"<< idx[j] << "\t" << w[idx[j]] << endl;
		
		for(int j=0;j<alphabet.size();j++)
		{
			if (w[idx[j]] > 0)
			{		
                string color;
                if (colors.find(alphabet[idx[j]]) != colors.end()) {color = colors[alphabet[idx[j]]];}
                else {color = "0.5 0.5 0.5";}
                 	
				out << postscript_text(alphabet.substr(idx[j],1),x, y_pos, 1, w[idx[j]],color); 
				y_pos += w[idx[j]] * ystep ;
			}
		}
		for(int j=alphabet.size()-1;j>=0;j--)		
		{
			if (w[idx[j]] < 0)
			{
				y_neg += w[idx[j]] * fontsize * scaley;
                string color;
                if (colors.find(alphabet[idx[j]]) != colors.end()) {color = colors[alphabet[idx[j]]];}else {color = "0.5 0.5 0.5";}
				out << postscript_text(alphabet.substr(idx[j],1), x, y_neg, 1, -w[idx[j]],color); 
				//y_neg += w[idx[j]] * ystep * 0.05;
			}
		}
		//cout << i << "," << w[idx[0]] * absmax / max_scale << "," << w[idx[3]] * absmax / max_scale<< "," << score_cutoff << endl;
		//if (w[idx[3]] * absmax / max_scale > score_cutoff)	out << postscript_text("*",x+0.25*xstep,max(y_pos+0.5*ystep, y0+coord_size + int(height_pos) * ystep),0.6,0.6,"1 0 0",0);
		//if (w[idx[0]] * absmax / max_scale < -score_cutoff) out << postscript_text("*",x+0.25*xstep,min(y_neg-0.5*ystep, y0-int(height_neg) * ystep) ,0.6,0.6,"1 0 0",0);
		
	}

	out << "showpage" << endl;
	out << "%%EOF" << endl;
	
	out.close();
						
}


// vertical kmer, color each bases differently	  
// height: total height, can be negative
string postscript_kmer(string seq, double x, double y, int fontsize, double scaley, double width, double height, map<char,string> colormap, double rotate)
{
	// length of the sequence
	int L = seq.size();
	
	// determine height for each base
	if (height < 0) 
	{
		y = y + height * fontsize * scaley;
		height = - height / L;
	}
	else height = height / L;
	
	string res,color;
	
	for(int i = L-1;i >= 0;i--)
	{
		if (colormap.find(seq[i]) == colormap.end()) color = "0.5 0.5 0.5";
        else color = colormap[seq[i]];
		//cout << i << endl;
		res += "\n"
   		 "gsave \n "
    	 "  "+to_string(x)+" "+to_string(y + fontsize * scaley * (L-i-1) * height)+" moveto \n"
         "  "+to_string(width) +" "+to_string(height*0.95) + " scale \n "
	     "  "+to_string(rotate)+ " rotate\n"
         "  "+color + " setrgbcolor ("+seq[i]+") dup stringwidth pop 2 div neg 0 rmoveto show  \n"
         "grestore \n"
	     "\n";
	}	
	return res;
}


void postscript_logo_from_PKA_output(string infile, string outfile, map<char,string> colors, int seqLen, double score_cutoff, int startPos, int fontsize, int cScore,string y_label, double max_scale)
{
	int cSeq = 0;
	int cPos = 1; 
	int cStat = 3;
 	cScore -= 1;
	
	int L = seqLen;
	vector<bool> significant;
	for(int i=0;i<L;i++) significant.push_back(false);
	
	// determine the max / min scores
	double maxv = -1;	
	double minv = 1e100;
	
 	ifstream fin(infile.c_str());
	string line;
	vector<string> flds;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		double score = stof(flds[cScore]);
		int pos = stoi(flds[cPos]);
		if (pos < 0) pos += 1;
		pos = pos + startPos - 1;
		if (score > score_cutoff) significant[pos-1] = true;
		if(stof(flds[cStat]) < 0 ) 
		{
			if (cScore != cStat) score = -score;
		}
		if (score > maxv) maxv = score;
		if (score < minv) minv = score;

	}
	fin.close();
	
	// fontsize
	//int fontsize = 20;
	
	double scalex = 0.75;
	double scaley = 0.75;
	

	double xstep = fontsize * scalex;
	double ystep = fontsize * scaley;
	double coord_size = fontsize * scaley * 2;
	
	double absmax = max(maxv,-minv);
				
	double height_pos = maxv  / absmax * max_scale ;
	double height_neg = max(0.0,-minv/absmax * max_scale) ;	
	
	//cout << minv << endl;
	//cout << height_pos << endl;
	//cout << height_neg << endl;
	

	
	// origin in the plot: left, center
	int x0 = 5 * xstep;
	int y0 = (height_neg + 1.5) * ystep;
	
	// ps file header, define plot window
	string header = ""
		"%!PS-Adobe-3.0 EPSF-3.0            \n"
		"%%Title: Sequence Logo : Logo      \n"
		"%%Creator: PKA 					\n"
		"%%CreationDate: "+current_time()+" \n"
		"%%BoundingBox:   0  0  "+to_string(L * xstep + x0 * 1.5)+" "+to_string( (height_pos+height_neg + 5) * ystep )+" \n"
		"%%Pages: 0                         \n"
		"%%DocumentFonts:                   \n"
		"%%EndComments                      \n"
		"\n"
		"/Helvetica-Bold findfont "+to_string(fontsize)+" scalefont setfont\n";
	
	ofstream out(outfile.c_str());
	out << header << endl;
	
	// draw y axis, +
	out << postscript_line(x0-xstep*1.5, y0+coord_size,  x0-xstep*1.5, y0 + coord_size + height_pos * ystep);
    //	out << postscript_text("0", x0-xstep, y0+ystep/2, 0.6,0.6, "0 0 0", 0);
	double y_max = height_pos * absmax / max_scale;
	out << postscript_text(to_string_with_precision(y_max), x0 - xstep*1.5 , y0+coord_size + (0.4+height_pos) * ystep, 0.6,0.6, "0 0 0", 0);
	// draw y axis, -
	if(minv < 0)
	{
		out << postscript_line(x0-xstep*1.5, y0,  x0-xstep*1.5, y0 - height_neg * ystep);
		//out << postscript_text("0", x0-xstep*1.5, y0-ystep/2, 1, 1, "0 0 0", 0);
		double ymax = height_neg* absmax / max_scale;
	    out << postscript_text(to_string_with_precision(ymax), x0-xstep*1.5, y0-(height_neg+1) * ystep , 0.6,0.6, "0 0 0", 0);
	    out << postscript_text(y_label, x0-xstep*3, y0+coord_size/2+(height_pos - height_neg)*ystep/2, 0.6, 0.6, "0 0 0", 90);
	}
	// yaxis label
	else out << postscript_text(y_label, x0-xstep*3, y0+coord_size+(height_pos - height_neg)*ystep/2, 0.6, 0.6, "0 0 0", 90);

	// x axis
	for(int i = 0 ; i < L; i++) 
	{
		int pos = i - startPos + 2;
		if (pos < 1) pos -= 1;
		string color_sig = "0.5 0.5 0.5";
		if (significant[i]) color_sig = "1 0 0";
		out << postscript_text(to_string(pos),x0 + xstep * (i+0.3) ,y0+ystep,0.6,0.6,color_sig,90);
	}
	
	
	fin.open(infile.c_str());
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		double score = stof(flds[cScore]) ; 
		if(stof(flds[cStat]) < 0 ) 
		{
			if (cScore != cStat) score = -score;
		}
		int pos = stoi(flds[cPos]);
		if (pos < 0) pos += 1;
		pos = pos + startPos - 1;
		//cout << flds[cSeq] << endl;
		if (score > 0 ) out << postscript_kmer(flds[cSeq], x0+ xstep * (pos-1), y0+coord_size, fontsize, scaley, 1, score/ absmax * max_scale, colors, 0);
		else if (score < 0) out << postscript_kmer(flds[cSeq], x0+ xstep * (pos-1), y0, fontsize, scaley, 1, score/ absmax * max_scale, colors, 0);
		//cout << pos << "," <<  score_cutoff << "," << score<< endl;
		//if (score > score_cutoff)	out << postscript_text("*",x0+ xstep * (pos-0.75),y0+coord_size + int(height_pos) * ystep,0.6,0.6,"1 0 0",0);
		//else if (score <  -score_cutoff) out << postscript_text("*",x0+ xstep * (pos-0.75),y0-int(height_neg) * ystep ,0.6,0.6,"1 0 0",0);
			
	}
	fin.close();
	
	out << "showpage" << endl;
	out << "%%EOF" << endl;
	
	out.close();
			
}


// write pwm in meme format
void write_pwm_in_meme_format(boost::numeric::ublas::matrix<double> pwm, string motifname, string filename)
{
    ofstream out;
    out.open(filename.c_str());
    out << "MEME version 4.4\n\n";
    out << "ALPHABET= ACGT\n\n";
    out << "strands: + -\n\n";
    out << "Background letter frequencies\n";
    out << "A 0.25 C 0.25 G 0.25 T 0.25" << "\n\n";
    out << "MOTIF "<<motifname<<"\n\n";
    out << "letter-probability matrix: alength= 4 w= "<<pwm.size2()<<" nsites= 1000 E= 0"<<endl;
    for (int i = 0; i < pwm.size2(); ++ i)
    {
        for (int j = 0; j < pwm.size1(); ++ j)
        {
            out << pwm(j,i) << "\t";
        }
        out << endl;
    }
    out << endl;
    out.close();
}


// build position weight matrix from a set of sequences of the same length
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_seqs(vector<string> seqs)
{

    // length of sequences
    int L = seqs[0].size();

    // initialize an empty matrix
    boost::numeric::ublas::matrix<double> pwm (4,L);
    for (int i = 0; i < 4; ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<char,int> letter2pos;
    letter2pos['A'] = 0;
    letter2pos['C'] = 1;
    letter2pos['G'] = 2;
    letter2pos['T'] = 3;
 
    for(int k=0;k<seqs.size();k++)
    {
        for( int i =0;i<L;i++)
        {
        	pwm(letter2pos[seqs[k][i]],i) += 1;
        }
    }
	
	pwm = pwm / seqs.size();
	
    return pwm;
}


// build position weight matrix for a positional kmer, can include flanking sequence
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_kmer(vector<string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos)
{
    //position: 0-based
    // need to back-calculate position = position + startPos
    position = position + startPos; 
    
    //cout << kmer << "\t" << position;     //debug

    // length of sequences
    int seqLen = seqs[0].size();
    int k = kmer.size();

    // begin and end of region in the sequence
    int start = max(0,position-d);
    int end = min(seqLen-1,position + k + d - 1);

    // initialize an empty matrix
    int L = end - start + 1;
    boost::numeric::ublas::matrix<double> pwm (4,L);
    for (int i = 0; i < pwm.size1(); ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<char,int> letter2pos;
    letter2pos['A'] = 0;
    letter2pos['C'] = 1;
    letter2pos['G'] = 2;
    letter2pos['T'] = 3;
 
    // expand kmer
    vector<string> dkmersexp = expand_degenerate_kmer(kmer,iupac); 

    // find sequence with motif match and update matrix       
    int total_positive_seq = 0; // total number of sequences contain the motif
    for(int k=0;k<seqs.size();k++)
    {
        for( int i =0;i<dkmersexp.size();i++)
        {
            if (seqs[k].substr(position,k) == dkmersexp[i])
            {
                // update matrix
                for( int j=start;j<=end;j++)
                {
                    pwm(letter2pos[seqs[k][j]],j-start) += 1;
                }                
                total_positive_seq ++;
                break;
            }
        }
    }

    //cout << "matrix done" << endl;


    // divide
    pwm /= total_positive_seq;

    //debug
    //cout << "after division" << endl;
    //cout << kmer << "," << total_positive_seq << endl;
    //print_matrix(pwm);
    //    cout << kmer << "\t" << total_positive_seq << endl;

    return pwm;
}

//PKA : create logo for a single kmer
void create_logo_for_kmer(vector<string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos,string output)
{
    boost::numeric::ublas::matrix<double> pwm;
    pwm = create_position_weight_matrix_from_kmer(seqs, kmer, position,iupac,d,startPos);
    
    string meme_filename = output+"_"+kmer+"_"+to_string(position)+".meme"; // file name use 1-based position
    write_pwm_in_meme_format(pwm,kmer,meme_filename);

    string cmd = "ceqlogo -i "+meme_filename+" -o "+meme_filename+".eps -t "+meme_filename+" -x '' -b 2 -c 1 -w 10 -h 5";  
    system(cmd.c_str());
    cmd = "ps2pdf -dEPSCrop "+meme_filename+".eps "+meme_filename+".pdf";
    system(cmd.c_str());
}

// create logos for top 1 kmer passing pCutoff_B
void create_logo_for_topN_sig_kmer_per_position(vector<string> seqs, string filename,int d, double pCutoff,bool Bonferroni,int startPos,string output)
{
    ifstream fin;
    fin.open(filename.c_str());
    string line;
    vector<string> flds;

    map<char,string> iupac = define_IUPAC();


    while(fin)
    {  
        getline(fin,line);
        if (line.length() == 0)
            continue;
		
		if(line[0] == '#') continue;
        flds = string_split(line);
        //cout << flds[9] <<","<< pCutoff << endl;
        if(stod(flds[9]) < pCutoff)
        {
            create_logo_for_kmer(seqs,flds[0],stoi(flds[2]),iupac,d,startPos,output); // change back to 0-based coordinates
        }
    }    
    fin.close();
}

/*
vector<int> count_one_kmer_in_all_seqs_regex(map<string,string> seqs, boost::regex pattern, int k, int seq_len_minus_k_plus_1, int k_plus_shift)
{
    vector<int> counts;
    for(int j=0;j<seq_len_minus_k_plus_1;j++) counts.push_back(0);
    for(map<string,string>::iterator it=seqs.begin();it!=seqs.end();it++)
    {   
        string seq = (*it).second+"@@@@@@@@@@@@@@";
         // for each position + shift window
        for( int m=0;m<seq_len_minus_k_plus_1;m++)
        {   
            if (boost::regex_search(seq.substr(m,k_plus_shift),pattern)) counts[m]++;
        }
    }
    return counts;
}
*/

/*
// find significant degenerate kmers with shift
array<int,2> find_significant_degenerate_kmer_with_shift(vector<string> kmers, map<string,string> seqs1, map<string,string> seqs2, int shift, double pCutoff, double pCutoff_B,float pseudo, int startPos, int nTest, string output)
{
    map<char,string> define_iupac = define_IUPAC();
    //map<char,string> interpret_iupac = interpret_IUPAC();

    int seq_len = seqs1.begin()->second.size();

    int nSeq1 = seqs1.size();
    int nSeq2 = seqs2.size();

    // significant kmers
    array<int,2> nSig  = {0,0};

    ofstream fout;
    fout.open(output.c_str());
    fout.precision(3);

    for( int i=0;i<kmers.size();i++)
    { 
        // each kmer could be of different size
        int k = kmers[i].size();
        // convert each kmer to pattern
        string regex_format = degenerate_kmer_to_regex(kmers[i],define_iupac);
        boost::regex pattern(regex_format); 
        vector<int> counts1 = count_one_kmer_in_all_seqs_regex(seqs1,pattern, k, seq_len - k + 1, k + shift);
        vector<int> counts2 = count_one_kmer_in_all_seqs_regex(seqs2,pattern, k, seq_len - k + 1, k + shift);

        //
        int total_counts1 = sum(counts1);

        for( int m=0; m < seq_len - k + 1;m++)
        {   
            //if (counts1[m] == 0) continue;
            double f2 = float(counts2[m]+pseudo)/nSeq2;
            double p = binom_test(nSeq1,counts1[m],f2);
            if (p < pCutoff)
            {
                nSig[0]++;
                double f1 = float(counts1[m])/nSeq1;
                double expected = nSeq1 * f2;
                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
                fout  << kmers[i] << "\t" << regex_format << "\t"  << m-startPos <<"\t"<< p << "\t"<< min(1.0,p*nTest)<< "\t"<< z << "\t" << f1/f2 << "\t" << f1 << "\t" << f2 << "\t" << local_r << endl;
                if (p*nTest<pCutoff) nSig[1]++;
            }
        }
    }

    return nSig;
}
*/

// doesn't work with shift+degenerate
// for each kmer, count its freq at each position(start), allowing shifts (to the right)
map<string,vector<int> > count_all_kmer_in_seqs(vector<string> kmers, vector<string> seqs, int shift)
{
    // kmers: vector of all kmers
    // 
    map<string, vector<int> > data;

    int i,j,k;
	int m;
    string name,seq;
    vector<int> counts;
    set<int> found;

    int seq_len = seqs[0].size();

    for (i=0;i<kmers.size();i++)
    {  
        // kmers could be of different size
        k = kmers[i].size();
        // initialize counts
        for(j=0;j<seq_len-k+1;j++) counts.push_back(0);
        // loop through each sequencs, count this kmer
        for(int x=0;x<seqs.size();x++)
        {   
            seq = seqs[x];

            // find all matches
            found = findall(seq,kmers[i]);
            // add shifts
            vector<int> to_be_inserted;
            for (set<int>::iterator it=found.begin(); it!=found.end(); ++it)
            {   
                for (m=1;m<shift+1;m++)
                {   
                    if( *it > m) to_be_inserted.push_back(*it-m);
                }
            }
            for (m=0;m<to_be_inserted.size();m++)  found.insert(to_be_inserted[m]);
            for (set<int>::iterator it=found.begin(); it!=found.end(); ++it) counts[*it]++;
        }
        data[kmers[i]] = counts;
        counts.clear();
    }
    return data;
}   

//PKA
void print_kmer_positional_profile(map<string,vector<int> > data)
{   
    for (map<string,vector<int> >::iterator it=data.begin();it!=data.end();it++)
    {
        cout << (*it).first;
        for( int i=0;i<(*it).second.size();i++) cout << ',' << (*it).second[i];
        cout << endl;
    }
}


map<string,vector<int> > degenerate_kmer_counts(vector<string> dkmers,map<string,vector<int> > data, map<char,string> define_iupac)
{
    map<string,vector<int> > new_data;
    vector<int> tmp;
   // for each degenerate kmer, combine element kmer counts
    for( int i=0;i<dkmers.size();i++)
    {  
        vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);
        new_data[dkmers[i]] = data[dkmersexp[0]];
//        cout << i << "\t" << dkmers[i] << "\t" << dkmersexp.size() << endl;
        for( int j=1;j<dkmersexp.size();j++)
        {
  //          cout << i << "," << j << endl;
            tmp = sum(new_data[dkmers[i]],data[dkmersexp[j]]);
            new_data[dkmers[i]] = tmp;
    //        cout << "done" << endl;
        }
    }
    return new_data;
}

int find_significant_kmer_from_one_seq_set(
	vector<string>seqs1, 
	map<string,double> probs_kmer,vector<string>kmers, 
	vector<string> dkmers, 
	int min_shift,
	int max_shift,
	bool degenerate,
	double pCutoff, 
	bool Bonferroni/*=true*/,
	int startPos,
	int nTest, 
	string outfile, 
	string output_count_file)
{
    int nSeq1 = seqs1.size();
    int nSig = 0;

    ofstream outstream;
    outstream.open(outfile.c_str());

    ofstream outcounts;
    outcounts.open(output_count_file.c_str());


	for(int shift = min_shift; shift <= max_shift; shift ++ )
	{
	    // kmer counts in foreground
	    // note that each count vector could be of different length due to kmer length difference
	    message("counting exact kmers in foreground sequences...shift="+to_string(shift));
	    map<string, vector<int> > data1 = count_all_kmer_in_seqs(kmers, seqs1, shift);

	    // compute p-value for exact kmers
	    message("computing p-values for exact kmers...");
	    for (map<string,vector<int> >::iterator it=data1.begin(); it!=data1.end(); ++it) 
	    {
	        // counts in foreground
	        vector<int> counts1 = it->second;

	        // total counts 
	        int total_counts1 = sum(counts1);

	        double f2 = min(1.0,probs_kmer[it->first]*(1+shift));


	        // for each position
	        for( int m=0;m<counts1.size();m++)
	        { 
	            if (counts1[m] == 0 && f2 ==0) continue; 
	            // compare background p estimate from shuffling and markov model
	            double p = binom_test(nSeq1,counts1[m],f2);
	            if (p < pCutoff)
	            {
					if(p == 0) p = 1e-16;
	                double corrected_p = min(1.0,p*nTest);					
					if(Bonferroni && corrected_p > pCutoff) continue;
	                nSig++;
	                double f1 = float(counts1[m])/nSeq1;
	                double expected = nSeq1 * f2;
	                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
	                // local enrichment
	                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
										
	                outstream << it->first << "\t" << m-startPos << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	                outcounts << it->first << ":" << m-startPos;
	                for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
	                outcounts << endl;
	            }
	        }
	    }

	    // compute p-value for non-exact kmers, only work without shift
	    if (degenerate)
	    {   

	        message( to_string( nSig)+ " significant exact kmers identified");
	        message( "computing p-values for degenerate kmers...");

	        // map from degeneate bases to all allowed bases, e.g. R => AG
	        map<char,string> define_iupac = define_IUPAC();
	        // for output purpurse only, e.g. will use g to replace H = [ACT] = non-G
	        // map<char,string> interpret_iupac = interpret_IUPAC();

	        // for each degenerate kmer, combine element kmer counts
	        for( int i=0;i<dkmers.size();i++)
	        {
	            // if already in exact kmers, i.e. contain no degenerate bases, skip
	            if ( find(kmers.begin(), kmers.end(), dkmers[i])!=kmers.end() ) continue;

	            // expand a degenerate kmer to all possible element/exact kmers
	            vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);

	            // for output only, also generate regex version of the degenerate kmer
	            string dkmersexp_readable = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // initialize the count at each position with the first element kmer
	            vector<int> counts1  = data1[dkmersexp[0]];
	            // sum over the rest element kmers
	            for( int j=1;j<dkmersexp.size();j++) counts1 = sum(counts1,data1[dkmersexp[j]]);

	            // calculate probs from markov model: sum of element prob
	            double f2 = probs_kmer[dkmersexp[0]];
	            for( int j=1;j<dkmersexp.size();j++) f2 += probs_kmer[dkmersexp[j]]; 
	            //f2 = f2 * (1+shift); // can't have shift on degenerate kmers

	            // for output only
	            string regex_format = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // total counts
	            int total_counts1 = sum(counts1);

	            // calculate  p-value 
	            for( int m=0;m<counts1.size();m++)
	            { 
	                if (counts1[m] == 0 && f2==0) continue;
	                double p = binom_test(nSeq1,counts1[m],f2);
	                if (p < pCutoff)
	                {
						if(p == 0) p = 1e-16;
	                    double corrected_p = min(1.0,p*nTest);
						if(Bonferroni && corrected_p > pCutoff) continue;
						
	                    nSig++;
	                    double f1 = float(counts1[m])/nSeq1;
	                    double expected = nSeq1 * f2;
	                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
	                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
												
	                    outstream  << dkmers[i] << "\t" << m-startPos << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	        
                       outcounts << dkmers[i] << ":" << m-startPos;
                       for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
                       outcounts << endl;
	                }
	            }
	        }
	    }
	}
	

    outstream.close();

    outcounts.close();

    return nSig;
} // end of function


// two file comparison, not allow shift and degenerate at the same time
int find_significant_kmer_from_two_seq_sets(vector<string>seqs1, vector<string>seqs2, vector<string>kmers, vector<string> dkmers, int min_shift,int max_shift,bool degenerate,double pCutoff, bool Bonferroni/*=true*/,double pseudo,int startPos,int nTest, string outfile,string output_count_file)
{
    int nSeq1 = seqs1.size();
    int nSeq2 = seqs2.size();
    int nSig = 0;



    ofstream outstream;
    outstream.open(outfile.c_str());

    ofstream outcounts;
    outcounts.open(output_count_file.c_str());

	for(int shift = min_shift; shift <= max_shift; shift ++ )
	{
	    // kmer counts in foreground
	    // note that each count vector could be of different length due to kmer length difference
	    message("counting exact kmers in foreground sequences...shift="+to_string(shift));
	    map<string, vector<int> > data1 = count_all_kmer_in_seqs(kmers, seqs1, shift);

	    // kmer counts in foreground
	    message("counting exact kmers in background sequences...shift="+to_string(shift));
	    map<string, vector<int> > data2 = count_all_kmer_in_seqs(kmers, seqs2, shift);
		
	    // compute p-value for exact kmers
	    message("computing p-values for exact kmers...");
	    for (map<string,vector<int> >::iterator it=data1.begin(); it!=data1.end(); ++it) 
	    {
	        // counts in foreground
	        vector<int> counts1 = it->second;
	        // counts in background
	        vector<int> counts2 = data2[it->first];

	        // total counts 
	        int total_counts1 = sum(counts1);

	        // for each position
	        for( int m=0;m<counts1.size();m++)
	        { 
	            if (counts1[m] == 0 && counts2[m] == 0) continue; 
	            double f2 = float(counts2[m]+pseudo)/nSeq2;
	            // compare background p estimate from shuffling and markov model
	            double p = binom_test(nSeq1,counts1[m],f2);
	            if (p < pCutoff)
	            {
					if(p == 0) p = 1e-16;
	                double corrected_p = min(1.0,p*nTest);					
					if(Bonferroni && corrected_p > pCutoff) continue;
	                nSig++;
	                double f1 = float(counts1[m])/nSeq1;
	                double expected = nSeq1 * f2;
	                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
	                // local enrichment
	                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
					
	                outstream << it->first << "\t"  << m-startPos + 2 << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	               
	                    outcounts << it->first << ":" << m-startPos +2;
	                    for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
	                    outcounts << endl;
	               

	            }
	        }
	    }

	    // compute p-value for non-exact kmers, only work without shift
	    if (degenerate)
	    {   

	        message(to_string(nSig) + " significant exact kmers identified");
	        message("computing p-values for degenerate kmers...");

	        // map from degeneate bases to all allowed bases, e.g. R => AG
	        map<char,string> define_iupac = define_IUPAC();
	        // for output purpurse only, e.g. will use g to replace H = [ACT] = non-G
	        // map<char,string> interpret_iupac = interpret_IUPAC();

	        // for each degenerate kmer, combine element kmer counts
	        for( int i=0;i<dkmers.size();i++)
	        {
	            // if already in exact kmers, i.e. contain no degenerate bases, skip
	            if ( find(kmers.begin(), kmers.end(), dkmers[i])!=kmers.end() ) continue;

	            // expand a degenerate kmer to all possible element/exact kmers
	            vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);

	            // for output only, also generate regex version of the degenerate kmer
	            string dkmersexp_readable = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // initialize the count at each position with the first element kmer
	            vector<int> counts1  = data1[dkmersexp[0]];
	            vector<int> counts2  = data2[dkmersexp[0]];

	            // sum over the rest element kmers
	            for( int j=1;j<dkmersexp.size();j++)
	            {
	                counts1 = sum(counts1,data1[dkmersexp[j]]);
	                counts2 = sum(counts2,data2[dkmersexp[j]]);
	            }

	            // for output only
	            string regex_format = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // total counts
	            int total_counts1 = sum(counts1);

	            // calculate  p-value 
	            for( int m=0;m<counts1.size();m++)
	            { 
	                if (counts1[m] == 0 && counts2[m]==0) continue;
	                double f2 = float(counts2[m]+pseudo)/nSeq2;
	                double p = binom_test(nSeq1,counts1[m],f2);
	                if (p < pCutoff)
	                {
						if(p == 0) p = 1e-16;
	                    double corrected_p = min(1.0,p*nTest); 
						if(Bonferroni && corrected_p > pCutoff) continue;
	                    nSig++;
	                    double f1 = float(counts1[m])/nSeq1;
	                    double expected = nSeq1 * f2;
	                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
	                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
						
						
	                    outstream  << dkmers[i] << "\t" << m-startPos + 2 << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r  << endl;
	                        outcounts << dkmers[i] << ":" << m-startPos + 2;
	                        for( int x=0;x<counts1.size();x++) outcounts << "\t" << double(counts1[x])/nSeq1;
	                        outcounts << endl;

	                }
	            }
	        }
	    }

	}

    outstream.close();
    outcounts.close();

    return nSig;
} // end of function



void plot_frequency_for_significant_kmer(string inputfile, string outputfile)
{
    string script = 
    "pdf('"+outputfile+"',width=10,height=5) \n"
    "con  <- file('"+inputfile+"', open = 'r') \n"
    "while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) { \n"
    "   myVector <- (strsplit(oneLine, '\\t')) \n"
    "   x=myVector[[1]] \n"
    "   name=x[1] \n"
    "   freq=as.numeric(x[2:length(x)])\n"
    "   plot(freq,type='h',xlab='position',ylab='frequency',main=name) \n"
    "} \n"
    "close(con) \n"
    "dev.off() \n";

    R_run(script);
    
} 
   



char complement(char ch){
        switch (ch)
        {
        case 'A':return 'T';
        case 'C':return 'G';
        case 'G':return 'C';
        case 'T':return 'A';
        default:return 'N';
        }
}

string reverseComplement(string seq){
        int L = seq.length();
        string rc (L,'0');
        for( int i=0;i<L; i++)
        {
                rc[L-i-1] = complement(seq[i]);
        }
        return rc;
}

// check if file is fasta by looking at the first character
bool is_fasta(string filename)
{
	ifstream fin(filename.c_str());
	string line;
	getline(fin,line);
	return line[0] == '>';
}

void ReadOneSeqFromFasta(ifstream& infile, string& name, string& seq){
  // read one sequence from fasta file
  getline(infile,name); // read identifier line
  name = name.substr(1);// remove leading '>'
  seq = "";		// initialize sequence
  string str;		
  while(infile.peek() != '>' && infile.good())
  {// before next '>' and before hitting the end of the file
        getline(infile,str);
        seq.append(str);
    }
  seq = to_upper(seq);
  // trim
  seq.erase(seq.find_last_not_of(" \n\r\t")+1);
} 

//read all sequences in a fasta file
map<string,string> ReadFasta(string filename){
  ifstream fin(filename.c_str());
  
  map<string,string> seqs;
  string name,seq;
  while(fin.good())
  {
    ReadOneSeqFromFasta(fin,name,seq);
    seqs[name] = seq;
  }
  fin.close();
  return seqs;
}

void ReadFastaToVectors(string filename, vector<string> &names, vector<string> &seqs)
{
  ifstream fin(filename.c_str());
  string name,seq;
  while(fin.good())
  {
    ReadOneSeqFromFasta(fin,name,seq);
    names.push_back(name);
	seqs.push_back(seq);
  }
  fin.close();
}

void WriteFasta(map<string,string> seqs, string filename){
    ofstream fout;
    fout.open(filename.c_str());
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
    {
        fout << ">" << it->first << endl << it->second << endl;
    }
}

vector<string> first_n_bases(vector<string> seqs,int n){
    // take the first n bases of each sequences
    // if the sequences is shorter than n, discard it
    vector<string> res;
    for (int i=0;i<seqs.size();i++)
    {   
        if(seqs[i].size() >= n)
        {
            res.push_back(seqs[i].substr(0,n));
        }
    }
    return res;
}

vector<string> last_n_bases(vector<string> seqs,int n){
    // take the last n bases of each sequences
    // if the sequences is shorter than n, discard it
    vector<string> res;
    for (int i=0;i<seqs.size();i++)
    {   
        if(seqs[i].size() >= n)
        {
            res.push_back(seqs[i].substr(seqs[i].size()-n,n));
        }
    }
    return res;
}

// start from position a and ends at b
// if a or b < 1: distance from end
// examples
// from position 10 to the end: a=10, b=0
// the last 10 letters except the last one: a=-10,b=-1
// first 10: a=1;b=10
// last 10: a=-9,b=0;
vector<string> sub_sequences(vector<string> seqs, int a, int b)
{
    vector<string> res;
    int start,end;
    for (int i=0;i<seqs.size();i++)
    {  
		int L = seqs[i].size();
		if (a>0) start = a;
        else start = L + a -1;
		if (b>0) end = b;
		else end = L + b - 1;  

		if( (start > 0 && start <= end) && (end>0 && end < L))      	
        {  
            res.push_back(seqs[i].substr(start-1,end-start+1));
        }
    }

	if(res.size()<1)
	{
		message("ERROR: sub_sequences: no sequences survive the triming!");
		exit(1);
	}

    return res;
}
// convert fasta file to a letter matrix
void fasta_to_letter_matrix(string input, string output){
	ofstream out(output.c_str());
	map<string,string> seqs = ReadFasta(input);
	
	// sequence length, assume fixed
	int L = int(seqs.begin()->second.size());
	
	//header
	out << "SeqID";
	for(int i=0;i<L;i++) out << "\t" << "k1:p" << i+1 ;
	out << endl;
	
    for (map<string,string>::iterator it=seqs.begin(); it!=seqs.end(); ++it)
	{
	    out << it->first;
		for(int i=0;i<L;i++) out << "\t" << it->second[i] ;
		out << endl;
	}
	out.close();
}

// convert fasta file to a letter matrix, no header
void tab_seq_to_letter_matrix(string input, string output, int k_min, int k_max, int col){
	col = col - 1;
		
	ifstream in(input.c_str());
	ofstream out(output.c_str());

	string line;
	vector<string> flds;
	
	while(in.good())
	{
		getline(in,line);
		if(line.length()==0) continue;
		if(line[0] == '#') continue;
		flds = string_split(line);
		// output other columns first
		for (int i=0;i<flds.size();i++) 
		{
			if(i != col) out << flds[i] << "\t";
		}
		// split sequence as letters
		for(int k = k_min;k<=k_max;k++)
			for(int i=0;i<=flds[col].size()-k;i++) 
				out << flds[col].substr(i,k) << "\t" ;
		out << endl;
	}
	in.close();
	out.close();
}

/**************** ushuffle *************/

//uShuffle
//http://digital.cs.usu.edu/~mjiang/ushuffle/
// shuffle sequence preserving k-let
string shuffle_seq_preserving_k_let(string str,int k){
    const char * c = str.c_str();
    char *t = new char[str.size() + 1];
    shuffle(c,t, str.size(), k);
    t[str.size()] = '\0';
    string res(t);
    return res;
}

//
vector<string> shuffle_seqs_preserving_k_let(vector<string> seqs, int N, int k){
    // obtain a time-based seed
    srand(time(NULL));
    vector<string> seqs2;
    for (int i=0;i<seqs.size();i++)
    {
        for( int i =0;i<N;i++) seqs2.push_back( shuffle_seq_preserving_k_let(seqs[i],k));
    }
    return seqs2;
}



/*******  seq_match   *******/


void mismatches(map<string,string>& mutant,map<string,int>& dist, string motif, int n, set<char> alphabet){
  set<char>::iterator it;
  if (mutant.count(motif) == 0)
  {
    mutant[motif] = "";
    dist[motif]=n;
  }
  if(n==0){return;}
  for( int i=0;i<motif.length();i++)
  {
      string str=motif;
      set<char> ab = alphabet;
      ab.erase(str[i]);
      for (it = ab.begin(); it!=ab.end(); it++)
      {
         str[i] = *it;
         //cout << "mutate "<<motif<<" to "<<str<<endl;
         if (mutant.count(str) >0)
         {
           if(dist[str] >= n)
           {
             //cout << mutant[str] <<endl;
             continue;
           }
         }
         
         //mutated to a new sequence
           //cout <<"new mutation"<<endl;
           mutant[str] = mutant[motif];
           mutant[str].push_back(','); 
           mutant[str].push_back(motif[i]);
           mutant[str].append(to_string(i+1));
           mutant[str].push_back(str[i]);
           dist[str]=n;
           //cout << "tag="<<mutant[str]<<" dist="<<n<<endl;
   
         if (n>1)
         {
           //cout << "subproc" <<endl;
           mismatches(mutant,dist,str,n-1,alphabet);
         }
      }

  }
}

map<string,string> ExpandMotifs(map<string,string>& motifs, int nmismatch, bool rc, set<char> alphabet) { 
  map<string,string> expandedmotifs;
  // generate mismatched motifs
  map<string,string> mutants;
  map<string,int> dist;
  map<string,string>::iterator it;   // iterator for motifs
  map<string,int>::iterator it2; // iterator for mutants
  string name,seq,tmp;
  //cout<<"input motifs"<<endl;
  //PrintMap(motifs);
  for(it=motifs.begin();it!=motifs.end();it++)
  {
    name = (*it).first;
    seq = (*it).second;
   
    mismatches(mutants,dist,seq,nmismatch,alphabet);
    //cout << mutants.size()<<" mutants identified" <<endl;
    //PrintMap(mutants);
    // add mutants to motifs
    for(it2=dist.begin();it2!=dist.end();it2++)
    {
      string tmp = name;
      tmp.append(",").append((*it2).first);
      tmp.append(mutants[(*it2).first]);
      expandedmotifs[tmp] = (*it2).first;
      //cout << name <<","<<tmp<<","<<expandedmotifs[tmp]<<endl;
    }
    // clear the mutants list
    mutants.clear();
    dist.clear();
  }
  //PrintMap(expandedmotifs);
  //cout << expandedmotifs.size() <<" expanded motifs"<<endl;
  //cout <<"add reverse complement"<<endl;
  map<string,string> expandedmotifs_rc = expandedmotifs;
  if (rc)
  {
    for(it=expandedmotifs.begin();it!=expandedmotifs.end();it++)
    {
      name = (*it).first;
      expandedmotifs_rc[name.append(",rc")] = reverseComplement((*it).second);
    }
  }   
 
  return expandedmotifs_rc;
 
}



array<int,2> match(string motiffile, string seqfile, string outfile, int nmismatch, bool rc, set<char> alphabet) {
  int nsite = 0; // total number of matches to report
  int nseq = 0;
  ifstream fmotif, fseq;
  ofstream fout;
  
  // load motifs
  map<string,string> motifs = ReadFasta(motiffile);
  cout <<"["<<current_time()<<"] "<<motifs.size()<< " motifs loaded from "<<motiffile<<endl;
  
  // expand motifs
  map<string,string> expandedmotifs = ExpandMotifs(motifs,nmismatch,rc,alphabet);
  cout <<"["<<current_time()<<"] "<<expandedmotifs.size()<< " motifs after expanding (mismatch/reverse-complement)"<<endl;
  
  //PrintMap(expandedmotifs);
  
  // searching motifs in each sequence
  fseq.open(seqfile.c_str());
  fout.open(outfile.c_str());
  
  string seqname,seq,motifname,motif;
  while(fseq.good())
  {
    // read one sequence
    ReadOneSeqFromFasta(fseq,seqname,seq);
    nseq = nseq + 1;

    cout.flush();
    // iterate over motifs
    map<string,string>::iterator it;
    for(it=expandedmotifs.begin();it!=expandedmotifs.end();it++)
    {
      motifname = (*it).first;
      motif = (*it).second;
      //cout << "searching for "<<motifname<<":"<< motif <<endl;
      set<int> found = findall(seq,motif);
      for (set<int>::iterator it=found.begin(); it!=found.end(); ++it)
      {
        fout <<seqname<<"\t"<< *it << "\t"<< motifname <<"\t"<<motif<<endl; 
      }
      nsite = nsite + found.size();
    }
    cout <<"\r["<<current_time()<<"] " << nsite << " sites found in "<< nseq << " sequences             "  ;
  }
  
  cout << endl; 
  fseq.close();
  fout.close();
  
  array<int,2> res = {{nsite, nseq}};
  
  return res;
}


int tab2bed_galaxy(string infile, string outfile){
  //hg18_chr6_122208322_122209078_+     635     5ss,A7C,G8T-rc  AAGTACCTG
  //hg18_chr6_122208322_122209078_+     553     5ss,C1G,G3A     GAAGTAAGT
  ifstream fin;
  ofstream fout;
  fin.open(infile.c_str());
  fout.open(outfile.c_str());
  string line;
  vector<string> flds;
  vector<string> pos;
  vector<string> nm;
  while(fin)
  {
    getline(fin,line);
    if (line.length() == 0)
      continue;
    flds = string_split(line);
    pos = string_split(flds[0],"_");
    if (pos.size() < 5)
    {
      cout << "\n!! incorrect sequence name format!\n make sure sequence name looks like: hg18_chr6_122208322_122209078_+\n\n";
      return 0;
    }
    if (pos.size() > 5)
    {// something like chr1_random, skip
      continue;
    }
    string chr = pos[1];
    int start = atoi(pos[2].c_str());
    int end = atoi(pos[3].c_str());
    int match_start = atoi(flds[1].c_str());
    int motifLen = flds[3].length();
    // check if match on the other strand
    string strandness = "sense";
    if (flds[2].find("rc") !=string::npos)
      strandness = "antisense";
    string strand = pos[4];
    if (strand== "+")
    {
      start = start + match_start;
      if (strandness == "antisense")
      {
        strand = "-";
      }
    }
    else//sequence on the - strand of the genome
    {
      start = end - match_start - motifLen;
      if (strandness == "antisense")
      {
        strand = "+";
      }
    }
    end = start + motifLen;
    // number of mismatches
    nm = string_split(flds[2],",");
    int score = nm.size()-2;
    if (strandness == "antisense") {score--;}

    fout << chr <<"\t"<<start<<"\t" <<end<<"\t"<<flds[3]<< "\t"<<score <<"\t"<<strand<< "\t"<<strandness<<"\t"<<flds[2]<<"\t"<<flds[0]<<"\t"<<flds[1]<<endl;

  }
  fin.close();
  fout.close();
  return 1;//return 1 if successfully created bed file
}


int tab2bed_bedtools(string infile, string outfile){
  //chrX:20597309-20645164(-)       7533    u1,CAGGTAAGT    CAGGTAAGT
  //chr17:70312707-70951085(+)      486494  u1,CAGGTAAGT,rc ACTTACCTG
  //can be without strand

  ifstream fin;
  ofstream fout;
  fin.open(infile.c_str());
  fout.open(outfile.c_str());
  string line,strand,chr;
  vector<string> flds,flds0,flds1,flds2,pos,nm;
  while(fin)
  {
    getline(fin,line);
    if (line.length() == 0)
      continue;
    flds = string_split(line);
    flds0 = string_split(flds[0],":");
    chr = flds0[0];
    strand = "+";
    // seee if there is ( 
    if (flds0[1].find("(") == string::npos)
    { 
      pos = string_split(flds0[1],"-");
    }
    else
    {
      flds1 = string_split(flds0[1],")");
      flds2 = string_split(flds1[0],"(");
      strand = flds2[1];
      pos = string_split(flds2[0],"-");
    }
    int start = atoi(pos[0].c_str());
    int end = atoi(pos[1].c_str());
    int match_start = atoi(flds[1].c_str());
    int motifLen = flds[3].length();
    // check if match on the other strand
    string strandness = "sense";
    if (flds[2].find("rc") !=string::npos)
      strandness = "antisense";
    
    if (strand== "+")
    {
      start = start + match_start;
      if (strandness == "antisense") 
      {
        strand = "-";
      }    
    }
    else//sequence on the - strand of the genome
    {
      start = end - match_start - motifLen;
      if (strandness == "antisense") 
      {
        strand = "+";
      }
    }
    end = start + motifLen;
    // number of mismatches
    nm = string_split(flds[2],",");
    int score = nm.size()-2;
    if (strandness == "antisense") {score--;} 
   
    fout << chr <<"\t"<<start<<"\t" <<end<<"\t"<<flds[3]<< "\t"<<score <<"\t"<<strand<< "\t"<<strandness<<"\t"<<flds[2]<<"\t"<<flds[0]<<"\t"<<flds[1]<<endl;
    
  }
  fin.close();
  fout.close();
  return 1;//return 1 if successfully created bed file
}


