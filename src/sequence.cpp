#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <array>        // std::array
#include <algorithm>

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

// protein sequence alignment scoring using BLOSUM

// load pwm from file, contain only the matrix
/*
A	R	N	D	C	Q	E	G	H	I	L	K	M	F	P	S	T	W	Y	V	B	Z	X	*
4	-1	-2	-2	0	-1	-1	0	-2	-1	-1	-1	-1	-2	-1	1	0	-3	-2	0	-2	-1	0	-4
-1	5	0	-2	-3	1	0	-2	0	-3	-2	2	-1	-3	-2	-1	-1	-3	-2	-3	-1	0	-1	-4
-2	0	6	1	-3	0	0	0	1	-3	-3	0	-2	-3	-2	1	0	-4	-2	-3	3	0	-1	-4
-2	-2	1	6	-3	0	2	-1	-1	-3	-4	-1	-3	-3	-1	0	-1	-4	-3	-3	4	1	-1	-4
0	-3	-3	-3	9	-3	-4	-3	-3	-1	-1	-3	-1	-2	-3	-1	-1	-2	-2	-1	-3	-3	-2	-4
-1	1	0	0	-3	5	2	-2	0	-3	-2	1	0	-3	-1	0	-1	-2	-1	-2	0	3	-1	-4
-1	0	0	2	-4	2	5	-2	0	-3	-3	1	-2	-3	-1	0	-1	-3	-2	-2	1	4	-1	-4
0	-2	0	-1	-3	-2	-2	6	-2	-4	-4	-2	-3	-3	-2	0	-2	-2	-3	-3	-1	-2	-1	-4
-2	0	1	-1	-3	0	0	-2	8	-3	-3	-1	-2	-1	-2	-1	-2	-2	2	-3	0	0	-1	-4
-1	-3	-3	-3	-1	-3	-3	-4	-3	4	2	-3	1	0	-3	-2	-1	-3	-1	3	-3	-3	-1	-4
-1	-2	-3	-4	-1	-2	-3	-4	-3	2	4	-2	2	0	-3	-2	-1	-2	-1	1	-4	-3	-1	-4
-1	2	0	-1	-3	1	1	-2	-1	-3	-2	5	-1	-3	-1	0	-1	-3	-2	-2	0	1	-1	-4
-1	-1	-2	-3	-1	0	-2	-3	-2	1	2	-1	5	0	-2	-1	-1	-1	-1	1	-3	-1	-1	-4
-2	-3	-3	-3	-2	-3	-3	-3	-1	0	0	-3	0	6	-4	-2	-2	1	3	-1	-3	-3	-1	-4
-1	-2	-2	-1	-3	-1	-1	-2	-2	-3	-3	-1	-2	-4	7	-1	-1	-4	-3	-2	-2	-1	-2	-4
1	-1	1	0	-1	0	0	0	-1	-2	-2	0	-1	-2	-1	4	1	-3	-2	-2	0	0	0	-4
0	-1	0	-1	-1	-1	-1	-2	-2	-1	-1	-1	-1	-2	-1	1	5	-2	-2	0	-1	-1	0	-4
-3	-3	-4	-4	-2	-2	-3	-2	-2	-3	-2	-3	-1	1	-4	-3	-2	11	2	-3	-4	-3	-2	-4
-2	-2	-2	-3	-2	-1	-2	-3	2	-1	-1	-2	-1	3	-3	-2	-2	2	7	-1	-3	-2	-1	-4
0	-3	-3	-3	-1	-2	-2	-3	-3	3	1	-2	1	-1	-2	-2	0	-3	-1	4	-3	-2	-1	-4
-2	-1	3	4	-3	0	1	-1	0	-3	-4	0	-3	-3	-2	0	-1	-4	-3	-3	4	1	-1	-4
-1	0	0	1	-3	3	4	-2	0	-3	-3	1	-1	-3	-1	0	-1	-3	-2	-2	1	4	-1	-4
0	-1	-1	-1	-2	-1	-1	-1	-1	-1	-1	-1	-1	-1	-2	0	0	-2	-1	-1	-1	-1	-1	-4
-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	-4	1
*/
map<char,map<char,int>> LoadScoreMat(string filename){
	
	map<char,map<char,int>> scoreMat;
		
 	ifstream fin(filename.c_str());
	string line;
	string alphabet;
	vector<string> flds;
	
	// read line 1
	getline(fin,line);
	if(line.length() == 0){
		message("ERROR! line 1 empty in scoring file: "+filename);
		system_run("touch exit_with_error");exit(1);
	}
	flds = string_split(line);
	int N=flds.size();
	for(int i=0;i<N;i++) alphabet[i]=flds[i][0];

	int i=0;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		if (N != flds.size())
		{
			message("ERROR! rows have different number of columns in scoring file: "+filename);
			system_run("touch exit_with_error");exit(1);
		}
		for(int j=0;j<N;j++) {
			scoreMat[alphabet[i]][alphabet[j]]=stoi(flds[j]);
			//cout << alphabet[i] << "\t" << alphabet[j] << "\t" << scoreMat[alphabet[i]][alphabet[j]] << endl;
		}
		i = i + 1;
		if(i > N)
		{
			message("ERROR! scoring file has too many rows: "+filename);
			system_run("touch exit_with_error");exit(1);
		}
	}
	fin.close();

	return scoreMat;
}

// score two sequences for alignment score
//
int seqAlignmentScore(map<char,map<char,int>> M, string seq1, string seq2)
{
	if(seq1.size() != seq2.size()){
		message("ERROR! the two sequences are different in size: \n"+seq1+"\n"+seq2);
		system_run("touch exit_with_error");exit(1);
	} 
	int score=0;
	for(int i=0;i<seq1.size();i++) score += M[seq1[i]][seq2[i]];
	return score;
}



// return position of the char in str that is not found in alphabet
// if no such char, return -1
int bad_char(string str, string alphabet)
{
	int i,j;
	for (  i=0;i<str.size();i++)
	{
		for (  j=0;j<alphabet.size();j++)
		{
			if (str[i] == alphabet[j]) break;
		}
		if (j == alphabet.size()) 
			{ 
				//message(str+"_"+to_string(i)+"_"+str[i]+"_");
				return i; 
			}
	}
	return -1; // if not found
}

vector<int> bad_char(vector<string> str, string alphabet)
{
	vector<int> res = {-1,-1};
	int i,j;
	for (  i=0;i<str.size();i++)
	{
		j = bad_char(str[i],alphabet);
		if (j > -1)
		{
			res = {i,j};
			return res;
		}
	}
	return res;
}


// initialize a pwm from a sequence
boost::numeric::ublas::matrix<double> initialize_pwm_from_one_seq(string seq)
{
    // nucleotide to position
    map<char,int> residual2pos;
    residual2pos['A'] = 0;
    residual2pos['C'] = 1;
    residual2pos['G'] = 2;
    residual2pos['T'] = 3;
	
    boost::numeric::ublas::matrix<double> pwm (4,seq.size());
    for (int i = 0; i < pwm.size1(); ++ i) 
		for (int j = 0; j < pwm.size2(); ++ j) 
			pwm(i,j) = 0;
    for (int i = 0; i < pwm.size2(); ++ i) pwm(residual2pos[seq[i]],i) = 1.0;
  
    return pwm;
}

// score a sequence by pwm
double score_by_pwm(boost::numeric::ublas::matrix<double> pwm, map<char,int> residual2pos, string seq)
{
	double score = 0;
	for (int i = 0; i < pwm.size2(); ++ i) score += pwm(residual2pos[seq[i]],i);
	return score;
}

// find best match by pwm in a seq
// pos: position of best match
// return best match score
double best_score_by_pwm(boost::numeric::ublas::matrix<double> pwm, map<char,int> residual2pos, string seq, int &pos)
{
	double max_score = -100000;
	double score;
	int motif_size = int(pwm.size2());
	for(int i=0;i< int(seq.size())-motif_size;i++)
	{
		//cout << i << "," << seq.size()-motif_size << "," << motif_size << endl;
		score =  score_by_pwm(pwm,residual2pos,seq.substr(i,motif_size));
		if (score > max_score)
		{
			max_score = score;
			pos = i;
		}
	}
	return max_score;
}

// update pwm from sequences
// scores: normalized to 0-1 ()
// weight from each sequence: pwm score of the match * sequence score
boost::numeric::ublas::matrix<double> update_pwm_from_seqs(vector<string> seqs, vector<double> scores, boost::numeric::ublas::matrix<double> pwm,  map<char,int> residual2pos)
{
	// initialize a zero matrix
    boost::numeric::ublas::matrix<double> pwm2 (4,pwm.size2());
    for (int i = 0; i < pwm2.size1(); ++ i) 
		for (int j = 0; j < pwm2.size2(); ++ j) 
			pwm2(i,j) = 0;
	
	// update from each sequence
	int max_pos;
	double max_score;
	double composite_score;
	for(int i=0;i<seqs.size();i++)
	{
		max_score = best_score_by_pwm(pwm,residual2pos,seqs[i],max_pos);
		composite_score = max_score * scores[i];
		for (int j = 0; j < pwm2.size2(); ++ j) pwm2(residual2pos[seqs[i][max_pos+j]],j) += composite_score; 
	}

	// normalize so that each column sum up to 1
	for(int c=0;c<pwm2.size2();c++) // for each column
	{
		double row_sum = 0;
		for(int r=0;r<pwm2.size1();r++)
			row_sum += pwm2(r,c);
		for(int r=0;r<pwm2.size1();r++)
			pwm2(r,c) = pwm2(r,c)/row_sum;
	}
	
    return pwm2;
}


/*

// smith-waterman scoring without gap
// a and b should be of the same size
// mismatch_penalty: negative number
double smith_waterman_score(string a, string b, double mismatch_penalty)
{
	double score = 0;
	double score_max = 0;
	for(int i=0;i<a.size();i++)
	{
		if(a[i] == b[i]) score++;
		else score += mismatch_penalty;
		if (score > score_max) score_max = score;
	}
	return score_max;
}

void find_best_smith_waterman_match(string motif, string seq, double mismatch_penalty, int &pos, double &score)
{
	pos = 0;
	score = -10000000;
	for(int i = 0;i<seq.size()-motif.size();i++)
	{
		double s = smith_waterman_score(motif, seq, mismatch_penalty);
		if (s > score) 
		{
			score = s;
			pos = i;
		}
	}
}

*/

/********* for regression motif analysis ********/

// for each motif, count frequency in each target sequence

vector<int> motif_counts_in_seqs(string motif, vector<string> seqs)
{
	vector<int> counts;
	for( int i=0;i<seqs.size();i++)
	{		
		counts.push_back(count_non_overlap(seqs[i],motif));
	}
	return counts;
}
/*******************/

// test if a kmer's presence is associated with sequence scores/weights
// do not require position info
array<double,4>  kmer_rank_test(string kmer, vector<string> seqs, vector<double> ranks)
{
	vector<double> ranks_with_kmer;
	array<double,2> utest;
	for(int i=0;i<seqs.size();i++)
	{
		if(seqs[i].find(kmer) != string::npos) ranks_with_kmer.push_back(ranks[i]);
	}
	if (ranks_with_kmer.size()>2 && ranks_with_kmer.size() < seqs.size()-1 ) // cannot have too few or too many 
	{
		utest = Mann_Whitney_U_test(ranks_with_kmer, seqs.size()); // z, p
	} else
	{
		utest = {{0,1}};
		
	}
	array<double,4> res = {{double(seqs.size()),double(ranks_with_kmer.size()),-utest[0],utest[0]/fabs(utest[0])*log10(utest[1])}};

	// debug
	//cout << "debug" << endl;
	//cout << seqs.size() << endl;
	//cout << seqs[0] << endl;
	//cout << kmer << endl;
	//cout << ranks_with_kmer.size() << endl;
	
	return res;
}


void kmer_cdf(string kmer, vector<string> seqs, vector<double> scores, string outputfile)
{
	ofstream fout(outputfile.c_str());
	for(int i=0;i<seqs.size();i++)
	{
		if(seqs[i].find(kmer) != string::npos) fout << "1\t" << scores[i]  << endl;
		else fout << "0\t" << scores[i] << endl;
	}
	fout.close();
	
	string script = 
    "pdf('"+outputfile+".pdf',width=5,height=5) \n"
	"par(cex=1)\n"
    "x  = read.table('"+outputfile+"', header=F) \n"
	"pos = x[x[,1]==1,2]\n"
	"neg = x[x[,1]==0,2]\n"
    "median_p = round(median(pos),digits=3)\n"
	"median_n = round(median(neg),digits=3)\n"
    "if(median_p < median_n){\n"
	"    kstest = ks.test(pos,neg,alternative='greater')\n"
    "} else {\n"
    "    kstest = ks.test(pos,neg,alternative='less')\n "
    "}\n"
	"p = kstest$p.value\n"
	"ecdf_pos = ecdf(pos)\n"
	"ecdf_neg = ecdf(neg)\n"
	"xtick = sort(c(pos[(length(pos)*0.01):(length(pos)*0.99)],neg[(length(neg)*0.01):(length(neg)*0.99)]))\n"
	"plot(xtick,ecdf_pos(xtick),type='l',col='blue', bty='n',main=paste('ks test p=',format(p,digits=2),'\nmedian dif = ',round(median_p-median_n,digits=3),'=',median_p,'-',median_n,sep=''),xlab='score', ylab='cdf') \n"	
	"lines(xtick,ecdf_neg(xtick),col='red')\n"
	"legend('topleft',c(paste(length(pos),' with "+kmer+"',sep=''),paste(length(neg),' without "+kmer+"',sep='')),col=c('blue','red'),box.lty=0,lty=1,cex=1)\n"
		
    "dev.off() \n"
	"\n";

	R_run(script);		
}


// kmer can be degenerate, such as CNNC
void positional_kmer_cdf(string kmer, vector<string> seqs, vector<double> scores, string outputfile)
{
	// get all positional kmers, including 
	vector<positional_kmer> pkmers = positional_kmer_vector_from_string(kmer,define_IUPAC());
	if(pkmers.size()<1){
		cerr << "ERROR: no positional kmer recognized: "+kmer << endl;
		return;
	}
	
	ofstream fout(outputfile.c_str());
	for(int i=0;i<seqs.size();i++)
	{
		if(seq_has_any_of_positional_kmer(seqs[i],pkmers)) fout << "1\t" << scores[i]  << endl;
		else fout << "0\t" << scores[i] << endl;
	}
	fout.close();
	
	string script = 
    "pdf('"+outputfile+".pdf',width=5,height=5) \n"
	"par(cex=1)\n"
    "x  = read.table('"+outputfile+"', header=F) \n"
	"pos = x[x[,1]==1,2]\n"
	"neg = x[x[,1]==0,2]\n"
    "median_p = round(median(pos),digits=1)\n"
	"median_n = round(median(neg),digits=1)\n"
    "if(median_p < median_n){\n"
	"    kstest = ks.test(pos,neg,alternative='greater')\n"
    "} else {\n"
    "    kstest = ks.test(pos,neg,alternative='less')\n "
    "}\n"
	"p = kstest$p.value\n"
	//"ecdf_pos = ecdf(pos)\n"
	//"ecdf_neg = ecdf(neg)\n"
	//"xtick = sort(c(pos[(length(pos)*0.05):(length(pos)*0.95)],neg[(length(neg)*0.5):(length(neg)*0.95)]))\n"
	//"plot(xtick,ecdf_pos(xtick),type='l',col='blue', bty='n',main=paste('ks test p=',format(p,digits=2),'\nmedian dif: ',median_p,'-',median_n,'=',round(median_p-median_n,digits=1),sep=''),xlab='score', ylab='cdf') \n"	
	//"lines(xtick,ecdf_neg(xtick),col='red')\n"
	"plot(ecdf(pos),xlim=c(min(quantile(pos,0.025),quantile(neg,0.025)),max(quantile(pos,0.975),quantile(neg,0.975))),do.points = FALSE,verticals=TRUE,col='blue',bty='l',main=paste('ks test p=',format(p,digits=2),'\\nmedian dif: ',median_p,'-',median_n,'=',round(median_p-median_n,digits=1),sep=''),xlab='score', ylab='CDF')\n"
	"lines(ecdf(neg),do.points = FALSE,verticals=TRUE,col='red',)\n"
	"legend('topleft',c(paste(length(pos),' with "+kmer+"',sep=''),paste(length(neg),' without "+kmer+"',sep='')),text.col=c('blue','red'),box.lty=0,cex=1)\n"	
    "dev.off() \n"
	"\n";

	R_run(script);		
}


string reverse(string str)
{
	string res = str;
	std::reverse(res.begin(), res.end());
	return res;
}

// seq can not have residuals not in valid_residuals
bool valid_sequence(string seq, string valid_residuals)
{
	for (int i=0;i<seq.size();i++)
	{
		bool valid = false;
		for(int j=0;j<valid_residuals.size();j++)
		{
			if (seq[i] == valid_residuals[j])
			{
				valid = true;
				break;
			}
		}
		if (valid == false) return false;
	}
	return true;
}

// homopolymer: find the longest runs of any residual of residuals in s
// return length and start position (0-based)
array<int,2> find_longest_run(string s, string residuals)
{
	int max_L = 0;
	int start=-1;
	int L=0;
	// starting at each position
	for(int i =0;i<s.size();i++)
	{
		if (residuals.find(s[i])!=std::string::npos) 
		{
			int j;
			for(j=i+1;j<s.size();j++)
			{
				if (residuals.find(s[j])==std::string::npos) break;
			}
			L = j - i;
			if(L > max_L)
			{
				max_L = L;
				start = i;
			}
		}
	}
	array<int,2> res={max_L,start};
	return res;
}

// DNA sequence feature, output to screen
// 
void sequence_feature(string s, bool di, bool tri, bool p1, bool p2, bool h2, bool h3, bool h4)
{

	// string header = "Seq\tAn.L\tAn.P\tCn.L\tCn.P\tGn.L\tGn.P\tTn.L\tTn.P\tRn.L\tRn.P\tYn.L\tYn.P\tWn.L\tWn.P\tSn.L\tSn.P\tA\tC\tG\tT\tW\tS\tR\tY\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT\tAAA\tCAA\tGAA\tTAA\tACA\tAGA\tATA\tCCA\tCGA\tCTA\tGCA\tGGA\tGTA\tTCA\tTGA\tTTA\tAAC\tAAG\tAAT\tCAC\tCAG\tCAT\tGAC\tGAG\tGAT\tTAC\tTAG\tTAT\tACC\tACG\tACT\tAGC\tAGG\tAGT\tATC\tATG\tATT\tCCC\tCCG\tCCT\tCGC\tCGG\tCGT\tCTC\tCTG\tCTT\tGCC\tGCG\tGCT\tGGC\tGGG\tGGT\tGTC\tGTG\tGTT\tTCC\tTCG\tTCT\tTGC\tTGG\tTGT\tTTC\tTTG\tTTT";

	// column 1: sequence
	cout << s ;

    // column 2-9: longest runs and center position
    vector<string> runs = {"A","C", "G", "T","AG","TC","AT","GC"};

    for(int i=0;i<runs.size();i++)
    {  
        array<int,2> res = find_longest_run(s,runs[i]);
        cout << "\t" << res[0] << "\t" << res[0]/2 + res[1] ;
    }


	// column 10-17: mono + A/T + G/C + A/G + T/C
	int nA = count(s,"A");
	int nC = count(s,"C");
	int nG = count(s,"G");
	int nT = count(s,"T");
	
	cout << "\t" << nA << "\t" << nC << "\t" << nG << "\t" << nT << "\t" << nA+nT << "\t" << nG+nC << "\t" << nA+nG << "\t" << nT+nC;	


	if(di)
	{	
		// column 18-33: di
		array<string,16> di = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};

		for (int i=0;i<di.size();i++) cout << "\t" << count(s,di[i]) ;
	}

	if(tri)
	{
		vector<string> tri = generate_kmers(3,"ACGT");
		for (int i=0;i<tri.size();i++) cout << "\t" << count(s,tri[i]) ;
	}
	//header = header +"\t"+to_string(tri);


	if(p1)
	{
		for(int i=0;i<s.size();i++)
		{
			if(s[i] == 'A') cout << "\t1\t0\t0\t0";
			else if (s[i] == 'C') cout << "\t0\t1\t0\t0";
            else if (s[i] == 'G') cout << "\t0\t0\t1\t0";
            else cout << "\t0\t0\t0\t1";
		}		
	}

    if(p2)
    {
		array<string,16> di = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
        for(int i=0;i< s.size()-1;i++)
        {
			for(int j=0;j<16;j++)
			{
				if(s.substr(i,2) == di[j]) cout << "\t1";
				else cout << "\t0";
			}
        }
    }


	// homopolymer
	if(h2)
	{
		array<string,4> homo = {"AA","CC","GG","TT"};
		for(int i=0;i< s.size()-1;i++)
        {   
            for(int j=0;j<4;j++)
            {  
                if(s.substr(i,2) == homo[j]) cout << "\t1";
                else cout << "\t0";
            }
        }
	}
    if(h3)
    {   
        array<string,4> homo = {"AAA","CCC","GGG","TTT"};
        for(int i=0;i< s.size()-2;i++)
        {      
            for(int j=0;j<4;j++)
            {  
                if(s.substr(i,3) == homo[j]) cout << "\t1";
                else cout << "\t0";
            }
        }
    }
    if(h4)
    {   
        array<string,4> homo = {"AAAA","CCCC","GGGG","TTTT"};
        for(int i=0;i< s.size()-3;i++)
        {      
            for(int j=0;j<4;j++)
            {  
                if(s.substr(i,4) == homo[j]) cout << "\t1";
                else cout << "\t0";
            }
        }
    }

	cout << endl;

	//cerr << header << endl;
}

// features specified via strings
// freq_feat: A,CG,...
// posi_feat: A.-3,CG.9,...
// col: column of sequence, 0 based
// start: 1 based
void sequence_feature(string infile, string outfile, string freq_feat, string posi_feat, int col_seq/*=1*/, int col_id/*=1*/,int start/*=1*/)
{
	vector<string> freq_feats = string_split(freq_feat,",");
	vector<string> posi_feats = string_split(posi_feat,",");
	vector<string> feat;
	vector<int> posi;
	for(int i=0;i<posi_feats.size();i++)
	{
		vector<string> tmp = string_split(posi_feats[i],".");
		feat.push_back(tmp[0]);
		int p = stoi(tmp[1]);
		if (p<0) p = p + 1;
		p = p + start - 2;
		posi.push_back(p);
	}
	
	ofstream fout(outfile.c_str());
	
	fout << "SeqID\t" << to_string(freq_feats) << "\t" << to_string(posi_feats) << endl;
	
	ifstream fin(infile.c_str());
	
    string line;
    vector<string> flds;

    while(fin)
    {  
        getline(fin,line);
        if (line.length() == 0)
            continue;

        line.erase(line.find_last_not_of(" \n\r\t")+1);
		
		flds = string_split(line,"\t");

		string seq = to_upper(flds[col_seq]);
		string id = flds[col_id];
		
		fout << id ;
		
		// count freq
		for(int i=0;i<freq_feats.size();i++)
		{
			fout << "\t" << count(seq,freq_feats[i]);
		}

		//presence of position motif
		for(int i=0;i<feat.size();i++)
		{
			if(seq.substr(posi[i],feat[i].size()) == feat[i]) fout <<"\t1" ;
			else fout << "\t0";
		}
		
		fout << endl;

	}
	fin.close();
	fout.close();
}

void to_upper(vector<string> &seqs)
{
	for(int i=0;i<seqs.size();i++) seqs[i] = to_upper(seqs[i]);
}

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


vector<int> filter_sequences_by_alphabet(vector<string> &seqs, string alphabet)
{
    vector<int> removed;

    for(int i=seqs.size()-1;i>=0;i--)
    {  
        if(bad_char(seqs[i],alphabet) != -1)
        {  
            seqs.erase(seqs.begin()+i);
            removed.push_back(i);
        }
    }
    return removed;
}

// any of a list positional kmer is present in sequence
bool seq_has_any_of_positional_kmer(string seq, vector<positional_kmer> pkmers)
{
	//cout << pkmers.size() << endl;

	for (int i=0;i<pkmers.size();i++)
	{
		//cout << pkmers[i].as_string() << endl;
		if (pkmers[i].is_present_in(seq))
		    return true;
	}
	return false;
}

// filter sequence by positional kmer
// seqs: original sequence and negative after filering
// positives: positive after filtering
// value: bool vector indicating positives
vector<bool> filter_sequences_by_kmer(vector<string> &seqs, vector<string> &positives, vector<positional_kmer> pkmers)
{
    vector<bool> is_positive;

    for(int i=seqs.size()-1;i>=0;i--)
    {  
        if(seq_has_any_of_positional_kmer(seqs[i],pkmers))
        {  
			positives.push_back(seqs[i]);
            seqs.erase(seqs.begin()+i);
            is_positive.push_back(true);
        }
		else is_positive.push_back(false);
    }
	// note that the order needs to be reversed
    reverse(is_positive.begin(),is_positive.end());
	return is_positive;
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
	int minCount/*=5*/,
	string alphabet/*=ACGT*/ )
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
		vector<string> exp_kmers;
		if (alphabet == "ACGT") exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		else exp_kmers = {kmers[i]};
		
		int k = kmers[i].size();
		for (int shift = min_shift; shift <= max_shift; shift ++)
		{
			for( int pos=0; pos < lSeq-k+1; pos ++) // at each position
			{
				// ranks of sequences with this kmer at this position
				vector<double> ranks;
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
//	message(to_string(lSeq)+" length");
	
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

				//cout << kmers[i] << "\t" << pos-startPos+2 << "\t" << weights1.size() << "," << weights2.size() << endl;

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
	int minCount/*=5*/,
	string alphabet/*=ACGT*/) 
{
	//message(alphabet);
	
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
        //vector<string> exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		
		vector<string> exp_kmers;
		if (alphabet == "ACGT") exp_kmers = expand_degenerate_kmer(kmers[i],define_iupac);
		else exp_kmers = {kmers[i]};
		
		
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
			
            // if this kmer is not present anywhere, skip it
            if (total_counts == 0) continue;

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
	                double z = (counts[x] - expected) / sqrt(expected*(1-f2));
					
					if(p == 0) 
					{
						p = p_value_for_z_score(z);
						if(p == 0) p = 1e-300;
					}

	                double corrected_p = min(1.0,p*nTest);	
					if (Bonferroni && corrected_p > pCutoff) continue;
					
	                nSig++;
					double f1 = float(counts[x])/nSeq;
	                double local_r = f1 / ((f2 * nSeq - f1)/(nSeq-1));// double(counts[x])/(total_counts-counts[x])*(nSeq-1);
					int position = x-startPos+2;
					if (position < 1) position -= 1;
	                outstream << kmers[i] << "\t" << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
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
	int minCount/*=5*/,
	string alphabet) 
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

// nucleotide plot from kpLogo2 weighted output
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
		"    residuals=unlist(strsplit(txt,'')) \n"
		"    thisy<-y \n"
		"    for(txtstr in 1:length(residuals)) { \n"
		"        text(x,thisy,residuals[txtstr],col=col[[residuals[txtstr]]],adj=0,srt=90,cex=cex) \n"
		"        thisy<-thisy+strheight(residuals[txtstr],cex=cex)*1.1 \n"
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

// use motif end position instead of start position
void use_end_position(string filename)
{
    ifstream fin(filename.c_str());
	string output = filename + ".tmp";
	ofstream fout(output.c_str());
    string line;
    vector<string> flds;

    while(fin)
    {
        getline(fin,line);
        if (line.length() == 0)
            continue;

        if(line[0] == '#') continue;

        flds = string_split(line);
		int oldpos = stoi(flds[1]);
		int pos  = stoi(flds[1])+int(flds[0].size())+stoi(flds[2])-1;
		if (oldpos < 0 && pos >=0) pos = pos + 1;
		fout << flds[0] << "\t" << to_string(pos);
		for(int i=2;i<flds.size();i++)
			fout << "\t" << flds[i];
		fout  << endl;
	}
	fin.close();
	fout.close();
	system_run("mv "+output+" "+filename);
}

// note input should be sorted by weight
vector<positional_kmer> build_model_from_kpLogo_output(string filename, int startPos)
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
double score_sequence_using_kpLogo_model(vector<positional_kmer> ranked_kmers,  string seq)
{
	double score = 0;
	for( int i=0;i<ranked_kmers.size();i++)
	{
		//cout << i << "\t"<< ranked_kmers[i].as_string() ;
		//cout << "\t" << seq.substr(ranked_kmers[i].pos-1,ranked_kmers[i].size) ;
		size_t found = seq.substr(ranked_kmers[i].pos,ranked_kmers[i].size).find(ranked_kmers[i].seq);
		if (found != std::string::npos) 
		{
			score += ranked_kmers[i].weight;
			//cout << "\t" << score;
		}
		//cout << endl;
	}
	//cout << "done with seq: " << seq << endl;
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
double score_sequence_using_kpLogo_model_use_group(vector<positional_kmer> ranked_kmers,  string seq)
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

// for both kpLogo and kpLogo2
void score_fasta_using_kpLogo_model(string seqfile, string outputfile, vector<positional_kmer> ranked_kmers)
{
    ifstream fin(seqfile.c_str());
	ofstream fout(outputfile.c_str());
  
    string name,seq;
    while(fin.good())
    {
		ReadOneSeqFromFasta(fin,name,seq);
     	double score = score_sequence_using_kpLogo_model(ranked_kmers, seq);
  		fout << name << "\t" << seq << "\t" << score << endl;		
    }
    fin.close();
	fout.close();
}

// sequence in column col, 1 based
void score_tabular_using_kpLogo_model(string tabfile, int col, string outputfile, vector<positional_kmer> ranked_kmers)
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
     	double score = score_sequence_using_kpLogo_model(ranked_kmers, flds[col]);
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

// read kpLogo output into two vectors
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

// read kpLogo2 output into two vectors
void read_significant_positional_kmer_from_kpLogo2_output(string inputfile, vector<string> &kmers, vector<int> &positions)
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
void significant_feature_matrix_kpLogo2(vector<string> seqs, vector<double> weights, vector<positional_kmer> ranked_kmers, string outfile)
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

// kpLogo: remove overlapping motifs
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

int count(string seq,string motif)
{
    // count the number of matches, allowing overlap 
    int n = 0;
    size_t pos = 0; //start at pos
    pos = seq.find(motif); // the first match
    while(pos != string::npos)
    {  
        n++;
        pos = seq.find(motif,pos+1); // start search for the next match from pos + 1
    }
    return n;
}

int count_non_overlap(string seq,string motif)
{
    // count the number of matches, not allowing overlap 
	int motif_size = motif.size();
    int n = 0;
    size_t pos = 0; //start at pos
    pos = seq.find(motif); // the first match
    while(pos != string::npos)
    {  
        n++;
        pos = seq.find(motif,pos+motif_size); // start search for the next match from pos + motif_size
    }
    return n;
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

// generate all possible combinations of residuals in alphabet of fixed length k, i.e. kmer
vector<string> generate_kmers(int k, string alphabet)
{
    vector<string> kmers;
    
    // first position, each residual in alphabet
    for(int i = 0; i< alphabet.size(); i++)
    {
        string s(1,alphabet[i]);
        kmers.push_back(s);
    }

    // for the rest k-1 positions, fill all possible residuals 
    for(int i = 0; i< k-1;i++)
    {
        int n = kmers.size();
        // for each current kmer, append all possible residuals
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


// for max_shift > 0: take the most strongest
boost::numeric::ublas::matrix<double> position_weight_matrix_from_kpLogo_output(string filename, string alphabet, int seqLen, int startPos, int cScore)
{
	int cSeq = 0;
	int cPos = 1; 
	int cStat = 3;
 	cScore -= 1;
	
	// initialize an empty matrix
    boost::numeric::ublas::matrix<double> pwm (alphabet.size(),seqLen);
    for (int i = 0; i < pwm.size1(); ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<string,int> residual2pos;
	for(int i=0;i<alphabet.size();i++)
	{
		residual2pos[alphabet.substr(i,1)] = i;
	}
 
 	ifstream fin(filename.c_str());
	string line;
	vector<string> flds;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		if(residual2pos.find(flds[cSeq]) != residual2pos.end())
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
			if(abs(pwm(residual2pos[flds[cSeq]],pos)) < abs(score)) // for multiple shift values
			{
				pwm(residual2pos[flds[cSeq]],pos) = score;
			}
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

// to avoid calculating information content, set nSeq=0
// if nSeq > 0, calculate information content and small sample correction
// if nSeq == 0, do not calculate information content
// if nSeq < 0, calculate info content without small sample correction
// small sample correction: subtract the following term from each value: (alphabet.size-1)/ln2/2/n
// descending: false. if true, will put important residual at the bottom
void generate_ps_logo_from_pwm(boost::numeric::ublas::matrix<double> pwm, string filename, string alphabet,vector<int> fixed_position, vector<string> fixed_residual, map<char,string> colors, double score_cutoff, int startPos, int fontsize, string y_label, double max_scale, int nSeq, int bottom_up)
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

	double small_sample_correction = 0;
	
	if (nSeq != 0) 
	{
		if(nSeq > 0)
		{
			small_sample_correction = (alphabet.size()-1.0)/log(2.0)/2.0/nSeq;
		}

		vector<double> info_content = matrix_column_information_content(pwm);

		for (int i=0;i<pwm.size2();i++)
		{
			info_content[i] -= small_sample_correction;
			if(info_content[i] < 0) info_content[i] = 0;
			for (int j=0;j<pwm.size1();j++)
			{
				pwm(j,i) *= info_content[i];
			}
		}
	}

	// debug
	//print_matrix(pwm);


	double maxv = matrix_max(pwm);
	double minv = matrix_min(pwm);
	vector<double> col_max = matrix_column_max(pwm);
	vector<double> col_min = matrix_column_min(pwm);
		
	// for fixed position
	double max_col_sum = max(matrix_column_sum(pwm,1));
	for (int i=0;i<fixed_position.size();i++)
	{
		for(int j=0;j<alphabet.size();j++)
		{
			if(fixed_residual[i] == alphabet.substr(j,1))
			{
				pwm(j,fixed_position[i]) = 1.1 * max_col_sum;
			}
			else pwm(j,fixed_position[i]) = 0 ;
		}
	}


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
		"%%Creator: kpLogo 					\n"
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
		// fixed position
		if (find (fixed_position.begin(), fixed_position.end(), i) != fixed_position.end() )  color_sig = "0 0 0";
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
		vector<size_t> idx = sort_indexes(w,bottom_up);
		//for(int j=0;j<4;j++) cout << i << "\t" << residuals[idx[j]] << "\t"<< idx[j] << "\t" << w[idx[j]] << endl;
		
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

void remove_kmers_overlapping_with_fixed_positions(string infile,vector<int> fixed_positions,int startPos)
{	
	string outputfile = infile+".remove_kmers_overlapping_with_fixed_positions";
	ofstream fout(outputfile.c_str());
	
 	ifstream fin(infile.c_str());
	string line;
	vector<string> flds;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		
		int kmer_start = stoi(flds[1]);
		if (kmer_start > 0) kmer_start = kmer_start + startPos - 2;
		else kmer_start = kmer_start + startPos - 1;
		int kmer_end = kmer_start + flds[0].size() - 1;
		bool no_overlap = true;
		for (int i=0;i<fixed_positions.size();i++)
		{
			if (kmer_start <= fixed_positions[i] && kmer_end >= fixed_positions[i])
			{
				no_overlap = false;
				break;
			}
		}
		if (no_overlap) fout << line << endl;
	}
	fin.close();
	fout.close();
	system_run("mv "+outputfile+" "+infile);
}
	

void postscript_logo_from_kpLogo_output(string infile, string outfile, vector<int> fixed_position, vector<string> fixed_residual, map<char,string> colors, int seqLen, double score_cutoff, int startPos, int fontsize, int cScore,string y_label, double max_scale)
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
		"%%Creator: kpLogo 					\n"
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
		if (find (fixed_position.begin(), fixed_position.end(), i) != fixed_position.end() )  color_sig = "0 0 0";
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
	
	// fixed positions
	for(int i=0;i<fixed_position.size();i++)
	{
		string s = string(fixed_residual[i]);
		//cout << s << endl;
		out << postscript_kmer(s, x0+ xstep * (fixed_position[i]), y0+coord_size, fontsize, scaley, 1, 1.1*maxv/ absmax * max_scale, colors, 0);
	}
	
	out << "showpage" << endl;
	out << "%%EOF" << endl;
	
	out.close();
			
}

// reverse a pwm columns
boost::numeric::ublas::matrix<double> reverse_pwm(boost::numeric::ublas::matrix<double> pwm)
{
	boost::numeric::ublas::matrix<double> pwm2 (pwm.size1(),pwm.size2());
	for (int i = 0; i < pwm.size1(); ++ i)
	    for (int j = 0; j < pwm.size2(); ++ j)
	        pwm2(i,j) = pwm(i,pwm.size2()-j-1);
	return pwm2;
}

// load pwm from file, contain only the matrix
boost::numeric::ublas::matrix<double> load_pwm_from_file(string filename, string alphabet){
	// determine the size of the matrix
 	ifstream fin(filename.c_str());
	string line;
	vector<string> flds;
	int i=0;
	int L=0;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		if (L==0) L=flds.size();
		else if (L != flds.size())
		{
			message("ERROR! rows have different number of columns in pwm file: "+filename);
			system_run("touch exit_with_error");exit(1);
		}
		i = i + 1;
				
		if(i > alphabet.size())
		{
			message("ERROR! pwm file has too many rows: "+filename);
			system_run("touch exit_with_error");exit(1);
		}
	}
	fin.close();
	if(i < alphabet.size())
	{
		message("ERROR! pwm file has too few rows: "+filename);
		system_run("touch exit_with_error");exit(1);
	}
	
	// initialize an empty matrix
	boost::numeric::ublas::matrix<double> pwm (alphabet.size(),L);
	for (int i = 0; i < alphabet.size(); ++ i)
	    for (int j = 0; j < pwm.size2(); ++ j)
	        pwm(i,j) = 0.0;
	
    // load from file
	fin.open(filename.c_str());
	i=0;
	while(fin)
	{
		getline(fin,line);
		if(line.length() == 0) continue;
		flds = string_split(line);
		for (int j=0;j<L;j++)
		{
			pwm(i,j) = stof(flds[j]);
		}
		i = i + 1;
	}
	fin.close();

	return pwm;
}



// write pwm in meme format
void write_pwm_in_meme_format(boost::numeric::ublas::matrix<double> pwm, string motifname, string filename)
{
    ofstream out;
    out.open(filename.c_str());
    out << "MEME version 4.4\n\n";
    out << "ALPHABET= ACGT\n\n";
    out << "strands: + -\n\n";
    out << "Background residual frequencies\n";
    out << "A 0.25 C 0.25 G 0.25 T 0.25" << "\n\n";
    out << "MOTIF "<<motifname<<"\n\n";
    out << "residual-probability matrix: alength= 4 w= "<<pwm.size2()<<" nsites= 1000 E= 0"<<endl;
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
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_seqs(vector<string> seqs,string alphabet)
{

    // length of sequences
    int L = seqs[0].size();

    // initialize an empty matrix
    boost::numeric::ublas::matrix<double> pwm (alphabet.size(),L);
    for (int i = 0; i < alphabet.size(); ++ i)
        for (int j = 0; j < pwm.size2(); ++ j)
            pwm(i,j) = 0.0;
  
    //cout << "matrix initiaze ok" << endl;

    // nucleotide to position
    map<char,int> residual2pos;
	for(int i=0;i<alphabet.size();i++) residual2pos[alphabet[i]] = i;
 
    for(int k=0;k<seqs.size();k++)
    {
        for( int i =0;i<L;i++)
        {
        	pwm(residual2pos[seqs[k][i]],i) += 1;
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
    map<char,int> residual2pos;
    residual2pos['A'] = 0;
    residual2pos['C'] = 1;
    residual2pos['G'] = 2;
    residual2pos['T'] = 3;
 
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
                    pwm(residual2pos[seqs[k][j]],j-start) += 1;
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

//kpLogo : create logo for a single kmer
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

//kpLogo
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
	string output_count_file,
	int minCount,
	string alphabet)
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
			if(f2 == 1) f2 = 0.9999;


	        // for each position
	        for( int m=0;m<counts1.size();m++)
	        { 
	            if (counts1[m] == 0 && f2 ==0) continue; 
				if (counts1[m] < minCount || counts1[m] > nSeq1 - minCount) continue;
	            // compare background p estimate from shuffling and markov model
	            double p = binom_test(nSeq1,counts1[m],f2);
	            if (p < pCutoff)
	            {
	                double f1 = float(counts1[m])/nSeq1;
	                double expected = nSeq1 * f2;
	                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
					if(p == 0) 
					{
						p = p_value_for_z_score(z);
						if(p == 0) p = 1e-300;
					}
	                double corrected_p = min(1.0,p*nTest);					
					if(Bonferroni && corrected_p > pCutoff) continue;
	                nSig++;

	                // local enrichment
	                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
					int position = m-startPos + 2;
					if(position < 1) position --;
	                outstream << it->first << "\t" << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	                outcounts << it->first << ":" << position;
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
	            // string dkmersexp_readable = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // initialize the count at each position with the first element kmer
	            vector<int> counts1  = data1[dkmersexp[0]];
	            // sum over the rest element kmers
	            for( int j=1;j<dkmersexp.size();j++) counts1 = sum(counts1,data1[dkmersexp[j]]);

	            // calculate probs from markov model: sum of element prob
	            double f2 = probs_kmer[dkmersexp[0]];
	            for( int j=1;j<dkmersexp.size();j++) f2 += probs_kmer[dkmersexp[j]]; 
	            //f2 = f2 * (1+shift); // can't have shift on degenerate kmers

	            // for output only
	            //string regex_format = degenerate_kmer_to_regex(dkmers[i],define_iupac);

	            // total counts
	            int total_counts1 = sum(counts1);

	            // calculate  p-value 
	            for( int m=0;m<counts1.size();m++)
	            { 
	                if (counts1[m] == 0 && f2==0) continue;
	                double p = binom_test(nSeq1,counts1[m],f2);
	                if (p < pCutoff)
	                {
	                    double f1 = float(counts1[m])/nSeq1;
	                    double expected = nSeq1 * f2;
	                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
						if(p == 0) 
						{
							p = p_value_for_z_score(z);
							if(p == 0) p = 1e-300;
						}
	                    double corrected_p = min(1.0,p*nTest);
						if(Bonferroni && corrected_p > pCutoff) continue;
						
	                    nSig++;

	                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
						
						int position = m-startPos+2;
						if(position < 1) position--;
												
	                    outstream  << dkmers[i] << "\t" << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	        
                       outcounts << dkmers[i] << ":" << position;
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
int find_significant_kmer_from_two_seq_sets(vector<string>seqs1, vector<string>seqs2, vector<string>kmers, vector<string> dkmers, int min_shift,int max_shift,bool degenerate,double pCutoff, bool Bonferroni/*=true*/,double pseudo,int startPos,int nTest, string outfile,string output_count_file,int minCount,string alphabet)
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
				if (counts1[m] < minCount || counts1[m] > nSeq1 - minCount) continue;
				
	            double f2 = float(counts2[m]+pseudo)/nSeq2;
				if (f2 >= 1) f2 = 1-pseudo/nSeq2;
				
	            // compare background p estimate from shuffling and markov model
	            double p = binom_test(nSeq1,counts1[m],f2);
	            if (p < pCutoff)
	            {
	                double f1 = float(counts1[m])/nSeq1;
	                double expected = nSeq1 * f2;
	                double z = (counts1[m] - expected) / sqrt(expected*(1-f2));
					if(p == 0) 
					{
						p = p_value_for_z_score(z);
						if(p == 0) p = 1e-300;
					}
	                double corrected_p = min(1.0,p*nTest);					
					if(Bonferroni && corrected_p > pCutoff) continue;
	                nSig++;

	                // local enrichment
	                double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
					
					int position = m-startPos+2;
					if (position < 1) position -= 1;
					
	                outstream << it->first << "\t"  << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r << endl;
	               
	                    outcounts << it->first << ":" << position;
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
	            //vector<string> dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);
				
				vector<string> dkmersexp;
				if (alphabet == "ACGT") dkmersexp = expand_degenerate_kmer(dkmers[i],define_iupac);
				else dkmersexp = {kmers[i]};

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
					if (counts1[m] < minCount || counts1[m] > nSeq1 - minCount) continue;

					double f2 = float(counts2[m]+pseudo)/nSeq2;
					if (f2 >= 1) f2 = 1-pseudo/nSeq2;

	                double p = binom_test(nSeq1,counts1[m],f2);
	                if (p < pCutoff)
	                {
	                    double f1 = float(counts1[m])/nSeq1;
	                    double expected = nSeq1 * f2;
	                    double z = (counts1[m] - expected) / sqrt(expected*(1-f2)); 
						if(p == 0) 
						{
							p = p_value_for_z_score(z);
							if(p == 0) p = 1e-300;
						}
	                    double corrected_p = min(1.0,p*nTest); 
						if(Bonferroni && corrected_p > pCutoff) continue;
	                    nSig++;

	                    double local_r = double(counts1[m]) / (total_counts1 - counts1[m]) * (counts1.size()-1);
						
		                int position = m-startPos+2;
			            if (position < 1) position -= 1;

	                    outstream  << dkmers[i] << "\t" << position << "\t" << shift << "\t" << z << "\t" << -log10(p) << "\t" << -log10(corrected_p) << "\t" << f1 << "\t" << f2 << "\t" << f1/f2 << "\t"  << local_r  << endl;
	                        outcounts << dkmers[i] << ":" << position;
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
	    str.erase(str.find_last_not_of(" \n\r\t")+1);
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
// 1-based
// if a or b < 1: distance from end
// examples
// from position 10 to the end: a=10, b=0
// the last 10 residuals except the last one: a=-10,b=-1
// first 10: a=1;b=10
// last 10: a=-9,b=0;
vector<string> sub_sequences(vector<string> seqs, int a, int b)
{
    vector<string> res;
    int start,end;
    for (int i=0;i<seqs.size();i++)
    {  
		int L = seqs[i].size();
		if (a>0) start = a - 1;
        else start = L + a -1;
		if (b>0) end = b - 1;
		else end = L + b - 1;  

		if( (start >= 0 && start <= end) && (end>=0 && end < L))      	
        {  
            res.push_back(seqs[i].substr(start,end-start+1));
        }
    }

	if(res.size()<1)
	{
		message("ERROR: sub_sequences: no sequences survive the triming!");
		system_run("touch exit_with_error");exit(1);
	}

    return res;
}
// convert fasta file to a residual matrix
void fasta_to_residual_matrix(string input, string output){
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

// convert fasta file to a residual matrix, no header
void tab_seq_to_residual_matrix(string input, string output, int k_min, int k_max, int col){
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
		// split sequence as residuals
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
        for( int j =0;j<N;j++) seqs2.push_back( shuffle_seq_preserving_k_let(seqs[i],k));
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
  WriteFasta(expandedmotifs,motiffile+".expanded.fa"); 
 
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


