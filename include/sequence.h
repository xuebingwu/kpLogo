#ifndef __SEQUENCE_H__
#define __SEQUENCE_H__

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <array>        // std::array

//#include <boost/regex.hpp>

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

#include <stdio.h>
#include <stdlib.h>
#include <algorithm>

#include "positional_kmer.h"

extern "C"{
#include "ushuffle.h"
}

using namespace std;


void use_end_position(string filename);


string postscript_line(double x1, double y1, double x2, double y2);

string postscript_text(string text, double x, double y, double width, double height, string color, double rotate=0);

void generate_ps_logo_from_pwm(boost::numeric::ublas::matrix<double> pwm, string filename,string alphabet, map<char,string> colors, double score_cutoff,int startPos=1, int fontsize=20, string ylabel="-log10(p)", double max_scale=6.0, int nSeq=0, bool bottom_up=false);

boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_seqs(vector<string> seqs, string alphabet);

boost::numeric::ublas::matrix<double> position_weight_matrix_from_PKA_output(string filename, string alphabet, int seqLen, int startPos=1,int cScore=5);

void postscript_logo_from_PKA_output(string infile, string outfile, map<char,string> colors,int seqLen, double score_cutoff, int startPos=1, int fontsize=20, int cScore=5,  string y_label="-log10(p)", double max_scale=6.0);

string postscript_kmer(string text, double x, double y, int fontsize, double scaley, double width, double height, map<char,string> colormap, double rotate);

bool is_fasta(string filename);
void ReadFastaToVectors(string filename, vector<string> &names, vector<string> &seqs);


vector<bool> filter_sequences_by_kmer(vector<string> &seqs, vector<string> &positives, vector<positional_kmer> pkmers);
bool seq_has_any_of_positional_kmer(string seq, vector<positional_kmer> pkmers);

vector<positional_kmer> build_model_from_PKA2_output(string filename, double pCutoff);
vector<positional_kmer> build_model_from_PKA_output(string filename, int startPos);
void save_model_to_file(vector<positional_kmer> ranked_kmers, string filename);

vector<positional_kmer> load_model_from_file(string filename);

double score_sequence_using_PKA_model(vector<positional_kmer> ranked_kmers, string seq);



void score_fasta_using_PKA_model(string seqfile, string outputfile, vector<positional_kmer> ranked_kmers);
void score_tabular_using_PKA_model(string tabfile, int col, string outputfile, vector<positional_kmer> ranked_kmers);


class paired_kmer
{
public:
	string seq1;
	string seq2;
	int dist; // distance between seq1 start and seq2 start
	int pos; // position of seq1
	int shift;
	double weight;
	paired_kmer();
	paired_kmer(const paired_kmer &a);
	paired_kmer(string seq1, string seq2, int dist=0, int shift=0, int pos=0, double weight=0);
	const paired_kmer &operator=(const paired_kmer &a);
	string as_string(string del="_", bool add_shift=true,bool add_pos=false, bool add_weight=false);
	bool equals(paired_kmer a, bool identical_pos=true, bool identical_shift=true);
	int len(); // total length = seq1 + seq2 + gap
	int gap();
};

vector<paired_kmer> build_paired_kmer_model(string filename);

double score_sequence_using_paired_kmer_model(vector<paired_kmer> model,  string seq);

vector<paired_kmer> generate_paired_kmers (
	string alphabet,
int seq1_len,
int seq2_len,
int max_dist,
int min_dist=1,
int max_shift=0,
int min_shift=0	
);

int find_significant_pairs_from_weighted_sequences(
	vector<string> seqs,
	vector<double> weights, 
	vector<paired_kmer> paired_kmers, 
	string outfile, 
	int nTest, 
	double pCutoff=0.05, 
	bool Bonferroni=true,
	int startPos=0,
	int minCount=5) ;

	
void plot_nucleotide_profile(string infile, string outfile, int lSeq, int col, int startPos);

void plot_most_significant_kmers(string infile, string outfile, int lSeq, int column, int startPos);

void load_weighted_sequences_to_vectors(string filename, vector<string> &seqs, vector<double> &weights,int cSeq=1,int cWeight=2);

vector<string> load_ranked_sequences_to_vectors(string filename, int cSeq=1);

int find_significant_kmer_from_ranked_sequences(vector<string> seqs, vector<string> kmers, string outfile, int nTest, double pCutoff=0.05,bool Bonferroni=true,int min_shift=0, int max_shift=0, int startPos=0, int minCount=5 );

int find_significant_kmer_from_weighted_sequences(vector<string> seqs,vector<double> weights, vector<string> kmers, string outfile, int nTest, double pCutoff=0.05, 	bool Bonferroni=true,int shift_min=0, int shift_max=2, int startPos=0, int minCount=5);

void save_feature_matrix(map<string,string> seqs, vector<string> kmers, string outfile, string label, bool append=false, int shift=0);

void read_significant_positional_kmer_from_file(string inputfile, vector<string> &kmers, vector<int> &positions);

map<string, string> seq_vector2map(vector<string> seqs);

void read_significant_positional_kmer_from_PKA2_output(string inputfile, vector<string> &kmers, vector<int> &positions);

void significant_feature_matrix(map<string,string> seqs, vector<string> kmers, vector<int> positions, string outfile, string label, bool append=false, int shift=0);

void significant_feature_matrix_PKA2(vector<string> seqs, vector<double> weights, vector<positional_kmer> ranked_kmers, string outfile);

vector<int> filter_sequences_by_size(vector<string> &seqs, int lSeq=0);

// convert fasta file to a letter matrix
void fasta_to_letter_matrix(string input, string output);

// convert fasta file to a letter matrix, no header
void tab_seq_to_letter_matrix(string input, string output, int k_min, int k_max, int col);

// PKA: remove overlapping motifs
// input: all significant motifs
int non_overlapping_sig_motifs(string inputfile, string outputfile);

// count substring in a map of sequences, i.e. from fasta
// used in generating markov model
int countSubstringInSeqs(vector<string>seqs, string sub);

// of not of identical length, return -1
int sequence_similarity(string a, string b);

// calculate and write pairwise similarity matrix of
void pairwise_sequence_similarity_matrix(vector<string> seqs, string filename);

	
set<int> findall(string seq, string motif);

// dict for IUPAC degenerate nucleotides
map<char,string> define_IUPAC();


// generate all possible combinations of letters in alphabet of fixed length k, i.e. kmer
vector<string> generate_kmers(int k, string alphabet);

// degenerate kmers based on IUPAC, DNA only
// remove those with terminal N, which is not a kmer, but (k-i)mer
vector<string> degenerate_kmer(int k, string alphabet="ACGTRYMKWSBDHVN");


// convert one degenerate kmer to regular expression, e.g. ARG -> A[AG]G
string degenerate_kmer_to_regex(string kmer,map<char,string> iupac);

// expand a degenerate kmer to all possible exact kmers
vector<string> expand_degenerate_kmer(string seq, map<char,string> iupac);

// implant a motif to a set of sequences
// motif can contain degenerate nucleotides
void implant_motif(map<string,string> &seqs, int position, string motif, double fraction);


// write pwm in meme format
void write_pwm_in_meme_format(boost::numeric::ublas::matrix<double> pwm, string motifname, string filename);


// build position weight matrix for a positional kmer, can include flanking sequence
boost::numeric::ublas::matrix<double> create_position_weight_matrix_from_kmer(vector<string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos);

//PKA : create logo for a single kmer
void create_logo_for_kmer(vector<string> seqs, string kmer, int position, map<char,string> iupac, int d, int startPos,string output);

// create logos for top 1 kmer passing pCutoff_B
void create_logo_for_topN_sig_kmer_per_position(vector<string> seqs, string filename,int d, int startPos,string output);


//vector<int> count_one_kmer_in_all_seqs_regex(map<string,string> seqs, boost::regex pattern, int k, int seq_len_minus_k_plus_1, int k_plus_shift);



// doesn't work with shift+degenerate
// for each kmer, count its freq at each position(start), allowing shifts (to the right)
map<string,vector<int> > count_all_kmer_in_seqs(vector<string> kmers, vector<string> seqs, int shift);  

//PKA
void print_kmer_positional_profile(map<string,vector<int> > data);


map<string,vector<int> > degenerate_kmer_counts(vector<string> dkmers,map<string,vector<int> > data, map<char,string> define_iupac);

int find_significant_kmer_from_one_seq_set(vector<string>seqs1, map<string,double> probs_kmer,vector<string>kmers, vector<string> dkmers, int min_shift,int max_shift, bool degenerate,double pCutoff, 	bool Bonferroni, int startPos,int nTest, string outfile, string output_count_file);

// two file comparison, not allow shift and degenerate at the same time
int find_significant_kmer_from_two_seq_sets(vector<string>seqs1, vector<string>seqs2, vector<string>kmers, vector<string> dkmers, int min_shift, int max_shift, bool degenerate,double pCutoff, 	bool Bonferroni,double pseudo,int startPos,int nTest, string outfile,string output_count_file);

int find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
	vector<string> seqs,
	vector<string> kmers, 
	string outfile, 
	int nTest, 
	double pCutoff=0.05, 
	bool Bonferroni=true,
	int min_shift=0, 
	int max_shift=2,
	int startPos=0,
	int minCount=5);

void plot_frequency_for_significant_kmer(string inputfile, string outputfile);
   

char complement(char ch);

string reverseComplement(string seq);

void ReadOneSeqFromFasta(ifstream& infile, string& name, string& seq);

map<string,string> ReadFasta(string filename);

void WriteFasta(map<string,string> seqs, string filename);

map<string,string> vector2map(vector<string> seqs);

void load_sequences_from_tabular(string filename, vector<string> &seqs, vector<double> &weights, int cSeq=1, int cWeight=2);

//uShuffle
//http://digital.cs.usu.edu/~mjiang/ushuffle/
// shuffle sequence preserving k-let
string shuffle_seq_preserving_k_let(string str,int k);

//
vector<string> shuffle_seqs_preserving_k_let(vector<string> seqs, int N, int k);

vector<string> first_n_bases(vector<string> seqs,int n);

vector<string> last_n_bases(vector<string> seqs,int n);

vector<string> sub_sequences(vector<string> seqs, int a, int b);

void mismatches(map<string,string>& mutant,map<string,int>& dist, string motif, int n, set<char> alphabet);

map<string,string> ExpandMotifs(map<string,string>& motifs, int nmismatch, bool rc, set<char> alphabet) ;


array<int,2> match(string motiffile, string seqfile, string outfile, int nmismatch, bool rc, set<char> alphabet);



int tab2bed_galaxy(string infile, string outfile);


int tab2bed_bedtools(string infile, string outfile);

#endif
