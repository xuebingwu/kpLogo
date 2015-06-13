#include "sequence.h"
#include "markov.h"
#include "container.h"
#include "utility.h"
#include "text.h"


void print_help()
{
    string txt = "\n"
    "PKA: Positional Kmer Analysis (version 0.3)\n"
    "\n"
    "   -by Xuebing Wu (wuxbl@wi.mit.edu), Bartel lab, Whitehead Institute\n"
    "\n"
    "   Identify statistically enriched/depleted short sequences of length k (kmer)\n"
	"   at every position in a set of aligned sequences, weighted or unweighted. \n"
	"   Degenerate nucleotides and small shift in positions can be allowed.\n"
    "\n"
    "Usage: PKA inpu_file [options]\n"
    "\n"
    "Options\n"
    "\n"
	"Type of analysis\n"
	"   (default)            activated by default. Identify significant kmers using Binomial test\n"
	"                        - input can be fasta/raw/tabular format\n"
	"   -ranked              Wilcoxon rank-sum test (i.e. Mann-Whitney U test) on ranked sequences \n"
	"                        - input can be fasta/raw/tabular format, but needs to be sorted\n"
	"   -weighted            Two-sample stutent's t test on weighted sequences\n"
	"                        - input can only be tabular format, does not need to be sorted\n"
	"   -predict prefix      use significant kmers from a previous run (-o prefix) to score input sequences\n"
    "Input\n"
    "   -alphabet STR        alphabet for generating kmers, default=ACGT, case insensitive\n"
    "                        note: 'dna' is equivalent to 'ACGT', 'U' will be converted to 'T'\n"
	"                              'protein' is equivalent to 'ACDEFGHIJKLMNOPQRSTUVWY'\n"
	"   -seq INT             for tabular input: sequences are in column INT. Default 1\n"
	"   -weight INT          for tabular input: weights are in column INT. Default 2\n"
    "   -region n1,n2        only consider subsequences from position n1 to n2 (start at 1). \n"
    "                        non-positive numbers interpreted as distance from the end\n"
    "   -select a,b,c,...    keep sequences contain any of specified kmers a, b, c, etc\n"
    "                        kmer format: seq:position:shift, e.g. CNNC:47:0 \n"
    "   -remove a,b,c,...    remove sequences contain any of specified kmers.\n"  
    "   -fix maxFreq         fix a position with a specific letter if it occurs in more than maxFreq of the sequences\n"  
	"                        fixed letters will be plotted as 1.1x of hight of the position with the highest total height\n"
    "Kmer counting\n"
    "   -k INT               use fixed kmer length INT\n"
	"   -max_k INT           consider all kmers of length 1,2,...,INT. default=4 \n"
    "   -shift INT           max shift (to right) allowed for kmer positions\n"
	"   -max_shift INT       consider shift from 0 to INT, default=0, i.e. no shift\n"	
	"   -degenerate STR      alphabet to use for degenerate kmers. Subset of all possible IUPAC DNA \n"
	"                        letters (ACGTRYMKWSBDHVN, equivalent to 'all'). One can use \n"
	"                        ACGTN to search gapped-kmers. Only work for DNA/RNA sequences\n"
	"   -gapped              allowing gapped kmer, equivalent to '-degenerate ACGTN' \n"
	"   -pair                also test all possible pairs of positional monomers\n"  			 	
    "Statistics & output\n"
    "   -o STR               prefix for all output files, default=PKA\n"
	"   -minCount NUM        minimum number of sequences to have this kmer to include in output\n"
	"                        if smaller than 1, treat as fraction of sequences (default=5)\n"
    "   -p FLOAT             p-value cut-off, default=1.01 (i.e. output all possible kmers)\n"
    "   -pc FLOAT            Bonferroni corrected p-value cut-off, default=0.05\n"
    "   -FDR                 adjust p value by FDR method ( default is Bonferroni correction)\n"
    "   -startPos INT        re-number position INT (1,2,3,..) as 1. The position before it will be -1\n"
    "   -last_letter         use a kmer's last letter position as the kmer's position. Default is first letter\n"
    "   -pseudo FLOAT        pseudocount added to background counts. default=1e-9. Ignored by -markov\n"
	"   -fontsize INT        font size for plotting sequence logos, default 20\n"
	"   -colorblind          use colorblind friendly color scheme\n"
    "   -email STR           send email notification to this address\n"
    "   -subject STR         email subject (quoted, default='PKA job done')\n"
    "   -content STR         email content (quoted, default='PKA job done')\n"
	"   -small_sample        correct for small sample size\n"
    "   -plot STR            which statistics to plot: p: raw p-value (default), b: Bonferroni corrected p, f: FDR, s: statisitcs\n"
    "   -bottom_up           stack letters from bottom to top, starting with the most significant one\n"
    "   -save                save shuffled sequences (*.shuffled.input) or the learned markov model (*.markov.model)\n"
	"                        in -predict mode, save feature matrix (*.feat.mat)\n"
    "Background model for unweighted & unranked sequences (ignore if using -ranked or -weighted)\n"
	"   (default)            compare to the same kmer at other positions \n"
	"   -bgfile FILE         background sequence file\n"	
    "   -markov INT          N-th order markov model trained from input or background (with -bgfile)\n"
    "                        N=0,1,or 2. Default N=1: first order captures upto di-nucleotide bias\n"
    "   -shuffle N,M         shuffle input N times, preserving M-nucleotide frequency\n"
    "   -no_bg_trim          no background sequence trimming (-region). valid with -markov and -bgfile\n"
    "\n\n";
	cout << txt;
}

 
int main(int argc, char* argv[]) {

	/*
	string commandline = argv[0];
	for (int i = 1; i < argc; i++) {
		commandline += " " + string(argv[i]);
	}
	message("The command line you have entered is:");
	message(commandline);
	*/

    ///////////////////////////////////////////////////////////////
    //       part 1: parameters and default settings             //
	///////////////////////////////////////////////////////////////

	// type of analysis
	string analysis = "default"; // or ranked or weighted
	
	// prediction mode: -predict prefix
	string prefix;
	
	// default output prefix
    string output = "PKA"; 
	
    string alphabet = "ACGT";   // default DNA
	
	// kmer length, shift, degenerate bases
    int k=0;    // k
	int min_k = 1;
    int max_k = 0;  // upto max_k
	
    int shift=0;
	int min_shift = 0;
	int max_shift = 0;
	
    bool degenerate = false;    // no degenerate bases allowed by default
	
	string degenerate_alphabet = "ACGTRYMKWSBDHVN"; 
	// if degenerate allowed, default to use all IUPAC DNA bases. to look at gappmer only, use ACGTN
	
	/*
	string protein = "ACDEFGHIJKLMNOPQRSTUVWY";
	string dna = "ACGT";
	string DNA_gap = "ACGTN";
	string DNA_all = "ACGTRYMKWSBDHVN";
	*/
	
	// tabular input format
	int cSeq=0; // -seq, sequence in this column, first column is 0
	int cWeight = -1; // -weight, -1 means no weight
    int first=1; //    no 3' trim
    int last=0; //     no 5' trim
    string select_pkmers;  // select sequences with some pkmers 
    string remove_pkmers;  // remove sequences with some pkmers
    bool no_bg_trim = false;    // only valid when -b and -markov used
	
	// statistics
    double pseudo = 1e-9;   // pseudocounts, only used when -b or -shift used
    double pCutoff=1.01;    // binomial p cutoff
	double pCutoff_corrected = 0.05;
	bool Bonferroni = true;
	int small_sample_correction = -1; // -1 or 1
	
	// output
    int startPos= 1; // coordinates
	int fontsize=40;
	string plot = "p"; // or b or f or s
    bool last_letter = false;
	bool bottom_up = false;

	// color blind
	bool colorblind = false;
	map<char,string> colors;
    colors['A'] = "0.05 0.75 0.05";
    colors['C'] = "0 0 1";
    colors['G'] = "1 0.75 0";
    colors['T'] = "1 0 0";
    colors['N'] = "0.5 0.5 0.5";
    colors['Y'] = "1 1 0";
    colors['R'] = "0 1 1";
    colors['D'] = "0 0.5 0.5";
    colors['E'] = "0 0 0.5";
    colors['F'] = "0 0.5 0";
    colors['H'] = "0.5 0 0";
    colors['I'] = "1 0.5 0.5";
    colors['J'] = "0.5 1 0.5";
    colors['K'] = "0.5 0.5 1";
    colors['L'] = "1 0 0.5";
    colors['M'] = "0 1 0.5";
    colors['O'] = "0 0.5 1";
    colors['P'] = "1 0.5 0";
    colors['Q'] = "0.5 1 0";
    colors['S'] = "0.5 0 1";   
	colors['U'] = "1 0.5 1";
    colors['V'] = "0.5 0.75 1";
    colors['W'] = "1 0 1"; 
    colors['X'] = "0.5 0 0.5";
    colors['Z'] = "0.5 0.5 0";
    colors['B'] = "0.5 0.5 0.75";


    map<char,string> colorblind_colors = colors;
    colorblind_colors['A'] = "0 0.620 0.451";
    colorblind_colors['C'] = "0 0.450 0.698";
    colorblind_colors['G'] = "0.941 0.894 0.259";
    colorblind_colors['T'] = "0.835 0.369 0";
    colorblind_colors['N'] = "0.5 0.5 0.5";
	
		
	// background 
	bool local = true;
    int shuffle_N = 0;  // shuffle N times of each input sequence
    int preserve = 2;   // preserving dinucleotide
    int markov_order = -1;  // order of markov model, -1 means not used
	string mmmodel_str;
	
	double minCount = 5.0; // minimum number of sequences to have this kmer to be reported in output
	double maxFrac = 0.75; // if a letter is present in more than maxFrac * 100% of the foreground, fix it
	// for fixed positions, ignore it during p value calculation
	// for visualization, assign the p-value to be the most significant one * 1.1, color will be black?
	vector<int> fixed_position;
	vector<string> fixed_letter;
	
	bool pair = false; // test all pairs of monomers
	
	bool build_model = false;

    // email
    string email = "";
    string subject = "PKA job done";
    string content = "PKA job done";

	bool save_to_file = false;	
    string seqfile1,seqfile2,str;

	///////////////////////////////////////////////////////////////
    /*******        part 2: get commandline arguments   */
	///////////////////////////////////////////////////////////////

	// total number of arguments
    if (argc < 2) 
	{
		print_help(); // if only type PKA, print help and exit
		exit(1);
	}

	// the first argument needs to be input sequence file
    seqfile1 = argv[1]; 
    if (seqfile1 == "-h" || seqfile1 == "--help" ) 
	{
		print_help(); // PKA -h or PKA --help: print help and exit
		exit(0);
	}
	else if (seqfile1[0] == '-')
	{
		print_help();
		message("ERROR: the first argument needs to be input file name! ");
		exit(1);
	}

    // other arguments
    for (int i = 2; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-bgfile") {  // background
                seqfile2 = argv[i + 1];
				local = false;
                i=i+1;
            } else if (str == "-ranked") {
                analysis = "ranked";
				cWeight = -1;
            } else if (str == "-weighted") {
                analysis = "weighted";
				if(cWeight < 0) cWeight = 1;
			} else if (str == "-predict") {
				prefix = argv[i + 1];
				i=i+1;
            } else if (str == "-o") {   // output prefix
                output = argv[i + 1];
                i=i+1;
            } else if (str == "-k") {   // fixed k
                k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_k") {    // max k
                max_k = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-shift") {
                shift = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-max_shift") {
                max_shift = atoi(argv[i + 1]);
                i=i+1;
            } else if (str == "-alphabet") {
                alphabet = argv[i + 1];
				if(alphabet == "ACGU") alphabet="ACGT";
                i=i+1;
            } else if (str == "-degenerate") {
                degenerate = true;
				degenerate_alphabet = argv[i+1];
				if (degenerate_alphabet == "all") degenerate_alphabet = "ACGTRYMKWSBDHVN";
				i=i+1;
            } else if (str == "-gapped") {
                degenerate = true;
				degenerate_alphabet = "ACGTN";
			} else if (str == "-weight") {
				cWeight = atoi(argv[i + 1]) - 1;
				i=i+1;
			} else if (str == "-seq") {
				cSeq = atoi(argv[i + 1]) - 1;
				i=i+1;
            } else if (str == "-no_bg_trim") {
                no_bg_trim = true;
            } else if (str == "-select") { 
                select_pkmers = argv[i + 1];
                i=i+1;
            } else if (str == "-remove") {   
                remove_pkmers = argv[i + 1];
                i=i+1; 
            } else if (str == "-pseudo") {
                pseudo = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-p") {
                pCutoff = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-pc") {
                pCutoff_corrected = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-FDR") {
                Bonferroni = false;
            } else if (str == "-last_letter") {
                last_letter = true;
            } else if (str == "-startPos") {
                startPos = atoi(argv[i + 1]);
                i=i+1;
			} else if (str == "-fontsize") {
				fontsize = atoi(argv[i + 1]);
				i=i+1;
            } else if (str == "-colorblind") {
                colorblind = true;
            } else if (str == "-bottom_up") {
                bottom_up = true;
            } else if (str == "-small_sample") {
				small_sample_correction = 1;
            } else if (str == "-markov") {
				str = argv[i+1];
				if(str.find(',') == std::string::npos)
				{
					markov_order = atoi(argv[i + 1]);
				} else {
					mmmodel_str = str;
				}
				local = false;
                i=i+1;
            } else if (str == "-shuffle") {
                string s(argv[i+1]);
                vector<string> ss = string_split(s,",");
                shuffle_N = atoi(ss[0].c_str());
                preserve = atoi(ss[1].c_str());
				local = false;
                i=i+1;
            } else if (str == "-region") {
                string s(argv[i+1]);
                vector<string> ss = string_split(s,",");
                first = atoi(ss[0].c_str());
                last = atoi(ss[1].c_str());
                i=i+1;
            } else if (str == "-minCount") {
                minCount = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-fix") {
                maxFrac = atof(argv[i + 1]);
                i=i+1;
            } else if (str == "-save") {
                save_to_file = true;
            } else if (str == "-plot") {
                plot = argv[i + 1];
                i=i+1;
            } else if (str == "-email") {
                email = argv[i + 1];
                i=i+1;            
			} else if (str == "-subject") {
                subject = argv[i + 1];
                i=i+1;
            } else if (str == "-content") {
                content = argv[i + 1];
                i=i+1;
            } else {
                message("ERROR: Unknown options: "+str);
                print_help();
				exit(1);
            }
        }
    }
	

	if(colorblind) colors = colorblind_colors;
	
	if(boost::algorithm::to_lower_copy(alphabet) == "dna") 
		alphabet="ACGT";
	else if (boost::algorithm::to_lower_copy(alphabet) == "protein") 
		alphabet = "ACDEFGHIJKLMNOPQRSTUVWY";

    // determine the length of k-mer
    if (k == 0 && max_k == 0) // neither is specified, do upto 4-mers
    {
        max_k = 4;
    }
    else if ( k > 0 && max_k == 0) // 
    {  
        min_k = k;
        max_k = k; // start from k to k
    } // else do 1..max_k
	
	// determine shift
	if(max_shift == 0 && shift > 0) // -shift but not -max_shift
	{
		min_shift = shift;
		max_shift = shift;
	}
	
    if(analysis != "weighted")
		{
			cWeight = -1;
		}
	
	// both shift and degenerate for unranked unweighted
	if(analysis == "default" && degenerate == true && max_shift > 0 && local == false) 
	{
		message("Warnining: for unweighted sequences, only under -local mode one can use both shift and degenerate bases. Switch to -local. Ignore -markov / -bgfile / -shuffle");
		local = true;
	}

    // determine background model
    if(analysis == "default" && seqfile2.size()==0 && local == false){ // no file specified using -b
        if (markov_order > -1 || mmmodel_str.size()>0) // markov model
        {
            shuffle_N = 0; // ignore -shuffle
            if (markov_order>2)
            {
                message("ERROR: markov_order can only be 0, 1, or 2! ");
                exit(0);
            }
        } else if (shuffle_N == 0) // no background specified, use -markov 1
        {
            markov_order = 1;
        }
    }


    // print out parameters used
        message("Summary: ");
        message("   Input       :   " + seqfile1);
    if(first != 1 || last != 0)
	{
    	message("   Subsequence :   " + to_string(first) + ","+to_string(last));
    } 
    if(local)
	{
		message("   Background  :   local");
	}
	else if(seqfile2.size()>0)
	{
    	message("   Background  :   " + seqfile2);
    	if(markov_order > -1) message( "                   " + \
			to_string( markov_order) +  " order markov model");
    } else 
	{
    	if(markov_order > -1) message("   Markov model:   order " + to_string( markov_order));
        else
	    {
    	message( "   Background  :   shuffle input");
    	message( "                   preserving " + to_string(preserve) + "-nuceotide frequency");
        }
    }
        message("   Output      :   " +output +".*");  
        message("   alphabet    :   " + alphabet);
        message("   min_kmer    :   " + to_string(min_k) );
        message("   max_kmer    :   " + to_string(max_k) );
        message("   min_shift   :   " + to_string(min_shift) );
        message("   max_shift   :   " + to_string(max_shift) );
		if (degenerate)
        message("   degenerate  :   " + degenerate_alphabet );
        message("   p           :   " + to_string(pCutoff ));
        message("   start at    :   " + to_string(startPos));

///////////////////////////////////////////////////////////////
    /***********    part 3: process input                 */
///////////////////////////////////////////////////////////////
		
	vector<string> names, seqs1,seqs2;
	vector<double> weights;
	if (is_fasta(seqfile1)) 
	{
		if(analysis == "weighted")
		{
			message("ERROR: You used -weighted option which requires input to be tabular instead of fasta");
			exit(1);
		}
		message("Input file is FASTA format");
		ReadFastaToVectors(seqfile1, names, seqs1);
	}
	else 
	{
		load_sequences_from_tabular(seqfile1,seqs1,weights,cSeq,cWeight);
	}

    if (seqs1.size() == 0)
    {
        message("No sequence present in file: " + seqfile1);
        exit(1);
    }
    // show the number of sequences loaded
    message(to_string(seqs1.size()) + " sequences loaded from " + seqfile1);

    // trim if -first or -last specified, note that too short sequences will be discarded
    if (first != 1 || last != 0) 
	{
		seqs1 = sub_sequences(seqs1,first,last);
		message(to_string(seqs1.size()) + 
		" sequences remain after trimming");
	}

    //debug
    //WriteFasta(vector2map(seqs1),"trimmed.fa");

    // replace U with T
    if (alphabet == "ACGT")
    {
        for (int i=0;i<seqs1.size();i++)
        {
            replace(seqs1[i].begin(),seqs1[i].end(),'U','T');
        }
    }

    markov_model markov;        // background markov model
    map<string,double> kmer_probs;  // kmer probability from markov model

	if(analysis == "default" && local == false)
	{
	    // if background not available, use markov model or shuffle the foreground
	    if (seqfile2.size() == 0) // no background file specified using -b
	    {
	        if (markov_order < 0 && mmmodel_str.size()==0) // shuffle
	        {
	            // shuffle sequence preserving m-let frequency, m specified by -preserve
	            seqs2 = shuffle_seqs_preserving_k_let(seqs1,shuffle_N,preserve);

	            message(to_string(seqs2.size()) + " shuffled sequences generated, " + to_string( shuffle_N ) + " from each input sequence");

	            if (save_to_file) // save shuffled sequences
	            {        
	                WriteFasta(vector2map(seqs2),output+".shuffled.input");
	                message("Shuffled sequences saved to: "+output+".shuffled.input");
	            }
	        } else if(markov_order > -1)// generate markov model from input
	        {
	            markov = markov_model(markov_order,alphabet,seqs1);
	        } else { // generate markov model from string
				markov = markov_model(alphabet,mmmodel_str);
			}
	    }
	    // else load background sequence
	    else 
	    {
			if(is_fasta(seqfile2)) ReadFastaToVectors(seqfile2, names, seqs2);
			else load_sequences_from_tabular(seqfile2,seqs2,weights,cSeq,cWeight);
	
		    // replace U with T
		    if (alphabet == "ACGT")
		    {
		        for (int i=0;i<seqs2.size();i++)
		        {
		            replace(seqs2[i].begin(),seqs2[i].end(),'U','T');
		        }
		    }
			
	        if (seqs2.size() == 0)
	        {
	            message("No sequence present in file: " + seqfile2);
	            exit(1);
	        }
	        // show the number of sequences loaded
	        message(to_string(seqs2.size()) + " sequences loaded from " + seqfile2);
        
	        if (no_bg_trim == false || markov_order < 0)
	        {
	            // trim if -first or -last specified, note that too short sequences will be discarded
	            if (first != 1 || last != 0) 
				{
					seqs2 = sub_sequences(seqs2,first,last);
	            	message( to_string(seqs2.size()) \
					+ " background sequences remain after trimming");
				}
	        }
    
	        if (markov_order > -1)
	        {
	            markov = markov_model(markov_order,alphabet,seqs2);
	        }
	    }

	    if( (markov_order > -1 || mmmodel_str.size()>0)  && save_to_file) 
	    {   
	        markov.print(output+".markov.model");
	        message("Markov model saved to: "+output+".markov.model");
	    }
	}

    // total number of sequences loaded or after trimming
    // note that too short sequences after trimming will be discarded
    int nSeq1 = seqs1.size();
    int nSeq2;
    int seq_len1 = seqs1[0].size();
    int seq_len2;

    // make sure all sequences have the same size
    // sequences in foreground
	vector<int> removed = filter_sequences_by_size(seqs1);
	nSeq1 = seqs1.size();
	if (removed.size() > 0) 
	{
		message(to_string(nSeq1)+" sequences left after filtering by size");
		if(analysis == "weighted")
		{
			for(int i=0;i<removed.size();i++)
			{
				weights.erase(weights.begin()+removed[i]);
			}
			// save the data
			//ofstream tmpout("tmp.txt");
			//for(int i=0;i<seqs1.size();i++)
			//	tmpout << seqs1[i] << "\t" << weights[i] << endl;
			//tmpout.close();
		}
	}
        
	if(seqs1.size()>1) 
	{
		if(seqs2.size()>1) 
		{
			filter_sequences_by_size(seqs2,seqs1[0].size());
			if(seqs2.size()<2)
			{
				message("ERROR: less than 2 background sequences left after filtering by size");
				exit(1);
			}
		}
	}
    else
	{
		message("ERROR: less than 2 input sequences left after filtering by size");
		exit(1);
	}

  // filter sequences
	if (select_pkmers.size()>0)
	{
		message("selecting sequences with positional kmers: "+select_pkmers);
		vector<positional_kmer> select_pkmers_vector = positional_kmer_vector_from_string(select_pkmers,define_IUPAC());
		vector<string> positives;
    	vector<bool> is_positive = filter_sequences_by_kmer(seqs1, positives, select_pkmers_vector);
		if(seqs1.size()>0) // some sequene was removed
		{
			// remove weight
			if(analysis == "weighted")
			{
				for(int i=is_positive.size()-1;i>=0;i--)
				{
					if(is_positive[i] == false) weights.erase(weights.begin()+i);
				}
			}
		}
        seqs1 = positives;
        message(to_string(seqs1.size())+" sequences left after filtering by kmer");
	}
    if (remove_pkmers.size()>0)
    {
        message("removing sequences with positional kmers: "+remove_pkmers);
        vector<positional_kmer> remove_pkmers_vector = positional_kmer_vector_from_string(remove_pkmers,define_IUPAC());
        vector<string> positives;
        vector<bool> is_positive = filter_sequences_by_kmer(seqs1, positives, remove_pkmers_vector);
        if(positives.size()>0) // some sequene was removed
        {   
            // remove weight
            if(analysis == "weighted")
            {   
                for(int i=is_positive.size()-1;i>=0;i--)
                {  
                    if(is_positive[i] == true) weights.erase(weights.begin()+i);
                }
            }
        }
        message(to_string(seqs1.size())+" sequences left after filtering by kmer");
    }



	message("making frequency logo...");
    boost::numeric::ublas::matrix<double> pwm2 = create_position_weight_matrix_from_seqs(seqs1,alphabet);
    generate_ps_logo_from_pwm(pwm2, output+".freq.eps",alphabet,fixed_position, fixed_letter,colors,1,startPos,fontsize,"Frequency",sqrt(seq_len1)/3.0,0,bottom_up);
    system_run("ps2pdf -dEPSCrop "+output+".freq.eps "+output+".freq.pdf");
    system_run("convert "+output+".freq.eps "+output+".freq.png");

    message("making information content logo...");
    generate_ps_logo_from_pwm(pwm2, output+".info.eps",alphabet,fixed_position, fixed_letter,colors,1,startPos,fontsize,"Bits",sqrt(seq_len1)/3.0,small_sample_correction * seqs1.size(),bottom_up);
    system_run("ps2pdf -dEPSCrop "+output+".info.eps "+output+".info.pdf");
    system_run("convert "+output+".info.eps "+output+".info.png");
	
	// extract fixed position and letter
	for (int i=0;i<pwm2.size1();i++)
	{
		for (int j=0;j<pwm2.size2();j++)
		{
			if ( pwm2(i,j) > maxFrac )
			{
				fixed_position.push_back(j);
				fixed_letter.push_back(alphabet.substr(i,1));
			}
		}
	}
	if(fixed_position.size()>0)
	{
		message("fixed position: "+to_string(fixed_position));
		message("fixed letter  : "+to_string(fixed_letter));
	}

///////////////////////////////////////////////////////////////
//         part 4:  prediction mode
///////////////////////////////////////////////////////////////
	
	if(prefix.size() > 0) 
	{
		message("=== prediction mode ===");
		/**/
		message("building model from file: "+prefix+".pass.p.cutoff.txt");		
		vector<positional_kmer> ranked_kmers = build_model_from_PKA_output(prefix+".pass.p.cutoff.txt",startPos);
		message("- " + to_string(ranked_kmers.size())+" positional kmers included in the model");
	
		// save the model to file
		// save_model_to_file(ranked_kmers, output+".model.txt");	
	
		/**/
		message("scoring input sequences using the model...");
		string scoreFile = prefix+".score";
	    ofstream out1(scoreFile.c_str());
		//load_weighted_sequences_to_vectors(inputfile,seqs,weights,cSeq,cWeight);
		for (int i=0;i<seqs1.size();i++)
		{
			//cout << i << seqs1[i] << endl;
			double score = score_sequence_using_PKA_model(ranked_kmers, seqs1[i]);
			if(weights.size()>0)
				out1 << seqs1[i] << "\t" << weights[i] << "\t" << score << endl;	
			else
                out1 << seqs1[i] << "\t"  << score << endl;	
		}
		out1.close();
		message("- done");
		
		// writing feature matrix
		significant_feature_matrix_PKA2(seqs1, weights, ranked_kmers, output+".feat.mat");
	
		if (pair == false) return 0;
		
		message("loading paired monomers from file: "+prefix+".significant.pair");
		vector<paired_kmer> paired_kmer_model = build_paired_kmer_model(prefix+".significant.pair");
		message("- " + to_string(paired_kmer_model.size())+" pairs included in the model");

		message("scoring input sequences using the model...");
		string pairedScoreFile = prefix+".pair.score";
	    ofstream out(pairedScoreFile.c_str());
		for (int i=0;i<seqs1.size();i++)
		{
			double score = score_sequence_using_paired_kmer_model(paired_kmer_model, seqs1[i]);
			out << seqs1[i] << "\t" << weights[i] << "\t" << score << endl;		
		}	
		out.close();
		message("- done");
		
		
		return 0;
	}
	


/*
// debug
implant_motif(seqs1, 10, "G", 0.1);
implant_motif(seqs1, 20, "G", 0.5);
implant_motif(seqs1, 30, "CNC", 0.1);
implant_motif(seqs1, 40, "CNC", 0.5);
implant_motif(seqs1, 50, "GHG", 0.1);
WriteFasta(seqs1,"implanted.fa");
*/

    //debug
    //PrintMap(seqs1);
	

///////////////////////////////////////////////////////////////
    /********   part 5: generate kmers      */
///////////////////////////////////////////////////////////////
	
	// total number of tests to be performed, ~ n_kmer * seq_len * shift
	// to be used in multiple testing correction
	int nTest = 0; 
	
	// generate all exact kmers	
	vector<string> kmers = generate_kmers(min_k, alphabet);
	if (degenerate == false) nTest = kmers.size() * (seq_len1 - min_k + 1);
	for (k = min_k+1; k <= max_k; k++)
	{
		vector<string> tmp = generate_kmers(k, alphabet);	
		kmers += tmp;
		if (degenerate == false) nTest += tmp.size() * (seq_len1 - k + 1);
	}
	message(to_string(kmers.size()) +  " exact kmers to be tested" );
		
	vector<string> dkmers;
	// generate degenerate kmers, which will also include exact kmers
	if (degenerate)
	{
		dkmers = degenerate_kmer(min_k,degenerate_alphabet);
		nTest = dkmers.size() * (seq_len1 - min_k + 1);
		for (k = min_k+1; k <= max_k; k++)
		{
			vector<string> tmp = degenerate_kmer(k,degenerate_alphabet);	
			dkmers += tmp;
			nTest += tmp.size() * (seq_len1 - k + 1);
		}		
		message(to_string(dkmers.size()) + " k-mers allowing degenerate bases");
	}
	
	// multiple by shift
	nTest = nTest * (max_shift - min_shift + 1);
	message(to_string(nTest) + " tests (kmer x position x shifts) will be performed");


///////////////////////////////////////////////////////////////
//			part 6: search for significant kmers	
///////////////////////////////////////////////////////////////
	
	if(minCount < 1) minCount = seqs1.size() * minCount;
	
    // number of significant kmers
    int nSig = 0; // based on uncorrected and corrected p-value

	// final output file
	string out = output+".pass.p.cutoff.txt";
	
    // tmp output file for significant kmers, unsorted
    string outtmp = out+".tmp";

    // tmp output file for significnt kmers' frequency, for B significant kmers
    string output_freq = output+".Bonferroni.significant.frequency.txt";

    string header = "#kmer\tposition\tshift\tstatistics\tp-value\tcorrected.p";

	if(analysis != "default")
	{
		// output file
		
		if (analysis == "ranked") 
		{
			if(degenerate) nSig = find_significant_kmer_from_ranked_sequences(
				seqs1, dkmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
			else nSig = find_significant_kmer_from_ranked_sequences(
				seqs1, kmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
		}
		else // weighted sequecnes 
		{
			header += "\tn1\tmean1\tstd1\tn2\tmean2\tstd2";
			if(degenerate) nSig = find_significant_kmer_from_weighted_sequences(
				seqs1, weights, dkmers, outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift, startPos,minCount);
			else nSig = find_significant_kmer_from_weighted_sequences(
				seqs1, weights, kmers, outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift, startPos,minCount);
		}
	
		message( to_string (nSig) +  " positional kmers with p <"+to_string(pCutoff));
	
		if (pair)
		{
			// pairs
			int seq1_len = 1;
			int seq2_len = 1;
			int dist_min = 1;
			int dist_max = seq_len1 - seq1_len - seq2_len;
			vector<paired_kmer> paired_kmers = generate_paired_kmers ( 
				alphabet, seq1_len, seq2_len, dist_max,dist_min, max_shift,min_shift);
			message(to_string(paired_kmers.size())+" paired_kmers to test in total");
			// total number of tests
			int nTest1 = 0;
			for(int i=dist_min;i<= dist_max; i++)
			{
				nTest1 += pow(alphabet.size(),seq1_len) * pow(alphabet.size(),seq2_len)*(seq_len1 - \
					i - seq1_len - seq2_len +1);
			}
			message(to_string(nTest1)+" tests (kmer x position) to perform");
			int nSig2 = find_significant_pairs_from_weighted_sequences(
				seqs1,weights, paired_kmers, output+".significant.pair", 100, pCutoff,Bonferroni, startPos,minCount);
	
			message(to_string(nSig2)+" significant paired_kmers identified");
		}
	}
	else
	{
        header += "\tfrac.obs\tfrac.exp\tobs/exp\tlocal.r";
		if(local)
		{
			if(degenerate) nSig = find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
				seqs1, dkmers,outtmp, nTest, pCutoff, Bonferroni, min_shift, max_shift,startPos,minCount);
			else nSig = find_significant_degenerate_shift_kmer_from_one_set_unweighted_sequences(
				seqs1, kmers,outtmp, nTest, pCutoff,Bonferroni, min_shift, max_shift,startPos,minCount);
		}
	    else if (seqs2.size()>0){
	        // find significant kmers by comparing two set of sequences
	        nSig = find_significant_kmer_from_two_seq_sets(
				seqs1,seqs2,kmers,dkmers,min_shift,max_shift,degenerate,
			pCutoff,Bonferroni,pseudo,startPos,nTest,outtmp,output_freq);
	    }else{
	        // find significant kmers using markov model as background
			kmer_probs = markov.probs(kmers);
	        nSig = find_significant_kmer_from_one_seq_set(
				seqs1,kmer_probs,kmers,dkmers,min_shift,max_shift,degenerate,
			pCutoff,Bonferroni, startPos,nTest,outtmp,output_freq); 
	    }

	    message( to_string (nSig) +  " significant positional kmers identified in total"); 

	}
	
	if (nSig == 0) return 0;
	
///////////////////////////////////////////////////////////////
//			part 7: post-processing	
///////////////////////////////////////////////////////////////

	// use the last letter's position as the motif position
	if(last_letter) use_end_position(outtmp);
	
	if(Bonferroni == false || plot == "f")
	{
		////   calculating FDR
    	message("adjusting p-values using FDR...");
    	string script = 
        	"x=read.table('"+ outtmp +"',header=F)   \n"
			"adjusted = p.adjust(10^(-x[,5]),method='fdr',n="+ to_string(nTest) +")\n"
			"sub = adjusted < "+to_string(pCutoff)+"\n"
        	"x = cbind(x[sub,],-log10(adjusted[sub])) \n"
        	"write.table(x,file='"+outtmp+"',sep='\\t',col.names=F,row.names=F,quote=FALSE) ";
    	R_run(script);
		int n = count_lines(outtmp);
		message(to_string(n)+" significant kmers with FDR < "+to_string(pCutoff));
	}
	
	// sort output by p.value
	system_run("sort -k5,5gr "+outtmp+" > "+out);

	// most significant at each position
    //system_run(" cat "+out +" | sort -k2,2g -k5,5gr  > "+outtmp);
    //remove_duplicates(outtmp,output+".most.significant.each.position.txt",2,1,"");	
	
	/**/
    system_run(" cat "+out +" | awk '$4>0' | sort -k2,2g -k5,5gr -k4,4gr  > "+outtmp);
    remove_duplicates(outtmp,outtmp+".most.enriched",2,1,"");	
    system_run(" cat "+out +" | awk '$4<0' | sort -k2,2g -k5,5gr -k4,4g  > "+outtmp);
    remove_duplicates(outtmp,outtmp+".most.depleted",2,1,"");	
	system_run(" cat "+outtmp+".most.enriched "+outtmp+".most.depleted > "+output+".most.significant.each.position.txt");
	/**/
	
    //remove intermediate files
    system_run("rm "+outtmp+"*");
	
    message("plotting the most significant kmer at each position...");

	// column to plot, 1 based
	int cScore = 6; // bonferroni
	string ylabel = "-log10(p)";
	double score_cutoff = -log10(pCutoff_corrected);

	if (plot == "f") // fdr
	{
	    if (analysis == "ranked") cScore = 7;
		else if (analysis == "weighted") cScore = 13;
		else cScore = 11;
	}
	else if (plot == "p") 
	{
		cScore = 5;
		score_cutoff = -log10(pCutoff_corrected/nTest);
	}
	else if (plot == "s")
	{
		cScore = 4;
		ylabel = "test statistic";
		score_cutoff = 1e300;
	}
    //string plotfilename = output+".most.significant.each.position.pdf";

    //plot_most_significant_kmers(output+".most.significant.each.position.txt", output+".most.significant.each.position.pdf", seq_len1, cScore,startPos);
	
	postscript_logo_from_PKA_output(output+".most.significant.each.position.txt", output+".most.significant.each.position.eps",colors, seq_len1, score_cutoff, startPos, fontsize,cScore,ylabel,sqrt(seq_len1)/1.5);	

	// if monomer is included in the analysis
	if(min_k < 2) 
	{
		message("plotting single nucleotide profile...");
		//plot_nucleotide_profile( out,  output+".nucleotide.profile.pdf",  seq_len1, cScore, startPos);
		
		boost::numeric::ublas::matrix<double> pwm = position_weight_matrix_from_PKA_output(out,alphabet, seq_len1, startPos, cScore);

		//print_matrix(pwm);

	    generate_ps_logo_from_pwm(pwm, output+".eps",alphabet,fixed_position, fixed_letter,colors,score_cutoff,startPos,fontsize,ylabel,sqrt(seq_len1)/1.5,0,bottom_up);	
	}

		
	system_run("ps2pdf -dEPSCrop "+output+".eps "+output+".pdf");
	system_run("ps2pdf -dEPSCrop "+output+".most.significant.each.position.eps "+output+".most.significant.each.position.pdf");
	
	// merge pdf
	system_run("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile="+output+".all.pdf "+output+"*.pdf");
	
	// ps to png
	system_run("convert "+output+".eps "+output+".png");
	system_run("convert "+output+".most.significant.each.position.eps "+output+".most.significant.each.position.png");
	
	// add header to data file
	insert_header(out,header);
	insert_header(output+".most.significant.each.position.txt",header);
    message("Done!");

	// send email once done
	if(email.size()>1) system_run("echo '"+content+"' | mailx -s '"+subject+"' "+email);

    return 0;
} 

