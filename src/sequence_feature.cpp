#include "utility.h"
#include "iostream"
#include "fstream"

#include "sequence.h"


using namespace std;

void help()
{
    string str =
	"\n"
    "Generate sequence feature matrix\n"
	"	- Xuebing Wu (wuxb07@gmail.com)\n"
    "\n"
    "Usage: sequence_feature -i input -c n > output \n"
    "\n"
    "Options:\n"
    "\n"
    "   -i  input   input file\n"
	"   -cID     	id column number\n"
	"   -cSeq     	sequence column number\n"	
    "   -di         include di-nucletode freq\n"
    "   -tri        include tri-nucletode freq\n"
    "   -p1         include positional mono-nucleotode\n"
    "   -p2         include positional di-nucleotide\n"
	"   -h2         include 2 nt homopolymer\n"
    "   -h3         include 3 nt homopolymer\n"
    "   -h4         include 4 nt homopolymer\n"
	"   -start      start position, 1-based\n"
    "   -freq str   frequency feature      \n"
    "   -posi str   position feature      \n"
    "\n";

    cerr << str;

    exit(0);
}



int main(int argc, char* argv[]) {

    // default
    string input="input";
	int col=0; // column of seq, 0 based, but as commandline option, 1 based
	int col_id=0; // column of id, 0 based
	bool di = false;
	bool tri = false;
    bool p1 = false;
    bool p2 = false;
    bool h2 = false;
    bool h3 = false;
    bool h4 = false;
	int seqL = 0;
	int start = 1; // position of letter 1
	
	string freq_feat = "";
	string posi_feat = "";
	string output = "output";

	if (argc < 1) help();
	
    // parse arguments
    string str;
    for (int i = 1; i < argc; i++) { 
        if (i != argc) { 
            str=argv[i];
            if (str == "-i") {  
                input = argv[i + 1];
                i=i+1;
	        } else if (str == "-o") {  
				output = argv[i + 1];
				i=i+1;
			} else if (str == "-cSeq"){
				col = stoi(argv[i+1]) - 1;
				i++;
				if(col<0)message("invalid value for -c:"+to_string(col));
			} else if (str == "-cID"){
				col_id = stoi(argv[i+1]) - 1;
				i++;
				if(col_id<0)message("invalid value for -c:"+to_string(col_id));
			} else if (str == "-di"){
				di = true;
            } else if (str == "-tri"){
                tri = true;
            } else if (str == "-p1"){
                p1 = true;
            } else if (str == "-p2"){
                p2 = true;
            } else if (str == "-h2"){
                h2 = true;
            } else if (str == "-h3"){
                h3 = true;
            } else if (str == "-h4"){
                h4 = true;
            } else if (str == "-freq"){
                freq_feat = argv[i+1];
				i++;
            } else if (str == "-posi"){
                posi_feat = argv[i+1];
				i++;
            } else if (str == "-start"){
                start = stoi(argv[i+1]);
				i++;
            } else if (str == "-l"){
                seqL = stoi(argv[i+1]) ;
                i++;
            } else if (str == "-h" || str == "--help") { 
                help();
            } else {
                message("unknown option: "+str);   
                help();        
            }
        }
    }   

	if(freq_feat.size()>0 || posi_feat.size()>0)
	{
		sequence_feature(input, output, freq_feat, posi_feat, col, col_id, start);
		return 1;
	}

	string header = "Seq\tAn.L\tAn.P\tCn.L\tCn.P\tGn.L\tGn.P\tTn.L\tTn.P\tRn.L\tRn.P\tYn.L\tYn.P\tWn.L\tWn.P\tSn.L\tSn.P\tA\tC\tG\tT\tW\tS\tR\tY";
	if(di)	header += "\tAA\tAC\tAG\tAT\tCA\tCC\tCG\tCT\tGA\tGC\tGG\tGT\tTA\tTC\tTG\tTT";
	if(tri) header += "\tAAA\tCAA\tGAA\tTAA\tACA\tAGA\tATA\tCCA\tCGA\tCTA\tGCA\tGGA\tGTA\tTCA\tTGA\tTTA\tAAC\tAAG\tAAT\tCAC\tCAG\tCAT\tGAC\tGAG\tGAT\tTAC\tTAG\tTAT\tACC\tACG\tACT\tAGC\tAGG\tAGT\tATC\tATG\tATT\tCCC\tCCG\tCCT\tCGC\tCGG\tCGT\tCTC\tCTG\tCTT\tGCC\tGCG\tGCT\tGGC\tGGG\tGGT\tGTC\tGTG\tGTT\tTCC\tTCG\tTCT\tTGC\tTGG\tTGT\tTTC\tTTG\tTTT";
    if(p1)
	{
		if(seqL < 2) message("invalid sequence length -l:"+to_string(seqL));
		else
		{
			for(int i =0 ;i<seqL;i++)
			{
				int pos = i+2 - start;
				if(pos <= 0 ) pos = pos - 1;
				header += "\tA."+to_string(pos)+"\tC."+to_string(pos)+"\tG."+to_string(pos)+"\tT."+to_string(pos);
			}
		}
	}
	if(p2)
    {   
        if(seqL < 2) message("invalid sequence length -l:"+to_string(seqL));
        else
        {  
	        array<string,16> di = {"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
            for(int i =0 ;i< seqL-1;i++)
            {
                int pos = i+2 - start;
                if(pos <= 0 ) pos = pos - 1;
				for(int j=0;j<16;j++)
                header += "\t"+di[j]+"."+to_string(pos);
            }
        }
    }

    if(h2)
    {
        if(seqL < 2) message("invalid sequence length -l:"+to_string(seqL));
        else
        {
            array<string,4> homo = {"AA","CC","GG","TT"};
            for(int i =0 ;i< seqL-1;i++)
            {
                int pos = i+2 - start;
                if(pos <= 0 ) pos = pos - 1;
                for(int j=0;j<4;j++)
                header += "\t"+homo[j]+"."+to_string(pos);
            }
        }
    }

    if(h3)
    {   
        if(seqL < 3) message("invalid sequence length -l:"+to_string(seqL));
        else 
        {   
            array<string,4> homo = {"AAA","CCC","GGG","TTT"};
            for(int i =0 ;i< seqL-2;i++)
            {  
                int pos = i+2 - start;
                if(pos <= 0 ) pos = pos - 1;
                for(int j=0;j<4;j++)
                header += "\t"+homo[j]+"."+to_string(pos);
            }
        }
    }

    if(h4)
    {   
        if(seqL < 4) message("invalid sequence length -l:"+to_string(seqL));
        else 
        {   
            array<string,4> homo = {"AAAA","CCCC","GGGG","TTTT"};
            for(int i =0 ;i< seqL-3;i++)
            {  
                int pos = i+2 - start;
                if(pos <= 0 ) pos = pos - 1;
                for(int j=0;j<4;j++)
                header += "\t"+homo[j]+"."+to_string(pos);
            }
        }
    }
	cout << header << endl;

	ifstream fin(input.c_str());
	
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

		sequence_feature(flds[col],di,tri,p1,p2,h2,h3,h4);
	}
	fin.close();

    return 0;
}
