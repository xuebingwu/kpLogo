#include "markov.h"
#include "sequence.h"
#include "utility.h"

using namespace std;


// empty constructor
markov_model::markov_model()
{
}

// construct a markov model from a string
// string format: A:24,C:23,G:23,T:30,AC:2,AG:3,....
markov_model::markov_model(string alphabet, string model)
{
	this->alphabet = alphabet;
	vector<string> flds = string_split(model,",");

	//determine the order of the model based on length
	int n0 = alphabet.size();
	int n1 = n0 + n0*n0;
	int n2 = n1 + n0*n0*n0;

	// at least specify frequency of each letter
	if(flds.size() == n0) this->order = 0;
	else if (flds.size() == n1) this->order = 1;
	else if (flds.size() == n2) this->order = 2;
	else
	{
		message("ERROR: markov model format error: order "+model);
		system_run("touch exit_with_error");exit(1);
	}

//	cout << model << endl;
//	cout << this->order << endl;


	// initialize the model
	vector<string> di,tri;
	for(int i=0;i<alphabet.size();i++) this->f1[alphabet[i]] = -1;
	if(this->order >0)
	{
        di = generate_kmers(2,alphabet);
        for (int i = 0;i<di.size();i++) this->f2[di[i]] = -1;
	}
    if(this->order >1)
    {   
        tri = generate_kmers(3,alphabet);
        for (int i = 0;i<tri.size();i++) this->f3[tri[i]] = -1;
    }

	//
	for(int i=0;i<flds.size();i++)
	{
		vector<string> tmp = string_split(flds[i],":");
		if(tmp.size() != 2)
		{ 
			message("ERROR: : markov model format error: "+flds[i]);
			system_run("touch exit_with_error");exit(1);
		}  
		double freq = stof(tmp[1]);

//		cout << flds[i] << "," << freq << endl;

		if(freq < 0)
		{
			message("ERROR: markov model error: negative value: "+ tmp[1]);
			system_run("touch exit_with_error");exit(1);
		}
		if(tmp[0].size() == 1)
		{
			size_t found = alphabet.find(tmp[0]);
			if(found == std::string::npos)
			{ 
				message("ERROR: : markov model format error: "+tmp[0]+" not in alphabet "+alphabet);
				system_run("touch exit_with_error");exit(1);
			}  	
			this->f1[alphabet[found]] = freq;
			//cout << alphabet[found] << "->" << this->f1[alphabet[found]] << endl;
		}
		else if((tmp[0].size() == 2) && (this->order > 0))
		{
			if(this->f2.find(tmp[0]) == f2.end())
            {   
                message("ERROR: : markov model key error: "+tmp[0]+" is not two-letter of "+alphabet);
                system_run("touch exit_with_error");exit(1);
            }
			this->f2[tmp[0]] = freq;
		}
        else if((tmp[0].size() == 3) && (this->order > 1))
        {
            if(this->f3.find(tmp[0]) == f3.end())
            {
                message("ERROR: : markov model key error: "+tmp[0]+" is not three-letter of "+alphabet);
                system_run("touch exit_with_error");exit(1);
            }
            this->f3[tmp[0]] = freq;
        }
	}
	
	// normalization and build model
	int L = alphabet.size(); // number of letters

    // zero order
    double total_letter_count = 0;
    for (int i = 0;i<L;i++)
    {
        total_letter_count += this->f1[alphabet[i]];
    }
    // divide by total to get frequency
    for (int i = 0;i<L;i++) this->f1[alphabet[i]] /= total_letter_count;

    if(this->order>0)
    {
        // first order
        double total_di_count = 0;
		for (int i=0;i<di.size();i++)
        {
            total_di_count += this->f2[di[i]];
        }
        // prob of each di
        for (int i = 0;i<di.size();i++) this->f2[di[i]] /= total_di_count;

        //for (int i = 0;i<di.size();i++) cout << di[i] << "\t" << p2[di[i]] << endl;

        // prob of each f2 given f1
        for (int i =0;i<L;i++)
        {
            for(int j =0;j<L;j++)
            {
                this->first_order[alphabet[i]][alphabet[j]] = this->f2[string(1,alphabet[i])+string(1,alphabet[j])] / this->f1[alphabet[i]] ;
                //cout << alphabet[i] << "->" << alphabet[j] << "\t" << first_order[alphabet[i]][alphabet[j]] << endl;
            }
        }

        if(this->order == 2 )
        {
            // second order
            double total_tri_count = 0; // note doesn't equal to total sequence length * 3, consider non-triplet sequences 
            for (int i = 0;i<tri.size();i++)
            {
                total_tri_count += this->f3[tri[i]];
            }
            // prob of each tri
            for (int i = 0;i<tri.size();i++) this->f3[tri[i]] /=  total_tri_count;

            // prob of each f3 given f2
            for (int i =0;i<di.size();i++)
            {
                for(int j =0;j<alphabet.size();j++)
                {
                    this->second_order[di[i]][alphabet[j]] = this->f3[di[i]+string(1,alphabet[j])] / this->f2[di[i]] ;
                }
            }
        }
    }	
}

// construct a n-th order markov model from sequences 
markov_model::markov_model(int n, string alphabet, vector<string> seqs)
{
    this->order = n; // n = 0,1,2
    this->alphabet = alphabet;
    int L = alphabet.size(); // number of letters
    // zero order
    vector<int> letter_count;
    int total_letter_count = 0;
    for (int i = 0;i<L;i++) 
    {
        letter_count.push_back(  countSubstringInSeqs(seqs,string(1,alphabet[i])) ); 
        total_letter_count += letter_count[i];  
    }
    // divide by total to get frequency
    for (int i = 0;i<L;i++) f1[alphabet[i]] = double(letter_count[i])/total_letter_count;
    
    //for (int i = 0;i<L;i++) cout << alphabet[i] << "\t" << p1[alphabet[i]] << endl;
    
    if(n>0)
    {
        // first order
        // count all dinucleotides
        vector<string> di = generate_kmers(2,alphabet);
        int total_di_count = 0;
        for (int i = 0;i<di.size();i++)
        {
            f2[di[i]] = countSubstringInSeqs(seqs,di[i]);
            total_di_count += f2[di[i]];
        } 
        // prob of each di
        for (int i = 0;i<di.size();i++) f2[di[i]] = f2[di[i]] / total_di_count;

        //for (int i = 0;i<di.size();i++) cout << di[i] << "\t" << p2[di[i]] << endl;

        // prob of each f2 given f1
        for (int i =0;i<L;i++)
        {
            for(int j =0;j<L;j++)
            {
                first_order[alphabet[i]][alphabet[j]] = f2[string(1,alphabet[i])+string(1,alphabet[j])] / f1[alphabet[i]] ; 
                //cout << alphabet[i] << "->" << alphabet[j] << "\t" << first_order[alphabet[i]][alphabet[j]] << endl;
            }
        }
        if(n==2)
        {
            // second order
            // count all trinucleotides
            vector<string> tri = generate_kmers(3,alphabet);
            int total_tri_count = 0; // note doesn't equal to total sequence length * 3, consider non-triplet sequences 
            for (int i = 0;i<tri.size();i++)
            {
                f3[tri[i]] = countSubstringInSeqs(seqs,tri[i]);
                total_tri_count += f3[tri[i]];
            }
            // prob of each tri
            for (int i = 0;i<tri.size();i++) f3[tri[i]] = f3[tri[i]] / total_tri_count;

            // prob of each f3 given f2
            for (int i =0;i<di.size();i++)
            {   
                for(int j =0;j<alphabet.size();j++)
                {
                    second_order[di[i]][alphabet[j]] = f3[di[i]+string(1,alphabet[j])] / f2[di[i]] ;
                }
            }
        }
    }
}

// print the markov model that is compatible with MEME suite
void markov_model::print(string filename)
{
    ofstream out;
    out.open(filename.c_str());
    out << "# order\t" << order << endl;
    out << "# alphabet\t" << alphabet << endl;
    out << "# mono-nucleotide-frequency" << endl;
    for (int i=0;i<f1.size();i++) out << alphabet[i] << "\t" << f1[alphabet[i]] << endl;
    if (order > 0)
    {
        out << "# di-nucleotide-frequency" << endl;
        for(map<string,double>::iterator it=f2.begin();it!=f2.end();it++)
        {
            out << it->first << "\t" << it->second << endl;;
        }
        out << "# transitional probability, first order" << endl;
        for (int i =0;i<alphabet.size();i++)
        {
            for(int j =0;j<alphabet.size();j++)
            {  
                out << "# " <<alphabet[i] << "->" << alphabet[j] << "\t" << first_order[alphabet[i]][alphabet[j]] << endl;
            }
        }
        if (order > 1)
        {
            out << "# tri-nucleotide-frequency" << endl;
            for(map<string,double>::iterator it=f3.begin();it!=f3.end();it++)
            {  
                out << it->first << "\t" << it->second << endl;;
            }
            out << "# transitional probability, second order" << endl;
            for(map<string,double>::iterator it=f2.begin();it!=f2.end();it++)
            {   
                for(int j =0;j<alphabet.size();j++)
                {
                    out << "# " << it->first << "->" << alphabet[j] << "\t" << second_order[it->first][alphabet[j]] << endl;
                }
            }
        }
    }
    out.close();
}

// score all kmers by a markov model, kmers can be of variable lengths, but need to be exact without degenerate bases
map<string,double> markov_model::probs(vector<string> kmers)
{
    map<string,double> p;

    // order 0
    // product of mononucleotide frequency
    if(order == 0)
    {
        for (int i=0;i<kmers.size();i++)
        {
            p[kmers[i]] = f1[kmers[i][0]];
            for (int j=1;j<kmers[i].size();j++) p[kmers[i]] = p[kmers[i]] * f1[kmers[i][j]];
        }
        return p;
    }

    // order 1
    if(order == 1)
    {
        for (int i=0;i<kmers.size();i++)
        {
            int k = kmers[i].size();
            if (k == 1) p[kmers[i]] = f1[kmers[i][0]]; 
            else if (k==2) p[kmers[i]] = f2[kmers[i]];
            else
            {
                p[kmers[i]] = f2[kmers[i].substr(0,2)];
                for (int j=2;j<k;j++) p[kmers[i]] = p[kmers[i]] * first_order[kmers[i][j-1]][kmers[i][j]];
            } 
        }
        return p;
    }

    // order 2
    if(order == 2)
    {
        for (int i=0;i<kmers.size();i++)
        {   
            int k = kmers[i].size();
            if (k == 1) p[kmers[i]] = f1[kmers[i][0]];
            else if (k==2) p[kmers[i]] = f2[kmers[i]];
            else if (k==3) p[kmers[i]] = f3[kmers[i]];
            else
            {
                p[kmers[i]] = f3[kmers[i].substr(0,3)];
                for (int j=3;j<k;j++) p[kmers[i]] = p[kmers[i]] * second_order[kmers[i].substr(j-2,2)][kmers[i][j]];
            }
        }
        return p;
    }

    // other orders
    cerr << "Not supported: order = " << order << endl;
    system_run("touch exit_with_error");exit(1);
}

