#include <boost/math/distributions/binomial.hpp>
  using boost::math::binomial;
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/normal.hpp> // for normal_distribution
  using boost::math::normal; // typedef provides default type is double.

#include <boost/math/distributions/hypergeometric.hpp>
using boost::math::hypergeometric; 

#include <array>
#include <vector>
#include "container.h"
#include "stat.h"
  
  using namespace std;
  
double p_value_for_z_score(double z)
	{
		normal s;
		return cdf(s, -fabs(z));
	} 

// note both enrichment and depletion
double binom_test(int trials,int success,double success_fraction)
{
    binomial binom(trials, success_fraction);
    double p = 1.0 - cdf(binom,success);
    if (p > 0.5) p = 1.0 - p; //depletion
    return p;
}

double hypergeometric_test(unsigned k, unsigned r, unsigned n, unsigned N)
{	// k: positive samples
	// n: total samples
	// r: total positive
	// N: total 
	hypergeometric hg(r,n,N);
	double p = 1.0 - cdf(hg,k);
	return p;
}

// normal approximation is ok for n1 and n2 > 8
// see paper: Carine A. Bellera et al, Normal Approximations to the Distributions of the Wilcoxon
// Statistics: Accurate to What N? Graphical Insights, Journal of Statistics Education, Volume 18, Number 2, (2010)
// http://www.amstat.org/publications/jse/v18n2/bellera.pdf
array<double,2> Mann_Whitney_U_test(vector<int> ranks, int N)
{
	double n1 = ranks.size();
	double n2 = N - n1;
	double rank_sum = ranks[0];
	for (int i=1;i<n1;i++) rank_sum += ranks[i];
	double n1n2 = n1*n2;
	double U = n1n2+n1*(n1+1)/2-rank_sum;
	double mU = n1n2/2.0;
	double sigmaU = sqrt(n1n2*(N+1)/12.0);
	double z = (U-mU)/sigmaU;
	normal s;
	//cout << n1n2 << "\t" << n1n2*(N+1)/12.0 << "\t" << sigmaU << endl;
	double p = cdf(s, -fabs(z));
	array<double,2> res = {{z,p}};
	return res;
}


array<double,2> two_samples_t_test_equal_sd(double Sm1,double Sd1, unsigned Sn1, double Sm2, double Sd2,unsigned Sn2)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm</a>
   //
   using namespace std;
   using namespace boost::math;

   // Now we can calculate and output some stats:
   //
   // Degrees of freedom:
   double v = Sn1 + Sn2 - 2;
   // Pooled variance:
   double sp = sqrt(((Sn1-1) * Sd1 * Sd1 + (Sn2-1) * Sd2 * Sd2) / v);
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / (sp * sqrt(1.0 / Sn1 + 1.0 / Sn2));
   //
   // Define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));

	array<double,2> res = {{t_stat,q}};

	return res;
}

array<double,2> two_samples_t_test_unequal_sd(double Sm1,double Sd1, unsigned Sn1,double Sm2, double Sd2,unsigned Sn2)
{
   //
   // Sm1 = Sample Mean 1.
   // Sd1 = Sample Standard Deviation 1.
   // Sn1 = Sample Size 1.
   // Sm2 = Sample Mean 2.
   // Sd2 = Sample Standard Deviation 2.
   // Sn2 = Sample Size 2.
   // alpha = Significance Level.
   //
   // A Students t test applied to two sets of data.
   // We are testing the null hypothesis that the two
   // samples have the same mean and that any difference
   // if due to chance.
   // See <a href="http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm">http://www.itl.nist.gov/div898/handbook/eda/section3/eda353.htm</a>
   //
   using namespace std;
   using namespace boost::math;

   //
   // Degrees of freedom:
   double v = Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2;
   v *= v;
   double t1 = Sd1 * Sd1 / Sn1;
   t1 *= t1;
   t1 /=  (Sn1 - 1);
   double t2 = Sd2 * Sd2 / Sn2;
   t2 *= t2;
   t2 /= (Sn2 - 1);
   v /= (t1 + t2);
   // t-statistic:
   double t_stat = (Sm1 - Sm2) / sqrt(Sd1 * Sd1 / Sn1 + Sd2 * Sd2 / Sn2);
   //
   // Define our distribution, and get the probability:
   //
   students_t dist(v);
   double q = cdf(complement(dist, fabs(t_stat)));

	array<double,2> res = {{t_stat,q}};
	return res;
}

// return 6 double
// t.statistics, p-value, Sm1, Sd1,  Sm2, Sd2
array<double,6> t_test(vector<double> x, vector<double> y, bool equal_var)
{
	unsigned Sn1 = x.size();
	double Sm1 = mean(x);
	double Sd1 = sqrt(var(x));

    unsigned Sn2 = y.size();
    double Sm2 = mean(y); 
    double Sd2 = sqrt(var(y)); 
	
	//debug cout << Sn1 << "," << Sm1 << "," << Sd1 << endl;
	//debug cout << Sn2 << "," << Sm2 << "," << Sd2 << endl;	

	array<double,6> res = {{-1,-1,Sm1, Sd1,  Sm2, Sd2}};
	
	array<double,2> ttest;
	if (equal_var) ttest = two_samples_t_test_equal_sd(Sm1,Sd1,Sn1,Sm2,Sd2,Sn2);
	else ttest = two_samples_t_test_unequal_sd(Sm1,Sd1,Sn1,Sm2,Sd2,Sn2);	

	res[0] = ttest[0];
	res[1] = ttest[1];
	
	return res;
}
