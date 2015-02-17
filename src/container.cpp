#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <iomanip>
#include <array>        // std::array

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <algorithm>

#include "container.h"

// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;


// print a matrix
void print_matrix(boost::numeric::ublas::matrix<double> m)
{
    for (unsigned i = 0; i < m.size1(); ++ i)
    {
        for (unsigned j = 0; j < m.size2(); ++ j)
        {
            cout << m(i,j) << "\t";
        }
        cout << endl;
    }
}


double matrix_min(boost::numeric::ublas::matrix<double> x)
{
	double v = 1e100;
    for (int i = 0; i < x.size1(); ++ i)
        for (int j = 0; j < x.size2(); ++ j)
			if(x(i,j) < v)
				v = x(i,j);
	return v;
}

double matrix_max(boost::numeric::ublas::matrix<double> x)
{
	double v = -1e100;
    for (int i = 0; i < x.size1(); ++ i)
        for (int j = 0; j < x.size2(); ++ j)
			if(x(i,j) > v)
				v = x(i,j);
	return v;
}

// max of each column
vector<double> matrix_column_max(boost::numeric::ublas::matrix<double> x)
{
	vector<double> res;
    for (int j = 0; j < x.size2(); ++ j)
	{
		res.push_back(x(0,j));
		for (int i = 1; i < x.size1(); ++ i)
		{
			if(x(i,j) > res[j]) res[j] = x(i,j);
		}
	}
	return res;
}

// min of each column
vector<double> matrix_column_min(boost::numeric::ublas::matrix<double> x)
{
	vector<double> res;
    for (int j = 0; j < x.size2(); ++ j)
	{
		res.push_back(x(0,j));
		for (int i = 1; i < x.size1(); ++ i)
		{
			if(x(i,j) < res[j]) res[j] = x(i,j);
		}
	}
	return res;
}

// min of each column
vector<double> matrix_column_information_content(boost::numeric::ublas::matrix<double> x)
{
	vector<double> res;
    for (int j = 0; j < x.size2(); ++ j)
	{
		res.push_back(log2(x.size1()));
		for (int i = 0; i < x.size1(); ++ i)
		{
			if( x(i,j) > 0 && x(i,j)<= 1) res[j] += x(i,j) * log2(x(i,j));
			else if (x(i,j) < 0 || x(i,j)>1)
			{
				cerr << "ERROR: matrix_column_information_content, element out of range [0,1]" << endl;
				exit(1);
			}
		}
	}
	return res;
}
    
// positive:
// >0: only sum positive
// <0: only sum negative
// =0: both 
vector<double> matrix_column_sum(boost::numeric::ublas::matrix<double> x, int positive)
{
	vector<double> res;
    for (int i = 0; i < x.size2(); ++ i) // column
	{
		res.push_back(0);
        for (int j = 0; j < x.size1(); ++ j) // row
		{
			if( positive > 0)
			{
				if (x(j,i)>0) res[i] += x(j,i);
			}
		    else if (positive <0)
			{
				if (x(j,i)<0) res[i] += x(j,i);
			}
			else
				res[i] += x(j,i);
			//cout << i << "\t" << res[i] << endl;
		}
	}
	return res;
}
