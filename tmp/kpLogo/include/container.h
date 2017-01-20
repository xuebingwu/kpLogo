#ifndef __CONTAINER_H__
#define __CONTAINER_H__

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



// for position matrix
#include <boost/numeric/ublas/matrix.hpp>

using namespace std;


// merge two vectors
template <typename T>
std::vector<T> operator+(const std::vector<T> &A, const std::vector<T> &B)
{
    std::vector<T> AB;
    AB.reserve( A.size() + B.size() );                // preallocate memory
    AB.insert( AB.end(), A.begin(), A.end() );        // add A;
    AB.insert( AB.end(), B.begin(), B.end() );        // add B;
    return AB;
}

template <typename T>
std::vector<T> &operator+=(std::vector<T> &A, const std::vector<T> &B)
{
    A.reserve( A.size() + B.size() );                // preallocate memory without erase original data
    A.insert( A.end(), B.begin(), B.end() );         // add B;
    return A;                                        // here A could be named AB
}

// sort return index
template <typename T>
vector<size_t> sort_indexes(const vector<T> &v, int ascending = 1) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;	

  // return index without sorting
  if (ascending==0) return idx;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  if(ascending < 0) reverse(idx.begin(),idx.end());

  return idx;
}


// return the rank
template <typename T>
vector<double> get_ranks(vector<T> &v, bool descending=false)
{
	vector<size_t> s = sort_indexes(v);

	vector<double> rk(v.size(),0) ;
	for(int i=0;i<s.size();i++) 
	{
		rk[s[i]] = i+1;
	}
	// ties
	int tie_start = 0;
	for(int i=1;i<s.size();i++)
	{
		if(v[s[i]] != v[s[i-1]]) // a change in value
			{
				// tie from tie_start to i-1
				double new_rank = (tie_start + i + 1)/2.0;
				for(int j= tie_start;j<i;j++)
				{
					rk[s[j]] = new_rank;
				}
				tie_start = i; 
			}
			else if (i==s.size()-1)
			{
				double new_rank = (tie_start + i + 2)/2.0;
				for(int j= tie_start;j<=i;j++)
				{
					rk[s[j]] = new_rank;
				}
			}
	}
	return rk;
}
	
	

vector<double> matrix_column_information_content(boost::numeric::ublas::matrix<double> x);
vector<double> matrix_column_sum(boost::numeric::ublas::matrix<double> x, int positive=0);
vector<double> matrix_column_min(boost::numeric::ublas::matrix<double> x);
vector<double> matrix_column_max(boost::numeric::ublas::matrix<double> x);

// print a matrix
void print_matrix(boost::numeric::ublas::matrix<double> m);

void save_matrix(boost::numeric::ublas::matrix<double> m, string filename);

double matrix_min(boost::numeric::ublas::matrix<double> x);

double matrix_max(boost::numeric::ublas::matrix<double> x);

// sum of two vectors
template<typename T>
vector<T> sum(vector<T> a, vector<T> b)
{
    vector<T> res;
    for (T i =0;i< a.size();i++) res.push_back(a[i]+b[i]);
    return res;
}



template<class key,class val>
void PrintMap(map<key,val> m)
{// print out a map
  for (typename map<key,val>::iterator it = m.begin();it!=m.end(); it++)
  {
    cout << (*it).first << "\t=>\t" <<(*it).second<<endl;
  }
}

#endif

