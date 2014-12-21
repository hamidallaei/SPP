#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>

using namespace std;

template <class T>
class Stat{
public:
	vector<T> data;
	double mean, std, error, variance;
	void Compute();
	void Reset();
	void Add_Data(T input);
};

template <class T> void Stat<T>::Compute()
{
	double sum = accumulate(data.begin(), data.end(), 0.0);
	mean = sum / data.size();
	double sq_sum = inner_product(data.begin(), data.end(), data.begin(), 0.0);
	double sq_mean = sq_sum / data.size();
	variance = sq_mean - mean*mean;
	std = sqrt(variance);
	error = sqrt(variance / data.size());
}

template <class T> void Stat<T>::Reset()
{
	data.clear();
}

template <class T> void Stat<T>::Add_Data(T input)
{
	data.push_back(input);
}

