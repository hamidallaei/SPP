#ifndef _STATISTICS_
#define _STATISTICS_

#include <iostream>
#include <numeric>
#include <vector>
#include <cmath>
#include <boost/algorithm/string.hpp>

using namespace std;

template <class T>
class Stat{
public:
	vector<T> data;
	double mean, mean_square, std, error, variance, min, max, corr_len;
	void Compute();
	void Correlation();
	void Find_Correlation_Length();
	void Shift_Average();
	void Reset();
	void Add_Data(T input);
	void Histogram(const int& num_bins, const string& info);
	void Periodic_Transform(const double& value);
	template <class Tp>  friend std::ostream& operator<<(std::ostream&, Stat<Tp>&);
};

template <class Tp> std::ostream& operator<<(std::ostream& os, Stat<Tp>& s)
{
	for (int i = 0; i < s.data.size(); i++)
		os << s.data[i] << endl;
	return (os);
}

template <class T> void Stat<T>::Compute()
{
	min = data[0];
	max = min;
	for (int i = 0; i < data.size(); i++)
	{
		if (max < data[i])
			max = data[i];
		if (min > data[i])
			min = data[i];
	}
	double sum = accumulate(data.begin(), data.end(), 0.0);
	mean = sum / data.size();
	double sq_sum = inner_product(data.begin(), data.end(), data.begin(), 0.0);
	mean_square = sq_sum / data.size();
	variance = mean_square - mean*mean;
	std = sqrt(variance);
	Find_Correlation_Length();
}

template <class T> void Stat<T>::Correlation()
{
	for (int tau = 0; tau < data.size()/2; tau++)
	{
		T c = 0;
		for (int i = 0; i < (data.size() - tau); i++)
			c += (data[i] - mean)*(data[i+tau]-mean) / (data.size() - tau);
		c /= variance;
		corr_len = tau;
		cout << tau << "\t" << c << endl;
	}
}

template <class T> void Stat<T>::Find_Correlation_Length()
{
	for (int tau = 0; tau < data.size()/2; tau++)
	{
		T c = 0;
		for (int i = 0; i < (data.size() - tau); i++)
			c += (data[i] - mean)*(data[i+tau]-mean) / (data.size() - tau);
		c /= variance;
		corr_len = tau;
		if (c < 0.01)
			tau = data.size();
//		cout << tau << "\t" << c << endl;
	}
	error = sqrt(variance*corr_len / data.size());
}

template <class T> void Stat<T>::Shift_Average()
{
	Compute();
	for (int i = 0; i < data.size(); i++)
		data[i] -= mean;
}

template <class T> void Stat<T>::Reset()
{
	data.clear();
}

template <class T> void Stat<T>::Add_Data(T input)
{
	data.push_back(input);
}

template <class T> void Stat<T>::Histogram(const int& num_bins, const string& info)
{
	double p[num_bins] = {0};
	double bin_width = (max - min) / num_bins;
	for (int i = 0; i < data.size(); i++)
	{
		int index = floor(round((data[i] - min) / bin_width));
		p[index] += 1.0 / (data.size()*bin_width);
	}
	stringstream address("");
	address << info;
	ofstream out_file(address.str().c_str());
	for (int i = 0; i < num_bins; i++)
	{
		out_file << (min + i*bin_width) << "\t" << p[i] << endl;
	}
	double px = 0, py = 0, polarization;
	for (int i = 0; i < num_bins; i++)
	{
		px += bin_width*cos(min + i*bin_width)*p[i];
		py += bin_width*sin(min + i*bin_width)*p[i];
	}
	polarization = sqrt(px*px+py*py);
	cout << "Polarization: " << polarization << endl;
}

template <class T> void Stat<T>::Periodic_Transform(const double& value)
{
	for (int i = 0; i < data.size(); i++)
	{
		data[i] -= 2*value*floor((data[i]+value) / (2*value));
	}
}

#endif
