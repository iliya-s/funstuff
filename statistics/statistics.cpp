#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include "statistics.h"

using namespace std;

int main()
{
    //read data
    ifstream efile("E_loc.bin", ios::binary);
    efile.seekg(0, ios::end);
    long size = efile.tellg();
    efile.seekg(0, ios::beg);
    vector<double> e(size/sizeof(double));
    efile.read((char *)e.data(), size);
    efile.close();

    ifstream tfile("tau.bin", ios::binary);
    vector<double> t(size/sizeof(double));
    tfile.read((char *)t.data(), size);
    tfile.close();
    
    cout << "############################### Brute Force ###############################" << endl;
    //sample average
    double e_bar = average(e,t);
    //cout << "average:\t" << e_bar << endl;

    //variance
    double var = variance(e,t);
    //cout << "variance:\t" << var << endl;

    //correlation function
    vector<double> c;
    corr_func(c,e,t);
    //autocorrelation time
    double t_corr = corr_time(c);
    cout << "autocorrelation time:\t" << t_corr << endl;
    write_corr_func(c);
    

    //blocking
    cout << "############################### Blocking ###############################" << endl;
    vector<double> r_x, block_size;
    block(block_size, r_x, e, t);
    //autocorrelation time
    t_corr = corr_time(e.size(), block_size, r_x);
    cout << "autocorrelation time:\t" << t_corr << endl;
    write_block(block_size, r_x);
}


