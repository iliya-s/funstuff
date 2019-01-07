#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Statistics.h"

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
    
    //sample average
    double e_bar = Average(e,t);
    cout << "sample average:\t" << e_bar << endl;

    //variance
    double var = Variance(e,t);
    cout << "sample variance:\t" << var << endl;
    cout << "sample Neff:\t" << Neff(t) << endl;

    cout << "############################### Brute Force ###############################" << endl;
    //correlation function
    vector<double> c;
    CorrFunc(c,e,t);
    //autocorrelation time
    double t_corr = CorrTime(c);
    cout << "autocorrelation time:\t" << t_corr << endl;
    writeCorrFunc(c);
    

    //blocking
    cout << "############################### Blocking ###############################" << endl;
    vector<double> r_x, block_size;
    Block(block_size, r_x, e, t);
    //autocorrelation time
    t_corr = CorrTime(e.size(), block_size, r_x);
    cout << "autocorrelation time:\t" << t_corr << endl;
    writeBlock(block_size, r_x);
    cout << endl;

    cout << "variance of average:\t" << t_corr * var / Neff(t) << endl;

    cout << "############################### Sampling every tau ###############################" << endl;
    int tau = (int) t_corr;
    double E = 0.0;
    vector<double> test_e, test_t;
    for (int i = 0; i < e.size(); i += tau)
    {
        test_e.push_back(e[i]);
        test_t.push_back(t[i]);
    }
    //sample average
    double test_e_bar = Average(test_e,test_t);
    cout << "average:\t" << test_e_bar << endl;

    //variance
    double test_var = Variance(test_e,test_t);
    cout << "variance:\t" << test_var << endl;

    //correlation function
    vector<double> test_c;
    CorrFunc(test_c,test_e,test_t);
    //autocorrelation time
    double test_t_corr = CorrTime(test_c);
    cout << "autocorrelation time:\t" << test_t_corr << endl;
    cout << "Neff:\t" << Neff(test_t) << endl;
    cout << endl;

    cout << "variance of average:\t" << test_t_corr * test_var / Neff(test_t) << endl;
}


