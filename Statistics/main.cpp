#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include <chrono>
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
    double e_bar = average(e,t);
    cout << "sample average:\t" << e_bar << endl;

    //variance
    double var = variance(e,t);
    double Neff = neff(t);
    cout << "sample variance:\t" << var << endl;
    cout << "sample Neff:\t" << Neff << endl;

    cout << "############################### Brute Force ###############################" << endl;
    auto start = std::chrono::system_clock::now();
    //correlation function
    std::vector<double> c;
    corrFunc(c,e,t);
    //autocorrelation time
    double t_corr = corrTime(c);
    cout << "autocorrelation time:\t" << t_corr << endl;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "elapsed time: " << elapsed.count() << endl;
    writeCorrFunc(c);
    cout << endl;
    cout << "variance of average:\t" << t_corr * var / Neff << endl;

    //blocking
    cout << "############################### Blocking ###############################" << endl;
    start = std::chrono::system_clock::now();
    vector<double> r_x, block_size;
    block(block_size, r_x, e, t);
    //autocorrelation time
    t_corr = corrTime(Neff, block_size, r_x);
    cout << "autocorrelation time:\t" << t_corr << endl;
    end = std::chrono::system_clock::now();
    elapsed = end - start;
    cout << "elapsed time: " << elapsed.count() << endl;
    writeBlock(block_size, r_x);
    cout << endl;

    cout << "variance of average:\t" << t_corr * var / Neff << endl;

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
    double test_e_bar = average(test_e,test_t);
    cout << "average:\t" << test_e_bar << endl;

    //variance
    double test_var = variance(test_e,test_t);
    cout << "variance:\t" << test_var << endl;

    //correlation function
    vector<double> test_c;
    corrFunc(test_c,test_e,test_t);
    //blocking
    vector<double> B, R;
    block(B, R, test_e, test_t);
    //autocorrelation time
    double test_t_corr = corrTime(test_c);
    cout << "integrated autocorrelation time:\t" << test_t_corr << endl;
    test_t_corr = corrTime(neff(test_t), B,R);
    cout << "blocking autocorrelation time:\t" << test_t_corr << endl;
    cout << "Neff:\t" << neff(test_t) << endl;
    cout << endl;

    cout << "variance of average:\t" << test_t_corr * test_var / neff(test_t) << endl;
}
