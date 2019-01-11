#ifndef STATS_HEADER_H
#define STATS_HEADER_H
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

//various functions for calculating statistics of serially correlated data, functions overloaded for weighted and unweighted data sets


//average function
//	weighted data
double average(std::vector<double> &x, std::vector<double> &w);
//	unweighted data
double average(std::vector<double> &x);


//calculates effective sample size for weighted data sets. For unweighted data, will just return the size of the data set
double neff(std::vector<double> &w);


//variance function, takes advantage of bessel's correction
//	weighted data
double variance(std::vector<double> &x, std::vector<double> &w);
//	unweighted data
double variance(std::vector<double> &x);


//correlation function: given a weighted/unweighted data set, calculates C(t) = <(x(i)-x_bar)(x(i+t)-x_bar)> / <(x(i)-x_bar)^2>
//Input: x, w if weighted; x if unweighted
//Output: c
void corrFunc(std::vector<double> &c, std::vector<double> &x, std::vector<double> &w);
void corrFunc(std::vector<double> &c, std::vector<double> &x);


//autocorrelation time: given correlation function, calculates autocorrelation time: t = 1 + 2 \sum C(i)
double corrTime(std::vector<double> &c);


//writes correlation function to text file
void writeCorrFunc(std::vector<double> &c);


//blocking function, given a wighted or unweighted data set, calculates autocorrelation time vs. block size
//Input: x, w if weighted; x if unweighted
//Output: b_size - block size per iteration, r_t - autocorrelation time per iteration
void block(std::vector<double> &b_size, std::vector<double> &r_t, std::vector<double> &x, std::vector<double> &w);
void block(std::vector<double> &b_size, std::vector<double> &r_t, std::vector<double> &x);


//autocorrelation time: given blocking data, finds autocorrelation time based on the criteria: (block size)^3 > 2 * (number of original data points) * (autocorrelation time)^2
double corrTime(double n_original, std::vector<double> &b_size, std::vector<double> &r_t);


//writes blocking data to file
void writeBlock(std::vector<double> &b_size, std::vector<double> &r_t);


//class wrapper
class Statistics
{
  public:
    //data
    std::vector<double> X, W;
    //outputs
    double avg, n = -1.0, var;
    std::vector<double> C;
    std::vector<double> B, R;
    double t_corr, t_block;
    
    //append data point    
    int push_back(double x, double w)
    {
      X.push_back(x);
      W.push_back(w);
      return 2;
    }
    int push_back(double x)
    {
      X.push_back(x);
      return 1;
    }

    //write data to file
    void WriteData()
    {
      if (X.size() == 0)
        std::cout << "No data to write" << std::endl;
      else
      {
        std::ofstream xdata("X.bin", std::ios::binary);
        xdata.write((char *)&X[0], X.size() * sizeof(double));
        xdata.close();
        if (X.size() == W.size())
        {
          std::ofstream wdata("W.bin", std::ios::binary);
          wdata.write((char *)&W[0], W.size() * sizeof(double));
          wdata.close();
        }
      }
    }

    //calculates average of data
    double Average()
    {
      if (X.size() == 0)
      {
        std::cout << "No data to average" << std::endl;
        return 0.0;
      }
      else
      {
        if (X.size() == W.size())
          avg = average(X, W);
        else
          avg = average(X);
        return avg;
      }
    }
        
    //calculates effective number of data points
    double Neff()
    {
      if (X.size() == W.size())
        n = neff(W);
      else
        n = (double) X.size();
      return n;
    }

    //calculates variance of data set
    double Variance()
    {
      if (X.size() == 0)
      {
        std::cout << "No data to analyze" << std::endl;
        return 0.0;
      }
      else
      {
        if (X.size() == W.size())
          var = variance(X, W);
        else
          var = variance(X);
        return var;
      }
    }

    //calculates time correlation function
    void CorrFunc()
    {
      if (X.size() == 0)
        std::cout << "No data to analyze" << std::endl;
      else
      {
        if (X.size() == W.size())
          corrFunc(C, X, W);
        else
          corrFunc(C, X);
      }
    }

    //Integrates correlation function to calculate autocorrelation time
    double IntCorrTime() //"Integrated autocorrelation time" - brute force calc
    {
      if (C.size() == 0)
      {
        std::cout << "Run CorrFunc() before integrating for autocorrelation time" << std::endl;
        t_corr = 0.0;
        return t_corr;
      }
      else
      {
        t_corr = corrTime(C);
        return t_corr;
      }
    }
    
    //Writes Correlation function to file
    void WriteCorrFunc()
    {
      if (C.size() == 0)
        std::cout << "Correlation function is empty" << std::endl;
      else
        writeCorrFunc(C);
    }
    
    //reblocking analysis
    void Block()
    {
      if (X.size() == 0)
        std::cout << "No data to analyze" << std::endl;
      else
      {
        if (X.size() == W.size())
          block(B, R, X, W);
        else
          block(B, R, X);
      }
    }
    
    //find autocorrelation time by finding saturated block length
    double BlockCorrTime()
    {
      if (B.size() == 0)
      {
        std::cout << "Run Block() before checking for autocorrelation time convergence" << std::endl;   
        return 0.0;
      }
      else
      {
        if (n == -1)
            Neff();
        t_block = corrTime(n, B, R);
        return t_block;
      }
    }

    //Write block to file
    void WriteBlock()
    {
      if (B.size() == 0 && R.size() == 0)
        std::cout << "Blocking data is empty" << std::endl;
      else 
        writeBlock(B, R);
    }
};
#endif
