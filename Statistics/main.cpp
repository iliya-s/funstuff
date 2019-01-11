#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Statistics.h"


int main()
{
    Statistics Stats;
    Stats.Block();
    Stats.BlockCorrTime();
    Stats.WriteData();
    Stats.push_back(1, 2);
}


