#!/bin/bash

./main 0.0 | tee 0.0.out
./main 0.001 | tee 0.001.out
./main 0.0001 | tee 0.0001.out
./main 0.00001 | tee 0.00001.out
