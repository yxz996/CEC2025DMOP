This repository is created for the CEC2025 competition  of dynamic multi-objective optimisation, There are two tracks of competition, briefly described as follows:  
Track 1: Dynamic Unconstrained Multi-Objective Optimisation  
Track 2: Dynamic Constrained Multi-Objective Optimisation  

All the benchmark functions have been implemented in MATLAB code, your competition results can be submitted as a brief technical report.   
Please send your results directly to Dr Xiaozhong Yu (xzyu@ smail.xtu.edu.cn)  

More details can be found in "https://zoujuan1.github.io/#page-top"



# Search space

| Problems | $x_i=1$ | $x_i=2,..n$ |
| -------- | ------- | ----------- |
| DP1      | [0,1]   | [-1,1]      |
| DP2      | [0,1]   | [-1,1]      |
| DP3      | [0,1]   | [0,1]       |
| DP4      | [0,1]   | [-1,1]      |
| DP5      | [0,1]   | [-1,1]      |
| DP6      | [0,1]   | [-1,1]      |
| DP7      | [0,1]   | [-1,1]      |
| DP8      | [0,1]   | [-1,1]      |
| DP9      | [0,1]   | [-1,1]      |
| DP10     | [0,1]   | [0,1]       |
| DC1      | [0,1]   | [-1,1]      |
| DC2      | [0,1]   | [-1,1]      |
| DC3      | [0,1]   | [-1,1]      |
| DC4      | [0,1]   | [-1,1]      |
| DC5      | [0,1]   | [0,1]       |
| DC6      | [0,1]   | [-1,1]      |
| DC7      | [0,1]   | [0,1]       |
| DC8      | [0,1]   | [-1,1]      |
| DC9      | [0,1]   | [-1,1]      |
| DC10     | [0,1]   | [-1,1]      |



# Parameter Setting

- Population Size:100
- Number of variables: 10
- Frequency of change $\tau_t$: 10 (fast changing environments), 20 (slow changing environments).
- Severity of change $n_t$: 5 (severe changing environments), 10 (moderate changing environments).
- Number of changes: 30.
- Number of independent runs: 20
- Stopping criterion: a maximum number of 100(30$\tau_t$+50) fitness evaluations, where 500 fitness evaluations are given before the first environmental change occurs.
- Metrics: MIGD、MHV[1-2]

# Result Submission

It is expected that competition results can be submitted in tables in a format exemplified in Table 1. However, other ways of result presentation are also acceptable. Please do make sure your result is of high readability for submission, and multiple types of results shown in Table are clearly recorded, including the mean and standard deviation of the MIGD/MHV values of each test instance.

| Problem | $(\tau_t,n_t)$                        | MIGD(mean(std.))   | MHV(mean(std.))    |
| ------- | ------------------------------------- | ------------------ | ------------------ |
| DP1     | 10,5 <br> 10,10 <br> 20,5  <br> 20,10 | 1.234E-2(1.234E-3) | 1.234E-2(1.234E-3) |
| DP2     |                                       |                    |                    |
| ……      |                                       |                    |                    |
| DP10    | 10,5 <br> 10,10 <br> 20,5 <br> 20,10  |                    |                    |

# Reference

[1]Hu Y, Zou J, Zheng J, et al. A new framework of change response for dynamic multi-objective optimization[J]. Expert Systems with Applications, 2024, 248: 123344.

[2]Jiang S, Zou J, Yang S, et al. Evolutionary dynamic multi-objective optimisation: A survey[J]. ACM Computing Surveys, 2022, 55(4): 1-47.
