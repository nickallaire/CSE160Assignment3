# CSE 160: Introduction to Parallel Computing Assignment 3

Project Overview: The goal of this assignment is to parallelize provided serial code that multplies two matrices together through the process of master-worker. The matrix multiplication is performed by sending blocks to each worker to perform a "smaller" matrix multiply and then send the results back to the master who stores the results. Bench marking is also performed with message passing to learn about how message size affects the bandwidth and latency for unidirectional and bidirectional message sends.

A complete description of the assignment can be found in PR3.pdf
