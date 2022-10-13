# CLR
Draft and Code for "Conditional likelihood ratio test with many weak instruments"

This repository contains the code to create all the tables in figures in the latest draft of "Conditional likelihood ratio test with many weak instruments". There are also additional simulations based on Angrist and Frandsen (2022).

For critical value function output, run critical_val_fn.m

For empirical rejection frequency table, run size_table_output.m

For the power graphs, run power_graphs.m

The empirical simulation file has two codes to run: replication_empirical_mKLM_CLR_MCLR_tHLI_JAK.m and power_empirical_mKLM_CLR_MCLR_tHLI_JAK.m. The first generates size tables, and the second generates power curves for the simulated data, which is based on the Angrist and krueger (1991) dataset. The STATA file creates the simulated dataset from the original dataset.

