# Model.Microbiome
 
This package creates and subsamples simulated microbial communities. It first generates a reference community and some intrinsic statistics. Then it subsamples to create a "raw" count distribution, simulating a sequencing event. Then it runs sequence normalization using user-supplied methods. Then metrics are collected on each of the normalization methods provided.

Updates in development:
1. The current implementation is very demanding on RAM because it holds a lot of data in RAM while running the simulation. In the next update the function Benchmark.MM will be disaggregated into several functions that allow the user to run the base community generation separate from all the metric calculation steps. The purpose of this is to reduce the chance that the computer hangs up during the simulation because it runs out of RAM.
2. Integrate HMSC as the go-to framework for testing the relationship of the taxa to the environment and to each other.
3. Include explicit taxon interaction parameters in the taxon equations.
