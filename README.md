# Quantum Annealing Single-qubit Assessment (QASA)
QASA is a protocol for evaluating qubit performance across a quantum annealing device. Parameters such as effective temperature (β), input field bias (b), saturation (γ), and flux noise (η) are calculated by fitting the output statistics of each individual qubit. In particular, `collect_spin_table.py` applies a range of h values to each spin and counts the number of times the spin in the z basis is measured as -1 for each h value. `reconstruct_parameters.jl` uses maximum likelihood estimation to fit these empirical statistics to a parameterized model of the qubit as a mixture of quantum Gibbs distributions.

## Installation
QASA requires the python package `dwave-cloud-client v0.8` for connecting to D-Wave API servers. The julia packages `JuMP v0.21`, `Ipopt v0.6`, `DataFrames v0.22`, and `CSV v0.8` are needed as well.


## Basic Usage

### Data Collection
The script `collect_spin_table.py` is used to collect the QASA dataset. For example, the following command will collect statistics on spins 304, 305, 306, and 307 of the DW_2000Q machine for h values ranging from -1 to 1 with step size of 0.025. 1,000,000 samples are taken for each h value, an annealing time of 5 µs is used, every 100 samples a spin reversal transform is performed and if the call to D-Wave's API servers takes longer than 300 seconds, then the call will be resubmitted. The data will be outputted into  `qasa_data/spin_table.csv`, which will be created if it does not already exist. Note that all command line options have default values and are optional with the exception of the directory (-d) and the QPU profile (-p). Also note that if no spins are given with -ss, then the protocol will perform the data collection for the entire chip.
```
python collect_spin_table.py -p DW_2000Q -hr 1.0 -hs 0.025 -d qasa_data -nr 1000000 -ss 304 305 306 307 -at 5 -srtr 100 -to 300
```

In this example, the output csv file will have the column headers `h`, `samples`, `spin_304`, `spin_305`, `spin_306`, `spin_307`.

### Parameter Reconstruction
The script `reconstruct_parameters.jl` performs the maximum likelihood estimation to learn the parameters: effective temperature (β), input field bias (b), saturation (γ), and flux noise (η). The following command will use the dataset `qasa_data/spin_table.csv` and divide the computation among 4 threads. The learned parameters will be outputted in `qasa_data/model_parameters.csv`, which will have the column headers `spin`, `β`, `b`, `γ`, and `η`.
```
 julia -p 4 reconstruct_parameters.jl -d qasa_data
```

## License
QASA is provided under a BSD-ish license with a "modifications must be indicated" clause.  See the `LICENSE.md` file for the full text.
This package is part of the Hybrid Quantum-Classical Computing suite, known internally as LA-CC-16-032.
