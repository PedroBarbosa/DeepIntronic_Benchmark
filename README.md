# Benchmarking computational methods to predict deep intronic variation

Repository with the data and code to reproduce the evaluation of computational methods to predict functional intronic variation described in the paper:
https://academic.oup.com/gigascience/article/doi/10.1093/gigascience/giad085/7329463

## Datasets
`data` folder contains all the datasets in the hg19 version used for all comparisons. VCF files are annotated with predictions from the tools included in the study.\
`data_hg38` folder contains the main datasets in the hg38 version (via liftover). 

## Code
`intronic_benchmark.ipynb` contains the code to reproduce most of the analysis. Main results were generated with [VETA](https://github.com/PedroBarbosa/VETA), a tool to evaluate the performance of variant effect prediction tools in multiple contexts.\
`scripts` folder contains the post-processing scripts used to create main figures. These scripts take as input results generated by the `intronic_benchmark.ipynb` notebook.

## Reproducing the results
To fully reproduce the results we recommend using Docker. In the command line simply run:
```
git clone https://github.com/PedroBarbosa/DeepIntronic_Benchmark
cd DeepIntronic_Benchmark
./build_docker.sh
./run_docker.sh
```

This will create a Docker image with all software (and exact versions) that will run the `intronic_benchmark.ipynb` jupyter notebook. All results will be saved in the `out` folder.

## Contact
psbarbosa@ciencias.ulisboa.pt
