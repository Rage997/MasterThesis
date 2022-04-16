
# Installation

(Export env with conda env export --name invasive-species --file environment.yml)

conda env create--file  environment.yml


# Implementation Problems:

1. Run R code with my data

2. Can't load the whole dataset into memory. 

```numpy.core._exceptions.MemoryError: Unable to allocate 215. GiB for an array with shape (4510, 23191, 276) and data type float64```

I should try to use a sparse matrix.


# TODO

Using Igor's starting points

1. Use 5 years rather than 10
2. 