# Supplementary files for "Generation of Tree-Child phylogenetic networks"
## By G. Cardona and J.C. Pons

If you have `pipenv` installed, the requirements can be installed by executing `pipenv install`. Once installed, the jupyter 
notebooks can be run with `pipenv run jupyter-notebook`, and the scripts with `pipenv run python script.py`.

The included files are:

* `generating_TC.py`: Python module that generates BTC networks.
* `bounds.py`: Python script that computes the upper bound for the number of BTC networks. It runs without parameters and computes these bounds up to 10 leaves.
* `compute_offspring.py`: Python script that computes the offspring of a set of networks. Usage: `compute_offspring.py infile outfile taxon`. `infile` is a text file containing the eNewick strings of the networks whose offspring has to be computed; `outfile` text file where the eNewick strings of the offspring are written; `taxon` identifier for the new taxon to be used for the generated networks.
* `count_offspring.py`: Python script that computes the cardinality of the offspring of a set of networks. Usage: `count_offspring.py infile`. `infile` is a text file containing the eNewick strings of the networks whose offspring has to be counted.
* `nets_n.txt`: Text file containing the eNewick strings of the BTC networks with `n` leaves (for `n`=1,2,3,4). These files can be used as input for the scripts above.
* `demo.py`: Jupyter notebook that demonstrates the sequential and random generation of BTC networks with arbitrary number of leaves.

