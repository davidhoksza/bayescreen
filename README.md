# Bayescreen
Easy-to-use command line virtual screening tool utilizing machine learning to build activity-based model given a set of input compounds.
The model is built over features of molecular fragments features and uses Naive Bayes classifier coupled with likelihood ratio scoring to prioritize the compound library to be screened.

## Installation

### Requirements

- Python 2.7 or 3.5
- RDkit
- PaDEL (only when you want to use it to generate fragments features; otherwise RDKit will be used)

### Download Bayescreen

```
git clone https://github.com/davidhoksza/bayescreen.git
```

## Usage

Bayescreen consists of three utilities:

* Build a model of activity based on features of fragments on known actives and inactives
* Screen a compound library against a model
* Analyze a model

To carry out the screening, one needs a set of known active and inactive molecules and a set of molecules (either in SDF or SMILES format) to screen. The following commands show an example usage of Bayescreen to screen the [CDK2](https://en.wikipedia.org/wiki/Cyclin-dependent_kinase_2) enzyme. The dataset is part of the collections from the [lbvs-environment](https://github.com/skodapetr/lbvs-environment) project, but for convenience the data of the CDK2 dataset are distribued together with Bayescreen. An archive with the SDF files with the molecules are available in the test directory. To unzip them, run:

```
 tar xzvf test/cdk.tar.gz -C test
 ```

### Model construction
To build a model based on known actives stored in `test/a_cdk_train.sdf` and inactives stored in `test/i_cdk_train.sdf` using RDKit for generating fragment features and store it in `model.bm` run:

```
python build_model.py -a test/a_cdk_train.sdf -i test/i_cdk_train.sdf -m test/model.bm -c
```

The `model.bm` is a [JSON](https://en.wikipedia.org/wiki/JSON) file which stores the probabilities used to build the Naive Bayes classifier, imputation and normalization values for every feature and other informations about the model such as which feature generator was used (RDKit or PaDEL) or binning information.

If the fragments and corresponding features are needed, remove the `-c` (`--clean`) option.

### Screening

Screening of the compounds ranks the compounds based on their probability of being active against the same target as the `a_train.sdf` molecules, i.e. the molecules on which the model was built. 

In order to also evaluate the perfromance of Bayescreen we use the split to actives and inactives provided by the XXX. Thus we separately screen the `a_test.sdf` molecules and `i_test.sdf` corresponding to the active and inactive molecules, respectively (obviously, this division will not be known in real VS campaign). both of the

To screen the compounds run:

```
python screen.py -m test/model.bm -d test/a_cdk_test.sdf -o test/a_test.out
python screen.py -m test/model.bm -d test/i_cdk_test.sdf -o test/i_test.out
```

The `x_test.out` contains the screening library sorted by decreasing probability of those molecules being active (average probability ratios of their fragments feature vectors).

To evaluate the performance, we provide a utility which utilizes RDKit to obtain the area under the ROC curve (AUC) and enrichment factor (EF). To get the performance, run:

```
python eval.py -a test/a_test.out -i test/i_test.out
```

### Model analysis

To analyze the model run:

```
python analyze.py -m model.bm -f
```

The output includes the following information:

* General model information (list features, used fragment types, feature generator)
* Importance of feature values in decreasing order. If the `-f` option is turned on, all the values are present, otherwise only the top 50 are listed.

### Contributing

We welcome any bug reports, enhancement requests, and other contributions. To submit a bug report or enhancement request, please use the GitHub issues tracker or 
contact directly David Hoksza.



