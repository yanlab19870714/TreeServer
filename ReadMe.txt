#################################################
Tested with the Following Version
#################################################

MPI - 3.2
GCC - 5.4.0
HADOOP - 2.7.5

Boost C++ Libraries should be installed

#################################################
How to Compile :
#################################################

There are two files to compile:

1. we use our own put command to put a local file to HDFS. It will write one file per data column
- Go insde the "put" directory, you will find the cpp file which you can make to create that program (make clean;make)
- rename to "put" so that you can run "./put local_file HDFS_path"

2. the cpp file to run treeServer program is "run.cpp" in the root root directory
- make clean;make
- run the application on an uploaded dataset
- see below for detailed usage

#################################################
How to Train the Decision Tree Classifier :
#################################################
------ command ------
mpiexec -n [number_of_processes] ./run [job_configon_hdfs] [tree_config_local] ( [test_file_local] [test_meta_local] )

1. to run a treeServer program, you need 2 configuration files (job-config and tree-config), which should be passed as parameter argument

See the following sample config files and files specified inside those (i.e meta.csv) to understand how to use
(a) job-config (job.config): general system parameter, including the HDFS paths for training data and meta dat: [train_file_on_hdfs] [meta_file_on_hdfs]
(b) tree-config (tree2.config, tree.config): a list of training tasks, one for each decision tree or random forest

Note that we require a meta data where one row is for each data column. Refer to csv.h for the explanation comments on the metafile format

2. to calculate test accuracy, 2 optional configuration files (test.csv and test_meta.csv) can be added as argument to the program.

Here is a sample command to run:
mpiexec -n 3 ./run job.config data/tree2.config data/bank_test.csv data/bank_meta.csv

(note that test-metafile could be different from training-metafile, as the y-column may be missing)

3. In "config.h", you can switch on the following to print and save the output trees
bool SAVE_TREE = true
bool PRINT_TREE = false; // -> true

#################################################
Input File Checklist :
#################################################

* On HDFS (needed by all workers):
[job_config] (upload using hadoop fs -put)
[train_file_on_hdfs] (upload using our put cpp program)
[meta_file_on_hdfs] (upload using hadoop fs -put)

* Local to Master:
[tree_config]
[test_file_local]
[test_meta_local]

Note that testing is currently single-threaded, but can be easily parallelized with more coding
Models are output to the "job" folder at the local disk

#################################################
Toy Data :
#################################################
folders:
data - for classification
data_regression - for regression
