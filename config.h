//########################################################################
//## Copyright 2019 Da Yan http://www.cs.uab.edu/yanda
//##
//## Licensed under the Apache License, Version 2.0 (the "License");
//## you may not use this file except in compliance with the License.
//## You may obtain a copy of the License at
//##
//## //http://www.apache.org/licenses/LICENSE-2.0
//##
//## Unless required by applicable law or agreed to in writing, software
//## distributed under the License is distributed on an "AS IS" BASIS,
//## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//## See the License for the specific language governing permissions and
//## limitations under the License.
//########################################################################

//########################################################################
//## Contributors
//## CHOWDHURY, Md Mashiur Rahman    (Mashiur)
//## YAN, Da    (Daniel)
//########################################################################

#ifndef CONFIG_H
#define CONFIG_H

#include <iostream>
#include "matrix.h"
#include "tree_obj.h"
#include "ydhdfs.h"

using namespace std;

const string DECISION_TREE = "DT";
const string RANDOM_FOREST = "RF";
const string EXTRA_TREES = "ET";

//the configurations for a tree to be built
//(a job can build many trees with different hyperparameters)

void print_job_config() {
    cout<< "================= CONFIGURATION =================" << endl;
    
    cout<< "BRUTE_FORCE_MAX_ITEM = " << BRUTE_FORCE_MAX_ITEM << endl;
    cout<< "NO_WARNING = " << NO_WARNING << endl;
    cout<< "MAX_REPORT_DEPTH = " << MAX_REPORT_DEPTH << endl;
    cout<< "num_compers = " << num_compers << endl;
    cout<< "|subtree_D| = " << subtree_D << endl;
    cout<< "Y column index = " << y_index << endl;
    cout<< "Train file = " << train_file << endl;
    cout<< "Meta file = " << meta_file << endl;
    cout<< "Job output directory = " << job_dir << endl;
    
    cout<< endl;

#ifdef DEBUG_LOG
    fout<< "================= CONFIGURATION =================" << endl;

    fout<< "BRUTE_FORCE_MAX_ITEM = " << BRUTE_FORCE_MAX_ITEM << endl;
    fout<< "NO_WARNING = " << NO_WARNING << endl;
    fout<< "MAX_REPORT_DEPTH = " << MAX_REPORT_DEPTH << endl;
    fout<< "num_compers = " << num_compers << endl;
    fout<< "|subtree_D| = " << subtree_D << endl;
    fout<< "Y column index = " << y_index << endl;
    fout<< "Train file = " << train_file << endl;
    fout<< "Meta file = " << meta_file << endl;
    fout<< "Job output directory = " << job_dir << endl;

    fout<< endl;
#endif
}

class TreeConfig {
public:
    string type; // DT, RF
    double column_sample = 1; // percentage of columns to consider at each tree, 100% for a single tree
    bool sample_col_each_node = false;// sample column for each node //todo:: add serialization !!!!!!!!!!!!

    int IMPURITY_FUNC = -1;
    int MAX_TREE_DEPTH = INT_MAX;
    int MIN_SAMPLE_LEAF = 1;
    int num_trees = 1;
    int n_columns; // |sample_column|

    // depends on col_sample_percentage
    // computed by MASTER
    vector<vector<int>> column_distribution;
    vector<TreeNode*> rootList;
    Matrix test_set;
};

ibinstream & operator<<(ibinstream & m, const TreeConfig & treeConfig) {
	m << treeConfig.type;
	m << treeConfig.column_sample;
	m << treeConfig.sample_col_each_node;
    m << treeConfig.IMPURITY_FUNC;
    m << treeConfig.MAX_TREE_DEPTH;
    m << treeConfig.MIN_SAMPLE_LEAF;
    m << treeConfig.n_columns;

    return m;
}

obinstream & operator>>(obinstream & m, TreeConfig & treeConfig) {
	m >> treeConfig.type;
	m >> treeConfig.column_sample;
	m >> treeConfig.sample_col_each_node;
    m >> treeConfig.IMPURITY_FUNC;
    m >> treeConfig.MAX_TREE_DEPTH;
    m >> treeConfig.MIN_SAMPLE_LEAF;
    m >> treeConfig.n_columns;

    return m;
}

void set_config(TreeConfig src, TreeConfig & dest, int tree_index) { // for serialization
    dest.type = src.type;
    dest.column_sample = src.column_sample;
    dest.sample_col_each_node = src.sample_col_each_node;
    dest.num_trees = src.num_trees;

    // serialize/deserialize
    dest.IMPURITY_FUNC = src.IMPURITY_FUNC;
    dest.MAX_TREE_DEPTH = src.MAX_TREE_DEPTH;
    dest.MIN_SAMPLE_LEAF = src.MIN_SAMPLE_LEAF;
    dest.n_columns = src.column_distribution[tree_index].size();
}

void load_job_config(const char* job_file) {
	hdfsFS fs = getHdfsFS();
	hdfsFile in = getRHandle(job_file, fs);
	LineReader reader(fs, in);
    string line, value;

    //# get BRUTE_FORCE_MAX_ITEM
    reader.readLine();
    line = reader.getLine();
    while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream max_item(line);
    getline(max_item, value, '#');
    BRUTE_FORCE_MAX_ITEM = atoi(value.c_str());

    //# get num_compers
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream n_threads(line);
    getline(n_threads, value, '#');
    num_compers = atoi(value.c_str());

    //# get num_samples
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream n_samples(line);
    getline(n_samples, value, '#');
    _n_samples = atoi(value.c_str());

    //# get subtree_D
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
	istringstream s_d(line);
	getline(s_d, value, '#');
	subtree_D = atoi(value.c_str());

    //# get NO_WARNING
	reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream wrn(line);
    getline(wrn, value, '#');
    NO_WARNING = atoi(value.c_str());

    //# get train data file path
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream train_f(line);
    getline(train_f, value, '#');
    istringstream(value) >> train_file; //trim

    //# get meta file path
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream meta_f(line);
    getline(meta_f, value, '#');
    istringstream(value) >> meta_file; //trim

    //# get y_index
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream y_idx(line);
    getline(y_idx, value, '#');
    y_index = atoi(value.c_str());

    //# get job dir path
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream j_dir(line);
    getline(j_dir, value, '#');
    istringstream(value) >> job_dir; //trim

    //# get SAVE_TREE
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream save_tree(line);
    getline(save_tree, value, '#');
    SAVE_TREE = atoi(value.c_str());

    //# get BFS_PRIORITY_THRESHOLD
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream bfs_th(line);
    getline(bfs_th, value, '#');
    BFS_PRIORITY_THRESHOLD = atoi(value.c_str());

    //# get ACTIVE_TREE_THRESHOLD
    reader.readLine();
	line = reader.getLine();
	while(line.length() == 0 || line[0] == '#') line = reader.getLine(); // skip blank lines or lines starting with #
    istringstream at_th(line);
    getline(at_th, value, '#');
    ACTIVE_TREES_THRESHOLD = atoi(value.c_str());

    hdfsCloseFile(fs, in);
    hdfsDisconnect(fs);
}

void load_local_job_config(const char* job_file) {
    ifstream fin(job_file);
    string line, value;

    //# get BRUTE_FORCE_MAX_ITEM
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream max_item(line);
    getline(max_item, value, '#');
    BRUTE_FORCE_MAX_ITEM = atoi(value.c_str());

    //# get num_compers
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream n_threads(line);
    getline(n_threads, value, '#');
    num_compers = atoi(value.c_str());

    //# get num_samples
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream n_samples(line);
    getline(n_samples, value, '#');
    _n_samples = atoi(value.c_str());

    //# get subtree_D
	getline(fin, line);
	while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
	istringstream s_d(line);
	getline(s_d, value, '#');
	subtree_D = atoi(value.c_str());

    //# get NO_WARNING
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream wrn(line);
    getline(wrn, value, '#');
    NO_WARNING = atoi(value.c_str());

    //# get train data file path
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream train_f(line);
    getline(train_f, value, '#');
    istringstream(value) >> train_file; //trim

    //# get meta file path
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream meta_f(line);
    getline(meta_f, value, '#');
    istringstream(value) >> meta_file; //trim

    //# get y_index
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream y_idx(line);
    getline(y_idx, value, '#');
    y_index = atoi(value.c_str());

    //# get job dir path
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream j_dir(line);
    getline(j_dir, value, '#');
    istringstream(value) >> job_dir; //trim

    //# get SAVE_TREE
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream save_tree(line);
    getline(save_tree, value, '#');
    SAVE_TREE = atoi(value.c_str());

    //# get BFS_PRIORITY_THRESHOLD
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream bfs_th(line);
    getline(bfs_th, value, '#');
    BFS_PRIORITY_THRESHOLD = atoi(value.c_str());

    //# get ACTIVE_TREE_THRESHOLD
    getline(fin, line);
    while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #
    istringstream at_th(line);
    getline(at_th, value, '#');
    ACTIVE_TREES_THRESHOLD = atoi(value.c_str());

    fin.close();
}

void load_config(const char* tree_config_fname, vector<TreeConfig> &configList) {
    ifstream fin(tree_config_fname);
    string line, value;

    while(getline(fin, line)) {
        //# get MAX_TREE_DEPTH
        while(line.length() == 0 || line[0] == '#'){
            if(!getline(fin, line)) return; // skip blank lines or lines starting with #
        }

        TreeConfig t_config;

        //  Tree type RF / DT
        istringstream tree_type(line);
        getline(tree_type, value, '#');
        t_config.type = value;
        boost::trim(t_config.type);

#ifdef DEBUG_LOG
        fout << "type found = " << t_config.type << endl;
#endif


        if(t_config.type == RANDOM_FOREST || t_config.type == EXTRA_TREES) {
            //# get sample_percentage
            getline(fin, line);
            while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

            istringstream column_sample(line);
            getline(column_sample, value, '#');
            t_config.column_sample = atof(value.c_str());

#ifdef DEBUG_LOG
            fout << "column sample found = " << t_config.column_sample << endl;
#endif

            //# get sample_column_each_node
            getline(fin, line);
            while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

            istringstream s_each_node(line);
            getline(s_each_node, value, '#');
            t_config.sample_col_each_node = atoi(value.c_str());

#ifdef DEBUG_LOG
            fout << "sample column each node found = " << t_config.sample_col_each_node << endl;
#endif

            //# get num_trees
            getline(fin, line);
            while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

            istringstream num_trees(line);
            getline(num_trees, value, '#');
            t_config.num_trees = atoi(value.c_str());

#ifdef DEBUG_LOG
            fout << "num tress found = " << t_config.num_trees << endl;
#endif
        }

        //# get MAX_TREE_DEPTH
        getline(fin, line);
        while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

        istringstream depth(line);
        getline(depth, value, '#');
        t_config.MAX_TREE_DEPTH = atoi(value.c_str());

        //# get IMPURITY_FUNC
        getline(fin, line);
        while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

        istringstream func(line);
        getline(func, value, '#');

        t_config.IMPURITY_FUNC = atoi(value.c_str());

        //# get MIN_SAMPLE_LEAF
        getline(fin, line);
        while(line.length() == 0 || line[0] == '#') getline(fin, line); // skip blank lines or lines starting with #

        istringstream msl(line);
        getline(msl, value, '#');
        t_config.MIN_SAMPLE_LEAF = atoi(value.c_str());

        if(t_config.MIN_SAMPLE_LEAF < 0) {
            cout<<__LINE__ << " [ERROR] : MIN_SAMPLE_LEAF expected as configuration parameter" << endl;
            exit(-1);
        }

        configList.push_back(t_config);
    }

    fin.close();
}

#endif
