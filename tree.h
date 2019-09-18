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

#ifndef TREE_H
#define TREE_H

#include "splitter.h"
#include "tree_obj.h"
#include <climits>
#include <ctime>
#include "ioser.h"

void free_tree(TreeNode* root);
void print_tree(TreeNode* root, int tree_depth);

//////////////////// used by both Task_Slave, Task_Master /////////////////

template <class T>
class Best_Split {
public:
    bool isOrdinal = 3;
    double best_impurity; // best impurity found among the scanned columns
    T split_value; // "resp" need this information, master use this to create "TreeNode*"
    set<T> S;
    set<T> S1; // "resp" need this information, master use this to create "TreeNode*"

    size_t left_D; // so that master can decide plan type in advance
    size_t right_D; // so that master can decide plan type in advance
};

void delete_best_split(int best_column_index, void* best_split) {
#ifdef ASSERT
    assert(best_column_index >= 0);
#endif

    if(best_split == NULL) { //no "best" found yet to replace
        return;
    }

    Column* column = cserver.X.get_column(best_column_index);

    int data_type = column->data_type;

    if(data_type == ELEM_BOOL) {
        Best_Split<bool>* b_split = (Best_Split<bool>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_SHORT) {
        Best_Split<short>* b_split = (Best_Split<short>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_INT) {
        Best_Split<int>* b_split = (Best_Split<int>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_FLOAT) {
        Best_Split<float>* b_split = (Best_Split<float>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_DOUBLE) {
        Best_Split<double>* b_split = (Best_Split<double>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_CHAR) {
        Best_Split<char>* b_split = (Best_Split<char>*) best_split;
        delete b_split;
    } else if(data_type == ELEM_STRING) {
        Best_Split<string>* b_split = (Best_Split<string>*) best_split;
        delete b_split;
    } else {
        cout << "Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
        exit(-1);
    }

}

template <class T>
ibinstream& operator<<(ibinstream& m, const Best_Split<T>* best_split) {

    if(best_split->isOrdinal) {
        m << best_split->isOrdinal;
        m << best_split->best_impurity;
        m << best_split->split_value;
        m << best_split->left_D;
        m << best_split->right_D;
    } else {
        m << best_split->isOrdinal;
        m << best_split->best_impurity;
        m << best_split->S;
        m << best_split->S1;
        m << best_split->left_D;
        m << best_split->right_D;
    }

    return m;
}

//best_split must bed "new-ed" already
template <class T>
obinstream& operator>>(obinstream& m, Best_Split<T>* & best_split) {
#ifdef ASSERT
    assert(best_split != NULL);
#endif

    m >> best_split->isOrdinal;
    if(best_split->isOrdinal) {
        m >> best_split->best_impurity;
        m >> best_split->split_value;
        m >> best_split->left_D;
        m >> best_split->right_D;
    } else {
        m >> best_split->best_impurity;
        m >> best_split->S;
        m >> best_split->S1;
        m >> best_split->left_D;
        m >> best_split->right_D;
    }
    return m;
}

//////////////////// Task Slave ///////////////////////

class Candidate_Rows {
public:
    vector<size_t> indexes;
    bool stop_splitting = false; // best_left_rows_stop_condition
    string freq_y;
};

ibinstream& operator<<(ibinstream& m, const Candidate_Rows & candidate_rows) {
    m << candidate_rows.indexes;
    m << candidate_rows.stop_splitting;
    m << candidate_rows.freq_y;

    return m;
}

obinstream& operator>>(obinstream& m, Candidate_Rows & candidate_rows) {
    m >> candidate_rows.indexes;
    m >> candidate_rows.stop_splitting;
    m >> candidate_rows.freq_y;

    return m;
}

class Task_Slave {
public:
    int task_id;
    char task_type; // "sub_tree_task", "column_split_task", "data_serve_task"
    TreeConfig tree_config; // only for subtree_task and col_split_task
    Candidate_Rows candidate_rows;
};

class Data_Serve_Task : public Task_Slave {
public:
    int dest_id; // only or sub-tree task, to return the fetched column value to dest_id which is the slave to build subtree
    vector<int> column_indices; // local columns, copied from req_s2s_columns::column_indices

    Data_Serve_Task() {
        task_type = TASK_DATA_SERVE;
    }
};

class Task_Slave_Col_Split : public Task_Slave {
public:
    bool is_root_task = false;
    vector<int> column_indices; // local columns

    int best_column_idx = -1;

    void* best_split; // Best<T>* type, T depends on "best_column_index"

    int n_cols = 0; // "row_indices" serve for "n_cols",
    Candidate_Rows best_left_rows;
    Candidate_Rows best_right_rows;

    int parent_request = 0;

    Task_Slave_Col_Split() {
        task_type = TASK_COL_SPLIT;
    }
};

class Task_Slave_Subtree : public Task_Slave {
public:
    Matrix* matrix; // remote columns would be fetched and stored here
    TreeNode* root;
    size_t n_req = 0;
    size_t n_met = 0;
    vector<vector<int>> mac2colList; // mac2colList[mac_ID] = list of requested columns, populate from "plan" message
    int level = 0;

    vector<int> column_indices; // candidate column indices while building the sub_tree in key_slave, see usage in build_tree(..)

    Task_Slave_Subtree() {
        task_type = TASK_SUB_TREE;
        matrix = new Matrix;
    }

    ~Task_Slave_Subtree() {
        delete matrix;
    }
};


// !!! for debug only !!!
// [start, end) should've been sorted by Xi (recorded in SplitResult)
void print_internal_node_content(SplitResult* result, vector<size_t>::iterator &start,
                 vector<size_t>::iterator &end, Matrix &dataset, size_t y_index) {

    Column* Xi = dataset.get_column(result->column_idx);
    Column* Y = dataset.get_column(y_index);

    cout << "|D| = " << (end - start) << ", attr " << result->column_idx << " is the best" << endl;

    size_t count = 0;
    // printing index \t sample_index \t attribute_value  \t  y_value
    cout<< "No." << "\t" << "idx" << "\t" << "Xi" << "\t" << "y" << endl;
    for(auto it = start; it != end; it++){
        if(it == result->pos) {
            cout<< count << "*" << "\t" << *it << "*" << "\t" << Xi->string_value(*it) << "*" << "\t" << Y->string_value(*it) << "*" << endl;
        } else {
            cout<< count << "\t" << *it << "\t" << Xi->string_value(*it) << "\t" << Y->string_value(*it) << endl;
        }
        count++;
    }
}

// this function is to decide when to stop splitting given [start, end)
// param: start_iterator, end_iterator (exclusive), Column Y
// output: stop ?
// leaf_label is an output: majority label
template <class T>
bool stop_splitting_categorical(Matrix & dataset, vector<size_t>::iterator start, vector<size_t>::iterator end,
                    TreeConfig &treeConfig, T &leaf_label) { // checking whether all Y_values are same

    Column* Y = dataset.get_column(y_index);

#ifdef ASSERT
    assert(start != end);
#endif

    map<T, size_t> freq; // freq[y_label] = frequency of elements with y_label in [start, end)

    for(auto it = start; it != end; it++) {
        T value;
        Y->get(*it, &value);

        auto it2 = freq.find(value);
        if(it2 == freq.end()) {
            freq[value] = 1;
        } else {
            it2->second++;
        }
    }

    size_t max_freq = 0;

    for(auto it = freq.begin(); it != freq.end(); it++) {
        if(max_freq < it->second) {
            leaf_label = it->first;
            max_freq = it->second;
        }
    }

    return freq.size() == treeConfig.MIN_SAMPLE_LEAF; // stop splitting if there is only 1 item in Y[start, end)
}

// leaf_label is an output: average Y-value
template <class T>
bool stop_splitting_ordinal(Matrix & dataset, vector<size_t>::iterator start, vector<size_t>::iterator end,
                    TreeConfig & treeConfig, T & leaf_label) { // checking whether all Y_values are same
    Column* Y = dataset.get_column(y_index);

    double average = 0.0;

    for(auto it = start; it != end; it++) {
        T value;
        Y->get(*it, &value);
        average += value;
    }

    average = average / (end - start);
    leaf_label = average;

    return (end - start) == treeConfig.MIN_SAMPLE_LEAF;
}

template <class T>
TreeNode* create_leaf(vector<size_t>::iterator start, vector<size_t>::iterator end, T label)
{
    LeafNode<T>* node = new LeafNode<T>();
    node->start = start;
    node->end = end;
    node->size = end - start;
    node->column_index = -1; //-1 means leaf
    node->label = label;
    return node;
}

TreeNode* create_leaf_wrappper(Column* Y, vector<size_t>::iterator start,
                               vector<size_t>::iterator end, string label) {
    int y_type = Y->data_type;

    if(y_type == ELEM_BOOL) {
        bool node_label = boost::lexical_cast<bool>(label);
        return create_leaf<bool>(start, end, node_label);
    } else if(y_type == ELEM_SHORT) {
        short node_label = boost::lexical_cast<short>(label);
        return create_leaf<short>(start, end, node_label);
    } else if(y_type == ELEM_INT) {
        int node_label = boost::lexical_cast<int>(label);
        return create_leaf<int>(start, end, node_label);
    } else if(y_type == ELEM_FLOAT) {
        float node_label = boost::lexical_cast<float>(label);
        return create_leaf<float>(start, end, node_label);
    } else if(y_type == ELEM_DOUBLE) {
        double node_label = boost::lexical_cast<double>(label);
        return create_leaf<double>(start, end, node_label);
    } else if(y_type == ELEM_CHAR) {
        char node_label = boost::lexical_cast<char>(label);
        return create_leaf<char>(start, end, node_label);
    } else if(y_type == ELEM_STRING) {
        return create_leaf<string>(start, end, label);
    } else {
        cout << "ERROR: wrong type, type found = " << y_type << " File = " << __FILE__
             << ", Line = " << __LINE__ << endl;
        exit(-1);
    }
}

template <class T>
void set_categorical_node(TreeNode* &node, SplitResult *result, string node_label) {
    CategoricalSplitResult<T>* splitResult = (CategoricalSplitResult<T>*) result;
    CategoricalTreeNode<T>* treeNode = new CategoricalTreeNode<T>();

    node = treeNode;

    treeNode->label = node_label;
    treeNode->S = splitResult->S;
    treeNode->S1 = splitResult->S1;
}

TreeNode* build_tree(Task_Slave_Subtree* task, vector<size_t>::iterator start, vector<size_t>::iterator end, int tree_depth) {

    vector<int> & cols = task->column_indices;
    TreeConfig & treeConfig = task->tree_config;

    bool end_of_path = (tree_depth == treeConfig.MAX_TREE_DEPTH);

    Column* Y = task->matrix->get_column(y_index);
    int ytype = Y->data_type;
    bool stop;
    string node_label; // y value with highest frequency

    if(Y->is_ordinal) {
        if(ytype == ELEM_SHORT) {
            short label;
            stop = stop_splitting_ordinal<short>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<short>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_INT) {
            int label;
            stop = stop_splitting_ordinal<int>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<int>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_FLOAT) {
            float label;
            stop = stop_splitting_ordinal<float>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<float>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_DOUBLE) {
            double label;
            stop = stop_splitting_ordinal<double>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<double>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else {
            cout<<"File: " << __FILE__<<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
    } else {
        if(ytype == ELEM_BOOL) {
            bool label;
            stop = stop_splitting_categorical<bool>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<bool>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_SHORT) {
            short label;
            stop = stop_splitting_categorical<short>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<short>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_INT) {
            int label;
            stop = stop_splitting_categorical<int>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<int>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_FLOAT) {
            float label;
            stop = stop_splitting_categorical<float>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<float>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_DOUBLE) {
            double label;
            stop = stop_splitting_categorical<double>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<double>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_CHAR) {
            char label;
            stop = stop_splitting_categorical<char>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<char>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_STRING) {
            string label;
            stop = stop_splitting_categorical<string>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<string>(start, end, label);
                return node;
            }
            node_label = label;
        } else {
            cout<<"File: " << __FILE__<<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
    }

    //find the best split among all the features (Xi's)
    SplitResult* result = node_split(*task->matrix, start, end, cols, treeConfig); //sort [start, end) by the best Xi

    if(result->column_idx == -2) {
        TreeNode* leaf_node = create_leaf_wrappper(Y, start, end, node_label);
        return leaf_node;
    }

    // print_internal_node_content(result, start, end, train_set, y_index); //###### debug ######

    vector<size_t>::iterator middle = result->pos; // split row index left <, right >=

    //create current node and recurse on children
    Column* Xi = task->matrix->get_column(result->column_idx);

    TreeNode* curNode;

    if(Xi->is_ordinal) {
        //-----
        if(Xi->data_type == ELEM_BOOL) {
            OrdinalTreeNode<bool>* ordinalTreeNode = new OrdinalTreeNode<bool>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_SHORT) {
            OrdinalTreeNode<short>* ordinalTreeNode = new OrdinalTreeNode<short>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_INT) {
            OrdinalTreeNode<int>* ordinalTreeNode = new OrdinalTreeNode<int>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_FLOAT) {
            OrdinalTreeNode<float>* ordinalTreeNode = new OrdinalTreeNode<float>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_DOUBLE) {
            OrdinalTreeNode<double>* ordinalTreeNode = new OrdinalTreeNode<double>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_CHAR) {
            OrdinalTreeNode<char>* ordinalTreeNode = new OrdinalTreeNode<char>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_STRING) {
            OrdinalTreeNode<string>* ordinalTreeNode = new OrdinalTreeNode<string>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));
        } else {
            cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
        //------
        curNode->column_index = result->column_idx;
        curNode->start = start;
        curNode->end = end;
        curNode->size = end - start;
        curNode->left = build_tree(task, start, middle, tree_depth + 1); // < or less than attribute_val
        curNode->right = build_tree(task, middle, end, tree_depth + 1); // >= or greater than or equal attribute_val

        return curNode;
    } else {
        if(Xi->data_type == ELEM_BOOL) {
            set_categorical_node<bool>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_SHORT) {
            set_categorical_node<short>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_INT) {
            set_categorical_node<int>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_FLOAT) {
            set_categorical_node<float>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_DOUBLE) {
            set_categorical_node<double>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_CHAR) {
            set_categorical_node<char>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_STRING) {
            set_categorical_node<string>(curNode, result, node_label);
        } else {
            cout << "File = " << __FILE__ << ", Line = " << __LINE__ << "Wrong type and column index = "
                 << result->column_idx << endl;
            exit(-1);
        }

        curNode->column_index = result->column_idx;
        curNode->start = start;
        curNode->end = end;
        curNode->size = end - start;

        curNode->left = build_tree(task, start, middle, tree_depth + 1); // < or less than attribute_val
        curNode->right = build_tree(task, middle, end, tree_depth + 1); // >= or greater than or equal attribute_val
        return curNode;
    }
}

TreeNode* build_tree(Task_Slave_Subtree* task, vector<size_t>::iterator start, vector<size_t>::iterator end,
                     vector<int> & cols,
                     int tree_depth, vector<char> & progress, TreeConfig &treeConfig) { //with progress report
    
    bool end_of_path = (tree_depth == treeConfig.MAX_TREE_DEPTH);
    
    clock_t t = clock();
    //print progress
    if(tree_depth <= MAX_REPORT_DEPTH)
    {
        cout<<"[PROGRESS] training path: ";
        for(size_t i=0; i<progress.size(); i++) cout<<progress[i];
        cout<<" ......"<<endl;
    }

    Column* Y = cserver.X.get_column(y_index);
    int ytype = Y->data_type;
    bool stop;
    string node_label; // y value with highest frequency

    if(Y->is_ordinal) {
        if(ytype == ELEM_SHORT) {
            short label;
            stop = stop_splitting_ordinal<short>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<short>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_INT) {
            int label;
            stop = stop_splitting_ordinal<int>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<int>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_FLOAT) {
            float label;
            stop = stop_splitting_ordinal<float>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<float>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_DOUBLE) {
            double label;
            stop = stop_splitting_ordinal<double>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<double>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else {
            cout<<"File: " << __FILE__<<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
    } else {
        if(ytype == ELEM_BOOL) {
            bool label;
            stop = stop_splitting_categorical<bool>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<bool>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_SHORT) {
            short label;
            stop = stop_splitting_categorical<short>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<short>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_INT) {
            int label;
            stop = stop_splitting_categorical<int>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<int>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_FLOAT) {
            float label;
            stop = stop_splitting_categorical<float>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<float>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_DOUBLE) {
            double label;
            stop = stop_splitting_categorical<double>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<double>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_CHAR) {
            char label;
            stop = stop_splitting_categorical<char>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<char>(start, end, label);
                return node;
            }
            node_label = to_string(label);
        } else if(ytype == ELEM_STRING) {
            string label;
            stop = stop_splitting_categorical<string>(*(task->matrix), start, end, treeConfig, label);
            if(stop || end_of_path) {
                TreeNode* node = create_leaf<string>(start, end, label);
                return node;
            }
            node_label = label;
        } else {
            cout<<"File: " << __FILE__<<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
    }

    //find the best split among all the features (Xi's)
    SplitResult* result = node_split(*task->matrix, start, end, cols, treeConfig); //sort [start, end) by the best Xi

    // print_internal_node_content(result, start, end, train_set, y_index); //###### debug ######

    if(result->column_idx == -2) {
        TreeNode* leaf_node = create_leaf_wrappper(Y, start, end, node_label);
        return leaf_node;
    }

    vector<size_t>::iterator middle = result->pos; // split row index left <, right >=

    //create current node and recurse on children
    Column* Xi = task->matrix->get_column(result->column_idx);

    TreeNode* curNode;

    if(Xi->is_ordinal) {
        //-----
        if(Xi->data_type == ELEM_BOOL) {
            OrdinalTreeNode<bool>* ordinalTreeNode = new OrdinalTreeNode<bool>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_SHORT) {
            OrdinalTreeNode<short>* ordinalTreeNode = new OrdinalTreeNode<short>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_INT) {
            OrdinalTreeNode<int>* ordinalTreeNode = new OrdinalTreeNode<int>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_FLOAT) {
            OrdinalTreeNode<float>* ordinalTreeNode = new OrdinalTreeNode<float>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_DOUBLE) {
            OrdinalTreeNode<double>* ordinalTreeNode = new OrdinalTreeNode<double>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_CHAR) {
            OrdinalTreeNode<char>* ordinalTreeNode = new OrdinalTreeNode<char>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));

        } else if(Xi->data_type == ELEM_STRING) {
            OrdinalTreeNode<string>* ordinalTreeNode = new OrdinalTreeNode<string>();
            ordinalTreeNode->label = node_label;
            curNode = ordinalTreeNode;
            Xi->get(*result->pos, &(ordinalTreeNode->split_value));
        } else {
            cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [ERROR] we don't consider any other attribute type"<<endl;
            exit(-1);
        }
        //------
        curNode->column_index = result->column_idx;
        curNode->start = start;
        curNode->end = end;
        curNode->size = end - start;
        progress.push_back('0'); // 0 = left, can also use 'l'
        curNode->left = build_tree(task, start, middle, cols, tree_depth + 1, progress, treeConfig); // < or less than attribute_val
        progress.pop_back();
        progress.push_back('1'); // 1 = right, can also use 'r'
        curNode->right = build_tree(task, middle, end, cols, tree_depth + 1, progress, treeConfig); // >= or greater than or equal attribute_val
        progress.pop_back();
        if(tree_depth <= MAX_REPORT_DEPTH)
        {
            for(size_t i=0; i<progress.size(); i++) cout<<progress[i];
            cout << ": " << (float)(clock() - t)/CLOCKS_PER_SEC << " seconds" << endl;
        }
        return curNode;
    } else {
        if(Xi->data_type == ELEM_BOOL) {
            set_categorical_node<bool>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_SHORT) {
            set_categorical_node<short>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_INT) {
            set_categorical_node<int>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_FLOAT) {
            set_categorical_node<float>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_DOUBLE) {
            set_categorical_node<double>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_CHAR) {
            set_categorical_node<char>(curNode, result, node_label);
        } else if(Xi->data_type == ELEM_STRING) {
            set_categorical_node<string>(curNode, result, node_label);
        } else {
            cout << "File = " << __FILE__ << ", Line = " << __LINE__ << "Wrong type " << endl;
            exit(-1);
        }

        curNode->column_index = result->column_idx;
        curNode->start = start;
        curNode->end = end;
        curNode->size = end - start;

        progress.push_back('0'); // 0 = left, can also use 'l'
        curNode->left = build_tree(task, start, middle, cols, tree_depth + 1, progress, treeConfig); // < or less than attribute_val
        progress.pop_back();
        progress.push_back('1'); // 1 = right, can also use 'r'
        curNode->right = build_tree(task, middle, end, cols, tree_depth + 1, progress, treeConfig); // >= or greater than or equal attribute_val
        progress.pop_back();
        if(tree_depth <= MAX_REPORT_DEPTH)
        {
            for(size_t i=0; i<progress.size(); i++) cout<<progress[i];
            cout << ": " << (float)(clock() - t)/CLOCKS_PER_SEC << " seconds" << endl;
        }
        return curNode;
    }
}

template<class T>
void print_categorical_node(TreeNode* root) {
    CategoricalTreeNode<T>* treeNode = (CategoricalTreeNode<T>*) root;
    cout<< " in {";

    for(auto it = treeNode->S1.begin(); it != treeNode->S1.end(); it++) {
        if(it != treeNode->S1.begin()) cout<<", ";
        cout<< *it;
    }

    cout<<"}";
}

// printing the tree starting from root to each leaf
// recursive function called repeatedly until the end
// param  : TreeNode* root, Matrix& dataset and output_col_index
void print_tree(TreeNode* root, int tree_depth, string indent, bool last = true) {
    size_t sample_size = root->size;

    cout << indent << "+-";

    indent += last ? "   " : "|  ";

    if(root->column_index == -1) {
        cout<<"Leaf: |D| = "<<sample_size;
        Column* Y = cserver.X.get_column(y_index);

        if(Y->data_type == ELEM_BOOL) {
            LeafNode<bool>* leaf_node = (LeafNode<bool>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_SHORT) {
            LeafNode<short>* leaf_node = (LeafNode<short>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_INT) {
            LeafNode<int>* leaf_node = (LeafNode<int>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_FLOAT) {
            LeafNode<float>* leaf_node = (LeafNode<float>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leaf_node = (LeafNode<double>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_CHAR) {
            LeafNode<char>* leaf_node = (LeafNode<char>*) root;
            cout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_STRING) {
            LeafNode<string>* leaf_node = (LeafNode<string>*) root;
            cout<< ", label=" << leaf_node->label;
        } else {
            cout<< "[ERROR] Unknown Attribute(y) of type " << Y->data_type << endl;
            exit(-1);
        }
        cout<< endl;
    } else { // regular node other than leaf
        cout<<"Internal: |D| = "<<sample_size<<", X["<<(root->column_index)<<"]";

        Column* Xi = cserver.X.get_column(root->column_index);

        if(Xi->is_ordinal) {
            if(Xi->data_type == ELEM_BOOL) {
                OrdinalTreeNode<bool>* ordinalTreeNode = (OrdinalTreeNode<bool>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_SHORT) {
                OrdinalTreeNode<short>* ordinalTreeNode = (OrdinalTreeNode<short>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_INT) {
                OrdinalTreeNode<int>* ordinalTreeNode = (OrdinalTreeNode<int>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_FLOAT) {
                OrdinalTreeNode<float>* ordinalTreeNode = (OrdinalTreeNode<float>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_DOUBLE) {
                OrdinalTreeNode<double>* ordinalTreeNode = (OrdinalTreeNode<double>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_CHAR) {
                OrdinalTreeNode<char>* ordinalTreeNode = (OrdinalTreeNode<char>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_STRING) {
                OrdinalTreeNode<string>* ordinalTreeNode = (OrdinalTreeNode<string>*) root;
                cout<< "<" << ordinalTreeNode->split_value << ", freq_y = " << ordinalTreeNode->label;
            } else {
                cout << "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown Attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }
        } else {

            if(Xi->data_type == ELEM_BOOL) {
                print_categorical_node<bool>(root);
            } else if(Xi->data_type == ELEM_SHORT) {
                print_categorical_node<short>(root);
            } else if(Xi->data_type == ELEM_INT) {
                print_categorical_node<int>(root);
            } else if(Xi->data_type == ELEM_FLOAT) {
                print_categorical_node<float>(root);
            } else if(Xi->data_type == ELEM_DOUBLE) {
                print_categorical_node<double>(root);
            } else if(Xi->data_type == ELEM_CHAR) {
                print_categorical_node<char>(root);
            } else if(Xi->data_type == ELEM_STRING) {
                print_categorical_node<string>(root);
            } else {
                cout << "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown Attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }

        }

        cout<< endl;
        tree_depth++;
        print_tree(root->left, tree_depth, indent, false);
        print_tree(root->right, tree_depth, indent, true);
    }
}

//============== serializing tree ==============

template<class T>
void write_categorical_node(TreeNode* root, ifbinstream &fout) {
    CategoricalTreeNode<T>* treeNode = (CategoricalTreeNode<T>*) root;
    fout << treeNode->label; // ### write CategoricalTreeNode::label
    fout << treeNode->S; // ### write CategoricalTreeNode::label
    fout << treeNode->S1; // ### write CategoricalTreeNode::label
}

void save_subtree(TreeNode* root, ifbinstream & fout)
{
    //for "dataset" we actually just need its meta_data
    //we do not save [start, end)
    
    fout << root->column_index; // ### write TreeNode::column_index
    if(root->column_index == -1) {
        Column* Y = cserver.X.get_column(y_index);
        
        if(Y->data_type == ELEM_BOOL) {
            LeafNode<bool>* leaf_node = (LeafNode<bool>*) root;
            fout << leaf_node->label; // ### write LeafNode::label
        } else if(Y->data_type == ELEM_SHORT) {
            LeafNode<short>* leaf_node = (LeafNode<short>*) root;
            fout << leaf_node->label;
        } else if(Y->data_type == ELEM_INT) {
            LeafNode<int>* leaf_node = (LeafNode<int>*) root;
            fout << leaf_node->label;
        } else if(Y->data_type == ELEM_FLOAT) {
            LeafNode<float>* leaf_node = (LeafNode<float>*) root;
            fout << leaf_node->label;
        } else if(Y->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leaf_node = (LeafNode<double>*) root;
            fout << leaf_node->label;
        } else if(Y->data_type == ELEM_CHAR) {
            LeafNode<char>* leaf_node = (LeafNode<char>*) root;
            fout << leaf_node->label;
        } else if(Y->data_type == ELEM_STRING) {
            LeafNode<string>* leaf_node = (LeafNode<string>*) root;
            fout << leaf_node->label;
        } else {
            cout<< "[ERROR] Unknown attribute(y) of type " << Y->data_type << endl;
            exit(-1);
        }
    } else { // regular node other than leaf
        Column* Xi = cserver.X.get_column(root->column_index);
        
        if(Xi->is_ordinal) {
            if(Xi->data_type == ELEM_BOOL) {
                OrdinalTreeNode<bool>* ordinalTreeNode = (OrdinalTreeNode<bool>*) root;
                fout << ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_SHORT) {
                OrdinalTreeNode<short>* ordinalTreeNode = (OrdinalTreeNode<short>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_INT) {
                OrdinalTreeNode<int>* ordinalTreeNode = (OrdinalTreeNode<int>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_FLOAT) {
                OrdinalTreeNode<float>* ordinalTreeNode = (OrdinalTreeNode<float>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_DOUBLE) {
                OrdinalTreeNode<double>* ordinalTreeNode = (OrdinalTreeNode<double>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_CHAR) {
                OrdinalTreeNode<char>* ordinalTreeNode = (OrdinalTreeNode<char>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else if(Xi->data_type == ELEM_STRING) {
                OrdinalTreeNode<string>* ordinalTreeNode = (OrdinalTreeNode<string>*) root;
                fout << ordinalTreeNode->split_value;
                fout << ordinalTreeNode->label;
            } else {
                cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }
        } else {

            if(Xi->data_type == ELEM_BOOL) {
                write_categorical_node<bool>(root, fout);
            } else if(Xi->data_type == ELEM_SHORT) {
                write_categorical_node<short>(root, fout);
            } else if(Xi->data_type == ELEM_INT) {
                write_categorical_node<int>(root, fout);
            } else if(Xi->data_type == ELEM_FLOAT) {
                write_categorical_node<float>(root, fout);
            } else if(Xi->data_type == ELEM_DOUBLE) {
                write_categorical_node<double>(root, fout);
            } else if(Xi->data_type == ELEM_CHAR) {
                write_categorical_node<char>(root, fout);
            } else if(Xi->data_type == ELEM_STRING) {
                write_categorical_node<string>(root, fout);
            } else {
                cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }
        }

        save_subtree(root->left, fout);
        save_subtree(root->right, fout);
    }
}

void save_forest(const char* path, vector<TreeNode*> forest)
{
    ifbinstream fout(path);
    fout << forest.size();

    for(size_t index = 0; index < forest.size(); index++) {
        save_subtree(forest[index], fout);
    }

    fout.close();
}

template <class T>
void load_categorical_node(TreeNode* &result, ofbinstream & fin) {
    CategoricalTreeNode<T>* treeNode = new CategoricalTreeNode<T>;
    fin >> treeNode->label; // ### write CategoricalTreeNode::label
    fin >> treeNode->S; // ### write CategoricalTreeNode::label
    fin >> treeNode->S1; // ### write CategoricalTreeNode::label
    result = treeNode;
}

TreeNode* load_subtree(Matrix &dataset, int y_index, ofbinstream & fin)
{
    //for "dataset" we actually just need its meta_data
    //we do not load [start, end)
    int col_idx;
    fin >> col_idx;

    TreeNode* result;

    if(col_idx == -1) {
        Column* Y = dataset.get_column(y_index);

        if(Y->data_type == ELEM_BOOL) {
            LeafNode<bool>* leaf_node = new LeafNode<bool>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_SHORT) {
            LeafNode<short>* leaf_node = new LeafNode<short>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_INT) {
            LeafNode<int>* leaf_node = new LeafNode<int>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_FLOAT) {
            LeafNode<float>* leaf_node = new LeafNode<float>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leaf_node = new LeafNode<double>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_CHAR) {
            LeafNode<char>* leaf_node = new LeafNode<char>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else if(Y->data_type == ELEM_STRING) {
            LeafNode<string>* leaf_node = new LeafNode<string>;
            fin >> leaf_node->label; // ### write LeafNode::label
            leaf_node->column_index = col_idx; // ### write TreeNode::column_index
            result = leaf_node;
        } else {
            cout<< "[ERROR] Unknown attribute(y) of type " << Y->data_type << endl;
            exit(-1);
        }
    } else { // regular node other than leaf
        Column* Xi = dataset.get_column(col_idx);

        if(Xi->is_ordinal) {
            if(Xi->data_type == ELEM_BOOL) {
                OrdinalTreeNode<bool>* ordinalTreeNode = new OrdinalTreeNode<bool>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_SHORT) {
                OrdinalTreeNode<short>* ordinalTreeNode = new OrdinalTreeNode<short>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_INT) {
                OrdinalTreeNode<int>* ordinalTreeNode = new OrdinalTreeNode<int>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_FLOAT) {
                OrdinalTreeNode<float>* ordinalTreeNode = new OrdinalTreeNode<float>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_DOUBLE) {
                OrdinalTreeNode<double>* ordinalTreeNode = new OrdinalTreeNode<double>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_CHAR) {
                OrdinalTreeNode<char>* ordinalTreeNode = new OrdinalTreeNode<char>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else if(Xi->data_type == ELEM_STRING) {
                OrdinalTreeNode<string>* ordinalTreeNode = new OrdinalTreeNode<string>;
                fin >> ordinalTreeNode->split_value; // ### write OrdinalTreeNode::split_value
                fin >> ordinalTreeNode->label;
                ordinalTreeNode->column_index = col_idx; // ### write TreeNode::column_index
                result = ordinalTreeNode;
            } else {
                cout<< "[ERROR] Unknown attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }
        } else {
            if(Xi->data_type == ELEM_BOOL) {
                load_categorical_node<bool>(result, fin);
            } else if(Xi->data_type == ELEM_SHORT) {
                load_categorical_node<short>(result, fin);
            } else if(Xi->data_type == ELEM_INT) {
                load_categorical_node<int>(result, fin);
            } else if(Xi->data_type == ELEM_FLOAT) {
                load_categorical_node<float>(result, fin);
            } else if(Xi->data_type == ELEM_DOUBLE) {
                load_categorical_node<double>(result, fin);
            } else if(Xi->data_type == ELEM_CHAR) {
                load_categorical_node<char>(result, fin);
            } else if(Xi->data_type == ELEM_STRING) {
                load_categorical_node<string>(result, fin);
            } else {
                cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }

            result->column_index = col_idx; // ### write TreeNode::column_index
        }
        //recurse on subtrees
        result->left = load_subtree(dataset, y_index, fin);
        result->right = load_subtree(dataset, y_index, fin);
    }

    return result;
}

vector<TreeNode*> load_forest(Matrix & dataset, int y_index, const char* path)
{
    vector<TreeNode*> forest;
    ofbinstream fin(path);

    size_t n_tree;
    fin >> n_tree;

    cout << "# of trees found: " << n_tree << endl;

    for(size_t index = 0; index < n_tree; index++) {
        forest.push_back(load_subtree(dataset, y_index, fin));
    }

    fin.close();

    return forest;
}

//============== visualizing tree ==============

template <class T>
void print_categorical_node_json(TreeNode* root, ofstream &fout) {
    CategoricalTreeNode<T>* treeNode = (CategoricalTreeNode<T>*) root;
    fout<< "@{";

    for(auto it = treeNode->S1.begin(); it != treeNode->S1.end(); it++) {
        if(it != treeNode->S1.begin()) fout<<", ";
        string s1_str = boost::lexical_cast<string>(*it);
        fout<< filterQuote(s1_str);
    }
    fout<<"}\",\n";
}

void print_json(TreeNode* root,  Matrix &dataset, int y_index, ofstream & fout, string indent, bool add_comma) {
    size_t sample_size = root->size;
    fout<<indent<<"{\n"<<indent<<"\"name\": \"|D| = "<<sample_size;
    
    if(root->column_index == -1) {
        Column* Y = dataset.get_column(y_index);
        
        if(Y->data_type == ELEM_BOOL) {
            LeafNode<bool>* leaf_node = (LeafNode<bool>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_SHORT) {
            LeafNode<short>* leaf_node = (LeafNode<short>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_INT) {
            LeafNode<int>* leaf_node = (LeafNode<int>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_FLOAT) {
            LeafNode<float>* leaf_node = (LeafNode<float>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leaf_node = (LeafNode<double>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_CHAR) {
            LeafNode<char>* leaf_node = (LeafNode<char>*) root;
            fout<< ", label=" << leaf_node->label;
        } else if(Y->data_type == ELEM_STRING) {
            LeafNode<string>* leaf_node = (LeafNode<string>*) root;
            fout<< ", label=" << filterQuote(leaf_node->label);
        } else {
            cout<< "[ERROR] Unknown Attribute(y) of type " << Y->data_type << endl;
            exit(-1);
        }
        fout<<"\"\n";
    } else { // regular node other than leaf
        fout<<", X["<<(root->column_index)<<"]";
        Column* Xi = dataset.get_column(root->column_index);
        
        if(Xi->is_ordinal) {
            if(Xi->data_type == ELEM_BOOL) {
                OrdinalTreeNode<bool>* ordinalTreeNode = (OrdinalTreeNode<bool>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_SHORT) {
                OrdinalTreeNode<short>* ordinalTreeNode = (OrdinalTreeNode<short>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_INT) {
                OrdinalTreeNode<int>* ordinalTreeNode = (OrdinalTreeNode<int>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_FLOAT) {
                OrdinalTreeNode<float>* ordinalTreeNode = (OrdinalTreeNode<float>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_DOUBLE) {
                OrdinalTreeNode<double>* ordinalTreeNode = (OrdinalTreeNode<double>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_CHAR) {
                OrdinalTreeNode<char>* ordinalTreeNode = (OrdinalTreeNode<char>*) root;
                fout<< "<" << ordinalTreeNode->split_value;
            } else if(Xi->data_type == ELEM_STRING) {
                OrdinalTreeNode<string>* ordinalTreeNode = (OrdinalTreeNode<string>*) root;
                fout<< "<" << filterQuote(ordinalTreeNode->split_value);
            } else {
                cout<< "[ERROR] Unknown Attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }
            fout<<"\",\n";
        } else {

            if(Xi->data_type == ELEM_BOOL) {
                print_categorical_node_json<bool>(root, fout);
            } else if(Xi->data_type == ELEM_SHORT) {
                print_categorical_node_json<short>(root, fout);
            } else if(Xi->data_type == ELEM_INT) {
                print_categorical_node_json<int>(root, fout);
            } else if(Xi->data_type == ELEM_FLOAT) {
                print_categorical_node_json<float>(root, fout);
            } else if(Xi->data_type == ELEM_DOUBLE) {
                print_categorical_node_json<double>(root, fout);
            } else if(Xi->data_type == ELEM_CHAR) {
                print_categorical_node_json<char>(root, fout);
            } else if(Xi->data_type == ELEM_STRING) {
                print_categorical_node_json<string>(root, fout);
            } else {
                cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << "[ERROR] Unknown Attribute(Xi) of type " << Xi->data_type << endl;
                exit(-1);
            }

        }
        
        fout<<indent<<"\"children\": [\n";
        print_json(root->left, dataset, y_index, fout, indent+"\t", true);
        print_json(root->right, dataset, y_index, fout, indent+"\t", false);
        fout<<indent<<"]\n";
    }
    if(add_comma) fout<<indent<<"},\n";
    else fout<<indent<<"}\n";
}

void export_json(TreeNode* root,  Matrix &dataset, int y_index, const char* fname)
{
    ofstream fout(fname);
    print_json(root, dataset, y_index, fout, "", false);
    fout.close();
}

//if you directly open index.html, you see blank page

//you have to start a web server with the following command
//python -m http.server

//then go to your browser to see the tree

#endif
