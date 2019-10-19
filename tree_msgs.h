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

#ifndef TREE_MSGS_H
#define TREE_MSGS_H

#include "tree.h"
#include "communication.h"
#include "mpi.h"
#include "master/utils.h"

void init_worker(int * argc, char*** argv)
{
    int provided;
    MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
    if(provided != MPI_THREAD_MULTIPLE)
    {
        printf("MPI do not Support Multiple thread\n");
        exit(0);
    }
    MPI_Comm_size(MPI_COMM_WORLD, &_num_workers);
    MPI_Comm_rank(MPI_COMM_WORLD, &_my_rank);
}

void worker_finalize()
{
    MPI_Finalize();
}

/////////////// Tree Node /////////////////////////////

void serialize_tree_node(ibinstream& m, const TreeNode* tree_node);
void deserialize_tree_node(obinstream &m, TreeNode* & tree_node);

template <>
ibinstream& operator<<(ibinstream &m, const TreeNode* tree_node) {
    serialize_tree_node(m, tree_node);
    return m;
}

template <>
obinstream& operator>>(obinstream &m, TreeNode* & tree_node) {
    deserialize_tree_node(m, tree_node);
    return m;
}

///////////////////// Ordinal Tree Node ///////////////////////////

template <class T>
ibinstream& operator<<(ibinstream &m, const OrdinalTreeNode<T>* tree_node) {
    m << tree_node->column_index;
    m << tree_node->size;
    m << tree_node->split_value;
    m << tree_node->label;

    serialize_tree_node(m, tree_node->left);
    serialize_tree_node(m, tree_node->right);

    return m;
}

template <class T>
obinstream& operator>>(obinstream &m, OrdinalTreeNode<T>* & tree_node) {
    //m >> tree_node->column_index;
    m >> tree_node->size;
    m >> tree_node->split_value;
    m >> tree_node->label;

    deserialize_tree_node(m, tree_node->left);
    deserialize_tree_node(m, tree_node->right);

    return m;
}

////////////////// Categorical Tree Node /////////////////////////////

template <class T>
ibinstream& operator<<(ibinstream &m, const CategoricalTreeNode<T>* tree_node) {
    m << tree_node->column_index;
    m << tree_node->size;
    m << tree_node->label;
    m << tree_node->S;
    m << tree_node->S1;

    serialize_tree_node(m, tree_node->left);
    serialize_tree_node(m, tree_node->right);

    return m;
}

template <class T>
obinstream& operator>>(obinstream &m, CategoricalTreeNode<T>* & tree_node) {
    //m >> tree_node->column_index;
    m >> tree_node->size;
    m >> tree_node->label;
    m >> tree_node->S;
    m >> tree_node->S1;

    deserialize_tree_node(m, tree_node->left);
    deserialize_tree_node(m, tree_node->right);

    return m;
}

//////////////// Leaf Tree Node ///////////////////////////

template <class T>
ibinstream& operator<<(ibinstream &m, const LeafNode<T>* tree_node) {
    m << tree_node->column_index;
    m << tree_node->size;
    m << tree_node->label;

    return m;
}

template <class T>
obinstream& operator>>(obinstream &m, LeafNode<T>* & tree_node) {
    //m >> tree_node->column_index;
    m >> tree_node->size;
    m >> tree_node->label;

    return m;
}

void serialize_tree_node(ibinstream& m, const TreeNode* tree_node) {
    Matrix & meta = cserver.X;

    if(tree_node->column_index == -1) {
        Column* column = meta.get_column(y_index);

        if(column->data_type == ELEM_BOOL) {
            LeafNode<bool>* leafNode = (LeafNode<bool>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_SHORT) {
            LeafNode<short>* leafNode = (LeafNode<short>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_INT) {
            LeafNode<int>* leafNode = (LeafNode<int>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_FLOAT) {
            LeafNode<float>* leafNode = (LeafNode<float>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leafNode = (LeafNode<double>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_CHAR) {
            LeafNode<char>* leafNode = (LeafNode<char>*) tree_node;
            m << leafNode;
        } else if(column->data_type == ELEM_STRING) {
            LeafNode<string>* leafNode = (LeafNode<string>*) tree_node;
            m << leafNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }

        return;
    }

    Column* column = meta.get_column(tree_node->column_index);

    if(column->is_ordinal) {
        if(column->data_type == ELEM_BOOL) {
            OrdinalTreeNode<bool>* ordinalTreeNode = (OrdinalTreeNode<bool>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_SHORT) {
            OrdinalTreeNode<short>* ordinalTreeNode = (OrdinalTreeNode<short>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_INT) {
            OrdinalTreeNode<int>* ordinalTreeNode = (OrdinalTreeNode<int>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_FLOAT) {
            OrdinalTreeNode<float>* ordinalTreeNode = (OrdinalTreeNode<float>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_DOUBLE) {
            OrdinalTreeNode<double>* ordinalTreeNode = (OrdinalTreeNode<double>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_CHAR) {
            OrdinalTreeNode<char>* ordinalTreeNode = (OrdinalTreeNode<char>*) tree_node;
            m << ordinalTreeNode;
        } else if(column->data_type == ELEM_STRING) {
            OrdinalTreeNode<string>* ordinalTreeNode = (OrdinalTreeNode<string>*) tree_node;
            m << ordinalTreeNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    } else {
        if(column->data_type == ELEM_BOOL) {
            CategoricalTreeNode<bool>* categoricalTreeNode = (CategoricalTreeNode<bool>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_SHORT) {
            CategoricalTreeNode<short>* categoricalTreeNode = (CategoricalTreeNode<short>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_INT) {
            CategoricalTreeNode<int>* categoricalTreeNode = (CategoricalTreeNode<int>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_FLOAT) {
            CategoricalTreeNode<float>* categoricalTreeNode = (CategoricalTreeNode<float>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_DOUBLE) {
            CategoricalTreeNode<double>* categoricalTreeNode = (CategoricalTreeNode<double>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_CHAR) {
            CategoricalTreeNode<char>* categoricalTreeNode = (CategoricalTreeNode<char>*) tree_node;
            m << categoricalTreeNode;
        } else if(column->data_type == ELEM_STRING) {
            CategoricalTreeNode<string>* categoricalTreeNode = (CategoricalTreeNode<string>*) tree_node;
            m << categoricalTreeNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    }
}

void deserialize_tree_node(obinstream &m, TreeNode* & tree_node) {
    Matrix & meta = cserver.X;

    int column_index;
    m >> column_index;

    if(column_index == -1) {
        Column* y_column = meta.get_column(y_index);

        if(y_column->data_type == ELEM_BOOL) {
            LeafNode<bool>* leafNode = new LeafNode<bool>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_SHORT) {
            LeafNode<short>* leafNode = new LeafNode<short>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_INT) {
            LeafNode<int>* leafNode = new LeafNode<int>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_FLOAT) {
            LeafNode<float>* leafNode = new LeafNode<float>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_DOUBLE) {
            LeafNode<double>* leafNode = new LeafNode<double>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_CHAR) {
            LeafNode<char>* leafNode = new LeafNode<char>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else if(y_column->data_type == ELEM_STRING) {
            LeafNode<string>* leafNode = new LeafNode<string>();
            leafNode->column_index = column_index;
            m >> leafNode;
            tree_node = leafNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }

        return;
    }

    Column* column = meta.get_column(column_index);

#ifdef ASSERT
    assert(column != nullptr);
#endif

    if(column->is_ordinal) {
        if(column->data_type == ELEM_BOOL) {
            OrdinalTreeNode<bool>* ordinalTreeNode = new OrdinalTreeNode<bool>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_SHORT) {
            OrdinalTreeNode<short>* ordinalTreeNode = new OrdinalTreeNode<short>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_INT) {
            OrdinalTreeNode<int>* ordinalTreeNode = new OrdinalTreeNode<int>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_FLOAT) {
            OrdinalTreeNode<float>* ordinalTreeNode = new OrdinalTreeNode<float>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_DOUBLE) {
            OrdinalTreeNode<double>* ordinalTreeNode = new OrdinalTreeNode<double>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_CHAR) {
            OrdinalTreeNode<char>* ordinalTreeNode = new OrdinalTreeNode<char>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else if(column->data_type == ELEM_STRING) {
            OrdinalTreeNode<string>* ordinalTreeNode = new OrdinalTreeNode<string>();
            ordinalTreeNode->column_index = column_index;
            m >> ordinalTreeNode;
            tree_node = ordinalTreeNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    } else {
        if(column->data_type == ELEM_BOOL) {
            CategoricalTreeNode<bool>* categoricalTreeNode = new CategoricalTreeNode<bool>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_SHORT) {
            CategoricalTreeNode<short>* categoricalTreeNode = new CategoricalTreeNode<short>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_INT) {
            CategoricalTreeNode<int>* categoricalTreeNode = new CategoricalTreeNode<int>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_FLOAT) {
            CategoricalTreeNode<float>* categoricalTreeNode = new CategoricalTreeNode<float>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_DOUBLE) {
            CategoricalTreeNode<double>* categoricalTreeNode = new CategoricalTreeNode<double>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_CHAR) {
            CategoricalTreeNode<char>* categoricalTreeNode = new CategoricalTreeNode<char>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else if(column->data_type == ELEM_STRING) {
            CategoricalTreeNode<string>* categoricalTreeNode = new CategoricalTreeNode<string>();
            categoricalTreeNode->column_index = column_index;
            m >> categoricalTreeNode;
            tree_node = categoricalTreeNode;
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    }
}

//////////////// slave's response to master //////////////////
class resp2master {
public:
    char message_type; // SUB_TREE_RESP or COL_SPLIT_RESP
    int task_id;

    int best_column_index = -1; // -1 means the current node is a leaf, -2 means can not split on column_split_resp
    // this is needed to cast resp2master* to column_split_resp<T>*
    // as <T> needs to be decided first
    // "best_column_index" is not serialized for subclass "subtree_resp"

    string node_label; // need to keep this information, categorical value missing in train_set, but present in test_set
};

ibinstream& operator<<(ibinstream & m, const resp2master & resp) {
    m << resp.message_type;
    m << resp.best_column_index;
    m << resp.task_id;
    m << resp.node_label;

    return m;
}

obinstream& operator>>(obinstream & m,  resp2master & resp) {
    //m >> resp.message_type >> resp.best_column_index; should be done outside
    m >> resp.task_id;
    m >> resp.node_label;

    return m;
}

class subtree_resp : public resp2master {
public:
    TreeNode* root;

    subtree_resp() {
        message_type = SUB_TREE_RESP;
    }
};

ibinstream& operator<<(ibinstream& m, const subtree_resp & resp) {
    m << resp.message_type;
    m << resp.task_id;
    m << resp.root;

    return m;
}

obinstream& operator>>(obinstream& m, subtree_resp & resp) {
    // m >> resp.message_type
    m >> resp.task_id;
    m >> resp.root;

    return m;
}

template<class T>
class column_split_resp : public resp2master { // best column split
public:
    int slave_id; // this is to track who holds the best column

    Best_Split<T>* best_split; // T depends on "best_column_index", pointed to Task_Slave_Col_Split::best_split

    column_split_resp() {
        message_type = COL_SPLIT_RESP;
    }
};

// to get frequent y_label from master (tree_depth == MAX_TREE_DEPTH)
class resp_s2m_label : public resp2master {
public :

    resp_s2m_label() {
        message_type = RESP_LABEL_FETCH;
    }

    friend ibinstream & operator<<(ibinstream & m, const resp_s2m_label & resp) {
        m << resp.message_type;
        m << resp.task_id;
        m << resp.node_label;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, resp_s2m_label & resp) {
        // m >> resp.message_type;
        m >> resp.task_id;
        m >> resp.node_label;

        return m;
    }
};

template <class T>
ibinstream& operator<<(ibinstream& m, const column_split_resp<T> & resp) {
    m << resp.message_type;
    m << resp.best_column_index;
    m << resp.task_id;
    m << resp.slave_id;

    if(resp.best_column_index != -1) {
        m << resp.best_split;
    }

    m << resp.node_label;

    return m;
}

//!!! need to set "best_column_index" outside first, will use inside
//!!! need to new "best_split" outside first, will use inside
template <class T>
obinstream& operator>>(obinstream& m, column_split_resp<T> & resp) {
    //m >> resp.message_type;
    //m >> resp.best_column_index;
    m >> resp.task_id;
    m >> resp.slave_id;

    if(resp.best_column_index != -1) { //!!! using "best_column_index"
        m >> resp.best_split;
    }

    m >> resp.node_label;

    return m;
}

//////////////// Task master //////////////////

class Load_Content {
public:
    int row_idx;
    int col_idx;
    size_t value;

    Load_Content(int r, int c, size_t val) {
        row_idx = r;
        col_idx = c;
        value = val;
    }
};

class Task_Master { // task table's task object
// "sub_tree_task" uses it directly
// "column_split_task" will use its subclass Task_Master_Col_Split
public:
    int root_task_id;// to track each tree progress
    char task_type;
    TreeConfig tree_config;
    vector<int> column_indices;

    int task_id;
    TreeNode* & node;
    // for subtree-task, will store the subtree
    // for column_split_task, will store the node

    int size; // left_D / right_D, only for printing

    // for load balancing, keep track of the cost to adjust before deleting the task
    vector<Load_Content> task_load;

    int level = 0;

    const clock_t begin_time = clock();

    Task_Master(TreeNode* & p_node) : node(p_node) {
        //input is the location to hook node/tree
        task_type = TASK_SUB_TREE;
    }

    ~Task_Master() {
        tree_progress.decrement(root_task_id);
        float time_to_complete = float(clock() - begin_time) / CLOCKS_PER_SEC;
        cout << "task = " << task_id << ", time to complete = " << time_to_complete << " seconds" << endl;
    }
};

class Task_Master_Col_Split : public Task_Master { // column_split task
public:
    int best_column_index = -2; // -2 means uninitialized/can not split on set/subset, since -1 means a leaf
    double best_impurity = -DBL_MAX;

    int best_slave_id;
    vector<int> slave_ids;
    // the IDs of those slaves that are assigned to process columns for this task
    // will send_delete-plans to those slaves_id != best_slave_id

    size_t n_met = 0; // # of plan responses got back from slaves
    // when all the responses got back (slave_ids.size() == n_met), will create "node"

    void* best_split = NULL;
    string node_label;
    // keep track of the best splitting condition
    // need to be void* since we do not know <T> for Best<T>
    // once finalized, will use it to create "TreeNode*" (node)
    // <T> depends on "best_column_index"

    Task_Master_Col_Split(TreeNode* & p_node) : Task_Master(p_node) {
        //input is the location to hook node/tree
        task_type = TASK_COL_SPLIT;
    }
};

class Leaf_Task : public Task_Master { // to fetch most frequent y_label
public:
    string node_label;

    Leaf_Task(TreeNode* & p_node) : Task_Master(p_node) {
        task_type = TASK_LEAF;
    }
};

/////////////// plan from master to slave ///////////////

TreeNode* dummy;

class plan {
// used directly for delete-plan and end-plan
// for task-plans, use the 2 subclasses defined right after
public :
    int root_task_id = -1; // to track tree progress
    char message_type; // 4 types of plan message described above
    int task_id; // "node_id" (sub_tree or column_split), ("delete_task_id" for "delete_plan", for "delete_plan" no entry in "task_table")

    int dest_id = 0; // ### need it since a plan first gets into plan_buffer, and when main-thread takes it, it needs to now the dest_id to put into the right sending queue
    int level = 0;

    // slave informations to fetch data from
    int parent_task_id = -1; // -1 for root, task_id to fetch "row_indices" from.
    int parent_slave_id = -1; // -1 for root, "slave" to fetch from
    bool is_left = 3; // left_rows or right_rows
    int total_parent_req = 0;

    int size = 0; // left_D / right_D / _n_sample

    TreeNode* & node;

    plan(TreeNode* & tree_node) : node(tree_node){}

    friend ibinstream & operator<<(ibinstream & m, const plan & p_plan) {
        m << p_plan.message_type;
        m << p_plan.task_id;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, plan & p_plan) {
        // m >> p_plan.message_type;
        m >> p_plan.task_id;

        return m;
    }
};

class subtree_plan : public plan { // sub_tree plan
public :

    subtree_plan(TreeNode* & tree_node) : plan(tree_node) {
        message_type = SUB_TREE_PLAN;
    }

    vector<vector<int>> mac2colList;
    vector<int> column_indices; // subset of columns for subtree
    TreeConfig tree_config; // tree configuration

    friend ibinstream & operator<<(ibinstream & m, const subtree_plan & p_plan) {
        m << p_plan.message_type;
        m << p_plan.task_id;
        m << p_plan.level;
        m << p_plan.parent_task_id;
        m << p_plan.parent_slave_id;
        m << p_plan.is_left;
        m << p_plan.mac2colList;
        m << p_plan.tree_config;
        m << p_plan.total_parent_req;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, subtree_plan & p_plan) {
        // m >> p_plan.message_type;
        m >> p_plan.task_id;
        m >> p_plan.level;
        m >> p_plan.parent_task_id;
        m >> p_plan.parent_slave_id;
        m >> p_plan.is_left;
        m >> p_plan.mac2colList;
        m >> p_plan.tree_config;
        m >> p_plan.total_parent_req;

        return m;
    }
};

class column_split_plan : public plan { // column_split plan
public :

    column_split_plan(TreeNode* & tree_node) : plan(tree_node) {
        message_type = COL_SPLIT_PLAN;
    }

    vector<int> column_indices; // partial columns assigned to be on a single machine
    TreeConfig tree_config; // tree configuration

    friend ibinstream & operator<<(ibinstream & m, const column_split_plan & p_plan) {
        m << p_plan.message_type;
        m << p_plan.task_id;
        m << p_plan.parent_task_id;
        m << p_plan.parent_slave_id;
        m << p_plan.is_left;
        m << p_plan.column_indices;
        m << p_plan.tree_config;
        m << p_plan.total_parent_req;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, column_split_plan & p_plan) {
        // m >> p_plan.message_type;
        m >> p_plan.task_id;
        m >> p_plan.parent_task_id;
        m >> p_plan.parent_slave_id;
        m >> p_plan.is_left;
        m >> p_plan.column_indices;
        m >> p_plan.tree_config;
        m >> p_plan.total_parent_req;

        return m;
    }
};

class plan_m2s_label : public plan {
public:

    int n_local_cols = 0;

    plan_m2s_label(TreeNode* & tree_node) : plan(tree_node) {
        message_type = LABEL_FETCH_PLAN;
    }

    friend ibinstream & operator<<(ibinstream & m, const plan_m2s_label & plan) {
        m << plan.message_type;
        m << plan.task_id;
        m << plan.parent_task_id;
        m << plan.is_left;
        m << plan.n_local_cols;
        m << plan.total_parent_req;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, plan_m2s_label & plan) {
        // m >> req.message_type
        m >> plan.task_id;
        m >> plan.parent_task_id;
        m >> plan.is_left;
        m >> plan.n_local_cols;
        m >> plan.total_parent_req;

        return m;
    }
};

/////////////// request slave to slave machine /////////

class req_s2s {
public:
    char message_type; // "columns_attribute_value_fetch only for sub_tree" or
    // "row_indices_fetch for both "column_split" and "sub_tree_task"
    int task_id;
    int parent_req = 0;
};

class req_s2s_columns : public req_s2s { // for sub_tree column fetch request
public:
    int parent_task_id = -1;
    int parent_slave_id = -1;
    bool is_left = false;
    vector<int> column_indices; // subset of columns to request for, copied to "data_serve_task"

    req_s2s_columns() {
        message_type = REQ_COL_FETCH;
    }

    friend ibinstream & operator<<(ibinstream &m, const req_s2s_columns & req) {
        m << req.message_type;
        m << req.task_id;
        m << req.parent_task_id;
        m << req.parent_slave_id;
        m << req.is_left;
        m << req.column_indices;
        m << req.parent_req;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, req_s2s_columns & req) {
        //m >> req.message_type;
        m >> req.task_id;
        m >> req.parent_task_id;
        m >> req.parent_slave_id;
        m >> req.is_left;
        m >> req.column_indices;
        m >> req.parent_req;

        return m;
    }
};

class req_s2s_row_indices : public req_s2s { // requesting row_indices from parent node
public:
    int parent_task_id; //  req_s2s_columns::parent_task_id
    bool is_left = false; // req_s2s_columns::is_left
    int n_local_cols = 0; // size_of req_s2s_columns::column_indices
    // FixMe:: n_local_cols to process for these rows, (m * n_local_columns) == 2 * total_column, m = # of slaves

    req_s2s_row_indices() {
        message_type = REQ_ROW_FETCH;
    }

    friend ibinstream & operator<<(ibinstream & m, const req_s2s_row_indices & req) {
        m << req.message_type;
        m << req.task_id;
        m << req.parent_task_id;
        m << req.is_left;
        m << req.n_local_cols;
        m << req.parent_req;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, req_s2s_row_indices & req) {
        // m >> req.message_type
        m >> req.task_id;
        m >> req.parent_task_id;
        m >> req.is_left;
        m >> req.n_local_cols;
        m >> req.parent_req;

        return m;
    }
};


template <class T>
void serialize_dense_attributes(ibinstream &m, Column *col) {
    DenseColumn<T>* column = (DenseColumn<T>*) col;
    //m << column->size;
    m << column->data_type;
    m << column->is_ordinal;
    m << column->is_dense;

    m << column->elem;
    m << column->missing_value;

}

template <class T>
void serialize_sparse_attributes(ibinstream &m, Column *col) {
    SparseColumn<T>* column = (SparseColumn<T>*) col;
    m << column->data_type;
    m << column->is_ordinal;
    m << column->is_dense;

    m << column->size; // for sparse column
    m << column->elem;
    m << column->row_index;
    m << column->missing_value;
    m << column->default_value;
}

template <class T>
DenseColumn<T>* deserialize_dense_attribute(obinstream &m, size_t data_type, bool is_ordinal) {
    DenseColumn<T>* column = new DenseColumn<T>(data_type, is_ordinal);

    m >> column->elem;
    m >> column->missing_value;

    column->size = column->elem.size();

    return column;
}

template <class T>
SparseColumn<T>* deserialize_sparse_attribute(obinstream &m, size_t data_type, bool is_ordinal) {
    SparseColumn<T>* column = new SparseColumn<T>(data_type, is_ordinal);

    m >> column->size;
    m >> column->elem;
    m >> column->row_index;
    m >> column->missing_value;
    m >> column->default_value;

    return column;
}

class resp_s2s {
public :
    char message_type; // columns fetch or row indices fetch
    int task_id; // same task_id for row_fetch_req, column_fetch_req but would be interpreted differently on message_type
};


class resp_s2s_row_indices : public resp_s2s {
public :
    Candidate_Rows candidate_rows;

    resp_s2s_row_indices() {
        message_type = RESP_ROW_FETCH;
    }

    friend ibinstream & operator<<(ibinstream & m, const resp_s2s_row_indices & resp) {
        m << resp.message_type;
        m << resp.task_id;
        m << resp.candidate_rows;

        return m;
    }

    friend obinstream & operator>>(obinstream & m, resp_s2s_row_indices & resp) {
        // m >> resp.message_type;
        m >> resp.task_id;
        m >> resp.candidate_rows;

        return m;
    }
};

class resp_s2s_columns : public resp_s2s { // response message from target slave to source slave
public:
    vector<Column*> cols;

    resp_s2s_columns() {
        message_type = RESP_COL_FETCH;
    }

    friend ibinstream& operator<<(ibinstream &m, const resp_s2s_columns & resp) {
        m << resp.message_type;
        m << resp.task_id;
        m << resp.cols.size();

        for(size_t i = 0; i < resp.cols.size(); i++) {
            Column* temp = resp.cols[i];

#ifdef ASSERT
            assert(temp != nullptr);
#endif
            if(temp->is_dense) {
                if(temp->data_type == ELEM_BOOL) {
                    serialize_dense_attributes<bool>(m, temp);
                } else if(temp->data_type == ELEM_SHORT) {
                    serialize_dense_attributes<short>(m, temp);
                } else if(temp->data_type == ELEM_INT) {
                    serialize_dense_attributes<int>(m, temp);
                } else if(temp->data_type == ELEM_FLOAT) {
                    serialize_dense_attributes<float>(m, temp);
                } else if(temp->data_type == ELEM_DOUBLE) {
                    serialize_dense_attributes<double>(m, temp);
                } else if(temp->data_type == ELEM_CHAR) {
                    serialize_dense_attributes<char>(m, temp);
                } else if(temp->data_type == ELEM_STRING) {
                    serialize_dense_attributes<string>(m, temp);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else {
                if(temp->data_type == ELEM_BOOL) {
                    serialize_sparse_attributes<bool>(m, temp);
                } else if (temp->data_type == ELEM_SHORT) {
                    serialize_sparse_attributes<short>(m, temp);
                } else if (temp->data_type == ELEM_INT) {
                    serialize_sparse_attributes<int>(m, temp);
                } else if (temp->data_type == ELEM_FLOAT) {
                    serialize_sparse_attributes<float>(m, temp);
                } else if (temp->data_type == ELEM_DOUBLE) {
                    serialize_sparse_attributes<double>(m, temp);
                } else if (temp->data_type == ELEM_CHAR) {
                    serialize_sparse_attributes<char>(m, temp);
                } else if (temp->data_type == ELEM_STRING) {
                    serialize_sparse_attributes<string>(m, temp);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            }
        }

        return m;

    }

    friend obinstream& operator>>(obinstream &m, resp_s2s_columns & resp) {
        // m >> resp.message_type
        size_t n_col;
        int data_type;
        bool is_ordinal;
        bool is_dense;

        m >> resp.task_id;
        m >> n_col;

        resp.cols.resize(n_col);

        for(size_t i = 0; i < n_col; i++) {
            m >> data_type;
            m >> is_ordinal;
            m >> is_dense;

            if(is_dense) {
                if(data_type == ELEM_BOOL) {
                    DenseColumn<bool>* column = deserialize_dense_attribute<bool>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_SHORT) {
                    DenseColumn<short>* column = deserialize_dense_attribute<short>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_INT) {
                    DenseColumn<int>* column = deserialize_dense_attribute<int>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_FLOAT) {
                    DenseColumn<float>* column = deserialize_dense_attribute<float>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_DOUBLE) {
                    DenseColumn<double>* column = deserialize_dense_attribute<double>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_CHAR) {
                    DenseColumn<char>* column = deserialize_dense_attribute<char>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if(data_type == ELEM_STRING) {
                    DenseColumn<string>* column = deserialize_dense_attribute<string>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }

            } else {
                if(data_type == ELEM_BOOL) {
                    SparseColumn<bool>* column = deserialize_sparse_attribute<bool>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_SHORT) {
                    SparseColumn<short>* column = deserialize_sparse_attribute<short>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_INT) {
                    SparseColumn<int>* column = deserialize_sparse_attribute<int>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_FLOAT) {
                    SparseColumn<float>* column = deserialize_sparse_attribute<float>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_DOUBLE) {
                    SparseColumn<double>* column = deserialize_sparse_attribute<double>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_CHAR) {
                    SparseColumn<char>* column = deserialize_sparse_attribute<char>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else if (data_type == ELEM_STRING) {
                    SparseColumn<string>* column = deserialize_sparse_attribute<string>(m, data_type, is_ordinal);
                    resp.cols[i] = column;
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            }
        }

        return m;
    }
};

template <class T>
Column* copy_dense_data(Column* src, vector<size_t> &row_indices) {
    DenseColumn<T>* source = (DenseColumn<T>*) src;
    DenseColumn<T>* target = new DenseColumn<T>(source->data_type, source->is_ordinal);
    target->missing_value = source->missing_value;

    for(size_t i = 0; i < row_indices.size(); i++) {
        T value;
        source->get(row_indices[i], &value);
        target->append(i, &value);
    }

    return target;
}

template <class T>
Column* copy_sparse_data(Column* src, vector<size_t> &row_indices) {
    SparseColumn<T>* source = (SparseColumn<T>*) src;
    SparseColumn<T>* target = new SparseColumn<T>(source->data_type, source->is_ordinal, source->default_value);
    target->missing_value = source->missing_value;

    for(size_t i = 0; i < row_indices.size(); i++) {
        T value;
        source->get(row_indices[i], &value);
        if(value != source->default_value) {
            target->append(i, &value);
        }
    }

    target->finish(row_indices.size());

    return target;
}


Column* getColumn(int col_idx, vector<size_t> & row_indices) {
    Column* src = cserver.X.get_column(col_idx);

    if(src->is_dense) {
        if(src->data_type == ELEM_BOOL) {
            return copy_dense_data<bool>(src, row_indices);
        } else if(src->data_type == ELEM_SHORT) {
            return copy_dense_data<short>(src, row_indices);
        } else if(src->data_type == ELEM_INT) {
            return copy_dense_data<int>(src, row_indices);
        } else if(src->data_type == ELEM_FLOAT) {
            return copy_dense_data<float>(src, row_indices);
        } else if(src->data_type == ELEM_DOUBLE) {
            return copy_dense_data<double>(src, row_indices);
        } else if(src->data_type == ELEM_CHAR) {
            return copy_dense_data<char>(src, row_indices);
        } else if(src->data_type == ELEM_STRING) {
            return copy_dense_data<string>(src, row_indices);
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    } else {
        if(src->data_type == ELEM_BOOL) {
            return copy_sparse_data<bool>(src, row_indices);
        } else if(src->data_type == ELEM_SHORT) {
            return copy_sparse_data<short>(src, row_indices);
        } else if(src->data_type == ELEM_INT) {
            return copy_sparse_data<int>(src, row_indices);
        } else if(src->data_type == ELEM_FLOAT) {
            return copy_sparse_data<float>(src, row_indices);
        } else if(src->data_type == ELEM_DOUBLE) {
            return copy_sparse_data<double>(src, row_indices);
        } else if(src->data_type == ELEM_CHAR) {
            return copy_sparse_data<char>(src, row_indices);
        } else if(src->data_type == ELEM_STRING) {
            return copy_sparse_data<string>(src, row_indices);
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
            exit(-1);
        }
    }
}

void create_leaf(Task_Master* task_master) {
#ifdef ASSERT
    assert(task_master->task_type == TASK_COL_SPLIT || task_master->task_type == TASK_LEAF);
#endif

    string node_label;
    int size;

    if(task_master->task_type == TASK_COL_SPLIT) {
        Task_Master_Col_Split* task = (Task_Master_Col_Split*) task_master;
        node_label = task->node_label;
        size = task->size;
    } else if (task_master->task_type == TASK_LEAF) {
        Leaf_Task* task = (Leaf_Task*) task_master;
        node_label = task->node_label;
        size = task->size;
    } else {
        cout << "[ERROR] wrong type = " << task_master->task_type << ", File = " << __FILE__
             << ", Line = " << __LINE__ << endl;
        exit(-1);
    }

    int data_type =  cserver.X.get_column(y_index)->data_type;

    if(data_type == ELEM_BOOL) {
        LeafNode<bool>* leafNode = new LeafNode<bool>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<bool>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_SHORT) {
        LeafNode<short>* leafNode = new LeafNode<short>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<short>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_INT) {
        LeafNode<int>* leafNode = new LeafNode<int>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<int>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_FLOAT) {
        LeafNode<float>* leafNode = new LeafNode<float>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<float>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_DOUBLE) {
        LeafNode<double>* leafNode = new LeafNode<double>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<double>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_CHAR) {
        LeafNode<char>* leafNode = new LeafNode<char>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<char>(node_label);

        task_master->node = leafNode;

    } else if (data_type == ELEM_STRING) {
        LeafNode<string>* leafNode = new LeafNode<string>;
        leafNode->column_index = -1;
        leafNode->size = size;
        leafNode->label = boost::lexical_cast<string>(node_label);

        task_master->node = leafNode;

    } else {
        cout << "ERROR: File " << __FILE__ << ", Line " << __LINE__ << endl;
        exit(-1);
    }
}

void set_node_label(Matrix & dataset, vector<size_t> & rows, TreeConfig & treeConfig, string & node_label) {

    Column* Y = dataset.get_column(y_index);
    vector<size_t>::iterator start = rows.begin();
    vector<size_t>::iterator end = rows.end();


    if(Y->data_type == ELEM_BOOL) { // only categorical
        bool label;

        if(!Y->is_ordinal) {
            stop_splitting_categorical<bool>(cserver.X, start, end, treeConfig, label);
        } else {
            cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_SHORT) {
        short label;

        if(Y->is_ordinal) {
            stop_splitting_ordinal_int<short>(cserver.X, start, end, treeConfig, label);
        } else {
            stop_splitting_categorical<short>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_INT) {
        int label;

        if(Y->is_ordinal) {
            stop_splitting_ordinal_int<int>(cserver.X, start, end, treeConfig, label);
        } else {
            stop_splitting_categorical<int>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_FLOAT) {
        float label;

        if(Y->is_ordinal) {
            stop_splitting_ordinal<float>(cserver.X, start, end, treeConfig, label);
        } else {
            stop_splitting_categorical<float>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);


    } else if(Y->data_type == ELEM_DOUBLE) {
        double label;

        if(Y->is_ordinal) {
            stop_splitting_ordinal<double>(cserver.X, start, end, treeConfig, label);
        } else {
            stop_splitting_categorical<double>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_CHAR) { // only categorical
        char label;

        if(!Y->is_ordinal) {
            stop_splitting_categorical<char>(cserver.X, start, end, treeConfig, label);
        } else {
            cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_STRING) { // only categorical
        string label;

        if(!Y->is_ordinal) {
            stop_splitting_categorical<string>(cserver.X, start, end, treeConfig, label);
        } else {
            cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        node_label = label;

    } else {
        cout << "Error, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
        exit(-1);
    }

}

#endif
