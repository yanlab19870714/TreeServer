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

#ifndef MASTERRECEIVER_H_
#define MASTERRECEIVER_H_

#include "../tree_msgs.h"
#include "../conmap2t.h"
#include "PlanQueue.h"
#include "../ReqQueue.h"
#include "../deque_p.h"
#include <unistd.h> //for usleep()
#include <thread>

using namespace std;

// The MasterRecver class is where master received processed task results, there are 2 kinds:

// T1-subtree, for which the received node will be hooked to tree
// delete the task from task table

// T2-colSplit, for which we will aggregate the received "local best"
// case 1: "local best" not received yet from every slave, update task (in task table) in terms of "global best" & progress
// case 2: "local best" all received:
//          - create tree-node
//          - create delete-plan and add to plan_queue (for sending)
//          - create lchild-plan and rchild-plan, add to plan buffer (for further processing by main thread)
//          - delete task from task table


// Note:
//- planQueue is where master sends plans for slave processing
//- planBuffer is where master's main-thread gets plans for processing

//!!!!!!!!!!!!!!!! should remove: leaf does not check non-y columns, should do no column split
// T should be the label type
// this is to find label majority to create leaf node

//column_id -> machine_list mapping
vector<vector<int>> mac_map; //mac_map[column_id] = list of machine IDs
//need to resize() after loading metadata !!!

vector<int> all_columns; // all column indices except Y, initialize in Worker::run(), should be synchrinized afterwards

void random_shuffle(int k, vector<int> & result) {
#ifdef ASSERT
    assert(all_columns.size() == (_num_columns - 1));
#endif

    int num_columns = all_columns.size();

    for(int index = 0; index < k; index++) {
        int random_pos = index + (rand() % (num_columns - index));
        std::swap(all_columns[index], all_columns[random_pos]);
    }

    for(int index = 0; index < k; index++) {
        result.push_back(all_columns[index]);
    }

#ifdef ASSERT
    vector<bool> check(_num_columns, false);
    check[y_index] = true;

    for(int index = 0; index < num_columns; index++) {
        check[all_columns[index]] = true;
    }

    bool is_true = true;
    for(int index = 0; index < num_columns; index++) {
        is_true = is_true && check[index];
    }

    if(!is_true) {
        cout << "column_indices are : " << endl;
        for(size_t index = 0; index < all_columns.size(); index++) {
            cout << "|" << all_columns[index] << "|, found = " << check[index] << ", File = " << __FILE__
                 << ", Line = " << __LINE__ << endl;
        }
    }

    assert(is_true);
#endif
}

//machine loads

#define RECV_IDX 0 //  data receiver thread (REQ_CHANNEL & RESP_CHANNEL)
#define SEND_IDX 1 // data sender thread (REQ_CHANNEL & RESP_CHANNEL)
#define COMPT_IDX 2 // computation load in the machine

vector<vector<size_t>> load_matrix; // thread locked, only considering above 3
mutex load_matrix_lock;

void init_load_matrix() {
    unique_lock<mutex> lck(load_matrix_lock);

    load_matrix.resize(_num_workers);

    for(int i = 0; i < _num_workers; i++) {
        load_matrix[i].resize(3, 0); // 0 = receive, 1 = send, 2 = computation
    }
}

void adjust_load_matrix(vector<Load_Content> & task_load) {
    unique_lock<mutex> lck(load_matrix_lock);

    for(size_t index = 0; index < task_load.size(); index++) {
        Load_Content & l = task_load[index];
        load_matrix[l.row_idx][l.col_idx] -= l.value;
    }
}

int min_col_split_worker(size_t D, vector<int> & machine_indices, int parent_slave,
                         vector<bool> & parent_seen_flag) {
    // load_matrix locked before coming here

    int min_rank = machine_indices[0]; // default
    size_t min_load = INT_MAX;

    for(size_t i = 0; i < machine_indices.size(); i++) {
        int worker_id = machine_indices[i];

#ifdef ASSERT
        assert(worker_id != MASTER_RANK);
#endif

        // TODO:: need to introduce tunable param here
        size_t current_load = load_matrix[worker_id][RECV_IDX] + load_matrix[worker_id][COMPT_IDX]; // column_split_resp const
        size_t load_add = 0;

        load_add += D; // computation

        if(parent_slave != -1 && parent_slave != worker_id && !parent_seen_flag[worker_id]) {
            load_add += D; // receive
        }

        if(min_load >= (current_load + load_add)) {
            min_load = current_load + load_add;
            min_rank = worker_id;
        }
    }

    return min_rank;
}

void set_col_split_assignment(size_t D, vector<vector<int>> & mac2colList,
                              vector<int> & cols, int parent_slave, vector<Load_Content> & task_load) {

    // cases
    // 1. slaves |cols| from master (ignored now)
    // 2. slaves fetch (receive) |parent_rows| or |D| from "Pa" and "Pa" send |parent_rows| to slaves
    // 3. slaves compute (|cols| * |D|)
    // 4. slaves send "resp" to MASTER (ignored now)

    unique_lock<mutex> lck(load_matrix_lock);

    vector<bool> parent_fetch_flag(_num_workers, false);

    for(size_t i = 0; i < cols.size(); i++) {
        vector<int> & machine_indices = mac_map[cols[i]];
        int min_rank = min_col_split_worker(D, machine_indices, parent_slave, parent_fetch_flag);

        if(parent_slave != -1 && parent_slave != min_rank && !parent_fetch_flag[min_rank]) {
            parent_fetch_flag[min_rank] = true;

            load_matrix[parent_slave][SEND_IDX] += D; // case 2
            load_matrix[min_rank][RECV_IDX] += D; // case 2

            Load_Content l1(parent_slave, SEND_IDX, D);
            Load_Content l2(min_rank, RECV_IDX, D);
            task_load.push_back(l1); // case 2
            task_load.push_back(l2); // case 2
        }

        load_matrix[min_rank][COMPT_IDX] += D; // case 3

        Load_Content l(min_rank, COMPT_IDX, D);
        task_load.push_back(l); // case 3

        mac2colList[min_rank].push_back(cols[i]);
    }
}

int get_min_data_server(int D, vector<int> & machine_indices, int parent_slave, vector<bool> & parent_fetch_flag) {
    // load_matrix locked before coming here

    int min_rank = machine_indices[0];
    size_t min_load = INT_MAX;

    // case
    // 1. RECV parent rows |D|
    // 2. SEND back 1 * alpha * |D| column values to KS

    for(size_t index = 0; index < machine_indices.size(); index++) {
        int worker_id = machine_indices[index];

#ifdef ASSERT
        assert(worker_id != MASTER_RANK);
#endif
        size_t current_load = std::max(load_matrix[worker_id][RECV_IDX], load_matrix[worker_id][SEND_IDX]); // todo improve later
        size_t load_to_add = 0;

        if(parent_slave != -1 && parent_slave != worker_id && !parent_fetch_flag[worker_id]) {
            load_to_add += D; // 1 * alpha * D
        }

        if(min_load >= (current_load + load_to_add)) {
            min_load = current_load + load_to_add;
            min_rank = worker_id;
        }
    }

    return min_rank;
}

void set_subtree_assignment(subtree_plan* s_plan, vector<Load_Content> & task_load) {
    vector<int> & cols = s_plan->column_indices;
    s_plan->mac2colList.resize(_num_workers);

    unique_lock<mutex> lck(load_matrix_lock);

    // cases
    // 1. slaves receive req_s2s_columns from KS  and KS sends subtree to master (ignored now)
    // 2. slaves fetch (receive) |sample_rows| from parent (send)
    // 3. slaves send |sample| * |cols| to KS  and KS receives
    // 4. computation load on KS = |D| * n_columns * log2(|D|)
    // 5. KS receive |sample_rows| from parent (send)

    // first calculating computation rank
    size_t min_load = INT_MAX;
    int computation_rank = 1;

    for(int worker_id = 0; worker_id < _num_workers; worker_id++) {
        if(worker_id == MASTER_RANK) {
            continue;
        }

        if(min_load >= load_matrix[worker_id][COMPT_IDX]) {
            min_load = load_matrix[worker_id][COMPT_IDX];
            computation_rank = worker_id;
        }
    }

    s_plan->dest_id = computation_rank;

    int n_columns = s_plan->tree_config.n_columns;

    size_t temp_load = s_plan->size * n_columns * log2(s_plan->size);
    Load_Content l(computation_rank, COMPT_IDX, temp_load);
    task_load.push_back(l); // case 4

    load_matrix[computation_rank][COMPT_IDX] += temp_load; // case 4

    // case 5
    if(s_plan->dest_id != s_plan->parent_slave_id && s_plan->parent_slave_id != -1) {
        Load_Content l1(computation_rank, RECV_IDX, s_plan->size);
        Load_Content l2(s_plan->parent_slave_id, SEND_IDX, s_plan->size);
        task_load.push_back(l1);
        task_load.push_back(l2);

        load_matrix[computation_rank][RECV_IDX] += s_plan->size;
        load_matrix[s_plan->parent_slave_id][SEND_IDX] += s_plan->size;
    }

    // adding |row| fetching load from "Pa" while doing column assignment
    // "Pa" can also be a candidate worker for column assignment
    // "SEND" and "RECV" load need to adjusted "only once" for assigned |cols| in a worker
    // following parent_fetch_flag is used to keep track of already assigned worker
    vector<bool> parent_fetch_flag(_num_workers, false);

    for(size_t i = 0; i < cols.size(); i++) {
        vector<int> & mac_indices  = mac_map[cols[i]];
        int min_rank = get_min_data_server(s_plan->size, mac_indices, s_plan->parent_slave_id, parent_fetch_flag);
        s_plan->mac2colList[min_rank].push_back(cols[i]);

        if(min_rank == computation_rank) { // case 5
            continue;
        }

        Load_Content l1(min_rank, SEND_IDX, s_plan->size);
        Load_Content l2(s_plan->dest_id, RECV_IDX, s_plan->size);
        task_load.push_back(l1);// case 3
        task_load.push_back(l2); // case 3

        load_matrix[min_rank][SEND_IDX] += s_plan->size;// case 3
        load_matrix[s_plan->dest_id][RECV_IDX] += s_plan->size; // case 3

        if(!parent_fetch_flag[min_rank]) {
            parent_fetch_flag[min_rank] = true;

            if(s_plan->parent_slave_id != -1 && s_plan->parent_slave_id != min_rank) { // for each worker
                Load_Content l1(min_rank, RECV_IDX, s_plan->size);
                Load_Content l2(s_plan->parent_slave_id, SEND_IDX, s_plan->size);
                task_load.push_back(l1); // case 2
                task_load.push_back(l2); // case 2

                load_matrix[min_rank][RECV_IDX] += s_plan->size;// case 2
                load_matrix[s_plan->parent_slave_id][SEND_IDX] += s_plan->size; // case 2
            }
        }
    }
}

void add_leaf_load(vector<Load_Content> & task_load, plan_m2s_label* l_plan) {
    // case
    // 1. MASTER sends frequent label fetch request to slave ("Pa") and slave ("Pa") receives[ignored now]
    // 2. "Pa" compute on |D| (left/right)
    // 3. "Pa" sends response to MASTER [ignored now]

    unique_lock<mutex> lck(load_matrix_lock);

    load_matrix[l_plan->dest_id][COMPT_IDX] += l_plan->size; // case 2

    Load_Content l(l_plan->dest_id, COMPT_IDX, l_plan->size);
    task_load.push_back(l);// case 2
}

class MasterRecver {
public:

    PlanQueue & planQueue; //to put task-delete msgs for sending
    thread main_thread;
    conmap2t<int, Task_Master*> & task_table; // to get pending tasks for update & delete
    deque_p<plan> & plan_buffer; //to insert child plans

    template <class T>
    void task_update_col_split_with_leaf(obinstream & m) {
        //### to get the response msg
        column_split_resp<T> resp;

        resp.best_column_index = -1;  //used by m >> resp
        m >> resp;

        //### response msg tells which task object to update
        //### get task object from task_table
        conmap2t_bucket<int, Task_Master*> & bucket = task_table.get_bucket(resp.task_id);
        bucket.lock();

        unordered_map<int, Task_Master*> & kvmap = bucket.get_map();
        auto it = kvmap.find(resp.task_id);
#ifdef ASSERT
        assert(it != kvmap.end());
#endif
        Task_Master_Col_Split* task_split = (Task_Master_Col_Split*) it->second;
        task_split->n_met++; //### track that one resp of the task is received

        //### should be leaf node, remove any "best" if exist
        if(task_split->best_column_index >= 0) {
            delete_best_split(task_split->best_column_index, task_split->best_split);
        }

        //### update task object's fields using resp
        task_split->best_column_index = -1;
        task_split->best_slave_id = resp.slave_id; //may not be useful...

        //### all resps of the task are received
        //### task object's "best" is finalized
        if(task_split->n_met == task_split->slave_ids.size()) {
            //### create leaf node
            LeafNode<T>* leafNode = new LeafNode<T>;
            leafNode->column_index = -1;
            leafNode->size = task_split->size;
            leafNode->label = boost::lexical_cast<T>(resp.node_label);
            //### hook it to the tree
            task_split->node = leafNode; //note that task_split->node is a pointer reference
            //### delete entry and task

            adjust_load_matrix(task_split->task_load);

            kvmap.erase(it);
            delete task_split;
        }

        bucket.unlock();
    }

    column_split_plan* create_empty_plan(Task_Master_Col_Split* task_split, bool is_left, int D, int parent_req) {

        column_split_plan* c_plan;

        if(is_left) {
            c_plan = new column_split_plan(task_split->node->left);
        } else {
            c_plan = new column_split_plan(task_split->node->right);
        }

        int new_task_id = ++task_id_counter;
        c_plan->task_id = new_task_id;

        c_plan->root_task_id = task_split->root_task_id;
        c_plan->tree_config = task_split->tree_config;
        c_plan->level = task_split->level + 1;
        c_plan->parent_slave_id = task_split->best_slave_id;
        c_plan->parent_task_id = task_split->task_id;
        c_plan->is_left = is_left;
        c_plan->size = D;
        c_plan->total_parent_req = parent_req;

        if(!task_split->tree_config.sample_col_each_node) {
            c_plan->column_indices = task_split->column_indices;// passing from parent task, root task get it form root_plan
        }

        return c_plan;

    }

    subtree_plan* create_subtree_plan(Task_Master_Col_Split* task_split, bool is_left,
                                      int sample_size, int task_id, int parent_req) {
        subtree_plan* s_plan;

        if(is_left) {
            s_plan = new subtree_plan(task_split->node->left);
        } else {
            s_plan = new subtree_plan(task_split->node->right);
        }

        // for slave to fetch task object
        s_plan->task_id = task_id;
        s_plan->root_task_id = task_split->root_task_id;
        s_plan->tree_config = task_split->tree_config;

        // for slave to fetch rows from parent node's slave
        s_plan->parent_task_id = task_split->task_id;
        s_plan->parent_slave_id = task_split->best_slave_id;
        s_plan->is_left = is_left;

        s_plan->size = sample_size; // only used for printing
        s_plan->level = task_split->level + 1;
        s_plan->total_parent_req = parent_req;

        if(!task_split->tree_config.sample_col_each_node) {
            s_plan->column_indices = task_split->column_indices;// passing from parent task, root task get it form root_plan
        }

        return s_plan;
    }


    plan_m2s_label* create_label_plan(Task_Master_Col_Split* task_split, int D, bool is_left, int parent_req) {

        plan_m2s_label* leaf_plan;
        int new_task_id = ++task_id_counter; // get a task ID for left child node

        if(is_left) {
            leaf_plan = new plan_m2s_label(task_split->node->left);
        } else {
            leaf_plan = new plan_m2s_label(task_split->node->right);
        }

        leaf_plan->parent_task_id = task_split->task_id;
        leaf_plan->root_task_id = task_split->root_task_id;
        leaf_plan->task_id = new_task_id;
        leaf_plan->is_left = is_left;
        leaf_plan->dest_id = task_split->best_slave_id;
        leaf_plan->size = D;
        leaf_plan->level = task_split->level + 1;
        leaf_plan->n_local_cols += task_split->tree_config.n_columns;
        leaf_plan->total_parent_req = parent_req;

        return leaf_plan;

    }

    template <class T>
    void create_node_and_plans(Task_Master_Col_Split* task_split) {

        Best_Split<T>* task_best_split = (Best_Split<T>*) task_split->best_split;

        if(task_best_split->isOrdinal) {
            OrdinalTreeNode<T>* treeNode = new OrdinalTreeNode<T>();
            treeNode->column_index = task_split->best_column_index;
            treeNode->split_value = task_best_split->split_value;
            treeNode->size = task_split->size; // just for print
            treeNode->label = task_split->node_label;

            task_split->node = treeNode;
        } else {
            CategoricalTreeNode<T>* treeNode = new CategoricalTreeNode<T>();
            treeNode->column_index = task_split->best_column_index;

            if(task_best_split->S.size() > 0) {
                treeNode->S.insert(task_best_split->S.begin(), task_best_split->S.end());
            }

            if(task_best_split->S1.size() > 0) {
                treeNode->S1.insert(task_best_split->S1.begin(), task_best_split->S1.end());
            }

            treeNode->size = task_split->size; // just for print
            treeNode->label = task_split->node_label;

            task_split->node = treeNode;
        }

        // ------ append delete-plans for slaves to delete their task objects

        for(int i = 0; i < task_split->slave_ids.size(); i++) {
            if(task_split->best_slave_id != task_split->slave_ids[i]) {
                plan* d_plan = new plan(dummy); // no need to use dummy, global one in tree.h as placeholder
                d_plan->message_type = TASK_DELETE_PLAN;
                d_plan->task_id = task_split->task_id;
                d_plan->dest_id = task_split->slave_ids[i];
                planQueue.add(d_plan);
            }
        }

        TreeConfig & treeConfig = task_split->tree_config;

        int total_parent_request = treeConfig.n_columns * 2;

        if(task_split->level < treeConfig.MAX_TREE_DEPTH - 1) {
            if(task_best_split->left_D < subtree_D && task_best_split->left_D > treeConfig.MIN_SAMPLE_LEAF) {
                total_parent_request++;
            }

            if(task_best_split->right_D < subtree_D && task_best_split->right_D > treeConfig.MIN_SAMPLE_LEAF) {
                total_parent_request++;
            }
        }

        if(task_split->level == treeConfig.MAX_TREE_DEPTH - 1) {

            plan_m2s_label* left_leaf_plan = create_label_plan(task_split, task_best_split->left_D,
                                                               true, total_parent_request);
            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> left-leaf-task " << left_leaf_plan->task_id << " (p="
                 << left_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> left-leaf-task " << left_leaf_plan->task_id << " (p="
                 << left_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(left_leaf_plan->root_task_id);

            if(left_leaf_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(left_leaf_plan);
            } else {
                plan_buffer.push_back(left_leaf_plan);
            }

            plan_m2s_label* right_leaf_plan = create_label_plan(task_split, task_best_split->right_D,
                                                                false, total_parent_request);

            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> right-leaf-task " << right_leaf_plan->task_id << " (p="
                 << right_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> right-leaf-task " << right_leaf_plan->task_id << " (p="
                 << right_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(right_leaf_plan->root_task_id);

            if(right_leaf_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(right_leaf_plan);
            } else {
                plan_buffer.push_back(right_leaf_plan);
            }

            return;
        }


        if(task_best_split->left_D <= treeConfig.MIN_SAMPLE_LEAF) { // left leaf plan

            plan_m2s_label* left_leaf_plan = create_label_plan(task_split, task_best_split->left_D,
                                                               true, total_parent_request);

            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> left-leaf-task " << left_leaf_plan->task_id << " (p="
                 << left_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> left-leaf-task " << left_leaf_plan->task_id << " (p="
                 << left_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(left_leaf_plan->root_task_id);

            if(left_leaf_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(left_leaf_plan);
            } else {
                plan_buffer.push_back(left_leaf_plan);
            }

        } else if (task_best_split->left_D < subtree_D) { // left subtree plan

            // --- subtree plan as left child, one left-plan
            int new_task_id = ++task_id_counter;

            subtree_plan* left_plan = create_subtree_plan(task_split, true,
                                                          task_best_split->left_D,
                                                          new_task_id, total_parent_request);

            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> left-task " << left_plan->task_id << " (p="
                 << left_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> left-task " << left_plan->task_id << " (p="
                 << left_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(left_plan->root_task_id);

            if(left_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(left_plan);
            } else {
                plan_buffer.push_back(left_plan);
            }

        } else { //  left col_split plan

            // --- col-split plan as left child, a bunch of left-cplans

            column_split_plan* left_plan = create_empty_plan(task_split, true,
                                                             task_best_split->left_D, total_parent_request);

            tree_progress.increment(left_plan->root_task_id);

            if(left_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(left_plan);
            } else {
                plan_buffer.push_back(left_plan);
            }

        }

        if(task_best_split->right_D <= treeConfig.MIN_SAMPLE_LEAF) { // right leaf plan

            plan_m2s_label* right_leaf_plan = create_label_plan(task_split, task_best_split->right_D,
                                                                false, total_parent_request);
            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> right-leaf-task " << right_leaf_plan->task_id << " (p="
                 << right_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> right-leaf-task " << right_leaf_plan->task_id << " (p="
                 << right_leaf_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(right_leaf_plan->root_task_id);

            if(right_leaf_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(right_leaf_plan);
            } else {
                plan_buffer.push_back(right_leaf_plan);
            }

        } else if (task_best_split->right_D < subtree_D) { // right subtree plan

            // --- subtree plan as right child, one right-plan
            int new_task_id = ++task_id_counter;

            subtree_plan* right_plan = create_subtree_plan(task_split, false,
                                                           task_best_split->right_D,
                                                           new_task_id, total_parent_request);
            //>>>>>>>>>>>>>>>>>>>>>>>
#ifdef DEBUG_LOG
            cout << "[MASTER] masterRecver: --> right-task " << right_plan->task_id << " (p="
                 << right_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;

            fout << "[MASTER] masterRecver: --> right-task " << right_plan->task_id << " (p="
                 << right_plan->parent_task_id<<") into plan_buffer"
                 << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
#endif

            tree_progress.increment(right_plan->root_task_id);

            if(right_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(right_plan);
            } else {
                plan_buffer.push_back(right_plan);
            }

        } else { // right col_split plan

            // --- col-split plan as right child, a bunch of right-cplans

            column_split_plan* right_plan = create_empty_plan(task_split, false,
                                                              task_best_split->right_D, total_parent_request);
            tree_progress.increment(right_plan->root_task_id);

            if(right_plan->size < BFS_PRIORITY_THRESHOLD) {
                plan_buffer.push_front(right_plan);
            } else {
                plan_buffer.push_back(right_plan);
            }

        }

    }

    void create_node_and_plans_wrapper(Task_Master_Col_Split* task_split) {

        Column* best_column = cserver.X.get_column(task_split->best_column_index);
        int data_type = best_column->data_type;

        if(data_type == ELEM_BOOL) {
            create_node_and_plans<bool>(task_split);
        } else if (data_type == ELEM_SHORT) {
            create_node_and_plans<short>(task_split);
        } else if(data_type == ELEM_INT) {
            create_node_and_plans<int>(task_split);
        } else if (data_type == ELEM_FLOAT) {
            create_node_and_plans<float>(task_split);
        } else if (data_type == ELEM_DOUBLE) {
            create_node_and_plans<double>(task_split);
        } else if (data_type == ELEM_CHAR) {
            create_node_and_plans<char>(task_split);
        } else if (data_type == ELEM_STRING) {
            create_node_and_plans<string>(task_split);
        } else {
            cout << "ERROR: Wrong type, File " << __FILE__ << ", Line " << __LINE__ << endl;
            exit(-1);
        }

    }

    template <class T>
    void task_update_col_split(obinstream & m, int best_column_index) {

        //### to get the response msg
        column_split_resp<T> resp;
        //create "best" and hook to resp object, needed by m >> resp
        Best_Split<T>* b_split = new Best_Split<T>;
        resp.best_split = b_split;
        //needed by m >> resp
        resp.best_column_index = best_column_index;
        m >> resp;

        //### to get the task indicated by the response msg
        conmap2t_bucket<int, Task_Master*> & bucket = task_table.get_bucket(resp.task_id);
        bucket.lock();

        unordered_map<int, Task_Master*> & kvmap = bucket.get_map();
        auto it = kvmap.find(resp.task_id);

#ifdef ASSERT
        assert(it != kvmap.end());
#endif

        Task_Master_Col_Split* task_split = (Task_Master_Col_Split*) it->second;
        task_split->n_met++; //### track that one resp of the task is received

        //######################## update task object from one received response ########################
        if(task_split->best_column_index == -1) { //todo:: this may need to be removed later
            //------ Case 1: existing task object already received a leaf-response
            // leaf node has been hooked to task_split->node by a previous leaf-response
            // we just need to erase the table entry and task object


            if(task_split->n_met == task_split->slave_ids.size()) {
                adjust_load_matrix(task_split->task_load);
                kvmap.erase(it);
                delete task_split;
            }

            bucket.unlock();

            // delete best_split object in received response, as it's not useful
            delete b_split;

            return; //this case just ends

        }

        //now existing task object is not leaf-responsed

        if ((task_split->best_impurity < resp.best_split->best_impurity) ||
               (task_split->best_impurity == resp.best_split->best_impurity
                 && task_split->best_column_index > resp.best_column_index)
               ) {
            //------ Case 2: find a better col-splitting from response
            // erase old task_split->best_split

            if(task_split->best_column_index >= 0) {
                delete_best_split(task_split->best_column_index, task_split->best_split);
            }

            // update with new best-split
            task_split->best_impurity = resp.best_split->best_impurity;
            task_split->best_column_index = resp.best_column_index;
            task_split->best_split = resp.best_split;
            task_split->best_slave_id = resp.slave_id;
            task_split->node_label = resp.node_label;

        } else {
            //------ Case 3: new response is not better
            // keep task_split->best_xxx as is
            // delete best_split object in received response, as it's not useful
            delete resp.best_split;
        }

        //######################## if all responses received ########################
        // - create node to hook, using the finalized best-split
        // - append delete-plans for slaves to delete their task objects
        // - create and append left & right child-plans
        // - remove table entry and task object

#ifdef ASSERT
        assert(task_split->n_met <= task_split->slave_ids.size());
#endif

        if(task_split->n_met == task_split->slave_ids.size()) { // if all responses received
            // ------ create node to hook, using the finalized best-split

            // get current task object's best-split

            create_node_and_plans_wrapper(task_split);

            kvmap.erase(it);
            bucket.unlock();

            adjust_load_matrix(task_split->task_load);

            // removing master task
            if(task_split->best_column_index >= 0) {
                delete_best_split(task_split->best_column_index, task_split->best_split);
            }

            delete task_split;

        } else { //not the last response
            bucket.unlock();
        }

    }

    void skip_response(obinstream & m) {
        resp2master resp;
        resp.message_type = COL_SPLIT_RESP;
        resp.best_column_index = -2;

        m >> resp;

        conmap2t_bucket<int, Task_Master*> & bucket = task_table.get_bucket(resp.task_id);
        bucket.lock();

        unordered_map<int, Task_Master*> & kvmap = bucket.get_map();
        auto it = kvmap.find(resp.task_id);

        Task_Master_Col_Split* task_split = (Task_Master_Col_Split*) it->second;
        task_split->n_met++;

        if(task_split->n_met == task_split->slave_ids.size()) {
            kvmap.erase(it);

            if(task_split->best_column_index == -2) { // can not split on this subset of columns
                //### create leaf node
                task_split->node_label = resp.node_label;
                create_leaf(task_split);

                // ------ append delete-plans for slaves to delete their task objects

                for(size_t i = 0; i < task_split->slave_ids.size(); i++) {
                    plan* d_plan = new plan(dummy); // no need to use dummy, global one in tree.h as placeholder
                    d_plan->message_type = TASK_DELETE_PLAN;
                    d_plan->task_id = task_split->task_id;
                    d_plan->dest_id = task_split->slave_ids[i];
                    planQueue.add(d_plan);
                }

            } else {
                // here comes mean best split is not leaf
                create_node_and_plans_wrapper(task_split);
            }

            adjust_load_matrix(task_split->task_load);

            if(task_split->best_column_index >= 0) {
                delete_best_split(task_split->best_column_index, task_split->best_split);
            }

            delete task_split;
        }

        bucket.unlock();
    }

    void thread_func(char * buf, int size) {
        obinstream m(buf, size);
        char resp_type;

        while (m.end() == false) { // get each msg in this batch
            m >> resp_type;

            if(resp_type == SUB_TREE_RESP) { // subtree response
                subtree_resp resp;
                m >> resp;

                // get the task object for update & delete
                conmap2t_bucket<int, Task_Master*> & bucket = task_table.get_bucket(resp.task_id);
                bucket.lock();

                unordered_map<int, Task_Master*> & kvmap = bucket.get_map();
                auto it = kvmap.find(resp.task_id);

                // step 1: hook subtree to right place
                Task_Master* subtree_task = (Task_Master*) it->second;
                subtree_task->node = resp.root;

                // load balancing
                adjust_load_matrix(subtree_task->task_load);

                // step 2: delete the task from task table
                kvmap.erase(it); // remove table entry (task_id, task*)
                delete subtree_task; // free task space

                bucket.unlock();

            } else if (resp_type == COL_SPLIT_RESP) {
                int best_column_index;
                m >> best_column_index;
                // converting "m" to "column_split_resp" is deferred inside task_update_col_splitXXX()

                if(best_column_index == -2) { // can not find split on this response
                    skip_response(m);
                } else if(best_column_index == -1) { //resp2master::best_column_index == -1 means leaf
                    Column* column = cserver.X.get_column(y_index);

                    if(column->data_type == ELEM_BOOL) {
                        task_update_col_split_with_leaf<bool>(m);
                    } else if(column->data_type == ELEM_SHORT) {
                        task_update_col_split_with_leaf<short>(m);
                    } else if(column->data_type == ELEM_INT) {
                        task_update_col_split_with_leaf<int>(m);
                    } else if(column->data_type == ELEM_FLOAT) {
                        task_update_col_split_with_leaf<float>(m);
                    } else if(column->data_type == ELEM_DOUBLE) {
                        task_update_col_split_with_leaf<double>(m);
                    } else if(column->data_type == ELEM_CHAR) {
                        task_update_col_split_with_leaf<char>(m);
                    } else if(column->data_type == ELEM_STRING) {
                        task_update_col_split_with_leaf<string>(m);
                    } else {
                        cout << "Wrong type : , File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                        exit(-1);
                    }

                } else { //resp2master::best_column_index != -1 means its a column to split into left & right child-nodes

                    Column* column = cserver.X.get_column(best_column_index);

                    if(column->data_type == ELEM_BOOL) {
                        task_update_col_split<bool>(m, best_column_index);
                    } else if(column->data_type == ELEM_SHORT) {
                        task_update_col_split<short>(m, best_column_index);
                    } else if(column->data_type == ELEM_INT) {
                        task_update_col_split<int>(m, best_column_index);
                    } else if(column->data_type == ELEM_FLOAT) {
                        task_update_col_split<float>(m, best_column_index);
                    } else if(column->data_type == ELEM_DOUBLE) {
                        task_update_col_split<double>(m, best_column_index);
                    } else if(column->data_type == ELEM_CHAR) {
                        task_update_col_split<char>(m, best_column_index);
                    } else if(column->data_type == ELEM_STRING) {
                        task_update_col_split<string>(m, best_column_index);
                    } else {
                        cout << "Wrong type : , File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                        exit(-1);
                    }

                }

            } else if (resp_type == RESP_LABEL_FETCH) {
                resp_s2m_label resp;
                m >> resp;

                conmap2t_bucket<int, Task_Master*> & bucket = task_table.get_bucket(resp.task_id);
                bucket.lock();

                unordered_map<int, Task_Master*> & kvmap = bucket.get_map();
                auto it = kvmap.find(resp.task_id);

                Leaf_Task* leaf_task = (Leaf_Task*) it->second;
                leaf_task->node_label = resp.node_label;

                create_leaf(leaf_task); // set leaf node here

                kvmap.erase(it);
                bucket.unlock();

                // load balancing
                adjust_load_matrix(leaf_task->task_load);

                // delete leaf_task
                delete leaf_task;

            } else {
                cout << "Wrong resp_type " << resp_type << " : , File = " << __FILE__
                     << ", Line = " << __LINE__ << endl;

                exit(-1);
            }
        }
    }

    void run()
    {
        bool first = true;
        thread t;
        //------
        while(global_end_label == false) //otherwise, thread terminates
        {
            int has_msg;
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, PLAN_CHANNEL, MPI_COMM_WORLD, &has_msg, &status);
            if(!has_msg) usleep(WAIT_TIME_WHEN_IDLE);
            else
            {
                int size;
                MPI_Get_count(&status, MPI_CHAR, &size); // get size of the msg-batch (# of bytes)
                char * buf = new char[size]; //space for receiving this msg-batch, space will be released by obinstream in thread_func(.)
                MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(!first) t.join(); //wait for previous CPU op to finish; t can be extended to a vector of threads later if necessary
                t = thread(&MasterRecver::thread_func, this, buf, size);
                first = false;
            }
        }
        if(!first) t.join();
    }

    MasterRecver(conmap2t<int, Task_Master*> & p_table, PlanQueue & pplanQueue,
                 deque_p<plan> & p_buffer) : task_table(p_table), planQueue(pplanQueue),
                                              plan_buffer(p_buffer) {

        main_thread = thread(&MasterRecver::run, this);

    }

    ~MasterRecver() {
        main_thread.join();
    }
};

#endif
