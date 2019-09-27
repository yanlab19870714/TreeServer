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

#ifndef COMPER_H
#define COMPER_H

#include <unistd.h> //for usleep()
#include <thread>
#include "../tree_msgs.h"
#include "../conque_p.h"
#include "../conmap2t.h"

using namespace std;

void set_stop_cond_and_indices(vector<size_t>::iterator start, vector<size_t>::iterator end, TreeConfig & treeConfig,
                               Candidate_Rows & best_rows) {

    Column* Y = cserver.X.get_column(y_index);

    string & node_label = best_rows.freq_y;
    bool & stop = best_rows.stop_splitting;
    best_rows.indexes.insert(best_rows.indexes.end(), start, end); // todo scope to improve

    //if should be leaf, then add resp (which send the "best", but means leaf-node)
    if(Y->data_type == ELEM_BOOL) { // only categorical
        bool label;

        if(!Y->is_ordinal) {
            stop = stop_splitting_categorical<bool>(cserver.X, start, end, treeConfig, label);
        } else {
            cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_SHORT) {
        short label;

        if(Y->is_ordinal) {
            stop = stop_splitting_ordinal<short>(cserver.X, start, end, treeConfig, label);
        } else {
            stop = stop_splitting_categorical<short>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_INT) {
        int label;

        if(Y->is_ordinal) {
            stop = stop_splitting_ordinal<int>(cserver.X, start, end, treeConfig, label);
        } else {
            stop = stop_splitting_categorical<int>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_FLOAT) {
        float label;

        if(Y->is_ordinal) {
            stop = stop_splitting_ordinal<float>(cserver.X, start, end, treeConfig, label);
        } else {
            stop = stop_splitting_categorical<float>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_DOUBLE) {
        double label;

        if(Y->is_ordinal) {
            stop = stop_splitting_ordinal<double>(cserver.X, start, end, treeConfig, label);
        } else {
            stop = stop_splitting_categorical<double>(cserver.X, start, end, treeConfig, label);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_CHAR) { // only categorical
        char label;

        if(!Y->is_ordinal) {
            stop = stop_splitting_categorical<char>(cserver.X, start, end, treeConfig, label);
        } else {
            cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        node_label = to_string(label);

    } else if(Y->data_type == ELEM_STRING) { // only categorical
        string label;

        if(!Y->is_ordinal) {
            stop = stop_splitting_categorical<string>(cserver.X, start, end, treeConfig, label);
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

void create_end_split_response(Task_Slave_Col_Split* task, conque_p<resp2master> & ready_buffer) {
    Column* Y = cserver.X.get_column(y_index);
    int y_data_type = Y->data_type;

    if(y_data_type == ELEM_BOOL) {
        column_split_resp<bool>* resp = new column_split_resp<bool>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_SHORT) {
        column_split_resp<short>* resp = new column_split_resp<short>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_INT) {
        column_split_resp<int >* resp = new column_split_resp<int>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_FLOAT) {
        column_split_resp<float>* resp = new column_split_resp<float>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_DOUBLE) {
        column_split_resp<double>* resp = new column_split_resp<double>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_CHAR) {
        column_split_resp<char>* resp = new column_split_resp<char>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else if(y_data_type == ELEM_STRING) {
        column_split_resp<string>* resp = new column_split_resp<string>;
        resp->best_column_index = -1;
        resp->task_id = task->task_id;
        resp->slave_id = _my_rank;
        resp->node_label = task->candidate_rows.freq_y;

        ready_buffer.enqueue(resp);
    } else {
        cout << "Error: Wrong type " << y_data_type << ", File = " << __FILE__ << ", Line = " << __LINE__ << endl;
        exit(-1);
    }

    if(task->best_column_idx >= 0) {
        delete_best_split(task->best_column_idx, task->best_split);
    }

    delete task;
}

class Comper {
public:
    thread main_thread;
    conque_p<Task_Slave> & task_buffer;
    conque_p<resp2master> & ready_buffer; //for sending subtree or "best"
    conmap2t<int, Task_Slave*> & slave_task_table;
    //after "best" is sent, move task back to slave_task_table
    //wait for delete (by delete-plan or by progress-counter)

    void build_subtree(Task_Slave* task_slave) {
        //take the task object, using its data to build subtree
        //put a resp containing the subtree, and remove the task

        //take the task object
        Task_Slave_Subtree* subtree_task = (Task_Slave_Subtree*) task_slave;

        if(subtree_task->candidate_rows.stop_splitting) {
            subtree_resp* resp = new subtree_resp;
            resp->task_id = subtree_task->task_id;

            Column* Y = cserver.X.get_column(y_index);
            vector<size_t> & rows = subtree_task->candidate_rows.indexes;
            subtree_task->root = create_leaf_wrapper(Y, rows.begin(), rows.end(), subtree_task->candidate_rows.freq_y);

            resp->root = subtree_task->root;
            ready_buffer.enqueue(resp);
            delete subtree_task;

        } else {
            //create resp
            subtree_resp* resp = new subtree_resp;
            resp->task_id = subtree_task->task_id;

            //get data from task object
            Matrix & task_data = *(subtree_task->matrix);

            task_data.col[y_index] = getColumn(y_index, subtree_task->candidate_rows.indexes);

            vector<size_t> rows;

            for(size_t i = 0; i < subtree_task->candidate_rows.indexes.size(); i++) {
                rows.push_back(i);
            }

            subtree_task->root = build_tree(subtree_task, rows.begin(), rows.end(), subtree_task->level);
            resp->root = subtree_task->root;
            ready_buffer.enqueue(resp);
            delete subtree_task;
        }

    }

    //create resp with "best"
    //add task back to slave-task-table, to be deleted or queried for rowIDs
    template <class T>
    void set_best_split(vector<size_t>::iterator start, vector<size_t>::iterator end, Column* column,
                        SplitResult* splitResult, Task_Slave_Col_Split* task_slave) {

        string & node_label = task_slave->candidate_rows.freq_y;
        Best_Split<T>* best_split = new Best_Split<T>(); // deleted in SubTreeSender communication thread

        //take "best" info from "splitResult" to best_split
        best_split->best_impurity = splitResult->best_impurity;

        if(column->is_ordinal) {
            OrdinalSplitResult<T>* result = (OrdinalSplitResult<T>*) splitResult;
            best_split->isOrdinal = true;
            best_split->split_value = result->attribute_value;
            best_split->left_D = splitResult->pos - start;
            best_split->right_D = end - splitResult->pos;
        } else {
            CategoricalSplitResult<T>* result = (CategoricalSplitResult<T>*) splitResult;
            best_split->isOrdinal = false;
            best_split->S = result->S;
            best_split->S1 = result->S1;
            best_split->left_D = splitResult->pos - start;
            best_split->right_D = end - splitResult->pos;
        }

        //give it to task object
        task_slave->best_split = best_split;
        task_slave->best_column_idx = splitResult->column_idx;

        //to set child-node info, will be needed to serve rowIDs later (if I am the best)
        set_stop_cond_and_indices(start, splitResult->pos, task_slave->tree_config, task_slave->best_left_rows);
        set_stop_cond_and_indices(splitResult->pos, end, task_slave->tree_config, task_slave->best_right_rows);

        //put task back to task-table
        conmap2t_bucket<int, Task_Slave*> & bucket = slave_task_table.get_bucket(task_slave->task_id);
        bucket.lock();
        bucket.insert(task_slave->task_id, task_slave);
        bucket.unlock();

        //create resp
        column_split_resp<T>* c_resp = new column_split_resp<T>;
        c_resp->task_id = task_slave->task_id;
        c_resp->slave_id = _my_rank; //track parent id for querying rowIDs
        c_resp->best_column_index = task_slave->best_column_idx;
        c_resp->best_split = (Best_Split<T>*) task_slave->best_split;
        c_resp->node_label = node_label;

        ready_buffer.enqueue(c_resp);
    }

    void stop_splitting(vector<size_t>::iterator start, vector<size_t>::iterator end,
                        TreeConfig & treeConfig, Task_Slave_Col_Split* task) {

        Column* Y = cserver.X.get_column(y_index);
        string & node_label = task->candidate_rows.freq_y;

        bool & stop = task->candidate_rows.stop_splitting;

        //if should be leaf, then add resp (which send the "best", but means leaf-node)
        if(Y->data_type == ELEM_BOOL) { // only categorical
            bool label;

            if(!Y->is_ordinal) {
                stop = stop_splitting_categorical<bool>(cserver.X, start, end, treeConfig, label);
            } else {
                cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<bool>* resp = new column_split_resp<bool>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);
            }

        } else if(Y->data_type == ELEM_SHORT) {
            short label;

            if(Y->is_ordinal) {
                stop = stop_splitting_ordinal<short>(cserver.X, start, end, treeConfig, label);
            } else {
                stop = stop_splitting_categorical<short>(cserver.X, start, end, treeConfig, label);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<short>* resp = new column_split_resp<short>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);
            }


        } else if(Y->data_type == ELEM_INT) {
            int label;

            if(Y->is_ordinal) {
                stop = stop_splitting_ordinal<int>(cserver.X, start, end, treeConfig, label);
            } else {
                stop = stop_splitting_categorical<int>(cserver.X, start, end, treeConfig, label);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<int>* resp = new column_split_resp<int>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);

            }


        } else if(Y->data_type == ELEM_FLOAT) {
            float label;

            if(Y->is_ordinal) {
                stop = stop_splitting_ordinal<float>(cserver.X, start, end, treeConfig, label);
            } else {
                stop = stop_splitting_categorical<float>(cserver.X, start, end, treeConfig, label);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<float>* resp = new column_split_resp<float>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);

            }


        } else if(Y->data_type == ELEM_DOUBLE) {
            double label;

            if(Y->is_ordinal) {
                stop = stop_splitting_ordinal<double>(cserver.X, start, end, treeConfig, label);
            } else {
                stop = stop_splitting_categorical<double>(cserver.X, start, end, treeConfig, label);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<double>* resp = new column_split_resp<double>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);

            }

        } else if(Y->data_type == ELEM_CHAR) { // only categorical
            char label;

            if(!Y->is_ordinal) {
                stop = stop_splitting_categorical<char>(cserver.X, start, end, treeConfig, label);
            } else {
                cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }

            node_label = to_string(label);

            if(stop) {
                column_split_resp<char>* resp = new column_split_resp<char>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);

            }

        } else if(Y->data_type == ELEM_STRING) { // only categorical
            string label;

            if(!Y->is_ordinal) {
                stop = stop_splitting_categorical<string>(cserver.X, start, end, treeConfig, label);
            } else {
                cout << "Error: Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }

            node_label = label;

            if(stop) {
                column_split_resp<string>* resp = new column_split_resp<string>;
                resp->best_column_index = -1;
                resp->task_id = task->task_id;
                resp->slave_id = _my_rank;
                resp->node_label = node_label;

                ready_buffer.enqueue(resp);
            }

        } else {
            cout << "Error, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }

        if(stop) { //if leaf, no need task and its best-split, GC
            if(task->best_column_idx >= 0) {
                delete_best_split(task->best_column_idx, task->best_split);
            }

            delete task;
        }

    }

    void find_best_split(Task_Slave* task_slave) { //take task object, compute the best split

        //take task object
        Task_Slave_Col_Split* task = (Task_Slave_Col_Split*) task_slave;

        if(task->candidate_rows.stop_splitting) {
            create_end_split_response(task, ready_buffer);
            return;
        }

        vector<int> & cols = task->column_indices;
        vector<size_t> & rows = task->candidate_rows.indexes;

        TreeConfig & treeConfig = task->tree_config;

        vector<size_t>::iterator start = rows.begin();
        vector<size_t>::iterator end = rows.end();

        if(task->is_root_task) {
            stop_splitting(start, end, treeConfig, task);

            if(task->candidate_rows.stop_splitting) { // leaf case already handled, return
                return;
            }
        }

        //find best-split
        SplitResult* splitResult = node_split(cserver.X, start, end, cols, treeConfig);

        if(splitResult->column_idx == -2) { // can not split on assigned columns

            // inserting because of delete_plan to come from MASTER
            conmap2t_bucket<int, Task_Slave*> & bucket = slave_task_table.get_bucket(task_slave->task_id);
            bucket.lock();
            bucket.insert(task_slave->task_id, task_slave);
            bucket.unlock();

            resp2master* resp = new resp2master;

            resp->message_type = COL_SPLIT_RESP;
            resp->best_column_index = -2;
            resp->task_id = task->task_id;
            resp->node_label = task->candidate_rows.freq_y;

            ready_buffer.enqueue(resp);

            delete splitResult;

            return;
        }

        Column* column = cserver.X.get_column(splitResult->column_idx);

        int data_type = column->data_type;

        //prepare resp containing best-split to send to master
        if(data_type == ELEM_BOOL) {
            set_best_split<bool>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_SHORT) {
            set_best_split<short>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_INT) {
            set_best_split<int>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_FLOAT) {
            set_best_split<float>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_DOUBLE) {
            set_best_split<double>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_CHAR) {
            set_best_split<char>(start, end, column, splitResult, task);
        } else if(data_type == ELEM_STRING) {
            set_best_split<string>(start, end, column, splitResult, task);
        } else {
            cout << "Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
            exit(-1);
        }


        free(splitResult);
    }

    void run() {

        while (global_end_label == false) {
            bool sth_found = false;

            Task_Slave* task_slave = task_buffer.dequeue();

            while (task_slave != NULL) {
                sth_found = true;

                cout<<"in comper::run(), task_slave->tree_config.type = "<<task_slave->tree_config.type<<endl;//@@@@@@@@@@@@@@@@@@@@@@@

                if(task_slave->task_type == TASK_SUB_TREE) {
                    build_subtree(task_slave);
                } else if (task_slave->task_type == TASK_COL_SPLIT) {
                    find_best_split(task_slave);
                } else {
                    cout << "ERROR : File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                    exit(-1);
                }

                task_slave = task_buffer.dequeue();

            }

            if(sth_found == false) {
                usleep(WAIT_TIME_WHEN_IDLE);
            } else {
                sth_found = false;
            }
        }

    }

    Comper(conque_p<Task_Slave> &t_buff, conque_p<resp2master> & r_buff, conmap2t<int, Task_Slave*> & t_table)
        : task_buffer(t_buff), ready_buffer(r_buff), slave_task_table(t_table) {
        main_thread = thread(&Comper::run, this);
    }

    ~Comper() {
        main_thread.join();
    }
};

#endif
