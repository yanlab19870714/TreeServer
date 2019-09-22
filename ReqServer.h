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

#ifndef REQSERVER_H_
#define REQSERVER_H_

#include "tree_msgs.h"
#include <unistd.h> //for usleep()
#include <thread>
#include "RespQueue.h"
#include "ReqQueue.h"
#include "conmap2t.h"

using namespace std;

// This is to run on slaves
// - task_table is needed to get the task object for update states
// - q_resp is needed to send back response of reqs (e.g., after processing)
class ReqServer {
public:
    conmap2t<int, Task_Slave*> & task_table;
    RespQueue & q_resp; //will create _num_workers responding threads
    ReqQueue & q_req;
    thread main_thread;

    // directly respond rows back, since this is the slave that contains parent node
    void serve_rows(req_s2s_row_indices & req, int src) {
        //### create response using current task ID
        resp_s2s_row_indices* resp = new resp_s2s_row_indices;
        resp->task_id = req.task_id;

        //### parent is at this slave, get the parent task
        conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(req.parent_task_id);
        bucket.lock();

        unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();
        auto it = kvmap.find(req.parent_task_id);

#ifdef ASSERT
        if(it == kvmap.end()) {
            cout << "task not found = |" << req.parent_task_id << "| rank = " << _my_rank << " File = "
                 << __FILE__ << ", Line = " << __LINE__ << endl;
        }
        assert(it != kvmap.end()); //parent should be at this slave
        assert(it->second->task_type == TASK_COL_SPLIT); //parent is at this slave
#endif

        Task_Slave_Col_Split* parent_task = (Task_Slave_Col_Split*) it->second; //parent task
        parent_task->parent_request = req.parent_req;
        parent_task->n_cols += req.n_local_cols; //track task progress: how man cols processed

        //### get to the right child, and get its rows
        if(req.is_left) {
            resp->candidate_rows = parent_task->best_left_rows;
        } else {
            resp->candidate_rows = parent_task->best_right_rows;
        }

        //### put fetched rowIDs to send
        q_resp.add(resp, src);

        //### remove the task if all columns of both left and right children are served
        if(parent_task->n_cols == parent_task->parent_request) {

            if(parent_task->best_column_idx >= 0) {
                delete_best_split(parent_task->best_column_idx, parent_task->best_split); //garbage collect parent task's best-split
            }

            delete parent_task;

            kvmap.erase(it);

        }

        bucket.unlock();
    }

    void serve_columns(req_s2s_columns & req, int src) { // 2 stages for col-spilt tasks, 1st is to fetch row-indices, though called serve_columns
        resp_s2s_columns* resp = new resp_s2s_columns;

        resp->task_id = req.task_id;
        resp->cols.resize(req.column_indices.size());

        if(req.parent_task_id == -1) { //skip the parent hop, go directly to subtree-building slave
            //### special case: no parent-slave, just root node
            //fetch entire column values (since rows should be all of them)

            //we need to copy entire column for sending
            //since q_resp will delete the column after sent
            vector <size_t> row_indices;

            for (size_t i = 0; i < _n_samples; i++) {
                row_indices.push_back(i);
            }

            for (size_t i = 0; i < req.column_indices.size(); i++) {
                resp->cols[i] = getColumn(req.column_indices[i], row_indices);
            }

            //send the fetched entire row (for requested columns) back
            q_resp.add(resp, src);

        } else if (req.parent_slave_id == _my_rank) { //avoid communication if parent-task-data is local

            conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(req.parent_task_id);

            bucket.lock();

            unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();
            auto it = kvmap.find(req.parent_task_id);

#ifdef ASSERT
            assert(it != kvmap.end());
#endif

            vector<int> & cols = req.column_indices;
            Candidate_Rows rows;

            Task_Slave_Col_Split* parent_task = (Task_Slave_Col_Split*) it->second;
            parent_task->parent_request = req.parent_req;

            if(req.is_left) {
                rows = parent_task->best_left_rows;
            } else {
                rows =  parent_task->best_right_rows;
            }

            parent_task->n_cols += req.column_indices.size();

            if(!rows.stop_splitting) {
                for (size_t i = 0; i < cols.size(); i++) {
                    resp->cols[i] = getColumn(cols[i], rows.indexes);
                }

                //send the fetched entire row (for requested columns) back
                q_resp.add(resp, src);
            }

            //### remove the task if all columns of both left and right children are served
            if(parent_task->n_cols == parent_task->parent_request) {

                if(parent_task->best_column_idx >= 0) {
                    delete_best_split(parent_task->best_column_idx, parent_task->best_split); //garbage collect parent task's best-split
                }

                delete parent_task;

                kvmap.erase(it);
            }

            bucket.unlock();

        } else { // can happen both in MASTER and SLAVE nodes
            // delegate column fetch request
            //### this is only for a subtree task
            // *** this is at a slave where the key-slave (who builds the subtree) is requesting for column values ***
            //there are 2 steps: a. fetch row IDs, b. fetch column values
            //we are at step "a"

            //now requesting data toward parent-slave, need to wait for resp by waiting in slave-task-table
            Data_Serve_Task* task = new Data_Serve_Task; //we create the task object
            task->task_id = req.task_id;
            task->column_indices.swap(req.column_indices); //take task column IDs from req
            task->dest_id = src; // src is the slave to build subtree, who require columns from this slave

            //get the bucket to insert task
            conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(task->task_id);

            bucket.lock();

            unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();
            auto it = kvmap.find(task->task_id);

#ifdef ASSERT
            assert(it == kvmap.end());
            // a task assigned to a slave to build subtree, then that slave should request for column data from this slave only once
#endif

            bucket.insert(task->task_id, task);

            bucket.unlock();

            //now send the rowID request to parent-slave
            req_s2s_row_indices* req_rows = new req_s2s_row_indices;
            req_rows->task_id = req.task_id; //to allow resp to indicate which task is in processing
            req_rows->parent_task_id = req.parent_task_id;
            req_rows->is_left = req.is_left;
            req_rows->n_local_cols += task->column_indices.size(); //for progress check by parent-slave (to GC)
            req_rows->parent_req = req.parent_req;

            q_req.add(req_rows, req.parent_slave_id);
        }
    }

    void thread_func(char * buf, int size, int src)
    {
        obinstream m(buf, size);

        char message_type;

        while(m.end() == false)
        {
            m >> message_type;

            if(message_type == REQ_COL_FETCH) { // only for sub_tree task
                //### this is only for a subtree task
                // *** this is at a slave where the key-slave (who builds the subtree) is requesting for column values ***
                //there are 2 steps: a. fetch row IDs, b. fetch column values
                //we are at step "a"

                //### special case: no parent-slave, just root node
                //fetch entire column values (since rows should be all of them)

                // can happen in MASTER, only for data serving

                req_s2s_columns req;
                m >> req;

                serve_columns(req, src);

            } else if(message_type == REQ_ROW_FETCH) { // row_indices serve for both sub_tree + column_split tasks
                // *** this is at parent-slave ***
                //### this is a row-fetch request, meaning that:
                //current slave holds the parent node (and left/right split data)
                //parent just send back child's rowIDs
#ifdef ASSERT
                assert(_my_rank != MASTER_RANK);
#endif

                req_s2s_row_indices req;
                m >> req;

                serve_rows(req, src);

            } else {
                cout << "ERROR : File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }
        }

    }
    //*/

    void run()
    {
        bool first = true;
        thread t;
        //------
        while(global_end_label == false) //otherwise, thread terminates
        {
            int has_msg;
            MPI_Status status;
            MPI_Iprobe(MPI_ANY_SOURCE, REQ_CHANNEL, MPI_COMM_WORLD, &has_msg, &status);
            if(!has_msg) usleep(WAIT_TIME_WHEN_IDLE);
            else
            {
                int size;
                MPI_Get_count(&status, MPI_CHAR, &size); // get size of the msg-batch (# of bytes)
                char * buf = new char[size]; //space for receiving this msg-batch, space will be released by obinstream in thread_func(.)
                MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                if(!first) t.join(); //wait for previous CPU op to finish; t can be extended to a vector of threads later if necessary
                t = thread(&ReqServer::thread_func, this, buf, size, status.MPI_SOURCE); //insert to q_resp[status.MPI_SOURCE]
                first = false;
            }
        }
        if(!first) t.join();
    }

    ReqServer(conmap2t<int, Task_Slave*> & t_table, RespQueue & respQueue, ReqQueue & reqQueue)
        : task_table(t_table), q_resp(respQueue), q_req(reqQueue) //get local_table from Worker
    {
        main_thread = thread(&ReqServer::run, this);
    }

    ~ReqServer()
    {
        main_thread.join();
    }
};

#endif
