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

//this is the server of managing vcache
#ifndef RESPSERVER_H_
#define RESPSERVER_H_

//it receives batches of resps, add new comers to vcache
//- if capacity is less than limit, just insert
//- otherwise, try to replace; if cannot work, then insert (overflow); finally, try again to trim to limit
//this is repeated till no msg-batch is probed, in which case it sleeps for "WAIT_TIME_WHEN_IDLE" usec and probe again

#include "tree_msgs.h"
#include <unistd.h> //for usleep()
#include <thread>
#include <vector>
#include "conque_p.h"
#include "conmap2t.h"
#include "ReqServer.h"

using namespace std;

// there are 2 cases:

// === 1. subtree task ===
// 1. parent-server sends rowIDs back [RESP_ROW_FETCH]
// 1.1: <special case: this column-data-server is key-slave itself>: collect columns; if all cols ready, move to task buffer
// 1.2: fetch columns and send to key-slave
// 2. column-data-server sends column-data back to key-slave (that builds subtree) [RESP_COL_FETCH]: collect columns; if all cols ready, move to task buffer

// === 2. col-split task ===
// 1. parent-server sends rowIDs back [RESP_ROW_FETCH]: now have both rows and columns needed, move to task buffer
// // the tasks will then be taken by compers to find best for each assigned column, place global best to send-buffer

class RespServer {
public:
	typedef conmap2t<int, Task_Slave*> TaskTable;
	RespQueue & q_resp; //put fetched column-data for sending (only used by subtree tasks)
	TaskTable & task_table; //has been waiting for data
	conque_p<Task_Slave> & task_buffer; //to put task here for compers to process
	thread main_thread;

	void thread_func(char * buf, int size, int src)
	{
		//insert received batch
		obinstream m(buf, size);
		char resp_message_type;

		while(m.end() == false)
		{
			m >> resp_message_type;

			if(resp_message_type == RESP_COL_FETCH) {
				resp_s2s_columns resp;
				m >> resp;

				// get task object from slave-task-table
				conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(resp.task_id);
				bucket.lock();

				unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();

				auto it = kvmap.find(resp.task_id);

#ifdef ASSERT
				assert(it != kvmap.end()); //it must be there, since it's the 2nd step
#endif

				Task_Slave* task = it->second;
				// task obtained
#ifdef ASSERT
				assert(task->task_type == TASK_SUB_TREE);
#endif
				Task_Slave_Subtree* subtree_task = (Task_Slave_Subtree*) task;
				vector<int> & col_list = subtree_task->mac2colList[src]; // which cols are responded back?
				vector<Column*> & task_cols = subtree_task->matrix->col; // take task-local data to put columns
				for(size_t col_idx = 0; col_idx < col_list.size(); col_idx++) {
					task_cols[col_list[col_idx]] = resp.cols[col_idx];
				}
				subtree_task->n_met++; // track progress (all columns recv-ed?)


				if(subtree_task->n_met == subtree_task->n_req) { //all columns recv-ed
					task_buffer.enqueue(subtree_task);
					kvmap.erase(it);
				}

				bucket.unlock();

			} else if (resp_message_type == RESP_ROW_FETCH) { // both col_split and data_serve (sub_tree)
				resp_s2s_row_indices resp;
				m >> resp; //rowIDs are sent back

				conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(resp.task_id);
				bucket.lock();

				unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();

				auto it = kvmap.find(resp.task_id);
#ifdef ASSERT
				assert(it != kvmap.end());
#endif
				if(it->second->task_type == TASK_COL_SPLIT) { // this is the column-serving slave
					// parent-server sends rowIDs back [RESP_ROW_FETCH]: now have both rows and columns needed, move to task buffer

					Task_Slave_Col_Split* task = (Task_Slave_Col_Split*) it->second;
					task->candidate_rows = resp.candidate_rows;

					kvmap.erase(it);
					bucket.unlock();

					task_buffer.enqueue(task);

				} else if (it->second->task_type == TASK_DATA_SERVE) { //this is step 2 of subtree task (fetch rowIDs, fetch column values*)
					Data_Serve_Task* task = (Data_Serve_Task*) it->second;
					task->candidate_rows = resp.candidate_rows;

					if(!resp.candidate_rows.stop_splitting) {
						//create resp
						resp_s2s_columns* resp_col = new resp_s2s_columns;
						resp_col->task_id = task->task_id;
						resp_col->cols.resize(task->column_indices.size());

						//get values for assigned columns
						for(size_t i = 0; i < task->column_indices.size(); i++) {
							resp_col->cols[i] = getColumn(task->column_indices[i], task->candidate_rows.indexes);
						}

						q_resp.add(resp_col, task->dest_id);
					}

					kvmap.erase(it);
					delete task;

					bucket.unlock();

				} else if (it->second->task_type == TASK_SUB_TREE) {
					// columns local but rows remote, on arrival of remote rows
					Task_Slave_Subtree* subtree_task = (Task_Slave_Subtree*) it->second;
					subtree_task->candidate_rows = resp.candidate_rows;

					vector<size_t> & rows = resp.candidate_rows.indexes; // rowIDs are sent back
					vector<int> & cols = subtree_task->mac2colList[_my_rank]; // assigned columns are local to this key-slave


					if(cols.size() > 0) {
						for(size_t j = 0; j < cols.size(); j++) {

							Column* column = getColumn(cols[j], rows);
							subtree_task->matrix->col[cols[j]] = column;
						}

						subtree_task->n_met++; // local columns
					}

					subtree_task->n_met++;// row fetch for parent subtree_task

					if((subtree_task->n_req == subtree_task->n_met) || subtree_task->candidate_rows.stop_splitting) {
						kvmap.erase(it);
						task_buffer.enqueue(subtree_task);
					}

					bucket.unlock();

				} else {
					cout << "ERROR: , File = " << __FILE__ << ", Line = " << __LINE__ << endl;
					exit(-1);
				}

			} else {
				cout << "ERROR: , File = " << __FILE__ << ", Line = " << __LINE__ << endl;
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
    		MPI_Iprobe(MPI_ANY_SOURCE, RESP_CHANNEL, MPI_COMM_WORLD, &has_msg, &status);
    		if(!has_msg) usleep(WAIT_TIME_WHEN_IDLE);
    		else
    		{
    			int size;
    			MPI_Get_count(&status, MPI_CHAR, &size); // get size of the msg-batch (# of bytes)
    			char * buf = new char[size]; //space for receiving this msg-batch, space will be released by obinstream in thread_func(.)
    			MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			if(!first) t.join(); //wait for previous CPU op to finish; t can be extended to a vector of threads later if necessary
    			t = thread(&RespServer::thread_func, this, buf, size, status.MPI_SOURCE);
    			first = false;
    		}
    	}
    	if(!first) t.join();
    }

    RespServer(TaskTable & tt, conque_p<Task_Slave> & tbf, RespQueue & resp_queue)
		: task_table(tt), task_buffer(tbf), q_resp(resp_queue) //get cache_table from Worker
    {
    	main_thread = thread(&RespServer::run, this);
    }

    ~RespServer()
	{
    	main_thread.join();
	}
};

#endif
