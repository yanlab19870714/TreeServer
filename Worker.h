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

#ifndef WORKER_H_
#define WORKER_H_

#include <iostream>
#include <atomic>
#include <cassert>
#include <algorithm>
#include <queue>
#include "ydhdfs.h"
#include "ReqServer.h"
#include "RespServer.h"
#include "ReqQueue.h"
#include "tree_msgs.h"
#include "conque_p.h"
#include "deque_p.h"
#include "conque.h"
#include "master/MasterRecver.h"
#include "conmap2t.h"
#include "slave/comper.h"
#include "slave/SlaveSender.h"
#include "master/PlanQueue.h"

using namespace std;

const int MODE_ALL_COLUMNS = 20; //a column is in every machine
const int MODE_REPLICATE = 21; //replication factor k: a column is in k machines

struct WorkerParams {
    // configuration parameters
	string job_file_path = "";
	string tree_file_path = "";
	string test_file_path = "";
	string test_meta_file = "";

    int column_assignment;
    int replicate;
};

struct min_heap_entry //for MODE_REPLICATE
{
	size_t count;
	int rank;

	bool operator<(const min_heap_entry& o) const
	{
		return count > o.count;
	}
};

class Worker
{
public:

    Matrix test_set; // per job 1 test data file for now
    //### for debug
    void print_mac_map()
    {
    	for(size_t i=0; i<mac_map.size(); i++)
		{
    		cout<<"col "<<i<<": ";
			vector<int> & macs = mac_map[i];
			for(size_t j=0; j<macs.size(); j++) cout<<macs[j]<<" ";
			cout<<endl;
		}
    }

	Worker()
    {
    	init_worker(NULL, NULL);

        if(_my_rank == MASTER_RANK) {
            init_load_matrix();
        }
	}

    ~Worker()
	{
#ifdef DEBUG_LOG
		fout.close();// todo debug remove later
#endif
		worker_finalize();
	}

    //only called by master
	void load_tree_configs(const WorkerParams & params, vector<TreeConfig> & configList) {
		load_config(params.tree_file_path.c_str(), configList); //then load tree configs
		_mkdir(job_dir.c_str());

#ifdef DEBUG_LOG
		string log_dir = job_dir + "/logs";
		_mkdir(log_dir.c_str());
        string log_file_path = log_dir + "/worker_" + to_string(_my_rank) + ".log";
        fout.open(log_file_path.c_str());

		print_job_config();
		fout<<configList.size()<<" tree-configs loaded"<<endl;
		fout<<endl;
#endif
	}

	// for master start ////////////////////////////

	subtree_plan* create_subtree_plan(TreeConfig & plan_config, TreeNode* & root, vector<int> & column_indices) {

		subtree_plan* p_plan = new subtree_plan(root);

		int task_id = ++task_id_counter;

		p_plan->task_id = task_id;
		p_plan->tree_config = plan_config;

		// for slave to fetch rows from parent node's slave
		p_plan->parent_task_id = -1;
		p_plan->parent_slave_id = -1;
		//p_plan->is_left = is_left;

		p_plan->size = _n_samples;
		p_plan->level = 0;
		p_plan->column_indices = column_indices;

		return p_plan;
	}

	void create_root_subtree_plans(queue<plan*> & root_plan_buffer, vector<TreeConfig> & configList) {

		for(size_t index = 0; index < configList.size(); index++) {
			TreeConfig & treeConfig = configList[index];

			if(treeConfig.type == DECISION_TREE) {
				TreeConfig p_config;
				set_config(treeConfig, p_config, 0);

				subtree_plan* subtreePlan = create_subtree_plan(p_config, treeConfig.rootList[0],
																  treeConfig.column_distribution[0]);
				subtreePlan->root_task_id = subtreePlan->task_id;

				root_plan_buffer.push(subtreePlan);

			} else if (treeConfig.type == RANDOM_FOREST) {
				int n_plans = treeConfig.column_distribution.size();

				for(int k = 0; k < n_plans; k++) {
					TreeConfig p_config;
					set_config(treeConfig, p_config, k);

					subtree_plan* subtreePlan = create_subtree_plan(p_config, treeConfig.rootList[k],
																	  treeConfig.column_distribution[k]);
					subtreePlan->root_task_id = subtreePlan->task_id;

					root_plan_buffer.push(subtreePlan);
				}
			} else {
				cout << "ERROR: wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
				exit(-1);
			}
		}
	}

	column_split_plan* create_empty_plan(TreeConfig & plan_config, TreeNode* & root, vector<int> & column_indices) {

		column_split_plan* empty_plan = new column_split_plan(root);

		int new_task_id = ++task_id_counter;

		empty_plan->task_id = new_task_id;
		empty_plan->tree_config = plan_config;
		empty_plan->level = 0;
		empty_plan->parent_slave_id = -1;
		empty_plan->parent_task_id = -1;
		empty_plan->column_indices = column_indices;
		//c_plan->is_left = 3;
		empty_plan->size = _n_samples;

		return empty_plan;
	}

	void create_root_col_split_plan(queue<plan*> & root_plan_buffer, vector<TreeConfig> & configList) {
        //todo:: right now just 1 root to build a tree, later need to load a list of roots (those of random/deep forest) from config-file

        for(size_t i = 0; i < configList.size(); i++) {
			TreeConfig & treeConfig = configList[i];

			if(treeConfig.type == DECISION_TREE) {
				TreeConfig p_config;
				set_config(treeConfig, p_config, 0);

				column_split_plan* empty_plan = create_empty_plan(p_config, treeConfig.rootList[0],
																  treeConfig.column_distribution[0]);
                empty_plan->root_task_id = empty_plan->task_id;

                root_plan_buffer.push(empty_plan);


			} else if (treeConfig.type == RANDOM_FOREST) {
				int n_plans = treeConfig.column_distribution.size();

				for(int k = 0; k < n_plans; k++) {
					TreeConfig p_config;
					set_config(treeConfig, p_config, k);

					column_split_plan* empty_plan = create_empty_plan(p_config, treeConfig.rootList[k],
																	  treeConfig.column_distribution[k]);
                    empty_plan->root_task_id = empty_plan->task_id;

                    root_plan_buffer.push(empty_plan);
                }
			} else {
				cout << "ERROR: wrong type, type found = " << treeConfig.type << " File = " << __FILE__
					 << ", Line = " << __LINE__ << endl;
				exit(-1);
			}
		}
    }

	void send_end_plan() {
		plan end_plan(dummy);
		end_plan.message_type = END_PLAN;
        //end_plan.task_id = -1; not used

		ibinstream m1;
		m1 << end_plan;

		for(int i = 0; i < _num_workers; i++) {
			if(i != MASTER_RANK) {
				send_ibinstream(m1, i, PLAN_CHANNEL);
			}
		}

		for(int i = 0; i < _num_workers; i++) {
			if(i != MASTER_RANK) {
				bool local_end; // temp container
				recv_data<bool>(i, STATUS_CHANNEL, local_end);
			}
		}

		global_end_label = true; // terminate all the threads of master, ReqServer, RespQueue, MasterRecver
	}

	void create_col_plans(Task_Master_Col_Split* task_master, column_split_plan* empty_plan,
								   PlanQueue & planQueue) {

		vector<vector<int>> mac2colList;
		mac2colList.resize(_num_workers);

		vector<int> & cols = empty_plan->column_indices;

		set_col_split_assignment(empty_plan->size, mac2colList, cols,
								 empty_plan->parent_slave_id, task_master->task_load);

		for(size_t worker_id = 0; worker_id < mac2colList.size(); worker_id++) {
			if (mac2colList[worker_id].size() > 0) {
				// to track the assigned slave for processing the task
				task_master->slave_ids.push_back(worker_id);
			}
		}

		//create cplans from "mac2colList"
		for(size_t worker_id = 0; worker_id < mac2colList.size(); worker_id++) {
			if(mac2colList[worker_id].size() > 0) {
				column_split_plan* c_plan = new column_split_plan(dummy);

				c_plan->task_id = empty_plan->task_id;
				c_plan->tree_config = empty_plan->tree_config;
				c_plan->level = empty_plan->level;
				c_plan->parent_slave_id = empty_plan->parent_slave_id;
				c_plan->parent_task_id = empty_plan->parent_task_id;
				c_plan->is_left = empty_plan->is_left; // ignore this field for root_plan
				c_plan->size = empty_plan->size;
				c_plan->dest_id = worker_id;
				c_plan->column_indices = mac2colList[worker_id];
				c_plan->total_parent_req = empty_plan->total_parent_req;

				planQueue.add(c_plan);

			}
		}
	}

    // creating or updating task object (for this plan) in master-task-table
    // note that a task may generate may cplans
    // updating means to append slaveID for tracking response progress
	void update_task_master(plan* tree_plan, conmap2t<int, Task_Master*> & task_table,
							PlanQueue & planQueue) { //tree_plan rename to the_plan


        // get the task object from task table, using the plan
		conmap2t_bucket<int, Task_Master*> & bkt = task_table.get_bucket(tree_plan->task_id);

		bkt.lock();

		if(tree_plan->message_type == SUB_TREE_PLAN) {
            //plan's fields has been set outside
            subtree_plan* s_plan = (subtree_plan*) tree_plan;

			Task_Master* task_master = new Task_Master(s_plan->node); //relay the node ref from plan
			task_master->task_id = s_plan->task_id;
            task_master->root_task_id = s_plan->root_task_id;
			task_master->tree_config = s_plan->tree_config;
			task_master->size = s_plan->size;
            task_master->level = s_plan->level;

            if(s_plan->tree_config.sample_col_each_node) {
				s_plan->column_indices.clear();
                random_shuffle(s_plan->tree_config.n_columns, s_plan->column_indices);
            }

			task_master->column_indices = s_plan->column_indices;

            // load balancing inside
            set_subtree_assignment(s_plan, task_master->task_load);

			bkt.insert(s_plan->task_id, task_master);

			bkt.unlock();

			planQueue.add(s_plan);

		} else if (tree_plan->message_type == COL_SPLIT_PLAN) { // COL_SPLIT_PLAN

            column_split_plan* empty_plan = (column_split_plan*) tree_plan;

			Task_Master_Col_Split* task_master_col_split = new Task_Master_Col_Split(empty_plan->node); //to take plan's treeNode ref

			task_master_col_split->task_id = empty_plan->task_id;
            task_master_col_split->root_task_id = empty_plan->root_task_id;
			task_master_col_split->tree_config = empty_plan->tree_config;
			task_master_col_split->size = empty_plan->size;
			task_master_col_split->level = empty_plan->level;
			task_master_col_split->n_met = 0;

			if(empty_plan->tree_config.sample_col_each_node) {
				empty_plan->column_indices.clear();
				random_shuffle(empty_plan->tree_config.n_columns, empty_plan->column_indices);
			}

			task_master_col_split->column_indices = empty_plan->column_indices;

			// inserting to task_table
			bkt.insert(empty_plan->task_id, task_master_col_split);

            bkt.unlock();

            create_col_plans(task_master_col_split, empty_plan, planQueue);

			delete empty_plan;

		} else if (tree_plan->message_type == LABEL_FETCH_PLAN) {
            plan_m2s_label* plan = (plan_m2s_label*) tree_plan;

            Leaf_Task* task = new Leaf_Task(plan->node);
            task->task_id = plan->task_id;
            task->root_task_id = plan->root_task_id;
            task->size = plan->size;
            task->level = plan->level;

            bkt.insert(task->task_id, task);

            bkt.unlock();

            // adding load to load_matrix
            add_leaf_load(task->task_load, plan);

			planQueue.add(plan);
        }
	}

	void run_master(PlanQueue & planQueue, queue<plan*> & root_plan_buffer, deque_p<plan> & plan_buffer,
						 conmap2t<int, Task_Master*> & task_table) {

		plan* temp_plan = plan_buffer.pop_front();

		while (temp_plan != NULL || tree_progress.active_trees() != 0
               || !root_plan_buffer.empty()) {

			if(tree_progress.active_trees() < ACTIVE_TREES_THRESHOLD && !root_plan_buffer.empty()) {

                plan* temp_root_plan = root_plan_buffer.front();
				root_plan_buffer.pop();

				tree_progress.increment(temp_root_plan->root_task_id);

                if(temp_root_plan->size < BFS_PRIORITY_THRESHOLD) {
                    plan_buffer.push_front(temp_root_plan);
                } else {
                    plan_buffer.push_back(temp_root_plan);
                }
			}

#ifdef DEBUG_LOG
            if(tree_progress.active_trees() > ACTIVE_TREES_THRESHOLD) {
                cout << "######################################################"<< endl;
                cout << "[ERROR] active_trees_count = " << tree_progress.active_trees() << endl;
                tree_progress.print_table();
                cout << "######################################################"<< endl;
            }
#endif

#ifdef ASSERT
            assert(tree_progress.active_trees() <= ACTIVE_TREES_THRESHOLD);
#endif

            if(temp_plan == NULL) {
                temp_plan = plan_buffer.pop_front();
            }

			if(temp_plan == NULL) {
                usleep(WAIT_TIME_WHEN_IDLE);
			} else {
                //todo:: make it two conditions:
                //todo:: for SUB_TREE_PLAN, just compute 2-level slave assignments, and then call planQueue.add(temp_plan, dest_ID); (also need to add task)

#ifdef ASSERT
                assert(temp_plan->message_type == SUB_TREE_PLAN
                       || temp_plan->message_type == COL_SPLIT_PLAN
                       || temp_plan->message_type == LABEL_FETCH_PLAN);
#endif

#ifdef DEBUG_LOG
					cout << "[MASTER] main-thread: pop task |" << temp_plan->task_id
						 << "| (p="<<temp_plan->parent_task_id<<") from plan_buffer "
						 << ", level = " << temp_plan->level
                         << ", root_task_id = " << temp_plan->root_task_id
						 << ", File " << __FILE__ << ", Line = " << __LINE__ << endl;

					fout << "[MASTER] main-thread: pop task |" << temp_plan->task_id
						 << "| (p="<<temp_plan->parent_task_id<<") from plan_buffer "
						 << ", level = " << temp_plan->level
						 << ", root_task_id = " << temp_plan->root_task_id
						 << ", File " << __FILE__ << ", Line = " << __LINE__ << endl;
#endif

                // creating or updating task object (for this plan) in master-task-table
                // note that a task may generate may cplans
                // updating means to append slaveID for tracking response progress
                update_task_master(temp_plan, task_table, planQueue);

			}

			temp_plan = plan_buffer.pop_front();
        }

		cout << "MASTER_END, All result found here, Now would send end plan, file = "
			 << __FILE__ << ", Line = " << __LINE__ << endl;

		send_end_plan();

	}

	// for master end ////////////////////////////

	// for slave start ////////////////////////////

	void process_end_plan_slave() {
		bool end_flag  = true;
		ibinstream m;
		m << end_flag;

		send_ibinstream(m, MASTER_RANK, STATUS_CHANNEL);
		global_end_label = true;
	}

    // data has been obtained, now create task into task buffer (for compers to find "best")
	void process_col_split_task(conque_p<Task_Slave> & task_buf, Candidate_Rows & candidate_rows,
								column_split_plan & c_plan) {

		Task_Slave_Col_Split* task = new Task_Slave_Col_Split;
		task->is_root_task = (c_plan.parent_slave_id == -1);
		task->task_id = c_plan.task_id;
		task->tree_config = c_plan.tree_config;
		task->candidate_rows = candidate_rows;
		task->column_indices.swap(c_plan.column_indices);

		task_buf.enqueue(task);
	}

	bool run_slave(ReqQueue & req_queue, conmap2t<int, Task_Slave*> & task_table,
						conque_p<Task_Slave> & task_buf, conque_p<resp2master> & send_buffer) {
		int has_msg;
		MPI_Status status;
		MPI_Iprobe(MASTER_RANK, PLAN_CHANNEL, MPI_COMM_WORLD, &has_msg, &status);

		if(!has_msg) {
			return false; // sleep
		}

		int size;
		MPI_Get_count(&status, MPI_CHAR, &size); // get size of the msg-batch (# of bytes)
		char * buf = new char[size]; //space for receiving this msg-batch, space will be released by obinstream in thread_func(.)
		MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		obinstream m(buf, size);

		char msg_type;

		while (m.end() == false) {
			m >> msg_type;

			if(msg_type == SUB_TREE_PLAN) { // on arrival of sub_tree plan

				subtree_plan p_plan(dummy);
				m >> p_plan;

				Task_Slave_Subtree* subtree_task = new Task_Slave_Subtree;
				subtree_task->task_id = p_plan.task_id;
				subtree_task->tree_config = p_plan.tree_config;

				subtree_task->level = p_plan.level;
				subtree_task->matrix->col.resize(_num_columns);
				subtree_task->mac2colList.swap(p_plan.mac2colList);

				vector<vector<int>> & mac_cols = subtree_task->mac2colList; // destination map (index/mac_id, list_of_col_indices)

				int num_slaves = 0;
                vector<int> & columns_to_consider = subtree_task->column_indices;

				for(size_t i = 0; i < mac_cols.size(); i++) {
					if(mac_cols[i].size() > 0) {
						num_slaves++;
                        for(size_t col_idx = 0; col_idx < mac_cols[i].size(); col_idx++) {
                            columns_to_consider.push_back(mac_cols[i][col_idx]);
                        }
                    }
				}

				subtree_task->n_req = num_slaves + 1; // extra 1 for parent row_indices

				conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(p_plan.task_id);
				bucket.lock();

				unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();

#ifdef ASSERT
				auto it = kvmap.find(p_plan.task_id);
				assert(it == kvmap.end());
#endif

				bucket.insert(p_plan.task_id, subtree_task);

				if (p_plan.parent_slave_id != -1 && p_plan.parent_slave_id != _my_rank) {
					req_s2s_row_indices* req = new req_s2s_row_indices;
					req->task_id = p_plan.task_id;
					req->parent_task_id = p_plan.parent_task_id;
					req->is_left = p_plan.is_left;
					req->parent_req = p_plan.total_parent_req;

					if(mac_cols[_my_rank].size() > 0) {
						req->n_local_cols += mac_cols[_my_rank].size() + 1;
					} else { // extra request of parent_sub_tree_task
						req->n_local_cols += 1;
					}

					req_queue.add(req, p_plan.parent_slave_id);

					bucket.unlock();

				} else { // (p_plan.parent_slave_id ==  _my_rank || p_plan.parent_slave_id == -1)

					Candidate_Rows & parent_rows = subtree_task->candidate_rows;

					if(p_plan.parent_slave_id == -1) {

						for(size_t i = 0; i < _n_samples; i++) {
							parent_rows.indexes.push_back(i);
						}

					} else { // p_plan.parent_slave_id ==  _my_rank

						conmap2t_bucket<int, Task_Slave*> & bkt = task_table.get_bucket(p_plan.parent_task_id);

						if(!(task_table.bucket_hash(p_plan.task_id) == task_table.bucket_hash(p_plan.parent_task_id))) {
						    bkt.lock();
						}

						unordered_map<int, Task_Slave*> & kvmap = bkt.get_map();

						auto it2 = kvmap.find(p_plan.parent_task_id);

#ifdef ASSERT
                        assert(it2 != kvmap.end());
#endif

						Task_Slave_Col_Split* parent_task = (Task_Slave_Col_Split*) it2->second;
						parent_task->parent_request = p_plan.total_parent_req;

						if(p_plan.is_left) {
							parent_rows = parent_task->best_left_rows;
						} else {
							parent_rows = parent_task->best_right_rows;
						}

						if(parent_rows.stop_splitting) { // KS ->master, no need to fetch remote columns
						    for(size_t count = 0; count < mac_cols.size(); count++) {
                                parent_task->n_cols += mac_cols[count].size();
                            }

						    parent_task->n_cols += 1;

						} else {
                            if(mac_cols[_my_rank].size() > 0) {
                                parent_task->n_cols += mac_cols[_my_rank].size() + 1;
                            } else {
                                parent_task->n_cols += 1;
                            }
                        }

						if(parent_task->n_cols == parent_task->parent_request) {
                            if(parent_task->best_column_idx >= 0) {
                                delete_best_split(parent_task->best_column_idx, parent_task->best_split); //garbage collect parent task's best-split
                            }

                            delete parent_task;
                            kvmap.erase(it2);
                        }

                        if(!(task_table.bucket_hash(p_plan.task_id) == task_table.bucket_hash(p_plan.parent_task_id))) {
                            bkt.unlock();
                        }

					}

					subtree_task->n_met++; // extra one request for parent_sub_tree_task rows

					if(mac_cols[_my_rank].size() > 0) {
						for(size_t j = 0; j < mac_cols[_my_rank].size(); j++) {
							Column* column = getColumn(mac_cols[_my_rank][j], parent_rows.indexes);
							subtree_task->matrix->col[mac_cols[_my_rank][j]] = column;
						}

						subtree_task->n_met++;
					}

					if((subtree_task->n_met == subtree_task->n_req) || parent_rows.stop_splitting) {
						auto it2 = kvmap.find(subtree_task->task_id);
						kvmap.erase(it2);
						task_buf.enqueue(subtree_task);
                        bucket.unlock();
                        continue;
					}

					bucket.unlock();
				}

				// "main_thread_slave" sending the data fetch request for sub_tree_task
				for(size_t i = 0; i < mac_cols.size(); i++) {

					if(i != _my_rank && mac_cols[i].size() > 0) {

						req_s2s_columns* req = new req_s2s_columns;
						req->task_id = p_plan.task_id;
						req->column_indices = mac_cols[i];
						req->parent_slave_id = p_plan.parent_slave_id;
						req->parent_task_id = p_plan.parent_task_id;
						req->is_left = p_plan.is_left;
						req->parent_req = p_plan.total_parent_req;

						req_queue.add(req, i);

					}

				}

			} else if (msg_type == COL_SPLIT_PLAN) { // on arrival of col_split_plan

				column_split_plan c_plan(dummy);
				m >> c_plan;

				if(c_plan.parent_slave_id == _my_rank) { // parent_slave is local, no need to rowID fetching

                    //this is a slave that master gives a few columns (rows to be obtained from parent task)

                    //get rowIDs from parent task (which is in this slave!!!)
					conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(c_plan.parent_task_id);
					bucket.lock();

					unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();

					auto it = kvmap.find(c_plan.parent_task_id);

#ifdef ASSERT
                    if(it == kvmap.end()) {
						cout << "parent task id = |" << c_plan.parent_task_id << "| not found, task_id = "
							 << c_plan.task_id << ", rank = " << _my_rank << " File = " << __FILE__
							 << ", Line = " << __LINE__ << endl;
					}

                    assert(it != kvmap.end());
#endif

					Task_Slave_Col_Split* parent_task = (Task_Slave_Col_Split*) it->second;
					parent_task->parent_request = c_plan.total_parent_req;

                    // from parent task (myself), get the child rowIDs
					Candidate_Rows rows;

					if(c_plan.is_left) {
						rows = parent_task->best_left_rows;
					} else {
						rows = parent_task->best_right_rows;
					}

                    //parent task updates progress
					parent_task->n_cols += c_plan.column_indices.size();

					if(parent_task->n_cols == parent_task->parent_request) {
						kvmap.erase(it);

						if(parent_task->best_column_idx >= 0) {
							delete_best_split(parent_task->best_column_idx, parent_task->best_split);
						}

						delete parent_task;
					}

					bucket.unlock();

                    // local data has been obtained, now create task into task buffer (for compers to find "best")
                    process_col_split_task(task_buf, rows, c_plan);


                } else if (c_plan.parent_slave_id == -1) { // entire column
                    //root node, no parent task

					Candidate_Rows candidate_rows;
					vector<size_t> & row_indices = candidate_rows.indexes;
					for(size_t i = 0; i < _n_samples; i++) {
						row_indices.push_back(i);
					}

					process_col_split_task(task_buf, candidate_rows, c_plan);

				} else { // fetch row_indices from remote, create entry task_table

					conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(c_plan.task_id);
					bucket.lock();

					Task_Slave_Col_Split* task = new Task_Slave_Col_Split;
					task->task_id = c_plan.task_id;
					task->tree_config = c_plan.tree_config;
					task->column_indices.swap(c_plan.column_indices);

					bucket.insert(task->task_id, task);

					bucket.unlock();

					req_s2s_row_indices* req = new req_s2s_row_indices;
					req->task_id = task->task_id;
					req->parent_task_id = c_plan.parent_task_id;
					req->is_left = c_plan.is_left;
					req->n_local_cols += task->column_indices.size();
					req->parent_req = c_plan.total_parent_req;

					req_queue.add(req, c_plan.parent_slave_id);

				}

			} else if (msg_type == TASK_DELETE_PLAN) {

				plan delete_plan(dummy);
				delete_plan.message_type = TASK_DELETE_PLAN; //todo:: may delete
				m >> delete_plan; //just to get the task ID

				conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(delete_plan.task_id);
				bucket.lock();

				unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();
				auto it = kvmap.find(delete_plan.task_id);

#ifdef ASSERT
                assert(it != kvmap.end());
#endif

				Task_Slave_Col_Split* task_slave_col_split = (Task_Slave_Col_Split*) it->second;

				//do the deletion
				kvmap.erase(it); //delete entry
				bucket.unlock();

				if(task_slave_col_split->best_column_idx >= 0) {
					delete_best_split(task_slave_col_split->best_column_idx, task_slave_col_split->best_split);
				}

				delete task_slave_col_split;//delete task

			} else if (msg_type == END_PLAN) {

				plan p_plan(dummy);
				p_plan.message_type = END_PLAN;
				m >> p_plan;

				process_end_plan_slave();

			} else if (msg_type == LABEL_FETCH_PLAN) {

                plan_m2s_label plan(dummy);
                m >> plan;

                conmap2t_bucket<int, Task_Slave*> & bucket = task_table.get_bucket(plan.parent_task_id);
                bucket.lock();

                unordered_map<int, Task_Slave*> & kvmap = bucket.get_map();
                auto it = kvmap.find(plan.parent_task_id);

#ifdef ASSERT
                assert(it != kvmap.end());
#endif

                Task_Slave_Col_Split* parent_task = (Task_Slave_Col_Split*) it->second;
				parent_task->parent_request = plan.total_parent_req;

                string label;
                TreeConfig & treeConfig = parent_task->tree_config;

                if (plan.is_left) {
                    vector<size_t> & rows = parent_task->best_left_rows.indexes;
                    set_node_label(cserver.X, rows, treeConfig, label);
                } else {
                    vector<size_t> & rows = parent_task->best_right_rows.indexes;
                    set_node_label(cserver.X, rows, treeConfig, label);
                }

				parent_task->n_cols += plan.n_local_cols;

                if(parent_task->n_cols == parent_task->parent_request) {

                    if(parent_task->best_column_idx >= 0) {
                        delete_best_split(parent_task->best_column_idx, parent_task->best_split);
                    }

                    delete parent_task;

                    kvmap.erase(it);
                }

                bucket.unlock();

                resp_s2m_label* resp = new resp_s2m_label;
                resp->task_id = plan.task_id;
                resp->node_label = label;

                send_buffer.enqueue(resp);

            } else {
				cout << "ERROR: File = " << __FILE__ << ", Line = " << __LINE__ << endl;
				exit(-1);
			}
		}

		return true;
	}

	//program entry point
    void run(const WorkerParams& params)
    {
    	assert(params.replicate < _num_workers); //otherwise, a machine will keep a column more than once
    	//note that only n-1 slaves are partitioning the data

		vector<TreeConfig> configList; //only used my master
		load_job_config(params.job_file_path.c_str()); //every worker loads job-config
		if(_my_rank == MASTER_RANK) load_tree_configs(params, configList);

        //========== load metadata
        cserver.load_meta(meta_file.c_str()); // HDFS path
        _num_columns = cserver.X.col.size();

		if(_my_rank == MASTER_RANK) {

			// get all non-y column IDs into "all_columns"
            for(int j = 0; j < _num_columns; j++) {
                if(j == y_index) {
                    continue;
                }
                all_columns.push_back(j);
            }

            // initialize columns of tree-configs
			for(size_t i = 0; i < configList.size(); i++) {
				TreeConfig & config = configList[i];

				if(config.type == DECISION_TREE) {
					config.column_distribution.resize(1);
					config.rootList.resize(1);

					config.column_distribution[0] = all_columns;

				} else if (config.type == RANDOM_FOREST) {

                    config.column_distribution.resize(config.num_trees);
                    config.rootList.resize(config.num_trees);

                    int num_columns = (_num_columns - 1) * config.column_sample; //number of sampled columns

                    for(int tree_index = 0; tree_index < config.num_trees; tree_index++) {
                    	//config.column_distribution[tree_index].clear();
                        random_shuffle(num_columns, config.column_distribution[tree_index]);
                    }

				} else {
					cout << "ERROR : Wrong type, type = " << config.type << ", File " << __FILE__
						 << ", Line " << __LINE__ << endl;
					exit(-1);
				}
			}

		}

		//========== check path + init
        if (_my_rank == MASTER_RANK)
        {
        	if (dirCheck(meta_file.c_str()) == -1) return;
            if (dirCheck(train_file.c_str()) == -1) return;
        }

        //========== build mac_map
        if(_my_rank == MASTER_RANK)
        {
        	mac_map.resize(_num_columns);
        	if(params.column_assignment == MODE_ALL_COLUMNS)
        	{
        		for(int i=0; i<_num_columns; i++) {
					for(int j=0; j<_num_workers; j++) {
						if(j == MASTER_RANK) {
							continue;
						}
						mac_map[i].push_back(j);
					}
				}
        	}
        	else if(params.column_assignment == MODE_REPLICATE)
        	{
        		//init min-heap
        		priority_queue<min_heap_entry> min_heap;
        		for(int i=0; i<_num_workers; i++)
				{
					if(i == MASTER_RANK) { //master is not serving data, to save bandwidth for control msgs
						continue;
					}

        			min_heap_entry en;
					en.count = 0;
					en.rank = i;
					min_heap.push(en);
				}
        		//assign:
        		//algo: each time assign k replicates of a column to lightest-loaded machines
        		for(int i=0; i<_num_columns; i++)
        		{
        			if(i == y_index)
        			{//y
        				for(int j=0; j<_num_workers; j++) {
							if(j == MASTER_RANK) {
								continue;
							}
                            mac_map[i].push_back(j);
                        }
        			}
        			else
        			{//Xi
        				vector<min_heap_entry> mins;
						for(int j=0; j<params.replicate; j++)
						{
							mins.push_back(min_heap.top());
							min_heap.pop();
						}
						for(size_t j=0; j<mins.size(); j++)
						{
							mac_map[i].push_back(mins[j].rank);
							mins[j].count++;
							min_heap.push(mins[j]);
						}
        			}
        		}
        	}
        	else
        	{
        		cout<< "ERROR:: File = " << __FILE__ << ", Line = " << __LINE__
        				<< ": this column-assignment mode is not implemented yet" << endl;
        		exit(-1);
        	}
        	//bcast-send
        	masterBcast(mac_map); // blocking operation
        }
        else
        {
        	slaveBcast(mac_map); // blocking operation
        }

        //========== load columns
		if(_my_rank != MASTER_RANK) {
			for(size_t i=0; i<mac_map.size(); i++)
			{
				vector<int> & macs = mac_map[i];
				auto it = find(macs.begin(), macs.end(), _my_rank);
				if(it == macs.end()) { // can happen, column 'i' does not need to be loaded in _my_rank
					continue;
				}
				cserver.load_column(train_file.c_str(), i);
			}
		}

		if(_my_rank == MASTER_RANK) {
			deque_p<plan> plan_buffer;
            queue<plan*> root_plan_buffer;
			conmap2t<int, Task_Master*> task_table; // only for master
			PlanQueue planQueue;
			MasterRecver master_receiver(task_table, planQueue, plan_buffer);

			if(_n_samples < subtree_D) {
				create_root_subtree_plans(root_plan_buffer, configList);
			} else {
				create_root_col_split_plan(root_plan_buffer, configList);
			}

			auto start = chrono::system_clock::now();
			run_master(planQueue, root_plan_buffer, plan_buffer, task_table);
			auto end = chrono::system_clock::now();
			auto train_time = chrono::duration_cast<chrono::milliseconds>(end - start).count();

			cout << "Training done... Now printing and deleting the tree, File = " << __FILE__ << ", Line = "
				 << __LINE__ << endl;

			//load test data (if specified)
            if(!params.test_file_path.empty()) {
                // loading test_set.csv and test_meta.csv
                load_csv(params.test_file_path.c_str(), params.test_meta_file.c_str(), test_set);
            }

            // print tree configs + accuracy results
			for(size_t i = 0; i < configList.size(); i++) {
				TreeConfig & treeConfig = configList[i];
                cout << i <<": " << treeConfig.type << " (testing and) saving ... " << endl;

                if(!params.test_file_path.empty()) {
                    if(treeConfig.IMPURITY_FUNC == IMPURITY_VARIANCE) {
                        double RMSE = calculate_RMSE(treeConfig.rootList, test_set, y_index, INT_MAX);
                        cout << "RMSE found for Task " << i  << " is : "
                             << RMSE << endl << endl;
                    } else if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY
                               || treeConfig.IMPURITY_FUNC == IMPURITY_GINI
                               || treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR){
                        cout << "Accuracy found for Task " << i  << " is : "
                             << get_accuracy(test_set, treeConfig.rootList, y_index, INT_MAX) << endl << endl;
                    } else {
                        cout << "ERROR : File = " << __FILE__ << ", Line = " << __LINE__ << endl;
                        exit(-1);
                    }
                }

				vector<TreeNode*> & rootList = treeConfig.rootList;

                string tree_dir = job_dir;

                if(SAVE_TREE) {
                    string tree_path = tree_dir + "/model_" + to_string(i);
					save_forest(tree_path.c_str(), rootList);
                }

				vector<vector<int>> & column_distribution = treeConfig.column_distribution;

				for(size_t root_idx = 0; root_idx < rootList.size(); root_idx++) {

					if(PRINT_TREE) {
						cout << "##################################################" << endl;
						cout << "config : " << treeConfig.type << ", tree index = " << root_idx << endl;

						cout << "##################################################" << endl;
						cout << "Candidate Columns = ";
						for(size_t columns_index = 0; columns_index < column_distribution[root_idx].size(); columns_index++) {
							cout << column_distribution[root_idx][columns_index] << ", ";
						}
						cout << endl;

						cout << "##################################################" << endl;
						print_tree(rootList[root_idx]);
						cout << endl << endl;
					}

#ifdef DEBUG_LOG
					cout << "trying to delete index = " << root_idx << endl;
#endif
					free_tree(rootList[root_idx], y_index, cserver.X);
#ifdef DEBUG_LOG
					cout << "deleted tree index = " << root_idx << endl;
#endif
				}
			}

			cout << "...................................................................." << endl;
			cout << "total train time = " << train_time << " milliseconds " << endl;

		} else { // slaves
			conque_p<Task_Slave> task_buf; // use by slave compers
			conque_p<resp2master> send_buffer;

			conmap2t<int, Task_Slave*> task_table;
			RespQueue resp_queue;
			ReqQueue req_queue;
			ReqServer req_server(task_table, resp_queue, req_queue); // run on master, slave, serving remote data for sub_tree, row_indices (subtree and col_split)

			SlaveSender slaveSender(send_buffer);
			RespServer resp_server(task_table, task_buf, resp_queue);
            // would run only on slaves
            // since RespServer is the recv-side of data serving

			vector<Comper*> compers(num_compers);
			for(int i = 0; i < num_compers; i++) {
				compers[i] = new Comper(task_buf, send_buffer, task_table);
			}

			while (global_end_label == false) {
				bool computing = run_slave(req_queue, task_table, task_buf, send_buffer);
				if(!computing) {
					usleep(STATUS_SYNC_TIME_GAP);
				}
			}

			for(size_t i = 0; i < compers.size(); i++) {
				delete compers[i];
			}

		}

#ifdef ASSERT
        if (_my_rank == MASTER_RANK) {
            bool isOk = true;

            for(size_t i = 0; i < load_matrix.size(); i++) {
                for(size_t j = 0; j < load_matrix[i].size(); j++) {
                    isOk = isOk && (load_matrix[i][j] == 0);
                }
            }

            if(!isOk) {
                for(size_t i = 0; i < load_matrix.size(); i++) {
                    for(size_t j = 0; j < load_matrix[i].size(); j++) {
                        cout << load_matrix[i][j] << ", ";
                    }
                    cout << endl;
                }
            }

            assert(isOk); // todo

            assert(tree_progress.active_trees() == 0);
        }
#endif
        //cout << "DEBUG: Terminating worker " << _my_rank << endl;
	}
};

#endif
