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

#ifndef PLANQUEUE_H_
#define PLANQUEUE_H_

// this is the plan sender on the master

#include "../tree_msgs.h"
#include "../conque_p.h"

#include <atomic>
#include <thread>
#include <unistd.h> //for usleep()

using namespace std;

// would be used only in master
class PlanQueue {
public:
    typedef conque_p<plan> Buffer;
    typedef vector<Buffer> Queue;

    Queue q;
    thread main_thread;

    void get_msgs(int i, ibinstream & m)
    {
        if(i == _my_rank) return;

        Buffer & buf = q[i];
        plan* temp = buf.dequeue(); //for fetching KeyT items
        while(temp != NULL){ //fetching till reach list-head
            if(temp->message_type == SUB_TREE_PLAN) {
                subtree_plan* s_plan = (subtree_plan*) temp;
                m << *s_plan;
                delete s_plan; //serialized for sending, plan object should be deleted
            } else if (temp->message_type == COL_SPLIT_PLAN) {
                column_split_plan* c_plan = (column_split_plan*) temp;
                m << *c_plan;
                delete c_plan; //serialized for sending, plan object should be deleted
            } else if (temp->message_type == TASK_DELETE_PLAN || temp->message_type == END_PLAN) {
                m << *temp;
                delete temp; //serialized for sending, plan object should be deleted
            } else if (temp->message_type == LABEL_FETCH_PLAN) {
                plan_m2s_label* plan = (plan_m2s_label*) temp;
                m << *plan;

                delete plan;
            } else {
                cout << "ERROR : File " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }

            temp = buf.dequeue();
        }
    }

    void thread_func() //managing requests to tgt_rank
    {
        int i = 0; //target worker to send
        bool sth_sent = false; //if sth is sent in one round, set it as true
        ibinstream* m0 = new ibinstream;
        ibinstream* m1 = new ibinstream;
        //m0 and m1 are alternating
        thread t(&PlanQueue::get_msgs, this, i, ref(*m0)); //assisting thread

        bool use_m0 = true; //tag for alternating
        /*clock_t last_tick = clock();*/
        while(global_end_label == false) //otherwise, thread terminates
        {
            t.join(); //m0 or m1 becomes ready to send
            int j = i+1;
            if(j == _num_workers) j = 0;
            if(use_m0) //even
            {
                //use m0, set m1
                t = thread(&PlanQueue::get_msgs, this, j, ref(*m1));
                if(m0->size() > 0)
                {
                    sth_sent = true;
                    //send reqs to tgt
                    MPI_Send(m0->get_buf(), m0->size(), MPI_CHAR, i, PLAN_CHANNEL, MPI_COMM_WORLD);
                    //------
                    delete m0;
                    m0 = new ibinstream;
                }
                use_m0 = false;
            }
            else //odd
            {
                //use m1, set m0
                t = thread(&PlanQueue::get_msgs, this, j, ref(*m0));
                if(m1->size() > 0)
                {
                    sth_sent = true;

                    MPI_Send(m1->get_buf(), m1->size(), MPI_CHAR, i, PLAN_CHANNEL, MPI_COMM_WORLD);
                    //------
                    delete m1;
                    m1 = new ibinstream;
                }
                use_m0 = true;
            }
            //------------------------
            i = j;
            if(j == 0)
            {
                if(!sth_sent) {
                    usleep(WAIT_TIME_WHEN_IDLE);
                } else {
                    sth_sent = false;
                }
            }
        }
        t.join();
        delete m0;
        delete m1;
    }



    PlanQueue()
    {
        q.resize(_num_workers);
        main_thread = thread(&PlanQueue::thread_func, this);
    }

    ~PlanQueue()
    {
        main_thread.join();
    }

    void add(plan* p_plan)
    {
#ifdef DEBUG_LOG
        if(p_plan->dest_id == MASTER_RANK) {
            cout << "[PlanQueue] task id = |" << p_plan->task_id << "|, plan_type = "
                 << p_plan->message_type << ", File = " << __FILE__ << ", Line = "
                 << __LINE__ << endl;
        }
#endif
#ifdef ASSERT
        assert(p_plan->dest_id != MASTER_RANK);
#endif
        Buffer & buf = q[p_plan->dest_id];
        buf.enqueue(p_plan);
    }
};

#endif