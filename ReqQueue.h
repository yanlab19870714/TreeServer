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

#ifndef REQQUEUE_H_
#define REQQUEUE_H_

#include "tree_msgs.h"
#include "conque_p.h"
#include <atomic>
#include <thread>
#include <unistd.h> //for usleep()

using namespace std;

// for data serving
// - REQ_ROW_FETCH
// - REQ_COL_FETCH

// This is to run on slaves
class ReqQueue {
public:
    typedef conque_p<req_s2s> Buffer;
    typedef vector<Buffer> Queue;

    Queue q;
    thread main_thread;

    void get_msgs(int i, ibinstream & m)
    {
        Buffer & buf = q[i];
        req_s2s* temp = buf.dequeue(); //fetch a request to send
        while(temp != NULL){ //fetching till reach list-head
            if(temp->message_type == REQ_COL_FETCH) {
                // this req is inserted by main-threads of key-task (who builds subtree) toward a slave that serves column values
                // only used by subtree task
                req_s2s_columns* req = (req_s2s_columns*) temp;
                m << *req;
                delete req;
            }
            else if(temp->message_type == REQ_ROW_FETCH) {
                // this req is inserted by ReqServer toward parent-task
                req_s2s_row_indices* req = (req_s2s_row_indices*) temp;
                m << *req;
                delete req;
            }
            else
            {
                cout << "[Error]:: File " << __FILE__ << " Line = " << __LINE__ << endl;
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
        thread t(&ReqQueue::get_msgs, this, i, ref(*m0)); //assisting thread

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
                t = thread(&ReqQueue::get_msgs, this, j, ref(*m1));
                if(m0->size() > 0)
                {
                    sth_sent = true;
                    //send reqs to tgt
                    MPI_Send(m0->get_buf(), m0->size(), MPI_CHAR, i, REQ_CHANNEL, MPI_COMM_WORLD);
                    //------
                    delete m0;
                    m0 = new ibinstream;
                }
                use_m0 = false;
            }
            else //odd
            {
                //use m1, set m0
                t = thread(&ReqQueue::get_msgs, this, j, ref(*m0));
                if(m1->size() > 0)
                {
                    sth_sent = true;

                    MPI_Send(m1->get_buf(), m1->size(), MPI_CHAR, i, REQ_CHANNEL, MPI_COMM_WORLD);
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

    ReqQueue()
    {
        q.resize(_num_workers);
        main_thread = thread(&ReqQueue::thread_func, this);
    }

    ~ReqQueue()
    {
        main_thread.join();
    }

    void add(req_s2s* req, int tgt)
    {
        Buffer & buf = q[tgt];
        buf.enqueue(req);
    }
};

#endif