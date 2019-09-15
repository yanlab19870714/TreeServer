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

#ifndef RESPQUEUE_H_
#define RESPQUEUE_H_

#include "conque_p.h"
#include "tree_msgs.h"
#include <atomic>
#include <thread>
#include <unistd.h> //for usleep()
using namespace std;


// T should be of type resp_msg_s2s
class RespQueue { //part of ReqServer
public:
    typedef conque_p<resp_s2s> Buffer;
    typedef vector<Buffer> Queue;

    Queue q;
    thread main_thread;

    void get_msgs(int i, ibinstream & m)
    {
        Buffer & buf = q[i];

        resp_s2s* temp = buf.dequeue(); //for fetching VertexT items

        while(temp != NULL) //fetching till reach list-head
        {
            if(temp->message_type == RESP_COL_FETCH) {
                //inserted by RespServer
                //- rowIDs is received from parent-task
                //- this response is processed by RespServer to insert fetched data to this queue
                resp_s2s_columns* resp = (resp_s2s_columns*) temp;
                m << *resp;

                //now that it's given to the send queue
                //we delete the resp object
                for(size_t i = 0; i < resp->cols.size(); i++) {
                    delete resp->cols[i];
                }

                delete resp;

            }
            else if(temp->message_type == RESP_ROW_FETCH)
            {
                //subtree column fetch
                //- this is the parent-slave, send back row indices
                //- target is to the slave serving column data
                resp_s2s_row_indices* resp = (resp_s2s_row_indices*) temp;

                m << *resp;

                delete resp;
            } else if (temp->message_type == RESP_LABEL_FETCH) {
                resp_s2m_label* resp = (resp_s2m_label*) temp;
                m << *resp;

                delete resp;
            } else {
                cout << "ERROR : File = " << __FILE__ << ", Line = " << __LINE__ << endl;
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
        thread t(&RespQueue::get_msgs, this, 0, ref(*m0)); //assisting thread
        bool use_m0 = true; //tag for alternating
        while(global_end_label == false) //otherwise, thread terminates
        {
            t.join();//m0 or m1 becomes ready to send
            int j = i+1;
            if(j == _num_workers) j = 0;
            if(use_m0) //even
            {
                //use m0, set m1
                t = thread(&RespQueue::get_msgs, this, j, ref(*m1));
                if(m0->size() > 0)
                {
                    sth_sent = true;
                    //send reqs to tgt
                    MPI_Send(m0->get_buf(), m0->size(), MPI_CHAR, i, RESP_CHANNEL, MPI_COMM_WORLD);
                    delete m0;
                    m0 = new ibinstream;
                }
                use_m0 = false;
            }
            else
            {
                //use m1, set m0
                t = thread(&RespQueue::get_msgs, this, j, ref(*m0));
                if(m1->size() > 0)
                {
                    sth_sent = true;
                    //send reqs to tgt
                    MPI_Send(m1->get_buf(), m1->size(), MPI_CHAR, i, RESP_CHANNEL, MPI_COMM_WORLD);
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
                if(!sth_sent) usleep(WAIT_TIME_WHEN_IDLE);
                else sth_sent = false;
            }
        }
        t.join();

        delete m0;
        delete m1;
    }

    RespQueue()
    {
        q.resize(_num_workers);
        main_thread = thread(&RespQueue::thread_func, this);
    }

    ~RespQueue()
    {
        main_thread.join();
    }

    void add(resp_s2s* resp, int tgt)
    {
        Buffer & buf = q[tgt];
        buf.enqueue(resp);
    }
};

#endif