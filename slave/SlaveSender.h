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

#ifndef SLAVESENDER_H_
#define SLAVESENDER_H_

#include <unistd.h> //for usleep()
#include <thread>
#include "../conque_p.h"
#include "../tree_msgs.h"

using namespace std;

// run only on slaves
class SlaveSender {
public:
    thread main_thread;
    conque_p<resp2master> & ready_buffer;

    void get_msgs(ibinstream & m) {
        resp2master* temp = ready_buffer.dequeue();

        while (temp != NULL) {

            if (temp->message_type == RESP_LABEL_FETCH) {
                resp_s2m_label* resp = (resp_s2m_label*) temp;
                m << *resp;

                delete resp;

            } else if(temp->message_type == SUB_TREE_RESP) {
                subtree_resp* resp = (subtree_resp*) temp;
                m << *resp;
                
                free_tree(resp->root, y_index, cserver.X);

                delete resp;
            } else if (temp->message_type == COL_SPLIT_RESP) {

                if(temp->best_column_index == -2) { // column_split_resp and can not split on assigned columns
                    m << *temp;
                    delete temp;
                } else {
                    Column* column;

                    //get column
                    if(temp->best_column_index == -1) {
                        column = cserver.X.get_column(y_index);
                    } else {
                        column = cserver.X.get_column(temp->best_column_index);
                    }

                    //serialize column for sending, and delete resp object
                    //no need to delete best-split, as it is managed by task object
                    if(column->data_type == ELEM_BOOL) {
                        column_split_resp<bool>* resp = (column_split_resp<bool>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_SHORT) {
                        column_split_resp<short>* resp = (column_split_resp<short>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_INT) {
                        column_split_resp<int>* resp = (column_split_resp<int>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_FLOAT) {
                        column_split_resp<float>* resp = (column_split_resp<float>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_DOUBLE) {
                        column_split_resp<double>* resp = (column_split_resp<double>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_CHAR) {
                        column_split_resp<char>* resp = (column_split_resp<char>*) temp;
                        m << *resp;

                        delete resp;
                    } else if(column->data_type == ELEM_STRING) {
                        column_split_resp<string>* resp = (column_split_resp<string>*) temp;
                        m << *resp;

                        delete resp;
                    } else {
                        cout << "Error : File " << __FILE__ << ", Line = " << __LINE__ << endl;
                        exit(-1);
                    }
                }

            } else {
                cout << "Error : File " << __FILE__ << ", Line = " << __LINE__ << endl;
                exit(-1);
            }

            temp = ready_buffer.dequeue();
        }
    }

    void run() {
        bool sth_sent = false;
        ibinstream* m0 = new ibinstream;
        ibinstream* m1 = new ibinstream;
        bool use_m0 = true;

        thread t = thread(&SlaveSender::get_msgs, this, ref(*m0));

        while (global_end_label == false) {
            t.join();

            if(use_m0) {
                t = thread(&SlaveSender::get_msgs, this, ref(*m1));

                if(m0->size() > 0) {
                    send_ibinstream(*m0, MASTER_RANK, PLAN_CHANNEL);
                    sth_sent = true;
                    delete m0;
                    m0 = new ibinstream;
                }
                use_m0 = false;
            } else {
                t = thread(&SlaveSender::get_msgs, this, ref(*m0));

                if(m1->size() > 0) {
                    send_ibinstream(*m1, MASTER_RANK, PLAN_CHANNEL);
                    sth_sent = true;
                    delete m1;
                    m1 = new ibinstream;
                }
                use_m0 = true;
            }

            if(sth_sent == false) {
                usleep(WAIT_TIME_WHEN_IDLE);
            } else {
                sth_sent = false;
            }
        }

        t.join();

        delete m0;
        delete m1;
    }

    SlaveSender(conque_p<resp2master> & r_buffer) : ready_buffer(r_buffer) {
        main_thread = thread(&SlaveSender::run, this);
    }

    ~SlaveSender() {
        main_thread.join();
    }
};

#endif