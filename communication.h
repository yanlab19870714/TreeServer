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

//Acknowledgements: this code is implemented by referencing pregel-mpi (https://code.google.com/p/pregel-mpi/) by Chuntao Hong.

#ifndef COMMUNICATION_H
#define COMMUNICATION_H

#include <mpi.h>
#include "timer.h"
#include "serialization.h"

//============================================
//binstream-level send/recv
void send_ibinstream(ibinstream& m, int dst, int tag)
{
    MPI_Send(m.get_buf(), m.size(), MPI_CHAR, dst, tag, MPI_COMM_WORLD);
}

obinstream recv_obinstream(int src, int tag)
{
	MPI_Status status;
	MPI_Probe(src, tag, MPI_COMM_WORLD, &status);
	int size;
	MPI_Get_count(&status, MPI_CHAR, &size); // get size of the msg-batch (# of bytes)
	char * buf = new char[size];
	MPI_Recv(buf, size, MPI_CHAR, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return obinstream(buf, size);
}

//============================================
//obj-level send/recv
template <class T>
void send_data(const T& data, int dst, int tag)
{
    ibinstream m;
    m << data;
    send_ibinstream(m, dst, tag);
}

template <class T>
T recv_data(int src, int tag)
{
    obinstream um = recv_obinstream(src, tag);
    T data;
    um >> data;
    return data;
}

template <class T>
void recv_data(int src, int tag, T& data)
{
    obinstream um = recv_obinstream(src, tag);
    um >> data;
}

//============================================
//all-to-all

template <class T>
void all_to_all(vector<vector<T*> > & to_exchange, int tag)
{
    StartTimer(COMMUNICATION_TIMER);
    int np = get_num_workers();
    int me = get_worker_id();
    for (int i = 0; i < np; i++)
    {
        int partner = (i - me + np) % np;
        if (me != partner)
        {
            if (me < partner)
            {
                StartTimer(SERIALIZATION_TIMER);
                //send
                ibinstream * m = new ibinstream;
                *m << to_exchange[partner];
                for(int k = 0; k < to_exchange[partner].size();k++)
                    delete to_exchange[partner][k];
                vector<T*>().swap(to_exchange[partner]);
                StopTimer(SERIALIZATION_TIMER);
                StartTimer(TRANSFER_TIMER);
                send_ibinstream(*m, partner, tag);
                delete m;
                StopTimer(TRANSFER_TIMER);
                //receive
                StartTimer(TRANSFER_TIMER);
                obinstream um = recv_obinstream(partner, tag);
                StopTimer(TRANSFER_TIMER);
                StartTimer(SERIALIZATION_TIMER);
                um >> to_exchange[partner];
                StopTimer(SERIALIZATION_TIMER);
            }
            else
            {
                StartTimer(TRANSFER_TIMER);
                //receive
                obinstream um = recv_obinstream(partner, tag);
                StopTimer(TRANSFER_TIMER);
                StartTimer(SERIALIZATION_TIMER);
                //send
                ibinstream * m = new ibinstream;
                *m << to_exchange[partner];
                for(int k = 0; k < to_exchange[partner].size();k++)
                    delete to_exchange[partner][k];
                vector<T*>().swap(to_exchange[partner]);
                StopTimer(SERIALIZATION_TIMER);
                StartTimer(TRANSFER_TIMER);
                send_ibinstream(*m, partner, tag);
                delete m;
                StopTimer(TRANSFER_TIMER);
                um >> to_exchange[partner];
            }
        }
    }
    StopTimer(COMMUNICATION_TIMER);
}

//================================================================
//bcast
template <class T>
void masterBcast(T& to_send)
{ //broadcast
    ibinstream m;
    m << to_send;
    int size = m.size();

    MPI_Bcast(&size, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

    char* sendbuf = m.get_buf();
    MPI_Bcast(sendbuf, size, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);
}

template <class T>
void slaveBcast(T& to_get)
{ //broadcast
    int size;

    MPI_Bcast(&size, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

    char* recvbuf = new char[size];
    MPI_Bcast(recvbuf, size, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

    obinstream um(recvbuf, size);
    um >> to_get;
}

#endif
