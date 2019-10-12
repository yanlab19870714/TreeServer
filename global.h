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

#ifndef GLOBAL_H
#define GLOBAL_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <set>
#include <chrono>
#include <atomic>
#include <mutex>
#include <condition_variable>
#include <sys/types.h>
#include <dirent.h>
#include <sys/stat.h>
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>
#include <cstring>
#include <unordered_set>
#include <unordered_map>
#include <fstream>
#include <iostream>

//#define ASSERT
//#define DEBUG_LOG

using namespace std;

#define MASTER_RANK 0
#define REQ_CHANNEL 201
#define RESP_CHANNEL 202
#define STATUS_CHANNEL 203
#define PLAN_CHANNEL 300


#define WAIT_TIME_WHEN_IDLE 100 //unit: usec, user-configurable, used by recv-er
#define STATUS_SYNC_TIME_GAP 100000 //unit: usec, used by Worker main-thread

//values that Column::data_type can take
const int ELEM_BOOL = 1; //binary
const int ELEM_SHORT = 2;
const int ELEM_INT = 3;
const int ELEM_FLOAT = 4;
const int ELEM_DOUBLE = 5;
const int ELEM_CHAR = 6;
const int ELEM_STRING = 7;

// constants related to impurity functions
const int IMPURITY_ENTROPY = 8;
const int IMPURITY_GINI = 9;
const int IMPURITY_CLASSIFICATION_ERROR = 10;
const int IMPURITY_VARIANCE = 11;

const double FEATURE_THRESHOLD = 1e-7;

bool global_end_label = false;

atomic<size_t> task_id_counter(0);

size_t ACTIVE_TREES_THRESHOLD = 50;

// plan from master to slaves
const char SUB_TREE_PLAN = 'a';
const char COL_SPLIT_PLAN = 'b';
const char TASK_DELETE_PLAN = 'c';
const char END_PLAN = 'd';
const char LABEL_FETCH_PLAN = 'e';

// response from slaves to master
const char SUB_TREE_RESP = 'f';
const char COL_SPLIT_RESP = 'g';

// master/slave task_type
const char TASK_SUB_TREE = 'h';
const char TASK_COL_SPLIT = 'i';
const char TASK_DATA_SERVE = 'j';
const char TASK_LEAF = 'k'; // only in master

// for data fetching
const char REQ_COL_FETCH = 'l'; // for sub_tree
const char REQ_ROW_FETCH = 'm'; // fetch row_indices from parent slave task

// for data fetching
const char RESP_COL_FETCH = 'n';
const char RESP_ROW_FETCH = 'o';
const char RESP_LABEL_FETCH = 'p';

//global job parameters:
int BRUTE_FORCE_MAX_ITEM = 10; //allowed size of categorical domain, beyond which the column will be ignored
bool NO_WARNING = true; //warning tag
bool SAVE_TREE = true;
bool PRINT_TREE = false;
int BFS_PRIORITY_THRESHOLD = 100;
int MAX_REPORT_DEPTH = 5; //do not report beyond this level; the report is only useful when building one tree using one thread
int num_compers = 1;
size_t subtree_D = 1000;
int y_index = -1;
string train_file;
string meta_file;
string job_dir = "job"; //to put tree outputs (.io for later loading, .json for visualization)


//========== for serialization (imported from G-thinker) =============

void _mkdir(const char *dir) {//taken from: http://nion.modprobe.de/blog/archives/357-Recursive-directory-creation.html
    char tmp[256];
    char *p = NULL;
    size_t len;
    
    snprintf(tmp, sizeof(tmp), "%s", dir);
    len = strlen(tmp);
    if(tmp[len - 1] == '/') tmp[len - 1] = '\0';
    for(p = tmp + 1; *p; p++)
        if(*p == '/') {
            *p = 0;
            mkdir(tmp, S_IRWXU);
            *p = '/';
        }
    mkdir(tmp, S_IRWXU);
}

void replaceAll(string & str, const string & from, const string & to) {
    if(from.empty()) return;
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length();
    }
}

string filterQuote(const string & str)
{
    string str1 = str;
    replaceAll(str1, "\"", "\\\"");
    return str1;
}

//========== for distributed computing =============

///worker info

int _my_rank;
int _num_workers;
int _num_columns;
size_t _n_samples;

inline int get_worker_id()
{
    return _my_rank;
}
inline int get_num_workers()
{
    return _num_workers;
}


vector<int> all_columns; // all column indices except Y, initialize in Worker::run(), should be synchrinized afterwards

void random_shuffle(int k, vector<int> & result) {

#ifdef ASSERT
	assert(all_columns.size() == (_num_columns - 1));
#endif

    vector<int> all_cols = all_columns;

    std::random_shuffle(all_cols.begin(), all_cols.end());

    for(int i=0; i<k; i++) result.push_back(all_cols[i]);
}

//@@@@@@@@@@@@@@@@@@@@ Da Yan debug todo:: remove
#ifdef DEBUG_LOG
ofstream fout; // todo debug remove later, note that many threads write to the same file; do this for simplicity but may have mixed text streams
#endif

#endif
