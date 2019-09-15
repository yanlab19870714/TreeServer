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

//====== this is a distributed example
//remember to
// 1. using our column-file put to upload [train_file_on_hdfs]
// 2. upload [job-config_on_hdfs] [meta_file_on_hdfs]

#include "csv.h"
#include "config.h"

using namespace std;

//* //testing Worker::run(.)
#include "Worker.h"
#include <iostream>
#include <cassert>

int main(int argc, char** argv){
    WorkerParams params;

    params.job_file_path = "job.config"; //contains the HDFS paths for training data and meta data
    params.tree_file_path = "tree2.config";
    params.test_file_path = "data/bank_test.csv";
    params.test_meta_file = "data/bank_meta.csv";

    if(argc == 3) {
        params.job_file_path = argv[1];
        params.tree_file_path = argv[2];
        params.test_file_path = "";
        params.test_meta_file = "";
    }
    else if(argc == 5) {
        params.job_file_path = argv[1];
        params.tree_file_path = argv[2];
        params.test_file_path = argv[3];
        params.test_meta_file = argv[4];
    }
    else
    {
    	cout<<"Wrong input argument number !!!"<<endl;
    	assert(false);
    }

    Worker worker;
    params.column_assignment = MODE_REPLICATE;
    params.replicate = 2;
    worker.run(params);
}
//*/
