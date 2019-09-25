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

//====== this is for loading a model for testing
//it loads a local-file model, a test dataset and its metafile

#include "tree.h"
#include "csv.h"
#include "config.h"
#include "tree_obj.h"

#include <iostream>
using namespace std;

//arg1: local file path for a forest model
//arg2: local test dataset path
//arg3: local metafile path
//arg4: local job-config path

int main(int argc, char** argv){
	PRINT_TREE = true;

    if(argc != 6)
    {
    	cout<<"arg1: local file path for a forest model"<<endl;
    	cout<<"arg2: local test dataset path"<<endl;
    	cout<<"arg3: local metafile path"<<endl;
    	cout<<"arg4: local job-config path (for y_index)"<<endl;
    	cout<<"arg5: is_classification? (0/1)"<<endl;
    	return;
    }

    Matrix test_set;
    load_csv(argv[2], argv[3], test_set);
    load_local_job_config(argv[4]);
    vector<TreeNode*> rootList = load_forest(test_set, y_index, argv[1]);


    if(*argv[5] == '0') {
		double RMSE = calculate_RMSE(rootList, test_set, y_index, INT_MAX);
		cout << "RMSE is : "
			 << RMSE << endl << endl;
	}
    else{
		cout << "Accuracy is : "
			 << get_accuracy(test_set, rootList, y_index, INT_MAX) << endl << endl;
	}

    for(size_t root_idx = 0; root_idx < rootList.size(); root_idx++) {

		if(PRINT_TREE) {
			cout << "##################################################" << endl;
			print_tree(rootList[root_idx], test_set);
			cout << endl << endl;
		}
		free_tree(rootList[root_idx], y_index, test_set);
	}

}
//*/

// ./run job/model_1 data_regression/test.csv data_regression/meta.csv data_regression/job.config 0
