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

#ifndef UTILS_H
#define UTILS_H

class Tree_Progress {
private:
    map<int, int> progress_table; // root_task_id->num_plans_in_progress
    size_t tree_build_count = 0;
    mutex lock;

public:

    size_t active_trees() {
        unique_lock<mutex> lck(lock);
        return progress_table.size();
    }

    void increment(int root_task_id) {
        unique_lock<mutex> lck(lock);

        auto it = progress_table.find(root_task_id);

        if(it == progress_table.end()) {
            progress_table[root_task_id] = 1;
        } else {
            it->second++;
        }
    }

    void decrement(int root_task_id) {
        unique_lock<mutex> lck(lock);

        auto it = progress_table.find(root_task_id);

#ifdef ASSERT
        assert(it != progress_table.end());
#endif

        it->second--;

        if(it->second == 0) {
        	cout << "tree = " << it->first << " built, num_trees built so far " << tree_build_count << endl;
            progress_table.erase(it);
            tree_build_count++;
        }
    }

    // todo, can be removed
    void print_table() {
        unique_lock<mutex> lck(lock);
        cout << "tree progress table is : " << endl;

        for(auto it = progress_table.begin(); it != progress_table.end(); it++) {
            cout << "key = " << it->first << ", counter = " << it->second << endl;
        }
    }
};

Tree_Progress tree_progress;
#endif
