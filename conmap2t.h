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

#ifndef CONMAP2T_H
#define CONMAP2T_H

//concurrent map for two threads

#define CONMAP2T_BUCKET_NUM 100 //used by TaskMap, should allow a comper and the response-processing thread to minimize collision

//idea: 2-level hashing
//1. id % CONMAP_BUCKET_NUM -> bucket_index
//2. bucket[bucket_index] -> give id, get content

#include "tree_msgs.h"
#include <vector>

using namespace std;

template <typename K, typename V> struct conmap2t_bucket
{
	typedef unordered_map<K, V> KVMap;
	mutex mtx;
	KVMap bucket;

	inline void lock()
	{
		mtx.lock();
	}

	inline void unlock()
	{
		mtx.unlock();
	}

	KVMap & get_map()
	{
		return bucket;
	}

	//returns true if inserted
	//false if an entry with this key alreqdy exists
	bool insert(K key, const V & val)
	{
		auto ret = bucket.insert(
			std::pair<K, V>(key, val)
		);
		return ret.second;
	}

	//returns whether deletion is successful
	bool erase(K key)
	{
		size_t num_erased = bucket.erase(key);
		return (num_erased == 1);
	}
};

template <typename K, typename V> struct conmap2t
{
private:
	bool empty() { // only for internal use
		bool isEmpty = true;

		for(size_t i = 0; i < CONMAP2T_BUCKET_NUM; i++) {
			bucket & current = pos(i);
			unordered_map<K, V> & kvmap = current.get_map();
			isEmpty = isEmpty & kvmap.empty();

			if(!kvmap.empty()) { // todo remove later
				for(auto it = kvmap.begin(); it != kvmap.end(); it++) {
                    if(it->second->task_type == TASK_COL_SPLIT) {
                        Task_Slave_Col_Split* task = (Task_Slave_Col_Split*) it->second;
                        cout << "id = " << it->first << ", type = " << it->second->task_type
                             << ", n_cols = " << task->n_cols
							 << ", parent_request = " << task->parent_request
                             << ", |size| = " << task->best_left_rows.indexes.size() + task->best_right_rows.indexes.size() << endl;
                    }
				}
			}
		}

		return isEmpty;
	}
public:
	typedef conmap2t_bucket<K, V> bucket;
	bucket* buckets;
	string name; // todo debug remove later

	conmap2t()
	{
		buckets = new bucket[CONMAP2T_BUCKET_NUM];
	}

	bucket & get_bucket(K key)
	{
		return buckets[key % CONMAP2T_BUCKET_NUM];
	}

	bucket & pos(size_t pos)
	{
		return buckets[pos];
	}

	size_t bucket_hash(K key) {
		return (key % CONMAP2T_BUCKET_NUM);
	}

	~conmap2t()
	{
		if(!empty()) { // todo debug remove later
#ifdef ASSERT
            assert(empty());
#endif
		}
		delete[] buckets;
	}
};

#endif
