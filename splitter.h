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

#ifndef SPLITTER_H
#define SPLITTER_H

#include "criterion.h"
#include "config.h"
#include "column_server.h"
#include <typeinfo>
#include <mutex>
#include <thread>
#include <random>

using namespace std;

// 1. design sketch
// 2. impurity(sketch) = score
// 3. split based on best score

// #################### column-split #######################
// -------- for regression --------

double impurity_improvement_fn(double left, double right, size_t n_l, size_t n_r);

template <class Tx, class Ty>
struct id_value {
    size_t id;
    Tx x_value;
    Ty y_value;

    friend bool operator< (const id_value &c1, const id_value &c2)
    {
        return c1.x_value < c2.x_value;
    }
};

struct Best {
    double impurity_improvement = -DBL_MAX;
    size_t offset;
    vector<size_t> S1;
    int column_idx = -2; // default / can not split condition
    vector<size_t>::iterator end;
};

struct category_item { //each X-label
    string category; //string used for any Xi-category type
    double sum_1 = 0.0; // this is like sum_1(var_sketch)
    double sum = 0.0; //double for any y-numerical type, like sum_y (var_sketch)
    double average; //!!! needs to be manually set after sum is finalized
    double sum_y2 = 0.0; // calculating (sum * sum)  for variance calculation later
    vector<size_t> indices; //  row indices of "category"
};

struct Category_comparator //for sorting X-labels by y-average
{
    bool operator()(const category_item item1, const category_item item2)
    {
        return item1.average < item2.average;
    }
};


struct var_sketch //sketch for variance
{//motivated by Google's PLANET paper
    double sum_y = 0.0;
    double sum_y2 = 0.0;
    double sum_1 = 0.0;
};

double var(var_sketch & sketch)
{//compute variance from sketch
    double MSS = (double) sketch.sum_y2 / sketch.sum_1; //mean sum of square: E[x^2]
    double mean = (double) sketch.sum_y / sketch.sum_1; //E[x]

    double mean2 = mean * mean; //square of mean

    return MSS - mean2;
}

template<class Tx, class Ty>
void set_sketch(var_sketch & sketch, vector<id_value<Tx, Ty>> &pairs)
{//compute sketch over rows [start, end), using labels y[row_idx]
#ifdef ASSERT
    assert(typeid(Ty) != typeid(char));
    assert(typeid(Ty) != typeid(string));
    assert(typeid(Ty) != typeid(bool));
#endif

    sketch.sum_y = 0.0;
    sketch.sum_y2 = 0.0;
    sketch.sum_1 = 0.0;

    double d_value;

    for(size_t i = 0; i < pairs.size(); i++) {
        d_value = pairs[i].y_value;

        sketch.sum_y += d_value;
        sketch.sum_y2 += (d_value * d_value);
        sketch.sum_1 += 1;
    }

}

template <class T>
void update_sketches(var_sketch & sketch_L, var_sketch & sketch_R, T value) //T should be numeric
{//incrementally update left and right child-nodes' sketches,
    //after moving "value" from right child to left child
#ifdef ASSERT
    assert(typeid(T) != typeid(char));
    assert(typeid(T) != typeid(string));
    assert(typeid(T) != typeid(bool));
#endif

    double d_value = value;
    sketch_L.sum_y += d_value;
    sketch_L.sum_y2 += (d_value * d_value);
    sketch_L.sum_1 += 1;

    sketch_R.sum_y -= d_value;
    sketch_R.sum_y2 -= (d_value * d_value);
    sketch_R.sum_1 -= 1;
}

// for batch sketch update for categorical Xi
void update_sketches(var_sketch & sketch_L, var_sketch & sketch_R, category_item &value)
{//incrementally update left and right child-nodes' batch sketches for  categorical data,
    //after moving "value" from right child to left child
    sketch_L.sum_y += value.sum;
    sketch_L.sum_y2 += value.sum_y2;
    sketch_L.sum_1 += value.sum_1;

    sketch_R.sum_y -= value.sum;
    sketch_R.sum_y2 -= value.sum_y2;
    sketch_R.sum_1 -= value.sum_1;
}

// #################### column-split for classification #######################

template <class T>
class classify_sketch {
public:
    map<T, size_t> freq; //freq[class] = # of times the class occurred
    size_t N = 0; //total number of elements
    //prob_class = freq[class] / N, used for computing impurity

    void add(T key) {
        auto it = freq.find(key);

        if(it == freq.end()) freq[key] = 1;
        else it->second++;

        N++;
    }

    void remove(T key) {
        auto it = freq.find(key);

#ifdef ASSERT
        assert(it->second > 0);
        assert(it != freq.end());
#endif
        it->second--;

        if(it->second == 0) freq.erase(it);
        N--;
    }

    void add(map<string, size_t> &f_y) {
        for(auto it = f_y.begin(); it != f_y.end(); it++) {
            T key = boost::lexical_cast<T>(it->first);
            auto it2 = freq.find(key);

            if(it2 == freq.end()) {
                freq[key] = it->second;
            } else {
                it2->second += it->second;
            }

            N += it->second;
        }
    }

    void remove(map<string, size_t> &f_y) {
        for(auto it = f_y.begin(); it != f_y.end(); it++) {
            T key = boost::lexical_cast<T>(it->first);
            auto it2 = freq.find(key);

#ifdef ASSERT
            assert(it2 != freq.end());
#endif
            it2->second -= it->second;

#ifdef ASSERT
            assert(it2->second >= 0);
#endif
            if(it2->second == 0) {
                freq.erase(it2);
            }

            N -= it->second;
        }
    }

    void set_ratio(vector<double> &target_ratios) {
        for(auto it = freq.begin(); it != freq.end(); it++) {
            double ratio = (double) it->second / N;
            target_ratios.push_back(ratio);
        }
    }

};

template <class Tx, class Ty>
void set_sketch(classify_sketch<Ty> &sketch, vector<id_value<Tx, Ty>> &pairs)
{//compute sketch over rows [start, end), using labels y[row_idx]
    for(size_t i = 0; i < pairs.size(); i++) {
        id_value<Tx, Ty> &pair =  pairs[i];
        sketch.add(pair.y_value);
    }
}

template <class T>
void update_sketches(classify_sketch<T> &sketch_L, classify_sketch<T> &sketch_R, T value) //T should be numeric
{//incrementally update left and right child-nodes' sketches,
    //after moving "value" from right child to left child
    sketch_L.add(value);
    sketch_R.remove(value);
}

// #################### brute-force #######################

// category item for brute force algorithm
struct Item {
    string category;
    vector<size_t> indices; // row indices
    map<string, size_t> freq_y; // y frequency
};


// updating the iterator [start, end) in order vector<Item> &items
// helper function to update original iterator indices [start, end)
void update_iterator(vector<size_t>::iterator start, vector<Item> &items) {
    // updating indices from [start, end) based on sorted order
    int offset = 0;

    for(size_t i = 0; i < items.size(); i++) {
        vector<size_t> &indices = items[i].indices;
        std::copy(indices.begin(), indices.end(), start + offset);
        offset += indices.size(); //  (end  - begin)
    }
}

// scan [start, end) to group items by Xi-label
// return a item vector "ordered_items" for enumerating S1
// "ordered_items" can be a vector of labels of any order: (Apple, Banana, ..) or (Banana, Apple,...)
template<class Tx, class Ty>
bool set_grouped_item(vector<id_value<Tx, Ty>> &pairs, vector<Item> &ordered_items,
                      Best &container) {
    // return is_skip? true to skip splitting Xi
    map<string, Item> category_map; // collecting unique categories

    for(size_t i = 0; i < pairs.size(); i++) {
        id_value<Tx, Ty> &pair  = pairs[i];

        // TODO::
        string x_value = boost::lexical_cast<string>(pair.x_value);
        // TODO::
        string y_value = boost::lexical_cast<string>(pair.y_value);

        auto it2 = category_map.find(x_value);

        if(it2 == category_map.end()) {
            if(category_map.size() > BRUTE_FORCE_MAX_ITEM) {
                if(!NO_WARNING) cout<< "[Warning] (Xi Categorical, Y Categorical), |category_item| = "
                    << category_map.size() << " exceeds the limit " << BRUTE_FORCE_MAX_ITEM << endl;
                return true;
            }

            Item item;
            item.category = x_value;
            item.indices.push_back(pair.id);
            item.freq_y[y_value] = 1;
            category_map[x_value] = item;
        } else {
            map<string, size_t> &freq_y = it2->second.freq_y;

            // creating y_freq map for each category
            auto it3 = freq_y.find(y_value);
            if(it3 == freq_y.end()) {
                freq_y[y_value] = 1;
            } else {
                it3->second++;
            }

            it2->second.indices.push_back(pair.id);
        }
    }

    // collecting category items in a vector
    for(auto it = category_map.begin(); it != category_map.end(); it++) {
        ordered_items.push_back(it->second);
    }

    return false;
}

// helper function called from get_impurity_improvement(..)
// to calculate impurity based on Y data-type
template <class T>
double impurity_improvement(vector<Item> &ordered_items, vector<size_t> &s1, TreeConfig &treeConfig) {
    //output:
    //impurity_improvement

    classify_sketch<T> left;
    classify_sketch<T> right;

    vector<bool> isS2(ordered_items.size(), true);

    for(size_t i = 0; i < s1.size(); i++) {
        left.add(ordered_items[s1[i]].freq_y);
        isS2[s1[i]] = false;
    }

    for(size_t i = 0; i < isS2.size(); i++) {
        if(isS2[i]) {
            right.add(ordered_items[i].freq_y);
        }
    }

    //compute pi
    vector<double> left_tgt_ratio;
    vector<double> right_tgt_ratio;
    left.set_ratio(left_tgt_ratio);
    right.set_ratio(right_tgt_ratio);

    double left_impurity;
    double right_impurity;
    if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY) {
        left_impurity = entropy(left_tgt_ratio);
        right_impurity = entropy(right_tgt_ratio);
    } else if (treeConfig.IMPURITY_FUNC == IMPURITY_GINI) {
        left_impurity = gini(left_tgt_ratio);
        right_impurity = gini(right_tgt_ratio);
    } else if (treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR) {
        left_impurity = classification_error(left_tgt_ratio);
        right_impurity = classification_error(right_tgt_ratio);
    } else {
        cout<< "[ERROR] for classification we only consider entropy, "
                "gini and classification error as impurity function " << endl;
        exit(-1);
    }

    return impurity_improvement_fn(left_impurity, right_impurity, left.N, right.N);
}

// calculate impurity improvement for a particular combination |S1| and |S-S1|
double get_impurity_improvement(Matrix & dataset, vector<Item> &ordered_items,
                                vector<size_t> &s1,
                                TreeConfig &treeConfig) {
    int type = dataset.get_column(y_index)->data_type;

    if(type == ELEM_BOOL) {
        return impurity_improvement<bool>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_SHORT) {
        return impurity_improvement<short>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_INT) {
        return impurity_improvement<int>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_FLOAT) {
        return impurity_improvement<float>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_DOUBLE) {
        return impurity_improvement<double>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_CHAR) {
        return impurity_improvement<char>(ordered_items, s1, treeConfig);
    } else if(type == ELEM_STRING) {
        return impurity_improvement<string>(ordered_items, s1, treeConfig);
    } else
    {
        cout<<"Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
        exit(-1);
    }
}

// helper function to calculate the offset from ordered_items and S1
size_t get_offset(vector<Item> &ordered_items, vector<size_t> &S1) {
    size_t offset = 0;

    for(size_t i = 0; i < S1.size(); i++) {
        offset += ordered_items[S1[i]].indices.size();
    }

    return offset;
}

// helper function of find_best_split_enum()
void do_enumerate(Matrix & dataset, size_t curr_val, vector<size_t> &temp_s1,
                  vector<Item> &ordered_items, Best &current,
                  TreeConfig &treeConfig) {
    //N: # of different Xi values
    //curr_val: latest (largest) element in enumerating set S1
    //temp_s1: current elements in S1

    // k = 1 .. floor(n/2), n is odd (no duplicate), n even (duplicate)
    size_t k_range = ordered_items.size() / 2;

    if(temp_s1.size() > 0)
        if((ordered_items.size() % 2 == 0) && (temp_s1.size() == k_range - 1) && (temp_s1[0] > 0))
            return;

    temp_s1.push_back(curr_val);

    double impurity = get_impurity_improvement(dataset, ordered_items, temp_s1, treeConfig);
    if(impurity > current.impurity_improvement) {
        current.impurity_improvement = impurity;
        current.S1 = temp_s1;
        current.offset = get_offset(ordered_items, temp_s1);
    }

    if(temp_s1.size() < k_range) {
        for(size_t i = curr_val + 1; i < ordered_items.size(); i++) {
            do_enumerate(dataset, i, temp_s1, ordered_items, current, treeConfig);
        }
    }

    temp_s1.pop_back();
}

// param: number of unique item in vector<Item>
// output: best_impurity_improvement, best_offset, best_S1
void find_best_split_enum(Matrix & dataset, vector<Item> &ordered_items, Best &current, TreeConfig &treeConfig) {
    size_t N = ordered_items.size();

    if(N <= 1) {
        cout << "[ERROR] Line : " << __LINE__ << ": N should be greater than 1, current N = " << N << endl;
        exit(-1);
    }

    current.impurity_improvement = -DBL_MAX;

    vector<size_t> temp_s1;// tree path stack

    for(size_t i = 0; i < N; i++) {
        do_enumerate(dataset, i, temp_s1, ordered_items, current, treeConfig); //i-root prefix tree enumeration: i must be in S1
    }
}

// updating the vector<Item> &ordered_item with the new S1
// the best S1 found while doing brute force impurity calculation
void update_category_order(vector<size_t>::iterator start, vector<Item> &ordered_item,
                           vector<size_t> &S1) {

    vector<Item> new_order;
    vector<bool> isS2(ordered_item.size(), true);

    for(size_t i = 0; i < S1.size(); i++) {
        Item &item = ordered_item[S1[i]];
        new_order.push_back(item);
        isS2[S1[i]] = false;
    }

    for(size_t i = 0; i < isS2.size(); i++) {
        if(isS2[i]) {
            Item &item = ordered_item[i];
            new_order.push_back(item);
        }
    }

    update_iterator(start, new_order);
}

// #################### node-split #######################

class SplitResult {
public:
    // column (attribute) index to split on
    int column_idx;
    // position (+ 1) in the sample index array to split for column_idx
    vector<size_t>::iterator pos;
    double best_impurity;
};

template <class T>
class OrdinalSplitResult : public SplitResult {
public:
    // attribute value on which to split column left and right
    T attribute_value;
    OrdinalSplitResult(){}
};

template <class T>
class CategoricalSplitResult : public SplitResult {
public:
    set<T> S; // unique categories between [start, end)
    set<T> S1;
};

void free(SplitResult* sr) {
    Column* column = cserver.X.get_column(sr->column_idx);
    int type = column->data_type;
    bool is_ordinal = column->is_ordinal;

    if(is_ordinal) {
        if(type == ELEM_BOOL) {
            OrdinalSplitResult<bool>* result = (OrdinalSplitResult<bool>*) sr;
            delete result;
        } else if(type == ELEM_SHORT) {
            OrdinalSplitResult<short>* result = (OrdinalSplitResult<short>*) sr;
            delete result;
        } else if(type == ELEM_INT) {
            OrdinalSplitResult<int>* result = (OrdinalSplitResult<int>*) sr;
            delete result;
        } else if(type == ELEM_FLOAT) {
            OrdinalSplitResult<float>* result = (OrdinalSplitResult<float>*) sr;
            delete result;
        } else if(type == ELEM_DOUBLE) {
            OrdinalSplitResult<double>* result = (OrdinalSplitResult<double>*) sr;
            delete result;
        } else if(type == ELEM_CHAR) {
            OrdinalSplitResult<char>* result = (OrdinalSplitResult<char>*) sr;
            delete result;
        } else if(type == ELEM_STRING) {
            OrdinalSplitResult<string>* result = (OrdinalSplitResult<string>*) sr;
            delete result;
        }
    }
    else {
    	if(type == ELEM_BOOL) {
    		CategoricalSplitResult<bool>* result = (CategoricalSplitResult<bool>*) sr;
			delete result;
		} else if(type == ELEM_SHORT) {
			CategoricalSplitResult<short>* result = (CategoricalSplitResult<short>*) sr;
			delete result;
		} else if(type == ELEM_INT) {
			CategoricalSplitResult<int>* result = (CategoricalSplitResult<int>*) sr;
			delete result;
		} else if(type == ELEM_FLOAT) {
			CategoricalSplitResult<float>* result = (CategoricalSplitResult<float>*) sr;
			delete result;
		} else if(type == ELEM_DOUBLE) {
			CategoricalSplitResult<double>* result = (CategoricalSplitResult<double>*) sr;
			delete result;
		} else if(type == ELEM_CHAR) {
			CategoricalSplitResult<char>* result = (CategoricalSplitResult<char>*) sr;
			delete result;
		} else if(type == ELEM_STRING) {
			CategoricalSplitResult<string>* result = (CategoricalSplitResult<string>*) sr;
			delete result;
		}
    }
}

double impurity_improvement_fn(double left, double right, size_t n_l, size_t n_r) {
    return (- left * n_l - right * n_r);
}


// find the best position for a column
// !! sort [start, end) first !!! X,Y both continuous (ordinal)
template <class Tx, class Ty>
void set_best_pos_classification(vector<id_value<Tx, Ty>> &pairs, Best &current, TreeConfig &treeConfig) {
    //output:
    //Best current

#ifdef ASSERT
    assert(typeid(Tx) != typeid(char));
    assert(typeid(Tx) != typeid(string));
#endif
    //ToDo :: ordered string is rare, can implement if needed

    current.impurity_improvement = -DBL_MAX;

    classify_sketch<Ty> left;
    classify_sketch<Ty> right;

    set_sketch<Tx, Ty>(right, pairs);

    for (size_t i = 0; i < pairs.size() - 1; i++) {
        id_value<Tx, Ty> &pair = pairs[i];
        update_sketches<Ty>(left, right, pair.y_value); //update with the current element

        if(pairs[i+1].x_value == pair.x_value) {
            continue;
        }

        //compute pi
        vector<double> left_tgt_ratio;
        vector<double> right_tgt_ratio;
        left.set_ratio(left_tgt_ratio);
        right.set_ratio(right_tgt_ratio);

        double left_impurity;
        double right_impurity;
        if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY) {
            left_impurity = entropy(left_tgt_ratio);
            right_impurity = entropy(right_tgt_ratio);
        } else if (treeConfig.IMPURITY_FUNC == IMPURITY_GINI) {
            left_impurity = gini(left_tgt_ratio);
            right_impurity = gini(right_tgt_ratio);
        } else if (treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR) {
            left_impurity = classification_error(left_tgt_ratio);
            right_impurity = classification_error(right_tgt_ratio);
        } else {
            cout<< "[ERROR] for classification we only consider entropy, "
                    "gini and classification error as impurity function " << endl;
            exit(-1);
        }

        double impurity_improvement = impurity_improvement_fn(left_impurity, right_impurity, left.N, right.N);

        if (current.impurity_improvement < impurity_improvement) {
            current.impurity_improvement = impurity_improvement;
            current.offset = i; // (it - start) ==  (i - 0)
        }
    }

}

// find the best position for a column
// !!! sort [start, end) first !!! X,Y both continuous (ordinal)
template <class Tx, class Ty>
void set_best_pos_regression(vector<id_value<Tx, Ty>> &pairs, Best &current) {
    //output:
    //best_impurity_improvement
    //best offest: the chosen split position from "start"

    // feature_value (Xi) of type numerical expected (ref: FEATURE_THRESHOLD)
#ifdef ASSERT
    assert(typeid(Tx) != typeid(char));
    assert(typeid(Tx) != typeid(string));
    assert(typeid(Tx) != typeid(bool));
#endif

    current.impurity_improvement = -DBL_MAX;

    var_sketch left;
    var_sketch right;

    set_sketch<Tx, Ty>(right, pairs);

    for (size_t i = 0; i < pairs.size() - 1; i++) {
        id_value<Tx, Ty> &pair = pairs[i];
        update_sketches<Ty>(left, right, pair.y_value); //update with the current element

        if(pairs[i+1].x_value == pair.x_value) {
            continue;
        }

        double var_L = var(left);
        double var_R = var(right);

        double impurity_improvement = impurity_improvement_fn(var_L, var_R, left.sum_1, right.sum_1);

        if (current.impurity_improvement < impurity_improvement) {
            current.impurity_improvement = impurity_improvement;
            current.offset = i; //like (it - start) == (i - 0);
        }
    }

}


// find the best position for a column (categorical value)
// !!! sort [start, end) first !!! Xi categorical, Y continuous (ordinal)
// need to update sketch on batch
template <class Tx, class Ty>
void set_best_pos_regression_breiman(vector<id_value<Tx, Ty>> &pairs, Best &current,
                                     vector<category_item> &category_items) {
    //output:
    //best_impurity_improvement
    //best offest: the chosen split position from "start"

    current.impurity_improvement = -DBL_MAX;

    var_sketch left;
    var_sketch right;
    set_sketch<Tx, Ty>(right, pairs);

    int total_indices = 0;

    for(size_t i = 0; i < category_items.size(); i++) {
        current.S1.push_back(i);
        update_sketches(left, right, category_items[i]);
        total_indices += category_items[i].indices.size();

        double var_L = var(left);
        double var_R = var(right);

        double impurity_improvement = impurity_improvement_fn(var_L, var_R, left.sum_1, right.sum_1);

        if (current.impurity_improvement < impurity_improvement) {
            current.impurity_improvement = impurity_improvement;
            current.offset = total_indices;// (end - start)
        }
    }

}

template <class Tx, class Ty>
void set_sorted_category_item(vector<id_value<Tx, Ty>> &pairs,
                              vector<category_item> &category_items) { // breiman, y should be numerical

    map<string, category_item> item_map;

    for(size_t i = 0; i < pairs.size(); i++) { //scan data to get X-label groups and their y-avg's accumulating values
        id_value<Tx, Ty> &pair = pairs[i];

        string attr_val = boost::lexical_cast<string>(pair.x_value); // TODO::

        double y_value = boost::lexical_cast<double>(pair.y_value); // TODO::

        auto it2 = item_map.find(attr_val);

        if(it2 == item_map.end()) {
            category_item ct;
            ct.category = attr_val;
            ct.sum_1 += 1;
            ct.sum += y_value;
            ct.sum_y2 += y_value * y_value;
            ct.indices.push_back(pair.id);
            item_map[attr_val] = ct;

        } else {
            it2->second.sum_1++;
            it2->second.sum += y_value;
            it2->second.sum_y2 += y_value * y_value;
            it2->second.indices.push_back(pair.id);
        }

    }
    //compute y-avg for each x-group
    for(auto it = item_map.begin(); it != item_map.end(); it++){
        it->second.average = it->second.sum / it->second.sum_1;
        //cout<< "key = " << it->second.category << ", avg = " << it->second.average << endl;
        category_items.push_back(it->second);
    }
    //sort the X-groups by y-avg
    Category_comparator ctg_comparator;
    std::sort(category_items.begin(), category_items.end(), ctg_comparator);
    //=================================
}

// setting the sorted indices (sorted based on category sort) from  [start, end)
void set_sorted_indices(vector<size_t>::iterator start, vector<category_item> &category_items) {
    int offset = 0;

    for(size_t i = 0; i < category_items.size(); i++) {
        vector<size_t> &indices = category_items[i].indices;
        std::copy(indices.begin(), indices.end(), start + offset);
        offset += indices.size(); //  (end  - begin)
    }

}

// wrapper function to do the sorting
// we need vector<category_item> &category_items in other parts of breiman's algorithm
template <class Tx, class Ty>
void row_sort_categorical(vector<size_t>::iterator start, vector<id_value<Tx, Ty>> &pairs,
                          vector<category_item> &category_items){
    set_sorted_category_item<Tx, Ty>(pairs, category_items);
    //set_sorted_indices(start, category_items); should be done only in best case, after scanning all the feature column
}

//form the final result after all columns are checked
//!!! return value needs to be deleted outside !!!
template <class Tx, class Ty>
SplitResult* create_split_result(Matrix & dataset, vector<size_t>::iterator start, Best &best,
                                 vector<id_value<Tx, Ty>> &pairs,
                                 vector<category_item> &category_items,
                                 vector<Item> &ordered_items) {
    //shuffle back row_idx array for recursive split later

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* Y = dataset.get_column(y_index);

    if(best_Xi->is_ordinal) {
#ifdef ASSERT
        assert(pairs.size() > 0);
#endif
        OrdinalSplitResult<Tx>* ordinalSplitResult = new OrdinalSplitResult<Tx>();
        ordinalSplitResult->pos = start + best.offset + 1;
        ordinalSplitResult->column_idx = best.column_idx;
        ordinalSplitResult->best_impurity = best.impurity_improvement;

        // split_value (0 + offset + 1)
        ordinalSplitResult->attribute_value = boost::lexical_cast<Tx>(pairs[best.offset + 1].x_value);

        return ordinalSplitResult;
    } else if (Y->is_ordinal) { // breiman's algorithm Xi categorical, Y ordinal (continuous)
#ifdef ASSERT
        assert(category_items.size() > 0);
#endif

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < category_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    } else {//if (!Y->is_ordinal) // best_Xi categorical and Y categorical brute force algorithm

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < ordered_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    }
}

//form the final result after all columns are checked
//!!! return value needs to be deleted outside !!!
template <class Tx, class Ty>
SplitResult* create_split_result_mid(Matrix & dataset, vector<size_t>::iterator start, Best &best, //take mid value for floating point numbers
                                 vector<id_value<Tx, Ty>> &pairs,
                                 vector<category_item> &category_items,
                                 vector<Item> &ordered_items) {
    //shuffle back row_idx array for recursive split later

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* Y = dataset.get_column(y_index);

    if(best_Xi->is_ordinal) {
#ifdef ASSERT
        assert(pairs.size() > 0);
#endif
        OrdinalSplitResult<Tx>* ordinalSplitResult = new OrdinalSplitResult<Tx>();
        ordinalSplitResult->pos = start + best.offset + 1;
        ordinalSplitResult->column_idx = best.column_idx;
        ordinalSplitResult->best_impurity = best.impurity_improvement;

        // split_value (0 + offset + 1)
        Tx mid = pairs[best.offset].x_value + pairs[best.offset + 1].x_value;
        ordinalSplitResult->attribute_value = mid/2;

        return ordinalSplitResult;
    } else if (Y->is_ordinal) { // breiman's algorithm Xi categorical, Y ordinal (continuous)
#ifdef ASSERT
        assert(category_items.size() > 0);
#endif

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < category_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    } else {//if (!Y->is_ordinal) // best_Xi categorical and Y categorical brute force algorithm

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < ordered_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    }
}

//form the final result after all columns are checked
//!!! return value needs to be deleted outside !!!
template <class Tx, class Ty>
SplitResult* create_split_result_mid_int(Matrix & dataset, vector<size_t>::iterator start, Best &best, //take mid value for floating point numbers
                                 vector<id_value<Tx, Ty>> &pairs,
                                 vector<category_item> &category_items,
                                 vector<Item> &ordered_items) {
    //shuffle back row_idx array for recursive split later

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* Y = dataset.get_column(y_index);

    if(best_Xi->is_ordinal) {
#ifdef ASSERT
        assert(pairs.size() > 0);
#endif
        OrdinalSplitResult<Tx>* ordinalSplitResult = new OrdinalSplitResult<Tx>();
        ordinalSplitResult->pos = start + best.offset + 1;
        ordinalSplitResult->column_idx = best.column_idx;
        ordinalSplitResult->best_impurity = best.impurity_improvement;

        // split_value (0 + offset + 1)
        Tx mid = pairs[best.offset].x_value + pairs[best.offset + 1].x_value;
        if(mid % 2 == 1) mid++; //round up
        ordinalSplitResult->attribute_value = mid/2;

        return ordinalSplitResult;
    } else if (Y->is_ordinal) { // breiman's algorithm Xi categorical, Y ordinal (continuous)
#ifdef ASSERT
        assert(category_items.size() > 0);
#endif

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < category_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(category_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    } else {//if (!Y->is_ordinal) // best_Xi categorical and Y categorical brute force algorithm

        CategoricalSplitResult<Tx>* splitResult = new CategoricalSplitResult<Tx>();
        splitResult->pos = start + best.offset;
        splitResult->column_idx = best.column_idx;
        splitResult->best_impurity = best.impurity_improvement;

        for(size_t i = 0; i < best.S1.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[best.S1[i]].category);
            splitResult->S1.insert(value);
        }

        for(size_t i = 0; i < ordered_items.size(); i++) {
            Tx value = boost::lexical_cast<Tx>(ordered_items[i].category);
            splitResult->S.insert(value);
        }

        return splitResult;
    }
}

template <class Tx, class Ty>
void align_missing(Matrix & dataset, int column_index, vector<size_t>::iterator start,
                   vector<size_t>::iterator &end, vector<id_value<Tx, Ty>> &pairs) {

    Column *Xi = dataset.get_column(column_index);
    Column *y =  dataset.get_column(y_index);

    if(!Xi->is_dense || !y->is_dense) {
        cout<< "File = " << __FILE__ << ", Line = " << __LINE__
            << ": We have not considered sparse column yet" << endl;
        exit(-1);
    }

    DenseColumn<Tx>* denseColumn = (DenseColumn<Tx>*) Xi;

    while (start != end) {
        Tx x_val;
        Xi->get(*start, &x_val);

        if(x_val == denseColumn->missing_value) {
            std::swap(*start, *(end - 1));
            end--;
        } else {
            Ty y_val;
            y->get(*start, &y_val);

            id_value<Tx, Ty> pair;
            pair.id = *start;
            pair.x_value = x_val;
            pair.y_value = y_val;
            pairs.push_back(pair);

            start++;
        }
    }
}

template <class Tx, class Ty>
void update_pair_indices(vector<size_t>::iterator start, vector<id_value<Tx, Ty>> &pairs) {
    for(size_t i = 0; i < pairs.size(); i++) {
        *start = pairs[i].id; // updating memory
        start++;
    }
}

template <class Tx, class Ty>
void update_best_split_ordinal_classification(Matrix &dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                               Best &current, Best &best) {
    vector<id_value<Tx, Ty>> pairs;

    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);

    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = " << current.column_idx << " (all are missing)" << endl;
        return;
    }

    std::sort(pairs.begin(), pairs.end());
    set_best_pos_classification<Tx, Ty>(pairs, current, treeConfig);

    //---
    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.end = current.end;
    }
}

template <class Tx, class Ty>
void update_best_split_ordinal_regression(Matrix &dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                                              Best &current, Best &best) {


    vector<id_value<Tx, Ty>> pairs;

    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);

    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = " << current.column_idx << " (all are missing)" << endl;
        return;
    }

    std::sort(pairs.begin(), pairs.end());
    set_best_pos_regression<Tx, Ty>(pairs, current);

    //---
    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.end = current.end;
    }
}

// X categorical, Y continuous
template <class Tx, class Ty>
void update_best_split_breiman(Matrix &dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                               Best &current, Best &best) {

    vector<id_value<Tx, Ty>> pairs;
    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);

    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = " << current.column_idx << " (all are missing)" << endl;
        return;
    }

    vector<category_item> category_items; // sorted category items (of Xi) with indices
    row_sort_categorical(start, pairs, category_items);

    if(category_items.size() <= 1) {
        if(!NO_WARNING) cout<< "[Warning] Skip the categorical feature with index = " << current.column_idx << " (only <=1 items, nothing to split)" << endl;
        return;
    }

    set_best_pos_regression_breiman<Tx, Ty>(pairs, current, category_items);

    //---
    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.S1 = current.S1;
        best.end = current.end;
    }
}

template <class Tx, class Ty>
void update_best_split_brute_force(Matrix & dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                                   Best &current, Best &best) {

    vector<id_value<Tx, Ty>> pairs;
    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);

    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = "
                            << current.column_idx<< " (all are missing)" << endl;
        return;
    }

    vector<Item> ordered_items;
    bool skip_Xi = set_grouped_item<Tx, Ty>(pairs, ordered_items, current);

    if(skip_Xi) {
        if(!NO_WARNING) cout<< "[Warning] Skip the categorical feature with index = "
                            << current.column_idx << " (too many items to enumerate)" << endl;
        return;
    }

    if(ordered_items.size() <= 1) {
        if(!NO_WARNING) cout<< "[Warning] Skip the categorical feature with index = "
                            << current.column_idx << " (only <=1 items, nothing to split)" << endl;
        return;
    }

    //---
    find_best_split_enum(dataset, ordered_items, current, treeConfig);

    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.S1 = current.S1;
        best.end = current.end;
    }
};

template<class Tx, class Ty>
SplitResult* create_split_result_wrapper(Matrix &dataset, vector<size_t>::iterator start,
                                         vector<size_t>::iterator &end, Best &best) {

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* y = dataset.get_column(y_index);

    vector<id_value<Tx, Ty>> pairs; // Xi ordinal
    vector<category_item> category_items; // breiman
    vector<Item> ordered_items; // brute force

    align_missing<Tx, Ty>(dataset, best.column_idx, start, end, pairs);

    if(best_Xi->is_ordinal) {
        std::sort(pairs.begin(), pairs.end());
        update_pair_indices<Tx, Ty>(start, pairs);

        return create_split_result<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else if (y->is_ordinal) { // breiman
        row_sort_categorical(start, pairs, category_items);
        set_sorted_indices(start, category_items);

        return create_split_result<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else //if (!y->is_ordinal)
    { // brute force
        set_grouped_item<Tx, Ty>(pairs, ordered_items, best);
        update_category_order(start, ordered_items, best.S1);

        return create_split_result<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    }
}

template<class Tx, class Ty>
SplitResult* create_split_result_mid_wrapper(Matrix &dataset, vector<size_t>::iterator start,
                                         vector<size_t>::iterator &end, Best &best) {

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* y = dataset.get_column(y_index);

    vector<id_value<Tx, Ty>> pairs; // Xi ordinal
    vector<category_item> category_items; // breiman
    vector<Item> ordered_items; // brute force

    align_missing<Tx, Ty>(dataset, best.column_idx, start, end, pairs);

    if(best_Xi->is_ordinal) {
        std::sort(pairs.begin(), pairs.end());
        update_pair_indices<Tx, Ty>(start, pairs);

        return create_split_result_mid<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else if (y->is_ordinal) { // breiman
        row_sort_categorical(start, pairs, category_items);
        set_sorted_indices(start, category_items);

        return create_split_result_mid<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else //if (!y->is_ordinal)
    { // brute force
        set_grouped_item<Tx, Ty>(pairs, ordered_items, best);
        update_category_order(start, ordered_items, best.S1);

        return create_split_result_mid<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    }
}

template<class Tx, class Ty>
SplitResult* create_split_result_mid_int_wrapper(Matrix &dataset, vector<size_t>::iterator start,
                                         vector<size_t>::iterator &end, Best &best) {

    Column* best_Xi = dataset.get_column(best.column_idx);
    Column* y = dataset.get_column(y_index);

    vector<id_value<Tx, Ty>> pairs; // Xi ordinal
    vector<category_item> category_items; // breiman
    vector<Item> ordered_items; // brute force

    align_missing<Tx, Ty>(dataset, best.column_idx, start, end, pairs);

    if(best_Xi->is_ordinal) {
        std::sort(pairs.begin(), pairs.end());
        update_pair_indices<Tx, Ty>(start, pairs);

        return create_split_result_mid_int<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else if (y->is_ordinal) { // breiman
        row_sort_categorical(start, pairs, category_items);
        set_sorted_indices(start, category_items);

        return create_split_result_mid_int<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    } else //if (!y->is_ordinal)
    { // brute force
        set_grouped_item<Tx, Ty>(pairs, ordered_items, best);
        update_category_order(start, ordered_items, best.S1);

        return create_split_result_mid_int<Tx, Ty>(dataset, start, best, pairs, category_items, ordered_items);
    }
}

void find_best_split_column(Matrix & dataset, vector<size_t>::iterator start, Best & current,
                            Best & best, TreeConfig & treeConfig) {
    //sort row-idx array based on Xi
    Column* Xi = dataset.get_column(current.column_idx);
    Column* y = dataset.get_column(y_index);

    if(Xi->is_ordinal)
    {
        if(treeConfig.IMPURITY_FUNC == IMPURITY_VARIANCE) {
            if(Xi->data_type ==  ELEM_SHORT) {
                if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_regression<short, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_regression<short, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_regression<short, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_regression<short, double>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_INT) {
                if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_regression<int, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_regression<int, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_regression<int, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_regression<int, double>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_FLOAT) {
                if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_regression<float, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_regression<float, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_regression<float, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_regression<float, double>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_DOUBLE) {
                if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_regression<double, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_regression<double, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_regression<double, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_regression<double, double>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            }
        } else if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY
                   || treeConfig.IMPURITY_FUNC == IMPURITY_GINI
                   || treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR) {

            if(Xi->data_type ==  ELEM_SHORT) {
                if(y->data_type ==  ELEM_BOOL) {
                    update_best_split_ordinal_classification<short, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_classification<short, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_classification<short, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_classification<short, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_classification<short, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_CHAR) {
                    update_best_split_ordinal_classification<short, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_STRING) {
                    update_best_split_ordinal_classification<short, string>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_INT) {
                if(y->data_type ==  ELEM_BOOL) {
                    update_best_split_ordinal_classification<int, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_classification<int, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_classification<int, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_classification<int, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_classification<int, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_CHAR) {
                    update_best_split_ordinal_classification<int, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_STRING) {
                    update_best_split_ordinal_classification<int, string>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_FLOAT) {
                if(y->data_type ==  ELEM_BOOL) {
                    update_best_split_ordinal_classification<float, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_classification<float, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_classification<float, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_classification<float, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_classification<float, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_CHAR) {
                    update_best_split_ordinal_classification<float, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_STRING) {
                    update_best_split_ordinal_classification<float, string>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            } else if(Xi->data_type ==  ELEM_DOUBLE) {
                if(y->data_type ==  ELEM_BOOL) {
                    update_best_split_ordinal_classification<double, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_SHORT) {
                    update_best_split_ordinal_classification<double, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_INT) {
                    update_best_split_ordinal_classification<double, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_FLOAT) {
                    update_best_split_ordinal_classification<double, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_DOUBLE) {
                    update_best_split_ordinal_classification<double, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_CHAR) {
                    update_best_split_ordinal_classification<double, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type ==  ELEM_STRING) {
                    update_best_split_ordinal_classification<double, string>(dataset, start, treeConfig, current, best);
                } else {
                    cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
                    exit(-1);
                }
            }
        }

    } else if (y->is_ordinal) { // is Xi is categorical, Y ordinal, Breiman's alogrithm

        if(treeConfig.IMPURITY_FUNC == IMPURITY_VARIANCE) {
            if(Xi->data_type ==  ELEM_BOOL) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<bool, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<bool, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<bool, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<bool, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_SHORT) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<short, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<short, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<short, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<short, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_INT) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<int, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<int, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<int, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<int, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_FLOAT) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<float, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<float, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<float, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<float, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_DOUBLE) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<double, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<double, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<double, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<double, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_CHAR) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<char, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<char, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<char, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<char, double>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type ==  ELEM_STRING) {
                if(y->data_type == ELEM_SHORT) {
                    update_best_split_breiman<string, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_breiman<string, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_breiman<string, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_breiman<string, double>(dataset, start, treeConfig, current, best);
                }
            } else {
                cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong Data type " << endl;
                exit(-1);
            }
        } else {
            cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ <<": [ERROR:Breiman] Currently we have only Variance (for regression) as impurity function" << endl;
            exit(-1);
        }

    } else if (!y->is_ordinal) { // Xi is categorical, Y categorical, Brute force classification

        if((treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY) || (treeConfig.IMPURITY_FUNC == IMPURITY_GINI)
           || (treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR)) {
            if(Xi->data_type == ELEM_BOOL) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<bool, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<bool, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<bool, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<bool, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<bool, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<bool, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<bool, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_SHORT) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<short, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<short, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<short, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<short, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<short, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<short, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<short, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_INT) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<int, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<int, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<int, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<int, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<int, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<int, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<int, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_FLOAT) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<float, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<float, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<float, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<float, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<float, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<float, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<float, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_DOUBLE) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<double, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<double, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<double, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<double, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<double, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<double, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<double, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_CHAR) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<char, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<char, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<char, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<char, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<char, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<char, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<char, string>(dataset, start, treeConfig, current, best);
                }
            } else if(Xi->data_type == ELEM_STRING) {
                if(y->data_type == ELEM_BOOL) {
                    update_best_split_brute_force<string, bool>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_SHORT) {
                    update_best_split_brute_force<string, short>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_INT) {
                    update_best_split_brute_force<string, int>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_FLOAT) {
                    update_best_split_brute_force<string, float>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_DOUBLE) {
                    update_best_split_brute_force<string, double>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_CHAR) {
                    update_best_split_brute_force<string, char>(dataset, start, treeConfig, current, best);
                } else if(y->data_type == ELEM_STRING) {
                    update_best_split_brute_force<string, string>(dataset, start, treeConfig, current, best);
                }
            }

        } else {
            cout<<"File: " << __FILE__<< ", Line " << __LINE__<< ": [ERROR] Currently we have only Entropy, Gini, Classification Error as impurity function for classification (Xi is categorical)" << endl;
            exit(-1);
        }
    }
}


// =========== newly added for ExtraTrees ==========
std::random_device rd;
std::mt19937 mt(rd());
double sample(double left, double right)
{
	std::uniform_real_distribution<double> dist(left, right);
	return dist(mt);
}
// =========== newly added for ExtraTrees ==========
template <class Tx, class Ty>
void set_random_pos_regression(vector<id_value<Tx, Ty>> &pairs, Best &current, double split_val) {
    //output:
    //best_impurity_improvement
    //best offest: the chosen split position from "start"
    // feature_value (Xi) of type numerical expected (ref: FEATURE_THRESHOLD)
#ifdef ASSERT
    assert(typeid(Tx) != typeid(char));
    assert(typeid(Tx) != typeid(string));
    assert(typeid(Tx) != typeid(bool));
#endif
    current.impurity_improvement = -DBL_MAX;
    var_sketch left;
    var_sketch right;
    set_sketch<Tx, Ty>(right, pairs);
    for (size_t i = 0; i < pairs.size() - 1; i++) {
        id_value<Tx, Ty> &pair = pairs[i];
        update_sketches<Ty>(left, right, pair.y_value); //update with the current element
        if(pairs[i+1].x_value == pair.x_value) continue;
        if(pairs[i+1].x_value >= split_val) // this is the only place to split and check
        {
        	double var_L = var(left);
			double var_R = var(right);
			double impurity_improvement = impurity_improvement_fn(var_L, var_R, left.sum_1, right.sum_1);
			if (current.impurity_improvement < impurity_improvement) {
				current.impurity_improvement = impurity_improvement;
				current.offset = i; //like (it - start) == (i - 0);
			}
			return; // once reach the place, return immediately
        }
    }
}
// =========== newly added for ExtraTrees ==========
template <class Tx, class Ty>
void set_random_pos_classification(vector<id_value<Tx, Ty>> &pairs, Best &current, TreeConfig &treeConfig, double split_val) {
    //output:
    //Best current
#ifdef ASSERT
    assert(typeid(Tx) != typeid(char));
    assert(typeid(Tx) != typeid(string));
#endif
    //ToDo :: ordered string is rare, can implement if needed
    current.impurity_improvement = -DBL_MAX;
    classify_sketch<Ty> left;
    classify_sketch<Ty> right;
    set_sketch<Tx, Ty>(right, pairs);
    for (size_t i = 0; i < pairs.size() - 1; i++) {
        id_value<Tx, Ty> &pair = pairs[i];
        update_sketches<Ty>(left, right, pair.y_value); //update with the current element
        if(pairs[i+1].x_value == pair.x_value) continue;
        if(pairs[i+1].x_value >= split_val) // this is the only place to split and check
        {
        	//compute pi
			vector<double> left_tgt_ratio;
			vector<double> right_tgt_ratio;
			left.set_ratio(left_tgt_ratio);
			right.set_ratio(right_tgt_ratio);
			double left_impurity;
			double right_impurity;
			if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY) {
				left_impurity = entropy(left_tgt_ratio);
				right_impurity = entropy(right_tgt_ratio);
			} else if (treeConfig.IMPURITY_FUNC == IMPURITY_GINI) {
				left_impurity = gini(left_tgt_ratio);
				right_impurity = gini(right_tgt_ratio);
			} else if (treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR) {
				left_impurity = classification_error(left_tgt_ratio);
				right_impurity = classification_error(right_tgt_ratio);
			} else {
				cout<< "[ERROR] for classification we only consider entropy, "
						"gini and classification error as impurity function " << endl;
				exit(-1);
			}
			double impurity_improvement = impurity_improvement_fn(left_impurity, right_impurity, left.N, right.N);
			if (current.impurity_improvement < impurity_improvement) {
				current.impurity_improvement = impurity_improvement;
				current.offset = i; // (it - start) ==  (i - 0)
			}
			return; // once reach the place, return immediately
        }
    }
}
// =========== newly added for ExtraTrees ==========
template <class Tx, class Ty>
void update_random_split_ordinal_regression(Matrix &dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                                              Best &current, Best &best) {
    vector<id_value<Tx, Ty>> pairs;
    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);
    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = " << current.column_idx << " (all are missing)" << endl;
        return;
    }
    std::sort(pairs.begin(), pairs.end());
    double split_val = sample(boost::lexical_cast<double>(pairs.front().x_value), boost::lexical_cast<double>(pairs.back().x_value));
    set_random_pos_regression<Tx, Ty>(pairs, current, split_val);
    //---
    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.end = current.end;
    }
}
// =========== newly added for ExtraTrees ==========
template <class Tx, class Ty>
void update_random_split_ordinal_classification(Matrix &dataset, vector<size_t>::iterator start, TreeConfig &treeConfig,
                               Best &current, Best &best) {
    vector<id_value<Tx, Ty>> pairs;
    align_missing<Tx, Ty>(dataset, current.column_idx, start, current.end, pairs);
    if(start == current.end) { // all data values are missing
        if(!NO_WARNING) cout<< "[Warning] Skip the feature with index = " << current.column_idx << " (all are missing)" << endl;
        return;
    }
    std::sort(pairs.begin(), pairs.end());
    double split_val = sample(boost::lexical_cast<double>(pairs.front().x_value), boost::lexical_cast<double>(pairs.back().x_value));
    set_random_pos_classification<Tx, Ty>(pairs, current, treeConfig, split_val);
    //---
    if(current.impurity_improvement > best.impurity_improvement
       || (current.impurity_improvement == best.impurity_improvement
           && current.column_idx < best.column_idx)) {
        best.impurity_improvement = current.impurity_improvement;
        best.column_idx = current.column_idx;
        best.offset = current.offset;
        best.end = current.end;
    }
}
// =========== newly added for ExtraTrees ==========
void find_random_split_column(Matrix & dataset, vector<size_t>::iterator start, Best & current,
                            Best & best, TreeConfig & treeConfig) {
	Column* Xi = dataset.get_column(current.column_idx);
	Column* y = dataset.get_column(y_index);
	if(Xi->is_ordinal)
	{
		if(treeConfig.IMPURITY_FUNC == IMPURITY_VARIANCE) {
			if(Xi->data_type ==  ELEM_SHORT) {
				if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_regression<short, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_regression<short, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_regression<short, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_regression<short, double>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_INT) {
				if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_regression<int, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_regression<int, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_regression<int, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_regression<int, double>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_FLOAT) {
				if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_regression<float, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_regression<float, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_regression<float, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_regression<float, double>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_DOUBLE) {
				if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_regression<double, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_regression<double, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_regression<double, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_regression<double, double>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			}
		} else if (treeConfig.IMPURITY_FUNC == IMPURITY_ENTROPY
				   || treeConfig.IMPURITY_FUNC == IMPURITY_GINI
				   || treeConfig.IMPURITY_FUNC == IMPURITY_CLASSIFICATION_ERROR) {
			if(Xi->data_type ==  ELEM_SHORT) {
				if(y->data_type ==  ELEM_BOOL) {
					update_random_split_ordinal_classification<short, bool>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_classification<short, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_classification<short, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_classification<short, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_classification<short, double>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_CHAR) {
					update_random_split_ordinal_classification<short, char>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_STRING) {
					update_random_split_ordinal_classification<short, string>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_INT) {
				if(y->data_type ==  ELEM_BOOL) {
					update_random_split_ordinal_classification<int, bool>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_classification<int, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_classification<int, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_classification<int, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_classification<int, double>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_CHAR) {
					update_random_split_ordinal_classification<int, char>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_STRING) {
					update_random_split_ordinal_classification<int, string>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_FLOAT) {
				if(y->data_type ==  ELEM_BOOL) {
					update_random_split_ordinal_classification<float, bool>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_classification<float, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_classification<float, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_classification<float, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_classification<float, double>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_CHAR) {
					update_random_split_ordinal_classification<float, char>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_STRING) {
					update_random_split_ordinal_classification<float, string>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			} else if(Xi->data_type ==  ELEM_DOUBLE) {
				if(y->data_type ==  ELEM_BOOL) {
					update_random_split_ordinal_classification<double, bool>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_SHORT) {
					update_random_split_ordinal_classification<double, short>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_INT) {
					update_random_split_ordinal_classification<double, int>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_FLOAT) {
					update_random_split_ordinal_classification<double, float>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_DOUBLE) {
					update_random_split_ordinal_classification<double, double>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_CHAR) {
					update_random_split_ordinal_classification<double, char>(dataset, start, treeConfig, current, best);
				} else if(y->data_type ==  ELEM_STRING) {
					update_random_split_ordinal_classification<double, string>(dataset, start, treeConfig, current, best);
				} else {
					cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": Wrong type " << endl;
					exit(-1);
				}
			}
		}
	} else {
		cout<<"File: " << __FILE__<< ", Line " << __LINE__<< ": [ERROR] Xi is categorical, not supported for completely random trees" << endl;
		exit(-1);
	}
}


SplitResult* node_split(Matrix & dataset, vector<size_t>::iterator start, vector<size_t>::iterator &end,
                        vector<int> & cols, TreeConfig &treeConfig) {
    // need to pass end param as reference, value might change (align_missing(.) case)
    Column* y = dataset.get_column(y_index);

    if((treeConfig.IMPURITY_FUNC == IMPURITY_VARIANCE)
       && !(y->data_type >= ELEM_SHORT && y->data_type <= ELEM_DOUBLE)) {
        cout<<"File: " << __FILE__<< ", Line = " <<__LINE__ << ": [ERROR] output column should be numeric for regression " << endl;
        exit(-1);
    }

    Best best;

    for(size_t col_idx = 0; col_idx < cols.size(); col_idx++) {
        if(cols[col_idx] == y_index) continue; // skipping target_column from splitting

        Best current;
        current.column_idx = cols[col_idx];
        current.end = end;

        if(treeConfig.type == EXTRA_TREES) find_random_split_column(dataset, start, current, best, treeConfig);
        else find_best_split_column(dataset, start, current, best, treeConfig);

    }

    if(best.column_idx == -2) {
        SplitResult* result = new SplitResult;
        result->column_idx = -2;
        result->best_impurity = -DBL_MAX;

        return result;
    }

    Column* best_Xi = dataset.get_column(best.column_idx);

    if(best_Xi->data_type == ELEM_BOOL) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_wrapper<bool, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_wrapper<bool, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_wrapper<bool, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_wrapper<bool, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_wrapper<bool, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_wrapper<bool, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_wrapper<bool, string>(dataset, start, end, best);
        }  else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_SHORT) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_mid_int_wrapper<short, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_mid_int_wrapper<short, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_mid_int_wrapper<short, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_mid_int_wrapper<short, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_mid_int_wrapper<short, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_mid_int_wrapper<short, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_mid_int_wrapper<short, string>(dataset, start, end, best);
        }  else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_INT) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_mid_int_wrapper<int, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_mid_int_wrapper<int, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_mid_int_wrapper<int, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_mid_int_wrapper<int, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_mid_int_wrapper<int, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_mid_int_wrapper<int, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_mid_int_wrapper<int, string>(dataset, start, end, best);
        } else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_FLOAT) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_mid_wrapper<float, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_mid_wrapper<float, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_mid_wrapper<float, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_mid_wrapper<float, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_mid_wrapper<float, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_mid_wrapper<float, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_mid_wrapper<float, string>(dataset, start, end, best);
        } else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_DOUBLE) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_mid_wrapper<double, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_mid_wrapper<double, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_mid_wrapper<double, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_mid_wrapper<double, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_mid_wrapper<double, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_mid_wrapper<double, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_mid_wrapper<double, string>(dataset, start, end, best);
        } else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_CHAR) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_mid_int_wrapper<char, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_mid_int_wrapper<char, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_mid_int_wrapper<char, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_mid_int_wrapper<char, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_mid_int_wrapper<char, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_mid_int_wrapper<char, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_mid_int_wrapper<char, string>(dataset, start, end, best);
        } else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    } else if(best_Xi->data_type == ELEM_STRING) {
        if(y->data_type == ELEM_BOOL) {
            return create_split_result_wrapper<string, bool>(dataset, start, end, best);
        } else if(y->data_type == ELEM_SHORT) {
            return create_split_result_wrapper<string, short>(dataset, start, end, best);
        } else if(y->data_type == ELEM_INT) {
            return create_split_result_wrapper<string, int>(dataset, start, end, best);
        } else if(y->data_type == ELEM_FLOAT) {
            return create_split_result_wrapper<string, float>(dataset, start, end, best);
        } else if(y->data_type == ELEM_DOUBLE) {
            return create_split_result_wrapper<string, double>(dataset, start, end, best);
        } else if(y->data_type == ELEM_CHAR) {
            return create_split_result_wrapper<string, char>(dataset, start, end, best);
        } else if(y->data_type == ELEM_STRING) {
            return create_split_result_wrapper<string, string>(dataset, start, end, best);
        } else {
            cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    }
    else {
        cout<<"File: " << __FILE__<<", Line "<< __LINE__ << ": Type Impossible!" << endl;
        exit(-1);
    }
}


#endif
