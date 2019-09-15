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

#ifndef CRITERION_H
#define CRITERION_H

#include "column.h"
//################### impurity functions ###################
//[Classification] link of the concepts of (1) entropy, Gini index, or classification error, and (2) Information Gain (IG):
//http://www.bogotobogo.com/python/scikit-learn/scikt_machine_learning_Decision_Tree_Learning_Informatioin_Gain_IG_Impurity_Entropy_Gini_Classification_Error.php
//we pass entropy, Gini index, or classification error as a function pointer to IG

double entropy(vector<double> target_ratios) { //input to info_gain()
    double sum = 0.0;

    for (size_t i = 0; i < target_ratios.size(); i++) {
        if (target_ratios[i] != 0 && target_ratios[i] != 1) { // to deal with NaN
            sum += target_ratios[i] * log2(target_ratios[i]);
        }
    }

    return (sum == 0) ? sum : (-1) * sum;
}

double gini(vector<double> target_ratios) { //input to info_gain()
    double square_sum = 0.0;

    for (size_t i = 0; i < target_ratios.size(); i++) {
        double square = target_ratios[i] * target_ratios[i];
        square_sum += square;
    }

    return (1 - square_sum);
}

double classification_error(vector<double> target_ratios) { //input to info_gain()
    return (1 - *max_element(target_ratios.begin(), target_ratios.end()));
}

//subfunction of info_gain
template <class T>
double classify_impurity(vector<size_t>::iterator start,
                        vector<size_t>::iterator end,
                        Column* Y, double (*f)(vector<double>)) {//used by info_gain() below

    //group by categories, and get counts in each category
    map<T, size_t> category_count_map;
    for (auto index = start; index != end; index++) {
        T key;
        Y->get(*index, &key);

        auto it = category_count_map.find(key);
        if (it == category_count_map.end()) {
            category_count_map[key] = 1;
        } else {
            it->second++;
        }
    }

    size_t size = end - start;

    //get the counts to compute impurity
    vector<double> target_ratios;
    for (auto it = category_count_map.begin();
         it != category_count_map.end(); it++) {
        double ratio = (double) it->second / size;
        target_ratios.push_back(ratio);
    }

    //call impurity function
    return (*f)(target_ratios);
}

//wrapper of information_gain()
double classification_impurity(vector<size_t>::iterator start,
                 vector<size_t>::iterator end,
                 Column* Y, double (*f)(vector<double>)) {//impurity function for classification

    int data_type = Y->data_type;

    if(data_type == ELEM_BOOL){
        return classify_impurity<bool>(start, end, Y, f);
    } else if(data_type == ELEM_SHORT){
        return classify_impurity<short>(start, end, Y, f);
    } else if(data_type == ELEM_INT){
        return classify_impurity<int>(start, end, Y, f);
    } else if(data_type == ELEM_FLOAT){
        return classify_impurity<float>(start, end, Y, f);
    } else if(data_type == ELEM_DOUBLE){
        return classify_impurity<double>(start, end, Y, f);
    } else if(data_type == ELEM_CHAR){
        return classify_impurity<char>(start, end, Y, f);
    } else if(data_type == ELEM_STRING){
        return classify_impurity<string>(start, end, Y, f);
    }
    else{
        cout<< "[ERROR] Unsupported type in info_gain() ! " << endl;
        exit(-1);
    }
}

//subfunction of var()
template<class T>
double variance(vector<size_t>::iterator start,
                vector<size_t>::iterator end,
                Column* Y) {//used by var() below
    //compute mean
    T sum = 0;
    for (auto index = start; index != end; index++) {
        T value;
        Y->get(*index, &value);
        sum += value;
    }

    size_t size = end - start;

    double mean = ((double)sum / size);
    //compute sum of squares
    double numerator_sum = 0.0;
    for (auto index = start; index != end; index++) {
        T value;
        Y->get(*index, &value);
        double x = value - mean;
        numerator_sum += (x * x);
    }
    //compute variance
    return (numerator_sum / size);
}

double var(vector<size_t>::iterator start,
           vector<size_t>::iterator end,
           Column* Y) {//impurity function for regression

    int data_type = Y->data_type;
    if(data_type == ELEM_SHORT){
        return variance<short>(start, end,  Y);
    } else if(data_type == ELEM_INT){
        return variance<int>(start, end, Y);
    } else if(data_type == ELEM_FLOAT){
        return variance<float>(start, end, Y);
    } else if(data_type == ELEM_DOUBLE){
        return variance<double>(start, end, Y);
    } else{
        cout<< "[ERROR] Column type is not numeric when calculating variance ! " << endl;
        exit(-1);
    }
}

#endif