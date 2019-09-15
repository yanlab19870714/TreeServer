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

#ifndef CSV_H
#define CSV_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include "matrix.h"
#include "config.h"

using namespace std;

//metadata format:
//- it's a csv file
//- one line per attribute

//format:
//data_type, is_ordinal, is_dense, (missing_value, default_value)

//Default value is marked with **
//- data_type = bool / short / int / float / double** / char / string
//- is_ordinal = true** / false
//- is_dense = true** / false
//- missing_value = user-specified, depending on data_type
//---- bool: default = 2
//---- short: default = SHRT_MIN
//---- int: default = INT_MIN
//---- float: default = FLT_MIN
//---- double: default = DBL_MIN
//---- char: default = 0
//---- string: default = ""
//- default_value = user-specified, depending on data_type
//---- bool: default = false
//---- short: default = 0
//---- int: default = 0
//---- float: default = 0.0
//---- double: default = 0.0
//---- char: default = 0
//---- string: default = ""

Column * create_dense_column(int data_type, bool is_ordinal)
{
    Column * column;
    if(data_type == ELEM_BOOL) column = new DenseColumn<bool>(ELEM_BOOL, is_ordinal);
    else if(data_type == ELEM_SHORT) column = new DenseColumn<short>(ELEM_SHORT, is_ordinal);
    else if(data_type == ELEM_INT) column = new DenseColumn<int>(ELEM_INT, is_ordinal);
    else if(data_type == ELEM_FLOAT) column = new DenseColumn<float>(ELEM_FLOAT, is_ordinal);
    else if(data_type == ELEM_DOUBLE) column = new DenseColumn<double>(ELEM_DOUBLE, is_ordinal);
    else if(data_type == ELEM_CHAR) column = new DenseColumn<char>(ELEM_CHAR, is_ordinal);
    else if(data_type == ELEM_STRING) column = new DenseColumn<string>(ELEM_STRING, is_ordinal);
    else {
		cout<<"File: " << __FILE__ <<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
		exit(-1);
	}
    return column;
}

Column * create_sparse_column(int data_type, bool is_ordinal)
{
    Column * column;
    if(data_type == ELEM_BOOL) column = new SparseColumn<bool>(ELEM_BOOL, is_ordinal, false);
    else if(data_type == ELEM_SHORT) column = new SparseColumn<short>(ELEM_SHORT, is_ordinal, 0);
    else if(data_type == ELEM_INT) column = new SparseColumn<int>(ELEM_INT, is_ordinal, 0);
    else if(data_type == ELEM_FLOAT) column = new SparseColumn<float>(ELEM_FLOAT, is_ordinal, 0.0);
    else if(data_type == ELEM_DOUBLE) column = new SparseColumn<double>(ELEM_DOUBLE, is_ordinal, 0.0);
    else if(data_type == ELEM_CHAR) column = new SparseColumn<char>(ELEM_CHAR, is_ordinal, 0);
    else if(data_type == ELEM_STRING) column = new SparseColumn<string>(ELEM_STRING, is_ordinal, "");
    else {
		cout<<"File: " << __FILE__ <<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
		exit(-1);
	}
    return column;
}

Column * create_sparse_column(int data_type, bool is_ordinal, char * default_value)
{
    Column * column;
    if(data_type == ELEM_BOOL)
    {
        SparseColumn<bool> * col = new SparseColumn<bool>(ELEM_BOOL, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_SHORT)
    {
        SparseColumn<short> * col = new SparseColumn<short>(ELEM_SHORT, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_INT)
    {
        SparseColumn<int> * col = new SparseColumn<int>(ELEM_INT, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_FLOAT)
    {
        SparseColumn<float> * col = new SparseColumn<float>(ELEM_FLOAT, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_DOUBLE)
    {
        SparseColumn<double> * col = new SparseColumn<double>(ELEM_DOUBLE, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_CHAR)
    {
        SparseColumn<char> * col = new SparseColumn<char>(ELEM_CHAR, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else if(data_type == ELEM_STRING)
    {
        SparseColumn<string> * col = new SparseColumn<string>(ELEM_STRING, is_ordinal);
        col->set_default_value(default_value);
        column = col;
    }
    else {
		cout<<"File: " << __FILE__ <<", Line: "<< __LINE__<<", [ERROR] we don't consider any other attribute type"<<endl;
		exit(-1);
	}
    return column;
}

void load_meta(const char* metafile, Matrix &X) {
    //load metadata to vectors in memory
    ifstream fin(metafile);
    //------
    char line[10000]; //a line in metafile
    string value; //an item separated by comma
    //------
    //mandatory files??fields
    int data_type;
    bool is_ordinal, is_dense;
    //tokenizer
    char * prev; //last position of split + 1
    char * pch; //current position of split
    //------
    while(fin.getline(line, 10000))
    {
        if(line[0] == '#') continue; //skip lines starting with '#'
        prev = line; //tokenizer: move to line head
        //============== read in field 1 data_type ==============
        pch = (char*) memchr (prev, ',', strlen(prev));
        if(pch == NULL)
        {
            cout<<"[ERROR] A line with no data_type field is detected !!!"<<endl;
            exit(-1);
        }
        else *pch = '\0'; //tokenizer: cut
        if(strcmp(prev, "bool") == 0) data_type = ELEM_BOOL;
        else if(strcmp(prev, "short") == 0) data_type = ELEM_SHORT;
        else if(strcmp(prev, "int") == 0) data_type = ELEM_INT;
        else if(strcmp(prev, "float") == 0) data_type = ELEM_FLOAT;
        else if(strcmp(prev, "double") == 0) data_type = ELEM_DOUBLE;
        else if(strcmp(prev, "char") == 0) data_type = ELEM_CHAR;
        else if(strcmp(prev, "string") == 0) data_type = ELEM_STRING;
        else
        {
            cout<<"[ERROR] Unsupported data type: "<<prev<<" !!!"<<endl;
            exit(-1);
        }
         //============== read in field 2 is_ordinal ==============
        prev = pch + 1; //tokenizer: move prev
        pch = (char*) memchr (prev, ',', strlen(prev)); //tokenizer: move to line head
        if(pch == NULL)
        {
            cout<<"[ERROR] A line with no is_ordinal field is detected !!!"<<endl;
            exit(-1);
        }
        else *pch = '\0'; //tokenizer: cut
        if(strcmp(prev, "true") == 0) is_ordinal = true;
        else if(strcmp(prev, "false") == 0) is_ordinal = false;
        else if(strcmp(prev, "") == 0) is_ordinal = true;
        else
        {
            cout<<"[ERROR] Unsupported is_ordinal field: "<<prev<<" !!!"<<endl;
            exit(-1);
        }
        //============== read in field 3 is_dense ==============
        prev = pch + 1; //tokenizer: move prev
        pch = (char*) memchr (prev, ',', strlen(prev)); //tokenizer: move to line head
        if(pch == NULL) //end of line
        {
            if(strlen(prev) == 0)//default dense
            {
                Column * column = create_dense_column(data_type, is_ordinal);
                X.append(column);
            }
            else
            {
                if(strcmp(prev, "true") == 0) is_dense = true;
                else if(strcmp(prev, "false") == 0) is_dense = false;
                else if(strcmp(prev, "") == 0) is_dense = true;//default dense
                else
                {
                    cout<<"[ERROR] Unsupported is_dense field: "<<prev<<" !!!"<<endl;
                    exit(-1);
                }
                Column * column;
                if(is_dense) column = create_dense_column(data_type, is_ordinal);
                else column = create_sparse_column(data_type, is_ordinal);
                X.append(column);
            }
            continue;
        }
        else
        {
            *pch = '\0'; //tokenizer: cut
            if(strcmp(prev, "true") == 0) is_dense = true;
            else if(strcmp(prev, "false") == 0) is_dense = false;
            else if(strcmp(prev, "") == 0) is_dense = true;
            else
            {
                cout<<"[ERROR] Unsupported is_dense field: "<<prev<<" !!!"<<endl;
                exit(-1);
            }
            //============== read in field 4 missing value ==============
            prev = pch + 1; //tokenizer: move prev
            pch = (char*) memchr (prev, ',', strlen(prev)); //tokenizer: move to line head
            if(pch == NULL) //end of line
            {
                Column * column;
                if(is_dense) column = create_dense_column(data_type, is_ordinal);
                else column = create_sparse_column(data_type, is_ordinal);
                X.append(column);
                if(strlen(prev) > 0) column->set_missing_value(prev);
                continue;
            }
            else
            {
                *pch = '\0'; //tokenizer: cut
                string miss = prev;
                //============== read in field 5 default sparse value ==============
                prev = pch + 1; //tokenizer: move prev
                pch = (char*) memchr (prev, ',', strlen(prev)); //tokenizer: move to line head
                if(pch != NULL) *pch = '\0'; //tokenizer: cut
                Column * column;
                if(is_dense)
                {
                    column = create_dense_column(data_type, is_ordinal);
                }
                else
                {
                    if(strlen(prev) > 0) column = create_sparse_column(data_type, is_ordinal, prev);
                    else column = create_sparse_column(data_type, is_ordinal);
                }
                if(miss.length() > 0) column->set_missing_value(miss);
                X.append(column);
            }
        }
    }

    fin.close();
}

void load_data(const char* datafile, Matrix &X) //if a value is empty string, use default value (for sparse column)
{
    ifstream fin(datafile);
    string line, value;
    //------
    int row_idx = 0;
    while(getline(fin, line))
    {
        istringstream sin(line);
        int col_idx = 0;
        while(getline(sin, value, ','))
        {
            if(value != "") // this means missing value, the column is assumed to be sparse
            {
                Column * column = X.get_column(col_idx);
                if(column->data_type == ELEM_BOOL)
                {
                    bool temp;
                    if(strcmp(value.c_str(), "true") == 0) temp = true;
                    else if(strcmp(value.c_str(), "false") == 0) temp = false;
                    else
                    {
                        cout << "[CSV ERROR] Boolean field gets unrecognized value: " << value << endl;
                        exit(-1);
                    }
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_SHORT)
                {
                    short temp = (short) atoi(value.c_str());
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_INT)
                {
                    int temp = atoi(value.c_str());
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_FLOAT)
                {
                    float temp = stof(value);
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_DOUBLE)
                {
                    double temp = stod(value);
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_CHAR)
                {
                    char temp = value.c_str()[0];
                    column->append(row_idx, &temp);
                }
                else if(column->data_type == ELEM_STRING)
                {
                    column->append(row_idx, &value);
                }
            } else {
                cout<< "[CSV ERROR] Sparse Data Error" << endl;
                exit(-1);
            }
            col_idx++;
        }
        row_idx++;
    }
    fin.close();
    //------
    vector<Column*> & columns = X.col;
    for(size_t i=0; i<columns.size(); i++)
    	if(columns[i]->is_dense == false)
    		columns[i]->finish(row_idx);
}

void load_csv(const char* datafile, const char* metafile, Matrix &X)
{
    load_meta(metafile, X);
    auto start = chrono::system_clock::now();
    load_data(datafile, X);
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout<< "Load time = " << elapsed.count() << " milliseconds" << endl;
}

//=============================================================================

// shuffling the indices and collecting number of split in a vector
// input : sample_size, n_split
// output: vector<vector<size_t>> indices, n_split number of data sets
void shuffle_indices(size_t sample_size, size_t n_split, vector<vector<size_t>> &shuffled_indices) {
    vector<size_t> indices;// collecting all the indices
    for(size_t i = 0; i < sample_size; i++) {
        indices.push_back(i);
    }
    
    random_shuffle(indices.begin(), indices.end()); // shuffling the indices
    
    size_t chunk = sample_size / n_split;
    size_t remain = sample_size % n_split;
    size_t offset = 0;
    
    for(size_t i = 0; i < remain; i++) {
        vector<size_t> temp;
        for(size_t j = 0; j <= chunk; j++) {
            temp.push_back(indices[offset + j]);
        }
        offset += chunk + 1;
        shuffled_indices.push_back(temp);
    }
    
    for(size_t i = remain; i < n_split; i++) {
        vector<size_t> temp;
        for(size_t j = 0; j < chunk; j++) {
            temp.push_back(indices[offset + j]);
        }
        offset += chunk;
        shuffled_indices.push_back(temp);
    }
}

// copying data (indices) from  data_set to dest Matrix
void set_data(Matrix &data_set, vector<size_t> &indices, Matrix &dest) {
    for(size_t i = 0; i < data_set.size; i++) { // traverse through columns
        Column* sourceColumn = data_set.get_column(i);
        Column* dest_column = dest.get_column(i);
        for(size_t j = 0; j < indices.size(); j++) {
            if(sourceColumn->data_type == ELEM_BOOL) {
                bool value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_SHORT) {
                short value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_INT) {
                int value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_FLOAT) {
                float value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_DOUBLE) {
                double value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_CHAR) {
                char value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            } else if(sourceColumn->data_type == ELEM_STRING) {
                string value;
                sourceColumn->get(indices[j], &value);
                dest_column->append(j, &value);
            }
        }
    }
}

// split data_set to n_split set
// input : Matrix &data_set, n_split
// output : vector<Matrix> collection of splitted sets
void split_data(Matrix &data_set, size_t n_split, vector<Matrix> &split_set) {
    size_t n_sample = data_set.get_column(0)->size;
    split_set.reserve(n_split);
    
    //create split datasets, and assign meta to them
    Matrix dummy;
    for(size_t i = 0; i < n_split; i++) {
        split_set.push_back(dummy);
        Matrix & dest = split_set.back();
        data_set.copy_meta(dest);
    }
    
    vector<vector<size_t>> shuffled_indices;
    shuffle_indices(n_sample, n_split, shuffled_indices);
    
    for(size_t i = 0; i < shuffled_indices.size(); i++) {
        Matrix &dest = split_set[i];
        set_data(data_set, shuffled_indices[i], dest);
    }
}

// wrapper function to load data, meta from csv file and split data to n_split chunks
void load_data(const char* datafile, const char* metafile, size_t n_split,
               vector<Matrix> &split_set, size_t &sample_size) {
    Matrix data_set;
    load_csv(datafile, metafile, data_set);
    sample_size = data_set.get_column(0)->size;
    split_data(data_set, n_split, split_set);
}

#endif
