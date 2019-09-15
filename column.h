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

#ifndef COLUMN_H
#define COLUMN_H

#include <vector>
#include <sstream>
#include <cassert>
#include <algorithm> // std::lower_bound (binary search with position returned)
#include <map>
#include <climits>
#include <cfloat>
#include "global.h"

using namespace std;

class Column { //one "whole column" of (1) a data table or (2) its partition (subset)
public:
	size_t size; //number of rows/elements in the column-vector
	int data_type; //what is the element type <T> in the column subclass
	bool is_ordinal; //whether the values have order; this is used to decide the node-splitting algorithm
    bool is_dense; // whether the column is dense or sparse
	virtual void get(size_t row_idx, void* output) = 0; //returns a pointer to the row-th element in the column
	virtual void append(size_t row_idx, void* value) = 0; //append a value to the column
	virtual void finish(size_t size) = 0; //number of elements
	//should never be called from the DenseColumn
    //should be called from SparseColumn, after a series of append:
    //- append, append, ..., append, finish
    //reason: the last appended elements may not be the last element

    // param: row_idx
    // output: string value of column[row_idx]
    string string_value(size_t row_idx) { //ToDo::
        int type = this->data_type;

        if(type == ELEM_BOOL){
            bool value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_SHORT){
            short value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_INT){
            int value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_FLOAT){
            float value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_DOUBLE){
            double value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_CHAR){
            char value;
            this->get(row_idx, &value);
            return boost::lexical_cast<string>(value);
        } else if(type == ELEM_STRING){
            string value;
            this->get(row_idx, &value);
            return value;
        } else
        {
            cout<<"Line "<< __LINE__ << ": [ERROR] unknown data type" << endl;
            exit(-1);
        }
    }

    template<class T>
    struct ColumnSort_comparator
    {
        Column* col;

        ColumnSort_comparator(Column* column)
        {
            col = column;
        }

        bool operator()(const size_t row1, const size_t row2)
        {
            T entry1;
            col->get(row1, &entry1);

            T entry2;
            col->get(row2, &entry2);

            return entry1 < entry2;
        }
    };

    template <class T>
    inline void row_sort(vector<size_t>::iterator it1,  vector<size_t>::iterator it2) { //it2 is exclusive
        ColumnSort_comparator<T> comp(this);
        sort(it1, it2, comp);
    }

    void row_sort(vector<size_t>::iterator it1, vector<size_t>::iterator it2) { //it2 is exclusive
#ifdef ASSERT
        assert(is_ordinal);
#endif

        if(data_type == ELEM_BOOL)
        {
            row_sort<bool>(it1, it2);
        }
        else if(data_type == ELEM_SHORT)
        {
            row_sort<short>(it1, it2);
        }
        else if(data_type == ELEM_INT)
        {
            row_sort<int>(it1, it2);
        }
        else if(data_type == ELEM_FLOAT)
        {
            row_sort<float>(it1, it2);
        }
        else if(data_type == ELEM_DOUBLE)
        {
            row_sort<double>(it1, it2);
        }
        else if(data_type == ELEM_CHAR)
        {
            row_sort<char>(it1, it2);
        }
        else if(data_type == ELEM_STRING)
        {
            row_sort<string>(it1, it2);
        }
    }

	virtual void print_missing_value() = 0; //used by csv.h load_meta()

	virtual void set_missing_value(string value) = 0; //used by csv.h load_meta()

	virtual Column* copy_meta() = 0;

	virtual ~Column() {};
};

template<class T>
class DenseColumn: public Column {
public:
	vector<T> elem; //dense vector of elements
	T missing_value;

	virtual ~DenseColumn() {};

	DenseColumn(int data_type, bool is_ordinal)
	{
		size = 0;
		this->data_type = data_type;
		this->is_ordinal = is_ordinal;
		is_dense = true;
		//------
		//set missing value to default ones
		if(data_type == ELEM_BOOL) *(bool*)&missing_value = 2;
		else if(data_type == ELEM_SHORT) *(short*)&missing_value = SHRT_MIN;
		else if(data_type == ELEM_INT) *(int*)&missing_value = INT_MIN;
		else if(data_type == ELEM_FLOAT) *(float*)&missing_value = FLT_MIN;
		else if(data_type == ELEM_DOUBLE) *(double*)&missing_value = DBL_MIN;
		else if(data_type == ELEM_CHAR) *(char*)&missing_value = 0;
		else if(data_type == ELEM_STRING) *(string*)&missing_value = "";
	}

	virtual void get(size_t row_idx, void* output) {
#ifdef ASSERT
		assert(row_idx >=0 && row_idx < elem.size()); ///### debugging (array boundary check), can be removed in real deployment
#endif
		T* ref = (T*) output;
		(*ref) = elem[row_idx];
	}

	virtual void append(size_t row_idx, void* value) { //needs to set row_idx properly; to be consistent with the base API
#ifdef ASSERT
		assert(elem.size() == row_idx); ///### debugging (row_idx check), can be removed in real deployment
#endif
		T* val = (T*) value;
		elem.push_back(*val);
		size++;
	}

	virtual void finish(size_t size) {
		cout<< "[ERROR] DenseColumn::finish(size) should not be called" << endl;
		exit(-1);
	}

	virtual void set_missing_value(string value)
	{
		missing_value = boost::lexical_cast<T>(value);
	}

	virtual void print_missing_value()
	{
		cout << missing_value;
	}

	virtual Column* copy_meta() {
		DenseColumn<T>* column = new DenseColumn<T>(data_type, is_ordinal);
		column->missing_value = missing_value;
		return column;
	}
};

template<class T>
class SparseColumn: public Column {
public:
	vector<T> elem;
	vector<size_t> row_index; //ordered
	T default_value; //value of an unspecified element
	T missing_value;

	virtual ~SparseColumn() {};

	SparseColumn(int data_type, bool is_ordinal)
	{
		this->data_type = data_type;
		this->is_ordinal = is_ordinal;
		is_dense = false;
		//------
		//set missing value to default ones
		if(data_type == ELEM_BOOL) *(bool*)&missing_value = 2;
		else if(data_type == ELEM_SHORT) *(short*)&missing_value = SHRT_MIN;
		else if(data_type == ELEM_INT) *(int*)&missing_value = INT_MIN;
		else if(data_type == ELEM_FLOAT) *(float*)&missing_value = FLT_MIN;
		else if(data_type == ELEM_DOUBLE) *(double*)&missing_value = DBL_MIN;
		else if(data_type == ELEM_CHAR) *(char*)&missing_value = 0;
		else if(data_type == ELEM_STRING) *(string*)&missing_value = "";
		//------
		//set default_value to default ones
		if(data_type == ELEM_BOOL) *(bool*)&default_value = false;
		else if(data_type == ELEM_SHORT) *(short*)&default_value = 0;
		else if(data_type == ELEM_INT) *(int*)&default_value = 0;
		else if(data_type == ELEM_FLOAT) *(float*)&default_value = 0.0;
		else if(data_type == ELEM_DOUBLE) *(double*)&default_value = 0.0;
		else if(data_type == ELEM_CHAR) *(char*)&default_value = 0;
		else if(data_type == ELEM_STRING) *(string*)&default_value = "";
	}

	SparseColumn(int data_type, bool is_ordinal, T default_value)
	{
		this->data_type = data_type;
		this->is_ordinal = is_ordinal;
		this->default_value = default_value;
		is_dense = false;
		//------
		//set missing value to default ones
		if(data_type == ELEM_BOOL) *(bool*)&missing_value = 2;
		else if(data_type == ELEM_SHORT) *(short*)&missing_value = SHRT_MIN;
		else if(data_type == ELEM_INT) *(int*)&missing_value = INT_MIN;
		else if(data_type == ELEM_FLOAT) *(float*)&missing_value = FLT_MIN;
		else if(data_type == ELEM_DOUBLE) *(double*)&missing_value = DBL_MIN;
		else if(data_type == ELEM_CHAR) *(char*)&missing_value = 0;
		else if(data_type == ELEM_STRING) *(string*)&missing_value = "";
		//------
		this->default_value = default_value;
	}

	virtual void get(size_t row_idx, void* output) {
		//binary search
        T* ref = (T*) output;
		auto it = lower_bound(row_index.begin(), row_index.end(), row_idx); //binary search
		if(it == row_index.end())
		{
			(*ref) = default_value; //no such element recorded, return default value
		}
		else if(*it > row_idx)
		{
			(*ref) = default_value; //no such element recorded, return default value
		}
		else //*it == row_idx
		{
			size_t pos = distance(row_index.begin(), it);
            (*ref) = elem[pos];
		}
	}

	virtual void append(size_t row_idx, void* value) {
        row_index.push_back(row_idx);
		T* val = (T*) value;
		elem.push_back(*val);
	}

	virtual void finish(size_t size) {
		this->size = size;
	}

	virtual void print_missing_value()
	{
		cout << missing_value;
	}

	virtual void set_missing_value(string value)
	{
		missing_value = boost::lexical_cast<T>(value);
	}

	void set_default_value(string value) //used by csv metadata load
	{
		default_value = boost::lexical_cast<T>(value);
	}

	virtual Column* copy_meta() {
		SparseColumn<T>* column = new SparseColumn<T>(data_type, is_ordinal);
		column->missing_value = missing_value;
		column->default_value = default_value;
        return column;
	}
};


#endif
