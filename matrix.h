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

#ifndef MATRIX_H
#define MATRIX_H

#include "column.h"

using namespace std;

class Matrix {
	//"new Column" + append;
	//"delete" is done by destructor
public:
	vector<Column*> col;
	size_t size = 0;

	~Matrix()
	{
		for(size_t i=0; i<col.size(); i++) delete col[i];
	}

	inline void append(Column* column)
	{
		col.push_back(column);
		size++;
	}

	inline Column* get_column(int col_idx)
	{
		return col[col_idx];
	}

	inline void get_entry(size_t i, size_t j, void* output) { //row i, column j
		return col[j]->get(i, output);
	}

    void copy_meta(Matrix &dest) {
        for(size_t i = 0; i < col.size(); i++) {
            Column* sc = get_column(i);
            Column* column = sc->copy_meta();
            dest.append(column);
        }
    }

    //=== for debugging
    void print()
    {
        if(col.size() == 0)
        {//empty matrix
            return;
        }

        size_t num_rows = col[0]->size;
        size_t num_cols = col.size();

        for(size_t i = 0; i < num_rows; i++){
            for(size_t j = 0; j < num_cols; j++){
                int type = col[j]->data_type;

                if(type == ELEM_BOOL)
                {
                    bool entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_SHORT)
                {
                    short entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_INT)
                {
                    int entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_FLOAT)
                {
                    float entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_DOUBLE)
                {
                    double entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_CHAR)
                {
                    char entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
                else if(type == ELEM_STRING)
                {
                    string entry;
                    get_entry(i, j, &entry);
                    if(j > 0) cout << ",";
                    cout << entry;
                }
            }
            cout << endl;
        }
    }
    
    void print_column(size_t j)
	{
		if(col.size() < j)
		{//no such column
			return;
		}

		size_t num_rows = col[j]->size;

		for(size_t i = 0; i < num_rows; i++){
			int type = col[j]->data_type;
			if(type == ELEM_BOOL)
			{
				bool entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_SHORT)
			{
				short entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_INT)
			{
				int entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_FLOAT)
			{
				float entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_DOUBLE)
			{
				double entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_CHAR)
			{
				char entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			else if(type == ELEM_STRING)
			{
				string entry;
				get_entry(i, j, &entry);
				cout << entry;
			}
			cout << endl;
		}
	}

    void print_meta()
    {
        if(col.size() == 0)
        {//empty matrix
            return;
        }
        
        size_t num_cols = col.size();
        
        for(size_t i = 0; i < num_cols; i++){
            int type = col[i]->data_type;
            if(type == ELEM_BOOL) cout << i <<": bool, ";
            else if(type == ELEM_SHORT) cout << i <<": short, ";
            else if(type == ELEM_INT) cout << i <<": int, ";
            else if(type == ELEM_FLOAT) cout << i <<": float, ";
            else if(type == ELEM_DOUBLE) cout << i <<": double, ";
            else if(type == ELEM_CHAR) cout << i <<": char, ";
            else if(type == ELEM_STRING) cout << i <<": string, ";
            //----
            if(col[i]->is_ordinal) cout << "ordinal, ";
            else cout << "categorical, ";
            //----
            if(col[i]->is_dense) cout << "dense, ";
            else
            {
                cout << "sparse (default = ";
                if(type == ELEM_BOOL) cout << ((SparseColumn<bool> *)col[i])->default_value;
                else if(type == ELEM_SHORT) cout << ((SparseColumn<short> *)col[i])->default_value;
                else if(type == ELEM_INT) cout << ((SparseColumn<int> *)col[i])->default_value;
                else if(type == ELEM_FLOAT) cout << ((SparseColumn<float> *)col[i])->default_value;
                else if(type == ELEM_DOUBLE) cout << ((SparseColumn<double> *)col[i])->default_value;
                else if(type == ELEM_CHAR) cout << ((SparseColumn<char> *)col[i])->default_value;
                else if(type == ELEM_STRING) cout << ((SparseColumn<string> *)col[i])->default_value;
                cout <<"), ";
            }
            cout<<"missing_value = ";
            col[i]->print_missing_value();
            cout<<endl;
        }
    }
};

#endif
