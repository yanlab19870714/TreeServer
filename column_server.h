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

#ifndef COLUMN_SERVER_H
#define COLUMN_SERVER_H

#include "csv.h"
#include "ydhdfs.h"

#include <fstream>

using namespace std;

class ColumnServer { //one "whole column" of (1) a data table or (2) its partition (subset)
public:
	Matrix X; //columns

	void load_meta(const char* metafile) {
		hdfsFS fs = getHdfsFS();
		hdfsFile in = getRHandle(metafile, fs);
		LineReader reader(fs, in);
		string value; //an item separated by comma
		//mandatory fields
		int data_type;
		bool is_ordinal, is_dense;
		//tokenizer
		char * prev; //last position of split + 1
		char * pch; //current position of split
		while(true)
		{
			reader.readLine();
			if (!reader.eof())
			{
				char * line = reader.getLine();
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
			else break;
		}
		hdfsCloseFile(fs, in);
		hdfsDisconnect(fs);
	}

	void load_local_meta(const char* metafile) {
		ifstream in(metafile);
		string value; //an item separated by comma
		//mandatory fields
		int data_type;
		bool is_ordinal, is_dense;
		//tokenizer
		char * prev; //last position of split + 1
		char * pch; //current position of split
		string line_str;
		while(getline(fin, line_str))
		{
			char* line = line_str.c_str();
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
		in.close();
	}

	void load_column(const char* dir, int col_idx) { //pass in directory name & column ID
		Column * column = X.get_column(col_idx);
		//------
		hdfsFS fs = getHdfsFS();
		string path = dir;
		path += "/column_" + to_string(col_idx);
		hdfsFile in = getRHandle(path.c_str(), fs);
		LineReader reader(fs, in);
		size_t row_idx = 0;
		while(true)
		{
			reader.readLine();
			if (!reader.eof())
			{
				char * value = reader.getLine();
				if(column->data_type == ELEM_BOOL)
				{
					bool temp;
					if(strcmp(value, "true") == 0) temp = true;
					else if(strcmp(value, "false") == 0) temp = false;
					else
					{
						cout << "[CSV ERROR] Boolean field gets unrecognized value: " << value << endl;
						exit(-1);
					}
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_SHORT)
				{
					short temp = (short) atoi(value);
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_INT)
				{
					int temp = atoi(value);
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_FLOAT)
				{
					float temp = stof(string(value));
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_DOUBLE)
				{
					double temp = stod(string(value));
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_CHAR)
				{
					char temp = value[0];
					column->append(row_idx, &temp);
				}
				else if(column->data_type == ELEM_STRING)
				{
					string val(value);
					column->append(row_idx, &val);
				}
				else
				{
					cout<< "File = " << __FILE__ << ", Line = " << __LINE__
					            << ": unknow data type" << endl;
					exit(-1);
				}
				row_idx++;
			}
			else break;
		}
		hdfsCloseFile(fs, in);
		hdfsDisconnect(fs);
		//------
		vector<Column*> & columns = X.col;
		for(size_t i=0; i<columns.size(); i++)
			if(columns[i]->is_dense == false)
				columns[i]->finish(row_idx);
	}

};

ColumnServer cserver;

#endif
