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

#ifndef TREE_OBJ_H
#define TREE_OBJ_H

#include "global.h"
#include "ioser.h"
#include "matrix.h"

#include <unordered_map>

using namespace std;

/////////////// Tree Node /////////////////////////////

class TreeNode {
public:
    int column_index; // column index [0......m-1] if there are m columns
    int size; // for serialize/deserialize, TODO:: may be removed
    vector<size_t>::iterator start; // start iterator of the data with which the node is formed
    vector<size_t>::iterator end; ///  end iterator of the  data with which the node is formed

    TreeNode* left = nullptr; //  left child pointer
    TreeNode* right = nullptr; //  right child pointer
};

///////////////////// Ordinal Tree Node ///////////////////////////

template <class T>
class OrdinalTreeNode : public TreeNode {
public:
    T split_value;
    string label;
};

////////////////// Categorical Tree Node /////////////////////////////

template <class T>
class CategoricalTreeNode : public TreeNode {
public:
    string label;// frequent y_label between [start, end)
    set<T> S;
    set<T> S1;
};

//////////////// Leaf Tree Node ///////////////////////////

template <class T>
class LeafNode : public TreeNode {
public:
    T label;
};

// predicting y for row_idx in test_set
// input: root (Tree_Root), test_set, y, row_idx
// output: predicted_y
template <class T>
void predict(TreeNode* root, size_t row_idx, T &predicted_y, Matrix & dataset, int MAX_TEST_DEPTH) {
    TreeNode* temp = root;

    int depth = 0;

    while (temp->column_index != -1) { // column_index == -1 means leaf node, see create_leaf function
        Column* column = dataset.get_column(temp->column_index);

        depth++;

        if(column->is_ordinal) {
            if(column->data_type == ELEM_SHORT) {
                short value;
                column->get(row_idx, &value);

                OrdinalTreeNode<short>* node = (OrdinalTreeNode<short>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_INT) {
                int value;
                column->get(row_idx, &value);

                OrdinalTreeNode<int>* node = (OrdinalTreeNode<int>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_FLOAT) {
                float value;
                column->get(row_idx, &value);

                OrdinalTreeNode<float>* node = (OrdinalTreeNode<float>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_DOUBLE) {
                double value;
                column->get(row_idx, &value);

                OrdinalTreeNode<double>* node = (OrdinalTreeNode<double>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_CHAR) {
                char value;
                column->get(row_idx, &value);

                OrdinalTreeNode<char>* node = (OrdinalTreeNode<char>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_STRING) {
                string value;
                column->get(row_idx, &value);

                OrdinalTreeNode<string>* node = (OrdinalTreeNode<string>*) temp;

                if(depth == MAX_TEST_DEPTH) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    return;
                }

                if(value >= node->split_value) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else {
                cout<< "File: " << __FILE__ << " Line : " << __LINE__ << "[Error] Unknown data type : "
                    << column->data_type << endl;
                exit(-1);
            }
        } else {
            // duplicate code in each case to avoid 49 cases
            if(column->data_type == ELEM_BOOL) {
                bool value;
                column->get(row_idx, &value);

                CategoricalTreeNode<bool>* node = (CategoricalTreeNode<bool>*) temp;

                if(node->S.find(value) == node->S.end() || (depth ==  MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_SHORT) {
                short value;
                column->get(row_idx, &value);

                CategoricalTreeNode<short>* node = (CategoricalTreeNode<short>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_INT) {
                int value;
                column->get(row_idx, &value);

                CategoricalTreeNode<int>* node = (CategoricalTreeNode<int>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_FLOAT) {
                float value;
                column->get(row_idx, &value);

                CategoricalTreeNode<float>* node = (CategoricalTreeNode<float>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_DOUBLE) {
                double value;
                column->get(row_idx, &value);

                CategoricalTreeNode<double>* node = (CategoricalTreeNode<double>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_CHAR) {
                char value;
                column->get(row_idx, &value);

                CategoricalTreeNode<char>* node = (CategoricalTreeNode<char>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else if(column->data_type == ELEM_STRING) {
                string value;
                column->get(row_idx, &value);

                CategoricalTreeNode<string>* node = (CategoricalTreeNode<string>*) temp;

                if(node->S.find(value) == node->S.end() || (depth == MAX_TEST_DEPTH)) {
                    predicted_y = boost::lexical_cast<T>(node->label);
                    if(!NO_WARNING) cout<< "File = " << __FILE__ << ", Line = " << __LINE__ << ", [Warning] value = " << value << " not found, return  max freq_label = "
                        << node->label << endl;
                    return;
                }

                if(node->S1.find(value) == node->S1.end()) {
                    temp = temp->right;
                } else {
                    temp = temp->left;
                }
            } else {
                cout<< "File: " << __FILE__ << " Line : " << __LINE__ << "[Error] Unknown data type : "
                    << column->data_type << endl;
                exit(-1);
            }

        }
    }

#ifdef ASSERT
    assert(temp->column_index == -1);
#endif

    LeafNode<T>* leafNode = (LeafNode<T>*) temp;
    predicted_y = leafNode->label;
}

template <class T>
double RMSE(vector<T> predicted_ys, Matrix & dataset, int y_index) {
    Column* y  = dataset.get_column(y_index);

    double RMSE = 0.0;
    for(size_t i = 0; i < predicted_ys.size(); i++) {
        T y_val;
        y->get(i, &y_val);
        double value = y_val - predicted_ys[i];
        RMSE += value * value;
    }

    RMSE /= predicted_ys.size();
    RMSE = sqrt(RMSE);
    return RMSE;
}

template <class T>
void free_internal_node(TreeNode* root, Column* Xi) {
    if (Xi->is_ordinal) {
        OrdinalTreeNode<T>* node = (OrdinalTreeNode<T>*) root;
        delete node;
    } else {
        CategoricalTreeNode<T>* node = (CategoricalTreeNode<T>*) root;
        delete node;
    }
}

void free_tree(TreeNode* root, int y_index, Matrix & meta) {

    if(root->column_index == -1) { //current node is a leaf
        Column* Y = meta.get_column(y_index);
        if(Y->data_type == ELEM_BOOL) {
            LeafNode<bool>* node = (LeafNode<bool>*) root;
            delete node;
        } else if(Y->data_type == ELEM_SHORT) {
            LeafNode<short>* node = (LeafNode<short>*) root;
            delete node;
        } else if(Y->data_type == ELEM_INT) {
            LeafNode<int>* node = (LeafNode<int>*) root;
            delete node;
        } else if(Y->data_type == ELEM_FLOAT) {
            LeafNode<float>* node = (LeafNode<float>*) root;
            delete node;
        } else if(Y->data_type == ELEM_DOUBLE) {
            LeafNode<double>* node = (LeafNode<double>*) root;
            delete node;
        } else if(Y->data_type == ELEM_CHAR) {
            LeafNode<char>* node = (LeafNode<char>*) root;
            delete node;
        } else if(Y->data_type == ELEM_STRING) {
            LeafNode<string>* node = (LeafNode<string> *) root;
            delete node;
        } else {
            cout<< "[ERROR] " << "y: Unknown data type " << Y->data_type << endl;
            exit(-1);
        }
    } else {
        free_tree(root->left, y_index, meta);
        free_tree(root->right, y_index, meta);

        Column* Xi = meta.get_column(root->column_index);
        if(Xi->data_type == ELEM_BOOL) {
            free_internal_node<bool>(root, Xi);
        } else if(Xi->data_type == ELEM_SHORT) {
            free_internal_node<short>(root, Xi);
        } else if(Xi->data_type == ELEM_INT) {
            free_internal_node<int>(root, Xi);
        } else if(Xi->data_type == ELEM_FLOAT) {
            free_internal_node<float>(root, Xi);
        } else if(Xi->data_type == ELEM_DOUBLE) {
            free_internal_node<double>(root, Xi);
        } else if(Xi->data_type == ELEM_CHAR) {
            free_internal_node<char>(root, Xi);
        } else if(Xi->data_type == ELEM_STRING) {
            free_internal_node<string>(root, Xi);
        } else {
            cout<< "[ERROR] " << "Xi: Unknown data type " << Xi->data_type << endl;
            exit(-1);
        }
    }
}

template <class T>
T get_majority(vector<T> predictions)
{
	unordered_map<T, size_t> histogram;
	for(int i=0; i<predictions.size(); i++)
	{
		T key = predictions[i];
		auto it = histogram.find(key);
		if(it == histogram.end()) histogram[key] = 1;
		else it->second ++;
	}
	//------
	T result;
	size_t max = 0;
	for(auto it = histogram.begin(); it!=histogram.end(); it++)
	{
		if(it->second > max)
		{
			max = it->second;
			result = it->first;
		}
	}
	return result;
}

template <class T>
double get_accuracy(Matrix & dataset, vector<TreeNode*> & rootList, int y_index, int MAX_TEST_DEPTH) {
    Column* y = dataset.get_column(y_index);
    size_t correct = 0;

    int D = y->size;

    for(size_t index = 0; index < D; index++) {
    	T ytrue;
    	y->get(index, &ytrue);
    	//------
    	vector<T> predictions(rootList.size());
    	for(int i=0; i<predictions.size(); i++)
    	{
    		T yhat; //this is needed to cope with vector<bool>
    		predict<T>(rootList[i], index, yhat, dataset, MAX_TEST_DEPTH);
    		predictions[i] = yhat;
    	}
    	//------
    	if(get_majority(predictions) == ytrue) correct++;
    }
    return correct/(double)D;
}


double get_accuracy(Matrix & dataset, vector<TreeNode*> & rootList, int y_index, int MAX_TEST_DEPTH) {
	Column* y = dataset.get_column(y_index);

	if(y->data_type == ELEM_BOOL) {
		return get_accuracy<bool>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_SHORT) {
		return get_accuracy<short>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_INT) {
		return get_accuracy<int>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_FLOAT) {
		return get_accuracy<float>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_DOUBLE) {
		return get_accuracy<double>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_CHAR) {
		return get_accuracy<char>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else if(y->data_type == ELEM_STRING) {
		return get_accuracy<string>(dataset, rootList, y_index, MAX_TEST_DEPTH);
	} else {
		cout << "Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
		exit(-1);
	}
}


// input: Tree_Root, dataset, y //todo: -------- remove
// output: predicted_y values
template <class T>
void predict(TreeNode* root, vector<T> &predicted_ys, Matrix & dataset, int MAX_TEST_DEPTH) {
    int n_sample = dataset.get_column(0)->size;
    for(size_t i = 0; i < n_sample; i++) {
        T predicted_y;
        predict<T>(root, i, predicted_y, dataset, MAX_TEST_DEPTH);
        predicted_ys.push_back(predicted_y);
    }
}

template <class T>
double calculate_RMSE(vector<TreeNode*> & rootList, Matrix & dataset, int y_index, int MAX_TEST_DEPTH) {

    Column* y = dataset.get_column(y_index);
    int D = y->size;
    int A = rootList.size();

    vector<double> avg_predictions(D, 0);

    for(int index=0; index<D; index++)
    {
    	T ytrue;
		y->get(index, &ytrue);
		//------
		double sum = 0.0;
		for(int i=0; i<A; i++)
		{
			T yhat;
			predict<T>(rootList[i], index, yhat, dataset, MAX_TEST_DEPTH);
			sum += yhat;
		}
		avg_predictions[index] = sum/A;
    }
    return RMSE(avg_predictions, dataset, y_index);
}

double calculate_RMSE(vector<TreeNode*> & rootList, Matrix & dataset, int y_index, int MAX_TEST_DEPTH) {
	Column* y = dataset.get_column(y_index);
	if(y->data_type == ELEM_SHORT) {
		return calculate_RMSE<short>(rootList, dataset, y_index, MAX_TEST_DEPTH);
	} else if (y->data_type == ELEM_INT) {
		return calculate_RMSE<int>(rootList, dataset, y_index, MAX_TEST_DEPTH);
	} else if (y->data_type == ELEM_FLOAT) {
		return calculate_RMSE<float>(rootList, dataset, y_index, MAX_TEST_DEPTH);
	} else if (y->data_type == ELEM_DOUBLE) {
		return calculate_RMSE<double>(rootList, dataset, y_index, MAX_TEST_DEPTH);
	} else {
		cout << "Wrong type, File = " << __FILE__ << ", Line = " << __LINE__ << endl;
		exit(-1);
	}
}

#endif
