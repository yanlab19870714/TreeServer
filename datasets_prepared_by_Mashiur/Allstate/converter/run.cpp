#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

string MISSING_VALUE_TRAIN = "?";
string MISSING_VALUE_TEST = "##";

int num_train_columns = 35;
int num_test_columns = 34;
vector<bool> missing_train(num_train_columns, false);
vector<bool> missing_test(num_test_columns, false);

vector<int> train_column_no_missing = {3, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34};
vector<int> test_column_no_missing = {3, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33};

vector<int> drop_columns = {0,1,2,4,5,6,7};

void convert_test_csv(int argc, char** argv) {
    if(argc < 3) {
        cout << "arg1 = input_file, arg2 = output_file requires" << endl;
        exit(-1);
    }

    string input_file = argv[1];
    string output_file = argv[2];

    string line, value;
    int row_index = 0;
    int col_idx = 0;
    bool once = false;

    ifstream fin(input_file);
    ofstream fout(output_file);

    while(getline(fin, line)) {
        if(row_index == 0) {
            row_index++;
            continue;
        }

        row_index++;
        col_idx = 0;
        once = false;

        istringstream sin(line);

        while (getline(sin, value, ',')) {
            boost::trim(value);

            vector<int>::iterator it = std::find(drop_columns.begin(), drop_columns.end(), col_idx);

            if(it != drop_columns.end()) {
                col_idx++;
                continue;
            }

            if(value.empty()) {
                missing_test[col_idx] = missing_test[col_idx] || true;

                vector<int>::iterator it2 = std::find(test_column_no_missing.begin(),
                                                      test_column_no_missing.end(), col_idx);
                assert(it2 == test_column_no_missing.end());

                value = MISSING_VALUE_TEST;
            }

            if(once) {
                fout << ",";
            }

            once = true;
            fout << value;

            col_idx++;
        }

        fout << endl;
    }

    cout << "train file written to path : " << output_file << endl;
    cout << "Missing value found in following columns : " << endl;

    for(int index = 0; index < missing_test.size(); index++) {
        if(missing_test[index]) {
            cout << index << endl;
        }
    }

    fin.close();
    fout.close();
}

void convert_train_csv(int argc, char** argv) {
    if(argc < 3) {
        cout << "arg1 = input_file, arg2 = output_file requires" << endl;
        exit(-1);
    }

    string input_file = argv[1];
    string output_file = argv[2];

    string line, value;
    int row_idx = 0;
    int col_idx = 0;
    bool once = false;

    ifstream fin(input_file);
    ofstream fout(output_file);

    while(getline(fin, line)) {
        if(row_idx == 0) {
            row_idx++;
            continue;
        }

        row_idx++;

        col_idx = 0;
        once = false;
        istringstream sin(line);

        while (getline(sin, value, ',')) {
            boost::trim(value);

            vector<int>::iterator it = std::find(drop_columns.begin(), drop_columns.end(), col_idx);

            if(it != drop_columns.end()) {
                col_idx++;
                continue;
            }

            if(value.empty()) {
                missing_train[col_idx] = missing_train[col_idx] || true;

                vector<int>::iterator it2 = std::find(train_column_no_missing.begin(),
                                                      train_column_no_missing.end(), col_idx);
                assert(it2 == train_column_no_missing.end());

                value = MISSING_VALUE_TRAIN;
            }

            if(once) {
                fout << ",";
            }

            once = true;
            fout << value;

            col_idx++;
        }

        fout << endl;
    }

    cout << "train file written to path : " << output_file << endl;
    cout << "Missing value found in following columns : " << endl;

    for(int index = 0; index < missing_train.size(); index++) {
        if(missing_train[index]) {
            cout << index << endl;
        }
    }

    fin.close();
    fout.close();
}

int main(int argc, char** argv) {
    //convert_train_csv(argc, argv);
    convert_test_csv(argc, argv);
    return 0;
}