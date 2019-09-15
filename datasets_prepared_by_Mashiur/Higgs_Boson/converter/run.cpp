#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

int num_attributes = 29;
vector<bool> missing_train(num_attributes, false);
vector<bool> missing_test(num_attributes, false);

void empty_check(int argc, char** argv) {
    if(argc < 2) {
        cout << "arg1 = input file requires" << endl;
        exit(-1);
    }

    string input_file = argv[1];

    ifstream fin(input_file);

    string line, value;
    int col_idx = 0;
    int row_idx = 0;

    while(getline(fin, line)) {
        col_idx = 0;

        istringstream sin(line);
        while(getline(sin, value, ',')) {
            boost::trim(value);

            if(value.empty()) {
                cout << "empty found row = " << row_idx << ", col_idx = " << col_idx << endl;
                exit(-1);
            }

            double numeric_val = stod(value);

            col_idx++;
        }

        row_idx++;
    }

    fin.close();

    cout << "data set found here : " << input_file << " is ok" << endl;
}

int main(int argc, char** argv) {
    empty_check(argc, argv);
    return 0;
}