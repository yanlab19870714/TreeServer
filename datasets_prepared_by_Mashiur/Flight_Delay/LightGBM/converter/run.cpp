#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>

using namespace std;

int NUM_COLS = 29;
string MISSING_VALUE = "?";
vector<bool> missing_cols(NUM_COLS, false);

#define MONTH 1
#define DAY_OF_MONTH 2
#define DAY_OF_WEEK 3
#define DEP_TIME 4
#define UNIQUE_CARRIER 8
#define DEP_DELAY 15 // would be converted to y_label(dep_delayed_15_min) 0/1, classification problem
#define ORIGIN 16
#define DEST 17
#define DISTANCE 18

vector<int> column_index = {MONTH, DAY_OF_MONTH, DAY_OF_WEEK, DEP_TIME, UNIQUE_CARRIER, DEP_DELAY,
                            ORIGIN, DEST, DISTANCE};

void convert_csv(int argc, char** argv) {
    if(argc < 3) {
        cout << "arg1 = input_file_path, arg2 = output_file_path requires " << endl;
        exit(-1);
    }

    string input_path = argv[1];
    string output_path = argv[2];

    ifstream fin(input_path);
    ofstream fout(output_path);

    string line,value;
    int row_index = 0;
    int col_index = 0;


    while(getline(fin, line)) { // skipping the header information

        if(row_index == 0) {
            row_index++;
            continue;
        }

        row_index++;

        istringstream sin(line);
        col_index = 0;

        while(getline(sin, value, ',')) {
            boost::trim(value);

            auto it = std::find(column_index.begin(), column_index.end(), col_index);

            if(it != column_index.end()) {
                if(col_index > MONTH) {
                    fout << ",";
                }

                if(value.empty()) {// it seems never executed this branch
                    missing_cols[col_index] = true;

                    assert((col_index == UNIQUE_CARRIER) || (col_index == ORIGIN) || (col_index == DEST));

                    value = MISSING_VALUE;
                }

                if(col_index == DEP_DELAY) {
                    try {
                        double dep_delay = stod(value);
                        fout << (dep_delay > 15); // 1/0
                    } catch(const std::invalid_argument& ia) {
                        cout << "could not convert value : " << value << " to double " << endl;
                        fout << false;
                    }

                } else {
                    fout << value;
                }
            }

            col_index++;
        }

        fout << endl;
    }

    cout << "Missing values found in following columns : " << endl;

    for(int counter = 0; counter < missing_cols.size(); counter++) {
        if(missing_cols[counter]) {
            cout << counter << endl;
        }
    }

    cout << "csv written to path : " << output_path << endl;

    fin.close();
    fout.close();
}

int main(int argc, char** argv) {
    convert_csv(argc, argv);
    return 0;
}