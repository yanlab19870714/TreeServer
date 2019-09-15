#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/algorithm/string.hpp>

using namespace std;

int NUM_COLS = 40;
string MISSING_VALUE = "?";
vector<bool> missing(NUM_COLS, false);

void convert_csv(int argc, char** argv) {
    if(argc < 3) {
        cout << "arg1 = input_file_path, arg2 = output_file_path requires " << endl;
        exit(-1);
    }

    string input_file = argv[1];
    string output_file = argv[2];

    ifstream fin(input_file);
    ofstream fout(output_file);
    string line, value;

    int col_idx = 0;

    while(getline(fin, line)) {
        istringstream sin(line);

        col_idx = 0;

        while(getline(sin, value, '\t')) { // (tab separated) https://labs.criteo.com/2013/12/download-terabyte-click-logs-2/
            boost::trim(value);

            if(col_idx > 0) {
                fout << ",";
            }

            if(value.empty()) {// https://labs.criteo.com/2013/12/download-terabyte-click-logs-2/
                missing[col_idx] = true;
                fout << MISSING_VALUE;
            } else {
                fout << value;
            }

            col_idx++;
        }

        fout << endl;
    }

    cout << "missing values found in column index : " << endl;

    for(int col = 0; col < missing.size(); col++) {
        if(missing[col]) {
            cout << col << endl;
        }
    }

    cout << "converted csv file written to path : " << output_file << endl;

    fin.close();
    fout.close();
}

int main(int argc, char** argv) {
    convert_csv(argc, argv);
    return 0;
}