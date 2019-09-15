#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using namespace std;

int FEATURE_DIMENSION = 700;

void create_meta(int argc, char** argv) {
    if(argc < 2) {
        cout << "arg1 = output_file_path requires " << endl;
        exit(-1);
    }

    string output_file = argv[1];

    ofstream fout(output_file);

    string line;
    line = "#data_type,is_ordinal,is_dense,missing_value,default_value";
    fout << line << endl;
    fout << "int,false,true" << endl;

    for(int counter = 0; counter < FEATURE_DIMENSION; counter++) {
        fout << "double,true,true" << endl;
    }

    fout << endl;

    fout.close();
}

void convert_csv(int argc, char** argv) {
    if(argc < 3) {
        cout << "arg1 = input_file_path, arg2 = output_file_path requires " << endl;
        exit(-1);
    }
    string file_path = argv[1];
    string output_path = argv[2];

    ifstream fin(file_path);
    ofstream fout(output_path);

    string line, value, f_value1, f_value2;
    int col_idx = 0;

    while(getline(fin, line)) {
        //cout << line << endl;//#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        col_idx = 0;
        vector<double> features(FEATURE_DIMENSION, 0.0);
        istringstream sin(line);

        while(getline(sin, value, ' ')) {
            //cout << value << endl;//#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            if(col_idx == 0) {
                fout << value << ",";
                col_idx++;
                continue;
            }

            if(col_idx == 1) {
                col_idx++;
                continue;
            }

            istringstream sin2(value);

            getline(sin2, f_value1, ':');
            getline(sin2, f_value2, ':');

            //cout << "f_value1 = " << f_value1 << endl;//#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            //cout << "f_value2 = " << f_value2 << endl;//#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            features[stoi(f_value1)] = stod(f_value2);

            col_idx++;
        }

        for(int counter = 0; counter < features.size(); counter++) {
            fout << features[counter];

            if(counter != features.size() - 1) {
                fout << ",";
            }
        }

        fout << endl;
    }

    fin.close();
    fout.close();

    cout << "result written to path : " << output_path << endl;
}

int main(int argc, char** argv) {
    //convert_csv(argc, argv);
    create_meta(argc, argv);
    return 0;
}