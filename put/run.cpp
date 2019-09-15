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

#include "../ydhdfs.h"
#include <fstream>
#include <vector>
#include <string>

using namespace std;

void put_csv(char* localpath, char* hdfspath)
{
    if (dirCheck(hdfspath, false) == -1) return;
    //------------------------------------------
    ifstream fin(localpath);
    string line, value;
    int col_num = 0;
    if(getline(fin, line))
    {
        istringstream sin(line);
        while(getline(sin, value, ',')) col_num++;
    }
    else
    {
        cout << "[ERROR] File" << localpath << " is empty !" <<endl;
    }
    //------------------------------------------
    hdfsFS fs = getHdfsFS();
    vector<hdfsFile> handles;
    vector<string> buffer(col_num);
    for(int i=0; i<col_num; i++)
    {
        string path = hdfspath;
        path += "/column_" + to_string(i);
        hdfsFile handle = getWHandle(path.c_str(), fs);
        handles.push_back(handle);
    }
    //------------------------------------------
    do {
        istringstream sin(line);
        int col_idx = 0;
        while(getline(sin, value, ','))
        {
            value += '\n';
            buffer[col_idx].append(value);
            col_idx++;
        }
        //if buffer size > 100000, flush all files
        if (buffer[0].size() > 100000){
            for (int i = 0; i < col_num; i++)
            {
                tSize numWritten = hdfsWrite(fs, handles[i], buffer[i].c_str(), buffer[i].size());
                if (numWritten == -1) {
                    fprintf(stderr, "Failed to write file!\n");
                    exit(-1);
                }
                else if (numWritten != buffer[i].size()) {
                fprintf(stderr, "File written but content size does not match!\n");
                exit(-1);
                }
                buffer[i].clear();
            }
        }
    } while (getline(fin, line));
    //------------------------------------------
    fin.close();

    //flush rest data in buffer
    for (int i = 0; i < col_num; i++)
    {
        if(buffer[i].size() > 0)
        {
            tSize numWritten = hdfsWrite(fs, handles[i], buffer[i].c_str(), buffer[i].size());
            if (numWritten == -1) {
                fprintf(stderr, "Failed to write file!\n");
                exit(-1);
            }
            else if (numWritten != buffer[i].size()) {
            fprintf(stderr, "File written but content size does not match!\n");
            exit(-1);
            }
        }
    }

    for(int i=0; i<col_num; i++)
        hdfsCloseFile(fs, handles[i]);
    hdfsDisconnect(fs);
}

int main(int argc, char** argv)
{
    if(argc != 3)
    {
        cout<<"Put a big CSV file to HDFS as column files, one file per column"<<endl;
        return -1;
    }
    put_csv(argv[1], argv[2]);
    return 0;
}
