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

#include "csv.h"
#include "config.h"

/* //testing ColumnServer::load_meta(.), hdfs path
#include "column_server.h"

int main(int argc, char** argv){
    ColumnServer server;
    server.load_meta(argv[1]);
    server.X.print_meta();
}
//*/

/* //testing ColumnServer::load_column(.), hdfs paths
#include "column_server.h"

int main(int argc, char** argv){
    ColumnServer server;
    server.load_meta("/user/mashiur/DT/bank_meta.csv");
    server.X.print_meta();
    cout<<endl<<endl;
    server.load_column("/user/mashiur/DT/bank_train", 1);
    cout<<"=== Column 1 loaded ==="<<endl;
    server.X.print_column(1);
}
//*/
