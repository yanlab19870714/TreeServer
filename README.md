## Distributed Training of Decision Trees and Random Forests
Following divide and conquer, this system treats each tree node to construct as a task, and parallelizes all tasks as much as possible following the idea of divide and conquer.

This allows it to use all CPU cores in a cluster to train tree models over big data. We require data to be kept on Hadoop Distributed File System for parallel loading. For more details on how to run the system, please read the file "ReadMe.txt".

### Configure the Running Environment
Please follow the documentation here to deploy: https://yanlab19870714.github.io/yanda/gthinker/deploy.html
We additionally require you to install the Boost C++ library. You may need to update the Makefiles to run on your platform.

### Contact
Da Yan: http://www.cs.uab.edu/yanda

Video Demo: https://www.youtube.com/watch?v=4DnLv_OFlIg

Email: yanda@uab.edu

### Contributors
CHOWDHURY, Md Mashiur Rahman    (Mashiur)

YAN, Da    (Daniel)
