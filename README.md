# iSpan: Parallel Identification of Strongly Connected Components with Spanning Trees

## Introduction
Detecting strongly connected components (SCCs) in a directed graph is crucial for understanding the structure of graphs. Most real-world graphs have one large SCC that contains the majority of the vertices, as well as many small SCCs whose sizes are reversely proportional to the frequency of their occurrences. For both types of SCCs, current approaches that rely on depth or breadth first search (DFS and BFS) face the challenges of both strict synchronization requirement and high computation cost. In this paper, we advocate a new paradigm of identifying SCCs with simple spanning trees, since SCC detection requires only the knowledge of connectivity among the vertices. We have developed a prototype called ISPAN, which consists of parallel, relaxed synchronization construction of spanning trees for detecting the large and small SCCs, combined with fast trims for small SCCs. We further scale ISPAN to distributed memory system by applying different distribution strategies to the data and task parallel jobs. The evaluations show that iSpan is able to significantly outperform current state-of-the-art DFS and BFS-based methods by average 18× and 4×, respectively


## Tutorial
You can find the source code in "src/", some useful scripts in "script/", and some test result under "result/".
More tutorials will be released soon.


## Reference
If you use iSpan in your project, please cite the following paper.
``
@inproceedings{ji2018s,
    title={iSpan: parallel identification of strongly connected components with spanning trees},
    author={Ji, Yuede and Liu, Hang and Huang, H Howie},
    booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage, and Analysis},
    pages={58},
    year={2018},
    organization={IEEE Press}
}

<!--- ## TODO
More related codes and files will be released soon.
* User guide
* Graph converter
* ...
-->
