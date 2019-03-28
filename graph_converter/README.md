## Graph Converter

#### Text Edge List to Binary CSR

Usage: 
For forward edges (outgoing):
g++ txt\_to\_bin\_fw\_int.cpp -o txt\_to\_bin\_fw

./txt\_to\_bin\_fw toy\_edge\_list.txt

After that, the "fw\_begin.bin" stores the begin position and "fw\_adjacent.bin" stores the adjacent lists


For backward edges (incoming):
g++ txt\_to\_bin\_bw\_int.cpp -o txt\_to\_bin\_bw

Similarly to forward edges.

Note: Be careful about the data type (int, unsigned int, long int) you are using, one can change in the source code.
