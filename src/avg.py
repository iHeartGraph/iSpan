import os
import os.path
import sys
import csv

file_all = "file_name.txt"

def print_head(out_file):
    fp = open(out_file, "w")
    content = "Graph,Trim,Elephant SCC,Mice SCC,Total Time"
    fp.writelines(content + "\n")
    fp.close()

def deal_graph(file_name, out_file):
    read_csv = csv.reader(open(file_name))

    graph_name = (file_name.split('result_')[1]).split('.')[0]
    print graph_name

    flag = 0
    content = ""
    for line in read_csv:
        if flag == 0:
            if len(line) > 0:
                one = line[0].split()
                if len(one) > 0:
                    if one[0] == "Average":
                        flag = 1
                        content = graph_name
        elif flag < 5:
            content += "," + line[1].strip()
            flag += 1

#    print content
    fp = open(out_file, "a")
    fp.writelines(content + "\n")
    fp.close()

def deal_all(file_name, out_file):
    fp = open(file_name, "r")
    data = fp.readlines()
#    print data
    if len(data) < 1:
        print "No csv result files\n"
    else:
        print_head(out_file)
        for one in data:
            file_in = one.strip()
#            print file_in
            deal_graph(file_in, out_file)


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "\nUsage: python avg.py <input_path> <output_file>\n"
        exit()

    print len(sys.argv), sys.argv[0], sys.argv[1], sys.argv[2]
    cmd = "ls " + sys.argv[1] +  "result_*.csv > " + file_all;
    os.system(cmd)
    deal_all(file_all, sys.argv[2])

