import os
import sys
import csv

ignore_file = ["edge_fw.csv", "edge_bw.csv", "edge_count.csv"]

def get_line_list(input_file):

    line_list = []

    with open(input_file, "r") as csvfile:
        for line in csvfile:
            line_list.append(line)

    return line_list

def write_list_to_csv(result_list, output_file):

    with open(output_file, "w") as csv_file:
        csv_write = csv.writer(csv_file, delimiter = ',')
        for line in result_list:
            csv_write.writerow(line)


def get_edge_count(input_folder):

    file_list = os.listdir(input_folder)
    result_fw = []
    result_bw = []
#    print(file_list)
    for file_one in file_list:
        if file_one in ignore_file:
            continue
        input_file = os.path.join(input_folder, file_one)

        line_list = get_line_list(input_file)

        edge_fw = []
        edge_bw = []

        for line in line_list:

            if line.startswith("fw_edge_bu"):
                edge_num = int(line.split(',')[1].strip())
                edge_fw.append(edge_num)

            elif line.startswith("bw_edge_bu"):
                edge_num = int(line.split(',')[1].strip())
                edge_bw.append(edge_num)

        result_fw.append([file_one[:-4]] + edge_fw)
        result_bw.append([file_one[:-4]] + edge_bw)

    return result_fw, result_bw
#    print result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python get_edge_count.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file edge_fw/bw.csv will be stored under ", input_folder)
    result_fw, result_bw = get_edge_count(input_folder)
    write_list_to_csv(result_fw, os.path.join(input_folder, "edge_fw.csv"))
    write_list_to_csv(result_bw, os.path.join(input_folder, "edge_bw.csv"))
