import os
import sys
import csv

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


def get_bfs_time(input_folder):

    file_list = os.listdir(input_folder)
    result_list = []
#    print(file_list)
    for file_one in file_list:
        if file_one == "bfs_time.csv":
            continue
        input_file = os.path.join(input_folder, file_one)
        time_bfs = [0.0] * 6

        frequency_list = [0] * 6
#        print time_bfs

        line_list = get_line_list(input_file)

        index = 0

        for line in line_list:

            if line.startswith("Runtime"):
                index = 0
            elif line.startswith("--->Switch"):
                index += 1
            elif line.startswith("upperbound"):
                index += 1
            elif line.startswith("Level-"):

                info = line.split()
                time_bfs[index] += float(info[4][0:-2])
                frequency_list[index] += 1

        for i in range(6):
            if frequency_list[i] != 0:
                time_bfs[i] = time_bfs[i] / frequency_list[i]

        result_list.append([file_one[:-4]] + time_bfs)

    return result_list
#    print result_list

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print "Usage: python get_bfs_time.py <input_folder>"
        exit(-1)
    input_folder = sys.argv[1]
    print("Output file bfs_time.csv will be stored under ", input_folder)
    result_list = get_bfs_time(input_folder)
    write_list_to_csv(result_list, os.path.join(input_folder, "bfs_time.csv"))
