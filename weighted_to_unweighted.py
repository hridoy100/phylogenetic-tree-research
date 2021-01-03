import re
DEBUG = False
fw = open("Datasets/un_weighted.input", "w")
with open("Datasets/weighted.input", "r") as f:
    for line in f:
        patterns = re.findall(r':\d+\.\d+', line)
        patterns += re.findall(r':\d+[eE][+-]?\d+', line)
        for pattern in patterns:
            # print("pattern: ", str(pattern))
            line = line.replace(pattern, "", 1)
        # print(line)
        fw.write(line)
fw.close()

input_file = open("Datasets/un_weighted.input", "r")
fw = open("Datasets/un_weighted_sample.input", "w")
input_rows = {}
for index, line in enumerate(input_file):
    # print("index: ", index)
    flag = 0
    patterns = re.findall(r'\)\d+\)*[\d+)*]*', line)
    # print("patters: ", patterns)
    for pattern in patterns:
        if DEBUG:
            print("pattern: ", str(pattern))
        ending_brackets = re.findall(r'\)', pattern)
        brackets = ""
        for bracket in ending_brackets:
            brackets += bracket
        line = line.replace(pattern, brackets, 1)
    if DEBUG:
        print("line", line)
    input_rows[index] = line  # strip removes \n
    fw.write(line)
fw.close()

if DEBUG:
    print(input_rows)
input_file.close()

fw = open("Datasets/un_weighted.input", "w")
with open("Datasets/un_weighted_sample.input", "r") as f:
    for line in f:
        patterns = re.findall(r':\d+', line)
        patterns += re.findall(r'[eE][+-]?\d+', line)
        print(patterns)
        for pattern in patterns:
            line = line.replace(pattern, "", 1)
        fw.write(line)
fw.close()