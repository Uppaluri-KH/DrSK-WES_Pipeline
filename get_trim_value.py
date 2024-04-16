import os
import sys

def test():
    print ("exec main..")
    sys.stderr.write('execution ok\n')
    return "execution ok"


i1 = sys.argv[1]
i2 = sys.argv[2]
def get_trim_values(input1="", input2=""  ):
    input1 = open(input1, "r")
    input2 = open(input2, "r")

    #get the total sequence length and stored in total_seq_len variable but here we are using only one input file
    #because all files from same sample will have same information regarding total length sequence
    total_seq_len = 0
    for seq_len in input1:
        if seq_len.startswith(">>END_MODULE"):
            break
        if seq_len.startswith("Sequence length"):
            data = seq_len.strip().split("\t")
            data = data[1].strip().split("-")
            total_seq_len += int(data[1])-1

    #get the content quality values i.e., pass, fail or warn
    cond = ""
    cond2 = []
    for rows in input1:
        if rows.startswith(">>Per base sequence content"):                  #check for the desired part in the txt file
            row = rows.strip().split(">>") 
            row = row[1].strip().split("\t")                                #get the base sequence content information either pass, warn or fail
            if row[1] == "pass" or row[1] == "warn" or row[1] == "fail":    #if the sequence content is pass then store in one variable for later use
                cond2.append(row[1])
            break                                                           # break the program here only so work with rest of the data to get the other information 

    #get the information from the read second file
    for record in input2:
        if record.startswith(">>Per base sequence content"):                  #check for the desired part in the txt file
            record = record.strip().split(">>") 
            row = record[1].strip().split("\t")                               #get the base sequence content information either pass, warn or fail
            if row[1] == "pass" or row[1] == "warn" or row[1] == "fail":    #if the sequence content is pass then store in one variable for later use
                cond2.append(row[1])
            break 
    if "fail" in cond2 or "warn" in cond2:
        cond = "fail"
    else:
        cond = "pass"
        
    head_trim = 0
    tail_trim = 0

    if cond == "pass":                                                      #if condition is pass then no need to trimming 
        head_trim = 0
        tail_trim = 0
    elif cond == "warn" or cond == "fail":                                  #if condition is warn or fail script will get the number to trim the sequence upto few bases
        records = input1.readlines()
        records2 = input2.readlines()
        header = records[0].strip().split("\t")

        #get the head trim value if it is needed ohterwise it will print 0 as returned value
        for rows, dt in zip(records[1:10], records2[1:10]):
            if rows.startswith(">>END_MODULE") and dt.startswith(">>END_MODULE"):
                break
            row = rows.strip().split("\t")
            dt = dt.strip().split("\t")
            if float(row[1]) > 30 or float(row[1]) < 20 or float(row[2]) < 20.0 or float(row[2]) > 30 or float(row[3]) < 20.0 or float(row[3]) > 30 or float(row[4]) > 30.0 or float(row[4]) < 20 or float(dt[1]) > 30.0 or float(dt[1]) < 20 or float(dt[2]) < 20.0 or float(dt[2]) > 30 or float(dt[3]) < 20.0 or float(dt[3]) > 30 or float(dt[4]) > 30.0 or float(dt[4]) < 20:
                head_trim = max(int(row[0]),int(dt[0]))

        #get the tail trim value if it is needed otherwise it will print 0 as returned value
        for rows, dt in zip(records[33:39][::-1], records2[33:39][::-1]):
            if rows.startswith(">>END_MODULE") and dt.startswith(">>END_MODULE"):
                break
            row = rows.strip().split("\t")
            dt = dt.strip().split("\t")
            if float(row[1]) > 30 or float(row[1]) < 20 or float(row[2]) < 20.0 or float(row[2]) > 30 or float(row[3]) < 20.0 or float(row[3]) > 30 or float(row[4]) > 30.0 or float(row[4]) < 20 or float(dt[1]) > 30.0 or float(dt[1]) < 20 or float(dt[2]) < 20.0 or float(dt[2]) > 30 or float(dt[3]) < 20.0 or float(dt[3]) > 30 or float(dt[4]) > 30.0 or float(dt[4]) < 20:
                tail_trim_val_R1, tail_trim_val_R2= row[0].strip().split("-")[0], dt[0].strip().split("-")[0]
                tail_trim = min(tail_trim_val_R1,tail_trim_val_R2)
                if tail_trim == 0:
                    tail_trim = 0
                else:
                    tail_trim = total_seq_len - int(tail_trim)
    if tail_trim == 0:
        print("HEADCROP:{}".format(head_trim))
    else:
        print("HEADCROP:{} CROP:{}".format(head_trim,tail_trim))
    # sys.stderr.write('{}, {}\n'.format(head_trim, tail_trim))
    return head_trim, tail_trim

if __name__ == "__main__":
    get_trim_values(i1,i2)
