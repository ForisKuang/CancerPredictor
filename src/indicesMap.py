import gzip, sys, csv
    
"""
Use case is that you pass in all 4 types of sequencing of a certain type of 
cancer into here and this section will validate whether or not you have the
correct number of arguments
"""
args = sys.argv
if (len(args) < 4):
    print("Usage ./indicesMap.py <number of files> <files> <output file>")

if (len(args) < args[1]):
    print("Usage ./indicesMap.py <number of files> <files> <output file>")

if (len(args) < args[1] + 1):
    args.append("result.csv")

resulting_map = {}

print("Running counts")
for i in range(args[1]):
    cur_file = args[i+1]
    
    lst = []
    lst.append(cur_file + "muse.tsv")
    lst.append(cur_file + "mutect.tsv")
    lst.append(cur_file + "varscan.tsv")
    lst.append(cur_file + "somaticsniper.tsv")
    for cur in lst:
        with open(cur, 'r') as f:
            count = 0
            for line in f:
                if count != 0:
                    vals = line.split('\t')
                    tuple_val = (vals[0], vals[5], vals[6])
                    if tuple_val in resulting_map:
                        resulting_map[tuple_val] += 1
                    else:
                        resulting_map[tuple_val] = 1
                count += 1
for k in resulting_map.keys():
    resulting_map[k] = resulting_map[k]/4
with open(args[len(args)-1], 'wb') as f2:
    w = csv.writer(f2)
    for k, v in resulting_map.iteritems():
        w.writerow([k[0], (k[1], k[2]), v])
print("Finished")
