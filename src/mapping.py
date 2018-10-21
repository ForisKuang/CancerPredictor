import gzip, sys, csv
    
"""
Use case is that you pass in all 4 types of sequencing of a certain type of 
cancer into here and this section will validate whether or not you have the
correct number of arguments
"""
args = sys.argv
if (len(args) < 4):
    print("Usage ./mapping <file 1 to parse> <file 2 to parse> <file 3 to parse> <file 4 to parse> <output file>")

if (len(args) < 5):
    args.append("result.csv")

resulting_map = {}

print("Running counts")
for i in range(4):
    with open(args[i + 1], 'r+') as f:
        count = 0
        for line in f:
            if count != 0:
                vals = line.split('\t')
                if vals[0] in resulting_map:
                    resulting_map[vals[0]] += 1
                else:
                    resulting_map[vals[0]] = 1
            count += 1

with open(args[5], 'wb') as f2:
    w = csv.writer(f2)
    for row in resulting_map.iteritems():
        w.writerow(row)
print("Finished")
