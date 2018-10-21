import gzip, csv, sys

# Get and verify command line arguments
args = sys.argv
if (len(args) < 2):
	print("Usage: ./parseToTSV <fileToParse> <optional output file name>")
	exit(-1)

# Add output filename if not specified
if (len(args) < 3):
	args.append(args[1] + ".tsv")

print("Opening file")

# Read file contents
with gzip.open(args[1]) as f:
	stuff = f.read()

print("Parsing file")

# Parse contents to find correct data
stuff = stuff.decode()
# stuff = stuff[stuff.find("Hugo"):]

print("Writing to tsv")

# Write tab-separated contents to output file
f2 = open(args[2], "w+")
f2.write(stuff)

print("Done!")
