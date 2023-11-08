import re
import sys

input_filename = sys.argv[1]
output_filename = sys.argv[2]

with open(input_filename, "r") as input_file, open(output_filename, "w") as output_file:
    for line in input_file:
        # Replace decimal values immediately following ')'
        modified_line = re.sub(r'\)(\d+)', ')', line)

        # Write the modified line to the output file
        output_file.write(modified_line)
