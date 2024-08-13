**python3 tree_variation.py -h**



usage: four_center_juli_august.py [-h] [-p] csv

This script processes a CSV file along with specific Newick tree files. 
It identifies the center branch and calculates the average UFBoot values 
for branches of the same Order in Groups 1 & 2 and Groups 3 & 4. 
The results are then used to write out tree variations.

For example, you can use the file 'four_0.01_1.fa_p4_0.05.satute.csv' located in the 'four_0.01' test folder.

Command: python3 tree_variations.py <csv_file>

positional arguments:
  csv                   Path to the CSV file

options:
  -h, --help            show this help message and exit
  -p, --print-p-value-output
                        Print p-value output

