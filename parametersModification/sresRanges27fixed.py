
import sys
import math
import os

num_fixed_rates = 32 	# 27 fixed rates and 5 critical fixed rates
num_rates = 37			# total number of rates
def write_new_ranges(pr, df, lower, higher):
	for i in range (num_fixed_rates):
		df.write("["+ pr[i] + "," + pr[i] + "]")
		df.write("\n")
	for i in range (num_fixed_rates, num_rates):
		df.write("[" + str(lower) + "," + str(higher) + "]")
		df.write("\n")

	

def create_sres_ranges_files(source_file, num_sets, lower, higher):
	sf = open (source_file, 'r')
	for i in range (num_sets):
		line = sf.next().strip()
		pr = line.split(',')
		dest_name = "range" + str(i) + ".txt"
		df = open (dest_name, "w")
		## Calculate the new ranges
		write_new_ranges(pr, df, lower, higher)
		df.close()
	sf.close()

def main ():
	create_sres_ranges_files('perfect.txt', 10, 0.1, 3);
main()
