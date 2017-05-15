import sys, os
from subprocess import call
import random
def main():
	fin = open("../good.txt", "r")
	fout = open("../params_100_compare.txt", "w")
	lines = fin.readlines()
	random_sets = []
	for i in range(100):
		random_sets.append(random.randint(0, len(lines) -1))
		
	for i in range(100):
		fout.write(lines[random_sets[i]])
	fin.close()
	fout.close()
	
	#commandsto = ["mv", "params_100_compare.txt", "../stochastic/stochastic/"]
	#commanddet = ["mv", "params_100_compare.txt", "../deterministic2013Version/deterministic"]
	#if 1 == call(commandsto):
		#exit(1)
	#if 1 == call(commanddet):
		#exit(1)
main()
