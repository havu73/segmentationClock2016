import sys
import os

def create_fun_files (num_files):
	for i in range (num_files):
		df = open ("sres" + str(i) + ".run", "w")
		df.write("#PBS -N wildtype \n")
		df.write("#PBS -l nodes=1:ppn=10 \n")
		df.write("#PBS -l mem=2GB \n")
		df.write("#PBS -l file=100MB \n")
		df.write("#PBS -M hvu@colgate.edu \n")
		df.write("#PBS -q biomath \n")
		df.write("#PBS -j oe \n")
		df.write("#PBS -o output.pbs-job" + str(i) + "\n")
		df.write("#PBS -l walltime=480:00:00 \n")
		df.write("cd $PBS_O_WORKDIR \n")
		df.write("time mpirun -np 5\
				/home/hvu/segmentationClock2017/stochastic/sres/sres\
				-p 20\
				-r /home/hvu/segmentationClock2017/stochastic/sres/ranges/range" + str(i) + ".txt\
				-f /home/hvu/segmentationClock2017/stochatic/stochastic/stochastic\
				-t /home/hvu/segmentationClock2017/stochastic/sres/output/wt1.txt -G 0 > wildtype1.txt\n")
		df.close()

def main ():
	create_fun_files(10)
main()

