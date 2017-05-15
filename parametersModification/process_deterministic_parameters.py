import sys
import math
import os

def process_deterministic_parameters(source_name, dest_name, kOnPD, kOnph11):
	sf = open (source_name, 'r')
	df = open (dest_name, 'w')
	for line in sf:
		pr = line.split(",")
		
		critpd = float(pr[27])
		critph11 = float(pr[28])
		
		kag1pn = kOnPD
		kag7pn = kOnPD
		
		kag1ph11 = kOnph11
		kag7ph11 = kOnph11
		kagdph11 = kOnph11
		
		kdg1pn = kag1pn * math.pow(critpd, 2.0)
		kdg7pn = kag7pn * math.pow(critpd, 2.0)
		
		kdg1ph11 = kag1ph11 * math.pow(critph11, 2.0)
		kdg7ph11 = kag7ph11 * math.pow(critph11, 2.0)
		kdgdph11 = kagdph11 * math.pow(critph11, 2.0)
		stringpar = [kag1pn, kag1ph11, kag7pn, kag7ph11, kagdph11, kdg1pn, kdg1ph11, kdg7pn, kdg7ph11, kdgdph11]
		parameters = ['{:.3f}'.format(x) for x in stringpar]
		df.write(",".join(pr[:27] + parameters))
		df.write("\n")	
	sf.close()
	df.close()

def main():
	path = os.getcwd()
	for i in range(6):
		des_name = path + '/stopar'+ str(i + 1) + '.txt'
		process_deterministic_parameters('good.txt', des_name, 0.1 * (i+1), 0.1 * (i +1))
	
main()
