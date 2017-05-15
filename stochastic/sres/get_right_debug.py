
num_sets = 4
f = open('sres_debug.txt', 'r')
fout = []
for i in range(num_sets):
	current_fout = open("sres_debug" + str(i) + ".txt", 'w')
	fout.append(current_fout)
line = f.readline()
for i in range(num_sets):
	line = f.readline()
	while (line != "=====\n" and line != ''):
		(fout[i]).write(line)
		line = f.readline()
	(fout[i]).close()
#lines = f.readlines()
#fout = open('sres_fixed.txt','w')
#lines = lines[(4890521): (4890521 + 2445521)]
#for line in lines:
	#fout.write(line)
#f.close()
#fout.close()
