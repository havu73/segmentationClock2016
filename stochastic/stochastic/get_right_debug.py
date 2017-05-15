
f = open('random.txt', 'r')
lines = f.readlines()
fout = open('random_fixed.txt','w')
lines = lines[(10056 + 5132): (10056 + 5132 + 9824 + 5132) ]
for line in lines:
	fout.write(line)
f.close()
fout.close()
