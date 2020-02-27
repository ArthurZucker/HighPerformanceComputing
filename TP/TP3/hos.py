from subprocess import Popen, PIPE, check_output
output2=check_output("ip addr | grep \'inet 134\'", shell=True)
command = 'nmap --open  -p 22 {} -oG - | awk '.format(output2[9:26].decode())    + '\'/Up$/{print $2}\''
#--open -sT
print(command)
output=check_output(command, shell=True)
output = output.decode()
file = "hostfile"
f = open(file,"w+")

tab = output.split('\n')
f.write("#"+str(len(tab)-1)+"\n")
for k in tab:
	if k != '':
		k = int(k[12:])
		if k<10:
			f.write("pc400"+str(k)+"\n")

		elif k<100:
			f.write("pc40"+str(k)+"\n")
		else:
			f.write("pc4"+str(k)+"\n")
