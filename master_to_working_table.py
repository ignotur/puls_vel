import string
from math import *
import sys
import subprocess

master = open ('master.dat', 'r')
ages   = open ('ages.dat', 'r')
work   = open ('working.dat', 'w')
names  = open ('names.dat', 'w')
tex    = open ('puls_list', 'w')
dist_list = open ('dist.txt', 'w')
dist_sch  = open ('sch_dist.txt', 'w') 
prmot_list= open ('prmot.txt', 'w') 

long_lines=[]
long_lines_ages=[]

#########################################
## Constant to transfrom proper motion ##
#########################################
alpha_G = radians(192.85948)
delta_G = radians(27.12825)
l_omega = radians(32.93192)
#########################################

#########################################
## Write a caption for the table
#########################################
tex.write("\\begin{table}\n")
tex.write("\\begin{tabular}{lccccccccc}\n")
tex.write("\\hline\n")
tex.write("JName & $\\mu_l$ & $\\mu_b$ & Ref & DM & $\\omega$ & $D_{low}$ & $D_{up}$ & P & $\dot P$  \\\\ \n")
tex.write("\\hline\n")
tex.write("\\hline\n")
########################################

lines      = master.readlines()
lines_ages = ages.readlines() 
#lines_new = string.split(lines)

for line in lines:
	lines_new = string.split(line)
	long_lines.append(lines_new)	
#	print line

for line in lines_ages:
	lines_new = string.split(line)
	long_lines_ages.append(lines_new)


long_lines.sort(key=lambda x: x[0])


for line in long_lines:
	names.write(line[0]+"\n")

num_lines = len (long_lines)
num_lines = num_lines - 1

for i in range (0, num_lines):
	name_1 = long_lines[i][0]
	name_2 = long_lines_ages[i][0]
	age = float(long_lines_ages[i][3])
	alpha  = radians(float(long_lines_ages[i][1]))
	delta  = radians(float(long_lines_ages[i][2]))

	l	      = radians(float(long_lines[i][1])) 
	b             = radians(float(long_lines[i][2]))
	mu_alpha      = float(long_lines[i][3])
	mu_alpha_lerr = float(long_lines[i][4])
	mu_alpha_rerr = float(long_lines[i][5])
	mu_delta      = float(long_lines[i][6])
	mu_delta_lerr = float(long_lines[i][7])
	mu_delta_rerr = float(long_lines[i][8])
	
#	print name_1, name_2
	if (name_1 != name_2):
		print "We have found a problem: ", name_1, " is not equal to ", name_2
		sys.exit("Catalogues problem")


	C1 = sin(delta_G) * cos(delta)	- cos(delta_G)*sin(delta)*cos(alpha - alpha_G)
	C2 = cos(delta_G) * sin(alpha - alpha_G)
	cosb = sqrt(C1*C1 + C2*C2)

#	print "Please compare -- ", cosb, "\t", cos(b)

	mu_l  = ( C1*mu_alpha + C2 * mu_delta) / cosb
	mu_b  = (-C2*mu_alpha + C1 * mu_delta) / cosb

	mu_l_lerr = abs(C1*C1*mu_alpha_lerr + C2*C2*mu_delta_lerr)/cosb/cosb;
	mu_l_rerr = abs(C1*C1*mu_alpha_rerr + C2*C2*mu_delta_rerr)/cosb/cosb;

	mu_b_lerr = abs(C2*C2*mu_alpha_lerr + C1*C1*mu_delta_lerr)/cosb/cosb;
	mu_b_rerr = abs(C2*C2*mu_alpha_rerr + C1*C1*mu_delta_rerr)/cosb/cosb;

#	print name_1,"\t",  mu_alpha, "\t",  mu_delta, "\t",  mu_l, "\t",  mu_b, "\t", mu_l_lerr, "\t", mu_b_lerr

	parral  = float(long_lines[i][9])
	parrale = abs(float(long_lines[i][10]))
	parrald = abs(float(long_lines[i][11]))


	DM = float(long_lines[i][18])	

	dlow  = float(long_lines[i][12])
	dlowe = float(long_lines[i][13])
	dup   = float(long_lines[i][14])
	dupe  = float(long_lines[i][15])

	per   = float(long_lines[i][16])
	dper  = float(long_lines[i][17])

	s14   = float(long_lines[i][19])
	
	ref   = long_lines[i][21]

	if (parral == -1.):
		p = subprocess.Popen(["./ask_conversion_NE2001.out", str(l), str(b), str(DM)], stdout=subprocess.PIPE)
		output = p.communicate()[0]
		dist_ne2001 = float(output)
		p = subprocess.Popen(["./ask_conversion_M2.out", str(l), str(b), str(DM)], stdout=subprocess.PIPE)
		output = p.communicate()[0]
		dist_m2 = float(output)
		print name_1, "\t", parral, "\t", DM, "\t", dist_ne2001, "\t", dist_m2

		dist_list.write(str(name_1)+"\t"+str(dist_ne2001)+"\t -1 \t -1\t "+ str(dlow)+"\t"+str(dlowe)+"\t"+str(dup)+"\t"+str(dupe)+"\t"+str(s14)+"\t"+str(degrees(l))+"\t"+str(degrees(b))+"\n")
		dist_sch.write(str(name_1)+"\t"+str(dist_m2)+"\t -1 \t -1\t "+ str(dlow)+"\t"+str(dlowe)+"\t"+str(dup)+"\t"+str(dupe)+"\t"+str(s14)+"\t"+str(degrees(l))+"\t"+str(degrees(b))+"\n")
	
	else:
		dist_list.write(str(name_1)+"\t"+str(parral)+"\t"+ str(parrale) + "\t" + str(parrald)+ "\t "+ str(dlow)+"\t"+str(dlowe)+"\t"+str(dup)+"\t"+str(dupe)+"\t"+str(s14)+"\t"+str(degrees(l))+"\t"+str(degrees(b))+"\n")
		dist_sch.write(str(name_1)+"\t"+str(parral)+"\t"+ str(parrale) + "\t" + str(parrald)+ "\t "+ str(dlow)+"\t"+str(dlowe)+"\t"+str(dup)+"\t"+str(dupe)+"\t"+str(s14)+"\t"+str(degrees(l))+"\t"+str(degrees(b))+"\n")
		
	prmot_list.write(str(mu_l)+"\t"+str(mu_l_lerr)+"\t"+str(mu_l_rerr)+"\t"+str(mu_b)+"\t"+str(mu_b_lerr)+"\t"+str(mu_b_rerr)+"\n")

	mu_l = round(mu_l, 2)
	mu_l_lerr = round(mu_l_lerr, 2)
	mu_l_rerr = round(mu_l_rerr, 2)

	mu_b = round(mu_b, 3)
	mu_b_lerr = round(mu_b_lerr, 2)
	mu_b_rerr = round(mu_b_rerr, 2)

	per = round(per, 5)
#	dper= round(dper, 3)


	tex.write(str(name_1)+"  &  $"+str(mu_l))
	if (mu_l_lerr == mu_l_rerr):
		tex.write("\pm "+str(mu_l_lerr)+" $  &  ")
	else:
		tex.write("^{+" + str(mu_l_lerr) + "}_{-" + str(mu_l_rerr) + "}$  & ")
	tex.write(" $" +str(mu_b))
	if (mu_b_lerr == mu_b_rerr):
		tex.write("\pm "+str(mu_b_lerr)+" $  &  ")
	else:
		tex.write("^{+" + str(mu_b_lerr) + "}_{-" + str(mu_b_rerr) + "}$  &  ")
	tex.write("\cite{" + ref + "}   & ")
	tex.write(str(DM))
	if (parral == -1):
		tex.write("  & * &  ")
	else:
		tex.write("  &  $" + str(parral) + "^{+" + str(parrale) + "}_{-" + str(parrald) + "}$  & ")

	if (dlow == -1):
		tex.write("  *  &  ")
	else:
		tex.write(" $"+str(dlow)+"\pm "+str(dlowe)+"$ & ")
	if (dup == -1):
		tex.write("  * &  ")
	else:
		tex.write(" $"+str(dup)+"\pm " + str(dupe)+"$ & ")
	tex.write(str(per)+" & " +str('%.3e' % dper) + " \\\\ \n")
	

#######################################
## Finishing our tex table
#######################################
tex.write("\\hline \n")
tex.write("\\end{tabular}\n")
tex.write("\\end{table}\n")


#for line in long_lines:
#	tex.write()

#print long_lines
