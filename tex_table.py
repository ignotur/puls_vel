#!/bin/python

import string
from decimal import *

f   = open('model.dat', 'r')
res = open('table.tex', 'w')

evid_energy_maxw = Decimal(f.readline()) # energy maxw
evid_flat_maxw   = Decimal(f.readline()) # flat maxw
val_maxw         = f.readline() # the actual maximum value
cred_maxw1	 = f.readline() # upper credential value
tmp=string.split (cred_maxw1)
cred_maxw1 = tmp[0]
cred_maxw2 = tmp[1]
#cred_maxw2	 = f.readline() # lower credential value

evid_energy_norm = Decimal(f.readline()) # energy norm
evid_flat_norm   = Decimal(f.readline()) # flat norm
val_norm         = f.readline() # the actual maximum value
cred_norm1	 = f.readline() # upper credential value
tmp=string.split (cred_norm1)
cred_norm1 = tmp[0]
cred_norm2 = tmp[1]
#cred_norm2	 = f.readline() # lower credential value

evid_energy_pach = Decimal(f.readline()) # energy pach
evid_flat_pach   = Decimal(f.readline()) # flat pach
val_pach         = f.readline() # the actual maximum value
cred_pach1	 = f.readline() # upper credential value
tmp=string.split (cred_pach1)
cred_pach1 = tmp[0]
cred_pach2 = tmp[1]
#cred_pach2	 = f.readline() # lower credential value

evid_energy_uni  = Decimal(f.readline()) # energy uni
evid_flat_uni    = Decimal(f.readline()) # flat uni
f.readline() # upper credential value
f.readline() # lower credential value

evid_energy_sum_maxw = Decimal(f.readline()) # energy maxw
evid_flat_sum_maxw   = Decimal(f.readline()) # flat maxw

str_complex = f.readline() ## the actual maximum value for the sum of maxwellian

arr=string.split (str_complex)

val_sum_maxw1 = round(float(arr[0]), 1)
val_sum_maxw2 = round(float(arr[1]), 1)
val_sum_maxw3 = round(float(arr[2]), 2)

str_complex = f.readline() ## credential interval for the sum of maxwellian
arr=string.split (str_complex)
cred_maxw11 = arr[0]
cred_maxw12 = arr[1]

str_complex = f.readline() ## credential interval for the sum of maxwellian
arr=string.split (str_complex)
cred_maxw21 = arr[0]
cred_maxw22 = arr[1]

str_complex = f.readline() ## credential interval for the sum of maxwellian
arr=string.split (str_complex)
cred_maxw31 = round(float(arr[0]), 2)
cred_maxw32 = round(float(arr[1]), 2)

print "& $ ", round(float(val_maxw.rstrip()), 1), "^{+", cred_maxw1.rstrip(),  "}_{-", cred_maxw2.rstrip(), "} $"  # the output of maxwellian 
print "& $ ", round(float(val_norm.rstrip()), 1), "^{+", cred_norm1.rstrip(),  "}_{-", cred_norm2.rstrip(), "} $"  # the output of normal
print "& $ ", round(float(val_pach.rstrip()), 1), "^{+", cred_pach1.rstrip(),  "}_{-", cred_pach2.rstrip(), "} $"  # the output of pachynski
print "& $ 1000 $ " #the output of uniform
print "& $ ",  val_sum_maxw1, "^{+", cred_maxw11.rstrip(),  "}_{-", cred_maxw12.rstrip(), "}, ", val_sum_maxw2, "^{+", cred_maxw21.rstrip(), " }_{-", cred_maxw22.rstrip(), "}, ", val_sum_maxw3, " ^{+", cred_maxw31, " }_{-", cred_maxw32, "} $"
print "\\\\"

evid11 = evid_flat_maxw/evid_flat_pach
if (evid11 > 0.01):
	evid11 = round(evid11, 2)
else:
	evid11 = '%.1E' % evid11
	evid11 = str(evid11).replace("E", "\\times 10^{")
	evid11 = evid11 + "}"

evid12 = evid_flat_norm/evid_flat_pach
if (evid12 > 0.01):
	evid12 = round(evid12, 2)
else:
	evid12 = '%.1E' % evid12
	evid12 = str(evid12).replace("E", "\\times 10^{")
	evid12 = evid12 + "}"

evid13 = evid_flat_pach/evid_flat_pach
if (evid13 > 0.01):
	evid13 = round(evid13, 2)
else:
	evid13 = '%.1E' % evid13
	evid13 = str(evid13).replace("E", "\\times 10^{")
	evid13 = evid13 + "}"

evid14 = evid_flat_uni/evid_flat_pach
if (evid14 > 0.01):
	evid14 = round(evid14, 2)
else:
	evid14 = '%.1E' % evid14
	evid14 = str(evid14).replace("E", "\\times 10^{")
	evid14 = evid14 + "}"

evid15 = evid_flat_sum_maxw/evid_flat_pach
if (evid15 > 0.01):
	evid15 = round(evid15, 2)
else:
	evid15 = '%.1E' % evid15
	evid15 = str(evid15).replace("E", "\\times 10^{")
	evid15 = evid15 + "}"



# $cred_sum_maxw3 $  \\\\" > table.tex
print "$\mathcal{E(\mathcal{M})_\mathrm{flat}}   $  & $ ", evid11,  " $ & $ ", evid12, " $ & $ ", evid13, " $ & $ ", evid14, " $ & $ ", evid15, " $ \\\\ "  

evid21 = evid_energy_maxw/evid_energy_pach
if (evid21 > 0.01):
	evid21 = round(evid21, 2)
else:
	evid21 = '%.1E' % evid21
	evid21 = str(evid21).replace("E", "\\times 10^{")
	evid21 = evid21 + "}"

evid22 = evid_energy_norm/evid_energy_pach
if (evid22 > 0.01):
	evid22 = round(evid22, 2)
else:
	evid22 = '%.1E' % evid22
	evid22 = str(evid22).replace("E", "\\times 10^{")
	evid22 = evid22 + "}"

evid23 = evid_energy_pach/evid_energy_pach
if (evid23 > 0.01):
	evid23 = round(evid23, 2)
else:
	evid23 = '%.1E' % evid23
	evid23 = str(evid23).replace("E", "\\times 10^{")
	evid23 = evid23 + "}"

evid24 = evid_energy_uni/evid_energy_pach
if (evid24 > 0.01):
	evid24 = round(evid24, 2)
else:
	evid24 = '%.1E' % evid24
	evid24 = str(evid24).replace("E", "\\times 10^{")
	evid24 = evid24 + "}"

evid25 = evid_energy_sum_maxw/evid_energy_pach
if (evid25 > 0.01):
	evid25 = round(evid25, 2)
else:
	evid25 = '%.1E' % evid25
	evid25 = str(evid25).replace("E", "\\times 10^{")
	evid25 = evid25 + "}"


print "$\mathcal{E(\mathcal{M})_\mathrm{energy}}   $  & $ ", evid21,  " $ & $ ", evid22, " $ & $ ", evid23, " $ & $ ", evid24, " $ & $ ", evid25, " $ \\\\ "  





#print  "$ " $evid_flat_norm   " $ & $ " $evid_flat_pach   " $ & $ " $evid_flat_uni   " $ & $ " $evid_flat_sum_maxw " $ \\\\" >> table.tex
#echo "$\mathcal{E(\mathcal{M})_\mathrm{energy}} $  & $ " $evid_energy_maxw " $ & $ " $evid_energy_norm " $ & $ " $evid_energy_pach " $ & $ " $evid_energy_uni " $ & $ " $evid_energy_sum_maxw " $ \\\\" >> table.tex 





