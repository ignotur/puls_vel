#!/usr/bin/python
import string
from math import *

file = open('prmot_pr_init.txt')

file_res = open('prmot_pr.txt', 'w')

alpha_G = radians(192.85948)
delta_G = radians(27.12825)
l_omega = radians(32.93192)

for line in file.readlines():

	# File should contain RAJD DecJD PMRA PMDec PML PMB with errors

	values = string.split(line)
	alpha = radians(float (values[4]))
	delta = radians(float (values[5]))
	mu_al = float (values[0])
	mu_dt = float (values[2])
	
	mu_al_err =  float (values[1])
	mu_dt_err =  float (values[3])

	mu_l_corr  = float(values[6])
	mu_b_corr  = float(values[8])

	mu_l_err = float(values[7])
	mu_b_err = float(values[9])


	C1 = sin(delta_G) * cos(delta)	- cos(delta_G)*sin(delta)*cos(alpha - alpha_G)
	C2 = cos(delta_G) * sin(alpha - alpha_G)
	cosb = sqrt(C1*C1 + C2*C2)

	mu_l  = ( C1*mu_al + C2 * mu_dt) / cosb
	mu_b  = (-C2*mu_al + C1 * mu_dt) / cosb


	print mu_l, mu_l_err, mu_l_err, mu_b, mu_b_err, mu_b_err
	file_res.write(str(mu_l)+"\t"+str(mu_l_err)+"\t"+str(mu_l_err)+"\t"+str(mu_b)+"\t"+str(mu_b_err)+"\t"+str(mu_b_err)+"\n")	


#	mu_l_err1  = abs(( C1*mu_al_err + C2 * mu_dt_err) / cosb)
#	mu_b_err1  = abs((-C2*mu_al_err + C1 * mu_dt_err) / cosb)
	

#	print mu_l_err, mu_l_err1, mu_b_err, mu_b_err1

#	if (abs(mu_l-mu_l_corr) > 6):
#		print mu_l, mu_l_corr, mu_b, mu_b_corr
	

#	print mu_l, mu_l_corr, mu_b, mu_b_corr

#	print values[1], mu_l, mu_b
	
#	v_l = 1e-3 * mu_l / 206265. * dist * 3e16 / 3.2e7
#	v_b = 1e-3 * mu_b / 206265. * dist * 3e16 / 3.2e7

	#v_z = v_b * sin(i)
#	v_x = v_l*cos(gb) * sin(gl)
#	v_y = v_l*cos(gb) * cos(gl)
	
#	v_z_sun = v_sun*sin(b_sun)
#	v_x_sun = v_sun*cos(b_sun)*sin(l_sun)
#	v_y_sun = v_sun*cos(b_sun)*cos(l_sun)

#	x_sun=0
#	y_sun=8.5
#	v_x_sun = -220
#	v_y_sun = 0

#	v_sun_pec = 15.6
#	l_sun = radians(48.8)
#	b_sun = radians(26.3)


#	v_sun_pec_x = v_sun_pec * cos(b_sun) * sin(l_sun)
#	v_sun_pec_y = - v_sun_pec * cos(b_sun) * cos(l_sun)
#	v_sun_pec_z = v_sun_pec * sin(b_sun)

#	theta_res = acos(y_sun*y/8.5/sqrt(x**2+y**2))

#	v_y_puls = -220*sin(theta_res)
#	v_x_puls = -220*cos(theta_res)

#	v_puls_res_x = v_x + v_x_sun + v_sun_pec_x - v_x_puls
#	v_puls_res_y = v_y + v_y_sun + v_sun_pec_y - v_y_puls
#	v_puls_res_z = v_b + v_sun_pec_z
#
#	vect_sun_pul_x = -x
#	vect_sun_pul_y = y_sun - y
#
#	norm_r = vect_sun_pul_x**2 + vect_sun_pul_y**2
#	v_orth_x = v_puls_res_x - (v_puls_res_x*vect_sun_pul_x + v_puls_res_y*vect_sun_pul_y)/norm_r*vect_sun_pul_x
#	v_orth_y = v_puls_res_y - (v_puls_res_x*vect_sun_pul_x + v_puls_res_y*vect_sun_pul_y)/norm_r*vect_sun_pul_y
#	v_orth_tot = sqrt(v_orth_x**2 + v_orth_y**2)	
#
#	v_l_res = sqrt(v_puls_res_x**2 + v_puls_res_y**2)

#	print "Before correction: ", v_l,  v_b
#	print values[1], v_l_res,  v_puls_res_z, v_l, v_b, sqrt(v_l**2+v_b**2), v_orth_tot

	
#	v = str("%.1f" % abs(v_l))
#	file_v_l.write(v)
#	file_v_l.write("\t")
#	v = str(x)
#	file_v_l.write(v)
#	file_v_l.write("\t")
#	v = str(8.5-y)
#	file_v_l.write(v)
#	file_v_l.write("\t")
#	file_v_l.write("\n");

#	v = str(z)
#	file_v_b.write(v)
#	file_v_b.write("\t")
#
#	v = str(v_b)
#	file_v_b.write(v)
#	file_v_b.write("\t")

#	if (z<0 and abs(z)>0.03):
#		z = z*(-1.)
#		v_b = v_b*(-1.);

#	print values[1], degrees(gb), dist, z, v_b
#	file_v_b.write(values[1])
#	file_v_b.write("\t")

#	v = str("%.1f" % v_b)
#	file_v_b.write(v)
#	file_v_b.write("\t")
#	v = str(z)
#	file_v_b.write(v)
#	file_v_b.write("\n");


#file_v_l.close();
#file_v_b.close();


