import aotensor as aotensor_mod

aotensor,Li,Lj,Lk,Lv=aotensor_mod.init_aotensor()

real_eps = 2.2204460492503131e-16

"""This module print the Ocean-Atmosphere tendencies tensor"""

for x in aotensor:
	print("aotensor["+str(x[0])+"]"+"["+str(x[1])+"]"+"["+str(x[2])+"]"+" = % .5E" % x[3])

