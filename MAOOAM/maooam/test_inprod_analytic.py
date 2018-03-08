import numpy as np
from inprod_analytic import *
from params_maooam import natm, noc

init_inprod()

real_eps = 2.2204460492503131e-16

"""This module print the coefficients computed in the inprod_analytic module"""

for i in range(0, natm):
    for j in range(0, natm):
        if(abs(atmos.a[i, j]) >= real_eps):
            print ("a["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % atmos.a[i, j])
        if(abs(atmos.c[i, j]) >= real_eps):
            print ("c["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % atmos.c[i, j])
        for k in range(0, natm):
            if(abs(atmos.b[i, j, k]) >= real_eps):
                    print (
                        "b["+str(i+1)+"]["+str(j+1)+"]["+str(k+1)+"] =%.5E"
                        % atmos.b[i, j, k])
            if(abs(atmos.g[i, j, k]) >= real_eps):

                print (
                    "g["+str(i+1)+"]["+str(j+1)+"]["+str(k+1)+"] = % .5E"
                    % atmos.g[i, j, k])
for i in range(0, natm):
    for j in range(0, noc):
        if(abs(atmos.d[i, j]) >= real_eps):
            print ("d["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % atmos.d[i, j])
        if(abs(atmos.s[i, j]) >= real_eps):
            print ("s["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % atmos.s[i, j])

for i in range(0, noc):
    for j in range(0, noc):
        if(abs(ocean.M[i, j]) >= real_eps):
            print ("M["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % ocean.M[i, j])
        if(abs(ocean.N[i, j]) >= real_eps):
            print ("N["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E" % ocean.N[i, j])
        for k in range(0, noc):
            if(abs(ocean.O[i, j, k]) >= real_eps):
                print (
                    "O["+str(i+1)+"]["+str(j+1)+"]["+str(k+1)+"] = % .5E"
                    % ocean.O[i, j, k])
            if(abs(ocean.C[i, j, k]) >= real_eps):
                print (
                    "C["+str(i+1)+"]["+str(j+1)+"]["+str(k+1)+"] = % .5E"
                    % ocean.C[i, j, k])
    for j in range(0, natm):
        if(abs(ocean.K[i, j]) >= real_eps):
            print (
                "K["+str(i+1)+"]"+"["+str(j+1)+"] = % .5E"
                % ocean.K[i, j])
        if(abs(ocean.W[i, j]) >= real_eps):
            print (
                "W["+str(i+1)+"]" + "["+str(j+1)+"] = % .5E"
                % ocean.W[i, j])

