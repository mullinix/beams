#!/usr/bin/python
import numpy as np
import sys

L = 0.45 # length of beam [m]
N = 3 # number of elements [1]
w = .002 # w: beam width (lateral) [m]
h = .03 # h: beam height (vertical) [m]
rho = 7850 # density of steel [kg/m^3]
m = L*w*h*rho # m: total beam mass [kg]
E = 210e9 # E: Young's modulus of steel [N/m^2]
nu = 3/10 # nu: Poisson's ratio
A = w*h # cross section area [m^2]
I = (w*h**3)/12 # moment of inertia, rectangular beam

l = L/N # equidistant node length

# element mass matrix
M = (rho*A*l/420)*np.array([[156,22*l,54,-13*l], \
							[22*l,4*l*l,13*l,-3*l*l], \
							[54,13*l,156,-22*l], \
							[-13*l,-3*l*l,-22*l,4*l*l]])

# element stiffness matrix
K = ((E*I)/(l**3))*np.array([[12,6*l,-12,6*l],\
							 [6*l,4*l*l,-6*l,2*l*l],\
							 [-12,-6*l,12,-6*l],\
							 [6*l,2*l*l,-6*l,4*l*l]])

# assembly block
# initialize global matrices
MM = np.zeros([2*(N+1),2*(N+1)])
KK = np.zeros([2*(N+1),2*(N+1)])
# loop to assemble
for i in xrange(N):
	j = 2*i
	k = j+4
	MM[j:k,j:k] = MM[j:k,j:k] + M
	KK[j:k,j:k] = KK[j:k,j:k] + K


print "Mass Matrix: "
print MM
print "Stiffness Matrix: "
print KK

eigs = np.linalg.eig(KK,MM)
print sqrt(eigs)
