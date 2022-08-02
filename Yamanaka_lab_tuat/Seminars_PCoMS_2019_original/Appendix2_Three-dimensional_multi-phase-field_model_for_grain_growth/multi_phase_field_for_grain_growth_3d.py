# Python program for three-dimensional multi-phase-field simulation of polycrystalline grain growrh
# This source code was developed by Yamanaka research group of Tokyo University of Agriculture and Technology
# August 2019

# import libraries
from numba import jit
import numpy as np
import random

# set parameters and physical values
nx = 64
ny = nx
nz = nx
number_of_grains = 20
dx = 0.5e-6
dt = 0.1
nsteps = 2000
pi = np.pi
sigma = 1.0
delta = 7.0 * dx
eee = 1.0e+7

# calculate phase-field parameters (a, W, and M^phi)
aaa = 2.0 / pi * np.sqrt(2.0*delta*sigma)
www = 4.0 * sigma/delta
pmobi = 1.0e-7

# declare arrays for phase-field parameters (a_ij, w_ij, and M_ij)
wij = np.zeros((number_of_grains,number_of_grains))
aij = np.zeros((number_of_grains,number_of_grains))
mij = np.zeros((number_of_grains,number_of_grains))

# declair arrays for the phase-field variables, grain ID, and number of grain
phi = np.zeros((number_of_grains,nx,ny,nz))
phi_new = np.zeros((number_of_grains,nx,ny,nz))
mf = np.zeros((15,nx,ny,nz),dtype = int)
nf = np.zeros((nx,ny,nz),dtype = int)

# declair arrays for dirving force
eij = np.zeros((number_of_grains,number_of_grains))

# calculate phase-field parameters (a_ij, w_ij, and M_ij) and driving force
for i in range(0,number_of_grains):
    for j in range(0,number_of_grains):
        wij[i,j] = www
        aij[i,j] = aaa
        mij[i,j] = pmobi
        eij[i,j] = 0.0
        if i == j:
            wij[i,j] = 0.0
            aij[i,j] = 0.0
            mij[i,j] = 0.0
        if i == 1 or j == 1:
            eij[i,j] = eee
        if i < j:
            eij[i,j] = -eij[i,j]

# function for cauclating number of grains and phase-field variables
# this function is accelerated by JIT compiler
@jit
def update(phi,phi_new,mf,nf,eij):
    for k in range(nz):
        for m in range(ny):
            for l in range(nx):
                l_p = l + 1
                l_m = l - 1
                m_p = m + 1
                m_m = m - 1
                k_p = k + 1
                k_m = k - 1
                if l_p > nx-1:
                    l_p = l_p - nx
                if l_m < 0:
                    l_m = l_m + nx
                if m_p > ny-1:
                    m_p = m_p - ny
                if m_m < 0:
                    m_m = m_m + ny
                if k_p > nz-1:
                    k_p = k_p - nz
                if k_m < 0:
                    k_m = k_m + nz
                n = 0
                for i in range(number_of_grains):
                    if phi[i,l,m,k] > 0.0 or (phi[i,l,m,k] == 0.0 and phi[i,l_p,m,k] > 0.0 or phi[i,l_m,m,k] > 0.0 or phi[i,l,m_p,k] > 0.0 or phi[i,l,m_m,k] > 0.0 or phi[i,l,m,k_p] > 0.0 or phi[i,l,m,k_m] > 0.0):
                        n += 1
                        mf[n-1,l,m,k] = i
                        nf[l,m,k] = n

    for n in range(nz):
        for m in range(ny):
            for l in range(nx):
                l_p = l + 1
                l_m = l - 1
                m_p = m + 1
                m_m = m - 1
                n_p = n + 1
                n_m = n - 1
                if l_p > nx-1:
                    l_p = l_p - nx
                if l_m < 0:
                    l_m = l_m + nx
                if m_p > ny-1:
                    m_p = m_p - ny
                if m_m < 0:
                    m_m = m_m + ny
                if n_p > nz-1:
                    n_p = n_p - nz
                if n_m < 0:
                    n_m = n_m + nz
                for n1 in range(nf[l,m,n]):
                    i = mf[n1,l,m,n]
                    dpi = 0.0
                    for n2 in range(nf[l,m,n]):
                        j = mf[n2,l,m,n]
                        ppp = 0.0
                        for n3 in range(nf[l,m,n]):
                             k = mf[n3,l,m,n]
                             ppp += (wij[i,k]-wij[j,k])*phi[k,l,m,n]+0.5*(aij[i,k]**2 - aij[j,k]**2)*(phi[k,l_p,m,n]+phi[k,l_m,m,n]+phi[k,l,m_p,n]+phi[k,l,m_m,n]+phi[k,l,m,n_p]+phi[k,l,m,n_m]-6.0*phi[k,l,m,n])/dx/dx
                        pee = phi[i,l,m,n]*phi[j,l,m,n]
                        dpi = dpi - 2.0 * mij[i,j] / float(nf[l,m,n]) * (ppp - 8.0/pi*np.sqrt(pee)*eij[i,j])
                    phi_new[i,l,m,n] = phi[i,l,m,n] + dpi *dt

    phi_new = np.where(phi_new <= 0.0,0.0,phi_new)
    phi_new = np.where(phi_new >= 1.0,1.0,phi_new)

    for k in range(nz):
        for m in range(ny):
            for l in range(nx):
                a = np.sum(phi_new[:,l,m,k])
                phi[:,l,m,k] = phi_new[:,l,m,k] / a

# function for generate VTK file to visualize results using ParaView software
@jit
def output(nstep,phi,nf,mf):
    f = open('p_3d%d.vtk' % nstep,'w')
    f.write('# vtk DataFile Version 3.0 \n')
    f.write('p%d.vtk \n' % nstep)
    f.write('ASCII \n')
    f.write('DATASET STRUCTURED_POINTS \n')
    f.write('DIMENSIONS %d %d %d \n' % (nx, ny, nz))
    f.write('ORIGIN 0.0 0.0 0.0 \n')
    f.write('SPACING 1.0 1.0 1.0 \n')
    f.write('POINT_DATA %d \n' % (nx*ny*nz))
    f.write('SCALARS grain_boundary float \n')
    f.write('LOOKUP_TABLE default \n')
    for k in range(0,nz):
        for m in range(0,ny):
            for l in range(0,nx):
                f.write('%f \n' % np.sum(phi[:,l,m,k]*phi[:,l,m,k]))
    f.write('SCALARS grain_ID int \n')
    f.write('LOOKUP_TABLE default \n')
    for k in range(0,nz):
        for m in range(0,ny):
            for l in range(0,nx):
                num_grain = 0
                max_phi = 0.
                for n in range(nf[l,m,k]):
                    n1 = mf[n,l,m,k]
                    if phi[n1,l,m,k] > max_phi:
                        max_phi = phi[n1,l,m,k]
                        num_grain = n1
                f.write('%d \n' % num_grain)
    f.close()

# initial distribution of phase-field variables and others
# initial grains are distributed using random numbers
phi[1,:,:,:] = 1.
nf[:,:,:] = 1
mf[1,:,:,:] = 1
for i in range(1,number_of_grains):
    center_x = random.randint(1, nx)
    center_y = random.randint(1, ny)
    center_z = random.randint(1, nz)
    print(center_x,center_y,center_z)
    for k in range(0,nz):
        for m in range(0,ny):
            for l in range(0,nx):
                x = (l-center_x)*dx
                y = (m-center_y)*dx
                z = (k-center_z)*dx
                r = np.sqrt(x*x+y*y+z*z)
                if r < 3.*dx:
                    nf[l,m,k] = 1
                    mf[1,l,m,k] = i
                    phi[i,l,m,k] = 1.

# output initial condition
output(0,phi,nf,mf)

# main for loop
for nstep in range(1,nsteps+1):
    update(phi,phi_new,mf,nf,eij)
    if nstep % 100 == 0:
        output(nstep,phi,nf,mf)
