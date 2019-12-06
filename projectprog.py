import numpy as np
import matplotlib.pyplot as plt

beta = np.arctan(15.)                 #Susy parameter
theta_w = np.arccos(80.379/91.1876)   #Weinberg Angle
m_W = 80.379    #[GeV]
m_Z =  91.1876 #[GeV]
m_tau = 1.77686 #[GeV]
m_stau = 101.0 #[GeV]
m_neutralino = 140.0 #[GeV]
e = -1. #Electric charge
g = e/np.sin(theta_w)       #yukawa coupling

N = np.array([ [ 9.50242935e-01, -9.10624471e-02, 2.58345402e-01, -1.48336266e-01 ],
               [-2.77486786e-01, -6.63819194e-01, 5.26808846e-01, -4.52567786e-01 ],
               [-6.82021090e-02,  9.23897924e-02, 6.92414485e-01,  7.12302449e-01 ],
               [-1.24048052e-01,  7.36556742e-01, 4.19871735e-01, -5.15561802e-01 ] ]) #Neutralino mixing matrix


L = np.array([ [1.77180902e-01, 9.84178301e-01],
[9.84178301e-01, -1.77180902e-01]  ]) #Left-right mixing

N_mark = np.zeros( (4,4) , float)

for j in range(4):
    for i in range(4):
        if i == 0:
            N_mark[j][i] =   N[j][0]*np.cos(theta_w) + N[j][1]*np.sin(theta_w)
        if i == 1:
            N_mark[j][i] = - N[j][0]*np.sin(theta_w) + N[j][1]*np.cos(theta_w)
        else:
            N_mark[j][i] = N[j][i]

rho1 = ( g * m_tau )/( 2.*m_W*np.cos(beta) )
rho2 = ( g * np.sin(theta_w)**2 )/np.cos(theta_w)
rho3 = 0.5*g/np.cos(theta_w)
k_mup_mu = m_stau**2 -0.5*(m_tau**2 + m_neutralino**2)
#Sfermion coefficients
alpha = 0
a = rho1*np.sum(N[:,2])*L[0,alpha] + L[1,alpha]*( e*np.sum(N_mark[:,0]) - rho2*np.sum(N_mark[:,1]) )
b = rho1*np.sum(N[:,2])*L[1,alpha] - L[0,alpha]*( e*np.sum(N_mark[:,0]) + \
rho3*np.sum(N_mark[:,1]) - rho2*np.sum(N_mark[:,1]) )

Gamma = 2*(k_mup_mu*(a**2+b**2) + 2*m_neutralino*m_tau*(a*b))/(16.*m_stau*np.pi)
lifetime = 1./Gamma

#print lifetime, "Brute force"
#print 10./500.*m_stau, "article"
for i in range(4):
    print N[0,i]**2 + N[1, i]**2 + N[2,i]**2 + N[3,i]**2
