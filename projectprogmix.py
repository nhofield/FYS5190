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

eps = 1.e-1


def vary(i_ind , j_ind):
    N = np.array([ [ 9.50242935e-01, -9.10624471e-02, 2.58345402e-01, -1.48336266e-01 ],
               [-2.77486786e-01, -6.63819194e-01, 5.26808846e-01, -4.52567786e-01 ],
               [-6.82021090e-02,  9.23897924e-02, 6.92414485e-01,  7.12302449e-01 ],
               [-1.24048052e-01,  7.36556742e-01, 4.19871735e-01, -5.15561802e-01 ] ]) #Neutralino mixing matrix


    L = np.array([ [1.77180902e-01, 9.84178301e-01],
[9.84178301e-01, -1.77180902e-01]  ]) #Left-right mixing

    N_mark = np.zeros( (4,4) , float)



    lifetimes = []
    for n in range( 100 ):
        pair = [i_ind, i_ind+1]
        if i_ind == 0:
            other_pair = [2,3]
        else:
            pass
        if i_ind == 1:
            other_pair = [0,3]
        else:
            pass
        if i_ind == 2:
            other_pair = [0,1]
        Sum = N[other_pair[0], j_ind]**2 + N[other_pair[1], j_ind]**2
        binoParam = np.linspace(eps, (np.sqrt(1. - Sum)-eps), 100)  #binofirst parameter
        #update N matrix
        N[i_ind, j_ind] = binoParam[n]


        N[i_ind+1, j_ind] = np.sqrt(1. - N[i_ind, j_ind]**2 - Sum )
        print "Checking equal 1 ",N[0,j_ind]**2 + N[1,j_ind]**2 + N[2,j_ind]**2 + N[3,j_ind]**2

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
        lifetimes.append(lifetime)
    return binoParam,lifetimes

for l in range(3):
    for n in range(4):
        binoParam, lifetimes1 = vary(l , n)
        plt.plot( binoParam , lifetimes1, label=r"$N_{%d%d}$"% (l,n) )
F = 12. #fontsize
plt.legend(fontsize=F)
plt.ylabel(r"$\tilde{\tau}$ Lifetime [ns]", fontsize=F)
plt.xlabel(r"$N_{ij}$", fontsize=F)
plt.grid(True)
plt.show()
