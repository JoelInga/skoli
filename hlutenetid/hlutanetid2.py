from scipy.optimize import fsolve,least_squares,minimize,Bounds,root_scalar,brentq
import numpy as np

P_tx  = 10 #transmitting power of each sensor 10mW or 10dBm
f_op = 900e6 #operating frequency of 900 MHz
P_re_min = -90 #minimum received power for correct demodulatino is -90dBm
G_t = 5 #antenna gains 5dB
G_r = 5

sand_percent = 0.5 #50% sand particle density
clay_percent = 0.15 #clay percent
p_b = 1.5 #gramms per cubic centimeter
p_s = 2.66 #gramms per cubic centimeter

case_1 = 0.05 # case one is water content in soil 5%
case_2 = 0.2 #case two is water content in soil 20%
epsilon_s = (1.01+0.44*p_s)**2 - 0.062
beta_1 = 1.2748-0.519*sand_percent-0.152*clay_percent
beta_2 = 1.33797 - 0.603*sand_percent-0.166*clay_percent
epsilon_winf = 4.9 #high-frequency limit of epsilon_fw_1
tau_w =1/(2*np.pi)*0.58*10**(-10) #relaxation time for water at room temperature
epsilon_w0 = 80.1 # static dielectric constant for water at room temperature
epsilon_0 = 8.854*10**(-12) # permittivity constant of free space

epsilon_fw_1 = epsilon_winf + (epsilon_w0-epsilon_winf)/(1+(2*np.pi*f_op*tau_w)**2) #real part
sigma_eff = 0.0467 + 0.2204*p_b-0.4111*sand_percent+0.6614*clay_percent

epsilon_fw_2 = lambda m_v:  (2*np.pi*f_op*tau_w*(epsilon_w0-epsilon_winf))/(1+(2*np.pi*f_op*tau_w)**2)+sigma_eff*(p_s-p_b)/(2*np.pi*epsilon_0*f_op*p_s*m_v)


real_dielectric = lambda m_v: 1.15*(1+p_b/p_s*epsilon_s**0.65+m_v**beta_1*epsilon_fw_1**0.65-m_v)**(1/0.65)-0.68
imag_dielectric = lambda m_v: (m_v**beta_2*epsilon_fw_2(m_v)**0.65)**(1/0.65)
mu = 4*np.pi*10**(-7) #magnetic permeability

print(epsilon_fw_1)
print(epsilon_fw_2(case_1))

alpha = lambda m_v: 2*np.pi*f_op*np.sqrt(mu*real_dielectric(m_v)*epsilon_0/2*(np.sqrt(1+(imag_dielectric(m_v)/real_dielectric(m_v))**2)-1))
beta = lambda m_v:  2*np.pi*f_op*np.sqrt(mu*real_dielectric(m_v)*epsilon_0/2*(np.sqrt(1+(imag_dielectric(m_v)/real_dielectric(m_v))**2)+1))

print('epsilon_1 for 5% water: ',real_dielectric(case_1))
print('epsilon_2 for 5% water: ', imag_dielectric(case_1))
print('epsilon_1 for 20% water: ',real_dielectric(case_2))
print('epsilon_2 for 20% water: ', imag_dielectric(case_2))
print('alpha for 5% water: ',alpha(case_1))
print('beta for 5% water: ', beta(case_1))
print('alpha for 20% water: ',alpha(case_2))
print('beta for 20% water: ',beta(case_2))

L_p_5 = lambda d: 6.4 + 20*np.log10(d)+20*np.log10(beta(case_1))+8.69*alpha(case_1)*d

L_p_20 = lambda d: 6.4 + 20*np.log10(d)+20*np.log10(beta(case_1))+8.69*alpha(case_2)*d

P_r_5 = lambda d: P_tx+G_r+G_t-L_p_5(d)

P_r_20 = lambda d: P_tx+G_r+G_t-L_p_20(d)


def obj_5(d):
    return P_r_5(d)+90
d_5 = brentq(lambda d:obj_5(d),0.000000001,10)

def obj_20(d):
    return P_r_20(d)+90
d_20 = brentq(lambda d:obj_20(d),0.000000001,10)
print('Transmission range for 5% volumetric water content', d_5)
print('Transmission range for 20% volumetric water content', d_20)