import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
import math
import csv
#--------------------------------------- useful functions from previous examples ----------------------------------
def coeff(x, f, sigma_f):
    al, bl, cl, dl, el, fl = 0, 0, 0, 0, 0, 0
    for i in range(len(f)):
        al += x[i]/(sigma_f[i]**2)
        bl += 1/(sigma_f[i]**2)
        cl += f[i]/(sigma_f[i]**2)
        dl += (x[i]**2)/(sigma_f[i]**2)
        el += (x[i]*f[i])/(sigma_f[i]**2)
        fl += (f[i]**2)/(sigma_f[i]**2)
    return al,bl,cl,dl,el,fl

def fit_2_params(A,B,C,D,E,F):
    slope = (E*B - C*A)/(D*B - A**2)
    slope_error = math.sqrt(B/(D*B - A**2))
    intercept = (D*C - E*A)/(D*B - A**2)
    intercept_error = math.sqrt(D/(D*B - A**2))
    return slope,slope_error,intercept,intercept_error

def fit_intercept(C, a, A, B):
    intercept = (C-a*A)/B
    intercept_error = math.sqrt(1/B)
    return intercept, intercept_error

def get_sm(x,f,sigma_f,slope,intercept):
    sm=0
    for i in range(len(f)):
        sm += ((f[i] - slope*x[i] - intercept)**2)/(sigma_f[i]**2)
    return sm

def myChi2(x,k):
    return x**(k/2-1)*np.exp(-x/2)/(2**(k/2)*gamma(k/2))

def get_pvalue(k,sm):
    pval=quad(lambda x,k: myChi2(x,k), sm, np.inf,args=(k))
    return pval[0]

def myGauss(xl,mul,s2l):
    return np.exp(-0.5*(((xl-mul)/s2l)**2))/np.sqrt(2*np.pi*s2l)

# -------------------------------------------------------------------------------------------------------------------

# making lists to store the data values coming in from the csv file (which I saved from the page on blackboard, which I named 'data.csv')
Na_v = []
Na_K = []
Na_K_u = []

Pt_v = []
Pt_K = []
Pt_K_u = []

Ag_v = []
Ag_K = []
Ag_K_u = []

K_v = []
K_K = []
K_K_u = []

Cs_v = []
Cs_K = []
Cs_K_u = []
# --------------------------------------------------------------------------------

# ---------- Open the .csv file and write the values to the lists above -----------------
with open('0423Ex1.csv', 'r') as csvfile:
    data_reader = csv.reader(csvfile, delimiter=' ')
    element = ''
    for row in data_reader:
        if(len(row)==0):    # if the row is empty, ignore it
            pass
        else:
            if(row[0] == 'Na'):
                element = 'Na'
            elif(row[0] == 'Pt'):
                element='Pt'
            elif(row[0] == 'Ag'):
                element='Ag'
            elif(row[0] == 'K'):
                element='K'
            elif(row[0] == 'Cs'):
                element='Cs'

            else:
                if(element == 'Na'):
                    Na_v.append(float(row[0]))
                    Na_K.append(float(row[1]))
                    Na_K_u.append(float(row[2]))
                elif(element=='Pt'):
                    Pt_v.append(float(row[0]))
                    Pt_K.append(float(row[1]))
                    Pt_K_u.append(float(row[2]))
                elif(element=='Ag'):
                    Ag_v.append(float(row[0]))
                    Ag_K.append(float(row[1]))
                    Ag_K_u.append(float(row[2]))
                elif(element=='K'):
                    K_v.append(float(row[0]))
                    K_K.append(float(row[1]))
                    K_K_u.append(float(row[2]))
                elif(element=='Cs'):
                    Cs_v.append(float(row[0]))
                    Cs_K.append(float(row[1]))
                    Cs_K_u.append(float(row[2]))
# -------------------------------------------------------------------------------------------------
W_true = 2.3
mu = 0  # used for later calculations w/ the Gaussian distribution (once a z-score has been determined)
s2 = 1  # also used in the guassian distribution (once a z-score has been determined)

# -------------------------------  Na calculations ----------------------------------------------------
alNa, blNa, clNa, dlNa, elNa, flNa = coeff(Na_v, Na_K, Na_K_u)
Na_slope, Na_slope_err, Na_intercept, Na_intercept_error = fit_2_params(alNa, blNa, clNa, dlNa, elNa, flNa)
# print('slope= {:.3f}'.format(Na_slope))
Na_Kmax_theory = np.array([(Na_slope*i + Na_intercept) for i in Na_v])
Sm_Na = get_sm(Na_v, Na_K, Na_K_u,Na_slope, Na_intercept)
p_val_Na = get_pvalue(len(Na_v)-2, Sm_Na)

# print('Intercept error = {:.3f}'.format(Na_intercept_error))
fm_b_Na = (-Na_intercept - W_true)/(Na_intercept_error)
# print('z score = ', fm_b_Na)

p_val_Na_intercept_tuple = quad(lambda x,mu,s2: myGauss(x,mu,s2), -fm_b_Na, fm_b_Na, args=(mu,s2))
p_val_Na_intercept = 1- abs(p_val_Na_intercept_tuple[0])
# print('Chi sq = ', Sm_Na)
# print('p-val = {:.3f}'.format(p_val_Na))
# print(p_val_Na_intercept)
# -------------------------------  end of Na calculations ----------------------------------------------------


# -------------------------------  Pt calculations ----------------------------------------------------
alPt, blPt, clPt, dlPt, elPt, flPt = coeff(Pt_v, Pt_K, Pt_K_u)
Pt_slope, Pt_slope_err, Pt_intercept, Pt_intercept_error = fit_2_params(alPt, blPt, clPt, dlPt, elPt, flPt)
# print('slope= {:.3f}'.format(Pt_slope))
Pt_Kmax_theory = np.array([(Pt_slope*i + Pt_intercept) for i in Pt_v])
Sm_Pt = get_sm(Pt_v, Pt_K, Pt_K_u,Pt_slope, Pt_intercept)
p_val_Pt = get_pvalue(len(Pt_v)-2, Sm_Pt)

# print('Intercept error = {:.3f}'.format(Pt_intercept_error))
fm_b_Pt = (-Pt_intercept - W_true)/(Pt_intercept_error)
# print('z score = ', fm_b_Pt)

p_val_Pt_intercept_tuple = quad(lambda x,mu,s2: myGauss(x,mu,s2), -fm_b_Pt, fm_b_Pt, args=(mu,s2))
p_val_Pt_intercept = 1-abs(p_val_Pt_intercept_tuple[0])
# print('Chi sq = ', Sm_Pt)
# print('p-val = {:.3f}'.format(p_val_Pt))
# print(p_val_Pt_intercept)
# -------------------------------  end of Pt calculations ----------------------------------------------------


# -------------------------------  Ag calculations ----------------------------------------------------
alAg, blAg, clAg, dlAg, elAg, flAg = coeff(Ag_v, Ag_K, Ag_K_u)
Ag_slope, Ag_slope_err, Ag_intercept, Ag_intercept_error = fit_2_params(alAg, blAg, clAg, dlAg, elAg, flAg)
# print('slope= {:.3f}'.format(Ag_slope))
Ag_Kmax_theory = np.array([(Ag_slope*i + Ag_intercept) for i in Ag_v])
Sm_Ag = get_sm(Ag_v, Ag_K, Ag_K_u, Ag_slope, Ag_intercept)
p_val_Ag = get_pvalue(len(Ag_v)-2, Sm_Ag)

# print('Intercept error = {:.3f}'.format(Ag_intercept_error))
fm_b_Ag = (-Ag_intercept - W_true)/(Ag_intercept_error)
# print('z score = ', fm_b_Ag)

p_val_Ag_intercept_tuple = quad(lambda x,mu,s2: myGauss(x,mu,s2), -fm_b_Ag, fm_b_Ag, args=(mu,s2))
p_val_Ag_intercept = 1- abs(p_val_Ag_intercept_tuple[0])
# print('Chi sq = ', Sm_Ag)
# print('p-val = {:.3f}'.format(p_val_Ag))
# print(p_val_Ag_intercept)
#-------------------------------  end of Ag calculations ----------------------------------------------------


# -------------------------------  K calculations ----------------------------------------------------
alK, blK, clK, dlK, elK, flK = coeff(K_v, K_K, K_K_u)
K_slope, K_slope_err, K_intercept, K_intercept_error = fit_2_params(alK, blK, clK, dlK, elK, flK)
# print('slope= {:.3f}'.format(K_slope))
K_Kmax_theory = np.array([(K_slope*i + K_intercept) for i in K_v])
Sm_K = get_sm(K_v, K_K, K_K_u, K_slope, K_intercept)
p_val_K = get_pvalue(len(K_v)-2, Sm_K)

# print('Intercept error = {:.3f}'.format(K_intercept_error))
fm_b_K = (-K_intercept - W_true)/(K_intercept_error)
# print('z score = ', fm_b_K)

p_val_K_intercept_tuple = quad(lambda x,mu,s2: myGauss(x,mu,s2), -fm_b_K, fm_b_K, args=(mu,s2))
p_val_K_intercept = 1- abs(p_val_K_intercept_tuple[0])
# print('Chi sq = ', Sm_K)
# print('p-val = {:.3f}'.format(p_val_K))
# print(p_val_K_intercept)
#-------------------------------  end of K calculations ----------------------------------------------------


# -------------------------------  Cs calculations ----------------------------------------------------
alCs, blCs, clCs, dlCs, elCs, flCs = coeff(Cs_v, Cs_K, Cs_K_u)
Cs_slope, Cs_slope_err, Cs_intercept, Cs_intercept_error = fit_2_params(alCs, blCs, clCs, dlCs, elCs, flCs)
# print('slope= {:.3f}'.format(Cs_slope))
Cs_Kmax_theory = np.array([(Cs_slope*i + Cs_intercept) for i in Cs_v])
Sm_Cs = get_sm(Cs_v, Cs_K, Cs_K_u, Cs_slope, Cs_intercept)
p_val_Cs = get_pvalue(len(Cs_v)-2, Sm_Cs)

# print('Intercept error = {:.3f}'.format(Cs_intercept_error))
fm_b_Cs = (-Cs_intercept - W_true)/(Cs_intercept_error)
# print('z score = ', fm_b_Cs)

p_val_Cs_intercept_tuple = quad(lambda x,mu,s2: myGauss(x,mu,s2), -fm_b_Cs, fm_b_Cs, args=(mu,s2))
p_val_Cs_intercept = 1- abs(p_val_Cs_intercept_tuple[0])
# print('Chi sq = ', Sm_Cs)
# print('p-val = {:.3f}'.format(p_val_Cs))
# print(p_val_Cs_intercept)
#-------------------------------  end of Cs calculations ----------------------------------------------------

# fig, ((ax1, ax2),(ax3, ax4), (ax5, ax6)) = plt.subplots(nrows = 3, ncols=2)
ax1 = plt.subplot2grid(shape=(2,6), loc=(0,0), colspan=2)
ax2 = plt.subplot2grid((2,6), (0,2), colspan=2)
ax3 = plt.subplot2grid((2,6), (0,4), colspan=2)
ax4 = plt.subplot2grid((2,6), (1,1), colspan=2)
ax5 = plt.subplot2grid((2,6), (1,3), colspan=2)


h_avg = (1/5)*(Na_slope + Ag_slope + Cs_slope + K_slope + Pt_slope)
st_dev_h = math.sqrt(Na_slope**2 + Ag_slope**2 + Cs_slope**2 + K_slope**2 + Pt_slope**2)
h_true =  0.4135667696

h_Z_score = (h_avg - h_true)/st_dev_h
print('h Z-score = {:.3f}'.format(h_Z_score))
h_p_val_integral = quad(lambda x: myGauss(x, 0, 1), 0, h_Z_score)
h_p_val = 1-abs(h_p_val_integral[0])
print('h p-value = {:.4f}'.format(h_p_val))

ax1.set_title('K_max vs v [Na]')
ax1.legend()
ax1.errorbar(Na_v, Na_K, xerr=0, yerr=Na_K_u, fmt='o', label = 'Experimental values')
ax1.plot(Na_v, Na_Kmax_theory, label = 'Line of best fit')
handles_Na, labels_Na = ax1.get_legend_handles_labels()
slope_text_Na = patches.Patch(label = '\n a = {:.3f}, \nσ_a = {:.3f} \nb = -W = {:.3f}\nσ_b = {:.3f}\nSm = {:.3f}\nP-value_(a) ={:.4f}\nP-value_(b) = {:.4f}'.format(Na_slope, Na_slope_err, Na_intercept, Na_intercept_error, Sm_Na, p_val_Na, p_val_Na_intercept), color='white')
handles_Na.append(slope_text_Na)
ax1.legend(loc=2, prop={'size': 6}, handles=handles_Na)


ax2.set_title('K_max vs v [Pt]')
ax2.errorbar(Pt_v, Pt_K, xerr=0, yerr=Pt_K_u, fmt='o')
ax2.plot(Pt_v, Pt_Kmax_theory, label = 'Line of best fit')
handles_Pt, labels_Pt = ax2.get_legend_handles_labels()
slope_text_Pt = patches.Patch(label = '\n a = {:.3f}, \nσ_a = {:.3f} \nb = -W = {:.3f}\nσ_b = {:.3f}\nSm = {:.3f}\nP-value_(a) ={:.4f}\nP-value_(b) = {:.4f}'.format(Pt_slope, Pt_slope_err, Pt_intercept, Pt_intercept_error, Sm_Pt, p_val_Pt, p_val_Pt_intercept), color='white')
handles_Pt.append(slope_text_Pt)
ax2.legend(loc=2, prop={'size': 6}, handles=handles_Pt)

ax3.set_title('K_max vs v [Ag]')
ax3.errorbar(Ag_v, Ag_K, xerr=0, yerr=Ag_K_u, fmt='o')
ax3.plot(Ag_v, Ag_Kmax_theory, label = 'Line of best fit')
handles_Ag, labels_Ag = ax3.get_legend_handles_labels()
slope_text_Ag = patches.Patch(label = '\n a = {:.3f}, \nσ_a = {:.3f} \nb = -W = {:.3f}\nσ_b = {:.3f}\nSm = {:.3f}\nP-value_(a) ={:.4f}\nP-value_(b) = {:.4f}'.format(Ag_slope, Ag_slope_err, Ag_intercept, Ag_intercept_error, Sm_Ag, p_val_Ag, p_val_Ag_intercept), color='white')
handles_Ag.append(slope_text_Ag)
ax3.legend(loc=2, prop={'size': 6}, handles=handles_Ag)

ax4.set_title('K_max vs v [K]')
ax4.errorbar(K_v, K_K, xerr=0, yerr=K_K_u, fmt='o')
ax4.plot(K_v, K_Kmax_theory, label = 'Line of best fit')
handles_K, labels_K = ax4.get_legend_handles_labels()
slope_text_K = patches.Patch(label = '\n a = {:.3f}, \nσ_a = {:.3f} \nb = -W = {:.3f}\nσ_b = {:.3f}\nSm = {:.3f}\nP-value_(a) ={:.4f}\nP-value_(b) = {:.4f}'.format(K_slope, K_slope_err, K_intercept, K_intercept_error, Sm_K, p_val_K, p_val_K_intercept), color='white')
handles_K.append(slope_text_K)
ax4.legend(loc=2, prop={'size': 6}, handles=handles_K)

ax5.set_title('K_max vs v [Cs]')
ax5.errorbar(Cs_v, Cs_K, xerr=0, yerr=Cs_K_u, fmt='o')
ax5.plot(Cs_v, Cs_Kmax_theory, label = 'Line of best fit')
handles_Cs, labels_Cs = ax5.get_legend_handles_labels()
slope_text_Cs = patches.Patch(label = '\n a = {:.3f}, \nσ_a = {:.3f} \nb = -W = {:.3f}\nσ_b = {:.3f}\nSm = {:.3f}\nP-value_(a) ={:.4f}\nP-value_(b) = {:.4f}'.format(Cs_slope, Cs_slope_err, Cs_intercept, Cs_intercept_error, Sm_Cs, p_val_Cs, p_val_Cs_intercept), color='white')
handles_Cs.append(slope_text_Cs)
ax5.legend(loc=2, prop={'size': 6}, handles=handles_Cs)


plt.tight_layout()
plt.show()






























