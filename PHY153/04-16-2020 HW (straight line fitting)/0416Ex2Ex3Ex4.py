import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma
import math

#--------------------------------------- Ex2 ---------------------------------------------------
def coeff(x, f, sigma_f):
#…. Your code
    al, bl, cl, dl, el, fl = 0, 0, 0, 0, 0, 0
    for i in range(len(f)):
        al += x[i]/(sigma_f[i]**2)
        bl += 1/(sigma_f[i]**2)
        cl += f[i]/(sigma_f[i]**2)
        dl += (x[i]**2)/(sigma_f[i]**2)
        el += (x[i]*f[i])/(sigma_f[i]**2)
        fl += (f[i]**2)/(sigma_f[i]**2)
    return al,bl,cl,dl,el,fl

def fit(A,B,C,D,E,F):
    slope = (E*B - C*A)/(D*B - A**2)
    slope_error = math.sqrt(B/(D*B - A**2))
    intercept = (D*C - E*A)/(D*B - A**2)
    intercept_error = math.sqrt(D/(D*B - A**2))
    return slope,slope_error,intercept,intercept_error

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
#--------------------------------------- END OF Ex2 ---------------------------------------------------


#--------------------------------------- Ex3 ---------------------------------------------------
y = [3.4, 8.5, 11.7, 19.3, 22.9, 32.4]
yerr = [0.3, 0.8, 1.2, 1.9, 2.3, 3.2]
x = [0, 1, 2, 3, 4, 5]
xerr = [0,0,0,0,0,0]
A, B, C, D, E, F = coeff(x, y, yerr)
f_x = []
a, a_err, b, b_err = fit(A,B,C,D,E,F)  # a=slope, b=intercept
#print(a, 'x + ', b, sep='')
best_fit_line = [a*t + b for t in x]
#a = 5
#b = 3.4
gt_1_line = [5*t + 3 for t in x]
gt_2_line = [4*t + 3 for t in x]

k = len(y)-2

sm_1 = get_sm(x,gt_1_line,yerr,a,b)
p_val_1 = get_pvalue(k, sm_1)

sm_2 = get_sm(x,gt_2_line,yerr,a,b)
p_val_2 = get_pvalue(k, sm_2)

print('For the line 5t + 3:')
print('S_m =',sm_1)
print('p-value = {:.3f}'.format(p_val_1))
print('Probabilistically, the likelihood of yielding a sample of size ', len(y), ' whose S_m value is less than or equal to {:.3f} is {:.3f}. A probability this low suggests that the line is not a good fit'.format(sm_1,p_val_1), sep ='')
print()
print('For the line 4t + 3:')
print('S_m =',sm_2)
print('Probabilistically, the likelihood of yielding a sample of size ', len(y), ' whose S_m value is less than or equal to {:.3f} is {:.3f}. A probability this high suggests that the line is a good fit'.format(sm_2,p_val_2), sep ='')

print('Close the graph to see Example 4\'s result')
print('\n\n\n')


ax = plt.subplot()
ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o')
ax.plot(x, gt_1_line, label = '5*t + 3, p-val = {:.3f}'.format(p_val_1))
ax.plot(x, gt_2_line, label = '4*t + 3, p-val = {:.3f}'.format(p_val_2))
#ax.plot(x, best_fit_line)
ax.legend(loc='upper left')
#plt.show()
#--------------------------------------- END OF Ex3 ---------------------------------------------------



#--------------------------------------- Ex4 ---------------------------------------------------
print('EXAMPLE 4')
y = [0.92, 4.15, 9.78, 14.46, 17.26, 21.9]
yerr = [0.5, 1.0, 0.75, 1.35, 1.0, 1.5]
x = [0, 1, 2, 3, 4, 5]
xerr = [0,0,0,0,0,0]

A, B, C, D, E, F = coeff(x, y, yerr)
f_x = []
a, a_err, b, b_err = fit(A,B,C,D,E,F)  # a=slope, b=intercept
k = len(y)-2
sm = get_sm(x, y, yerr, a, b)
p_value = get_pvalue(k, sm)

print('S_m ={:.2f}'.format(sm))
print('p-value = {:.3f}'.format(p_value))
print('Probabilistically, the likelihood of yielding a sample of size ', len(y), ' whose S_m value is less than or equal to {:.3f} is {:.3f}. A probability this high suggests that the line is a good fit'.format(sm,p_value), sep ='')
best_fit_line = [a*t + b for t in x]

ax.clear()
ax.errorbar(x, y, xerr=xerr, yerr=yerr, fmt='o')
ax.plot(x,best_fit_line)
plt.show()