import matplotlib.pyplot as plt
import numpy as np
import csv
from scipy import optimize
from scipy.special import gamma
from scipy.integrate import quad
import math
from matplotlib.lines import Line2D
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

fig, ((ax1, ax7), (ax2, ax8), (ax3, ax9), (ax4, ax10), (ax5, ax11), (ax6, ax12)) = plt.subplots(6,2)
count = []
voltage = []

with open('0507_data.csv', 'r') as csvfile:
    data_reader = csv.reader(csvfile, delimiter=',')
    element = ''
    for row in data_reader:
        if(len(row)==0):    # if the row is empty, ignore it
            pass
        else:
            voltage.append(float(row[0]))
            count.append(float(row[1]))
# observation of a plot of this data seems to suggest the following:


# V_th_1 ~ 698, and V_th_2 ~ 827 (starting values; will change for each plot by shifting one array value)
voltage_border_1 = voltage.index(698)
voltage_border_2 = voltage.index(827)



# setting initial values for variables that we will re-use for each plot

# BEGINNING BY MAKING THE EDUCATED GUESS THAT V_th_1 ~ 698, AND V_th_2 ~ 827 (will go +/- 2 array values and compare p values in range II)
# ----------------------------------------------------------------------------------------------------------------------
voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax1 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax1.legend(handles=legend_elements, loc='upper left')
ax1.plot(voltage_1, range1_regression_line)
ax1.set_title('Original Guess')
ax1.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax1.plot(voltage_2, range2_regression_line)
ax1.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax1.plot(voltage_3, range3_regression_line)
ax1.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# Now shifting the leftmost voltage border (region 1) over to the left by one array value,
# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 -= 1

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax2 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
legend_properties = {'weight':'bold'}
ax2.legend(handles=legend_elements, loc='upper left', prop=legend_properties)
ax2.set_title('Vth-left shifted 1 array unit to the left')
ax2.plot(voltage_1, range1_regression_line)
ax2.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax2.plot(voltage_2, range2_regression_line)
ax2.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax2.plot(voltage_3, range3_regression_line)
ax2.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# Now shifting the leftmost voltage border (region 1) over to the left by two array values,
# So, in this iteration, range 1 is shorter by 2 values from the original guess and range 2 is 2 values larger
# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 -= 1 # only subtracting by 1 since it's already been subtracted by 1 in the last block

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax3 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
ax3.set_title('Vth-left shifted 2 array units to the left')
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax3.legend(handles=legend_elements, loc='upper left')
ax3.plot(voltage_1, range1_regression_line)
ax3.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax3.plot(voltage_2, range2_regression_line)
ax3.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax3.plot(voltage_3, range3_regression_line)
ax3.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# Now shifting the leftmost voltage border (region 1) over to the left by 3 array values,
# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 -= 1 # total subtraction of 3 from the original on the left border

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax4 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
ax4.set_title('Vth-left shifted 3 array units to the left')
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax4.legend(handles=legend_elements, loc='upper left')
ax4.plot(voltage_1, range1_regression_line)
ax4.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax4.plot(voltage_2, range2_regression_line)
ax4.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax4.plot(voltage_3, range3_regression_line)
ax4.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')

# With this plot, we see the p-value for region 2 drop dramatically, so we know that we can stop incrementing in this region
# Now we will turn our attention toward the second threshold region, by incrementing to on the rightmost border from our original guess values

# ----------------------------------------------------------------------------------------------------------------------

# Now shifting the leftmost voltage border (region 1) over to the left by two array values,
# So, in this iteration, range 1 is shorter by 2 values from the original guess and range 2 is 2 values larger
# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 += 4 # Now test values to the right of the border, starting with an increment of 1 from the original

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax5 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax5.legend(handles=legend_elements, loc='upper left')
ax5.set_title('Vth-left shifted 1 array unit to the right')
ax5.plot(voltage_1, range1_regression_line)
ax5.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax5.plot(voltage_2, range2_regression_line)
ax5.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax5.plot(voltage_3, range3_regression_line)
ax5.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 += 1 # Total shift of 2 units to the right

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax6 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax6.legend(handles=legend_elements, loc='upper left')
ax6.set_title('Vth-left shifted 2 array units to the right')
ax6.plot(voltage_1, range1_regression_line)
ax6.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax6.plot(voltage_2, range2_regression_line)
ax6.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax6.plot(voltage_3, range3_regression_line)
ax6.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# Here it is evident that the p value for range 2 is beginning to shrink again, so we conclude that our ideal range was 1 to the left of our original guess as that is when the p-value was the highest
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
voltage_border_1 += 3  # setting the shift back to the best leftmost position (associated w/ highest p-value so far of 0.8786)
voltage_border_2 -= 1  # Total shift of 1 unit to the left

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax7 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax7.legend(handles=legend_elements, loc='upper left')
ax7.set_title('Vth-left-> 1L, Vth-right-> 1L')
ax7.plot(voltage_1, range1_regression_line)
ax7.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax7.plot(voltage_2, range2_regression_line)
ax7.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax7.plot(voltage_3, range3_regression_line)
ax7.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
voltage_border_2 -= 1 # Total shift of 2 units to the left

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax8 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax8.legend(handles=legend_elements, loc='upper left')
ax8.set_title('Vth-left-> 1L, Vth-right-> 2L')
ax8.plot(voltage_1, range1_regression_line)
ax8.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax8.plot(voltage_2, range2_regression_line)
ax8.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax8.plot(voltage_3, range3_regression_line)
ax8.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
voltage_border_2 -= 1 # Total shift of 3 units to the left

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax9 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax9.legend(handles=legend_elements, loc='upper left')
ax9.set_title('Vth-left-> 1L, Vth-right-> 3L')
ax9.plot(voltage_1, range1_regression_line)
ax9.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax9.plot(voltage_2, range2_regression_line)
ax9.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax9.plot(voltage_3, range3_regression_line)
ax9.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
voltage_border_2 += 3 # Total shift of 1 units to the right

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax10 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax10.legend(handles=legend_elements, loc='upper left')
ax10.set_title('Vth-left-> 1L, Vth-right-> 1R')
ax10.plot(voltage_1, range1_regression_line)
ax10.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax10.plot(voltage_2, range2_regression_line)
ax10.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax10.plot(voltage_3, range3_regression_line)
ax10.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
voltage_border_2 += 1 # Total shift of 2 units to the right

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax11 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax11.legend(handles=legend_elements, loc='upper left')
ax11.set_title('Vth-left-> 1L, vth-right-> 2R')
ax11.plot(voltage_1, range1_regression_line)
ax11.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax11.plot(voltage_2, range2_regression_line)
ax11.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax11.plot(voltage_3, range3_regression_line)
ax11.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

# ----------------------------------------------------------------------------------------------------------------------
voltage_border_2 += 1 # Total shift of 3 units to the right

voltage_1 = np.array(voltage[0:voltage_border_1]) # voltage corresponding to range 1
voltage_2 = voltage[voltage_border_1:voltage_border_2]  # voltage corresponding to range 2
voltage_3 = voltage[voltage_border_2:]  # voltage corresponding to range 3

count_1 = np.array(count[0:voltage_border_1])  # count corresponding to range 1
count_1_err = np.array(np.sqrt(count[0:voltage_border_1])) # taking the square root to be the approximate uncertainty (assuming poisson statistics)
count_2 = count[voltage_border_1:voltage_border_2]  # count corresponding to range 2
count_2_err = np.sqrt(count[voltage_border_1:voltage_border_2])
count_3 = count[voltage_border_2:]  # count corresponding to range 3
count_3_err = np.sqrt(count[voltage_border_2:])

print('ax12 --------------------------------------------')
range1_a, range1_b, range1_c, range1_d, range1_e, range1_f = coeff(voltage_1, count_1, count_1_err)
range1_slope, range1_slope_err, range1_intercept, range1_intercept_err = fit_2_params(range1_a, range1_b, range1_c, range1_d, range1_e, range1_f)
range1_regression_line = [range1_intercept + range1_slope*i for i in voltage_1]
sm_range1 = get_sm(voltage_1, count_1, count_1_err, range1_slope, range1_intercept)
print('V1 = {:.3f}, V2 = {:.3f}'.format(voltage_1[0], voltage_1[-1]))
print('a1, b1, ndf1 = {:.3f}, {:.3f}, {:.3f}'.format(range1_slope, range1_intercept, len(voltage_1)-2))
print('uncertainty of a (slope): {:.3f}'.format(range1_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range1_intercept_err))
print('sm = {:.4f}'.format(sm_range1))
pval1 = get_pvalue(len(voltage_1)-2, sm_range1)
print('p value for range 1 = {:.4f}'.format(pval1))


range2_a, range2_b, range2_c, range2_d, range2_e, range2_f = coeff(voltage_2, count_2, count_2_err)
range2_slope, range2_slope_err, range2_intercept, range2_intercept_err = fit_2_params(range2_a, range2_b, range2_c, range2_d, range2_e, range2_f)
range2_regression_line = [range2_intercept + range2_slope*i for i in voltage_2]
sm_range2 = get_sm(voltage_2, count_2, count_2_err, range2_slope, range2_intercept)
print('V2 = {:.3f}, V3 = {:.3f}'.format(voltage_2[0], voltage_2[-1]))
print('a2, b2, ndf2 = {:.3f}, {:.3f}, {:.3f}'.format(range2_slope, range2_intercept, len(voltage_2)-2))
print('uncertainty of a (slope): {:.3f}'.format(range2_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range2_intercept_err))
print('sm_2 = {:.4f}'.format(sm_range2))
pval2 = get_pvalue(len(voltage_2)-2, sm_range2)
print('p value for range 2 = {:.4f}'.format(pval2))


range3_a, range3_b, range3_c, range3_d, range3_e, range3_f = coeff(voltage_3, count_3, count_3_err)
range3_slope, range3_slope_err, range3_intercept, range3_intercept_err = fit_2_params(range3_a, range3_b, range3_c, range3_d, range3_e, range3_f)
range3_regression_line = [range3_intercept + range3_slope*i for i in voltage_3]
sm_range3 = get_sm(voltage_3, count_3, count_3_err, range3_slope, range3_intercept)
print('V3 = {:.3f}, V4 = {:.3f}'.format(voltage_3[0], voltage_3[-1]))
print('a3, b3, ndf3 = {:.3f}, {:.3f}, {:.3f}'.format(range3_slope, range3_intercept, len(voltage_3)-2))
print('uncertainty of a (slope): {:.3f}'.format(range3_slope_err))
print('uncertainty of b (intercept): {:.3f}'.format(range3_intercept_err))
print('sm_3 = {:.4f}'.format(sm_range3))
pval3 = get_pvalue(len(voltage_3)-2, sm_range3)
print('p-value for range 3 = {:.4f}'.format(pval3))
print('------------------------------------------------')
print()
legend_elements = [Line2D([0], [0], color='green', lw=4, label='p-value = {:.3f}'.format(pval2))]
ax12.legend(handles=legend_elements, loc='upper left')
ax12.set_title('Vth-left-> 1L, vth-right-> 3R')
ax12.plot(voltage_1, range1_regression_line)
ax12.errorbar(voltage_1, count_1, xerr=0, yerr=count_1_err, color = 'red', fmt='o')
ax12.plot(voltage_2, range2_regression_line)
ax12.errorbar(voltage_2, count_2, xerr=0, yerr=count_2_err, color = 'green', fmt='o')
ax12.plot(voltage_3, range3_regression_line)
ax12.errorbar(voltage_3, count_3, xerr=0, yerr=count_3_err, color='blue', fmt='o')
# ----------------------------------------------------------------------------------------------------------------------

print('To reiterate our approach...')
print('We began by assessing when the p-value of the second range (shown in green on the plot) was the')
print('largest. Checking values from the left border, this occurred when Vth-left was shifted 1 unit to '
      '\nthe left. Then we tested the right side with the same approach, while we held our Vth-left value'
      '\nconstant at the value that we previously determined (shifted one unit left). For each border we '
      '\nassessed points on both the left and right sides to make sure every possible scenario was'
      '\naccounted for.')
print()
print('As it turns out, our original guess was almost correct. The highest p-value occurred in the second plot,\n'
      'where the left (first) threshold voltage was shifted to the left by one plot point. It is indicated by\n'
      'being the only plot whose legend is in bold.\n')
print('The middle region set associated with the highest p-value has the following formula: \ny = ax +b ---> Count(middle region) = 1.506V + 6998 (this can be seen in the "ax2" section of the code above, among other significant parameters)')
plt.tight_layout()
plt.show()


