'''
PROBABILITY DISTRIBUTION GUESSER
This program takes a list of numbers and uses its moments (mean, variance, skew, etc.)
    to guess which probability distribution the data follows.
It assumes the sample moments are approximately the population moments,
    so don't use it on small data sets.
I'm not sure how the results change in response to small changes in the data.
The input file should be just a list of floats separated by newlines.

Author: Evan E.
Last modified: 9/21/24
Moment formulas found on Wikipedia
'''

from math import sqrt, log, exp, gamma, pi

def mean(data_list, size):
    total = 0
    for i in range(size):
        total += data_list[i]
    return (total/size)

def variance(data_list, size, mean):
    total = 0
    for i in range(size):
        total += (data_list[i])**2
    var = (total/size) - (mean**2)
    return var

def skew(data_list, size, avg, var):
    altered_data = data_list.copy()
    var_factor =(sqrt(var))**3
    for i in range(size):
        altered_data[i] = (altered_data[i] - avg)**3
    return mean(altered_data, size)/var_factor

def kurtosis(data_list, size, avg, var):
    altered_data = data_list.copy()
    var_factor = var**2
    for i in range(size):
        altered_data[i] = (altered_data[i] - avg)**4
    return mean(altered_data, size)/var_factor

class data_set:
    '''This is a structure containg a set of data and its key parameters.'''
    def __init__(self, data_list, n):
        self.data = data_list
        self.size = n
        self.mean = mean(self.data, self.size)
        self.var = variance(self.data, self.size, self.mean)
        self.skew = skew(self.data, self.size, self.mean, self.var)
        self.kurt = kurtosis(self.data, self.size, self.mean, self.var)

def compare_moments(dset_obj):
    '''
    The data is given a moment score for each distribution:
    sum((sample moment - moment assuming distribution)^2)/(# of independent moments)
    A moment score closer to 0 indicates the relationship between the moments
    lines up more closely with that distribution.
    '''
    s_size = dset_obj.size
    m1 = dset_obj.mean
    m2 = dset_obj.var
    m3 = dset_obj.skew
    m4 = dset_obj.kurt
    #c_2, c_3, c_p, etc. represent the predicted m2, m3, and parameter using the current distribution
    print("The sample size is", s_size)
    print("Mean:", m1, "Variance:", m2, "Skewness:", m3, "Kurtosis:", m4) 
    #Continuous uniform: a = m1 - sqrt(3*m2), b = m1 + sqrt(3*m2)
    print("Continuous uniform:", ((m3)**2 + (m4 - 1.8)**2)/2)
    #Irwin-Hall [0, 1]; n = 2*m1
    print("Irwin-hall [0, 1]:", ((m2 - (m1/6))**2 + (m3)**2 +(m4 - 3 + 0.6/m1)**2)/3)
    #Bates [0, 1]; n = 1/(12*m2)
    print("Bates [0, 1]:", ((m1 - 0.5)**2 + (m3)**2 +(m4 - 3 + 14.4*m2)**2)/3)
    #Beta [0, 1]: alpha = (m1 - (m1**2) - m2)/(m2 + m1 - (m1**2)), beta = (alpha/m1) - alpha
    c_p = (m1 - (m1**2) - m2)/(m2 + m1 - (m1**2))
    c_n = (c_p/m1) - c_p
    if (c_p + c_n + 1 >= 0) and (c_p*c_n >= 0):
        c_3 = 2*(c_n - c_p)*sqrt(c_p + c_n + 1)/((c_p + c_n + 2)*sqrt(c_p*c_n))
        c_4 = 3 - 6/(c_p + c_n + 3) + 6*((c_p - c_n)**2)*(c_p + c_n + 1)/((c_p*c_n)*(c_p + c_n + 2)*(c_p + c_n + 3))
        print("Beta [0, 1]:", ((m3 - c_3)**2 + (m4 - c_4)**2)/2)
    #Exponential: 1/lambda = m1 = sqrt(m2)
    print("Exponential:", ((m2 - (m1**2))**2 + (m3 - 2)**2 +(m4 - 9)**2)/3)
    #Gamma: beta = m1/m2, alpha = (m1**2)/m2
    c_p = (m1**2)/m2
    print("Gamma:", ((m3 - 2/sqrt(c_p))**2 + (m4 - 3 - 6/c_p)**2)/2)
    #Inverse gamma with alpha > 4: 2 + (m1**2)/m2
    c_p = 2 + (m1**2)/m2
    c_3 = 4*sqrt(c_p - 2)/(c_p - 3)
    c_4 = 3 + 6*(5*c_p - 11)/((c_p - 3)*(c_p - 4))
    if (c_p > 4):
        print("Inverse gamma:", ((m3 - 2/sqrt(c_p))**2 + (m4 - 6/c_p)**2)/2)
    #Inverse Gaussian: m1 > 0 is a parameter, lambda = (m1**3)/m2
    c_p = m2/(m1**2)
    print("Inverse Gaussian:", ((m3 - 3*sqrt(c_p))**2 + (m4 - 3 - 15*c_p)**2)/2)
    #Type 1 Pareto with alpha > 4: alpha = 1 + sqrt((m1**2)/m2), s = m1*(1 - 1/alpha)
    c_p = 1 + sqrt((m1**2)/m2)
    if (c_p > 4):
        c_3 = 2*(1 + c_p)*sqrt(1 - 2/c_p)/(c_p - 3)
        c_4 = 3 + 6*((c_p**2) + c_p - 6 - 2/c_p)/((c_p - 3)*(c_p - 4))
        print("Type 1 Pareto:", ((m3 - c_3)**2 + (m4 - c_4)**2)/2)
    #Chi-squared: m1 = k
    if m1 >= 0:
        print("Chi-squared:", ((m2 - 2*m1)**2 + (m3 - (8/sqrt(m1)))**2 +(m4 - (12/m1))**2)/3)
    #Chi (not squared): k = m2 + (m1**2)
    c_p = m2 + (m1**2)
    c_1 = sqrt(2)*gamma(0.5 + c_p/2)/gamma(c_p/2)
    c_3 = c_1*(1 - 2*m2)/(m2*sqrt(m2))
    c_4 = 3 + 2*(1 - c_1*sqrt(m2)*c_3 - m2)/m2
    print("Chi:", ((m1 - c_1)**2 + (m3 - c_3)**2 +(m4 - c_4)**2)/3)
    #F-distribution with d2 > 8: d2 = 2*m1/(m1 - 1), (d2 - 2)/((d2 - 4)*m2/(m1**2) - 1)
    c_p = 2*m1/(m1 - 1)
    if (c_p > 8):
        c_n = (c_p - 4)*m2/(m1**2)
        c_n = (c_p - 2)/(c_n - 1)
        c_3 = (2*c_n + c_p - 2)*sqrt(8*c_p - 32)/((c_p - 6)*sqrt(c_n*(c_n + c_p - 2)))
        c_4 = 12*((5*c_p - 22) + (c_p - 4)*(c_p - 2)*(c_p - 2)/((c_n + c_p - 2)*c_n))
        c_4 = 3 + c_4/((c_p - 6)*(c_p - 8))
        print("Fisher-Snedecor F:", ((m3 - c_3)**2 +(m4 - c_4)**2)/2)
    #Birnbaum-Saunders: alpha = sqrt(2*((m1**2) - m2 + sqrt((m1**4) + 18*m2*(m1**2) - 3*(m2**2)))/(m2 - 5*(m1**2)))
        #beta = m1/((alpha**2)/2)
    try:
        c_p = 2*((m1**2) - m2 + sqrt((m1**4) + 18*m2*(m1**2) - 3*(m2**2)))/(m2 - 5*(m1**2)) #alpha^2
        c_3 = 4*(11*c_p + 6)*(sqrt(c_p))/((sqrt(5*c_p + 4))**3)
        c_4 = 3 + 6*c_p*(93*c_p + 40)/((5*c_p + 4)**2)
        print("Birnbaum-Saunders:", ((m3 - c_3)**2 +(m4 - c_4)**2)/2)
    except ValueError:
        pass #this is to deal with all of the sqrt(-1) cases
    #Log-normal: sigma^2 = ln(1 + m2/(m1**2)), mu = ln(m1) - sigma^2/2
    c_p = log(1 + m2/(m1**2))
    c_3 = (exp(c_p) + 2)*sqrt(exp(c_p) - 1)
    c_4 = 3 + 3*exp(2*c_p) + 2*exp(3*c_p) + exp(4*c_p)
    print("Log-normal:", ((m3 - c_3)**2 +(m4 - c_4)**2)/2)
    #Normal: m1 and m2 are the parameters
    print("Normal:", ((m3)**2 + (m4 - 3)**2)/2)
    #Standard t-distribution with nu > 4: nu = 4 + 6/(m4 - 3)
    c_p = 4 + 6/(m4 - 3)
    if (c_p > 4):
        print("Standard t:", ((m1)**2 +(m2 - (c_p/(c_p - 2)))**2 +(m3)**2)/3)
    #Skew-normal (Xi = 0): delta = sqrt((pi/2)*(1 + (m2/(m1**2))), omega = m1/(delta*sqrt(2/pi))
    c_p = sqrt(1 + m2/(m1**2))   #delta*sqrt(2/pi)
    if (1 - c_p**2 >= 0):
        c_3 = (2 - pi/2)*((c_p)**4)/(sqrt((1 - c_p**2)**3))
        c_4 = 3 + 2*(pi-3)*((c_p)**4)/((1 - c_p**2)**2)
        print("Skew-normal (Xi = 0):", ((m3 - c_3)**2 +(m4 - c_4)**2)/2)
    #exGaussian (mu = 0): lambda = 1/m1, sigma^2 = m2 - (m1**2)
    c_p = m2 - (m1**2)
    c_p = (m1**2)/c_p
    if (c_p >= 0):
        c_3 = 2*(sqrt(c_p)**3)/(sqrt(1 + c_p)**3)
        c_4 = 3*(1 + 2*c_p + 3*(c_p**2))/((1 + c_p)**2)
        print("exGaussian (mu = 0):", ((m3 - c_3)**2 +(m4 - c_4)**2)/2)
    #Laplace: beta = sqrt(6*m2/(pi**2)), mu = m1 - beta*0.5772157
    print("Laplace:", ((m3)**2 + (m4 - 6)**2)/2)
    #Gumbel: m1 is a parameter, b = sqrt(m2/2)
    print("Gumbel:", ((m3 - 1.139547)**2 + (m4 - 5.4)**2)/2)
    #Logistic: m1 is a parameter, s = sqrt(3*m2/(pi**2))
    print("Logistic:", ((m3)**2 + (m4 - 4.2)**2)/2)                  
    #Discrete uniform (n values starting at k): n = sqrt(12*m2 + 1), k = m1 - (n - 1)/2
    n = sqrt(12*m2 + 1)
    c_4 = 3 - 1.2*(n*n + 1)/(n*n - 1)
    print("Discrete uniform:", ((m3)**2 + (m4 - c_4)**2)/2)
    #Geometric (trials until success): 1/p = m1
    c_p = 1/(m1)
    if (1 - c_p >= 0):
        c_2 = (m1**2)*(1 - c_p)
        c_3 = (2 - c_p)/(sqrt(1 - c_p))
        c_4 = 9 + 1/m2
        print("Geometric:", ((m2 - c_2)**2 + (m3 - c_3)**2 +(m4 - c_4)**2)/3)
    #Binomial: p = 1 - m2/m1, n = m1/p
    c_p = 1 - (m2/m1)
    c_n = m1/c_p
    c_3 = (1 - 2*c_p)/sqrt(m2)
    c_4 = (1 - 6*(m2/c_n))/m2
    print("Binomial:", ((m3 - c_3)**2 + (m4 - c_4)**2)/2)
    #Poisson: lambda = m1 = m2
    print("Poisson:", ((m2 - m1)**2 + (m3 - (1/sqrt(m2)))**2 +(m4 - (1/m2))**2)/3)

def process_file(file_string):
    file = open(file_string, 'r')
    data_list = file.readlines()
    size = len(data_list)
    i = 0
    for i in range(size):
        data_list[i] = float(data_list[i])
    file.close()
    return (data_list, size)

if __name__ == '__main__':
    file_string = input("Enter file: ")
    try:
        (measurements, num_points) = process_file(file_string)
        moments = data_set(measurements, num_points)
        compare_moments(moments)
    except TypeError:
        print("Incorrect file format; use floats separated by newlines.")
