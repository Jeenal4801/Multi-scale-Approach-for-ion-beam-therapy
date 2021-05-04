import numpy as np
import matplotlib.pyplot as plt
from random import randint
import math
import scipy. stats as ss
from math import exp


dose = np.arange(0,8,0.001)

# x0 ans x1 are the positive parameters of probability function
x0 = 0.35
x1 = 0.04
sig = 0.371
#Ng = 2500
Se = 150         # Stopping Potential
lamda=0.15

# A549, AG1522, Hela, NB1RGB, A1722, V79, CHO are the seven kinds of cells considered in this experiment

# Dn is the diameter of respective cell
Dn_A549= 9.6
Dn_AG1522= 13.4
Dn_Hela=16.7
Dn_NB1RGB=14.8
Dn_A1722=16.3
Dn_V79=10.6
Dn_CHO = 12.7

# ns is the number density of complex damage sites on chromatin of the respective cell
ns_A549=(1.2)*(10**(-3))
ns_AG1522=(4.2)*(10**(-4))
ns_Hela=(2.2)*(10**(-4))
ns_NB1RGB=(3.2)*(10**(-4))
ns_A1722=(2.4)*(10**(-4))
ns_V79=(7.2)*(10**(-4))
ns_CHO=(4.26)*(10**(-4))

# r is the radius of respective cell
rA549=4.8
rAG1522=6.7
rHela=8.35
rNB1RGB=7.4
rA1722=8.15
rV79=5.3
rCHO=6.36

# Se(x) is the stopping potential of the cell x 
# Se(x)H is the stopping potential of the cell x in hypoxic conditions
# Here x = 1,2,...,7 and these numbers corresponds to the above mentioned cells
Se1= 100
Se1H=25
Se2=122
Se2H=17
Se3=70
Se3H=15
Se4=54
Se4H=13
Se5=105
Se5H=78
Se6=120
Se6H=100
Se7=150
Se7H=72

#gradNumOfLesion = Ns*sig


# z is the  average length of ions’ traversing through a nucleus 
# Function to get z for respective cells
def getZ(Dn):
    return (math.pi*Dn)/4  #formula for getting z

z_A549= getZ(Dn_A549)
z_AG1522= getZ(Dn_AG1522)
z_Hela=getZ(Dn_Hela)
z_NB1RGB=getZ(Dn_NB1RGB)
z_A1722=getZ(Dn_A1722)
z_V79=getZ(Dn_V79)
z_CHO = getZ(Dn_CHO)

# Nbp is the number of base pairs 
# Function to get Nbp for respective cells
def getNbp(Ns,An,z):
    return (160*An*z*Ns)/(3*math.pi)  # Formula to find Nbp where An is the area of cell

Nbp_A549= getNbp(ns_A549,144,z_A549)
Nbp_AG1522= getNbp(ns_AG1522,100,z_AG1522)
Nbp_Hela=getNbp(ns_Hela,219,z_Hela)
Nbp_NB1RGB=getNbp(ns_NB1RGB,172,z_NB1RGB)
Nbp_A1722=getNbp(ns_A1722,209,z_A1722)
Nbp_V79=getNbp(ns_V79,88,z_V79)
Nbp_CHO = getNbp(ns_CHO,127,z_CHO)


# Ng is the genome size
# Function to get Ng for respective cells
def getNg(Nbp):
    return Nbp/3.33  # Formula of Ng

Ng_A549= getNg(Nbp_A549)*1000
Ng_AG1522= getNg(Nbp_AG1522)*1000
Ng_Hela=getNg(Nbp_Hela)*1000
Ng_NB1RGB=getNg(Nbp_NB1RGB)*1000
Ng_A1722=getNg(Nbp_A1722)*1000
Ng_V79=getNg(Nbp_V79)*1000
Ng_CHO = getNg(Nbp_CHO)*1000

# Function to find the number of lethal lesions i.e. Yl for respective cells
def Yield  (Ng,Se):
    return (math.pi*sig*Ng*dose)/(16*Se)  # Using equation 15

YlA549 = Yield(Ng_A549,Se1)
YlAG1522 = Yield(Ng_AG1522,Se2)
YlHela = Yield(Ng_Hela,Se3)
YlNB1RGB = Yield(Ng_NB1RGB,Se4)
YlA1722 = Yield(Ng_A1722,Se5)
YlV79 = Yield(Ng_V79,Se6)
YlCHO = Yield(Ng_CHO,Se7)

# N_r is the average probability for SSBs
# N_r=randint(1,10)
N_r=np.random.uniform(0.04,0.08)

# Probability of producing lethal damage
Pl_r= (1- (ss.poisson.pmf(0, N_r) + ss.poisson.pmf(1, N_r) + ss.poisson.pmf(2, N_r)))
Pl_r=Pl_r*lamda

#actual cell survival curve
def expected_survival(n):
    
    cell_surv1=[]
    
    for i in n:
        cell_surv=math.exp(-i)
        cell_surv1.append(cell_surv)
        
    return cell_surv1


# CALCULATING THE PROBABILITY OF CELL SURVIVAL
# Pie represents the probability
Pie=[]

for i in YlCHO:
    
    if(i<(x0/x1)):
        Pie_surv=math.exp(-(1-x0)*i - x1*i*i) # Probability of cell survival from equation 5
        Pie.append(Pie_surv)
        
    else:
        Pie_surv=math.exp(-i) # Probability of cell survival when X is not introduced. (Equation 2)
        Pie.append(Pie_surv)

Pie_expected = expected_survival(YlCHO)
plt.yscale('log')
plt.plot(dose,Pie,color='orange',label=' Actual cell survival curve ')
plt.show()
plt.yscale('log')
plt.plot(dose,Pie,label=' actual Cell survival curve ')
plt.plot(dose,Pie_expected,color ='g',label ='Expected survival curve in absence of repair factor of')
plt.xlabel('dose')
plt.ylabel('Probability of cell survival ')
plt.legend()
plt.show()

# CALCULATING THE PROBABILITY OF DAMAGE
damage=[]

for i in Pie:
    
    val=1-i # Probability of damage
    damage.append(val)

plt.xlabel('dose')
plt.ylabel('Probability of damage')
plt.plot(dose,damage)
plt.show()
