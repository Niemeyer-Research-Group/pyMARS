import cantera as ct

Soln=ct.Solution('gri30.cti')

"""
for i, name in enumerate(Soln.reactions()):
	if i is 33:
		print (i, name)
"""

rxn=Soln.reaction(33)

Arr=rxn.rate
print Arr
A=Arr.pre_exponential_factor*1000
b=Arr.temperature_exponent
E=Arr.activation_energy		#J/kmol
print A, b, E

T=Soln.T
R=ct.gas_constant	#J/kmol-k
order=  b+(E/(R*T))

print(order)

print A**(-1*order)

print dir(Soln)

# k_f = A T^b \exp{-\tfrac{E}{RT}}


#needs to be 2.08e19
