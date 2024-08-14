from sympy import symbols, cos, sqrt, diff, Function, simplify, atan2, ccode
# from sympy.printing import print_ccode


# Define the symbols
p, q, m, qa, pc, e, theta, eta = symbols('p q m qa pc e theta eta')


# Define the function g(theta) as given in the provided definition
# g_theta = Function('g')(theta)
coslode = cos(theta)
coslode2 = cos(theta)**2
num = 4 * (1 - e * e) * coslode2  + (2 * e - 1) * (2 * e - 1);
den = 2 * (1 - e * e) * coslode + (2 * e - 1) * sqrt(4 * (1 - e * e) * coslode2 + 5 * e * e - 4 * e);
# g_theta_expr = (2 * (1 - e**2) * cos(theta) + (2 * e - 1) * sqrt(4 * (1 - e**2) * cos(theta)**2 + 5 * e**2 - 4 * e)) / \
#                (4 * (1 - e**2) * cos(theta)**2 + (2 * e - 1)**2)

g_theta_expr = num / den

# Define the yield function F
# q * pow(1 + q / qa, m)  - eta * (p - pc) / g(theta, e);
# F = (q / qa)**m * g_theta - eta * (p - pc)

# F = q * (1 + q / qa)**m  - eta * (p - pc) / g_theta_expr
F = q * (1 + q / qa)**m * g_theta_expr  - eta * (p - pc) 

# Now we calculate the derivatives needed
# Derivative of F with respect to stress (Voigt notation for sigma)
sigma = symbols('sigma(0) sigma(1) sigma(2) sigma(3) sigma(4) sigma(5)')
# Mean stress p in terms of principal stresses (Voigt notation)
p_voigt = -sum(sigma[:3]) / 3
# Stress deviator q in terms of principal stresses (Voigt notation)
q_voigt = sqrt((sigma[0] - sigma[1])**2 + (sigma[1] - sigma[2])**2 + (sigma[0] - sigma[2])**2 + 6*(sigma[3]**2 + sigma[4]**2 + sigma[5]**2)) / sqrt(2)


# Define symbols for the stress components in Voigt notation
# sigma_11, sigma_22, sigma_33, sigma_12, sigma_23, sigma_13 = symbols('sigma_11 sigma_22 sigma_33 sigma_12 sigma_23 sigma_13')

# Lode angle calculation in terms of principal stresses
# Since we are using Voigt notation, we need to express the principal stresses (sigma_1, sigma_2, sigma_3) in terms of sigma_11, sigma_22, sigma_33, sigma_12, sigma_23, sigma_13.
# This requires finding the eigenvalues of the stress tensor, which is not trivial without numerical values.
# However, for the lode angle, we can use the invariants of the stress tensor which can be expressed in Voigt notation directly.

sigma_11 = sigma[0]
sigma_22 = sigma[1]
sigma_33 = sigma[2]
sigma_23 = sigma[3]
sigma_13 = sigma[4]
sigma_12 = sigma[5]


# First invariant I1 (sum of principal stresses)
I1 = sigma_11 + sigma_22 + sigma_33

# Second invariant J2 (equivalent to the von Mises stress)
J2 = ((sigma_11 - sigma_22)**2 + (sigma_22 - sigma_33)**2 + (sigma_33 - sigma_11)**2)/6 + sigma_12**2 + sigma_23**2 + sigma_13**2

# Third invariant J3
J3 = (sigma_11 + sigma_22 + sigma_33)**3/27 \
     - (sigma_11 + sigma_22 + sigma_33)/6*((sigma_11 - sigma_22)**2 + (sigma_22 - sigma_33)**2 + (sigma_33 - sigma_11)**2) \
     + sigma_11*sigma_22*sigma_33 - sigma_11*sigma_23**2 - sigma_22*sigma_13**2 - sigma_33*sigma_12**2 + 2*sigma_12*sigma_23*sigma_13

# Lode angle theta
theta_expr = simplify(atan2(sqrt(3) * (J3 * 2 * sqrt(3)), sqrt(J2**3)))

# Substitute I1, J2, and J3 with the expressions in terms of Voigt notation
theta_voigt = theta_expr.subs({
    I1: sigma_11 + sigma_22 + sigma_33,
    J2: J2,
    J3: J3
})



# Substitute p and q in terms of Voigt notation into F
F_voigt = F.subs({p: p_voigt, q: q_voigt, theta: theta_voigt})

# Derivative of F with respect to the Voigt stress components
dFd_sigma = [diff(F_voigt, s) for s in sigma]

# Derivative of F with respect to eta
dFd_eta = diff(F_voigt, eta)

# Substitute g(theta) into the derivative expressions to get the final forms
dFd_sigma = [(derivative) for derivative in dFd_sigma]
dFd_eta = (dFd_eta)

print(f"{F_voigt=}\n")
print(f"{dFd_sigma=}\n")
print(f"{dFd_eta=}\n")

print("\n\n")


fid = open("RMC_Expr.txt","w")


def printshow(s):
	print(s)
	fid.write(s+"\n")


printshow(ccode(F_voigt, "F_voigt"))
printshow("\n\n")

for i in range(6):
    printshow(ccode(dFd_sigma[i], f"dFd_sigma[{i}]"))

printshow("\n\n")
printshow(ccode(dFd_eta, "dFd_eta"))


fid.close()