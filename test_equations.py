import locuspoly_construct
from sage.all import GF, random_prime
from main_equations import kernel_poly_parameters, isogeny_codomain_parameters, locus_polynomial5

"*************************************************************************************"
"***************** Universal Locus Polynomials Psi_l for l = 3 and 5 *****************"
"*************************************************************************************"

print("*******************************************************************************")
print("************ The locus polynomials that can be universally used  **************")
print("*******************************************************************************\n")
exec(open("locuspoly_construct.py").read())


print("*******************************************************************************")
print("****** Working with a field whose characteristic has a bit length ~ 128  ******")
print("*******************************************************************************\n")
p128 = random_prime(2**128 - 1, lbound=2**127)
F128 = GF(p128)


"******************** Working with the prime l = 3  **********************************"

init_para_3 = kernel_poly_parameters(F128, 3)
param_3 = isogeny_codomain_parameters(init_para_3, 3)

print("# !!! The parameters used in the construction for subgroups of order 3 !!! #")
print("----------------------------------------------------------------------------\n")
print(f"The domain E:")
print(f"{param_3[0]}\n")
print(f"The kernel coefficients:")
print(f"{param_3[1]}\n")
print("The embedding field element:")
print(f"{param_3[2]}\n")
print("The subgroup identifier:")
print(f"{param_3[3]}\n")
print("The isogeny codomain:")
print(f"{param_3[4]}\n")


"*********************** Working with the prime l = 5 *******************************"

param_5 = locus_polynomial5(F128, 5)


print("# !!! The parameters used in the construction for subgroups of order 5 !!! #")
print("----------------------------------------------------------------------------\n")
print(f"The domain E:")
print(f"{param_5[0]}\n")
print(f"The kernel coefficients:")
print(f"{param_5[1]}\n")
print("The embedding field element:")
print(f"{param_5[2]}\n")
print("The subgroup identifier:")
print(f"{param_5[3]}\n")
print("The isogeny codomain:")
print(f"{param_5[4]}\n")


print("\n*******************************************************************************")
print("********** Testing if the locus polynomial satisfies the parameters  **********")
print("*******************************************************************************\n")
Psi_5ab = param_5[5]
lamda = param_5[3]
print(f"The evaluation of Psi_5 at the subgroup identifier lambda is:")
print(f"{Psi_5ab(lamda)}")
if Psi_5ab(lamda) == 0:
    print(f"As expected the evaluation at the unique identifier is equal to zero!!!")
else:
    print("The constructed locus polynomial does not satisfy elliptic curve parameters.")

