from sage.all import EllipticCurve, PolynomialRing, QQ, random_prime, GF, var

"""
To make the construction of the locus polynomial universal:
   - a and b should be indeterminates representing generic elliptic-curve coefficients, 
     rather than specific values.
   - θ should be an indeterminate representing a basis generator of Fq^d over Fq.
Consequently, Psi_l belongs to Fq^d[a,b,z][x].
The actual verifications on elliptic curves over finite fields are performed 
in the file main_equations.py. In this present file, the elliptic-curve coefficients 
a and b, as well as the embedded finite-field element θ (represented by z), 
are treated as indeterminates and appear in the resulting polynomials.
GF(q) is a constructor that requires a specified finite-field cardinality,
which we used in main_equations.py.
"""

"Definition of the base and embedding field parameters"
l = 5
p = random_prime(2^128 - 1, lbound = 5) #characteristic of the field
Fq = GF(p)
R = PolynomialRing(Fq, 'x')
x = R.gen()

d = (l-1) // 2
f = R.irreducible_element(d)
#print(f, f.is_irreducible()), test works
Fqd = Fq.extension(f, 'θ')
θ = Fqd.gen()


AB_θ = PolynomialRing(Fqd, names=['a', 'b', 'θ'])
a, b, z = AB_θ.gens()
R_ab = PolynomialRing(AB_θ, 'x')
#x = R.gen()

def affine_relations(a, b, θ, T, C):
    """
    The functions alpha, beta, and gamma are from the paper, Section 3.3.2.
    Q is the reduction of the degree-6 relation (intermediate polynomial in the paper)
    F_5(C) = C^6 + 20aC^4 + 160bC^3 - 80a^2C^2 - 128abC - 80b^2 to degree < 3 modulo R.
    """    

    #θ = var('θ')    
    alpha = (
        -θ**2 * T**2
        - (4*a*θ**2 + 18) * θ * T
        + (a**2 * θ**2 - 26*a + 28*b*θ) * θ**2
        - 36
    )
    beta = (
        12*θ * T**2
        + (28*a*θ**2 - 28*b*θ**3 + 36) * T
        + 8*a**2 * θ**3
        - 8*a*b*θ**4
        + 12*a*θ
        - 4*b*θ**2
    )
    gamma = 8*b*θ**2 * T + 16*a*b*θ**3 - 16*b**2 * θ**4 + 24*b*θ

    R_CT = θ*C**3 + 6*C**2 - (6*T + 2*a*θ)*C - 4*b*θ
    Q_CT = alpha*C**2 + beta*C + gamma
    
    
    return R_CT, Q_CT


def locuspoly_5():
    """
    Eliminate C from R(C, T) and Q(C, T) using the Sylvester resultant in C 
    The returned polynomial is a cofactor in QQ[a, b, θ][T].
    In this universal computation, we have integers and a, b, and θ in the locus.
    """
    #θ = var('θ') 
    Coeff = PolynomialRing(QQ, names=("a", "b", "θ"))
    a, b, θ = Coeff.gens()
    RT = PolynomialRing(Coeff, "T")
    T = RT.gen()
    RC = PolynomialRing(RT, "C")
    C = RC.gen()

    R_CT, Q_CT = affine_relations(a, b, θ, T, C)
    resultant = R_CT.resultant(Q_CT)#The result can be set in Fqd

    content = Coeff(16 * b**2 * θ**8)
    cofactor = []
    for coefficient in resultant.list():
        quotient, remainder = Coeff(coefficient).quo_rem(content)
        if remainder != 0:
            raise ValueError("16 b^2 θ^8 does not divide Res_C(R, Q)")
        cofactor.append(quotient)

    tilde_Psi_5 = RT(cofactor)  
  

    return tilde_Psi_5



def locuspoly_3():
    """
        For l = 3, the locus is just a normalization of the 3-division polynomial.
        In this case the θ in our paper is equal to 1.
    """

    Coeff = PolynomialRing(QQ, names=("a", "b"))
    a, b = Coeff.gens()
    RT = PolynomialRing(Coeff, "T")
    T = RT.gen()
    
    psi_3 = 3*T**4 + 6*a*T**2 + 12*b*T - a**2
    Psi_3 = psi_3 // 3

    return Psi_3


if __name__ == "__main__":
    print("------------------------------------------------------------------")
    print("# !!! Locus polynomial for cyclic subgroups of prime order 5 !!! #")
    print("------------------------------------------------------------------\n")
    tilde_Psi_5 = locuspoly_5()    
    print(f"{tilde_Psi_5}\n")
    
    Psi_5 = tilde_Psi_5 // 5
    print("The corresponding normalized polynomial is:")
    print(f"{Psi_5}\n")

    print("------------------------------------------------------------------")
    print("# !!! Locus polynomial for cyclic subgroups of prime order 3 !!! #")
    print("------------------------------------------------------------------")
    print("It corresponds to the normalization of the 3-division polynomial:")
    print(f"{locuspoly_3()}\n")
