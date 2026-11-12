from tokenize import String
from sage.all import EllipticCurve, GF, Mod, PolynomialRing, is_prime


# defining a finite field 
p128 = 340282366920938463463374607431768211507
F128 = GF(p128)

def kernel_poly_parameters(Fq, order_l):
    """
    Generates an elliptic curve E with l+1 subgroup of order l
    Give the parameters associated to the kernel polynomial
    Return lamda_G an identifier for a subgroup G of order l in E
    """
    
    "ELliptic curve of order divisible by order_l"
    while True:
        delta = 0
        while delta == 0:
            a, b = Fq.random_element(), Fq.random_element()
            delta = 4*(a**3) - 27*(b**2)

        E = EllipticCurve(Fq, [a, b]) # y^2 = x^3 + ax + b over Fq
        if Mod(E.cardinality(), order_l) == 0:
            break    

    "Subgroup G of order l"
    while True:
        Q = E.random_point()
        P = (E.cardinality() // order_l)*Q
        if P != E(0):
            break
    G = [i*P for i in range(order_l)]

    "Isogeny of kernel G"
    phi = E.isogeny(G)
    Kx = phi.kernel_polynomial()    

    "Constructing lamda_G, the scalar parameter of G using finite field extension"  
    Kx_coeffs = Kx.list()[::-1] # highest-degree-first order  
    d = Kx.degree()
    assert d == (order_l-1) // 2   

    R = PolynomialRing(Fq, 'x') 
    x = R.gen()
    f = R.irreducible_element(d)
    Fqd = Fq.extension(f, 'θ')
    θ = Fqd.gen()
    assert f(θ) == 0 

    E_ext = E.base_extend(Fqd)
    Bc = [θ**i for i in range(d)] #Canonical Basis
    # Pack elementary symmetric coeffs: Kx = x^d - s1 x^{d-1} + s2 x^{d-2} - ...
    # Sage stores Kx = x^d + c1 x^{d-1} + c2 x^{d-2} + ..., so s_k = (-1)^k c_k
    lamda = sum(((-1)**(i+1)) * Kx_coeffs[i+1] * θ**i for i in range(d))
    #print(Kx_coeffs, lamda)

    return [E, Kx_coeffs, θ, lamda]


def isogeny_codomain_parameters(init_para, order_l):
    """
    Given an elliptic curve constructed from 'kernel_poly_parameters'
       - Compute the cododmain using appropriate equations
       - Test if E and the codomain E' = E/G are isogenous
    init_para should be equal to kernel_poly_parameters(Fq, order_l)
    """

    E = init_para[0]    
    Kx_coeffs = init_para[1]
    θ, lamda = init_para[2], init_para[3]
    Fq = E.base_field()

    "Definition of the codomain curve isogenous to E"
    if len(Kx_coeffs) > 3:
        c1, c2, c3 = Kx_coeffs[1:4]
    elif len(Kx_coeffs) == 3:
        c1, c2 = Kx_coeffs[1:3]
        c3 = 0
    else:
        c1 = Kx_coeffs[1]
        c2 = c3 = 0

    l = order_l
    a, b = E.a4(), E.a6()
    s1, s2, s3 = -c1, c2, -c3 #To respect elementary symmetric functions signs
    a2 = -30*s1**2 + 60*s2 - (5*l - 6)*a 
    b2 = -70*(s1**3 - 3*s1*s2 + 3*s3) - 42*a*s1 - (14*l - 15)*b 
    E2 = EllipticCurve(Fq, [a2, b2])
    #print(E.cardinality(), E2.cardinality(), E.cardinality() - E2.cardinality())
    #print(E.is_isogenous(E2)) #it works

    return [E, Kx_coeffs, θ, lamda, E2]


"*********************************************************************"
"**************** On the Cyclic-subgroup Locus Polynomial ************"
"*********************************************************************"

"******* Construction/Verification of Psi_5, locus for l = 5 *********"

def locus_polynomial5(Fq, order_l):

    """
    Given E/Fq and G of order order_l verifies the Locus Psi_5
    """

    init_para_0 = kernel_poly_parameters(Fq, order_l)  
    init_para = isogeny_codomain_parameters(init_para_0, order_l)
    E = init_para[0] 
    a, b = E.a4(), E.a6()
    E2 = init_para[4]
    Kx_coeffs = init_para[1]
    c1 = Kx_coeffs[1]

    R = PolynomialRing(Fq, 'x')
    x = R.gen()
    #The intermediate polynomial
    F5x = x**6 + 20*a*x**4 + 160*b*x**3 - 80*a**2*x**2 - 128*a*b*x - 80*b**2
    s1 = -c1
    #print(F5x(s1)) # The verification works
    assert F5x(s1) == 0

    "Definition of the locus polynomial Psi_5"
    θ = init_para[2]
    t, d = θ, order_l // 2 
    Fqd = t.parent()
    Rd = R.base_extend(Fqd)
    x = Rd.gen()

    #This polynomial is computed from the code in locuspoly_construct file
    Psi_5t = (5*x**6 + 38*a*t*x**5 + (55*a**2*t**2 + 100*a + 420*b*t)*x**4
             + (-60*a**3*t**3 - 240*a**2*t + 560*a*b*t**2 - 640*b**2*t**3 + 800*b)*x**3
             + (35*a**4*t**4 - 200*a**3*t**2 + 280*a**2*b*t**3 - 400*a**2 + 160*a*b**2*t**4
                - 320*a*b*t + 720*b**2*t**2)*x**2
             + (-10*a**5*t**5 + 80*a**4*t**3 - 80*a**3*b*t**4 + 160*a**3*t - 80*a**2*b**2*t**5
                - 480*a**2*b*t**2 + 1280*a*b**2*t**3 - 640*a*b - 960*b**3*t**4 + 1200*b**2*t)*x
             + a**6*t**6 - 60*a**5*t**4 + 100*a**4*b*t**5 - 80*a**4*t**2 - 32*a**3*b**2*t**6
             - 320*a**2*b**2*t**4 - 64*a**2*b*t + 640*a*b**3*t**5 - 800*a*b**2*t**2
             - 256*b**4*t**6 + 640*b**3*t**3 - 400*b**2) 
    Psi_5ab = Psi_5t // 5

    lamda = init_para[3]
    #print(Psi_5(lamda)) # The verification works
    assert Psi_5ab(lamda) == 0            

    return [E, Kx_coeffs, θ, lamda, E2, Psi_5ab]


"*********************************************************************"
"****** Trick for reducing the density of the Locus Polynomial *******"
"*********************************************************************"    


"**************** Reduction of Psi_5 through char(Fq) *****************"    

def locus_reduction_usingchar(l_order):

    """
       Psi_5 can be made less dense using appropriate char(Fq)
    """
    
    int_term = {
                 "six": [1], "five": [38], 
                 "four": [11, 20, 84],
                 "three": [12, 48, 112, 128, 160], 
                 "two": [7, 40, 56, 80, 32, 64, 144],
                 "one": [2, 16, 16, 32, 16, 96, 256, 128, 192, 240],
                 "zero": [1, 12, 20, 16, 32, 64, 64, 128, 160, 256, 128, 80]
                 }

    #Looking for the larger integer of int_term
    LI = 0
    terms_tab = []
    for i in range(len(int_term)) :
        dic_elt = list(int_term.items())[i]
        term_i = dic_elt[1] 
        #terms_tab.append(term_i)
        terms_tab.extend(elt for elt in term_i if elt not in terms_tab) #terms_tab+term_i
        for j in term_i:
            if j > LI:
                LI = j
    #terms_tab = terms_tab.sort()
    print(terms_tab)
    shortest_form = []
    good_char = 2
    for l in range(2, LI + 1):
        if is_prime(l):
            reduced_term = []
            for tab_elt in terms_tab:
                if Mod(tab_elt, l) == 0:
                    reduced_term.append(tab_elt)
            if len(reduced_term) >= len(shortest_form):
                shortest_form = reduced_term
                good_char = l

    
    print(f"The shotest form and the corresponding char are:")
    print(f"{shortest_form, good_char}")                    

    return "End Seach"


"******** Reduction of Psi_l throug elliptic curve coeff selection *************"     

def locus_reduction_usingcoeff(Psi_l):
    """
        The terms of Psi_l can be made less complex by curve selection: 
          - Count how many terms has in them, a and b, eliminate the one who appears more 
          - Select a curve such that a = 0 or b = 0.
    """

    #Convert Psi_l in a chain and count the occurences of a and b
    chain_poly = String(Psi_l)
    num_a = 0
    num_b = 0
    for i in len(chain_poly):
        if chain_poly[i] == "a":
            num_a = num_a+1
        elif chain_poly[i] == "b":
            num_b = num_b+1
    print(num_a, num_b)

    if num_a <= num_b:
        print("The coefficient a has less occurences.") 
        return "It is better to select a curve with b = 0."
    else:
        print("The coefficient b has less occurences.")
        return "It is better to select a curve with a = 0."


