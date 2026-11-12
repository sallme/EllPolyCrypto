# Overview
This Python/SageMath repository contains the code and functions associated with our paper entitled 
"EℓℓPoly: On ℓ-Polynomials for Isogeny-Based Cryptography and zk-SNARKs."

The mathematical constructions and related algorithms are developed for cryptographic applications, 
particularly in isogeny-based cryptography and zero-knowledge proofs. To validate the constructions 
presented in the paper, we implement the associated algorithms and test their outputs, including the 
locus polynomials whose roots correspond to the proposed identifiers of cyclic subgroups of order ℓ 
on an elliptic curve.

Using the techniques described in the paper, including finite-field embeddings and properties of elliptic 
curves, the file 'locuspoly_construct.py' provides functions for constructing the locus polynomials. 
Since low-density polynomials are often desirable in cryptographic applications, we focus primarily on the 
cases ℓ=3 and ℓ=5. The resulting universal locus polynomials are then reused and tested in 'main_equations.py', 
where concrete elliptic-curve coefficients are defined.

In accordance with the scope of the paper, the implementation focuses on constructing the proposed polynomials 
and verifying the corresponding cyclic-subgroup and isogeny properties. We therefore do not implement the 
complete R1CS system or integrate it into a zk-SNARK platform.

Finally, the file 'test_equations.py' contains tests for the main functions and displays their outputs. Readers 
can run this file to examine examples involving different elliptic-curve parameters and the corresponding locus 
polynomials. As is customary in isogeny-based cryptography, the tests use finite fields of large characteristic p.
The tests verify not only the construction of the locus polynomials but also that they vanish at the corresponding identifiers of cyclic subgroups of order ℓ. They are performed over different finite fields and elliptic curves, consistently confirming the expected vanishing property established by the mathematical construction. The file 'Test_Screenshots.pdf' provides screenshots illustrating these tests.