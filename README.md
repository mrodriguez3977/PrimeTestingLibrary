# PrimeTestingLibrary
This is a library used to determine if a number is prime or not using the miller-rabin algorithm, and to factor a number suspected to not be prime using pollard rho factorization.

You can change the number of iterations for the miller rabin test in order to be more precise, the more passes, the less chance of the output being wrong, since the miller rabin is a probabilistic test, you should use as many passes as you can in order to prevent false positives. 
