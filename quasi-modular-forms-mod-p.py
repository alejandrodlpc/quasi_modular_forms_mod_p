#!/usr/bin/env python
# coding: utf-8

# # Quasi-modular Forms mod $p$

# Here we make $q$-series calculations of modular forms and their reductions modulo some prime $p$.

# ## 0. Notations and sets

# We fix the following notations and sets:
# 1. $B[[q]] = \mathbb{Q}[[q]]$ is the power series ring in the variable $q$.
# 2. $S[P,Q,R] = \mathbb{Q}[x,y,z]$ is the polynomial ring in the variables $x,y,z$.

# In[1]:


B.<q> = PowerSeriesRing(QQ,default_prec=300)
S.<P, Q, R> = PolynomialRing(QQ, 3)


# ## 1. $E_2$, $E_4$, $E_6$ and $\Delta$

# Below, we define the Eisenstein series $E_2$, $E_4$ and $E_6$. Then, using these we define $\Delta$. Finally, we define $E_{p-1}$ for a prime $p\geq5$ as a polynomial $A$ in $E_4$ and $E_6$.

# In[2]:


# Set precision
prec = 1000

# Define the graded ring of modular forms of level 1
M = ModularFormsRing(1)

# The graded ring is generated as an algebra by the Eisenstein series E4 and E6
E2 = 1-24*sum(sigma(n,1)*q^n for n in range(1,prec+1))
E4 = M.q_expansion_basis(4,prec+1)[0]
E6 = M.q_expansion_basis(6,prec+1)[0]

# Using E4 and E6, we define the Delta function
Delta = (E4^3 - E6^2)/1728;

# Print E4, E6 and Delta up to q^r where
r = 5
print('E_2 =',E2.add_bigoh(r))
print('E_4 =',E4.add_bigoh(r))
print('E_6 =',E6.add_bigoh(r))
print('Delta =',Delta.add_bigoh(r))


# We define a general Eisenstein series $E_k$.
# 
# The function `Ek` has:
# - input: positive integers `k`, the weight of $E_k$, and `r`, the precision of the $q$-expansion
# - output: the $q$ expansion of $E_k$ up to $q^{r-1}$

# In[3]:


def Ek(k,r):
    q_expansion = 1-(2*k/bernoulli(k))*sum(sigma(n,k-1)*q^n for n in range(1,r+2))
    return q_expansion.add_bigoh(r)


# In[4]:


Ek(14,6)


# Next, for a prime $p\geq 5$, we determine the polynomial $A_p\in\mathbb{Q}[x,y]$ for which $A_p(E_4,E_6)=E_{p-1}$.
# 
# The function `Ap` has:
# - input: a prime $p$
# - output: a polynomial in $\mathbb{Q}[x,y]$

# ## 2. Basis of $M^{\mathrm{qm}}_{\leq k}(\mathbb{Q})$

# Given a weight $k$, we compute a monomial basis of $M^{\mathrm{qm}}_{\leq k}(\mathbb{Q})$ which is a set of monomials of the form $E_2^a E_4^b E_6^c$ where $2a+4b+6c\leq k$.
# 
# The function `monomial_exponent` has:
# 
# - input: a weight `k`
# - output: a list of triples `[a,b,c]` with $2a+4b+6c \leq k$

# In[5]:


def monomial_exponents(k):
    mon = []
    for a in range(0,k):
        if 2*a<=k:
            for b in range(0,k):
                if 4*b <= k - 2*a:
                    for c in range(0,k):
                        if 6*c <= k-2*a-4*b:
                            mon.append([a,b,c])
                        else: pass
                else: pass
        else: pass
#    mon.pop(0)
    return mon


# For example,

# In[6]:


k = 6
monomial_exponents(k)


# Given a triple $\{a,b,c\}$ of nonnegative integers, we compute $E_2^a E_4^b E_6^b$.
# 
# The function `qmform_from_exp` has:
# - input: a triple `exponents` of the form `[a,b,c]` where $a,b,c\geq 0$
# - output: the $q$-series $E_2^a E_4^b E_6^c$

# In[7]:


def qmform_from_exp(exponents):
    return E2^(exponents[0]) * E4^(exponents[1]) * E6^(exponents[2])


# For example:

# In[8]:


e = [2,3,1]
print('E_2 ^',e[0],'* E_4 ^',e[1],'* E_6 ^',e[2],'has q-series:')
print(qmform_from_exp(e).add_bigoh(6))


# We can write the monomial generating set of $M^{\mathrm{qm}}_{\leq k}(\mathbb{Q})$ as a matrix whose columns are the coefficients of the $q$-expansion of each element of the basis.
# 
# The `basis_matrix` has:
# 
# - input: a weight `k`
# - output: a matrix $(a_{ij})$, where $a_{ij}$ is the $i$th $q$-expansion coefficient of the $j$th basis element.

# In[9]:


def basis_matrix(k):
    mons = monomial_exponents(k)
    basis_length = len(mons)
    basis_elements = [qmform_from_exp(u).add_bigoh(basis_length) for u in mons]
    return Matrix([[m.padded_list()[i] for i in range(0,basis_length)] for m in basis_elements]).transpose()


# For example

# In[10]:


# Fix a weight
k = 6

print(basis_matrix(k))


# Here we check that the `basis_matrix(k)` is invertible for even $k$ up to 20.

# In[11]:


import time
r = 10 # number of even weights
# For r = 5, runtime ~ 0.05 sec
# For r = 10, runtime ~ 0.9 sec
# For r = 15, runtime ~ 10 sec
# For r = 20, runtime ~ 137 sec


# Record the start time
start_time = time.time()

for k in [2*i for i in range(1, r+1)]:
    det_for_k = basis_matrix(k).det()
    print('the monomial generating set of weights <=',k,'is a basis:',det_for_k!=0)

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('----------------------------------------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# Checking the determinants modulo a large prime is much quicker though not completely reliable. Below, we check the determinant of `basis_matrix(k)` modulo $p$ where $p$ is the $(k/2)$th prime (this choice of prime was made empirically by computing the factorizations of the determinants and observing that, roughly, new prime factor appeared in the factorization everytime the weight increased by two).

# In[12]:


# The prime considered is the kth prime, where k is the weight

r = 5
# For r = 5, runtime ~ 0.06 sec
# For r = 10, runtime ~ 0.9 sec
# For r = 15, runtime ~ 7.5 sec
# For r = 20, runtime ~ 45 sec
# For r = 25, runtime ~ 201 sec

# Record the start time
start_time = time.time()

for k in [2*i for i in range(1, r+1)]:
    p = nth_prime(2*k)
    A = basis_matrix(k).change_ring(GF(p))
    det_modp = det(A)
    print('for weight',k,', the monomial generating set is a basis mod p =',p,'-->',det_modp !=0)

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('----------------------------------------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# We can write the generating set of $M^{\mathrm{qm}}_{\leq k}(\mathbb{Q})$ as a matrix whose columns are the coefficients of the $q$-expansion of each element of the monomial generating set.
# 
# The `basis_matrix` has:
# 
# - input: a weight `k` and a positive integer `rows`
# - output: a `rows`$\times N$ matrix $(a_{ij})$, where $N$ is the number of monomial generators and $a_{ij}$ is the $i$th $q$-expansion coefficient of the $j$th monomial.

# In[13]:


def generating_matrix(k,rows):
    mons = monomial_exponents(k)
    basis_length = len(mons)
    basis_elements = [qmform_from_exp(u).add_bigoh(rows) for u in mons]
    return Matrix([[m.padded_list()[i] for i in range(0,rows)] for m in basis_elements]).transpose()


# In[14]:


k = 6
rows = 10

print(generating_matrix(k,rows))


# Given an even weight $2k$, we have $$\#\{(a,b,c):2a+4b+6c\leq2k\}=\left\lceil \frac{2k^3+21k^2+66k+30}{72} \right\rceil$$ The function `number_of_monomials` gives this formula.

# In[15]:


def number_of_monomials(w):
    if w%2 ==0:
        k = w/2
    else:
        k = (w-1)/2
    return ceil((2*k^3+21*k^2+66*k+30)/72)


# In[16]:


k = 96
print('For weight =',k,',')
print('using monomial_exponents gives us',len(monomial_exponents(k)),'many monomials,')
print('using the above formula gives us',number_of_monomials(k),'many monomials.')


# Notice that the function `number_of_monomials` is much faster than computing `len(monomial_exponents(k))`

# In[17]:


w = 400

repetitions = 1

import timeit
from math import ceil

# Define a lambda function to include the calculation
calculation_1 = lambda: number_of_monomials(w)
calculation_2 = lambda: len(monomial_exponents(w))

execution_time_1 = timeit.timeit(calculation_1, number=repetitions)
execution_time_2 = timeit.timeit(calculation_2, number=repetitions)

print('weight =',w)
print("Execution time using number_of_monomials:", execution_time_1, "seconds")
print("Execution time using len(monomial_exponents):", execution_time_2, "seconds")


# ## 3. Generating quasi-modular forms from $\{E_2^a E_4^b E_6^c:a,b,c\geq0\}$

# Given a quasi-modular form $f\in M^{\mathrm{qm}(\mathbb{Q})}$ as a $q$-series, we express it as a $\mathbb{Q}$-linear combination of the monomial generating set computed in $\S 2$.
# 
# The function `linear_combination_with_prec` has:
# 
# - input: a positive even weight `k`, a positive integer `pr` and a $q$-series `f` that is a sum of quasi-modular forms of weight bounded by $k$,
# 
# - output: a list of numbers $\{a_1,\ldots,a_n\}$ such that $f=a_1 v_1+\cdots+a_n v_n$ where $\{v_1,\ldots,v_n\}$ is the ordered generating set computed with `monomial_basis`.

# In[18]:


def linear_combination_with_prec(k,pr,f):
    mons = monomial_exponents(k)
    basis_length = len(mons)
    if f.common_prec(f) < pr:
        F = f.lift_to_precision(pr)
    else:
        F = f.add_bigoh(pr)
    list_of_coefficients = F.padded_list(pr)
    v = vector([c for c in list_of_coefficients[:pr]])
    A = generating_matrix(k,pr)
    try:
        sol = A.solve_right(v)
    except ValueError:
        return "Error"
    return sol 


# In[19]:


k=12
mons = monomial_exponents(k)
print(len(mons))

pr = 29

f = Ek(k,pr);
w = linear_combination_with_prec(12,pr,f)

print(w)

#pol = sum(w[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)

#print(pol)

#print(f-pol(E2,E4,E6))


# Given a quasi-modular form $f\in M^{\mathrm{qm}(\mathbb{Q})}$ as a $q$-series, we express it as a $\mathbb{Q}$-linear combination of the monomial basis computed in $\S 2$.
# 
# The function `linear_comb` has:
# 
# - input: a positive even weight `k` and a $q$-series `f` that is a sum of quasi-modular forms of weight bounded by $k$,
# 
# - output: a list of numbers $\{a_1,\ldots,a_n\}$ such that $f=a_1 v_1+\cdots+a_n v_n$ where $\{v_1,\ldots,v_n\}$ is the ordered basis computed with `monomial_basis`.

# In[20]:


def linear_comb(k,f):
    basis_length = number_of_monomials(k)
    F = f.add_bigoh(basis_length)
    list_of_coefficients = F.padded_list(basis_length)
    v = vector([c for c in list_of_coefficients[:basis_length]])
    A = basis_matrix(k)
    return list(A^(-1)*v)


# Given a linear combination of monomials from `linear_comb`, we can check whether the coefficients are $p$-integral.
# 
# The function `p_integral_linear_comb` has:
# - input: a positive even weight `k`, $q$-series `f` that is a sum of quasi-modular forms of weight bounded by $k$, and a prime `p`
# 
# - output: a list of numbers $\{d_1,\ldots,d_n\}$ where $d_1,\ldots,d_n$ are the reduction modulo $p$ of the denominators of the coefficients of the linear combination given by `linear_comb(k,f)`. 

# In[21]:


def p_integral_linear_comb(k,f,p):
    return [denominator(u)%p for u in linear_comb(k,f)]


# Given a positive even weight $k$, we use `linear_comb` to produce the polynomial in $E_2,E_4$ and $E_6$ that equals $f$.
# 
# The function `associated_poly` has:
# - input: a positive even weight `k` and a $q$-series `f` that is a sum of quasi-modular forms of weight bounded by $k$,
# - output: a polynomial in $P$, $Q$ and $R$ that equals `f` when evaluated at $P=E_2$, $Q=E_4$ and $R=E_6$.

# In[22]:


def associated_poly(k,f):
    mons = monomial_exponents(k)
    v = linear_comb(k,f)
    return sum(v[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)


# In[23]:


def associated_poly_with_prec(k,pr,f):
    mons = monomial_exponents(k)
    try:
        v = linear_combination_with_prec(k,pr,f)
    except ValueError:
        return "Error"
    else:
        return sum(v[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)


# For example:

# In[24]:


p = 19 # fixed prime number
k = p-1 # weight of the form
r = 100 # max order when printing q-expansions

f = Ek(k,prec);
pol = associated_poly(k,f)

print('Let:')
print('f = ',f.add_bigoh(10))
print('---------------------------')
print('f is equal to the polynomial:')
print(pol)
print('---------------------------')
print('The difference f-polynomial is:')
print(f-pol(E2,E4,E6))


# We combine the above functions into the single function `fancy_linear_comb` that has:
# - inputs: `k` a highest weight, `f` a $q$-series, `p` a prime, `r` a precision for $q$-series
# - ouput: a verbose analysis of `f`

# In[25]:


def fancy_linear_comb(k,f,p,r):
    mons = monomial_exponents(k)
    basis_length = len(mons)
    #basis_elements = [qmform_from_exp(u).add_bigoh(basis_length) for u in mons]
    v = linear_comb(k,f)
    associated_poly = sum(v[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)
    denominators_mod_p = p_integral_linear_comb(k,f,p)
    
    if 0 in denominators_mod_p:
        p_integral_check = False
    else:
        p_integral_check = True

    print('Consider the form:')
    print('f =',f.add_bigoh(r))
    print('-----------------------------------------------')
    print('The monomial basis for weights <=',k,"is:")
    print([P^e[0]*Q^e[1]*R^e[2] for e in mons])
    print('-----------------------------------------------')
    print('f is equal to the following linear combination of the above monomial basis:')
    print(associated_poly)
    print('-----------------------------------------------')
    print('The coefficients of the linear combination are:')
    print(v)
    print('-----------------------------------------------')
    print('The denominators modulo',p,'of the denominators of these coefficients are:')
    print(denominators_mod_p)
    print('-----------------------------------------------')
    print('The linear combination is',p,'- integral:',p_integral_check)


# For example: we compute the Eisenstein series $E_{p-1}$ as a linear combination of monomials of the form $E_2^a E_4^b E_6^c$.

# In[26]:


p = 13 # fixed prime number
k = p-1 # weight of the form
r = 5 # max order when printing q-expansions

basis_length = len(monomial_exponents(k))

f = Ek(k,basis_length); #Define the form

fancy_linear_comb(k,f,p,r)


# ## $\S$4. Finding Graded Components

# Given a prime $p$ and a set $X$ of triples $(a,b,c)$, we partition this set into
# 
# $$X = \bigsqcup_{i=0}^{p-2} X_i\quad\text{with}\quad X_i :=\{(a,b,c)\in X \mid 2a+4b+6c\equiv i \pmod{p-1}\}$$
# 
# The function `graded_exponents_mod_p` takes
# - input: A set `X` of triples and a prime number `p`
# - output: A set $X_0,X_1,\ldots,X_{p-2}$ of subsets of `X` where $(a,b,c)\in X_i \Longleftrightarrow 2a+4b+6c\equiv i\pmod{p-1}$

# In[27]:


def graded_exponents_mod_p(X,p):
    parts = []
    for u in [l for l in range(0,p-2) if l%2 == 0]:
        parts_u = [w for w in X if (2*w[0]+4*w[1]+6*w[2]-u)%(p-1) == 0]
        parts.append(parts_u)
    return parts


# In[28]:


p = 13
k = 12 # bound for weight

mons = monomial_exponents(k)

print('-----------------------------------------------------------------------')
print('The set of monomials of weight <=',k,'is')
print([P^e[0]*Q^e[1]*R^e[2] for e in mons])

print('-----------------------------------------------------------------------')
    
g_exp = graded_exponents_mod_p(mons,p)
print('The set of monomials partitioned according to their weights modulo p-1:')

for ge in g_exp:
    print('k =',2*g_exp.index(ge),'-->',[P^t[0]*Q^t[1]*R^t[2] for t in ge])
print('-----------------------------------------------------------------------')


# In[128]:


p = 13
k = 8 # bound for weight

mons = monomial_exponents(k)

print('-----------------------------------------------------------------------')
print('The set of monomials of weight <=',k,'is')
print([P^e[0]*Q^e[1]*R^e[2] for e in mons])

print('-----------------------------------------------------------------------')
    
g_exp = graded_exponents_mod_p(mons,p)

print('The set of monomials partitioned according to their weights modulo p-1:')

for ge in g_exp:
    print('k =',2*g_exp.index(ge),'-->',[P^t[0]*Q^t[1]*R^t[2] for t in ge])
print('-----------------------------------------------------------------------')


# ## $\S$5. The Quasi-modular Forms $U_a$

# In this section we define the $q$-series $U^{(\beta)}_a$ where $a>1$ and $\beta$ is an ordered list of $(a-1)$ many 0's and 1's. 

# First, we convert an ordered list of 0's and 1's into a string of inequalities.
# 
# The function `vector_to_ineqs` has:
# - input: a list `vec` of 0's and 1's.
# - output: a list of strings with "$\geq$" in the place of the 0's and "$>$" in the place of the 1's.
# 
# For example, the vector `[0,0,1,1,0]` is mapped to `['>=','>=','>','>','>=']`

# In[30]:


inequality_signs = ['>=','>']

def vector_to_ineqs(vec):
    list_of_ineqs = []
    for v in vec:
        if v==0:
            list_of_ineqs.append('>=')
        elif v==1:
            list_of_ineqs.append('>')
        else:
            pass
    return list_of_ineqs


# For example,

# In[31]:


vec = [0,0,0,1,1,0,1,0,1,1]
print('beta =',vec)
print('associated list of inequalities =',vector_to_ineqs(vec))


# Given a partition $[\lambda_1,\ldots,\lambda_a]$ of $n$ of length $a$, we check whether this partition satisfies the list of inequalities given by a vector $\beta$ with $(a-1)$ entries of 0's and 1's. More precisely:
# 
# The function `check_partition_vs_ineqs` has:
# - input: a `partition` of $n$ of length $a$, and a list `choice_of_ineqs` of $a-1$ 0's and 1's
# - output: `True` or `False` depending on whether the partition $\{\lambda_1,\ldots,\lambda_a\}$ satisfies the given list of inequalities.

# In[32]:


def check_partition_vs_ineqs(partition,choice_of_ineqs):
    import itertools
    partition_str = [str(p) for p in partition]
    list_of_zipped_strings = [x for x in itertools.chain.from_iterable(itertools.zip_longest(partition_str,list(choice_of_ineqs))) if x]
    joined_string = "".join(list_of_zipped_strings)
    return eval(joined_string)


# For example:

# In[33]:


n = 11
a = 4

import random
vec = [random.choice([0, 1]) for _ in range(a-1)]
choice_of_ineqs = vector_to_ineqs(vec)

parts = Partitions(n,length=a)

print('list of inequalities:',choice_of_ineqs)

for p in parts:
    print('The partition',p,'satisfies the inequalites',choice_of_ineqs,':',check_partition_vs_ineqs(p,choice_of_ineqs))


# Given a number $n$, a length $a$ and a list $\beta$ of 0's and 1's, we compute all the partitions of $n$ of size $a$ that satisfy the inequalities associated to $\beta$.
# 
# The function `partitions_from_ineqs` has:
# - input: a number `n`, a length `a` and a list `vec` of $a-1$ 0's and 1's
# - output: the complete list of partitions of $n$ of length $a$ that satisfy the inequalities associated to `vec`.

# In[34]:


def partitions_from_ineqs(n,a,vec):
    parts = Partitions(n,length=a)
    choice_of_ineqs = vector_to_ineqs(vec)
    return [p for p in parts if check_partition_vs_ineqs(p,choice_of_ineqs)]


# For example:

# In[35]:


n = 11
a = 4

import random
vec = [random.choice([0, 1]) for _ in range(a-1)]
choice_of_ineqs = vector_to_ineqs(vec)

print('list of inequalities:',choice_of_ineqs)
print('partitions of',n,'of length',a,'that satisfy',choice_of_ineqs,':')
for p in partitions_from_ineqs(n,a,vec):
    print(p)
    
print(list(partitions_from_ineqs(n,a,vec)[0]))


# The $q$-series $U^{(\beta)}_a$ is defined as
# $$U^{(\beta)}_a (q) = \sum_{k} \frac{q^{|k|}}{(1-q^k)^2}$$
# where we are using multindex notation, i.e.
# $$k=(k_1,\ldots,k_a),\quad |k|=k_1+\cdots+k_a,\quad (1-q^k)^2=(1-q^{k_1})^2\cdots(1-q^{k_a})^2$$

# The function `U` has
# - input: `a` an positive integer, `vec` a list of $a-1$ 0's and 1's and `r` a positive integer
# - output: the $q$-series $U^{\beta}_a (q)$ up to $O(q^r)$.

# In[36]:


B.<q> = PowerSeriesRing(QQ,default_prec=200)

def U_term(part_index):
    return q^(sum(i for i in part_index))*(prod((1-q^u)^(-2) for u in part_index))

def U(a,vec,r):
    choice_of_ineqs = vector_to_ineqs(vec)
    index_set = []
    for n in range(0,r+1):
        for p in partitions_from_ineqs(n,a,vec):
            index_set.append(list(p))
    return sum(U_term(k) for k in index_set).add_bigoh(r)


# For example:

# In[37]:


a = 4
vec = [0,1,0]
r = 8

choice_of_ineqs = vector_to_ineqs(vec)
print('list of inequalities:',choice_of_ineqs)

index_set = []
for n in range(0,r+1):
    print('-----------------------------------------------------------------------------')
    print('Partitions of',n,'of length',a,'satisfying',choice_of_ineqs,'is:')
    for p in partitions_from_ineqs(n,a,vec):
        print(p,'-->',U_term(p).add_bigoh(r+2))
        index_set.append(list(p))
print('-----------------------------------------------------------------------------')       
print('The q-series U is:',U(a,vec,r))


# In[38]:


a = 4
vec = [0,0,0]
r = 20
f = U(a,vec,r)
pol = associated_poly(2*a,f)

print('Size of monomial basis:',len(monomial_exponents(2*a)))
print('---------------------------')
print('Let:')
print('f = ',f)
print('---------------------------')
print('The associated polynomial is:')
print('F =',pol)
print('---------------------------')
print('The difference f-F is:')
print(f-pol(E2,E4,E6))


# In[ ]:


k = 0 --> [1]
k = 2 --> [P] 5/172032
k = 4 --> [Q, P^2] 37/5529600, 37/2211840
k = 6 --> [R, P*Q, P^3] 1/870912 , 1/331776 , 5/1990656
k = 8 --> [Q^2, P*R, P^2*Q, P^4]
k = 10 --> []


# In[135]:


f = -27859/464486400+37/5529600*E4+37/2211840*E2^2
f.add_bigoh(5)


# To test the quasi-modularity of $U^{(\beta)}_a$ for all binary lists $\beta$, we have the following function that lists all possible such $\beta$.

# In[39]:


def generate_binary_lists(length):
    all_binary_lists = []
    total_possibilities = 2**length
    for i in range(total_possibilities):
        binary_str = bin(i)[2:].zfill(length)  # Convert to binary and pad with zeros
        binary_list = [int(bit) for bit in binary_str]
        all_binary_lists.append(binary_list)
    
    return all_binary_lists


# For example,

# In[40]:


a = 4
generate_binary_lists(a-1)


# Given $a$, there are $2^{a-1}$ different $U^{(\beta)}_a$. For each of these, we compute the associated polynomial and the difference between the evaluated associated polynomial and $U^{(\beta)}_a$.

# In[41]:


a = 4
all_vectors = generate_binary_lists(a-1)
mon_size = len(monomial_exponents(2*a))
r = mon_size +10

print('Size of monomial generating set:',mon_size)

for vec in all_vectors:
    f = U(a,vec,r)
    pol = associated_poly(2*a,f)
    print('-----------------------------')
    print('For U _',a,'^',vec,'=',f.add_bigoh(f.valuation()+1),', then U - associated_poly(U) =',(f-pol(E2,E4,E6)).add_bigoh(mon_size+1))


# In[42]:


def check_all_U_for_fixed_a(a):
    all_vectors = generate_binary_lists(a-1)
    mon_size = len(monomial_exponents(2*a))
    r = mon_size + 1
    quasi_modular_forms = []
    for vec in all_vectors:
        f = U(a,vec,r)
        pol = associated_poly(2*a,f)
        if (f-pol(E2,E4,E6)).truncate(mon_size+1)==0:
            quasi_modular_forms.append(vec)
        else:
            pass
    return quasi_modular_forms

def check_all_U_for_fixed_a_verbose(a):
    all_vectors = generate_binary_lists(a-1)
    mon_size = len(monomial_exponents(2*a))
    r = mon_size + 1
    print('Size of monomial generating set:',mon_size)
    print('--------------------------------------')
    for vec in all_vectors:
        f = U(a,vec,r)
        pol = associated_poly(2*a,f)
        print('For U _',a,'^',vec,'=',f.add_bigoh(f.valuation()+1),', then U - associated_poly(U) =',(f-pol(E2,E4,E6)).add_bigoh(mon_size+1))
        print('--------------------------------------')
    return ''


# In[43]:


a = 4

print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')
for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

print('------------------------------------------------')
print('Long Answer:')
print('--------------------------------------')
print(check_all_U_for_fixed_a_verbose(a))


# In[44]:


a = 5

print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')
for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

print('------------------------------------------------')
print('Long Answer:')
print('--------------------------------------')
print(check_all_U_for_fixed_a_verbose(a))


# In[45]:


#a = 6

#print('Short Answer:')
#print('--------------------------------------')
#print('Inequalities for which U _',a,'is quasi-modular:')
#for x in check_all_U_for_fixed_a(a):
#    print(vector_to_ineqs(x))
#
#print('------------------------------------------------')
#print('Long Answer:')
#print('--------------------------------------')
#print(check_all_U_for_fixed_a_verbose(a))


# In[46]:


# <1 sec for a=2
a = 2
import time


print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
start_time = time.time()

for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('--------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[47]:


# <1 sec for a=3
a = 3
import time


print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
start_time = time.time()

for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('--------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[48]:


# <1 sec for a=4
a = 4
import time


print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
start_time = time.time()

for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('--------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[49]:


# <1 sec for a=5
a = 5
import time


print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
start_time = time.time()

for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('--------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[50]:


# ~4.5 sec for a=6
a = 6
import time


print('Short Answer:')
print('--------------------------------------')
print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
start_time = time.time()

for x in check_all_U_for_fixed_a(a):
    print(vector_to_ineqs(x))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('--------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[51]:


# ~29 sec for a=7
#a = 7
#import time

#print('Short Answer:')
#print('--------------------------------------')
#print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
#start_time = time.time()

#for x in check_all_U_for_fixed_a(a):
#    print(vector_to_ineqs(x))

# Record the end time
#end_time = time.time()

# Calculate and print the elapsed time
#elapsed_time = end_time - start_time
#print('--------------------------------------')
#print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[52]:


# ~260 sec for a=8
#a = 8
#import time

#print('Short Answer:')
#print('--------------------------------------')
#print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
#start_time = time.time()

#for x in check_all_U_for_fixed_a(a):
#    print(vector_to_ineqs(x))

# Record the end time
#end_time = time.time()

# Calculate and print the elapsed time
#elapsed_time = end_time - start_time
#print('--------------------------------------')
#print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[53]:


# ~2972 sec ~ 50min for a=9
#a = 9
#import time


#print('Short Answer:')
#print('--------------------------------------')
#print('Inequalities for which U _',a,'is quasi-modular:')

# Record the start time
#start_time = time.time()

#for x in check_all_U_for_fixed_a(a):
#    print(vector_to_ineqs(x))

# Record the end time
#end_time = time.time()

# Calculate and print the elapsed time
#elapsed_time = end_time - start_time
#print('--------------------------------------')
#print(f"Elapsed time: {elapsed_time:.6f} seconds")


# ---
# 
# ## Scratch

# In[54]:


print('1 =',1)
print('E2 =',E2.add_bigoh(5))
print('E4 =',E4.add_bigoh(5))
print('E2^2 =',(E2*E2).add_bigoh(5))
print('E6 =',E6.add_bigoh(5))
print('E2*E4 =',(E2)*(E4).add_bigoh(5))
print('E2^3 =',(E2^3).add_bigoh(5))


# In[106]:


k = 10; print('Diagonal Basis for Quasimodular forms of mixed weight <=',k)
d = number_of_monomials(k); print('Dimension of the space =',d)

mons_exp = monomial_exponents(k)
mons = [qmform_from_exp(u) for u in mons_exp]

A = basis_matrix(k)
all_primes_in_all_denominators = []

def flatten(x):
    return [item for sublist in x if sublist for item in sublist]

print('----------------------------------------------------')

for r in range(0,d):
    v = vector([0]*d)
    v[r] = 1
    sol = A.solve_right(v)
    q_series = sum(sol[i]*mons[i] for i in range(0,d))
    dens = set([denominator(a) for a in q_series.padded_list()])
    primes_in_denominator = []
    
    for w in dens:
        primes_in_denominator.append(prime_factors(w))
    
    flattened_primes_in_denominators = flatten(primes_in_denominator)
    if flattened_primes_in_denominators == []:
        all_primes_in_dens = []
    else:
        all_primes_in_dens = list(set(flattened_primes_in_denominators));
        all_primes_in_all_denominators.append(all_primes_in_dens)
    #print(q_series.add_bigoh(d+1))
    #print('denominators =',dens)
    #print('prime factors of the denominators =',all_primes_in_dens)
    #print('----------------------------------------------------')
    
print('all prime factors in all denominators:',list(set(flatten(all_primes_in_all_denominators))))


# In[110]:


def flatten(x):
    return [item for sublist in x if sublist for item in sublist]

def primes_in_den_weight_k(k):
    d = number_of_monomials(k);
    mons_exp = monomial_exponents(k)
    mons = [qmform_from_exp(u) for u in mons_exp]

    A = basis_matrix(k)
    all_primes_in_all_denominators = []

    for r in range(0,d):
        v = vector([0]*d)
        v[r] = 1
        sol = A.solve_right(v)
        q_series = sum(sol[i]*mons[i] for i in range(0,d))
        dens = set([denominator(a) for a in q_series.padded_list()])
        primes_in_denominator = []
    
        for w in dens:
            primes_in_denominator.append(prime_factors(w))
    
        flattened_primes_in_denominators = flatten(primes_in_denominator)
        if flattened_primes_in_denominators == []:
            all_primes_in_dens = []
        else:
            all_primes_in_dens = list(set(flattened_primes_in_denominators));
            all_primes_in_all_denominators.append(all_primes_in_dens)
    return list(set(flatten(all_primes_in_all_denominators)))


# In[111]:


primes_in_den_weight_k(16)


# In[127]:


max_weight = 22
list_of_primes = []

for w in range(2,max_weight+1,2):
    list_of_primes.append(primes_in_den_weight_k(w))

final_list = sorted(list(set(flatten(list_of_primes))))
print(final_list)


# In[58]:


def theta(g):
    return q * derivative(g,q,1)

def iter_theta(g,n):
    for _ in range(1, n+1):
        g = theta(g)
    return g


# In[59]:


k = 8; print('weight <=',k)
r = number_of_monomials(k); print('upper bound on the dimension =',r)
p = 13; print('prime =',p)

print('------------------------------------------------------------')

f = E2^3 + E2*E6 + E2*E4 + E2; print('f = ')
print(f.add_bigoh(r+1))

print('------------------------------------------------------------')

fbar = f.change_ring(GF(p)); print('f modulo',p,'=')
print(fbar.add_bigoh(r+1))

print('------------------------------------------------------------')

for n in range(1,3):
    nth_der_of_f = iter_theta(f,n)
    nth_der_of_f_modp = nth_der_of_f.change_ring(GF(11))
    print('derivative #',n,'of fbar is:')
    print(nth_der_of_f_modp.add_bigoh(r+1))
    print('------------------------------------------------------------')


# In[ ]:


def Up(g,p):
    B.<q> = PowerSeriesRing(QQ,default_prec=200)
    coef_list = g.padded_list()
    new_coef_list = coef_list[0::p]
    return B(new_coef_list)

def Tn(k,g,ell):
    return Up(g,ell) + ell^(k-1) * (g.V(ell))


# In[ ]:


def det_mod_p(A,p):
    A_modp = A.change_ring(GF(p))
    return det(A_modp)


# In[60]:


k = 8; print('weight <=',k)
r = number_of_monomials(k); print('upper bound on the dimension =',r)
p = 3; print('prime =',p)

print('------------------------------------------------------------')

f = E2^3 + E2*E6 + E2*E4 + E2; print('f = ')
print(f.add_bigoh(r+1))

print('------------------------------------------------------------')

print('U_p(f) =')
print(Up(f,p).add_bigoh(r+1))

print('------------------------------------------------------------')

print('V_p(f) =')
print(f.V(p).add_bigoh(r+1))
print('------------------------------------------------------------')

print('T_5(f) =')
print(Tn(8,f,5).add_bigoh(r+1))


# In[ ]:


r = 10

# Record the start time
start_time = time.time()

for k in [2*i for i in range(1, r+1)]:
    A = basis_matrix(k)
    i = 1
    while det_mod_p(A,nth_prime(i)) == 0:
        i += 1
    print('for weights <=',k,', the first prime with nonzero determinant is p =',nth_prime(i),', with deteriminant =',det_mod_p(A,nth_prime(i)))

# Record the end time
end_time = time.time()

# Calculate and print the elapsed time
elapsed_time = end_time - start_time
print('----------------------------------------------------------------------')
print(f"Elapsed time: {elapsed_time:.6f} seconds")


# In[ ]:


f = (1/3)*q
f.change_ring(GF(5))


# In[ ]:


k = 8; print('weight <=',k)
r = number_of_monomials(k); print('upper bound on the dimension =',r)
p = 5; print('prime =',p)

print('------------------------------------------------------------')

f = E2^3 + E2*E6 + E2*E4 + E2; print('f = ')
print(f.add_bigoh(r+1))
print('f = ',associated_poly(k,f))

print('------------------------------------------------------------')

f = fbar = f.change_ring(GF(p)); print('f modulo',p,'=')
print(fbar.add_bigoh(r+1))


print('------------------------------------------------------------')


list_of_coefficients = fbar.padded_list(r)
v = vector([c for c in list_of_coefficients])
A = basis_matrix(k).change_ring(GF(p))
try:
    sol = A.solve_right(v)
    print(sol)
except ValueError:
    print("Error")
    
print('------------------------------------------------------------')

solQQ = [u.lift() for u in sol]
mons = monomial_exponents(k)
poly = sum(solQQ[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)

print(poly)


# In[ ]:


def associated_poly(k,f):
    mons = monomial_exponents(k)
    v = linear_comb(k,f)
    return sum(v[mons.index(e)]*P^e[0]*Q^e[1]*R^e[2] for e in mons)


# In[ ]:





# In[ ]:


a = 3
vec = [0,0]
r = 20
f = U(a,vec,r)

p = 7; print('Let p =',p) # fixed prime number
k = 2*a # weight of the form

print('-----------------------------------------------------------------------')

fancy_linear_comb(k,f,p,r)

mons = monomial_exponents(k)

print('-----------------------------------------------------------------------')
    
g_exp = graded_exponents_mod_p(mons,p)
print('The set of monomials partitioned according to their weights modulo p-1:')

for ge in g_exp:
    print('k =',2*g_exp.index(ge),'-->',[P^t[0]*Q^t[1]*R^t[2] for t in ge])
print('-----------------------------------------------------------------------')

poly = associated_poly(k,f)

for ge in g_exp:
    graded_poly = sum(poly.coefficient({P:m[0], Q:m[1], R:m[2]})*P^(m[0])*Q^(m[1])*R^(m[2]) for m in ge)
    print('k =',2*g_exp.index(ge),'-->',graded_poly)
print('-----------------------------------------------------------------------')

for ge in g_exp:
    graded_poly = sum(poly.coefficient({P:m[0], Q:m[1], R:m[2]})*P^(m[0])*Q^(m[1])*R^(m[2]) for m in ge)
    print('k =',2*g_exp.index(ge),'-->',graded_poly(E2,E4,E6).add_bigoh(a+1))
print('-----------------------------------------------------------------------')


# In[ ]:


640%7


# In[ ]:


a = 4
vec = [0]*(a-1)
r = 20
f = U(a,vec,r)

p = 7; print('Let p =',p) # fixed prime number
k = 2*a # weight of the form

print('-----------------------------------------------------------------------')

fancy_linear_comb(k,f,p,r)

mons = monomial_exponents(k)

print('-----------------------------------------------------------------------')
    
g_exp = graded_exponents_mod_p(mons,p)
print('The set of monomials partitioned according to their weights modulo p-1:')

for ge in g_exp:
    print('k =',2*g_exp.index(ge),'-->',[P^t[0]*Q^t[1]*R^t[2] for t in ge])
print('-----------------------------------------------------------------------')

poly = associated_poly(k,f)

for i in (0..a):
    f_i = poly.coefficient(P^i)
    print('f_',i,'=',f_i)
print('-----------------------------------------------------------------------')

for i in (0..a):
    fi = poly.coefficient(P^i)(E2,E4,E6)
    print('f_',i,'=',B(fi).add_bigoh(5))
print('-----------------------------------------------------------------------')


# In[ ]:


A = basis_matrix(4)
print(A,'\n')

print('det A=',factor(det(A)))

Abar = Matrix([[1/det(A),1,1,1],[0,240,-24,-48],[0,2160,-72,432],[0,6720,-96,3264]])

v = vector([0,0,1,0])

sol = Abar^(-1)*v

print(Abar^(-1)*v)

q_series = sol[0]*1/det(A) + sol[1]*E4 + sol[2]*E2 + sol[3]*E2^2

print(q_series.add_bigoh(20))

#B = Matrix([[1,1/240,1/24,1/28],[0,1,-1,-1],[0,9,3,9],[0,28,-4,68]])

#print(Abar^(-1))

#print(B)

#print('det B =',factor(det(B)))

#print(B^(-1))


# In[ ]:


2^15 * 3^5 * 5


# In[ ]:


4064256-7/1440-5/48+1/144


# In[ ]:


for k in (1..4):
    A = basis_matrix(2*k)
    d = det(A)
    print('for weights <=',2*k,'the determinant is:',factor(abs(d)))


# In[ ]:


def ff(x,*y):
    return print(x,y[0],y[1],y[2])


# In[ ]:


def qseries_matrix(r,li):
    return Matrix([g.padded_list(r) for g in li]).transpose()


# In[ ]:


k = 4
i = 1 # canonical basis vector index
r = number_of_monomials(k)

print('weights up to:',k,', the space has dimension',r)
print('------------------------------------------------------')

B.<q> = PowerSeriesRing(QQ,default_prec=300)
d = 1/det(basis_matrix(k))

mons = monomial_exponents(k) # exponents of monomial basis
mons.pop(0) # remove the constant monomial
monomial_basis = [B(d)]+[qmform_from_exp(m) for m in mons]

A = qseries_matrix(len(monomial_basis),monomial_basis)

zeros = [0]*r
zeros[i-1] = 1
v = vector(zeros) # ith canonical basis vector

sol = A^(-1)*v

print('linear combination of modified monomial basis required to generate q ^',i-1,'+ ... is:\n')
print(sol)
print('\n-----------------------------------------------------')

diag_q_series = sum(sol[i]*monomial_basis[i] for i in (0..r-1))

print('resulting q-series:\n')
print(diag_q_series.add_bigoh(r+1))


# In[ ]:


B.<q> = PowerSeriesRing(QQ,default_prec=300)

for l in (1..7):
    k = 2*l
    r = number_of_monomials(k)

    print('------------------------------------------------------------------')
    print('Diagonal basis for QM forms of weights <=',k,'of dimension',r,':\n')
    d = 1/det(basis_matrix(k))

    mons = monomial_exponents(k) # exponents of monomial basis
    mons.pop(0) # remove the constant monomial
    monomial_basis = [B(d)]+[qmform_from_exp(m) for m in mons]

    A = qseries_matrix(len(monomial_basis),monomial_basis)

    for i in (1..r):
        zeros = [0]*r
        zeros[i-1] = 1
        v = vector(zeros) # ith canonical basis vector

        sol = A^(-1)*v

        diag_q_series = sum(sol[i]*monomial_basis[i] for i in (0..r-1))

        print(diag_q_series.add_bigoh(r+2))


# In[ ]:




