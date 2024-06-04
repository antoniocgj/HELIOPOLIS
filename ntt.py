from univariate import *
from math import log2, floor

def ntt( primitive_root, values ):
    assert(len(values) & (len(values) - 1) == 0), "cannot compute ntt of non-power-of-two sequence"
    if len(values) <= 1:
        return values

    half = len(values) // 2

    odds = ntt(primitive_root**2, values[1::2])
    evens = ntt(primitive_root**2, values[::2])

    return [evens[i % half] + (primitive_root**i) * odds[i % half] for i in range(len(values))]


def _ntt_rk(p, w_n, depth, m, ntt_max_depth):
    n = len(p)
    if(depth == ntt_max_depth):
        p2 = [0]*n
        for i in range(n):
            for j in range(n):
                p2[i] += p[j]*(w_n**(i*j))
        return p2
    w_m = (w_n**(n//m))
    h = [0]*m
    for i in range(m):
        g = []
        for j in range(n//m):
            g.append(p[m*j + i])
        h[i] = _ntt_rk(g, (w_n**m), depth + 1, m, ntt_max_depth)
    
    p2 = [0]*n
    for k1 in range(n//m):
        for k2 in range(m):
            for i in range(m):
                p2[k1 + (n//m)*k2] += h[i][k1] * ((w_n**(i*k1)) * (w_m**(i*k2)))
    return p2


def _ntt_rk2(p, ws, w_power, depth, maskMod, m, ntt_max_depth):
  n = len(p)
  if n == 1: return p
  if(depth == ntt_max_depth or (n%m)):
    p2 = [FieldElement(0)]*n
    for i in range(n):
        for j in range(n):
            p2[i] += p[j]*ws[(w_power*i*j)&maskMod]
    return p2
  h = [0]*m
  for i in range(m):
    g = []
    for j in range(n//m):
      g.append(p[m*j + i])
    h[i] = _ntt_rk2(g, ws, m*w_power, depth + 1, maskMod, m, ntt_max_depth)
  p2 = [FieldElement(0)]*n
  for k1 in range(n//m):
    for k2 in range(m):
      for i in range(m):
        p2[k1 + (n//m)*k2] += h[i][k1] * (ws[(w_power*i*k1 + w_power*(n//m)*i*k2)&maskMod])
  return p2

next_power_of_2 = lambda x: 1<<int(log2(x))

def shallow_ntt( primitive_root, values ):
    assert(len(values) & (len(values) - 1) == 0), "cannot compute ntt of non-power-of-two sequence"
    if len(values) <= 1:
        return values

    ntt_max_depth = 2
    # radix
    m = next_power_of_2(floor((len(values)/(ntt_max_depth - 1))**(1/((ntt_max_depth - 1) + 1))))
    ws = [FieldElement(1)] + [0]*(len(values) - 1)
    for i in range(1, len(values)):
        ws[i] = ws[i - 1]*primitive_root
    return _ntt_rk2(values, ws, 1, 1, len(values) - 1, m, ntt_max_depth)

def shallow_expand( polynomial, generator, order ):
    values = shallow_ntt(generator, polynomial.coefficients + [FieldElement(0)] * (order - len(polynomial.coefficients)))
    return values
 
def intt( primitive_root, values ):
    assert(len(values) & (len(values) - 1) == 0), "cannot compute intt of non-power-of-two sequence"

    if len(values) == 1:
        return values

    ninv = FieldElement(len(values))**-1

    transformed_values = shallow_ntt(primitive_root**-1, values)
    return [ninv*tv for tv in transformed_values]

def fast_multiply( lhs, rhs, primitive_root, root_order ):
    assert(primitive_root**root_order == 1), "supplied root does not have supplied order"
    assert(primitive_root**(root_order//2) != 1), "supplied root is not primitive root of supplied order"

    if lhs.is_zero() or rhs.is_zero():
        return Polynomial([])

    root = primitive_root
    order = root_order
    degree = lhs.degree() + rhs.degree()

    if degree < 8:
        return lhs * rhs

    while degree < order // 2:
        root = root**2
        order = order // 2

    lhs_coefficients = lhs.coefficients[:(lhs.degree()+1)]
    while len(lhs_coefficients) < order:
        lhs_coefficients += [0]
    rhs_coefficients = rhs.coefficients[:(rhs.degree()+1)]
    while len(rhs_coefficients) < order:
        rhs_coefficients += [0]

    lhs_codeword = ntt(root, lhs_coefficients)
    rhs_codeword = ntt(root, rhs_coefficients)

    hadamard_product = [l * r for (l, r) in zip(lhs_codeword, rhs_codeword)]

    product_coefficients = intt(root, hadamard_product)
    return Polynomial(product_coefficients[0:(degree+1)])

def fast_zerofier( domain, primitive_root, root_order ):
    assert(primitive_root**root_order == 1), "supplied root does not have supplied order"
    assert(primitive_root**(root_order//2) != 1), "supplied root is not primitive root of supplied order"

    if len(domain) == 0:
        return Polynomial([])

    if len(domain) == 1:
        return Polynomial([-domain[0], 1])

    half = len(domain) // 2

    left = fast_zerofier(domain[:half], primitive_root, root_order)
    right = fast_zerofier(domain[half:], primitive_root, root_order)
    return fast_multiply(left, right, primitive_root, root_order)

def fast_evaluate( polynomial, domain, primitive_root, root_order ):
    assert(primitive_root**root_order == 1), "supplied root does not have supplied order"
    assert(primitive_root**(root_order//2) != 1), "supplied root is not primitive root of supplied order"

    if len(domain) == 0:
        return []

    if len(domain) == 1:
        return [polynomial.evaluate(domain[0])]

    half = len(domain) // 2

    left_zerofier = fast_zerofier(domain[:half], primitive_root, root_order)
    right_zerofier = fast_zerofier(domain[half:], primitive_root, root_order)

    left = fast_evaluate(polynomial % left_zerofier, domain[:half], primitive_root, root_order)
    right = fast_evaluate(polynomial % right_zerofier, domain[half:], primitive_root, root_order)

    return left + right

def fast_interpolate( domain, values, primitive_root, root_order ):
    assert(primitive_root**root_order == 1), "supplied root does not have supplied order"
    assert(primitive_root**(root_order//2) != 1), "supplied root is not primitive root of supplied order"
    assert(len(domain) == len(values)), "cannot interpolate over domain of different length than values list"

    if len(domain) == 0:
        return Polynomial([])

    if len(domain) == 1:
        return Polynomial([values[0]])

    half = len(domain) // 2

    left_zerofier = fast_zerofier(domain[:half], primitive_root, root_order)
    right_zerofier = fast_zerofier(domain[half:], primitive_root, root_order)

    left_offset = fast_evaluate(right_zerofier, domain[:half], primitive_root, root_order)
    right_offset = fast_evaluate(left_zerofier, domain[half:], primitive_root, root_order)

    if not all(not (v==0) for v in left_offset):
        print("left_offset:", " ".join(str(v) for v in left_offset))

    left_targets = [n / d for (n,d) in zip(values[:half], left_offset)]
    right_targets = [n / d for (n,d) in zip(values[half:], right_offset)]

    left_interpolant = fast_interpolate(domain[:half], left_targets, primitive_root, root_order)
    right_interpolant = fast_interpolate(domain[half:], right_targets, primitive_root, root_order)

    return left_interpolant * right_zerofier + right_interpolant * left_zerofier

def fast_coset_evaluate( polynomial, offset, generator, order ):
    scaled_polynomial = polynomial.scale(offset)
    values = ntt(generator, scaled_polynomial.coefficients + [0] * (order - len(polynomial.coefficients)))
    return values

def fast_coset_divide( lhs, rhs, offset, primitive_root, root_order ): # clean division only!
    assert(primitive_root**root_order == 1), "supplied root does not have supplied order"
    assert(primitive_root**(root_order//2) != 1), "supplied root is not primitive root of supplied order"
    assert(not rhs.is_zero()), "cannot divide by zero polynomial"

    if lhs.is_zero():
        return Polynomial([])

    assert(rhs.degree() <= lhs.degree()), "cannot divide by polynomial of larger degree"

    root = primitive_root
    order = root_order
    degree = max(lhs.degree(),rhs.degree())

    if degree < 8:
        return lhs / rhs

    while degree < order // 2:
        root = root**2
        order = order // 2

    scaled_lhs = lhs.scale(offset)
    scaled_rhs = rhs.scale(offset)
    
    lhs_coefficients = scaled_lhs.coefficients[:(lhs.degree()+1)]
    while len(lhs_coefficients) < order:
        lhs_coefficients += [0]
    rhs_coefficients = scaled_rhs.coefficients[:(rhs.degree()+1)]
    while len(rhs_coefficients) < order:
        rhs_coefficients += [0]

    lhs_codeword = ntt(root, lhs_coefficients)
    rhs_codeword = ntt(root, rhs_coefficients)

    quotient_codeword = [l / r for (l, r) in zip(lhs_codeword, rhs_codeword)]
    scaled_quotient_coefficients = intt(root, quotient_codeword)
    scaled_quotient = Polynomial(scaled_quotient_coefficients[:(lhs.degree() - rhs.degree() + 1)])

    return scaled_quotient.scale(offset**-1)

