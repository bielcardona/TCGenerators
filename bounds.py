from phylonetwork import memoize_function


@memoize_function
def p0(n, h, k):
    if h > n-1:
        return 0
    prod = 1
    for i in range(2*n-2*h-2*k, 2*n-2*h, 2):
        prod *= i
    return prod


@memoize_function
def p1(n, h, k):
    return k * p0(n, h, k-1)


@memoize_function
def p(n, h, k):
    return p0(n, h, k)+p1(n, h, k)


@memoize_function
def FH1(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1) * (2*n+h-r)


@memoize_function
def FH2(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1) * (2*n+h-r) * (2*n+h-r-1) // 2


@memoize_function
def FT1(n, h, rm1):
    r = rm1 + 1
    if r <= 0:
        return 0
    return p(n, h, r-1) * (2*n+h-r)


@memoize_function
def FT2(n, h, rm1):
    r = rm1 + 1
    if r <= 1:
        return 0
    return p(n, h, r-2) * (r-1) * 3*h


@memoize_function
def FT3A(n, h, rm1):
    r = rm1+1
    if r <= 2:
        return 0
    return (r-1)*(r-2)*p1(n, h, r-3)*(2*n-2*h-2*r+6)


@memoize_function
def FT3B(n, h, rm1):
    r = rm1+1
    if r <= 2:
        return 0
    return (r-1)*(r-2)*p0(n, h, r-3)*(2*n-2*h-2*r+4)


@memoize_function
def FT4(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1)*(r-1)


@memoize_function
def FT(n, h, rm1):
    return FT1(n, h, rm1)+FT2(n, h, rm1)+FT3A(n, h, rm1)+FT3B(n, h, rm1)+FT4(n, h, rm1)


@memoize_function
def FH(n, h, rm1):
    return FH1(n, h, rm1)+FH2(n, h, rm1)


@memoize_function
def BTCh(n, h):
    if n < 1:
        return 0
    if n == 1:
        if h == 0:
            return 1
        return 0
    return (sum([BTCh(n-1, h-(r-1))*FT(n-1, h-(r-1), r-1) for r in range(1, h+2)]) +
            sum([BTCh(n - 1, h - r) * FH(n - 1, h - r, r - 1) for r in range(1, h + 1)]))


@memoize_function
def BTC(n):
    if n == 1:
        return 1
    return sum([BTCh(n, h) for h in range(0, n)])


if __name__ == '__main__':
    for N in range(1, 21):
        print(f"|BTC({N})| <= {BTC(N)}")
