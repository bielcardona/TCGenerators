def p0(n, h, k):
    if h > n-1:
        return 0
    prod = 1
    for i in range(2*n-2*h-2*k, 2*n-2*h, 2):
        prod *= i
    return prod


def p1(n, h, k):
    return k * p0(n, h, k-1)


def p(n, h, k):
    return p0(n, h, k)+p1(n, h, k)


def FH1(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1) * (2*n+h-r)


def FH2(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1) * (2*n+h-r) * (2*n+h-r-1) // 2


def FT1(n, h, rm1):
    r = rm1 + 1
    if r <= 0:
        return 0
    return p(n, h, r-1) * (2*n+h-r)


def FT2(n, h, rm1):
    r = rm1 + 1
    if r <= 1:
        return 0
    return p(n, h, r-2) * (r-1) * 3*h


def FT3A(n, h, rm1):
    r = rm1+1
    if r <= 2:
        return 0
    return (r-1)*(r-2)*p1(n, h, r-3)*(2*n-2*h-2*r+6)


def FT3B(n, h, rm1):
    r = rm1+1
    if r <= 2:
        return 0
    return (r-1)*(r-2)*p0(n, h, r-3)*(2*n-2*h-2*r+4)


def FT4(n, h, rm1):
    r = rm1+1
    return p(n, h, r-1)*(r-1)


def FT(n, h, rm1):
    return FT1(n, h, rm1)+FT2(n, h, rm1)+FT3A(n, h, rm1)+FT3B(n, h, rm1)+FT4(n, h, rm1)


def FH(n, h, rm1):
    return FH1(n, h, rm1)+FH2(n, h, rm1)


def BTCh(n, h):
    if n < 1:
        return 0
    if n == 1:
        if h == 0:
            return 1
        return 0
    return (sum([BTCh(n-1, h-(r-1))*FT(n-1, h-(r-1), r-1) for r in range(1, h+2)]) +
            sum([BTCh(n - 1, h - r) * FH(n - 1, h - r, r - 1) for r in range(1, h + 1)]))


def BTC(n):
    if n == 1:
        return 1
    return sum([BTCh(n, h) for h in range(0, n)])


if __name__ == '__main__':
    for N in range(1, 11):
        print(f"BTC({N}) <= {BTC(N)}")
