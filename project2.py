#Task 1
""" Rebecca Selvaggini """
def kmer_hist(S, k) :
    d = {}
    for j in range(len(S)-k+1):
        kmer = S[j:j+k]
        count = d.get(kmer,0)
        d[kmer] = count +1
    m = max(d.values())
    h = [0]*(m+1)
    mfkmers = []
    for kmer in d.keys() :
        i = d[kmer]
        h[i] = h[i]+1
        if i == m :
            mfkmers.append(kmer)
    return h, mfkmers


#Task 3

mod, base = 1001, 5

def extend_by_one(h, c):
    return (base * h + ord(c)) % mod

def remove_left(h, c, bstar):
    return (h - (ord(c) * bstar) % mod) % mod

def hash_value(s):
    h = 0
    for c in s: h = extend_by_one(h, c)    
    return h

#Version1

# def kmer_search(S, L):
#     n, k, m = len(S), len(L[1]), len(L)
#     freq = [0]*m
#     pos = [-1]*m
#     hp = [hash_value(L[i]) for i in range(m)]
#     hs = hash_value(S[:k])
#     bstar = pow(base, k-1, mod)
#     for i in range(n-k+1):
#         for j in range(m):
#             if hp[j] == hs and S[i:i+k] == L[j]: 
#                 freq[j] += 1
#                 if pos[j] == -1 : pos[j] = i
#         if i < n-k:
#             hs = remove_left(hs, S[i], bstar)
#             hs = extend_by_one(hs, S[i+k])
#     mf = max(freq)
#     i = freq.index(mf)
#     if freq.count(mf) == 1 : 
#         return pos[i], mf
#     else :
#         if mf == 0 : return None, 0   
#         while max(freq) == mf: 
#             i = freq.index(mf)
#             freq.remove(mf)
#             j = freq.index(mf)
#             p = min(pos[i], pos[j])
#         return p, mf

#Version2

def kmer_search(S, L):
    n, k, m = len(S), len(L[1]), len(L)
    freq = 0
    pos = -1
    hp = [hash_value(L[i]) for i in range(m)]
    hs = [hash_value(S[:k])]
    bstar = pow(base, k-1, mod)
    for i in range(0, n-k):
        h = remove_left(hs[i], S[i], bstar)
        h = extend_by_one(h, S[i+k])
        hs.append(h)
    for j in range(m):
        f = 0
        p = -1
        for i in range(n-k+1):   
            if hp[j] == hs[i] and S[i:i+k] == L[j]:
                f += 1
                if p == -1 : p = i
        if f > freq : 
            pos = p
            freq = f
        if f == freq and pos > p:
            pos = p
            freq = f
    if freq == 0 :
        return None, 0
    else:
        return pos, freq
