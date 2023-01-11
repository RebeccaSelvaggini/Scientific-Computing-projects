def complete_shuffle(k):
    " " " Rebecca Selvaggini " " "
    L = [*range(1, 2*k +1)]
    s = [0]*(2*k)
    f = [0]*(2*k)
    for m in range(1, k+1) :
        tmp = []
        for i in range(0,m):
            tmp.append(L[i+m])
            tmp.append(L[i])
        L = tmp + L[2*m:]
        f[L[0]-1] = f[L[0]-1] + 1
        if s[L[0]-1] == 0 :
            s[L[0]-1] = m
    return L[0:3], s[L[0]-1], f[L[0]-1]


