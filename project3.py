import numpy as np
import numpy.linalg as la

"""Rebecca Selvaggini"""

def gs_step_1d(uh, fh): 
    N = len(fh) -1
    h = 4/N
    uh1 = [fh[0]] 
    for i in range(1, N) :
        uh1.append(((h**2)*fh[i]  + uh[i+1] + uh1[i-1])/(2+h**2))
    uh1.append(fh[N])
    uh1[:] = np.array(uh1)
    res = la.norm(uh1-uh, np.inf)
    uh[:] = uh1
    return res


def Auf(uh, fh): #compute fh - A*uh
    N = len(uh)-1
    h = 4/N
    v = np.zeros_like(uh)
    v[0] = fh[0]
    for i in range(1,len(uh)-1):
        v[i] = fh[i] + (uh[i-1] - (2+h**2)*uh[i]  +uh[i+1])/(h**2)
    v[len(fh)-1] = fh[len(fh)-1]
    return v


def I_2h(v2h): #prolongation operator in 1d
    v = np.zeros(2*len(v2h)-1)
    for i in range((len(v)-1) // 2):
        v[2*i] = v2h[i]
        v[2*i+1] = (v2h[i] + v2h[i+1]) / 2
    v[-1] = v2h[-1]
    return v

def I_h(vh): #full weighting restriction in 1d
    N = len(vh)-1
    l = N // 2 + 1
    v2h = np.zeros(l)
    for i in range(1,l-1):
      v2h[i] = (vh[2*i-1]+ 2*vh[2*i] + vh[2*i+1])/4
    return v2h 


def Injection(vh):
    N =len(vh)
    if N == 4:
        return np.array([vh[0], vh[2],vh[3]])
    v2h = np.zeros(N//2 +1)
    for i in range( N//2):
        v2h[i] = vh[2*i]
    v2h[-1] = vh[-1]
    return np.array(v2h)


def twogrid_step(uh,fh): #two grid correction scheme
    N = len(uh)-1
    h = 4/N
    gs_step_1d(uh, fh)
    rh = np.zeros_like(uh)
    for i in range(1,N):
        rh[i] = fh[i] + (uh[i-1] - (2 + h**2)*uh[i] + uh[i+1])/(h**2)
    r2h = I_h(rh)
    e2h = np.zeros_like(r2h)
    for i in range(5):
        gs_step_1d(e2h,r2h)
    eh = I_2h(e2h)
    uh[:] += eh
    res = gs_step_1d(uh, fh)
    return res


def v_cycle_step_1d(uh, fh, alpha1, alpha2) :
    N = len(uh)-1
    h = 4.0/N
    if N==2 :
        uh[:] = np.array([fh[0], ((h**2)*fh[1]+fh[0]+fh[2])/(2+h**2) ,fh[2]])
        return 0.0
    for i in range(alpha1):
        gs_step_1d(uh, fh)
    tmp = Auf(uh,fh)
    f2h = I_h(tmp).copy()
    u2h = np.zeros_like(f2h)
    v_cycle_step_1d(u2h,f2h,alpha1,alpha2) 
    uh[:] += I_2h(u2h)    
    for i in range(alpha2):
        res = gs_step_1d(uh, fh)
    return res



def full_mg_1d(uh, fh, alpha1, alpha2, nu):
    N = len(fh)-1
    h = 4.0/N
    if N==2 :
        uh[:] = np.array([fh[0], ((h**2)*fh[1]+fh[0]+fh[2])/(2+h**2) ,fh[2]])
        return 0.0
    else:
        f2h = Injection(fh)
        u2h = np.zeros_like(f2h)
        full_mg_1d(u2h, f2h, alpha1, alpha2, nu) 
        uh[:] = I_2h(u2h)
        for i in range(nu):
            res = v_cycle_step_1d(uh,fh,alpha1,alpha2)
        return res
        
    
#2d multigrid

def gs_step_2d(uh, fh):
    N = len(fh) -1
    h = 4/N
    uh1 = fh.copy()
    res = 0
    for i in range(1, N) :
        for j in range(1,N):
            uh1[i,j] = (h**2*fh[i,j] + uh1[i-1,j] + uh[i+1,j] + uh1[i,j-1] + uh[i,j+1]) / (4+h**2)
    res = la.norm((uh1-uh).flat, np.inf)
    uh[:] = uh1
    return res


def Prolongation2d(v2h):
    N = len(v2h)-1
    vh = np.array([np.zeros(2*N +1) for i in range(2*N+1)])
    for i in range(N+1):
        for j in range(N+1):
            vh[2*i, 2*j] = v2h[i,j]
    for i in range(N):
        for j in range(N+1):
            vh[2*i+1, 2*j] =(v2h[i,j]+v2h[i+1,j])/2
    for i in range(N+1):
        for j in range(N):
            vh[2*i, 2*j+1] = (v2h[i,j]+v2h[i,j+1])/2
    for i in range(N):
        for j in range(N):
            vh[2*i+1, 2*j+1] = (v2h[i,j]+v2h[i+1,j]+v2h[i,j+1]+v2h[i+1,j+1])/4
    return vh

def g(x, y): return (x**2 + y**2) / 10

def set_boundary(vh):
    N = len(vh)-1
    b = g(*np.meshgrid(np.linspace(-2,2,N+1),np.linspace(-2,2,N+1)))
    vh[0][:] = b[0][:]
    vh[N][:] = b[N][:]
    vh[:][0] = b[:][0]
    vh[:][N] = b[:][N]
    return vh
    

def Restriction2d(vh):
    N = len(vh) -1
    v2h = np.array([np.zeros(N//2+1) for i in range(N//2+1)])
    for i in range(1,N//2):
        for j in range(1,N//2):
            v2h[i,j] = (vh[2*i-1,2*j-1]+vh[2*i-1,2*j+1]+vh[2*i+1,2*j-1] + vh[2*i+1,2*j+1] + 2*(vh[2*i,2*j-1]+vh[2*i,2*j+1]+vh[2*i-1,2*j]+vh[2*i+1,2*j]) + 4*vh[2*i,2*j])/16
    return v2h

def Injection2d(vh):
    N = len(vh)-1
    k = N//2 +1
    v2h = np.array([np.zeros(k) for i in range(k)])
    for i in range(k):
        for j in range(k):
            v2h[i,j] = vh[2*i, 2*j]
    return v2h


def v_cycle_step_2d(uh, fh, alpha1, alpha2):
    N = len(uh)-1
    h = 4.0/N
    if N==2 :
        tmp = fh.copy()
        tmp[1,1] = (fh[0,1] + fh[2,1] +4*fh[1,1] + fh[1,0] + fh[1,2])/8
        uh[:] = tmp
        return 0.0
    for i in range(alpha1):
        gs_step_2d(uh, fh)
    tmp = np.zeros_like(uh)
    for i in range(1,N):
        for j in range(1,N):
            tmp[i][j] = fh[i][j] + (-(4 + h**2)*uh[i][j] + uh[i-1][j] + uh[i+1][j]+ uh[i][j+1]+ uh[i][j-1])/(h**2)
    f2h = Restriction2d(tmp).copy()
    u2h = np.zeros_like(f2h)
    v_cycle_step_2d(u2h,f2h,alpha1,alpha2)
    uh[:] += Prolongation2d(u2h)  
    for i in range(alpha2):
        res = gs_step_2d(uh, fh)
    return res    



def full_mg_2d(uh, fh, alpha1, alpha2, nu):
    N = len(fh)-1
    h = 4.0/N
    if N==2 :
        tmp = fh.copy()
        tmp[1,1] = (fh[0,1] + fh[2,1] +4*fh[1,1] + fh[1,0] + fh[1,2])/8
        uh[:] = tmp
        return 0.0
    else:
        f2h = Injection2d(fh)
        u2h = np.zeros_like(f2h)
        full_mg_2d(u2h, f2h, alpha1, alpha2, nu) 
        uh[:] = Prolongation2d(u2h)
        uh[:] = set_boundary(uh)
        for i in range(nu):
            res = v_cycle_step_2d(uh,fh,alpha1,alpha2)
        return res