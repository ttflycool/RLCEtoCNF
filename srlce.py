import reedsolo as rs
import numpy as np 
import numpy.matlib 
import scipy.linalg


m=3
n=7  # length of total message+ecc
k=4
w=2
t=1
nsym =n-k  # length of ecc


prim = rs.find_prime_polys(c_exp=m, fast_primes=True, single=True)
rs.init_tables(prim=prim, generator=2, c_exp=m)
g0=rs.rs_generator_poly(nsym)
def g(g0,n):#生成多项式读取扩张
    g1=np.frombuffer(g0,dtype=np.uint8)
    g2=np.pad(g1,(0,n-len(g1)),'constant')
    g2=np.array(g2)
    return g2

def GRS(g):#rs码生成矩阵Обобщенные коды Рида-Соломона
    g=g(g0,n)
    GRS=g
    for i in range(1,k):
        g=np.roll(g,1)
        GRS=np.concatenate((GRS,g))
    GRS=np.array(GRS.reshape(k,n))
    return GRS

def V(n):
    a=np.random.randint(1,n-2,n)
    v=np.diag(a)
    v=np.array(v)
    return v

def GFmatrixMul(A, B):
    res=[[0] * len(B[0]) for i in range(len(A))]
    for i in range(len(A)):
            for j in range(len(B[0])):
                for m in range(len(B)):
                     x=rs.gf_mul(A[i][m],B[m][j])
                     res[i][j]=rs.gf_add(res[i][j],x)
    return res


def Gs(g,n):#广义里德所罗门码生成矩阵 обобщенным кодом Рида-Соломона.
    A=GRS(g)
    B=V(n)
    Gs=GFmatrixMul(A,B)
    Gs=np.array(Gs)
    return Gs

def R(w,k):#вектор-столбец
    x=np.random.randint(1,n,size=(k,w))
    x=np.array(x)
    return x

def G1(g,n,w,k): #插入列 G1=(GS,R), Случайно вставить вектор-столбец R=(r0,r1,...,rw−1)
    m=Gs(g,n)
    for j in range(0,w):
        x=np.random.randint(1,n)
        m=np.insert(m,x,R(w,k)[:,j],axis=1)
    return m

def getMatrixA(k,w):#稀疏矩阵A A-(n+w)×(n+w) невырожденная матрица:
    rE=np.random.randint(1,k,size=(1,4*w+19))
    i=0
    j=0
    IA=np.matlib.eye(n-w,dtype=int)
    while i<w:
        det=0
        while det==0:
            det=rE[0,j]*rE[0,j+3]-rE[0,j+1]*rE[0,j+2]
            while det==0:
                j=j+1
            else:
                    A=np.zeros(shape=(2,2),dtype=int)
                    A[0][0]=rE[0,j]
                    A[0][1]=rE[0,j+1]
                    A[1][0]=rE[0,j+2]
                    A[1][1]=rE[0,j+3]
                    j=j+4
                    i=i+1
                    IA=scipy.linalg.block_diag(IA,A)        
    return IA

def MatrixG2(g,n,k,w):#G2 	G2-k×(n+w) матрица ,G2=G1A.
    a=G1(g,n,w,k)
    b=getMatrixA(k,w)
    G2=GFmatrixMul(a, b)
    G2=np.array(G2)
    return G2


def permutation(n,w):# P-матрица перестановок (n+w)×(n+w)
    x=np.matlib.eye(n+w,dtype=int)
    x=np.random.permutation(x)
    return x

def G3(g,n,k,w):#G2*P
    g3=GFmatrixMul(MatrixG2(g,n,k,w),permutation(n,w))
    g3=np.array(g3)
    return g3

def Snonsingula(n,k):#生成非奇异矩阵s S-произвольно выбранная невырожденная матрица размера k×k
    while True:
        x=np.random.randint(0,n-1,(k,k)) 
        B=np.linalg.det(x) 
        if B!= 0:
            break
    return x

def Gpub(g,n,k,w):
    #SAGP=SG3
    Gpub=GFmatrixMul(Snonsingula(n,k), G3(g,n,k,w))
    Gpub=np.array(Gpub)
    return Gpub



A=Gpub(g,n,k,w)
print("GPUBisok!")
