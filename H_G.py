import fieldmath
import numpy as np
import srlce  
from srlce import A,prim,k,n,w,t
from randomerror import rerror


A=A.tolist()
modprim=int(prim)

# A=[[0,0,1, 2, 1, 1, 1, 7, 0],[0, 0, 4, 5, 3, 5, 0, 0, 6],
# [5, 6, 0, 5, 2, 2, 0, 5, 1],[3,0, 0, 2, 4, 0, 0, 2, 5]]

def GPUBI(A):#生成矩阵的标准化
    f = fieldmath.BinaryField(modprim)#prim
    numRows = len(A)#行
    numCols = len(A[0])#列
    mat=fieldmath.Matrix(numRows, numCols, f)
    for x in range (numRows):
        for y in range (numCols):
            mat.set(x, y, A[x][y])
    mat.reduced_row_echelon_form()
    GPUB=np.zeros((k,n+w),dtype = np.int)
    xx=np.array(mat.values)
    return xx  

def GE(A,k):
    GE=GPUBI(A)[...,k:]
    return GE

def H(A,n,w,k):#校验矩阵
    I=np.eye(n+w-k,dtype =int)
    GE1=GE(A,k).T
    H=np.concatenate((GE1,I),axis=1)
    H=np.array(H)
    return H


def c(msg,A,n,w,t,rerror):
    c=srlce.GFmatrixMul(msg, A)
    c=np.array(c)
    x=c
    e=rerror
    for i in range(len(c[0])):
            x[0,i]=c[0,i]^e[i]
    return x

def hc(A,n,w,k): #计算中要用到的b参数即H*CT
    ct=c(msg,A,n,w,t,rerror).T#c的转置
    x=srlce.GFmatrixMul(H(A,n,w,k), ct)
    x=np.array(x)
    return x

msg='2'*k
msg=np.array(np.fromiter(msg, dtype=int)).reshape(1,k)
G=GPUBI(A)#最后的高斯消元后的系统形态的公钥
h=H(A,n,w,k)#校验矩阵
b= hc(A,n,w,k)
rc=c(msg,A,n,w,t,rerror)
print("生成矩阵\n",G)
print('校验矩阵\n',h)
print('输入的消息\n',msg)
print('随机错误向量\n',rerror)
print("接收到的消息C\n",rc)
print("H*cT=b\n",b)





