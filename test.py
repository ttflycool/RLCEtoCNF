import numpy as np 


def Snonsingula(n,k):#生成非奇异矩阵s S-произвольно выбранная невырожденная матрица размера k×k
    while True:
        x=np.random.randint(0,n-1,(k,k)) 
        B=np.linalg.det(x) 
        if B!= 0:
            break
    return x

x=np.random.randint(1,5,size=(3,7))
aa=Snonsingula(7,4)
print(aa)
print(x)