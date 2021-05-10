import numpy as np
import math
from Tocnf import y
from srlce import k, m, n,t, w

def output(temp):
    with open('out.cnf', 'a') as f:
        f.write('\n' + str(temp) + " 0")


def ttog2(num,m):#把一个数转化为m位二进制形式
    l=[]
    for i in range (1,m+1):
        num, remainder = divmod(num, 2)
        l.append(remainder)
    l=np.array(l)
    return l
    
def mxyqtocnf(x,y,z,q):#q=x*y or  y*z or z*x
    cause0=str(x) + " " + str(-q) 
    cause1=str(y) + " " + str(-q)
    cause2=str(z) + " " + str(-q)
    cause3=str(-x) + " " + str(-y) + " " + str(-z) + " " + str(q)
    output(cause0)
    output(cause1)
    output(cause2)
    output(cause3)



def zzzqxortocnf (x,y,z,q):# z0 xor z1 xor z xor q xor 1=1
        cause0 = str(x) + " " + str(y) + " " + str(z) + " " + str(q)
        cause1 = str(-x) + " " + str(-y) + " " + str(-z) + " " + str(-q)
        cause2 = str(-x) + " " + str(-y) + " " + str(z) + " " + str(q)
        cause3 = str(-x) + " " + str(y) + " " + str(-z) + " " + str(q)
        cause4 = str(-x) + " " + str(y) + " " + str(z) + " " + str(-q)
        cause5 = str(x) + " " + str(-y) + " " + str(-z) + " " + str(q)
        cause6 = str(x) + " " + str(-y) + " " + str(z) + " " + str(-q)
        cause7 = str(x) + " " + str(y) + " " + str(-z) + " " + str(-q)
        output(cause0)
        output(cause1)
        output(cause2)
        output(cause3)
        output(cause4)
        output(cause5)
        output(cause6)
        output(cause7)

def xyztocnf(x,y,z):#末尾 x xor y xor z=1
            cause0 = str(-x) + " " + str(-y) + " " + str(z)
            cause1 = str(-x) + " " + str(y) + " " + str(-z)
            cause2 = str(x) + " " + str(-y) + " " + str(-z)
            cause3 = str(x) + " " + str(y) + " " + str(z)
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)

def xytoq(x,y,z):#提供的q   x and y =q
    cause0 =str(x) + " " + str(-z)
    cause1 =str(y) + " " + str(-z)
    cause2 = str(-x) + " " + str(-y) + " " + str(z)
    output(cause0)
    output(cause1)
    output(cause2)

def zqtocnf(z,q):#z=q
    cause0 = str(z) + " " + str(-q)
    cause1 = str(-z) + " " + str(q)
    output(cause0)
    output(cause1)

def txyqxor(x,y,q,t):
        if t== 1:
            cause0 = str(-x) + " " + str(-y) + " " + str(q)
            cause1 = str(-x) + " " + str(y) + " " + str(-q)
            cause2 = str(x) + " " + str(-y) + " " + str(-q)
            cause3 = str(x) + " " + str(y) + " " + str(q)
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)
        else:
            cause0 = str(x) + " " + str(y) + " " + str(q)
            cause1 = str(-x) + " " + str(y) + " " + str(q)
            cause2 = str(x) + " " + str(-y) + " " + str(q)
            cause3 = str(x) + " " + str(y) + " " + str(-q)
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)

def tqz(q,t):
    if t==1:
        cause1=str(q)
        output(cause1)
    else:
        cause2=str(-q)
        output(cause2)
def txy(x,y,t):
    if t == 1:
        cause0 = str(x) + " " + str(y)
        cause1 = str(-x) + " " + str(-y)
        output(cause0)
        output(cause1)
    else:
        cause0 = str(x) + " " + str(-y)
        cause1 = str(-x) + " " + str(y)
        output(cause0)
        output(cause1)

def treetocnf():
    l=m+1
    for s in range (0,l):#深度
        h=l-1-s#第几层，由下开始向上
        sumq=(2**h)*(l-h)#每层q变量数量
        sumz=(2**h)*(l+1-h)#每层z变量数量
        lengthz=s+2#代表z1的位数
        if s==0:#底层
            z0=np.arange(y+1,y+1+(2**l))#z0此处为x的表
            z1=np.arange(z0[-1]+1,z0[-1]+1+sumz)#sumz这里表示z1有多少个
            q=np.arange(z1[-1]+1,z1[-1]+1+sumq)#q的储存表
            i=0
            qi=0
            for j in range (0,sumz,2):
                xyztocnf(z0[i],z0[i+1],z1[j])
                xytoq(z0[i],z0[i+1],q[qi])#x0 xor x1 xor q0=1
                zqtocnf(z1[j+1],q[qi])
                i=i+2
                qi=qi+1

        if s>0 and s<(l-1) :       #其它层
            lastnum=q[-1]+1 
            z0=z1
            z1=np.arange(lastnum,lastnum+sumz)
            q=np.arange(z1[-1]+1,z1[-1]+1+sumq)
            i=0     #用来对z0计数
            qi=0    #对q计数 
            for j in range (0,sumz,lengthz):#z1中循环，即高位的z循环  
                for p in range (0,lengthz):#p此处代表位数
                    if p==0:
                        xyztocnf(z0[i],z0[i+lengthz-1],z1[j])
                        i=i+1

                    if p==1:
                        zzzqxortocnf(z0[i],z0[i+lengthz-1],q[qi],z1[j+1])
                        xytoq(z0[i-1],z0[i+lengthz-2],q[qi]) 
                        i=i+1   
                        qi=qi+1                    

                    if p>1 and p<(lengthz-1):#中间的z1
                        zzzqxortocnf(z0[i],z0[i+lengthz-1],q[qi],z1[j+p])
                        mxyqtocnf(z0[i-1],z0[i+lengthz-2],q[qi-1],q[qi])
                        i=i+1
                        qi=qi+1

                    if p==(lengthz-1): #最高次的z=q
                        zqtocnf(z1[j+lengthz-1],q[qi])
                        mxyqtocnf(z0[i-1],z0[i+lengthz-2],q[qi],q[qi-1])
                        i=i+lengthz-1
                        qi=qi+1
        if s==l-1:
            z0=z1
            z1=ttog2(t,lengthz)
            q=np.arange(q[-1]+1,q[-1]+1+sumq)
            print("the sum x=",q[-1])
            #输出z1[0]
            txy(z0[0],z0[lengthz-1],z1[0])
            #z1[1]
            txyqxor(z0[1],z0[lengthz],q[0],z1[1])
            xytoq(z0[0],z0[lengthz-1],q[0])
            #z1[lengthz-1]
            tqz(q[-1],z1[-1])
            mxyqtocnf(z0[lengthz-2],z0[-1],q[-2],q[-1])
            for i in range(2,lengthz-1):
               txyqxor(z0[i],z0[i+lengthz-1],q[i-1],z1[i])
               mxyqtocnf(z0[i-1],z0[i+lengthz-2],q[i-1],q[i-2])
        

treetocnf()
print("100% \nok!")

