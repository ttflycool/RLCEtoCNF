import numpy as np
import reedsolo as rs
from H_G import b, h, rc
from srlce import A, k, m, n, prim, t, w


def output(temp):
    with open('out.cnf', 'a') as f:
        f.write('\n' + str(temp) + " 0")

def toxor(D, b):  # D是array，b是ture或flase即0，1
    l = len(D)
    if l == 1:
        if b == 1:
            output(D[0])
        else:
            output(-D[0])
    if l == 2:
        if b == 1:
            cause0 = str(D[0]) + " " + str(D[1])
            cause1 = str(-D[0]) + " " + str(-D[1])
            output(cause0)
            output(cause1)
        else:
            cause0 = str(D[0]) +" " + str(-D[1])
            cause1 = str(-D[0]) + " " + str(D[1])
            output(cause0)
            output(cause1)
    if l == 3:
        if b == 1:
            cause0 = str(-D[0]) + " " + str(-D[1]) + " " + str(D[2])
            cause1 = str(-D[0]) + " " + str(D[1]) + " " + str(-D[2])
            cause2 = str(D[0]) + " " + str(-D[1]) + " " + str(-D[2])
            cause3 = str(D[0]) + " " + str(D[1]) + " " + str(D[2])
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)
        else:
            cause0 = str(D[0]) + " " + str(D[1]) + " " + str(D[2])
            cause1 = str(-D[0]) + " " + str(D[1]) + " " + str(D[2])
            cause2 = str(D[0]) + " " + str(-D[1]) + " " + str(D[2])
            cause3 = str(D[0]) + " " + str(D[1]) + " " + str(-D[2])
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)
    if l == 4:
        if b == 1:
            cause0 = str(D[0]) + " " + str(D[1]) + " " + str(D[2]) + " " + str(D[3])
            cause1 = str(-D[0]) + " " + str(-D[1]) + " " + str(-D[2]) + " " + str(-D[3])
            cause2 = str(-D[0]) + " " + str(-D[1]) + " " + str(D[2]) + " " + str(D[3])
            cause3 = str(-D[0]) + " " + str(D[1]) + " " + str(-D[2]) + " " + str(D[3])
            cause4 = str(-D[0]) + " " + str(D[1]) + " " + str(D[2]) + " " + str(-D[3])
            cause5 = str(D[0]) + " " + str(-D[1]) + " " + str(-D[2]) + " " + str(D[3])
            cause6 = str(D[0]) + " " + str(-D[1]) + " " + str(D[2]) + " " + str(-D[3])
            cause7 = str(D[0]) + " " + str(D[1]) + " " + str(-D[2]) + " " + str(-D[3])
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)
            output(cause4)
            output(cause5)
            output(cause6)
            output(cause7)
        else:
            cause0 = str(-D[0]) + " " + str(-D[1]) + " " + str(-D[2]) + " " + str(D[3])
            cause1 = str(-D[0]) + " " + str(-D[1]) + " " + str(D[2]) + " " + str(-D[3])
            cause2 = str(-D[0]) + " " + str(D[1]) + " " + str(-D[2]) + " " + str(-D[3])
            cause3 = str(-D[0]) + " " + str(D[1]) + " " + str(D[2]) + " " + str(D[3])
            cause4 = str(D[0]) + " " + str(-D[1]) + " " + str(-D[2]) + " " + str(-D[3])
            cause5 = str(D[0]) + " " + str(-D[1]) + " " + str(D[2]) + " " + str(D[3])
            cause6 = str(D[0]) + " " + str(D[1]) + " " + str(-D[2]) + " " + str(D[3])
            cause7 = str(D[0]) + " " + str(D[1]) + " " + str(D[2]) + " " + str(-D[3])
            output(cause0)
            output(cause1)
            output(cause2)
            output(cause3)
            output(cause4)
            output(cause5)
            output(cause6)
            output(cause7)


def gptog2(num,m):#把一个数转化为m位二进制形式
    l=[]
    for i in range (1,m+1):
        num, remainder = divmod(num, 2)
        l.append(remainder)
    l=np.array(l)
    l=l.reshape(l.shape[0],1)
    return l



def binb():#hcT=b ，n+w-k 
    x=gptog2(b[0][0],m)
    for i in range (1,n+w-k):
        b1=gptog2(b[i][0],m)
        x=np.concatenate((x,b1),axis=0)
    x=x.flatten()
    return x


def binA():#Ay.T=b1 
    for i in range(0,n+w-k):#读取h矩阵的行
        hline=h[i]
        lenth=len(hline)
        for j in range(0,lenth):#读取每行中的一个元素
            for p in range(0,m):#每个元素扩充为m列
                b=rs.gf_mul(hline[j],2**p)
                a1=gptog2(b,m)
                if p==0 and j==0:
                    a=a1
                else:
                    a=np.concatenate((a,a1),axis=1)
        if i==0:
            Ablock=a 
        else:
            Ablock=np.concatenate((Ablock,a),axis=0)
    print("ABLOCK\n",Ablock)
    return  Ablock


def numblock():#hA矩阵转化为序号，代表第多少个未知数1=e1,2=e12,...
    HA=binA()
    numrow=HA.shape[0]
    numcol=HA.shape[1]
    for i in range(0,numrow):
        for j in range(0,numcol):
            if HA[i,j]==1:
                HA[i,j]=j+1
    return HA

def extocnf(yy):#x=e_11 or e_12...e_1m,y为空的新变量
    x=1
    for i in range (0,n+w):
        list0=[]
        for j in range (0,m):
            list0.append(x)
            x=x+1
        list1=list(list0)
        list1.append(yy)
        list1[0]=-list1[0]
        list1=[str(z) for z in list1]
        list1=" ".join(list1)
        output(list1)
        list0.append(-yy)
        list0=[str(z) for z in list0]
        list0=" ".join(list0)
        output(list0)
        yy=yy+1
    return yy  
    

def blockxortocnf():#生成cnf文件
    a=numblock()
    b=binb()
    b=b.flatten()#将b变为一维数组b[0]......以此来对应h[0]...
    numrow=a.shape[0]
    y=np.arange(m*(n+w)+1,10000000,1)
    county=0
    for i in range(0,numrow):#一共有m*(n+w)行的式子
        temp=a[i]
        temp=temp[temp != 0]#去除所有0的h的一行
        l=len(temp)#长度
        if l==0:
            continue
        if 0<l<5:
            b1=b[i]#b=[b11,b12,b13,....,b21,b22,....]
            toxor(temp,b1)
        if 4<l<7:
            block1=temp[0:3]#第一列的123个e11,e12,e13,y1 []左闭右开
            y1=y[county]
            block1=np.append(block1,y1)
            b1=b[i]
            toxor(block1,b1)
            block3=temp[3:l]#第二个xor块
            block3=np.append(block3,y1)
            toxor(block3,0)
            county=county+1
        if l>6:
            num, remainder = divmod(l-3, 2)#num表示中间有几个2
            block1=temp[0:3]#第一列的123个e11,e12,e13,y1
            y1=y[county]
            block1=np.append(block1,y1)
            b1=b[i]
            toxor(block1,b1)
            for j in range (1,num):
                s=1+j*2
                block2=temp[s:s+2]#j j+1,  e4 e5 y1 y2
                y2=y[county:county+2]
                county=county+1
                block2=np.append(block2,y2)
                toxor(block2,0)

            block3=temp[-(remainder+2):]
            y3=y[county]
            block3=np.append(block3,y3)
            toxor(block3,0)
            county=county+1
    yy=y[county]+1#空的新变量
    lastnum=extocnf(yy)
    return yy


y=blockxortocnf()
print("75%")
