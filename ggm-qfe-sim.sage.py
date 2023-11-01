import sympy as sy
import numpy as np
from sage.all import *
import copy
import time
import re

#设置一个结构体，多项式及其系数
class polys:
    p = ''
    coeff = ''
    def __init__(self,p,coeff):
        self.p = p
        self.coeff = coeff

class monomials:
    monomial = ''
    coeff = ''
    h = ''
    def __init__(self,monomial,coeff,h):
        self.monomial = monomial
        self.coeff = coeff
        self.h = h


#读取文件并分析文件
def read(name):
    G1_poly = []
    G2_poly = []
    list1 = [[],[],[],[],[],[],[],[]]
    count = 0
    file_object = open(name,'r')
    while True:
        line = file_object.readline()
        if line:
            for j in line.split(';'):
                b = j.find('in')
                if b >= 0:
                    tmp = j.find(':')
                    if tmp >=0:
                        list1[count] = list1[count] + j[tmp+1:b].replace(' ','').split(',')
                    else:
                        list1[count] = list1[count] + j[0:b].replace(' ','').split(',')
                else:
                    list1[count] = list1[count]+ ['no']
            count = count + 1
        else:
            break
    file_object.close()

    x = list1[0]
    y = list1[1]
    q = list1[2]
    variable = list1[3]
    public = list1[4]
    enc = list1[5]
    keygen = list1[6]
    offset = list1[7]

    tmp = ''
    tmp1 = ''
    param = []
    x_var = []
    y_var = []
    q_var = []
    offset_poly = []
    for i in range(len(x)):
        param.append(x[i])
        tmp = tmp + x[i]
        tmp = tmp + ','
    for j in range(len(y)):
        param.append(y[j])
        tmp = tmp + y[j]
        tmp = tmp + ','
    for k in range(len(q)-1):
        param.append(q[k])
        tmp = tmp + q[k]
        tmp = tmp + ','
    tmp = tmp + list1[2][len(q)-1]
    param.append(q[len(q)-1])
    for i in range(len(variable)-1):
        tmp1 =tmp1 + variable[i]
        tmp1 =tmp1 + ','
    tmp1 = tmp1 + list1[3][len(variable)-1]
    #生成多项式变量
    R_QQ = PolynomialRing(QQ,len(x)+len(y)+len(q),tmp)
    T_R = PolynomialRing(R_QQ,len(variable),tmp1)
    RT = R_QQ.gens()
    TT = T_R.gens()
    for i in range(len(param)):
        globals()[param[i]] = RT[i] 
        if i < len(x):
            x_var.append(globals()[param[i]])
        if i >= len(x) and i < (len(x)+len(y)):
            y_var.append(globals()[param[i]])
        if i >= (len(x)+len(y)):
            q_var.append(globals()[param[i]])

    for i in range(len(variable)):
        globals()[variable[i]] = TT[i]
    
    for i in public:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G2_poly.append(eval(m))
    
    for i in enc:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G2_poly.append(eval(m))

    for i in keygen:
        if i[-1] == '1':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G1_poly.append(eval(m))
        if i[-1] == '2':
            m = re.findall(r'\[(.*?)\]',i)[0]
            G2_poly.append(eval(m))

    if offset[0] != 'no':
        for i in offset:
            offset_poly.append(eval(i))
    return G1_poly,G2_poly,offset_poly

#将G1中的多项式组合和G2中的多项式组合进行pairing.
def parametric_completion(G1_poly,G2_poly,offset_poly):
    GT_poly = []
    G1_poly.append(G1_poly[0] * 0 + 1)
    G2_poly.append(G2_poly[0] * 0 + 1)
    for i in G1_poly:
        for j in G2_poly:
            poly = i * j
            GT_poly.append(poly)
    if offset_poly != []:
        GT_poly.remove(1)
        for i in offset_poly:
            GT_poly.append(i)
    return GT_poly

#提取GT中的变量组合和相关系数,并返回。
def merge(GT_poly):
    monomial = []
    coeff = []
    dict_merge = {}#用来合并同一变量组合的多项式的
    dict_count = {}#用来累计次数的
    dict_coeff = {}#用来累计系数的
    for i in GT_poly:
        monomial.append(i.monomials())
        coeff.append(i.coefficients())
#处理平方项问题，如果发现有大于2的项，则在后面添加一个新的只包含它的多项式，这样的话就不会被删掉了
    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            str_var = str(monomial[i][j])
            for k in str_var.split('*'):
                if '_' in k and '^' in k:
                    monomial.append([monomial[i][j]])
                    coeff.append([coeff[i][j]])
#设置一个字典，首先先初始化
    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            dict_merge[monomial[i][j]] = 0
            dict_count[monomial[i][j]] = 0
            dict_coeff[monomial[i][j]] = []
    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            dict_count[monomial[i][j]] = dict_count[monomial[i][j]] + 1

    set_monomial = []
    for i in range(len(monomial)):
        tmp = []
        for j in range(len(monomial[i])):
            tmp.append(monomials(monomial[i][j],var('h'+str(i))*coeff[i][j],var('h'+str(i))))
        set_monomial.append(tmp)
    
    #两轮后,留下可以配对的变量组合
    round = 0
    while(round < 2):
        for i in range(len(monomial)):
            flag1 = 0
            for j in range(len(monomial[i])):
                if dict_count[monomial[i][j]] == 1 and (coeff[i][j] == 1 or coeff[i][j] == -1):
                    flag1 = 1
                    break
            if flag1 == 1:
                for k in range(len(monomial[i])):
                    set_monomial[i][k] = monomials(monomial[i][k],0,0)
                    dict_count[monomial[i][k]] = dict_count[monomial[i][k]] - 1
        round = round + 1

    for i in range(len(set_monomial)):
        for j in range(len(set_monomial[i])):
            dict_merge[set_monomial[i][j].monomial] = dict_merge[set_monomial[i][j].monomial]+set_monomial[i][j].coeff
            dict_coeff[set_monomial[i][j].monomial].append(set_monomial[i][j].h)
    #print(dict_merge) 是用来合并多项式的
    #print(dict_coeff) 是用来存放多项式的系数的
    return dict_merge,dict_coeff

def verify(dict_merge,dict_coeff):    
    solve_left_q = []#存放模拟器可以模拟的多项式
    solve_left_xy = []#存放模拟器不能模拟的多项式
    right_q = []#存放模拟器可以模拟的系数
    right_xy = []#解出余下的系数值
    solve_left_sub = []#将解出来的系数代入多项式中
    for i in dict_merge.keys():
        a = str(dict_merge[i]).replace("q_ij",'')
        if (('x' not in a) and ('y' not in a)):
            if dict_merge[i] != 0 :
                solve_left_q.append(dict_merge[i])
                right_q = right_q + dict_coeff[i]
        else:
            if dict_merge[i] != 0 :
                solve_left_xy.append(dict_merge[i])
                right_xy = right_xy + dict_coeff[i]
    #先求带q的多项式的系数
    Kernel_q = sy.solve(solve_left_q,right_q)
    #判断
    if Kernel_q == []:
        return "FAIL!"
    print("Kernel_q:")
    print(Kernel_q)
    right_tmp = list(set(right_xy) - set(right_q))
    #代入
    for i in solve_left_xy:
        solve_left_sub.append(i.subs(Kernel_q))
    #再次求解系数
    Kernel = sy.solve(solve_left_sub,right_tmp)
    #判断
    print("Kernel_xy:")
    print(Kernel)
    #判断
    if right_tmp!= [] and Kernel == []:
        return "FAIL!"
    #判断
    flag = 1
    for i in Kernel.values():
        if i == 0:
            continue
        if ('x' in str(i) and 'y' not in str(i)) or ('x' in str(i) and 'y' not in str(i)):
            flag = 0
        if 'x_i*y_j' in str(i):
            a = str(i).replace("*x_i*y_j",'')
            if 'q_ij' in a:
                flag = 1
            else:
                for j in Kernel_q.values():
                    if str(j) in (a+'/q_ij'):
                        flag = 1
                        break
        if(flag == 0):
            break
    if flag == 0:
        return "FAIL!"
    #再次代入系数，这次如果都是0，则可以通过验证，否则验证失败
    solve_verify_0 = []
    for i in solve_left_sub:
        solve_verify_0.append(i.subs(Kernel))
    print("solve_verify_0:")
    print(solve_verify_0)
    for i in solve_verify_0:
        if i != 0:
            return "FAIL!" 
    return "PASS!"
 
#启动函数
def run(choose):
    G1_poly, G2_poly,offset_poly = read(choose+'.txt')
    print("=================Read Finished!========================")
    start = time.time()
    GT_poly = parametric_completion(G1_poly,G2_poly,offset_poly)
    print("========Monomial_combination Finished!=================")
    dict_merge,dict_coeff = merge(GT_poly)
    print("==============Merge Finished!==========================")
    result = verify(dict_merge,dict_coeff)
    print("==============Verify Finished!=========================")
    print(result)
    print("=========================Time==========================")
    end = time.time()
    print(end-start)

#测试
if __name__ == '__main__':
    choose = input("Enter a number: ")
    if choose == 'Wee20-new':
        run(choose)
    if choose == 'BCFG17':
        run(choose)
    if choose == 'RPB+19':
        run(choose)
    if choose == 'Wee20':
        run(choose)
    if choose == 'GQ-1-new':
        run(choose)
    if choose == 'GQ-1':
        run(choose)
    if choose == 'GQ-2':
        run(choose)  



        




        
