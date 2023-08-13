#先不考虑offset的问题，这个之后可以再进行考虑
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
    test_set = []
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
    #分解单独的x和y项
    poly = vector(x_var)*vector(y_var)*vector(q_var)
    poly_str = str(poly[0])
    if poly_str.strip('q_ij*x_i*y_j') == '':
        test_set.append('q_ij*x_i*y_j')
    else:
        tmp = poly_str.split('q_ij*x_i*y_j')[0][-1].split('+')
        test_set.append(tmp)
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
    return G1_poly,G2_poly,test_set,offset_poly

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
#提取GT中的变量组合.
def monomial_combination(GT_poly):
    monomial = []
    poly = []
    base = []
    dict = {}
    for i in GT_poly:
        monomial.append(i.monomials())
#处理平方项问题，如果发现有大于2的项，则在后面添加一个新的只包含它的多项式，这样的话就不会被删掉了
    for i in range(len(monomial)):
        for j in range(len(monomial[i])):
            str_var = str(monomial[i][j])
            for k in str_var.split('*'):
                if '_' in k and '^' in k:
                    monomial.append([monomial[i][j]])
    r = copy.deepcopy(monomial)
#设置一个字典，删除只出现了一次的项
    for i in range(len(r)):
        for j in range(len(r[i])):
            dict[r[i][j]] = 0
        
    round = 0
    for i in range(len(r)):
        for j in range(len(r[i])):
            dict[r[i][j]] = dict[r[i][j]] + 1
            if dict[r[i][j]] > round:
                round = dict[r[i][j]]
#================round rotation====================
    g = 0
    while(g < round):
        for i in range(len(r)):
            flag = 0
            for j in range(len(r[i])):
                if dict[r[i][j]] < 2:
                    for k in range(len(r[i])):
                        dict[r[i][k]] = dict[r[i][k]]-1
                    flag = 1
                    break
            if flag == 1:
                r[i] = []
        g = g + 1
    #得到多项式和基
    Poly_simplify = []
    for i in range(len(GT_poly)):
        if r[i] != []:
            Poly_simplify.append(GT_poly[i])
    for i in range(len(GT_poly),len(r)):
        if r[i] != []:
            Poly_simplify.append(r[i][0])
    
    for i in range(len(r)):
        for j in range(len(r[i])):
            if r[i][j] not in base:
                base.append(r[i][j])
    if base == []:
        return "FAIL!"
    return Poly_simplify

def merge(Poly_simplify):
    poly_monomials = []
    poly_coeff = []
    Merge = []#将系数赋值
    right = []
    coeff_tmp = []
    solve_left = []
    dict = {}
    for i in Poly_simplify:
        poly_monomials.append(i.monomials())
        poly_coeff.append(i.coefficients())
    #将系数赋值，并合并相同项的系数。
    for i in range(len(poly_monomials)):
        for j in range(len(poly_monomials[i])):
            tmp = polys(poly_monomials[i][j],var('h'+str(i))*poly_coeff[i][j])
            Merge.append(tmp)
        right.append(var('h'+str(i)))
    #得到要解的方程式和系数
    for i in Merge:
        dict[i.p] = []
    for i in Merge:
        dict[i.p].append(i.coeff)
    
    for value in dict.values():
        coeff_tmp.append(value)

    for i in range(len(coeff_tmp)):
        tmp = 0
        for j in range(len(coeff_tmp[i])):
            tmp = tmp + coeff_tmp[i][j]
        solve_left.append(tmp)
    return solve_left,right

#验证函数
def verify(solve_left,right,result,test_set):
    flag = 0
    for key in result.keys():
        if 'q_ij*x_i*y_j' in str(result[key]):
            right.remove(key)
            break
    Kernel = sy.solve(solve_left,right)
    print("Kernel:%s"%Kernel)
    tmp_result = []
    for i in Kernel.values():
        tmp_result.append(str(i))
    for i in tmp_result:
        if 'q_ij*x_i*y_j' in i:
            flag =1
    if flag == 0:
        return 'FAIL!'
    for i in tmp_result:
        for j in test_set:
            i = i.replace(j,'')
            i = i.replace("q_ij",'').replace("*",'')
            if 'x_i' in i or 'y_j' in i:
                return "FAIL!"
    return "PASS!"

#启动函数
def run(choose):
    G1_poly, G2_poly,test_set,offset_poly = read(choose+'.txt')
    print("=================Read Finished!========================:")
    start = time.time()
    GT_poly = parametric_completion(G1_poly,G2_poly,offset_poly)
    print("========Parametric_completion Finished!================:")
    print("GT_poly:%s"%GT_poly)
    Poly_simplify = monomial_combination(GT_poly)
    print("========Monomial_combination Finished!=================:")
    print("Poly_simplify:%s"%Poly_simplify)
    if Poly_simplify == []:
        return "FAIL!"
    print("==============Merge Finished!==========================:")
    solve_left,right = merge(Poly_simplify)
    ans = sy.solve(solve_left,right)
    if ans == []:
        return "FAIL!"
    Result = verify(solve_left,right,ans,test_set)
    print("==============Verify Finished!=========================:")
    print(Result)
    print("=========================Time==========================:")
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



        
