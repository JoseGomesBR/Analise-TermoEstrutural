    # -*- coding: utf-8 -*-
"""
Código que calcula deslocamentos, deformações e tensoes em um elemento
bidimensional (chapa) utilizando o programa GMSH e tem estrutura 
semelhante ao arquivo de mesmo nome versao 1.0

OBSERVAÇÕES:
>Restrições, forças e propriedades devem ser impostas em 'Physical Groups'(PG).

>A espessura é dada como a terceira dimensão no programa GMSH

>As dimensões da chapa no GMSH DEVEM ser em METROS!

>Imposição de restrição: desloc x;y (x e y são os valores dos deslocamentos
 na direção x e y, respectivamente). ~PG de linha ou ponto~

>Imposição de força: fd Fmin;Fmax E/OU fp Fx;Fy (fd = força distribuida(TENSÃO), 
 fp = força pontual e Fx;Fy é a força na direção x;y). ~PG de linha ou ponto~

>Imposição das propriedades do material: Prop: E = 200 alfa = 12e-6 k = 20 v= 0.3
por exemplo. SEMPRE coloque 'Prop:' para definir a superficie


            VERSÃO FINAL!!!!

######LEMBRAR DE FAZER TESTE DE CISALHAMENTO COM FORÇA PONTUAL#######

@author: Jose
"""
import time 
import numpy as np
"-----------------------------------------------------------"
inicio = time.time()
arquivo = open(r"c:\users\jose\desktop\GMSH\Testes\Estrutura\Validacao\chapa_tensao.msh")
lista = arquivo.readlines()

nomes = []
a = 0
for i in range(len(lista)):
    termo = lista[i]
    if termo == '$PhysicalNames\n':
        a = 1
    elif termo == '$EndPhysicalNames\n':
        a = 0
    if a == 1:
        nomes.append(lista[i])
nomes.pop(0)
nomes.pop(0)

"-----Propriedades do material-----"
dados = []
local = ''
restric = []
forca = [] #vetor ['tag',Fx, Fy]-pontual ['tag',Fmax,Fmin]-distrib
tag_fp = [] #tag das forças pontuais no arquivo .msh
tag_fd = [] #tag das forças distrib. no arquivo .msh

for i in range(len(nomes)):
    linha = nomes[i].split()
    if linha[2] == '"Prop:':
        dados = linha[3:]
        local = str(linha[1])
    if linha[2] == '"desloc':
        for char in '"':
            item = linha[-1]
            item = item.replace(char,'')
        valor = item.split(';')
        valor.insert(0,linha[1])       
        restric.append(valor)
    if linha[2] == '"fp':
        k = linha[1]
        tag_fp.append(k)
        for char in '"':
            item = linha[-1]
            item = item.replace(char,'')
        valor = item.split(';')
        valor[0] = float(valor[0]) ; valor[1] = float(valor[1])
        valor.insert(0,linha[1])       
        forca.append(valor)
    if linha[2] == '"fd':
        k = linha[1]
        tag_fd.append(k)
        for char in '"':
            item = linha[-1]
            item = item.replace(char,'')
        valor = item.split(';')
        valor[0] = float(valor[0]) ; valor[1] = float(valor[1])
        valor.insert(0,linha[1])
        forca.append(valor)

for simb in '"':
    item = dados[-1]
    item = item.replace(simb,'')
    dados[-1] = item

for i in range(len(dados)):
    if dados[i] == 'E':
        E = float(dados[i+2]) #Modulo de elasticidade
    if dados[i] == 'v': 
        v = float(dados[i+2]) #Coeficiente de Poisson
    if dados[i] == 'ro':
        ro = float(dados[i+2]) #Densidade do material

g = 1000 #9.80665
"----Nós----"
nos = [] # lista que deve receber a posicao dos nos
a = 0
for i in range(len(lista)):
    termo = lista[i]
    if termo == '$Nodes\n':
        a = 1
    if termo == '$EndNodes\n':
        a = 0
    if a == 1:
        nos.append(lista[i])
        
nos.pop(0) #Retira '$Nodes\n'
nos.pop(0) #Retira quantidade de nós

xn = [] #lista de coordenadas dos nos em X
yn = [] #lista de coordenadas dos nos em Y
zn = []

for j in range(len(nos)):
    linha = nos[j].split()
    xn.append(float(linha[1]))
    yn.append(float(linha[2]))
    zn.append(float(linha[3]))
t = zn[0] #Espessura do elemento

Npts = len(xn) #Quantidade de pontos/nós na malha
"---------------------------------------------------------------"
Elem = [] #Lista que recebe valores iniciais da matriz IEN
b = 0
for i in range(len(lista)):
    termo = lista[i]    
    if termo == '$Elements\n':
        b = 1
    if termo == '$EndElements\n':
        b = 0
    if b == 1:
        Elem.append(lista[i])
     
Elem.pop(0) #Retida o $Elements
Elem.pop(0) #Retira o numero de elementos da lista

arquivo.close()
"----------------------------------------------------"
triang = []
for j in range(len(Elem)): #Gerando lista de strings
    linha = Elem[j].split()
    if linha[3] == local: #Selecionando os elementos triangulares da superficie toda
        triang.append(linha)

elem=len(triang) #Quantidade de elementos bidimensionais na malha
IEN= np.zeros((len(triang),3), dtype = np.int)        

for k in range(len(triang)):
    termo = triang[k]
    for n in range(3):
        IEN[k][n] = termo[-n-1]

for j in range(len(IEN)): #Mudando posicao das colunas na matriz ien
    q = IEN[j][1] - 1
    IEN[j][1] = IEN[j][0] - 1
    IEN[j][0] = q
    IEN[j][2] -= 1 

"Matriz de conectividade ien para 2 GDL(Graus De Liberdade)"
ien = np.zeros((elem,6), dtype=float)
for i in range(elem):
    k = 0
    for j in range(6):
        if j%2 == 0: #Coluna Par
            ien[i][j] = IEN[i][k]*2
        if j%2 != 0: #Coluna IMPAR
            ien[i][j] = IEN[i][k]*2 + 1
            k += 1

"----------------Deslocamentos conhecidos------------------"
nos_restric = [] #[tag_PG,no1,no2,...]
desloc = [] #Lista com nós e seus deslocamentos exigidos
for i in range(len(Elem)):
    linha = Elem[i].split()
    for j in range(len(restric)):
        if linha[3] == restric[j][0]:
            a = [linha[3]]
            if linha[-1] == linha[-2]:
                a.append(int(linha[-1]) - 1)
            elif linha[-1] != linha[-2]:
                a.append(int(linha[-1]) - 1)
                a.append(int(linha[-2]) - 1)
            nos_restric.append(a)

for i in range(len(restric)):
    termo = restric[i]
    for j in range(len(nos_restric)):
        termo2 = nos_restric[j]
        if termo[0] == termo2[0]:
            desloc.append([termo2[1],termo[1],termo[2]])
            if len(termo2) > 2:
                desloc.append([termo2[2],termo[1],termo[2]])
bc = []
for i in range(len(desloc)):
    noX = desloc[i][0]*2 ; noY = desloc[i][0]*2 +1
    if desloc[i][1] != '':
        bc.append([noX,float(desloc[i][1])])
    if desloc[i][2] != '':
        bc.append([noY,float(desloc[i][2])])
"--------------Força aplicadas----------------"
F = np.zeros((2*Npts,1))
nos_forca = [] #vetor ['tag',no1, no2,...]
for i in range(len(forca)):
    nos_forca.append([forca[i][0]])
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == nos_forca[-1][0]: #Verifica o ultimo termo add
            termo = nos_forca[-1]
            a = int(linha[-1]) - 1
            b = int(linha[-2]) - 1
            if a != b :
                if b in termo:
                    termo.append(a)
                else:
                    termo.append(b)
                    termo.append(a)
            elif a == b :
                termo.append(a)
    nos_forca[-1] = termo

for i in range(len(tag_fp)): #Só é valido para força pontual 'fp'
    pos = nos_forca[i][0]
    no = int(nos_forca[i][1])
    for j in range(len(forca)):
        f = forca[j][0]
        if pos == f:
            F[2*no][0] += forca[j][1] #Fx
            F[2*no+1][0] += forca[j][2] #Fy


ElemNos = [] #[[tag_elem, nos...],[tag_elem, nos...],...]
for j in range(len(tag_fd)):
    ElemNos.append([])
    PGCC = ElemNos[j]
    for i in range(len(Elem)):
        linha = Elem[i].split()
        if linha[3] == tag_fd[j]:
            nline = linha[4] #tag do elemento de linha        
            if nline not in PGCC:
                PGCC.append(nline)
                PGCC.append(int(linha[5]) - 1)
                PGCC.append(int(linha[6]) - 1)
            else:
                if int(linha[5])-1 not in PGCC:
                    PGCC.append(int(linha[5])-1)
                else:  #if int(linha[6])-1 not in PGCC
                    PGCC.append(int(linha[6])-1)

for i in range(len(ElemNos)):
    termo = ElemNos[i]; tagElem = forca[i][0]
    lado = []
    for j in range(len(termo)):
        if type(termo[j]) == str:
            lado.append([termo[j]])
        else:
            lado[-1].append(termo[j])

    for j in range(len(lado)):
        face=[]
        for m in lado[j]:
            if type(m) == int:
                face.append(m)    

        elementos = len(face) - 1
        vetorL = []
        L = 0
        P = np.zeros((elementos+1, 1), dtype = float) #Vetor forca normal em cada nó [N/m^2]
        V = np.zeros((elementos+1, 1), dtype = float) #Vetor forca tangente
        Px = np.zeros((elementos+1, 1), dtype = float) #[N]
        Vx = np.zeros((elementos+1, 1), dtype = float)
        Py = np.zeros((elementos+1, 1), dtype = float) #[N]
        Vy = np.zeros((elementos+1, 1), dtype = float)
        r = []; s = [] #componente cosseno e seno, respectivamente
        for m in range(len(forca)):
            if forca[m][0] == tagElem:
                P[0:] = forca[m][1]
                V[0:] = forca[m][2]
        
        for m in range(elementos):
            noi = face[m]; noj = face[m+1]
            dx = xn[noi] - xn[noj]; dy = yn[noj] - yn[noi]
            l = np.sqrt(dx**2+dy**2)
            cos = dy/l ; sen = dx/l
            r.append(cos); s.append(sen)
            L += float(l)
            vetorL.append(l)

        for m in range(elementos+1):
            if m == 0:
                Px[m][0] += vetorL[m]*(2*P[m][0] + P[m+1][0])*r[m]
                Py[m][0] += vetorL[m]*(2*P[m][0] + P[m+1][0])*s[m]
                Vx[m][0] += vetorL[m]*(2*V[m][0] + V[m+1][0])*(-s[m])
                Vy[m][0] += vetorL[m]*(2*V[m][0] + V[m+1][0])*r[m]
            elif m == elementos: #ultimo da lista
                Px[m][0] += vetorL[m-1]*(2*P[m][0] + P[m-1][0])*r[m-1]
                Py[m][0] += vetorL[m-1]*(2*P[m][0] + P[m-1][0])*s[m-1]
                Vx[m][0] += vetorL[m-1]*(2*V[m][0] + V[m-1][0])*(-s[m-1])
                Vy[m][0] += vetorL[m-1]*(2*V[m][0] + V[m-1][0])*r[m-1]
            else:
                Px[m][0] = vetorL[m-1]*(2*P[m][0] + P[m-1][0])*r[m-1] + vetorL[m]*(2*P[m][0] + P[m+1])*r[m]
                Py[m][0] = vetorL[m-1]*(2*P[m][0] + P[m-1][0])*s[m-1] + vetorL[m]*(2*P[m][0] + P[m+1])*s[m]
                Vx[m][0] = vetorL[m-1]*(2*V[m][0] + V[m-1][0])*(-s[m-1]) + vetorL[m]*(2*V[m][0] + V[m+1])*(-s[m])
                Vy[m][0] = vetorL[m-1]*(2*V[m][0] + V[m-1][0])*r[m-1] + vetorL[m]*(2*V[m][0] + V[m+1])*r[m]

        Px = (t/6)*Px; Py = (t/6)*Py
        Vx = (t/6)*Vx; Vy = (t/6)*Vy
        xinicio = xn[face[0]]; xfim = xn[face[-1]]
        yinicio = yn[face[0]]; yfim = yn[face[-1]]
        if xinicio - xfim > 0 or yinicio - yfim > 0:
            face.reverse()

        for m in range(len(face)):
            pos = face[m]
            F[pos*2][0] += Px[m][0] + Vx[m][0]
            F[pos*2+1][0] += Py[m][0]+ Vy[m][0]

fb = -g*ro
"---------------------------------------------------------"
G = np.zeros((3,6)) #Matriz auxiliar na Matriz Deslocamento-Deformação
G[0][1] = G[1][5] = G[2][2] = G[2][4] = 1

D = np.zeros((3,3)) #Matriz de Elasticidade
D[0][0] = D[1][1] = 1 ; D[1][0] = D[0][1] = v
D[2][2] = (1-v)/2
D *= E/(1-v**2)

"Matriz de Rigidez Global K"
K = np.zeros((2*Npts,2*Npts))
for e in range(elem):
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    
    "Calculo da area do triangulo" 
    matriz_A = np.ones((3,3), dtype=np.float32)
    matriz_A[0][1] = xn[a]
    matriz_A[1][1] = xn[b]
    matriz_A[2][1] = xn[c]
    matriz_A[0][2] = yn[a]
    matriz_A[1][2] = yn[b]
    matriz_A[2][2] = yn[c]
    detA = np.linalg.det(matriz_A); AREA = 0.5*detA
    
    NosElem = [a,b,c]
    A = np.zeros((6,6))
    j = 0
    for i in range(3):
        A[i+j][0] = A[i+j+1][3] = 1
        A[i+j][1] = A[i+j+1][4] = xn[NosElem[i]]
        A[i+j][2] = A[i+j+1][5] = yn[NosElem[i]]
        j += 1
        
    A_1 = np.linalg.inv(A)
    B = np.dot(G,A_1)
    B_t = np.matrix.transpose(B)
    prod1 = np.dot(B_t,D)
    prod2 = np.dot(prod1,B)

    k = prod2*AREA*t  #Unidade: N/m
    for i in range(6):
        ii = int(ien[e][i])
        for j in range(6):
            jj = int(ien[e][j])
            K[ii][jj] += k[i][j]
    
    incr = 1
    for j in range(3): #Loop referente a forca Peso
        i = int(ien[e][j + incr])
        F[i][0] += (t*fb*AREA)/3
        incr += 1
"-----------------------------------------------------"
for i in range(len(bc)):
    ind = bc[i][0]; valor = bc[i][1]
    K[ind][ind] *= 1e10
    F[ind][0] = K[ind][ind] * valor

invK = np.linalg.inv(K)
"Calculo do deslocamento"
U = np.dot(invK,F)
Uo = [] #Deslocamento em todos os nós

for i in range(len(U)):
    Uo.append(float(U[i]))

Uox = []; Uoy = []
for i in range(Npts):
    Uox.append(float(Uo[2*i]))
    Uoy.append(float(Uo[2*i+1]))
"---------------------------------"
"Calculo das Deformações e Tensões"
DEF = np.zeros((elem,3))
TEN = np.zeros((elem,3))
VM = [] #Tensão de von Mises
for e in range(elem):
    delta = np.zeros((6,1))

    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    
    NosElem = [a,b,c]
    A = np.zeros((6,6))
    j = 0
    for i in range(3):
        A[i+j][0] = A[i+j+1][3] = 1
        A[i+j][1] = A[i+j+1][4] = xn[NosElem[i]]
        A[i+j][2] = A[i+j+1][5] = yn[NosElem[i]]
        j += 1
    
    A_1 = np.linalg.inv(A)
    B = np.dot(G,A_1) #Matriz Deslocamento-Deformação 3x6 [1/m]
    
    for i in range(len(ien[0])):
        delta[i][0] = Uo[int(ien[e][i])]
    
    eps = np.dot(B,delta) #Deformação dentro do elemento 'e'
    sig = np.dot(D,eps) #Tensão dentro do elemento 'e'
    
    p = (sig[0][0])**2 - sig[0][0]*sig[1][0] + (sig[1][0])**2 + 3*(sig[2][0])**2
    vm = np.sqrt(p); VM.append(vm)
    
    for i in range(3):
        DEF[e][i] += eps[i][0]
        TEN[e][i] += sig[i][0]

TENx = []; TENy = []
EPSx = []; EPSy = []
EPSxy = []; TENxy =[]

for i in range(elem):
    TENx.append(TEN[i][0])
    TENy.append(TEN[i][1])
    TENxy.append(TEN[i][2])
    EPSx.append(DEF[i][0])
    EPSy.append(DEF[i][1])
    EPSxy.append(DEF[i][2])

x = [] #Nova posição dos nós na direção x  
y = [] #Nova posição dos nós na direção y

alfa= 1#50

for i in range(len(xn)):
    a = float(xn[i])
    b = float(yn[i])
    A = float(Uo[2*i])
    B = float(Uo[2*i+1]) 
    x.append(a+A*alfa)
    y.append(b+B*alfa)
"-------------------------- INICIO PARAVIEW ---------------------------"
arq = open('chapa_tensao_7928elem.vtk','w')

texto = '''# vtk DataFile Version 4.2
2D Análise Estrutural
ASCII               
DATASET UNSTRUCTURED_GRID
FIELD FieldData 1
NODES 1 2 int
'''
arq.writelines(texto) #Cabeçalho
arq.write(str(Npts))
arq.write(' ')
arq.write(str(elem))
arq.write("\n\n")
arq.write('POINTS ') #Coordenadas
arq.write(str(Npts))
arq.write(' double')

for i in range(Npts):
    a,b = str(round(x[i],5)),str(round(y[i],5))
    arq.write("\n")
    arq.write(a)
    arq.write(' ')
    arq.write(b)
    arq.write(' ')
    arq.write('0.0')

arq.write("\n\n")
arq.write('CELLS ') #Elementos
arq.write(str(elem))
arq.write(' ')
arq.write(str(4*elem))

for i in range(elem):
    linha = IEN[i]
    a,b,c = str(linha[0]),str(linha[1]),str(linha[2])
    arq.write("\n")
    arq.write(str(len(linha)))
    arq.write(' ')
    arq.write(a)
    arq.write(' ')
    arq.write(b)
    arq.write(' ')
    arq.write(c)

arq.write("\n\n")
arq.write('CELL_TYPES ') #Tipo de Elemento da Malha
arq.write(str(elem))

for i in range(elem):
    arq.write("\n")
    arq.write('5') #Este número representa o tipo de elemento da malha
                   #5 = Triangulo  
arq.write("\n\n")
arq.write('CELL_DATA ') #Vetor velocidade: no caso de solidos, v = 0
arq.write(str(elem))
arq.write('\n')
arq.write('VECTORS Deformacao[m/m] float')
for d in range(len(EPSx)):
    arq.write('\n')
    arq.write(str(EPSx[d]))
    arq.write(' ')
    arq.write(str(EPSy[d]))
    arq.write(' ')
    arq.write(str(EPSxy[d]))

arq.write('\n')
arq.write('VECTORS Tensao[Pa] float')
for i in range(len(TENx)):
    arq.write('\n')
    arq.write(str(round(TENx[i],5)))
    arq.write(' ')
    arq.write(str(round(TENy[i],5)))
    arq.write(' ')
    arq.write(str(round(TENxy[i],5)))
    
arq.write('\n')
arq.write('VECTORS TensaoVM[Pa] float')
for i in range(len(VM)):
    arq.write('\n')
    arq.write(str(round(VM[i],5)))
    arq.write(' ')
    arq.write('0.0')
    arq.write(' ')
    arq.write('0.0')
    
arq.write('\n\n')
arq.write('POINT_DATA ') #Vetor velocidade: no caso de solidos, v = 0
arq.write(str(Npts))
arq.write('\n')
arq.write('VECTORS Deslocamento[m] float') #O nome Deslocamento pode ser muodificado.
for v in range(len(Uox)):
    arq.write('\n')
    arq.write(str(round(Uox[v],9)))
    arq.write(' ')
    arq.write(str(round(Uoy[v],9)))
    arq.write(' ')
    arq.write('0.0')

arq.close()

fim = time.time()
tempo = fim - inicio
print('Elementos: ',elem,'   ','Nós: ',Npts,'\n')
print('Tempo de execução: ',tempo,' segundos.')
if tempo >= 60.0:
    minuto = tempo/60
    print('                  ',minuto,' minuto(s).')
if tempo >= 3600.0:
    hora = tempo/3600
    print('                  ',hora,' hora(s).')
