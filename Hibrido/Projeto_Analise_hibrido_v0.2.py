# -*- coding: utf-8 -*-
"""
        NOVA VERSÃO - DESCRIÇÃO AINDA SERÁ EDITADA!

Código que calcula a: Distribuição de temperatura
Deslocamentos 
Deformações (Normal elastica e Termica)
Tensões (Normal e Cisalhante)

Geração de malha no GMSH e plotagem dos resultados no Paraview.

OBSERVAÇÕES:
>Restrições, forças e propriedades devem ser impostas em 'Physical Groups'(PG).

>Sempre defina a superficie da chapa com PG de superficie 
 com o termo 'superficie'

>As dimensões da peça no GMSH devem ser em METROS!

>Imposição de temperatura: T 50 ~Physycal Group de linha~
>Imposição de geração de calor: pontual 50 ~Physical Group de ponto
                                uniforme 5e+7 ~Physical Group de superficie
                                
>Imposição de propriedades: Prop: alfa = 12e-6 k = 20 ou kx = 20 ky = 25
 (Condutividade térmica em W/m.K) ~PG de superficie~
 
>Imposição de restrição: desloc x;y (x e y são os valores dos deslocamentos
 na direção x e y, respectivamente). ~PG de linha ou ponto~

>Imposição de força: fd Fmin;Fmax E/OU fp Fx;Fy (fd = força distribuida, 
 fp = força pontual e Fx/Fy é a força na direção x/y). ~PG de linha ou ponto~
 
@author: Jose Gomes
"""

import numpy as np
import time

inicio = time.time()
arquivo = open(r"c:\users\jose\desktop\GMSH\testes\Hibrido\tensao_termica_sem_furo.msh")
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

"-----Dados de Entrada-----"
dados = []
local = ''
Qpontual=[]
Quniforme = []
Fluxoq = []
convec = []
restric = []
forca = [] #vetor ['tag',Fx, Fy]-pontual ['tag',Fmax,Fmin]-distrib
tag_fp = [] #tag das forças pontuais no arquivo .msh
tag_fd = [] #tag das forças distrib. no arquivo .msh
temp = [] #lista de tag e temperatura ['tag', 'valor']
for i in range(len(nomes)):
    linha = nomes[i].split()
    if linha[2] == '"Prop:':
        dados = linha[3:]
        local = str(linha[1])
    if linha[-2] == '"pontual': #Recolhe as info de fonte de calor PONTUAL
        Qpontual.append([linha[1],linha[-1]])
    elif linha[-2] == '"uniforme': #Recolhe as info de fonte de calor UNIFORME
        Quniforme.append([linha[1],linha[-1]])
    if linha[2] == '"q':
        Fluxoq.append([linha[1],linha[-1]])
    if linha[2] == '"h':
        convec.append([linha[1],linha[3],linha[5]])
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
    if linha[2] == '"T':
        for char in '"':
            item = linha[-1]
            item = item.replace(char,'')
        temp.append([linha[1],float(item)])

for simb in '"':
    item = dados[-1]
    item = item.replace(simb,'')
    dados[-1] = item

"-----Propriedades do material-----"

for i in range(len(dados)):
    if dados[i] == 'alfa': #Coeficiente de dilatação linear [1/K ou 1/ºC]
        alfa = float(dados[i+2])
    if dados[i] == 'E': #Modulo de Elasticidade ou Young [Gpa]
        E = float(dados[i+2])
    if dados[i] == 'v': #Coeficiente de Pisson (menor do que 1)
        v = float(dados[i+2])
    if dados[i] == 'k': #Coeficiente de Condutividade térmica [W/m.K]
        k = float(dados[i+2])
    if dados[i] == 'ro':
        ro = float(dados[i+2])
    if dados[i] == 'Ti':
        Ti = float(dados[i+2])
g = 9.80665 #Aceleração da gravidade
"----- Nós-----"
nos = [] # lista que deve receber a lista da posicao dos nos
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
t = zn[0] #Espessura do elemento [m]

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
"----------------Matriz de conectividade IEN----------------"
triang = []
for j in range(len(Elem)): #Gerando lista de strings
    linha = Elem[j].split()
    if linha[3] == local: #Selecionando os elementos triangulares da superficie toda
        triang.append(linha)

elem=len(triang) #Quantidade de elementos na malha
IEN= np.zeros((len(triang),3), dtype = int)        

for i in range(len(triang)):
    termo = triang[i]
    for n in range(3):
        IEN[i][n] = termo[-n-1]

for j in range(len(IEN)): #Mudando posicao das colunas na matriz ien
    q = IEN[j][1] - 1
    IEN[j][1] = IEN[j][0] - 1
    IEN[j][0] = q
    IEN[j][2] -= 1 
"Matriz de conectividade ien para 2 GDL(Graus De Liberdade)"
ien = np.zeros((elem,6), dtype=float)
for i in range(elem):
    m = 0
    for j in range(6):
        if j%2 == 0: #Coluna Par
            ien[i][j] = IEN[i][m]*2
        if j%2 != 0: #Coluna IMPAR
            ien[i][j] = IEN[i][m]*2 + 1
            m += 1

######################## INICIO PARTE TEMPERATURA ##################
'---------Fonte de Calor------------'
"---------------------Fonte de Calor------------------------"
pontos = [] #coordenada dos nós onde possui fonte ou sumidouro PONTUAL
NosUni = [] #Lista da lista dos nós por elemento onde possui fonte/sumidouro UNIFORME

for i in range(len(Qpontual)):
    for char in '"': #Retira o caractere " do numero
        item = Qpontual[i][1] 
        item = item.replace(char,'')
        Qpontual[i][1] = float(item)
    for j in range(len(Elem)):
        linha = Elem[j].split()
        linha2 = Qpontual[i]
        if linha[-3] in linha2:
            linha2[0] = int(linha[-1]) - 1 #Adiciono o nó na lista Qpontual
    a = xn[Qpontual[i][0]]
    b = yn[Qpontual[i][0]]
    pontos.append([a,b])

for i in range(len(Quniforme)):
    for char in '"': #Retira o caractere " do numero
        item = Quniforme[i][1] 
        item = item.replace(char,'')
        Quniforme[i][1] = float(item)
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == Quniforme[i][0]:
            a = int(linha[-2]) - 1
            b = int(linha[-1]) - 1
            c = int(linha[-3]) - 1
            NosUni.append([a,b,c])

"--------------Fluxo de Calor------------"
"Vetor de 'força'"
Ft = np.zeros((Npts,1)) #Vetor forca da Temperatura
Fq = [0]*Npts
for i in range(len(Fluxoq)):
    for char in '"': #Retira o caractere " do numero
        item = Fluxoq[i][1] 
        item = item.replace(char,'')
        Fluxoq[i][1] = float(item)
    termo = []; termo.append(Fluxoq[i][0])
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == Fluxoq[i][0]:
            b = int(linha[-1]) - 1
            a = int(linha[-2]) - 1
            if a not in termo:
                termo.append(a)
            if b not in termo:
                termo.append(b)
    for m in range(1,len(termo)-1):
        no1 = termo[m]; no2 = termo[m+1]
        dx = xn[no1] - xn[no2]; dy = yn[no1]-yn[no2]
        l=np.sqrt(dx**2 + dy**2)
        q = Fluxoq[i][1] 
        fq = q*l*t*0.5
        Fq[no1] += fq; Fq[no2] += fq

"------------termo convectivo--------------"
Fh = [0]*Npts
nos_Fh = []
elem_Fh = []
for i in range(len(convec)): #[tag,'valor de h', 'valor de Ta']
    for char in '"':
        item = convec[i][-1]
        item = item.replace(char,'')
        convec[i][-1] = float(item); convec[i][-2] = float(convec[i][-2])
    termo = []; termo.append(convec[i][0])
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == convec[i][0]:
            b = int(linha[-1]) - 1
            a = int(linha[-2]) - 1
            if a not in termo:
                termo.append(a)
            if b not in termo:
                termo.append(b)
            if b == termo[1]: # Condicao de convec numa superficie fechada
                termo.append(b)
    
    nos_Fh.append(termo)
    for m in range(1,len(termo)-1):
        no1 = termo[m]; no2 = termo[m+1]
        dx = xn[no1] - xn[no2]; dy = yn[no1]-yn[no2]
        l=np.sqrt(dx**2 + dy**2)
        h = convec[i][1] ; Ta = convec[i][2] 
        fh = h*Ta*l*t*0.5
        Fh[no1] += fh; Fh[no2] += fh
    
    termo2 = []; termo2.append(convec[i][0])
    for o in range(1,len(termo)-1):
        c = termo[o]; d = termo[o+1]
        for ii in range(elem):
            if c in IEN[ii] and d in IEN[ii]:
                termo2.append(ii)        
    elem_Fh.append(termo2)

for i in range(Npts):
    Ft[i][0] += Fq[i] + Fh[i]

'-----------Temperatura no contorno-----------'
contorno = []
TCC = [0]*Npts #Lista das temperaturas em toda a malha
tbc = []
for i in range(len(temp)):
    tag_PG = temp[i][0]
    contorno.append([tag_PG])
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == tag_PG:
            no1 = int(linha[-2]) -1 ; no2 = int(linha[-1])-1
            if no1 not in contorno[-1]:
                contorno[-1].append(no1)
            if no2 not in contorno[-1]:
                contorno[-1].append(no2)

for i in range(len(temp)):
    linha = temp[i]
    tag = temp[i][0]; valor = temp[i][1]
    for j in range(len(contorno)):
        if tag == contorno[j][0]:
            for l in range(1,len(contorno[j])):
                tbc.append([contorno[j][l], valor])

for i in range(len(tbc)):
    pivo = tbc[i]
    for j in range(len(tbc)):
        if pivo[0] == tbc[j][0] and i != j:
            T1 = pivo[1]; T2 = tbc[j][1]
            if T1 < T2:
                tbc[j][0] = 'null'
            if T2 < T1:
                pivo[0] = 'null' 

tbc = [i for i in tbc if i[0] != 'null']
for i in tbc:
    TCC[i[0]] = i[1]
"----------------------------------------------------"
"Matriz de Rigidez Global K"
Kt = np.zeros((Npts,Npts))

for e in range(elem):
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    
    "Calculo da area do triangulo A"
    matriz_A = np.ones((3,3), dtype=np.float)
    matriz_A[0][1] = xn[a]
    matriz_A[1][1] = xn[b]
    matriz_A[2][1] = xn[c]
    matriz_A[0][2] = yn[a]
    matriz_A[1][2] = yn[b]
    matriz_A[2][2] = yn[c]
    
    detA = np.linalg.det(matriz_A); AREA = 0.5*detA
    
    "Calculo Matriz local Kt"
    kox = np.zeros((3,3), dtype = np.float)
    koy = np.zeros((3,3), dtype = np.float)
    
    qi = xn[c] - xn[b]
    qj = xn[a] - xn[c]
    qk = xn[b] - xn[a]    
    
    ri = yn[b] - yn[c]
    rj = yn[c] - yn[a]
    rk = yn[a] - yn[b]
    
    kox[0][1] = kox[1][0] = qj*qi
    kox[2][0] = kox[0][2] = qk*qi
    kox[2][1] = kox[1][2] = qk*qj
    kox[0][0] = qi**2
    kox[1][1] = qj**2
    kox[2][2] = qk**2
    
    koy[0][1] = koy[1][0] = rj*ri
    koy[2][0] = koy[0][2] = rk*ri
    koy[2][1] = koy[1][2] = rk*rj
    koy[0][0] = ri**2
    koy[1][1] = rj**2
    koy[2][2] = rk**2
    
    kot = (t/(4*AREA))*k*(kox + koy) #Matriz de rigidez local final
    
    matrizh = np.zeros((3,3), dtype = int)
    matrizh[0][0] = matrizh[1][1] = matrizh[2][2] = 2
    matrizh[1][0] = matrizh[0][1] = matrizh[2][1] = matrizh[1][2] = matrizh[0][2] = matrizh[2][0] = 1
    for i in range(len(elem_Fh)):
        pivo = elem_Fh[i]
        if e in pivo:
            tag = pivo[0]
            h = convec[i][1]
            if a in nos_Fh[i] and b in nos_Fh[i]:
                l = np.sqrt(qk**2 + rk**2)
                corte = 2
            elif b in nos_Fh[i] and c in nos_Fh[i]:
                l = np.sqrt(qi**2 + ri**2)
                corte = 0
            elif c in nos_Fh[i] and a in nos_Fh[i]:
                l = np.sqrt(qj**2 + rj**2)
                corte = 1
            
            for j in range(3):
                for w in range(3):
                    if j == corte or w == corte:
                        matrizh[j][w] = 0
            kh = (h*t*l)/6 * matrizh
            kot += kh
    
    for i in range(3):
        ii = int(IEN[e][i])
        for j in range(3):
            jj = int(IEN[e][j])
            Kt[ii][jj] += kot[i][j]

    for i in range(len(Qpontual)):
        if a==Qpontual[i][0] or b==Qpontual[i][0] or c==Qpontual[i][0]:
            xo = pontos[i][0] ; yo = pontos[i][1]
            N = np.ones((3,1))
            si = xn[b]*yn[c] - xn[c]*yn[b]
            sj = xn[c]*yn[a] - xn[a]*yn[c]
            sk = xn[a]*yn[b] - xn[b]*yn[a]
            N[0][0] = (si + ri*xo + qi*yo)
            N[1][0] = (sj + rj*xo + qj*yo)
            N[2][0] = (sk + rk*xo + qk*yo)
            G = Qpontual[i][1]
            fu = (1/(2*AREA))*G*t*N
            for j in range(3):
                ii = int(IEN[e][j])
                Ft[ii][0] += fu[j][0]
    
    for i in range(len(NosUni)):
        for j in range(len(Quniforme)):
            termo=NosUni[i]
            if a == termo[0] and b == termo[1] and c == termo[2]:
                one = np.ones((3,1))
                G = Quniforme[j][1]
                fu = one * (G*AREA*t)/3
                for m in range(3):
                    ii = int(IEN[e][m])
                    Ft[ii][0] += fu[m][0]
"-------------------------------------------------------------"
"Condicoes de contorno"
for i in range(len(tbc)):
    ind = tbc[i][0]; valor= tbc[i][1]
    Kt[ind][ind] *= 1e10
    Ft[ind][0] = Kt[ind][ind] * valor

"Invertendo a matriz global K"
invKt = np.linalg.inv(Kt)

"Calculo do vetor das temperaturas"
T = np.dot(invKt, Ft)
To= []
for i in range(Npts):
    To.append(float(T[i]))
T = To

"Calculo do gradiente de temperatura (Grad) e do fluxo de Calor (flux)"
Grad = np.zeros((elem,2), dtype = float)
flux = np.zeros((elem,2), dtype = float)
for e in range(elem):
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    
    B = np.zeros((2,3), dtype = float)
    temp = np.zeros((3,1), dtype = float)
    "Calculo da area do triangulo A"
    matriz_A = np.ones((3,3), dtype=np.float)
    matriz_A[0][1] = xn[a]; matriz_A[0][2] = yn[a]
    matriz_A[1][1] = xn[b]; matriz_A[1][2] = yn[b]
    matriz_A[2][1] = xn[c]; matriz_A[2][2] = yn[c]
    
    detA = np.linalg.det(matriz_A)    
    qi = (xn[c] - xn[b])*(1/detA) ; ri = (yn[b] - yn[c])*(1/detA)
    qj = (xn[a] - xn[c])*(1/detA) ; rj = (yn[c] - yn[a])*(1/detA)
    qk = (xn[b] - xn[a])*(1/detA) ; rk = (yn[a] - yn[b])*(1/detA)
    
    B[0][0] = ri      ;B[0][1] = rj        ;B[0][2] = rk 
    B[1][0] = qi      ;B[1][1] = qj        ;B[1][2] = qk 
    temp[0][0] = T[a] ;temp[1][0] = T[b]   ;temp[2][0] = T[c]
    
    grad = np.dot(B,temp)
    q = -k*grad

    for i in range(2):
        Grad[e][i] = grad[i][0]
        flux[e][i] = q[i][0]

"Variação da temperatura (dT) em cada elemento"
dT = [] #Variação da temperatura no elemento
for e in range(elem):
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    dT1 = T[a] - Ti; dT2 = T[b] - Ti; dT3 = T[c] - Ti 
    Media = (dT1 + dT2 + dT3)/3
    dT.append(Media)

"Deformação inicial devido a variacao de temperatura (DEFt)"
DEFt = np.zeros((elem,3)) #[[EPSx],[EPSy],[GAMMAxy]]
for i in range(elem):
    DEFt[i][0] = alfa*dT[i] #9.36e-4
    DEFt[i][1] = alfa*dT[i]

##################### INICIO PARTE ESTRUTURA ###########################
G = np.zeros((3,6)) #Matriz auxiliar na Matriz Deslocamento-Deformação
G[0][1] = G[1][5] = G[2][2] = G[2][4] = 1

D = np.zeros((3,3)) #Matriz de Elasticidade (3x3)
D[0][0] = D[1][1] = 1 ; D[1][0] = D[0][1] = v
D[2][2] = (1-v)/2
D *= E/(1-v**2)

"----------------Deslocamentos conhecidos------------------"
nos_restric = []
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
Fu = np.zeros((2*Npts,1)) #Matriz Global Força
nos_forca = []

for i in range(len(forca)):
    nos_forca.append([forca[i][0]])
    for j in range(len(Elem)):
        linha = Elem[j].split()
        if linha[3] == nos_forca[-1][0]:
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
"Força pontual"
for i in range(len(tag_fp)): #Só é valido para força pontual 'fp'
    pos = nos_forca[i][0]
    no = int(nos_forca[i][1])
    for j in range(len(forca)):
        f = forca[j][0]
        if pos == f:
            Fu[2*no][0] += forca[j][1]
            Fu[2*no+1][0] += forca[j][2]

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
            Fu[pos*2][0] += Px[m][0] + Vx[m][0]
            Fu[pos*2+1][0] += Py[m][0]+ Vy[m][0]


fb = -g*ro
"---------------------------------------------------------"
"Matriz de Rigidez Global K"
Ku = np.zeros((2*Npts,2*Npts))
"Vetor Força referente a variação da temperatura dT"
TETA = np.zeros((2*Npts,1))

for e in range(elem):
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    epst = np.array(DEFt[e]).reshape(3,1)
    
    "Calculo da area do triangulo A"
    matriz_A = np.ones((3,3), dtype=np.float)
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
    prod3 = np.dot(prod1,epst) #matriz 6x1

    ku = prod2*AREA*t  #Unidade: N/m
    for i in range(6):
        ii = int(ien[e][i])
        for j in range(6):
            jj = int(ien[e][j])
            Ku[ii][jj] += ku[i][j]
        
    teta = t*AREA*prod3 #Força termica saindo em N (6x1)        
    incr = 1
    for j in range(3): #Loop referente a forca
        x = NosElem[j]
        TETA[2*x][0] += teta[2*j]
        TETA[2*x+1][0] += teta[2*j+1]

        i = int(ien[e][j + incr])
        Fu[i][0] += (t*fb*AREA)/3
        incr += 1
"----------------------------------------"
Fu += TETA

for i in range(len(bc)):
    ind = bc[i][0]; valor = bc[i][1]
    Ku[ind][ind] *= 1e10
    Fu[ind][0] = Ku[ind][ind] * valor

invKu = np.linalg.inv(Ku)
"Calculo do deslocamento"
U = np.dot(invKu,Fu)
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
    delta = np.zeros((6,1)) #Matriz deslocamento do elemento
    a = IEN[e][0]
    b = IEN[e][1]
    c = IEN[e][2]
    
    epst = np.array(DEFt[e]).reshape(3,1)
    NosElem = [a,b,c]
    A = np.zeros((6,6))
    j = 0
    for i in range(3):
        A[i+j][0] = A[i+j+1][3] = 1
        A[i+j][1] = A[i+j+1][4] = xn[NosElem[i]]
        A[i+j][2] = A[i+j+1][5] = yn[NosElem[i]]
        j += 1
    
    A_1 = np.linalg.inv(A)
    B = np.dot(G,A_1) #* 1e-3 #Matriz Deslocamento - Deformação (unidade 1/m)    
    for i in range(len(ien[0])):
        delta[i][0] = Uo[int(ien[e][i])]
    
    eps = np.dot(B,delta) #Deformação total dentro do elemento 'e'
    epsR = eps - epst #Vetor da deformação real
    sig = np.dot(D, epsR) #Tensão real dentro do elemento 'e'[Pa]
    p = (sig[0][0])**2 - sig[0][0]*sig[1][0] + (sig[1][0])**2 + 3*(sig[2][0])**2
    vm = np.sqrt(p); VM.append(vm)
    
    for i in range(3):
        DEF[e][i] += epsR[i][0]
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

fator = 1#e+7*float(max(Uo))  #Fator de deformação visivel

x = [] #Nova posição dos nós na direção x  
y = [] #Nova posição dos nós na direção y
for i in range(len(xn)):
    a = float(xn[i])
    b = float(yn[i])
    A = float(Uo[2*i])
    B = float(Uo[2*i+1]) 
    x.append(a+A*fator)
    y.append(b+B*fator)

"------------------INICIO-DO-PARAVIEW-----------------------------"
arq = open('TensaoTermica_sem_furo.vtk','w')

texto = '''# vtk DataFile Version 4.2
2D Análise Estrutural e Termico
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
arq.write('VECTORS Deformacao_Termica float')
for d in range(len(DEFt)):
    arq.write('\n')
    arq.write(str(DEFt[d][0]))
    arq.write(' ')
    arq.write(str(DEFt[d][1]))
    arq.write(' ')
    arq.write(str(DEFt[d][2]))

arq.write('\n')
arq.write('VECTORS Deformacao float')
for d in range(len(EPSx)):
    arq.write('\n')
    arq.write(str(EPSx[d]))
    arq.write(' ')
    arq.write(str(EPSy[d]))
    arq.write(' ')
    arq.write(str(EPSxy[d]))

arq.write('\n')
arq.write('VECTORS Tensao float')
for t in range(len(TENx)):
    arq.write('\n')
    arq.write(str(round(TENx[t],5)))
    arq.write(' ')
    arq.write(str(round(TENy[t],5)))
    arq.write(' ')
    arq.write(str(round(TENxy[t],5)))


arq.write('\n\n')
arq.write('VECTORS TensaoVM[Pa] float')
for i in range(len(VM)):
    arq.write('\n')
    arq.write(str(round(VM[i],5)))
    arq.write(' ')
    arq.write('0.0')
    arq.write(' ')
    arq.write('0.0')
    
arq.write('\n')
arq.write('VECTORS Gradiente float')
for d in range(len(Grad)):
    arq.write('\n')
    arq.write(str(Grad[d][0]))
    arq.write(' ')
    arq.write(str(Grad[d][1]))
    arq.write(' ')
    arq.write('0')

arq.write('\n\n')
arq.write('VECTORS Fluxo float')
for i in range(len(flux)):
    arq.write('\n')
    arq.write(str(flux[i][0]))
    arq.write(' ')
    arq.write(str(flux[i][1]))
    arq.write(' ')
    arq.write('0')


arq.write('\n\n')
arq.write('POINT_DATA ') #Vetor velocidade: no caso de solidos, v = 0
arq.write(str(Npts))
arq.write('\n')
arq.write('VECTORS Deslocamento float') #O nome Deslocamento pode ser muodificado.
for v in range(len(Uox)):
    arq.write('\n')
    arq.write(str(round(Uox[v],5))) #Em milimetros
    arq.write(' ')
    arq.write(str(round(Uoy[v],5))) #Em milimetros
    arq.write(' ')
    arq.write('0.0')

arq.write("\n\n")
arq.write('SCALARS Temperatura double') #O nome Temperatura pode ser modificado.
arq.write('\n')
arq.write('LOOKUP_TABLE default')

for temp in range(len(T)):
    arq.write('\n')
    arq.write(str(T[temp]))
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
#"""