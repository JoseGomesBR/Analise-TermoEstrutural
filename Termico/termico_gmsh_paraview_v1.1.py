"""
Codigo que calcula a DISTRIBUIÇÃO de temperatura e o seu GRADIENTE, utilizando
os programa GMSH e Paraview. 
Aqui as condições impostas são de temperaturas especificas nas bordas e
fonte/sumidouro de calor pontual e/ou uniforme.

OBSERVAÇÕES:
>As dimensões da chapa no GMSH devem ser em METROS!

>Imposição de Condições de Contorno: temperatura- T 50 ~Physycal Group de linha~
                                     FontePontual - pontual 50 ~Physical Group de ponto~
                                     FonteUniforme - uniforme 5e+7 ~Physical Group de superficie
                                     convecção - h 25 ta 70 ~P.G. de linha~
                                     Fluxo de calor - q 100 ~P.G. de linha~
                                     
>Imposição de propriedades: Prop: k = 20 (Condutividade térmica em W/m.K) ~PG de superficie~


VERSÃO FINAL!!!!!!!!! 
"""
import numpy as np
import time

"---------------------------------------------------------------------"
inicio = time.time()
arquivo = open(r"c:\users\jose\desktop\GMSH\Testes\Termica\Validacao\condicao_fonte_pontual.msh")
lista = arquivo.readlines()

nomes = []
a = 0
for i in range(len(lista)):
    termo = lista[i]
    if termo == '$PhysicalNames\n':
        a = 1
    if termo == '$EndPhysicalNames\n':
        a = 0
    if a == 1:
        nomes.append(lista[i])
nomes.pop(0)
nomes.pop(0)

"----- Dados de Entrada -------"
dados = []
local = ''
Qpontual=[]
Quniforme = []
Fluxoq = []
convec = []
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
    if linha[2] == '"T':
        for char in '"':
            item = linha[-1]
            item = item.replace(char,'')
        temp.append([linha[1],float(item)])

for simb in '"':
    item = dados[-1]
    item = item.replace(simb,'')
    dados[-1] = item

for i in range(len(dados)):
    if dados[i] == 'k': #Coeficiente de Condutividade térmica [W/m.K]
        k = float(dados[i+2])

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
"----Elementos-----"
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
"--------------------------------------------------------"
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
F = np.zeros((Npts,1))
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
    F[i][0] += Fq[i] + Fh[i]

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
"--------------------------FINAL DO ARQUIVO MSH--------------------------------"
"Matriz de rigidez global K"
K = np.zeros((Npts,Npts), dtype = np.float)

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
    
    "Matriz de rigidez local"
    kox = np.zeros((3,3), dtype = np.float)
    koy = np.zeros((3,3), dtype = np.float)
    
    qi = xn[c] - xn[b]; ri = yn[b] - yn[c]
    qj = xn[a] - xn[c]; rj = yn[c] - yn[a] 
    qk = xn[b] - xn[a]; rk = yn[a] - yn[b]
    
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

    ko = (t/(4*AREA))*k*(kox + koy) #Matriz de rigidez local final
    
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
            ko += kh

    for i in range(3):
        ii = int(IEN[e][i])
        for j in range(3):
            jj = int(IEN[e][j])
            K[ii][jj] += ko[i][j]
        
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
                F[ii][0] += fu[j][0]
    
    for i in range(len(NosUni)):
        for j in range(len(Quniforme)):
            termo=NosUni[i]
            if a == termo[0] and b == termo[1] and c == termo[2]:
                one = np.ones((3,1))
                G = Quniforme[j][1]
                fu = one * (G*AREA*t)/3
                for m in range(3):
                    ii = int(IEN[e][m])
                    F[ii][0] += fu[m][0]
"--------------------------------------------------------------"
"Condicoes de contorno"
for i in range(len(tbc)):
    ind = tbc[i][0]; valor= tbc[i][1]
    K[ind][ind] *= 1e10
    F[ind][0] = K[ind][ind] * valor

"Invertendo a matriz global K"
invK = np.linalg.inv(K)

"Calculo do vetor das temperaturas"
T = np.dot(invK, F)
To= []
for i in range(Npts):
    To.append(float(T[i]))
T = To

"Calculo do gradiente de T e do fluxo de Calor"
G = np.zeros((elem,2), dtype = float) #Gradiente da temperatura [ºC/m]
flux = np.zeros((elem,2), dtype = float) #Fluxo de calor [ ]
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
    
    B[0][0] =ri      ;B[0][1] =rj      ;B[0][2] =rk 
    B[1][0] =qi      ;B[1][1] =qj      ;B[1][2] =qk 
    temp[0][0] =T[a] ;temp[1][0] =T[b] ;temp[2][0] =T[c]
    
    g = np.dot(B,temp)
    q = -k*g
    
    for i in range(2):
        G[e][i] = g[i][0]
        flux[e][i] = q[i][0]

"------------------INICIO-DO-PARAVIEW-----------------------------"
arq = open('Fonte_pontual_teste.vtk','w')

texto = '''# vtk DataFile Version 4.2
2D Laplace Python
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
    a,b = str(round(xn[i],3)),str(round(yn[i],3))
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
arq.write("\n")
arq.write('CELL_DATA ') #Vetor velocidade: no caso de solidos, v = 0
arq.write(str(elem))
arq.write('\n')
arq.write('VECTORS Gradiente float')
for d in range(len(G)):
    arq.write('\n')
    arq.write(str(G[d][0]))
    arq.write(' ')
    arq.write(str(G[d][1]))
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

arq.write("\n\n")
arq.write('POINT_DATA ') #Vetor velocidade: no caso de solidos, v = 0
arq.write(str(Npts))
arq.write("\n")
arq.write('SCALARS Temperatura double') #O nome Temperatura pode ser modificado.
arq.write('\n')
arq.write('LOOKUP_TABLE default')
for temp in range(len(T)):
    arq.write('\n')
    arq.write(str(T[temp]))

arq.close()
fim = time.time()
tempo = fim - inicio
print("Executado!")
print('\nNós:',Npts,' ','Elementos:',elem)
print('Tempo de execução: ',tempo,' segundos.')
