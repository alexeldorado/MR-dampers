#-------------------------------------------------------------------------    
#1. Importando módulos
#------------------------------------------------------------------------- 

import numpy as np
import scipy.linalg as sc
import matplotlib.pyplot as plt
import pandas as pd
from scipy import linalg as la

#-------------------------------------------------------------------------    
#2. Declarando as funções
#------------------------------------------------------------------------- 

def importa_dados (xlsx): 
    # Recebe uma string do tipo 'nome_do_arquivo.xlsx'
    # Importa duas abas do arquivo .xls nomeadas como "nos" e "barras"
    # Retorna o conteúdo das abas em dois conjuntos de dados
    nos    = pd.read_excel(xlsx,'nos')
    barras = pd.read_excel(xlsx,'barras')    
    return nos, barras

def lista_gl (u,v,teta,nn):
    # Recebe as listas dos graus de liberdade u, v e teta e o número de nós (nn)
    # Concatena as listas e agrupa por nós
    # Retorna uma lista com os graus de Liberdade do sistema agrupados em função do seu nó
    a_gls = []
    for i in range (nn):
        gls_i  = [u[i],v[i],teta[i]]
        a_gls.append(gls_i)    
    return a_gls

def matriz_id_nos (noN,noF,nb):
    # Recebe as listas de Nós Iniciais (noN) e de Nós Finais (noF) e o número de barras (nb)
    # Retorna a Matriz de identificação dos nós( Nó inicial (N) e Nó Final (F) de cada barra)
    IDN      = np.zeros((2,nb))
    IDN[0,:] = noN
    IDN[1,:] = noF
    return IDN

def matriz_id_barras (IDN,nb):     
    # Recebe a Matriz id dos nós (IDN) e o número de barras (nb)
    # Retorna a Matriz de identificação das Barras em relação aos Graus de Liberdade
    IDB = np.zeros((6,nb)) #Graus de liberdade por barra, número de barras
    for i in range(3):
        IDB[i,:]   = IDN[0,:]*3-2+i
        IDB[i+3,:] = IDN[1,:]*3-2+i
    return IDB

def geometria (a, b, c, d):
    # Recebe a lista com a base(direção z) e a altura(x/y) das barras do pórtico
    # Retorna a lista com as áreas (A) e as Inércias (I) das barras
    Area = []
    I    = []
    for i in range (len(a)):
        Ai = a[i]*b[i] + c[i]*d[i]
        Area.append(Ai)        
        A1  = a[i]*b[i]
        A2  = c[i]*d[i]
        AT  = A1 + A2
        A1x = b[i]/2        
        A2x = (d[i]/2) + b[i]        
        xg  = (A1*A1x)/AT + (A2 *A2x)/AT        
        Ii = (((a[i]*(b[i]**3))/12)+(A1*((xg-A1x)**2))) + (((c[i]*(d[i]**3))/12)+(A2*((A2x-xg)**2)))
        I.append(Ii)   
    return Area, I

def L_e_coss(IDN,nb,cx,cy):
    # Recebe a Matriz id dos nós (IDN), o número de barras (nb) e as coordenadas Cx e Cy
    # Retorna a lista do Comprimento de cada barra e os cossenos diretores
    Lx   = np.zeros(nb)
    Ly   = np.zeros(nb)
    cosx = np.zeros(nb)
    cosy = np.zeros(nb)
    L    = np.zeros(nb)
    for n in range (nb):
        k1      = int(IDN[0,n] -1)  # Indexador da matriz IDN
        k2      = int(IDN[1,n] -1)  # Indexador da matriz IDN
        Lx[n]   = cx[k2] - cx[k1]
        Ly[n]   = cy[k2] - cy[k1]
        L[n]    = np.sqrt(Lx[n]**2 + Ly[n]**2)
        cosx[n] = Lx[n]/L[n]
        cosy[n] = Ly[n]/L[n]
    return L, cosx, cosy

def matrizes (ngl,nb,L,cosx,cosy,IDB,RHO,E,A,I):
    # Recebe:
    # ngl = Número de graus de liberdade
    # nb  = número de barras
    # IDN = Matriz id do Nó inicial (N) e Nó Final (F) de cada barra
    # IDB = Matriz id das Barras em relação aos Graus de Liberdade
    # cx  = lista com as coordenadas em X dos nós
    # cy  = lista com as coordenadas em Y dos nós
    # E   = Módulo de elasticidade (N/m2)
    # RHO = Massa específica (Kg/m³)
    # A   = Lista com as áreas (m2) das barras
    # I   = Lista com as Inércias (m4) das barras
    # Retorna as Matrizes de Rigidez (K) e de Massa (M)    
    
    # Criando as Matrizes 
    K    = np.zeros((ngl,ngl)) 
    M    = np.zeros((ngl,ngl))
    
    for i in range (nb):
    # Matriz de rigidez local da barra
        k = np.array([[E*A[i]/L[i], 0, 0, -E*A[i]/L[i],0 ,0 ],
                      [0, 12*E*I[i]/(L[i]**3), 6*E*I[i]/(L[i]**2), 0, -12*E*I[i]/(L[i]**3),6*E*I[i]/(L[i]**2)],
                      [0,6*E*I[i]/(L[i]**2), 4*E*I[i]/L[i], 0, -6*E*I[i]/(L[i]**2), 2*E*I[i]/L[i] ],
                      [-E*A[i]/L[i], 0, 0, E*A[i]/L[i],0 ,0 ],
                      [0, -12*E*I[i]/(L[i]**3), -6*E*I[i]/(L[i]**2), 0,12*E*I[i]/(L[i]**3),-6*E*I[i]/(L[i]**2)],
                      [0,6*E*I[i]/(L[i]**2), 2*E*I[i]/L[i], 0,-6*E*I[i]/(L[i]**2), 4*E*I[i]/L[i] ]])

    # Matriz de massa local da barra
        m = ((RHO*A[i]*L[i])/420)*np.array([[140, 0, 0, 70, 0, 0],
                                            [0, 156, 22*L[i], 0, 54, -13*L[i]],
                                            [0, 22*L[i], 4*(L[i]**2), 0, 13*L[i], -3*(L[i]**2)],
                                            [70, 0, 0, 140, 0, 0],
                                            [0, 54, 13*L[i], 0, 156, -22*L[i]],
                                            [0, -13*L[i], -3*(L[i]**2), 0, -22*L[i], 4*(L[i]**2)]])
    
    # Matriz de rotação
        tau = np.array([[cosx[i], cosy[i], 0, 0 ,0 ,0],
                        [-cosy[i], cosx[i],0, 0, 0, 0],
                        [0,0,1,0,0,0],                     
                        [0,0,0,cosx[i], cosy[i], 0],
                        [0, 0, 0,-cosy[i], cosx[i],0],
                        [0,0,0,0,0,1]])

    # Matrizes locais rotacionadas
        k_r = np.dot(np.dot(tau.T, k),tau)
        m_r = np.dot(np.dot(tau.T, m),tau)

    # Alocação das matrizes locais na matriz global
        k_rG = np.zeros((ngl,ngl))
        a1 = int(IDB[0,i]-1)
        a2 = int(IDB[2,i])
        a3 = int(IDB[3,i]-1)
        a4 = int(IDB[5,i])
        k_rG[a1:a2,a1:a2] = k_r[0:3,0:3]
        k_rG[a3:a4,a1:a2] = k_r[3:6,0:3]
        k_rG[a1:a2,a3:a4] = k_r[0:3,3:6]
        k_rG[a3:a4,a3:a4] = k_r[3:6,3:6]
        K += k_rG 

        m_rG = np.zeros((ngl,ngl))
        a1 = int(IDB[0,i]-1)
        a2 = int(IDB[2,i])
        a3 = int(IDB[3,i]-1)
        a4 = int(IDB[5,i])
        m_rG[a1:a2,a1:a2] = m_r[0:3,0:3]
        m_rG[a3:a4,a1:a2] = m_r[3:6,0:3]
        m_rG[a1:a2,a3:a4] = m_r[0:3,3:6]
        m_rG[a3:a4,a3:a4] = m_r[3:6,3:6]
        M += m_rG 
    return K,M

def remover_glr(K,M,u,v,teta,ur,vr,tetar):
    # Recebe:
    # K           = Matriz de rigidez
    # M           = Matriz de massa
    # u,v,teta    = Lista dos graus de liberdade
    # ur,vr,tetar = Lista dos graus de liberdade restringidos    
    # Retorna as Matrizes de Rigidez (Kf) e de Massa (Mf) sem os graus restringidos

    # Montar array com os Graus de Liberdade Restringidos 
    gl         = np.array(u+v+teta)
    id_glr     = np.array(ur+vr+tetar)
    glr        = np.trim_zeros(sorted(gl*id_glr))
    remover_gl = np.array(glr)-1

    # Deletar Linhas e Colunas restringidas das matrizes K e M
    Ki = np.delete(K, remover_gl,axis=0)
    Kf = np.delete(Ki, remover_gl,axis=1)  
    Mi = np.delete(M, remover_gl,axis=0)
    Mf = np.delete(Mi, remover_gl,axis=1)  
    return Kf, Mf

def EIG(Kf,Mf):
    # Recebe:
    # Kf  = Matriz de rigidez restringida
    # Mf  = Matriz de massa restringida
    # Soluciona o problema de Autovalores e Autovetores
    # Retorna as frequencias Naturais wk(rad/s) e fk(Hz)
    # Retorna a matriz com os modos de vibração (Phi)
    lamb,Phi  = sc.eig(Kf,Mf)
    index_eig = lamb.argsort()    # indexando em ordem crescente
    lamb      = lamb[index_eig]   # Aplicando indexador
    Phi       = Phi[:,index_eig]  # Aplicando indexador
    w2        = np.real(lamb)     # Extraindo apenas a parte real de lambda
    wk        = np.sqrt(w2)       # rad/s
    fk        = wk/2/np.pi        # Hz    
    return wk, fk, Phi

def amortecimento(zeta,wk,Mf,Kf):
    # Recebe:
    # zeta = Razão de amortecimento típica para estrutura em análise
    # wk   = Frequencias naturais(rad/s) do primeiro e segundo Modos
    # Kf   = Matriz de rigidez restringida
    # Mf   = Matriz de rigidez restringida
    # Retorna a Matriz de amortecimento de Rayleigh (Cf)    
    a  = zeta*2*wk[0]*wk[1]/(wk[0]+wk[1])  # parâmetro da matriz de massa
    b  = zeta*2/(wk[0]+wk[1])              # parâmetro da matriz de rigidez
    Cf = a*Mf + b*Kf
    return Cf

def indice_gl_horizontal (nn, nrestr):
    # Recebe o número de nós (nn) e a quantida de nós restringidos na base (nrestr)
    # Retorna uma Matriz com:
    # O numero do nó no desenho/ seu indice final (ja com nós restritos descontados) / Numero do gl horizontal
    
    espelhoNos = np.zeros ((nn,3))    
    
    for i in range(nrestr):
        espelhoNos [i,1] = 99999999 # indetificacao visual do no restrito
        espelhoNos [i,2] = 99999999 # indetificacao visual do no restrito     
    j = 1
    k = 0
    z = 0
    for i in range(len(espelhoNos)):
        espelhoNos [i,0] = j
        j += 1        
    for i in range(len(espelhoNos)-nrestr):    
        espelhoNos [i+nrestr,1] = k        
        espelhoNos [i+nrestr,2] = z       
        k += 1
        z += 3        
    return espelhoNos

def SIGMA (n,mMR,configuraMR,espelhoNos):
    # Recebe o numero de graus de liberdade ativos (n), o numero de amortecedores (mMR)
    # a matriz de Configuraçao (configuraMR), o espelho do problema (espelhoNos)
    # Retorna a Matriz Sigma de localizaçao das Forças
    Sigma = np.zeros((n,mMR))
    idGLMR = np.zeros((mMR,1))
    
    for m in range(mMR):    
        idGLMR [m,0] = int(espelhoNos[int(configuraMR[m,0])-1,2])
        Sigma[int(espelhoNos[int(configuraMR[m,0])-1,2]),m] = 1
        Sigma[int(espelhoNos[int(configuraMR[m,1])-1,2]),m] = -1
   
    return Sigma, idGLMR

def Ricattifunc (P,dynsys): 
    # Recebe:
    # P  = Variável desconhecida da Euqação de Ricatti
    # dynsys = Matrizes da estratégia LQR = matrizA, matrizB, matrizQ, matrizR
    # Retorna o resultado da Euqação de Ricatti, que deve ser zero.
    matrizP = P       
    matrizRinv = np.linalg.inv(matrizR)
    matriz_output = matrizP @ matrizAparaLQR - (1/2)*(matrizP @ matrizB @ matrizRinv @ matrizB.T @ matrizP) + matrizAparaLQR.T @ matrizP + 2*matrizQ
    return matriz_output 

def BoucwenNMR(x,t,Bwm,I): 
    # Recebe:
    # x   = Vetor estado (Ui,vi,zi,yi).T
    # t   = tempo atual (apenas para variações temporais de força) 
    # Bwm = Constantes do modelo de Boucwen modificado
    # I   = A corrente a ser aplicada nos MRs (A)
    # Retorna as derivadas de x

    alfa     = -826.67*I**3 + 905.14*I**2 + 412.52*I + 38.24
    c0       = -11.73*I**3 + 10.51*I**2 + 11.02*I + 0.59
    c1       = -54.40*I**3 + 57.03*I**2 + 64.57*I + 4.73     
    ya_dot1  = ((1/(c0+c1))*((alfa*x[2,0]) + (c0*x[1,0]) + (k0*(x[0,0]-x[3,0]))))
    z_dot1   = -beta*abs(x[1,0]-ya_dot1)*x[2,0]*(abs(x[2,0])**(nbw-1)) - gama*(x[1,0]-ya_dot1)*(abs(x[2,0])**nbw) + A*(x[1,0]-ya_dot1)    
    saida       = np.zeros((4,1)) 
    saida [0,0] = x[0,0]
    saida [1,0] = x[1,0]
    saida [2,0] = z_dot1
    saida [3,0] = ya_dot1   

    return saida

#-------------------------------------------------------------------------    
#3. Dados de Entrada  
#------------------------------------------------------------------------- 
nos, barras = importa_dados('dados_de_entradaEdificioOriginal.xlsx')
RHO         = 2500                           # Massa específica (Kg/m³)
fck         = 50                             # mPa
E           = 5600*((fck)**(1/2))*0.85*10**6 # Módulo de elasticidade (N/m2)
zeta        = 0.02                           # Razão de amortecimento típica para estrutura aporticada de concreto

#-------------------------------------------------------------------------     
#4. Matrizes auxiliares 
#------------------------------------------------------------------------- 

n_id     = list(nos['Nó'])           # Identificação dos Nós
nn       = len(list(nos['Cx']))      # número de nós
ngl      = len(list(nos['Cx']))*3    # número de graus de liberdade 
ug       = list(nos['u'])[0:nn]      # Identificação do Grau de Liberdade
vg       = list(nos['v'])[0:nn]      # Identificação do Grau de Liberdade
tetag    = list(nos['teta'])[0:nn]   # Identificação do Grau de Liberdade
ur       = list(nos['ur'])[0:nn]     # Identificação dos Graus de Liberdade restringidos
vr       = list(nos['vr'])[0:nn]     # Identificação dos Graus de Liberdade restringidos
tetar    = list(nos['tetar'])[0:nn]  # Identificação dos Graus de Liberdade restringidos
cx       = list(nos['Cx'])[0:nn]     # coordenadas em X dos nós
cy       = list(nos['Cy'])[0:nn]     # coordenadas em Y dos nós
cz       = list(nos['Cz'])[0:nn]     # coordenadas em Z dos nós
noN      = list(barras['noN'])       # Nós Iniciais
noF      = list(barras['noF'])       # Nós Finais
basea    = list(barras['a'])         # Base dos elementos do pórtico (direção z)
alturab  = list(barras['b'])         # Altura dos elementos do pórtico (direção x/y)
basec    = list(barras['c'])         # base segunda parte do pilar em L
alturad  = list(barras['d'])         # Altura segunda parte do pilar em L
nb       = len(noN)                  # Número de Barras
a_gls    = lista_gl(ug,vg,tetag,nn)  # Graus de Liberdade agrupados por nó
IDN      = matriz_id_nos(noN,noF,nb) # Matriz id do Nó inicial (N) e Nó Final (F) de cada barra
IDB      = matriz_id_barras(IDN,nb)  # Matriz id das Barras em relação aos Graus de Liberdade

# Aplicando as dimensoes otimizadas nos elementos 

# Importar o Vetor de Dimensos otimas

x0 = DimOTM

# Localizadores dos grupos

# Pilares Laterais ; Obs.: basec e alturab são fixos e valem zero
# Grupo 1:
# x0 = Lado = baseaa = alturab
G1barras = np.zeros((10,1))
G1barra0 = 1 
qtdGrupo = 0
while qtdGrupo < 10:
    G1barras[qtdGrupo,0] = G1barra0 -1      #Barra indexada
    G1barras[qtdGrupo+1,0] = G1barra0 +3 -1 #Barra indexada
    G1barra0 += 7 
    qtdGrupo += 2
# Grupo 3:
# x2 = Lado = baseaa = alturab
G3barras = np.zeros((10,1))
G3barra0 = 36 
qtdGrupo = 0
while qtdGrupo < 10:
    G3barras[qtdGrupo,0] = G3barra0 -1      #Barra indexada
    G3barras[qtdGrupo+1,0] = G3barra0 +3 -1 #Barra indexada
    G3barra0 += 7 
    qtdGrupo += 2
# Grupo 5:
# x4 = Lado = baseaa = alturab
G5barras = np.zeros((10,1))
G5barra0 = 71 
qtdGrupo = 0
while qtdGrupo < 10:
    G5barras[qtdGrupo,0] = G5barra0 -1      #Barra indexada
    G5barras[qtdGrupo+1,0] = G5barra0 +3 -1 #Barra indexada
    G5barra0 += 7 
    qtdGrupo += 2    
# Grupo 7:
# x6 = Lado = baseaa = alturab
G7barras = np.zeros((10,1))
G7barra0 = 106
qtdGrupo = 0
while qtdGrupo < 10:
    G7barras[qtdGrupo,0] = G7barra0 -1      #Barra indexada
    G7barras[qtdGrupo+1,0] = G7barra0 +3 -1 #Barra indexada
    G7barra0 += 7 
    qtdGrupo += 2  
# Grupo 9:
# x8 = Lado = baseaa = alturab
G9barras = np.zeros((10,1))
G9barra0 = 141 
qtdGrupo = 0
while qtdGrupo < 10:
    G9barras[qtdGrupo,0] = G9barra0 -1      #Barra indexada
    G9barras[qtdGrupo+1,0] = G9barra0 +3 -1 #Barra indexada
    G9barra0 += 7 
    qtdGrupo += 2
# Grupo 11:
# x10 = Lado = baseaa = alturab
G11barras = np.zeros((10,1))
G11barra0 = 176 
qtdGrupo = 0
while qtdGrupo < 10:
    G11barras[qtdGrupo,0] = G11barra0 -1      #Barra indexada
    G11barras[qtdGrupo+1,0] = G11barra0 +3 -1 #Barra indexada
    G11barra0 += 7 
    qtdGrupo += 2
# Grupo 13:
# x12 = Lado = baseaa = alturab
G13barras = np.zeros((10,1))
G13barra0 = 211 
qtdGrupo = 0
while qtdGrupo < 10:
    G13barras[qtdGrupo,0] = G13barra0 -1      #Barra indexada
    G13barras[qtdGrupo+1,0] = G13barra0 +3 -1 #Barra indexada
    G13barra0 += 7 
    qtdGrupo += 2    

# Pilares Centro ; Obs.: alturab e basec são fixos e valem 0,32
# Grupo 2:
# x1 = D = alturad ; basea = 32 + D
G2barras = np.zeros((10,1))
G2barra0 = 2 
qtdGrupo = 0
while qtdGrupo < 10:
    G2barras[qtdGrupo,0] = G2barra0 -1      #Barra indexada
    G2barras[qtdGrupo+1,0] = G2barra0 +1 -1 #Barra indexada
    G2barra0 += 7 
    qtdGrupo += 2
# Grupo 4:
# x3 = D = alturad ; basea = 32 + D
G4barras = np.zeros((10,1))
G4barra0 = 37 
qtdGrupo = 0
while qtdGrupo < 10:
    G4barras[qtdGrupo,0] = G4barra0 -1      #Barra indexada
    G4barras[qtdGrupo+1,0] = G4barra0 +1 -1 #Barra indexada
    G4barra0 += 7 
    qtdGrupo += 2
# Grupo 6:
# x5 = D = alturad ; basea = 32 + D
G6barras = np.zeros((10,1))
G6barra0 = 72 
qtdGrupo = 0
while qtdGrupo < 10:
    G6barras[qtdGrupo,0] = G6barra0 -1      #Barra indexada
    G6barras[qtdGrupo+1,0] = G6barra0 +1 -1 #Barra indexada
    G6barra0 += 7 
    qtdGrupo += 2
# Grupo 8:
# x7 = D = alturad ; basea = 32 + D
G8barras = np.zeros((10,1))
G8barra0 = 107 
qtdGrupo = 0
while qtdGrupo < 10:
    G8barras[qtdGrupo,0] = G8barra0 -1      #Barra indexada
    G8barras[qtdGrupo+1,0] = G8barra0 +1 -1 #Barra indexada
    G8barra0 += 7 
    qtdGrupo += 2
# Grupo 10:
# x9 = D = alturad ; basea = 32 + D
G10barras = np.zeros((10,1))
G10barra0 = 142 
qtdGrupo = 0
while qtdGrupo < 10:
    G10barras[qtdGrupo,0] = G10barra0 -1      #Barra indexada
    G10barras[qtdGrupo+1,0] = G10barra0 +1 -1 #Barra indexada
    G10barra0 += 7 
    qtdGrupo += 2
# Grupo 12:
# x11 = D = alturad ; basea = 32 + D
G12barras = np.zeros((10,1))
G12barra0 = 177 
qtdGrupo = 0
while qtdGrupo < 10:
    G12barras[qtdGrupo,0] = G12barra0 -1      #Barra indexada
    G12barras[qtdGrupo+1,0] = G12barra0 +1 -1 #Barra indexada
    G12barra0 += 7 
    qtdGrupo += 2
# Grupo 14:
# x13 = D = alturad ; basea = 32 + D
G14barras = np.zeros((10,1))
G14barra0 = 212 
qtdGrupo = 0
while qtdGrupo < 10:
    G14barras[qtdGrupo,0] = G14barra0 -1      #Barra indexada
    G14barras[qtdGrupo+1,0] = G14barra0 +1 -1 #Barra indexada
    G14barra0 += 7 
    qtdGrupo += 2

# Vigas Laterais/centrais 
# Grupo 15:
# x14 = basea
# x15 = alturab
G15barras = np.zeros((10,1))
G15barra0 = 5 
qtdGrupo = 0
while qtdGrupo < 10:
    G15barras[qtdGrupo,0] = G15barra0 -1      #Barra indexada
    G15barras[qtdGrupo+1,0] = G15barra0 +2 -1 #Barra indexada
    G15barra0 += 7 
    qtdGrupo += 2
# Grupo 16:
# x16 = basea
# x17 = alturab
G16barras = np.zeros((5,1))
G16barra0 = 6
qtdGrupo = 0
while qtdGrupo < 5:
    G16barras[qtdGrupo,0] = G16barra0 -1      #Barra indexada    
    G16barra0 += 7 
    qtdGrupo += 1
# Grupo 17:
# x18 = basea
# x19 = alturab
G17barras = np.zeros((10,1))
G17barra0 = 40 
qtdGrupo = 0
while qtdGrupo < 10:
    G17barras[qtdGrupo,0] = G17barra0 -1      #Barra indexada
    G17barras[qtdGrupo+1,0] = G17barra0 +2 -1 #Barra indexada
    G17barra0 += 7 
    qtdGrupo += 2
# Grupo 18:
# x20 = basea
# x21 = alturab
G18barras = np.zeros((5,1))
G18barra0 = 41
qtdGrupo = 0
while qtdGrupo < 5:
    G18barras[qtdGrupo,0] = G18barra0 -1      #Barra indexada    
    G18barra0 += 7 
    qtdGrupo += 1
# Grupo 19:
# x22 = basea
# x23 = alturab
G19barras = np.zeros((10,1))
G19barra0 = 75 
qtdGrupo = 0
while qtdGrupo < 10:
    G19barras[qtdGrupo,0] = G19barra0 -1      #Barra indexada
    G19barras[qtdGrupo+1,0] = G19barra0 +2 -1 #Barra indexada
    G19barra0 += 7 
    qtdGrupo += 2
# Grupo 20:
# x24 = basea
# x25 = alturab
G20barras = np.zeros((5,1))
G20barra0 = 76
qtdGrupo = 0
while qtdGrupo < 5:
    G20barras[qtdGrupo,0] = G20barra0 -1      #Barra indexada    
    G20barra0 += 7 
    qtdGrupo += 1
# Grupo 21:
# x26 = basea
# x27 = alturab
G21barras = np.zeros((10,1))
G21barra0 = 110 
qtdGrupo = 0
while qtdGrupo < 10:
    G21barras[qtdGrupo,0] = G21barra0 -1      #Barra indexada
    G21barras[qtdGrupo+1,0] = G21barra0 +2 -1 #Barra indexada
    G21barra0 += 7 
    qtdGrupo += 2
# Grupo 22:
# x28 = basea
# x29 = alturab
G22barras = np.zeros((5,1))
G22barra0 = 111
qtdGrupo = 0
while qtdGrupo < 5:
    G22barras[qtdGrupo,0] = G22barra0 -1      #Barra indexada    
    G22barra0 += 7 
    qtdGrupo += 1
# Grupo 23:
# x30 = basea
# x31 = alturab
G23barras = np.zeros((10,1))
G23barra0 = 145 
qtdGrupo = 0
while qtdGrupo < 10:
    G23barras[qtdGrupo,0] = G23barra0 -1      #Barra indexada
    G23barras[qtdGrupo+1,0] = G23barra0 +2 -1 #Barra indexada
    G23barra0 += 7 
    qtdGrupo += 2
# Grupo 24:
# x32 = basea
# x33 = alturab
G24barras = np.zeros((5,1))
G24barra0 = 146
qtdGrupo = 0
while qtdGrupo < 5:
    G24barras[qtdGrupo,0] = G24barra0 -1      #Barra indexada    
    G24barra0 += 7 
    qtdGrupo += 1
# Grupo 25:
# x34 = basea
# x35 = alturab
G25barras = np.zeros((10,1))
G25barra0 = 180 
qtdGrupo = 0
while qtdGrupo < 10:
    G25barras[qtdGrupo,0] = G25barra0 -1      #Barra indexada
    G25barras[qtdGrupo+1,0] = G25barra0 +2 -1 #Barra indexada
    G25barra0 += 7 
    qtdGrupo += 2
# Grupo 26:
# x36 = basea
# x37 = alturab
G26barras = np.zeros((5,1))
G26barra0 = 181
qtdGrupo = 0
while qtdGrupo < 5:
    G26barras[qtdGrupo,0] = G26barra0 -1      #Barra indexada    
    G26barra0 += 7 
    qtdGrupo += 1
# Grupo 27:
# x38 = basea
# x39 = alturab
G27barras = np.zeros((10,1))
G27barra0 = 215 
qtdGrupo = 0
while qtdGrupo < 10:
    G27barras[qtdGrupo,0] = G27barra0 -1      #Barra indexada
    G27barras[qtdGrupo+1,0] = G27barra0 +2 -1 #Barra indexada
    G27barra0 += 7 
    qtdGrupo += 2
# Grupo 28:
# x40 = basea
# x41 = alturab
G28barras = np.zeros((5,1))
G28barra0 = 216
qtdGrupo = 0
while qtdGrupo < 5:
    G28barras[qtdGrupo,0] = G28barra0 -1      #Barra indexada    
    G28barra0 += 7 
    qtdGrupo += 1

#Pilares Laterais
for i in range(len (G1barras)):
    basea   [int(G1barras[i])]  = x0[0]
    alturab [int(G1barras[i])]  = x0[0] 
for i in range(len (G3barras)):
    basea   [int(G3barras[i])]  = x0[2]
    alturab [int(G3barras[i])]  = x0[2] 
for i in range(len (G5barras)):
    basea   [int(G5barras[i])]  = x0[4]
    alturab [int(G5barras[i])]  = x0[4]     
for i in range(len (G7barras)):
    basea   [int(G7barras[i])]  = x0[6]
    alturab [int(G7barras[i])]  = x0[6] 
for i in range(len (G9barras)):
    basea   [int(G9barras[i])]  = x0[8]
    alturab [int(G9barras[i])]  = x0[8] 
for i in range(len (G11barras)):
    basea   [int(G11barras[i])]  = x0[10]
    alturab [int(G11barras[i])]  = x0[10] 
for i in range(len (G13barras)):
    basea   [int(G13barras[i])]  = x0[12]
    alturab [int(G13barras[i])]  = x0[12]   
#Pilares centro
for i in range(len (G2barras)):
    basea   [int(G2barras[i])]  = x0[1]  + 0.32
    alturad [int(G2barras[i])]  = x0[1] 
for i in range(len (G4barras)):
    basea   [int(G4barras[i])]  = x0[3] + 0.32
    alturad [int(G4barras[i])]  = x0[3] 
for i in range(len (G6barras)):
    basea   [int(G6barras[i])]  = x0[5] + 0.32
    alturad [int(G6barras[i])]  = x0[5] 
for i in range(len (G8barras)):
    basea   [int(G8barras[i])]  = x0[7] + 0.32
    alturad [int(G8barras[i])]  = x0[7] 
for i in range(len (G10barras)):
    basea   [int(G10barras[i])]  = x0[9] + 0.32
    alturad [int(G10barras[i])]  = x0[9] 
for i in range(len (G12barras)):
    basea   [int(G12barras[i])]  = x0[11] + 0.32
    alturad [int(G12barras[i])]  = x0[11] 
for i in range(len (G14barras)):
    basea   [int(G14barras[i])]  = x0[13] + 0.32
    alturad [int(G14barras[i])]  = x0[13]     
#Vigas
for i in range(len (G15barras)):
    basea   [int(G15barras[i])]  = x0[14]
    alturab [int(G15barras[i])]  = x0[15]
for i in range(len (G16barras)):
    basea   [int(G16barras[i])]  = x0[16]
    alturab [int(G16barras[i])]  = x0[17]   
for i in range(len (G17barras)):
    basea   [int(G17barras[i])]  = x0[18]
    alturab [int(G17barras[i])]  = x0[19]   
for i in range(len (G18barras)):
    basea   [int(G18barras[i])]  = x0[20]
    alturab [int(G18barras[i])]  = x0[21]   
for i in range(len (G19barras)):
    basea   [int(G19barras[i])]  = x0[22]
    alturab [int(G19barras[i])]  = x0[23]  
for i in range(len (G20barras)):
    basea   [int(G20barras[i])]  = x0[24]
    alturab [int(G20barras[i])]  = x0[25]   
for i in range(len (G21barras)):
    basea   [int(G21barras[i])]  = x0[26]
    alturab [int(G21barras[i])]  = x0[27]   
for i in range(len (G22barras)):
    basea   [int(G22barras[i])]  = x0[28]
    alturab [int(G22barras[i])]  = x0[29]  
for i in range(len (G23barras)):
    basea   [int(G23barras[i])]  = x0[30]
    alturab [int(G23barras[i])]  = x0[31]    
for i in range(len (G24barras)):
    basea   [int(G24barras[i])]  = x0[32]
    alturab [int(G24barras[i])]  = x0[33]   
for i in range(len (G25barras)):
    basea   [int(G25barras[i])]  = x0[34]
    alturab [int(G25barras[i])]  = x0[35]   
for i in range(len (G26barras)):
    basea   [int(G26barras[i])]  = x0[36]
    alturab [int(G26barras[i])]  = x0[37]   
for i in range(len (G27barras)):
    basea   [int(G27barras[i])]  = x0[38]
    alturab [int(G27barras[i])]  = x0[39]   
for i in range(len (G28barras)):
    basea   [int(G28barras[i])]  = x0[40]
    alturab [int(G28barras[i])]  = x0[41]

Area, I  = geometria(basea,alturab,basec,alturad)        # Lista com as áreas (m2) e as Inércias (m4) das barras

# Comprimento de cada barra e cossenos diretores
L, cosx, cosy = L_e_coss(IDN,nb,cx,cy)

#-------------------------------------------------------------------------        
#5. Matrizes de Rigidez e Massa
#-------------------------------------------------------------------------  

#5.1 Definindo as Matrizes de Rigidez (K) e Massa (M) Iniciais do pórtico
K,M   = matrizes(ngl,nb,L,cosx,cosy,IDB,RHO,E,Area,I)

#5.2 Removendo os graus de liberdade restringidos
K,M = remover_glr(K,M,ug,vg,tetag,ur,vr,tetar)

#-------------------------------------------------------------------------     
#6. Matriz de Massa das Lajes (Lumped)
#------------------------------------------------------------------------- 

#6.1 Definindo a Matriz de Massa Lumped
ML = np.zeros((len(M),len(M))) 

nostemp       = list(np.arange(1,(nn+1),1))
nos_externos  = sorted(list(np.arange(4,(nn+4),4)) + list(np.arange(1,(nn+1),4)))
relacao_gl_no = list(np.arange(0,(2*nn),2))  #de zero até o dobro de nós
gl_externos   = []

for i in range(len(nostemp)+1):
    if i in nos_externos:
        utemp = nostemp[i-1]+relacao_gl_no[i-1]
        vtemp = utemp + 1
        ttemp = utemp + 2
        gl_externos.append(utemp)
        gl_externos.append(vtemp)
        gl_externos.append(ttemp)

for i in range(len(ML)):
     for j in range(len(ML)):
        if i == j:
            if i+1 in gl_externos:
                ML[i,j] = 3526.215  #Kg
            else:
                ML[i,j] = 5356.4175 #Kg


#6.2 Somando a Massa das Lajes à Matriz de Massa (Consistente + Lumped)
M  += ML

#-------------------------------------------------------------------------        
#7. Frequências Naturais e Modos de Vibração
#------------------------------------------------------------------------- 
# Frequencias Naturais wk(rad/s) e fk(Hz) e Modos de vibração (Phi)
wk, fk, Phi = EIG(K,M)

#-------------------------------------------------------------------------        
#8. Matriz de Amortecimento
#------------------------------------------------------------------------- 
C = amortecimento(zeta,wk,M,K)

#-------------------------------------------------------------------------        
#9. Propriedades dinâmicas
#------------------------------------------------------------------------- 

#Identificando os graus de liberdade a serem carregados
espelhoNos = indice_gl_horizontal(nn,4) #Matriz de Localizacao
ncarregados = list(np.arange(5,142,4))  #Nós carregados
glcarregados = []                       #Lista vazia para gl carregados
for i in range(len(ncarregados)):
    temp  = (ncarregados[i]*3)-2        #regra de identificação do gl
    temp -= 13                          #Indexização do gl, removendo gl restringidos
    glcarregados.append(temp)

n   = ngl - 12                          # Número de Graus de Liberdade
mFe = len(glcarregados)                 # Número de Forças externas
mMR = 34                                 # Número de atuadores MR
Id  = np.eye(n)
Id2 = np.eye(mMR)
matrizR = Id2*10**-7

Minv = np.linalg.inv(M)

matrizAparaLQR = np.zeros((2*n,2*n))    
A12 = np.eye((n))
A21 = -1*Minv @ K 
A22 = -1*Minv @ C 
matrizAparaLQR [0:n,n:] = A12
matrizAparaLQR [n:,0:n] = A21
matrizAparaLQR [n:,n:]  = A22

configuraMR = np.zeros((mMR,2)) #Matriz mMR x 2 , com os pares: Nó controlado, Nó de fixacao, conforme desenho.
configuraMR[0,0], configuraMR[1,1] = 10,7
configuraMR[1,0], configuraMR[1,1] = 14,11
configuraMR[2,0], configuraMR[2,1] = 18,15
configuraMR[3,0], configuraMR[3,1] = 22,19
configuraMR[4,0], configuraMR[4,1] = 26,23
configuraMR[5,0], configuraMR[5,1] = 30,27
configuraMR[6,0], configuraMR[6,1] = 34,31
configuraMR[7,0], configuraMR[7,1] = 38,35
configuraMR[8,0], configuraMR[8,1] = 42,39
configuraMR[9,0], configuraMR[9,1] = 46,43
configuraMR[10,0], configuraMR[10,1] = 50,47
configuraMR[11,0], configuraMR[11,1] = 54,51
configuraMR[12,0], configuraMR[12,1] = 58,55
configuraMR[13,0], configuraMR[13,1] = 62,59
configuraMR[14,0], configuraMR[14,1] = 66,63
configuraMR[15,0], configuraMR[15,1] = 70,67
configuraMR[16,0], configuraMR[16,1] = 74,71
configuraMR[17,0], configuraMR[17,1] = 78,75
configuraMR[18,0], configuraMR[18,1] = 82,79
configuraMR[19,0], configuraMR[19,1] = 86,83
configuraMR[20,0], configuraMR[20,1] = 90,87
configuraMR[21,0], configuraMR[21,1] = 94,91
configuraMR[22,0], configuraMR[22,1] = 98,95
configuraMR[23,0], configuraMR[23,1] = 102,99
configuraMR[24,0], configuraMR[24,1] = 106,103
configuraMR[25,0], configuraMR[25,1] = 110,107
configuraMR[26,0], configuraMR[26,1] = 114,111
configuraMR[27,0], configuraMR[27,1] = 118,115
configuraMR[28,0], configuraMR[28,1] = 122,119
configuraMR[29,0], configuraMR[29,1] = 126,123
configuraMR[30,0], configuraMR[30,1] = 130,127
configuraMR[31,0], configuraMR[31,1] = 134,131
configuraMR[32,0], configuraMR[32,1] = 138,135
configuraMR[33,0], configuraMR[33,1] = 142,139

matrizB = np.zeros((2*n,mMR))
Sigma,idGLMR   = SIGMA(n,mMR,configuraMR,espelhoNos)
B21     = Minv @ Sigma
matrizB[n:,0:] = B21

matrizQ = np.zeros((2*n,2*n))
matrizQ [0:n,0:n] = K 

dynsys =  {"M":M,"K":K,"C":C,"Id":Id,"n":n,"Minv":Minv,"mFe":mFe,"mMR":mMR,"matrizAparaLQR":matrizAparaLQR, "matrizB":matrizB,"matrizQ":matrizQ, "matrizR":matrizR}

# Bouc Wen modificado
A    = 10.013
# beta = 3.044
# gama = 0.103

gama = 3.044
beta = 0.103

k0   = 1.121
f0   = 40
nbw  = 2
Bwm  =  {"A":A,"beta":beta,"gama":gama,"k0":k0,"f0":f0,"nbw":nbw}

#-------------------------------------------------------------------------        
#10. Força do vento
#------------------------------------------------------------------------- 
# Observações:
# Importar para o spyder a Matriz de Força gerada em outro arquivo

F = F0

#-------------------------------------------------------------------------        
#11. Estratégia de Controle
#-------------------------------------------------------------------------
# Resolvendo a Equação de Ricatti
P = la.solve_continuous_are(matrizAparaLQR, matrizB, 2*matrizQ, 2*matrizR)

# Montando o vetor Ganho para solução ótima LQR

matrizRinv = np.linalg.inv(matrizR)
Ganho      = (-1/2) * matrizRinv @ matrizB.T @ P 

# teste do ganho
deveserzero = Ricattifunc (P,dynsys)

#-------------------------------------------------------------------------        
#12. Resolvendo a Equação do Equilíbrio Dinâmico
#------------------------------------------------------------------------- 
# Montar arrays Necessários
dur     = 300                    # Duração da análise (s)
dt      = 0.0005                 # Passo de tempo (s)
tf      = int(dur/dt)            # Tamanho do vetor de tempo
t       = np.linspace(0,dur,tf)  # Vetor de tempo discretizado
n       = len(F[:,0])            # Graus de liberdade
Acc     = np.zeros((n,tf))       # Aceleração (m/s²)
v       = np.zeros((n,tf))       # Velocidade (m/s)
U       = np.zeros((n,tf))       # Deslocamento (m) 
Fc      = np.zeros((n,tf))       # Vetor Força de amortecimento (N)
Fcvec   = np.zeros((mMR,tf))     # Vetor temporário das Forças de amortecimento (N)
Ivec    = np.zeros((mMR,tf))     # Vetor de Correntes (A)
z       = np.zeros((mMR,tf))     # Vetor da variável histerética
y       = np.zeros((mMR,tf))     # Vetor da variável modificada
Fco     = np.zeros((mMR,tf))     # Vetor de Forças ótimas (N)
estado  = np.zeros((2*n,tf))     # Vetor de estado 
Imax    = 0.5                    # Corrente máxima nos amortecedores

# Determinar as constantes do método de Newmark
delta    = 0.5
alfa1    = 0.25
a10      = 1/(alfa1*(dt**2))
a11      = 1/(alfa1*dt)
a12      = (1/(2*alfa1))-1
a13      = delta/(dt*alfa1)
a14      = (delta/alfa1) - 1
a15      = (dt/2)*((delta/alfa1) - 2)
C1       = np.linalg.inv(a10*M + a13*C + K)
Acc[:,0] = np.dot(np.linalg.inv(M),(F[:,0]-np.dot(C,v[:,0])-np.dot(K,U[:,0]))) #aceleração no tempo zero

# Looping principal
for i in range(tf-1): 

    
    # Seleção da Voltagem    
    estado [0:n,i]   = U[:,i] 
    estado [n:2*n,i] = v[:,i]               
    Fco[:,i]=  Ganho @ estado [:,i] #Força ótima (N)

    for m in range(len(Fcvec)):        
        Ivec[m,i+1] = Imax*np.heaviside((Fco[m,i]-(1*Fcvec[m,i]))*(1*Fcvec[m,i]),(Ivec[m,i]/Imax)) 
    
    #RUNGE KUTTA 4    
    
    contaBw = 0    
    for m in range(len(Fcvec)):
        contaGL  = int(idGLMR[m,0])
        xi       = np.zeros((4,1))    
        xi [0,0] = U[contaGL,i]*1000
        xi [1,0] = v[contaGL,i]*1000
        xi [2,0] = z[contaBw,i]
        xi [3,0] = y[contaBw,i]
                
        k1 = dt*BoucwenNMR(xi,t[i], Bwm, Ivec[contaBw,i+1])
        k2 = dt*BoucwenNMR(xi + k1/2, t[i] + dt/2, Bwm, Ivec[contaBw,i+1])
        k3 = dt*BoucwenNMR(xi + k2/2, t[i] + dt/2, Bwm, Ivec[contaBw,i+1])
        k4 = dt*BoucwenNMR(xi + k3, t[i], Bwm, Ivec[contaBw,i+1])
            
        dx = (k1 + 2*k2 + 2*k3 + k4)/6
            
        xi = xi + dx  
        
        z[contaBw,i+1] = xi [2,0]
        y[contaBw,i+1] = xi [3,0]     
        
                
        alfa     = -826.67*Ivec[contaBw,i+1]**3 + 905.14*Ivec[contaBw,i+1]**2 + 412.52*Ivec[contaBw,i+1] + 38.24
        c0       = -11.73*Ivec[contaBw,i+1]**3 + 10.51*Ivec[contaBw,i+1]**2 + 11.02*Ivec[contaBw,i+1] + 0.59
        c1       = -54.40*Ivec[contaBw,i+1]**3 + 57.03*Ivec[contaBw,i+1]**2 + 64.57*Ivec[contaBw,i+1] + 4.73
        ya_dot1  = ((1/(c0+c1))*((alfa*z[contaBw,i+1]) + (c0*v[contaGL,i]*1000) + (k0*(U[contaGL,i]*1000-y[contaBw,i+1]))))
        fcm      = (c1*ya_dot1 + f0)                      
        Fcvec [m,i+1] = -1*fcm*20 # A troca de sinal da EED está aqui                               
        contaBw += 1    

    # Newmark    
    Fc[:,i+1] = Sigma @ Fcvec [:,i+1]  

    var1     = F[:,i+1]+Fc[:,i+1]+np.dot(M,(a10*U[:,i]+ a11*v[:,i] + a12*Acc[:,i]))+np.dot(C,(a13*U[:,i]+ a14*v[:,i] 
                                                                                    + a15*Acc[:,i]))
    U[:,i+1]   = np.dot(C1,var1)
    v[:,i+1]   = a13*(U[:,i+1] - U[:,i]) - a14*v[:,i] - a15*Acc[:,i]    
    Acc[:,i+1] = a10*(U[:,i+1] - U[:,i]) - a11*v[:,i] - a12*Acc[:,i] 
