# FEM_python

# Delete Row: np.delete(matrix, row(wanna delete), axis = 0(row))
# Delete Column : np.delete(matrix, column(wana delete), axis = 1(column))

import numpy as np

# ne = number of element

df = 2 # degree of freedom
numel = 2 # number of element
nne = df * numel # number of degree of freedom of elements
TotNe = 5 # Number of total element

E = 30*(10**6)
L = 12*10
I = 50

# Moment and Concetrated Load

F = np.array([[0, 0, 10000, 0, 0, 0, 10000, 0, 0, 0]])
# F = np.zeros((1,df*TotNe))

v = np.zeros((df*TotNe, 1))

## 1. Local Stiffness Matrix of Beam(Concentrated Load)

def LocalBeamStiffnessMat(E, I, L, nne):
    Stiff = np.zeros((nne, nne))

    for i in range(0,nne):
        coe = np.array([-1**(i), -1**(i), -1**(i), -1**(i)])
        coe = coe.reshape(1,nne)
        if i == 0:
            Stiff[i,:] = [12, 6*L, -12, 6*L] 
        elif i == 2:
            Stiff[i,:] = coe * [12, 6*L, -12, 6*L] 
        elif i == 1:
            Stiff[i,:] = [6*L, 4*(L**2), -6*L, 2*(L**2)]
        else:
            Stiff[i,:] = [6*L, 2*(L**2), -6*L, 4*(L**2)]
            
    return  Stiff

## 2. Global Stiffness Matrix of Beam(Concentrated Load)

def GlobalBeamStiffnessMat(df, TotNe):
    
    TotNeDof = df * TotNe

    GlobalStiff = np.zeros((TotNeDof, TotNeDof))

    for i in range(0, TotNeDof+1):
        if i == nne:
            GlobalStiff[i-4:i, i-4:i] = GlobalStiff[i-4:i, i-4:i] + LocalBeamStiffnessMat(E, I, L, nne)
        for j in range(1, TotNe-1):
            if i == nne+2*j:
                GlobalStiff[i-4:i, i-4:i] = GlobalStiff[i-4:i, i-4:i] + LocalBeamStiffnessMat(E, I, L, nne)

    return  GlobalStiff * E*I/(L**3)

GlobalK = GlobalBeamStiffnessMat(df, TotNe)

print(GlobalK)

## 3. Solve Problem

### 3-1 : Get Degree of Freedom

def DofMatrix(*args):  # 가변인자를 통해 미지의 변수 통제

    Sol_K = GlobalBeamStiffnessMat(df, TotNe)
    lista = list(range(df*TotNe-1, -1, -1)) # list를 통해서 변위가 0인 지점 리스트화
    args = list(args) # 변위가 0이 아닌 지점

    for i in range(0, int(len(args))):
        lista.remove(args[i])

    for i in lista:
        Sol_K = np.delete(Sol_K, i, axis = 0)
        Sol_K = np.delete(Sol_K, i, axis = 1)

    return Sol_K

Sol_K = DofMatrix(2, 6)

## 3-2 : Get Displacemnet 

def LocalDisp(*args):

    LocalPower = np.zeros((1, int(len(args))))

    # args = list(args)

    for i in args:
        LocalPower = F[0, args]

    LocalPower = LocalPower.reshape(int(len(args)),1)

    Displacement = np.dot(np.linalg.inv(DofMatrix(*args)), LocalPower)

    Displacement = Displacement.tolist()

    for i in range(0, int(len(args))):
        v[args[i],0] = v[args[i], 0] + Displacement[i]
    
    return v

v = LocalDisp(2, 6)

print(v)

## 3-3 Get the power of point

sol = np.dot(GlobalK, v)

print(sol)

