import numpy as np
from matplotlib import pyplot as plt

def plot(cdnt, dk, result, real):
    x = [i[0] for i in cdnt]
    y = [i[1] for i in cdnt]
    realX, realY = real
    resultX, resultY = result
    c1 = plt.Circle((x[0], y[0]), dk[0], fill = False)
    c2 = plt.Circle((x[1], y[1]), dk[1], fill = False)
    c3 = plt.Circle((x[2], y[2]), dk[2], fill = False)
    c4 = plt.Circle((x[3], y[3]), dk[3], fill = False)
    c5 = plt.Circle((x[4], y[4]), dk[4], fill = False)
    fig, ax = plt.subplots()   
    ax.add_patch(c1)
    ax.add_patch(c2)
    ax.add_patch(c3)
    ax.add_patch(c4)
    ax.add_patch(c5)
    # Receptores
    plt.scatter(x, y, c = "blue")
    # Emissor    
    plt.scatter(float(resultX), float(resultY), c = "brown")
    # Real
    plt.scatter(float(realX), (float(realY)), c = "green")
    plt.show()  

def calcDk(pk0, lk, pk):
    result = []    
    for i in range(len(pk0)):
        result.append(10**((pk0[i]-pk[i])/(10*lk[i])))
    return result

def euclidianDistance(real, calculated):
    real = np.array(real)
    calculated = np.array(calculated)
    return np.linalg.norm(real-calculated)

def linearCombination(cdnt, dk): 
    A = np.empty((0, 2))
    B = np.empty((0, 2))

    x = [i[0] for i in cdnt]
    y = [i[1] for i in cdnt]
    # z = [i[2] for i in cdnt]

    for i in range(1, len(x)):
        A = np.append(A, [[2*(x[i]-x[0]), 2*(y[i]-y[0])]], axis=0)
    #print(A)
        
    for i in range(1, len(x)):
        B = np.append(B, [[((dk[i]**2 - (x[i]**2 + y[i]**2)) - ((dk[0]**2 - (x[0]**2 + y[0]**2))))]])
    #print(B)
    return np.dot(np.dot(np.linalg.inv(np.dot(A.transpose(), A)), A.transpose()), B)

# coordenadas conhecidas do receptor K
coordinates = [[1.55, 17.63, 1.35],
               [-4.02, 0.00, 1.35],
               [-4.40, 9.60, 1.35],
               [9.27, 4.64, 1.35],
               [9.15, 12.00, 1.35]]

pk0 = [-26.0, -33.8, -29.8, -31.2, -33.0]
lk = [2.1, 1.8, 1.3, 1.4, 1.5]

# caso 1
#pk = [-48.4, -50.6, -32.2, -47.4, -46.3]

# caso 2
pk = [-46.9, -46.4, -41.2, -45.8, -48.7]

dk = calcDk(pk0, lk, pk)
real = [3, 3]
result = linearCombination(coordinates, dk)
print("Real [x y]: {}".format(real))
print("Resultado [x y]: {}".format(result))
print("Distância euclidiana: {}".format(euclidianDistance(real, result)))
plot(coordinates, dk, result, real)
