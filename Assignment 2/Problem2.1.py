import numpy as np
import matplotlib.pyplot as plt


# Parameters
Length = 100
Time = 100
rho = 0.5
q = 8
dt = 1 / (rho * 100)

# Fixed points
u1 = 0
u2 = 1/2 * ( -1 + q + (np.sqrt(-4 * q  + rho * (1 + q)**2) / np.sqrt(rho))) # Largest fixed point 5.5615
u3 = 1/2 * ( -1 + q - (np.sqrt(-4 * q  + rho * (1 + q)**2) / np.sqrt(rho))) # 1.4384

# Initial conditions
x0 = 20
u0 = u2

# Functions for the problem

def u(x): # This function is used to initialize the matrixSpaceTime
    global u0
    global x0
    return u0 / (1+np.exp(x - x0))

def term(matrix, t, x): # This function is used to calculate the time evolution of the matrix
    global rho
    global q

    # Boundary condition 0, if x == 0 then matrix[t-1,x-1] is out of bounds, thus = 0
    if x == 0:
        return rho * matrix[t-1,x] * (1 - matrix[t-1,x]/q) - matrix[t-1,x]/(1+matrix[t-1,x]) + (matrix[t-1,x+1]-2*matrix[t-1,x])
    
    # Boundary condition L, if x == Length-1 then matrix[t-1,x+1] is out of bounds, thus = 0
    if x == Length-1:
        return rho * matrix[t-1,x] * (1 - matrix[t-1,x]/q) - matrix[t-1,x]/(1+matrix[t-1,x]) + (matrix[t-1,x-1]-2*matrix[t-1,x])
    
    # General case, no boundary condition, all x's are in bounds
    return rho * matrix[t-1,x] * (1 - matrix[t-1,x]/q) - matrix[t-1,x]/(1+matrix[t-1,x]) + (matrix[t-1,x+1]+matrix[t-1,x-1]-2*matrix[t-1,x])


## Main code

# Initialize the matrix
matrixSpaceTime = np.zeros((Time,Length))


# Set the initial values 
for i in range(0,Length):
    matrixSpaceTime[0,i] = u(i+1)


for t in range(1,Time):
    for x in range(0,Length-1):
        matrixSpaceTime[t,x] = matrixSpaceTime[t-1, x] + dt * term(matrixSpaceTime, t, x)
    

# plot a heatmap of the matrix and use a fixed colorbar
# plt.imshow(matrixSpaceTime, cmap='hot', interpolation='nearest')
# plt.colorbar()
# plt.title('Heatmap of the matrix')
# # axis ratio 1:5
# plt.show()

showInits = [0,19,39,59,79,99]
plt.figure()

for i in range(len(showInits)):
    plt.plot(matrixSpaceTime[showInits[i],:], label=str(showInits[i]))

plt.legend()

plt.show()