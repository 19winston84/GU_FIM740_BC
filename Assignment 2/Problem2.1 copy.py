import numpy as np
import matplotlib.pyplot as plt


# Parameters
Length = 100
Time = 10000
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

def derivateDUDX(matrix):
    derivate = np.zeros((Time,Length))
    for t in range(0,Time):
        for x in range(0,Length-1):
            if x == 0:
                derivate[t,x] = (matrix[t,x+1] - matrix[t,x]) / (x+1-x)
            else:
                derivate[t,x] = (matrix[t,x+1] - matrix[t,x-1]) / (x+1-x)
    return derivate

def approxVelocity(matrix, t0, t1, us):
    global dt
    # find in the matrix[t1,:] and matrix[t0,:] the values of us
    # and compared the position of the values

    # Find the position of the values, we need to consider a small tolerance of 0.01
    
    pos0 = np.where(np.abs(matrix[t0,:] - us) < 0.01)[0][0]
    pos1 = np.where(np.abs(matrix[t1,:] - us) < 0.01)[0][0]

    
    
    print(pos1, pos0)

    # Return the velocity
    return (pos1 - pos0) / (t1 - t0) / dt
    


## Main code

# Initialize the matrix
matrixSpaceTime = np.zeros((Time,Length))
matrixSpaceTime2 = np.zeros((Time, Length))


# Set the initial values 
for i in range(0,Length):
    matrixSpaceTime[0,i] = u(i+1)
    matrixSpaceTime2[0,i] = u(i+1)


for t in range(1,Time):
    for x in range(0,Length-1):
        matrixSpaceTime[t,x] = matrixSpaceTime[t-1, x] + dt * term(matrixSpaceTime, t, x)
        matrixSpaceTime2[t,x] = matrixSpaceTime2[t-1, x] + dt * term(matrixSpaceTime2, t, x)
    

# Plot the matrix at different times in 999 steps
showInits = [0,999,1999,2999,3999,4999,5999,6999,7999,8999,9999]
showInits = [0]

derivatedDUDX = derivateDUDX(matrixSpaceTime)

# plt.figure(figsize=(12, 6))

# plt.subplot(1, 2, 1)
# for i in range(len(showInits)):
#     plt.plot(matrixSpaceTime[showInits[i],:], label='At time ' + str(showInits[i]))
# plt.legend()
# plt.xlabel('Position xi')
# plt.ylabel('u(x)')
# plt.grid()

# plt.subplot(1, 2, 2)
# # plot u vs du/dx
# for i in range(len(showInits)):
#     plt.plot(matrixSpaceTime[showInits[i],:], derivatedDUDX[showInits[i],:], label='At time ' + str(showInits[i]))
# plt.legend()
# plt.xlabel('u(x)')
# plt.ylabel('du/dx')
# plt.grid()

# plt.show()

print(approxVelocity(matrixSpaceTime, 1000, 2000, 3))

