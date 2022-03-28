import numpy as np
import matplotlib.pyplot as plt
# the dielectric constants
n_0 = 1
n_1 = 1.5
n_2 = 2.5
pi = np.pi
N = 5
c = 1
# the wavelength we want to remove
lambda_removed = 450
d_1 = lambda_removed/(4*n_1)
d_2 = lambda_removed/(4*n_2)
# the P matrix for the first type of glass

for lambdas in range(300, 800):
    phi_1 = n_1 * ((2 * pi) / lambdas) * d_1
    A_1 = np.exp((-1j * phi_1))
    D_1 = np.exp((1j * phi_1))
    P1_matrix = np.matrix([[A_1, 0], [0, D_1]])

    # the P matrix for the second type of glass
    phi_2 = n_2 * ((2 * pi) / lambdas) * d_2
    A_2 = np.exp((-1j * phi_2))
    D_2 = np.exp((1j * phi_2))
    P2_matrix = np.matrix([[A_2, 0], [0, D_2]])

    # the collision between layers matrix 1
    BC_1 = n_2 - n_1
    AD_1 = n_1 + n_2
    B12_matrix = np.matrix([[AD_1, BC_1], [BC_1, AD_1]])
    B12_constant = 1 / (2 * n_2)
    B12 = B12_constant * B12_matrix

    # the collision between layers 2 and 1
    BC_2 = n_1 - n_2
    AD_2 = n_1 + n_2
    B21_matrix = np.matrix([[AD_2, BC_2], [BC_2, AD_2]])
    B21_constant = 1 / (2 * n_1)
    B21 = B21_constant * B21_matrix

    # the collision between air and layer 1
    BC_0 = n_1 - n_0
    AD_0 = n_1 + n_0
    B01_matrix = np.matrix([[AD_0, BC_0], [BC_0, AD_0]])
    B01 = B21_constant * B01_matrix

    # the collision between layer 2 and air
    BC_02 = n_0 - n_2
    AD_02 = n_0 + n_2
    B20_constant = 1 / (2 * n_0)
    B20_matrix = np.matrix([[AD_02, BC_02], [BC_02, AD_02]])
    B20 = B20_constant * B20_matrix

    # the inverse of B21 matrix
    inverse_B21 = np.linalg.inv(B21)
    # creating the M matrix
    # the M matrix should be = B20 * B21^-1 * [B21 * P2 * B12 *P1]^N * B01, where N is the number of layers.
    element_1 = np.dot(B20, inverse_B21)
    subelement_1 = np.dot(B21, P2_matrix)
    subelement_2 = np.dot(subelement_1, B12)
    subelement_3 = np.dot(subelement_2, P1_matrix)
    element_2 = np.linalg.matrix_power(subelement_3, N)
    element_3 = np.dot(element_1, element_2)
    M_matrix = np.dot(element_3, B01)
    M22 = M_matrix.A[1][1]
    transmission = 1 / (M22)
    T = abs(transmission)
    def T2(c):
        return T
    c = c + 1
    plt.plot(lambdas, T2(c), 'b.')
    width = [450]
    if T < 0.5:
        width.append(lambdas)
        x1 = min(width)
        x2 = max(width)
gap = x2 - x1
print(gap)
plt.xlabel("wavelength (nm)")
plt.ylabel("transmission")
plt.show()

