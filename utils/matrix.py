import numpy as np

def matrix_multiply(matrix1, matrix2):
    # Convert inputs to numpy arrays if they aren't already
    m1 = np.array(matrix1)
    m2 = np.array(matrix2)
    
    # Use numpy's matrix multiplication
    result = np.matmul(m1, m2)
    
    # Convert back to list of lists for compatibility
    return result.tolist()
    