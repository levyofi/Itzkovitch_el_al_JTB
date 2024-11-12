import numpy as np

def semi_random_sample(matrix_size, num_samples):
    submat_size = matrix_size // int(np.sqrt(num_samples))
    submat_steps = int(matrix_size / submat_size)
    
    points = []
    for i in range(submat_steps): 
        for j in range(submat_steps):
            x_start, x_end = i * submat_size, (i + 1) * submat_size
            y_start, y_end = j * submat_size, (j + 1) * submat_size
    
            x = np.random.randint(x_start, x_end)
            y = np.random.randint(y_start, y_end)
            
            points.append((x, y))
            
    while len(points) < 1000:
        x = np.random.randint(0, matrix_size)
        y = np.random.randint(0, matrix_size)
        points.append((x, y))
        
    return points
