import numpy as np

def build_forcing_vector(b, nmax):
    """
    Given:
      b: a NumPy array of shape (n,) representing the external input B[k+1]
         (the forcing applied to the current state),
      nmax: an integer, the maximum delay (i.e., the number of blocks in the augmented state).

    Returns:
      U: a NumPy column vector of shape (n * nmax, 1) representing the forcing vector
         for the augmented state–space system. It is defined as:

         U = [ b ]
             [ 0 ]
             [ . ]
             [ 0 ]

         where b appears in the top n rows and the remaining (nmax-1)*n rows are zero.
    """
    # Determine the number of nodes (n)
    n = b.shape[0]

    # Create a zero vector of length n * nmax.
    U = np.zeros((n * nmax, 1))

    # Set the top n entries equal to b.
    # If b is a 1D array, reshape it to a column vector.
    U[:n, 0] = b.reshape(-1)

    return U

def build_augmented_matrix(alpha, D, t_res, lambda_decay):
    """
    Given:
      - alpha: an (n x n) NumPy array with entries α_{j→i}.
               (For non-existing edges, the entry should be 0.)
      - D:     an (n x n) NumPy array with the integer delays (in time steps)
               associated with each edge from node j to node i.
               (For non-existing edges, set D[i,j] = 0.)
      - t_res: an (n x n) NumPy array with the residence times for each edge.
               (For non-existing edges, set t_res[i,j] = 0.)
      - lambda_decay: a scalar decay constant.

    The function computes the decay factor matrix:
          E = exp(-lambda_decay * t_res)
    (elementwise) and then forms the matrix
          A_prod = alpha * E
    (elementwise multiplication).

    Next, it splits A_prod into N_max blocks according to D, where:
          N_max = max(D)
    and for each m = 1, 2, ..., N_max, the block A_m (of size n x n) is defined by:
          (A_m)[i, j] = A_prod[i, j] if D[i, j] == m,
                        0           otherwise.

    Finally, it constructs and returns the augmented matrix A_aug:

         A_aug = [ A_1   A_2   ...   A_Nmax ]
                 [ I_n   0     ...    0     ]
                 [  0    I_n   ...    0     ]
                 [ ...   ...   ...   ...    ]
                 [  0     0    ...   I_n    0 ]

    where I_n is the (n x n) identity matrix and there are (N_max - 1) companion rows.

    Parameters:
      alpha:        NumPy array of shape (n, n) with α_{j→i} values.
      D:            NumPy array of shape (n, n) with integer delays.
      t_res:        NumPy array of shape (n, n) with residence times.
      lambda_decay: Scalar decay constant.

    Returns:
      A_aug: NumPy array of shape (n * N_max, n * N_max) representing the augmented matrix.
    """
    n = alpha.shape[0]
    N_max = int(np.max(D))  # maximum delay among all entries

    # Compute the decay factor matrix:
    E = np.exp(-lambda_decay * t_res)

    # Compute the elementwise product:
    A_prod = alpha * E

    # Create a list to store the blocks A_1, A_2, ..., A_N_max.
    A_blocks = []
    for m in range(1, N_max + 1):
        # For each m, create a block that retains the value from A_prod only where D equals m.
        A_m = np.where(D == m, A_prod, 0)
        A_blocks.append(A_m)

    # Build the top block row by horizontally stacking the A_m blocks.
    top_row = np.hstack(A_blocks)  # shape: (n, n * N_max)

    # Build the companion (shift) rows.
    companion_rows = []
    for i in range(N_max - 1):
        # For each companion row, we need N_max blocks of size (n x n).
        # In companion row i, place the identity matrix I_n in the i-th block and zeros elsewhere.
        row_blocks = []
        for j in range(N_max):
            if j == i:
                row_blocks.append(np.eye(n))
            else:
                row_blocks.append(np.zeros((n, n)))
        companion_row = np.hstack(row_blocks)
        companion_rows.append(companion_row)

    # Vertically stack the top row and the companion rows.
    A_aug = np.vstack([top_row] + companion_rows)
    return A_aug

def build_augmented_matrix_2(alpha, E, D):
    """
    Given three matrices:
      - alpha: an (n x n) NumPy array with entries α_{j→i}. (Non-existing edges: 0.)
      - E:     an (n x n) NumPy array with entries exp(-λ * t_{j,i}). (Non-existing edges: 0.)
      - D:     an (n x n) NumPy array with the integer delays (in time steps)
               associated with each edge from node j to node i.
               (For non-existing edges, set D[i, j] = 0.)

    The function computes:
       A_prod = alpha * E   (elementwise multiplication)
    and then splits A_prod into N_max blocks according to the delays in D,
    where
       N_max = max(D)
    and for m = 1, 2, ..., N_max, defines the block A_m (of size n x n) by:

       (A_m)[i, j] = A_prod[i, j]   if D[i, j] == m,
                     0             otherwise.

    Then it constructs and returns the augmented matrix A_aug defined by:

         A_aug = [ A_1   A_2   ...   A_Nmax ]
                 [ I_n   0     ...    0     ]
                 [  0    I_n   ...    0     ]
                 [ ...   ...   ...   ...    ]
                 [  0     0    ...   I_n    0 ]

    where I_n is the (n x n) identity matrix and there are (N_max - 1) companion rows.

    Parameters:
      alpha: NumPy array of shape (n, n) with α_{j→i} values.
      E:     NumPy array of shape (n, n) with exp(-λ * t_{j,i}) values.
      D:     NumPy array of shape (n, n) with integer delays.

    Returns:
      A_aug: NumPy array of shape (n*N_max, n*N_max) representing the augmented matrix.
    """
    # Determine the number of nodes and maximum delay.
    n = alpha.shape[0]
    N_max = int(np.max(D))  # maximum delay among all entries

    # Compute the elementwise product:
    A_prod = alpha * E

    # Create a list to store the blocks A_1, A_2, ..., A_Nmax.
    A_blocks = []
    for m in range(1, N_max + 1):
        # For each m, create a block that has the value from A_prod only where D equals m.
        A_m = np.where(D == m, A_prod, 0)
        A_blocks.append(A_m)

    # Build the top block row by horizontally stacking the A_m blocks.
    top_row = np.hstack(A_blocks)  # shape: (n, n * N_max)

    # Build the companion (shift) block rows.
    # There will be N_max - 1 companion rows.
    companion_rows = []
    for i in range(N_max - 1):
        # For each companion row, we need N_max blocks of size (n x n).
        # In companion row i (starting from i=0 for the first companion row),
        # place the identity I_n in block column i and zeros in all other blocks.
        row_blocks = []
        for j in range(N_max):
            if j == i:
                row_blocks.append(np.eye(n))
            else:
                row_blocks.append(np.zeros((n, n)))
        companion_row = np.hstack(row_blocks)
        companion_rows.append(companion_row)

    # Vertically stack the top row and the companion rows.
    A_aug = np.vstack([top_row] + companion_rows)
    return A_aug

def build_rr_augmented_matrix(alpha, D, t_res, lambda_decay):
    """
    Builds a special augmented matrix that models the temporal evolution of inner sources
    along edges with delays.

    Inputs:
      alpha       : (n x n) NumPy array with the α_{j→i} values (0 if no edge).
      D           : (n x n) NumPy array with integer delays (in time steps) for each edge.
                    (For non-existing edges, D[i,j] should be 0.)
      t_res       : (n x n) NumPy array with the residence times t_res_{i,j} (same convention).
      lambda_decay: A scalar decay constant.

    For an edge from node j to node i with delay m = D[i,j] and residence time t_res[i,j]:
      - If m == 1, the only block factor is:
            factor = α_{j→i} * (1 - exp(-lambda_decay * t_res[i,j]))
      - If m > 1, then for each delay stage l = 1,..., m, the factor is:
            factor(l) = α_{j→i} * (1 - exp(-lambda_decay*(t_res[i,j]/m))) *
                        exp(-lambda_decay*(t_res[i,j]/m)*(l-1))
      - For l > m, the contribution is 0.

    The function creates block matrices A_l (for l = 1,..., N_max, where N_max = max(D))
    of size (n x n) with:
       (A_l)[i,j] = factor(l) if D[i,j] >= l, 0 otherwise.

    Then, it constructs the augmented matrix:

        A_aug = [ A_1   A_2   ...   A_{N_max} ]
                [ I_n   0     ...    0       ]
                [ 0    I_n    ...    0       ]
                [ ...   ...   ...   ...      ]
                [ 0     0     ...   I_n    0  ]

    which will be of size (n*N_max x n*N_max). This matrix is meant to be multiplied
    with a state vector S (of length n*N_max) that contains the source term factors of the nodes
    for all delay stages.

    Returns:
      A_aug: The special augmented matrix (NumPy array) of shape (n*N_max, n*N_max).
    """
    n = alpha.shape[0]
    N_max = int(np.max(D))  # maximum delay (in time steps)

    # Initialize a list for the block matrices A_l (l=1,...,N_max)
    A_blocks = []

    # Loop over each delay stage l.
    for l in range(1, N_max + 1):
        # Initialize an (n x n) matrix for this block.
        A_l = np.zeros((n, n))
        # Loop over all node indices i, j.
        for i in range(n):
            for j in range(n):
                m = int(D[i, j])  # delay for edge from j to i
                if m > 0 and l <= m:
                    if m == 1:
                        # Only one delay unit; use the simple factor.
                        factor = alpha[i, j] * (1 - np.exp(-lambda_decay * t_res[i, j]))
                    else:
                        # For m > 1, spread the residence time over m subintervals.
                        delta = t_res[i, j] / m
                        factor = alpha[i, j] * (1 - np.exp(-lambda_decay * delta)) \
                                 * np.exp(-lambda_decay * delta * (l - 1))
                    A_l[i, j] = factor
                else:
                    A_l[i, j] = 0.0
        A_blocks.append(A_l)

    # Build the top block row by horizontally stacking the A_l blocks.
    top_row = np.hstack(A_blocks)  # shape: (n, n * N_max)

    # Build the companion (shift) block rows: there will be (N_max - 1) rows.
    companion_rows = []
    for i in range(N_max - 1):
        row_blocks = []
        for j in range(N_max):
            if j == i:
                row_blocks.append(np.eye(n))
            else:
                row_blocks.append(np.zeros((n, n)))
        companion_row = np.hstack(row_blocks)
        companion_rows.append(companion_row)

    # Stack the top row and the companion rows vertically.
    A_aug = np.vstack([top_row] + companion_rows)
    return A_aug



def build_decay_average_augmented_matrix(D, t_res, lambda_decay):
    """
    Builds an augmented matrix that models the temporal evolution of the average species
    concentration inside each node, based solely on the delay and residence time data.

    Inputs:
      D           : (n x n) NumPy array with integer delays for each edge.
                    (For non-existing edges, the entry should be 0.)
                    In each node, we assume all nonzero delays are equal.
      t_res       : (n x n) NumPy array with residence times for each edge.
                    (For non-existing edges, the entry should be 0.)
                    In each node, assume all nonzero residence times are equal.
      lambda_decay: A scalar decay constant.

    For each node i (i=0,...,n-1):
      Let d_i = max(D[i, :])   and T_i = max(t_res[i, :])  (assumed nonzero if there is any edge)
      Define, for delay stage l = 1,2,...,d_i, the factor:

         factor_i(l) = (1/d_i) * exp(-lambda_decay * ((l-1)*T_i/d_i)) *
                        (1 - exp(-lambda_decay * (T_i/d_i))) / (lambda_decay*(T_i/d_i))

      For l > d_i, define factor_i(l) = 0.

    Then for each delay stage l = 1,..., N_max (where N_max = max_i d_i),
    define the (n x n) diagonal matrix A_l by:
         (A_l)[i,i] = factor_i(l)
         (A_l)[i,j] = 0 for i != j.

    The special augmented matrix is then constructed as:

         A_aug = [ A_1  A_2  ...  A_{N_max} ]
                 [ I_n   0    ...   0      ]
                 [  0   I_n   ...   0      ]
                 [ ...  ...   ...  ...     ]
                 [  0    0    ...  I_n   0  ]

    This matrix has size (n*N_max) x (n*N_max) and is intended to be used in an augmented
    equation that multiplies a state vector X (of length n*N_max) that contains the source term factors.

    Returns:
      A_aug : The special augmented matrix as a NumPy array.
    """
    n = D.shape[0]

    # For each node, obtain its delay and residence time.
    # We assume that for node i, if there is any nonzero delay then all nonzero entries in row i are equal.
    d_vec = np.zeros(n)  # delay per node
    T_vec = np.zeros(n)  # residence time per node
    for i in range(n):
        # For row i, take the maximum nonzero delay
        nonzero_delays = D[i, D[i, :] > 0]
        if nonzero_delays.size > 0:
            d_vec[i] = np.max(nonzero_delays)
            # Similarly, for t_res, take the maximum of nonzero entries.
            nonzero_tres = t_res[i, t_res[i, :] > 0]
            T_vec[i] = np.max(nonzero_tres)
        else:
            # If no incoming edge, we can set d_vec[i] = 1 and T_vec[i] = 0 (or leave 0)
            d_vec[i] = 1
            T_vec[i] = 0

    # Determine the overall maximum delay
    N_max = int(np.max(d_vec))

    # Build the block matrices A_l for l=1,..., N_max.
    A_blocks = []
    for l in range(1, N_max + 1):
        # For each node i, compute its factor at delay stage l.
        factors = np.zeros(n)
        for i in range(n):
            d_i = d_vec[i]
            T_i = T_vec[i]
            if l <= d_i and T_i > 0:
                delta = T_i / d_i
                # Compute the factor according to the formula:
                # factor = (1/d_i) * exp(-lambda_decay*(l-1)*delta) *
                #          (1 - exp(-lambda_decay*delta)) / (lambda_decay*delta)
                factors[i] = (1.0/d_i) * np.exp(-lambda_decay*(l-1)*delta) * \
                             (1 - np.exp(-lambda_decay*delta)) / (lambda_decay*delta)
            else:
                factors[i] = 0.0
        # Construct an (n x n) diagonal matrix with these factors.
        A_l = np.diag(factors)
        A_blocks.append(A_l)

    # Build the top block row by horizontally stacking A_blocks.
    top_row = np.hstack(A_blocks)  # shape: (n, n * N_max)

    # Build the companion rows (there will be N_max - 1 of them).
    companion_rows = []
    for i in range(N_max - 1):
        row_blocks = []
        for j in range(N_max):
            if j == i:
                row_blocks.append(np.eye(n))
            else:
                row_blocks.append(np.zeros((n, n)))
        companion_row = np.hstack(row_blocks)
        companion_rows.append(companion_row)

    # Stack the top row and the companion rows vertically.
    A_aug = np.vstack([top_row] + companion_rows)
    return A_aug

def build_rr_average_augmented_matrix(D, t_res, lambda_decay):
    """
    Builds a special augmented matrix that models the temporal evolution
    of the average species concentration inside each node due solely to
    the evolution of the node's source term.

    Inputs:
      D           : (n x n) NumPy array with integer delays (in time steps) for each edge.
                    For each node i, assume that all nonzero entries in row i are equal;
                    let d_i denote that common delay (if a node has no incoming edges, d_i can be set to 1).
      t_res       : (n x n) NumPy array with residence times for each edge.
                    For each node i, assume that all nonzero entries in row i are equal;
                    let T_i denote that common residence time.
      lambda_decay: Scalar decay constant.

    The procedure for each node i:
      1. Set m_i = d_i.
      2. Compute t_res_d_i = T_i / m_i.
      3. Compute decay factor: df_i = exp(-lambda_decay * t_res_d_i)
      4. Compute reaction-rate decay factor: rrdf_i = 1 - df_i.
      5. Compute average inlet factor: aif_i = rrdf_i / (lambda_decay * t_res_d_i)
      6. Compute average reaction rate factor: arrf_i = 1 - aif_i.

    Then for each delay stage l = 1,2,..., N_max (with N_max = max_i m_i):
      For node i:
         if l <= m_i:
            if l == 1:
               factor_i(1) = (1/m_i) * arrf_i
            else:  # l > 1
               factor_i(l) = (1/m_i) * rrdf_i * (df_i ** (l-1)) * aif_i
         else:
            factor_i(l) = 0.

    For each delay stage l, form an (n x n) diagonal matrix A_l with:
         (A_l)[i,i] = factor_i(l).

    Finally, construct the augmented matrix:
         A_aug = [ A_1  A_2  ...  A_N_max ]
                 [ I_n   0   ...   0      ]
                 [  0   I_n  ...   0      ]
                 [ ...  ...  ...  ...     ]
                 [  0    0   ...  I_n   0  ]
    (This matrix has size (n*N_max) x (n*N_max).)

    Returns:
      A_aug: The special augmented matrix as a NumPy array.
    """
    n = D.shape[0]

    # For each node i, determine d_i and T_i.
    d_vec = np.zeros(n)
    T_vec = np.zeros(n)
    for i in range(n):
        # Take the first nonzero delay in row i (or if none, set to 1)
        nonzero_D = D[i, D[i, :] > 0]
        if nonzero_D.size > 0:
            d_vec[i] = nonzero_D[0]  # assume all nonzero are equal
        else:
            d_vec[i] = 1  # default if no incoming edges

        nonzero_tres = t_res[i, t_res[i, :] > 0]
        if nonzero_tres.size > 0:
            T_vec[i] = nonzero_tres[0]  # assume all nonzero are equal
        else:
            T_vec[i] = 0  # or some default value

    # For each node i, set m_i = d_i.
    m_vec = d_vec.copy()

    # Overall maximum delay over all nodes.
    N_max = int(np.max(m_vec))

    # For each node i, compute the per-unit residence time and factors.
    # We will store these in arrays of length n.
    t_res_d = np.zeros(n)
    df = np.zeros(n)
    rrdf = np.zeros(n)
    aif = np.zeros(n)
    arrf = np.zeros(n)

    for i in range(n):
        m_i = m_vec[i]
        T_i = T_vec[i]
        if m_i > 0 and T_i > 0:
            t_res_d[i] = T_i / m_i
            df[i] = np.exp(-lambda_decay * t_res_d[i])
            rrdf[i] = 1 - df[i]
            # Avoid division by zero:
            if lambda_decay * t_res_d[i] != 0:
                aif[i] = rrdf[i] / (lambda_decay * t_res_d[i])
            else:
                aif[i] = 0
            arrf[i] = 1 - aif[i]
        else:
            t_res_d[i] = 0
            df[i] = 1
            rrdf[i] = 0
            aif[i] = 0
            arrf[i] = 0

    # Build the block matrices A_l for l = 1,..., N_max.
    A_blocks = []
    for l in range(1, N_max + 1):
        # For each node i, compute its factor for delay stage l.
        factors = np.zeros(n)
        for i in range(n):
            m_i = m_vec[i]
            if l <= m_i:
                if l == 1:
                    factors[i] = (1.0 / m_i) * arrf[i]
                else:
                    factors[i] = (1.0 / m_i) * rrdf[i] * (df[i] ** (l - 1)) * aif[i]
            else:
                factors[i] = 0.0
        # Construct the diagonal matrix for this delay stage.
        A_l = np.diag(factors)
        A_blocks.append(A_l)

    # Build the top block row by horizontally stacking the A_l blocks.
    top_row = np.hstack(A_blocks)  # shape: (n, n * N_max)

    # Build the companion (shift) rows. There are (N_max - 1) such rows.
    companion_rows = []
    for i in range(N_max - 1):
        row_blocks = []
        for j in range(N_max):
            if j == i:
                row_blocks.append(np.eye(n))
            else:
                row_blocks.append(np.zeros((n, n)))
        companion_row = np.hstack(row_blocks)
        companion_rows.append(companion_row)

    # Stack the top row and companion rows vertically.
    A_aug = np.vstack([top_row] + companion_rows)
    return A_aug