import numpy as np
from scipy.optimize import linear_sum_assignment
from sympy import primerange



class Microstructure:
    def __init__(self, img: np.ndarray, grains: np.ndarray):
        self.img = img
        self.grains = grains
class MicrostructureCompute:
    """
    A class to compute microstructure properties and relabel grains.
    """
    def __init__(self):
        pass

    def create_microstructure(self, img, primeList):
        """
        Create a Microstructure object from an image.

        Parameters:
        img (numpy.ndarray): image representing the microstructure.

        Returns:
        Microstructure: An instance of the Microstructure class.
        """
        # Relabel all grains at once so a new prime can't collide with an old label
        labels, inverse = np.unique(img, return_inverse=True)
        primes = [primeList.pop(0) for _ in labels]
        img = np.array(primes, dtype=int)[inverse].reshape(img.shape)

        grains = np.array(np.unique(img), dtype=int)
        return Microstructure(img, grains)


    def compute_overlaid_structure(self, M, N) -> Microstructure:
        """
        Compute the overlaid structure of two microstructures.

        Parameters:
        M (Microstructure): The starting microstructure.
        N (Microstructure): The target microstructure.

        Returns:
        Microstructure: A new Microstructure object representing the overlaid structure.
        """
        X_img = np.array(np.multiply(M.img, N.img), dtype=int)
        grains = np.array(np.unique(X_img), dtype=int)
        return Microstructure(X_img, grains)
        

def VolumeExchange(M_img: np.ndarray, N_img: np.ndarray, return_overlay: bool = False):
    """
    Calculate the volume exchange required to transform microstructure M into N.

    Parameters:
    M (np.ndarray): The starting microstructure image with labeled grains.
    N (np.ndarray): The target microstructure image with labeled grains.
    return_overlay (bool): If True, also return the overlaid structure and matched mask.

    Returns:
    float: The total volume exchange required.
    If return_overlay is True, returns (total_cost, X_img, matched_mask) where
    X_img is the overlaid structure (product of prime-labeled M and N) and
    matched_mask is True on the overlap area of optimally matched grain pairs.
    """

    # M_img = M_img.reshape((20,20))
    # N_img = N_img.reshape((20,20))
    
    # =================== Create MicrostructureCompute instances ===================
    microstructureCompute = MicrostructureCompute()
    primeList = list(primerange(10, 10000))
    M = microstructureCompute.create_microstructure(M_img, primeList)
    N = microstructureCompute.create_microstructure(N_img, primeList)
    X = microstructureCompute.compute_overlaid_structure(M, N)
    # ==============================================================================

    # ============================ Helper Functions ================================
    def getArea(X, g_x):
        return np.sum(X.img == g_x)/X.img.size

    def getIntersectionGrains(M, N, X, g_m, g_n):
        mMask = (M.img == g_m)
        nMask = (N.img == g_n)
        intersection = np.logical_and(mMask, nMask)
        intersection = np.multiply(intersection, X.img)
        unique_grains = np.unique(intersection)
        unique_grains = np.delete(unique_grains, np.where(unique_grains == 0))
        return unique_grains
    # ==============================================================================

    # ========================= Volume Exchange Calculation ========================

    # Initialize square cost matrix with zeros
    n_m = len(M.grains)
    n_n = len(N.grains)
    n = max(n_m, n_n)
    cost_matrix = np.zeros((n, n))
    
    # Fill cost matrix with negative overlap areas
    for i, g_m in enumerate(M.grains):
        for j, g_n in enumerate(N.grains):
            overlap_grains = getIntersectionGrains(M, N, X,g_m, g_n)
            overlap_area = sum(getArea(X, g_x) for g_x in overlap_grains)
            cost_matrix[i, j] = -overlap_area  # Negative because we minimize cost

    # Use Jonker-Volgenant-Castanon algorithm (linear assignment) to find optimal assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    total_negative_overlap = cost_matrix[row_ind, col_ind].sum()
    total_overlap = -total_negative_overlap
    
    # Calculate total cost
    total_cost = 1 - total_overlap
    if not return_overlay:
        return total_cost

    # Matched area: overlaid grains g_m * g_n for each assigned (non-dummy) pair
    matched_mask = np.zeros(X.img.shape, dtype=bool)
    for i, j in zip(row_ind, col_ind):
        if i < n_m and j < n_n:
            matched_mask |= (X.img == M.grains[i] * N.grains[j])
    return total_cost, X.img, matched_mask

