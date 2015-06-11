from scipy.spatial.distance import cdist
import random

def H(A,B):
    m, n = len(A), len(B)
    distmatrix = cdist(A,B)
    minsa = [min(distmatrix[i,:]) for i in xrange(m)]
    minsb = [min(distmatrix[:,i]) for i in xrange(n)]
    return max(max(minsa),max(minsb))
    
if __name__ == '__main__':
    X = [(random.uniform(0,10),random.uniform(0,10)) for i in xrange(50)]
    Y = [(random.uniform(0,10),random.uniform(0,10)) for i in xrange(50)]
    Z = [(random.uniform(0,20),random.uniform(0,20)) for i in xrange(50)]
    print H(X, Y), H(X,Z)
