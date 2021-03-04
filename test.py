import scipy.sparse.linalg as sla
import numpy as np 
from scipy.sparse import spdiags
from scipy import sparse
import time
import pickle
from scipy import sparse, io
# up=np.arange(2400.0)
# down=up
# data=np.stack((up,down))
# mat=spdiags(data,np.array([-1,1]),2400,2400)
with open('ham.pickle','rb') as f:
    z=pickle.load(f)
# z=sparse.dia_matrix(z)
# start=time.time()
z=(z+z.T)/2
z=np.real(z).tocsc()
# val,vec=LA.eigsh(z,k=1000,sigma=0,v0=)
# print(time.time()-start)
# print(val)
# io.savemat('H.mat',dict(mat=z))

class SpSolve(sla.LinearOperator):
    """
    Use spsolve to directly solve Ax=B
    Requires: scikit-umfpack and SuiteSparse
    https://scikit-umfpack.github.io/scikit-umfpack
    """
    def __init__(self,M):
        self.M=M
        self.shape = M.shape

    def __matvec(self,x):
        return sla.spsolve(self.M,x)

zop=SpSolve(z)
zop(np.random.rand(2600,1))