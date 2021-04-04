import kwant
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as sla
from scipy.sparse import csc_matrix
import time
import warnings

import kwant.linalg.mumps as mumps
from scipy.sparse import identity


class Params:
    '''
    generate bilayer hetereostructure
    example: params=Params(g=1,ED=433)
    
    '''
    def __init__(self,
        a=np.array([0.1,0.1,0.1]), # lattice constant (nm)
        L_Al=np.array([10,10,5]),   # size of Al (nm)
        L_FM=np.array([2,10,5]), # size of FM (nm)
        m_Al=1,    # effective mass in Al in m_e, where m_e is 0.511 MeV
        m_FM=0.3,   # effective mass in FM in  m_e, where m_e is 0.511 MeV
        EF_Al=11, # Chemical potential of Al in eV
        EF_FM=-0.5,
        h_exc=0,  # exchange field in eV, acts like Zeeman field
        g=5.6,   # Coupling strength in eV
        E_D=1,    # Debye temperature in 433 K
        T=0,    # Temperature in eV      
        U_D=0,     # variance of disorder in eV
        Delta_0=3.5e-4,    # the initial SC gap
        periodic_boundary_condition=1, # 1: periodic boundary condition
        verbose=True,   #verbose
        store_history=False
        ):
        assert L_Al[1]==L_FM[1] and L_Al[2]==L_FM[2],"The y,z dimension of Al ({},{}) not same as FM ({},{})".format(L_Al[1],L_Al[2],L_FM[1],L_FM[2])
        self.a=np.array(a)*1e-9*5.076e6    # eV^-1
        self.L_Al=np.array(L_Al)*1e-9*5.076e6    # eV^-1
        self.L_FM=np.array(L_FM)*1e-9*5.076e6   # eV^-1
        self.m_Al=m_Al*0.511e6    # eV
        self.m_FM=m_FM*0.511e6    # eV
        self.EF_Al=EF_Al        # eV
        self.EF_FM=EF_FM        # eV
        self.h_exc=h_exc        # eV
        self.g=g            # eV
        self.E_D=E_D*433*8.617333262e-5       # eV
        self.T=T            # eV
        self.U_D=U_D  # eV
        self.Delta_0=Delta_0    #initial value of SC gap
        self.periodic_boundary_condition=periodic_boundary_condition
        self.verbose=verbose


        self.N_Al=np.rint(self.L_Al/self.a).astype(int)    # dimension of sites in Al
        self.N_FM=np.rint(self.L_FM/self.a).astype(int)    # dimension of sites in FM
        self.s0=np.array([[1,0],[0,1]])
        self.sx=np.array([[0,1],[1,0]])
        self.sy=np.array([[0,-1j],[1j,0]])
        self.sz=np.array([[1,0],[0,-1]])
        self.b=2*np.pi/self.a    # reciprocal unit vector in eV
        self.t_Al=1/(2*self.m_Al*self.a**2)    # NN hopping of Al in eV
        self.t_FM=1/(2*self.m_FM*self.a**2)    # NN hopping of FM in eV
        self.m_int=np.sqrt(self.m_Al*self.m_FM)
        self.t_int=1/(2*self.m_int*self.a**2)    # NN hopping at the interface in eV

        self.disorder={}    # generate empty dict for disorder profile
        self.Delta={}
        self.Delta_array=[]
        # self.Delta=np.zeros((self.N[0],self.N[1]))  
        # self.wfall=[[] for _ in range(self.N[2])]
        # self.energyall=[[] for _ in range(self.N[2])]    
        self.wfall=[]
        self.energyall=[]
        self.store_history=store_history
        if store_history:
            self.wfall_history=[]
            self.energyall_history=[]

        self.make_system()
        self.generate_disorder()
        self.update_Delta()
        self.energy_kz()
        self.estimate_eig()



    def energy_kz(self):
        kindex=np.arange(self.N_Al[2])
        uz=(2*kindex-self.N_Al[2]+1)/(2*self.N_Al[2])
        self.kz_list=uz*self.b[2]    # eV^-1
        self.energy_Al=2*self.t_Al[2]*(1-np.cos(self.a[2]*self.kz_list))   # bandstructure of Al
        self.energy_FM=2*self.t_FM[2]*(1-np.cos(self.a[2]*self.kz_list))   # bandstructure of Al

    
    def estimate_eig(self):
        k_list={}
        for i in range(2):
            k_list[i]=(2*np.arange(self.N_Al[i])-self.N_Al[i]+1)/(2*self.N_Al[i])*self.b[i]
        [kx_grid,ky_grid]=np.meshgrid(k_list[0],k_list[1])
        self.estimate=[]
        for kz in self.kz_list:
            E_mesh=2*self.t_Al[0]*(1-np.cos(kx_grid*self.a[0]))+2*self.t_Al[1]*(1-np.cos(ky_grid*self.a[1]))+2*self.t_Al[2]*(1-np.cos(kz*self.a[2]))-self.EF_Al
            self.estimate.append(4*np.count_nonzero((E_mesh<self.E_D)*(E_mesh>-self.E_D)))        
        self.estimate=np.array(self.estimate)
        
    def Fermi_dist(self):
        if self.T==0:
            self.F_D=[np.heaviside(-energy,0) if len(energy)>0 else [] for energy in self.energyall]
        else:
            self.F_D=[1./(np.exp(np.array(energy)/self.T)+1) if len(energy)>0 else [] for energy in self.energyall]

    def get_disorder(self,pos):
        try:
            return self.disorder[pos]
        except:
            self.disorder[pos]=np.random.uniform(-self.U_D,self.U_D)
            return self.disorder[pos]

    def generate_disorder(self):
        if len(self.disorder)==0:
            for site in self.system.sites:
                self.disorder[site.pos]=np.random.uniform(-self.U_D,self.U_D)
        else:
            raise ValueError("Disorder is already generated")

    def update_Delta(self):
        if len(self.Delta_array)==0:
            for site in self.system.sites:
                self.Delta[site.pos]=self.Delta_0
        else:
            for site_index,site in enumerate(self.system.sites):
                self.Delta[site.pos]=self.Delta_array[site_index]
        
        self.Delta_mean=np.mean([delta for position,delta in self.Delta.items() if position[0]>0])
                
    def make_system(self):
        def shape(pos):
            x,y=pos
            in_x=-self.L_FM[0]<=x<=self.L_Al[0]
            in_y=0<=y<=self.L_Al[1]
            return in_x and in_y
        
        def onsite(site,kz_index):
            x,y=site.pos
            if x>self.a[0]:
                return (2*self.t_Al[0]+2*self.t_Al[1]+self.energy_Al[kz_index]-self.EF_Al+self.disorder[site.pos])*np.kron(self.sz,self.s0)+self.Delta[site.pos]*np.kron(self.sy,self.sy)
            elif x>0:
                return (self.t_Al[0]+self.t_int[0]+2*self.t_Al[1]+self.energy_Al[kz_index]-self.EF_Al+self.disorder[site.pos])*np.kron(self.sz,self.s0)+self.Delta[site.pos]*np.kron(self.sy,self.sy)
            elif x>-self.a[0]:
                return (self.t_FM[0]+self.t_int[0]+2*self.t_FM[1]+self.energy_FM[kz_index]-self.EF_FM)*np.kron(self.sz,self.s0)+self.h_exc*np.kron(self.sz,self.sz)
            else:
                return (2*self.t_FM[0]+2*self.t_FM[1]+self.energy_FM[kz_index]-self.EF_FM)*np.kron(self.sz,self.s0)+self.h_exc*np.kron(self.sz,self.sz)

        def hopping_x(site1,site2):
            x1,y1=site1.pos
            x2,y2=site2.pos
            if (x1+x2)/2>self.a[0]:
                return -self.t_Al[0]*np.kron(self.sz,self.s0)
            if (x1+x2)/2<0:
                return -self.t_FM[0]*np.kron(self.sz,self.s0)
            if 0<(x1+x2)/2<self.a[0]:
                return -self.t_int[0]*np.kron(self.sz,self.s0)
            raise ValueError("The hopping in x dimension between ({},{}) and ({},{}) is not defined".format(x1,y1,x2,y2))


        def hopping_y(site1,site2):
            x1,y1=site1.pos
            x2,y2=site2.pos
            if (x1+x2)/2>0:
                return -self.t_Al[1]*np.kron(self.sz,self.s0)
            if (x1+x2)/2<self.a[0]:
                return -self.t_FM[1]*np.kron(self.sz,self.s0)
            raise ValueError("The hopping in y dimension between ({},{}) and ({},{}) is not defined".format(x1,y1,x2,y2))


        lat=kwant.lattice.Monatomic(((self.a[0],0),(0,self.a[1])))
        syst=kwant.Builder()
        syst[lat.shape(shape,(0,0))]=onsite
        syst[kwant.builder.HoppingKind((1, 0), lat, lat)]=hopping_x
        syst[kwant.builder.HoppingKind((0, 1), lat, lat)]=hopping_y
        self.system=syst.finalized()
        return self.system

    def energyMF(self):
        start_time=time.time()
        k_req_max=np.maximum(self.estimate[:self.N_Al[2]//2],20)
        for kz_index,k_req in enumerate(k_req_max):
            k_req=np.max([40,2*(1.1*k_req//2).astype(int)])
            self.make_system()
            H_bdg=self.system.hamiltonian_submatrix(params=dict(kz_index=kz_index),sparse=True)
            H_bdg=csc_matrix(np.real((H_bdg+H_bdg.T.conj())/2))
            try:
                val,vec=sla.eigsh(H_bdg,k=k_req,sigma=0,v0=self.wfall[kz_index][:,0])
            except:
                val,vec=sla.eigsh(H_bdg,k=k_req,sigma=0)

            print('kz_index={},min(val)={:e},max(val)={:e},'.format(kz_index,np.min(np.abs(val)),np.max(np.abs(val))))
            while np.max(val)<self.E_D:
                print('k_req ({}) is too small, restart with k_req ({})'.format(k_req,k_req*2))
                k_req=k_req*2
                try:
                    val,vec=sla.eigsh(H_bdg,k=k_req,sigma=0,v0=self.wfall[kz_index][:,0])
                except:
                    val,vec=sla.eigsh(H_bdg,k=k_req,sigma=0)

            debyeindex=(np.abs(val)<=self.E_D)
            val=val[debyeindex]
            vec=vec[:,debyeindex]
            sort_index=np.argsort(val)
            val=val[sort_index]
            vec=vec[:,sort_index]
            print('len(val)={}'.format(val.shape[0]))

            self.energyall.append(val)
            self.wfall.append(vec)

        if self.store_history:
            self.energyall_history.append(self.energyall)
            self.wfall_history.append(self.wfall)
        print('Elapsed time is: {:.1f}s'.format(time.time()-start_time))

    def ave(self):
        self.Fermi_dist()
        summation=[]
        # Here only half of k_z is considered because E(kz) is even
        for kz_index,energy in enumerate(self.energyall):
            if len(energy)>0:
                wf=self.wfall[kz_index].reshape((len(self.system.sites),4,-1))
                summation.append((wf[:,3,:].conj()*wf[:,0,:]-wf[:,2,:].conj()*wf[:,1,:])@self.F_D[kz_index])

        ave=np.mean(summation,0)/(2)    
        self.Delta_array=self.g*ave
        self.update_Delta()

    def save(self):
        save_dict=self.__dict__
        del save_dict['system']
        with open("Lz{:.2f}g{:.2f}ED{:.2f}.pickle".format(self.L_Al[2],self.g,self.E_D/(433*8.617333262e-5)),"wb") as f:
            pickle.dump(save_dict,f)

def run():
    params=Params(L_Al=np.array([10,10,10]),L_FM=np.array([2,10,10]),U_D=0,Delta_0=4.8e-4)
    params.Delta_mean_list=[params.Delta_mean]
    params.energyMF()
    for i in range(1000):
        params.ave()
        print('-'*10+'Iteration: {}, Average Delta: {:e} eV'.format(i,params.Delta_mean)+'-'*10)
        params.Delta_mean_list.append(params.Delta_mean)
        if np.abs(params.Delta_mean_list[-1]-params.Delta_mean_list[-2])<1e-8:
            break
        params.energyMF()
        params.save()
    
    return params


    
# def sparse_diag(matrix, k, sigma, **kwargs):
#     """Call sla.eigsh with mumps support.

#     Please see scipy.sparse.linalg.eigsh for documentation.
#     """
#     class LuInv(sla.LinearOperator):
#         def __init__(self, A):
#             inst = mumps.MUMPSContext()
#             inst.analyze(A, ordering='pord') # other values: pord, metis, scotch
#             inst.factor(A)
#             self.solve = inst.solve
#             sla.LinearOperator.__init__(self, A.dtype, A.shape)

#         def _matvec(self, x):
#             return self.solve(x.astype(self.dtype))

#     opinv = LuInv(matrix - sigma * identity(matrix.shape[0]))
#     return sla.eigsh(matrix, k, sigma=sigma, OPinv=opinv, **kwargs)

if __name__=="__main__":
    run()
