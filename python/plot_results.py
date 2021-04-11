import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import loadmat
from scipy.interpolate import interp1d
import pickle


def show_plot_pickle(data,export=False):
    fig,ax=plt.subplots(2,2)
    ax[0,0].plot(data['Delta_mean_list'])
    ax[0,0].text(.5,1,'mean delta = {:f}meV'.format(1000*data['Delta_mean_list'][-1]),ha='center',va='bottom',transform=ax[0,0].transAxes)
    # ax[0,1].plot(data['htotlist'][0])
    gs = ax[1, 0].get_gridspec()
    ax[1,0].remove()
    ax[1,1].remove()
    axbig = fig.add_subplot(gs[1,:])
    axbig.scatter(np.array(list(data['Delta'].keys()))[:,0]/0.05076,np.array(list(data['Delta'].keys()))[:,1]/0.05076,c=np.array(list(data['Delta'].values())))
    ax[0,0].text(0,1,'Al:[{},{},{}]'.format(data['N_Al'][0],data['N_Al'][1],data['N_Al'][2]),transform=ax[0,0].transAxes,ha='right',va='bottom')
    if export:
        fig.savefig('Nkz{}.pdf'.format(data['N_Al'][2]),bbox_inches='tight')
    return 1000*data['Delta_mean_list'][-1]


if __name__=="__main__":
    i='Lz70.000g5.60ED1.00.pickle'
    with open('{}'.format(i),'rb') as f:
        z=pickle.load(f)
    data_fn.append(z)
    dim_list.append(z['N_Al'])
    delta_list.append(show_plot_pickle(z,True))
