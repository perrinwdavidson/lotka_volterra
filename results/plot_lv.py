# import packages -----------------------------------------
import numpy as np
import matplotlib.pyplot as plt

# define plotting function --------------------------------
def plotlk(filenamex, filenamey, filenameargs, ichange, titlename):
    '''
    plot lotka-volterra outputs
    ---------------------------
    inputs:
        [1] filenamex - string filename for prey data
        [2] filenamey - string filename for predator data
        [3] filenameargs - argument names for change
        [4] titlename - title of the graph
    outputs:
        [1] plt - plot of values
    '''
    # read in data:
    lk_x = np.loadtxt(filenamex)
    lk_y = np.loadtxt(filenamey)
    argnames = np.loadtxt(filenameargs)

    # get dimensions:
    numarg = np.shape(argnames)[0]
    numparams = np.shape(argnames)[1]

    # set variable change:
    varnames = np.array(['$x_0 =$ ', '$y_0 =$ ', '$\\alpha =$ ', '$\\beta =$ ', '$\\delta =$ ', '$\\gamma =$ '])

    # get constant parameters:
    constparams = np.delete(argnames[0, :], ichange)
    constparamnames = np.delete(varnames, ichange)

    # make subplot title:
    subtitle = ''
    for i in range(numparams - 1):
        subtitle0 = constparamnames[i] + str(round(constparams[i], 2)) + ', '
        subtitle = subtitle + subtitle0
    subtitle = subtitle[:len(subtitle)-2]

    # Plot:
    for i in range(numarg):
        plt.plot(lk_x[:, i], lk_y[:, i], '-', label=varnames[ichange] + str(argnames[i, ichange]))
    plt.xlabel('Prey')
    plt.ylabel('Predator')
    plt.suptitle(titlename, fontsize=16)
    plt.title(subtitle, fontsize=12)
    plt.legend()
    plt.savefig('plots/lk.png', dpi=300)
    plt.show()
    plt.close()

# plot ----------------------------------------------------
if __name__ == '__main__':

    # Set basepath:
    basepath = '/Users/perrindavidson/Research/uchicago/current/lotka_volterra/'
    directorypath = 'run/'

    # Set inputs:
    filenamex = basepath + directorypath + 'lk_prey.txt'
    filenamey = basepath + directorypath + 'lk_pred.txt'
    filenameargs = basepath + directorypath + 'args.txt'
    ichange = 1
    titlename = 'Lotka-Volterra Dynamics'

    # Plot:
    plotlk(filenamex, filenamey, filenameargs, ichange, titlename)
