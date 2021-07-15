import numpy as np
import matplotlib.pyplot as plt

np.random.seed(42)

if __name__ == '__main__':

    N = 100
    mu1,s1 = 5.0, 3.0
    data1 = (np.random.randn(N) *s1) + mu1

    N = 200
    mu2,s2 = 9.0, 0.5
    data2 = (np.random.randn(N) *s2) + mu2

    #plt.hist(data1,bins=8)
    #plt.hist(data2,bins=8)
    np.savetxt('dataset1.txt',data1)
    np.savetxt('dataset2.txt',data2)

    trueX =5
    fig, ax = plt.subplots(figsize=(8,5))
    ax.errorbar([1],data1.mean(),data1.std(),marker='s',capsize=4,
                label='Experiment 1')
    ax.errorbar([2],data2.mean(),data2.std(),marker='s',capsize=4,
                label='Experiment 2')

    ax.plot([0,5],[trueX,trueX],ls=':',zorder=0,color='k')
    ax.set_xticks([1,2])
    ax.set_xticklabels([])
    ax.tick_params(axis='both',direction='in',labelsize=14)
    ax.set_xlim(0,4)
    ax.set_ylabel('Measured '+r'$X$',size=16)
    ax.legend(frameon=False,loc='lower right')
    plt.savefig('experimental_data.png')
    plt.show()
