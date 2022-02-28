import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

plt.style.use('classic')

if __name__ == '__main__':

    
    kepler_url="https://exoplanetarchive.ipac.caltech.edu/cgi-bin/nstedAPI/nph-nstedAPI?table=cumulative&where=koi_disposition+like+'CONFIRMED'&format=csv"
    df=pd.read_csv(kepler_url) 

    fig, ax = plt.subplots(1,1,figsize=(10,7))
    fig.subplots_adjust(top=0.96,right=0.97,bottom=0.14)
    #ax=fig.add_subplot(111)
    
    ax.errorbar(df['koi_period'],df['koi_prad'],yerr=[-df['koi_prad_err2'].values,df['koi_prad_err1'].values],
                xerr=[-df['koi_period_err2'].values,df['koi_period_err1'].values],fmt='s',lw=1.2,ms=3.0,capsize=3,alpha=0.5,mec=u'#1f77b4',color=u'#1f77b4')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$P_{\rm p}[{\rm d}]$',size=26)
    ax.set_ylabel(r'$R_{\rm p}[R_{\oplus}]$',size=26)
    ax.tick_params(axis='both',width=2,length=4,labelsize=18)
    ax.tick_params(axis='both',which='both',width=2)

    fig.savefig('kepler_period-radius.png')
    #plt.show()

    index = df['koi_srad'] > 4
    print(df[index].shape)
    print(df[index]['koi_srad'])
    print(df[index])
    for item in df[index]:
        print(item)
