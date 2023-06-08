import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
import math
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import statistics

def get_histograms(type):
    RMSDs = []
    for i in range(1, 10):
        print('loading run ' + str(i))
        temp_RMSD=[]
        try:
            print(type)
            u = mda.Universe('/nfs/hfdata/epstein/CALHM1/analysis-results/NTH-atomistic/'+type + '/' + 'minim-prot.gro','/nfs/hfdata/epstein/CALHM1/analysis-results/NTH-atomistic/'+ type + '/' + f"r{i}" + '_skip10_prot.xtc')
            print(u)
            ref = u.select_atoms('protein and name CA')
            ref_u = mda.Universe('/nfs/hfdata/epstein/CALHM1/analysis-results/NTH-atomistic/'+type + '/' +'minim-prot.gro')
            aligner = align.AlignTraj(u, ref, select='protein and name CA', in_memory=True).run()
            print('The number of frames are: ' + str(len(u.trajectory)))
            rmsd_array = []
            for ts in u.trajectory[1:100]:
            # for ts in u.trajectory:
                ref = ref_u.select_atoms("protein and name CA")
                mob = u.select_atoms("protein and name CA")
                ligand_rmsd = rmsd(mob.positions, ref.positions)
                rmsd_array.append(ligand_rmsd)
            # temp_RMSD.append(rmsd_array)
            temp_RMSD.append(rmsd_array)
        except Exception as e:
            print((f"An error occurred while loading run {i}: {e}. Skipping to the next run..."))
            continue
        # print(temp_RMSD)
        temp_RMSD=np.array(temp_RMSD)
        temp_RMSD.flatten()
        RMSDs.append(temp_RMSD)
    RMSDs=np.array(RMSDs)
    RMSDs.flatten()

    reordered = []
    for i in range(len(RMSDs)):
        print(i)
        for j in RMSDs[i]:
            for k in j:
                reordered.append(k)
    return(reordered)

apo_hist=get_histograms(type='apo')
chol_hist=get_histograms(type='chol')
popc_hist=get_histograms(type='popc')

print('apo median = ' + str(statistics.median(apo_hist)))
print('chol median = ' + str(statistics.median(chol_hist)))
print('popc median = ' + str(statistics.median(popc_hist)))

print('apo mean = '+ str(statistics.mean(apo_hist)))
print('chol mean = '+ str(statistics.mean(chol_hist)))
print('popc mean = '+ str(statistics.mean(popc_hist)))

plt.figure()
plt.hist(apo_hist,label='apo',color='r',bins=50,histtype='step',density=True)
plt.hist(chol_hist,label='chol',color='g',bins=50,histtype='step',density=True)
plt.hist(popc_hist,label='popc',color='b',bins=50,histtype='step',density=True)
plt.legend(loc='upper left')
plt.xlim(4.5,7.0)
# plt.xlim(0,101)
plt.axvline(statistics.mean(apo_hist),ymin=0,ymax=1,color='r',linestyle='dashed')
plt.axvline(statistics.mean(chol_hist),ymin=0,ymax=1,color='g',linestyle='dashed')
plt.axvline(statistics.mean(popc_hist),ymin=0,ymax=1,color='b',linestyle='dashed')


plt.savefig('RMSD-histograms.png')
plt.savefig('RMSD-histograms.eps',format="eps",dpi=300)
plt.savefig('RMSD-histograms.svg',format="svg",dpi=300)
plt.show()
