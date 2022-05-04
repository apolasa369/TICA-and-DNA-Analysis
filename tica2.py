import matplotlib.pyplot as plt
import numpy as np
import mdshare # pip install mdshare
import pyemma # pip install pyemma
import sys
from matplotlib import cm

#pdb = mdshare.fetch('alanine-dipeptide-nowater.pdb', working_directory='data')
#files = mdshare.fetch('alanine-dipeptide-*-250ns-nowater.xtc', working_directory='data')
pdb = sys.argv[1] #required pdb file of protein system can be just CA
files = sys.argv[2] # trajectory (dcd) file or xtc file of trajectory can be used
fname = sys.argv[2]
print(pdb)
print(files)
name_split = pdb.split('.')
print(name_split)
prt_name = name_split[0].split('/') 
#print(name_split)
print(prt_name[-1])
rep = int(name_split[1])
chain_name = name_split[1]
label_dict = {"6al2": {0:"YidC Set 1",1:"YidC Set 2"} # Please modify lables respective to your systems
             ,"6al2TM": {0:"YidC ΔPD Set 1",1:"YidC ΔPD Set 2"}
             ,"6al2-woloop": {0:"YidC ΔC2 Set 1",1:"YidC ΔC2 Set 2"}
             ,"6al2-woloopTM": {0:"YidC ΔC2 ΔPD Set 1",1:"YidC ΔC2 ΔPD Set 2"}
            }

feat = pyemma.coordinates.featurizer(pdb)
feat.add_selection(feat.select_Ca())
#feat.add_selection(feat.select_Heavy())
#feat.add_backbone_torsions(periodic=False)
#feat.add_sidechain_torsions(which='chi1', cossin=True, periodic=False)
#feat.add_distances_ca(periodic=False)
# Feature - pairs
#pairs = feat.pairs(feat.select_Heavy())
#feat.add_distances(pairs, periodic=False)

data = pyemma.coordinates.load(files, features=feat)

print('We have {} features.'.format(feat.dimension()))

#pca = pyemma.coordinates.pca(data, dim=2)
tica = pyemma.coordinates.tica(data, lag=3, dim=2)

#pca_concatenated = np.concatenate(pca.get_output())
tica_concatenated = np.concatenate(tica.get_output())

#np.savetxt('%s.txt', % output_prefix, tica_concatenated, fmt='%1.3f', newline=", ")

#cls_pca = pyemma.coordinates.cluster_kmeans(pca, k=100, max_iter=50, stride=10)
#cls_tica = pyemma.coordinates.cluster_kmeans(tica, k=3, max_iter=50, stride=10)

#its_pca = pyemma.msm.its(
 #   cls_pca.dtrajs, lags=[1, 2, 5, 10, 20, 50], nits=4, errors='bayes')
#its_tica = pyemma.msm.its(
 #   cls_tica.dtrajs, lags=[1, 2, 5, 10, 20, 50], nits=4, errors='bayes')
def plot_map(title, output_prefix):
   fig, axes = plt.subplots(1, 1, figsize=(6, 6))
   title1 = plt.title(title, fontsize=30,y=1.08)
 #  pyemma.plots.plot_feature_histograms(pca_concatenated, ax=axes[0, 0])
 #  pyemma.plots.plot_feature_histograms(tica_concatenated, ax=axes[1, 0])
 #  axes[0, 0].set_title('PCA')
 #  axes[1, 0].set_title('TICA')
 #  pyemma.plots.plot_density(*pca_concatenated.T, ax=axes[0, 1], cbar=False, alpha=0.1)
 #  axes[0, 1].scatter(*cls_pca.clustercenters.T, s=15, c='C1')
 #  axes[0, 1].set_xlabel('PC 1')
 #  axes[0, 1].set_ylabel('PC 2')
   new_map = cm.summer_r
   plt.plot(*tica_concatenated.T, ".",alpha=0.1)
   #axes.scatter(*cls_tica.clustercenters.T, s=15, c='C1')
   axes.set_xlabel('IC 1')
   axes.set_ylabel('IC 2')
 #  pyemma.plots.plot_implied_timescales(its_pca, ax=axes[0, 2], units='ps')
 #  pyemma.plots.plot_implied_timescales(its_tica, ax=axes[1, 2], units='ps')
 #  axes[0, 2].set_ylim(1, 2000)
 #  axes[1, 2].set_ylim(1, 2000)
   fig.tight_layout()

   plt.savefig('%s.pdf' % output_prefix, dpi=300,bbox_inches='tight')
   plt.close('all')

plot_map(f"{label_dict[prt_name[-1]][rep]}",f"{fname[:-4]}")
np.savetxt(f"{fname[:-4]}.txt", tica_concatenated, fmt='%1.3f', newline="\n")
