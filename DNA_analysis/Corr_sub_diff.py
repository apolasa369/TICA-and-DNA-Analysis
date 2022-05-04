import sys
from itertools import combinations
import mdtraj as md
import numpy as np
import matplotlib
#matplotlib.rcParams['font.size'] = 9
#matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = ['Times New Roman']
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib import cm

label_dict = {"6al2": {0:"YidC Set 1",1:"YidC Set 2"}
             ,"6al2TM": {0:"YidC ΔPD Set 1",1:"YidC ΔPD Set 2"}
             ,"6al2-woloop": {0:"YidC ΔC2 Set 1",1:"YidC ΔC2 Set 2"}
             ,"6al2-woloopTM": {0:"YidC ΔC2 ΔPD Set 1",1:"YidC ΔC2 ΔPD Set 2"}
            }
def plot_map(correlation,res_list,r_dict, title, output_prefix):
    #M = np.ma.masked_invalid(correlation)
    M = np.array(correlation)
    res = np.array(res_list,dtype=int) 
    res_max = res.max()
    Z = np.empty((res_max+1,res_max+1))
    Z.fill(np.nan)  
    xx,yy = np.meshgrid(res,res)
    Z[xx,yy]=M
    ax = plt.subplots()[1]
    new_map = cm.summer_r #[(cm.Spectral(i)) for i in range(0,250)]
    #colors = [('white')] + [(cm.coolwarm(i)) for i in range(0,250)]

    #new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=0, vmax=2)

    fig = plt.gcf()
    fig.set_size_inches(6,5.4)
    fig.subplots_adjust(top=0.831,
                    bottom=0.152,
                    left=0.121,
                    right=0.977,
                    hspace=0.2,
                    wspace=0.2)
    ax.set_frame_on(False)
    ax.grid(False)

    plt.xticks(rotation=90)

    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)

    for t in ax.yaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)
    
    mn = np.nan
    for di in r_dict:
        res_start = r_dict[di]['start']
        res_stop = r_dict[di]['stop']
        res_min = r_dict[di]['min'] 
        mn = min(res_min,mn)
        res_max = r_dict[di]['max']
        res_mid = (res_start + res_stop)/2
        if (r_dict[di]['txt'] == "mid"):
            txt_mid = res_mid
        else:
            txt_mid = r_dict[di]['txt']
        res_len = res_stop - res_start
        cl = r_dict[di]['color']
        ax.annotate(di,xy=(res_mid,res_max-8),xytext=(txt_mid,res_max+8),xycoords="data"
                    ,ha="center",va="center"
                    ,fontsize=8
                    ,c=cl
                    ,arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            facecolor=cl)
                    ,bbox=dict(fc='none',ec=cl)
                   )
        ax.annotate(di,xy=(res_max-8,res_mid),xytext=(res_max+8,txt_mid),xycoords="data"
                    ,ha="center",va="center"
                    ,c=cl
                    ,fontsize=8
                    ,rotation=-90
                    ,arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3")
                    ,bbox=dict(fc='none',ec=cl)
                   )
        recx = matplotlib.patches.Rectangle((res_start,res_min),res_len,res_max-mn
                                           ,lw=1
                                           ,fc='none',ec=cl)
        recy = matplotlib.patches.Rectangle((res_min,res_start),res_max-mn,res_len
                                           ,lw=1
                                           ,fc='none',ec=cl)
        ax.add_patch(recx)
        ax.add_patch(recy)

    ax.set_xlim(res_min,res_max)
    ax.set_ylim(res_min,res_max)
    title1 = plt.title(title, fontsize=20,y=1.08)
    xl1 = plt.xlabel('Residue Index', fontsize=24)
    xl2 = plt.ylabel("Residue Index", fontsize=24)

    cbar = plt.colorbar(heatmap, orientation="vertical")
    plt.tight_layout()

    #plt.show()
    plt.savefig('%s.pdf' % output_prefix, dpi=300,bbox_inches='tight')
    #plt.savefig('%s.eps' % output_prefix, dpi=300,bbox_inches='tight')
    plt.savefig('%s.png' % output_prefix, dpi=300,bbox_inches='tight')
    plt.close('all')


def plot_map_nolabel(correlation,res_list,r_dict, title, output_prefix):
    #M = np.ma.masked_invalid(correlation)
    M = np.array(correlation)
    res = np.array(res_list,dtype=int) 
    res_max = res.max()
    Z = np.empty((res_max+1,res_max+1))
    Z.fill(np.nan)  
    xx,yy = np.meshgrid(res,res)
    Z[xx,yy]=M
    ax = plt.subplots()[1]
    
    new_map = cm.Purples
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=0., vmax=1.35)
    fig = plt.gcf()
    fig.set_size_inches(6,5.4)
    fig.subplots_adjust(top=0.831,
                    bottom=0.152,
                    left=0.121,
                    right=0.977,
                    hspace=0.2,
                    wspace=0.2)
    ax.set_frame_on(False)
    ax.grid(False)

    plt.xticks(rotation=90)

    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)

    for t in ax.yaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)
    
    ax.set_xlim(0,res_max+40)
    ax.set_ylim(0,res_max+40)
    title1 = plt.title(title, fontsize=30,y=1.08)
    xl1 = plt.xlabel('Residue Index', fontsize=24)
    xl2 = plt.ylabel("Residue Index", fontsize=24)

    cbar = plt.colorbar(heatmap, orientation="vertical")
    plt.tight_layout()

    # Turn heatmap on, labels and titles off and plot only heatmap
    #heatmap.set_visible(True)
    plt.axis('off')
    cbar.remove()
    #cbar.outline.set_visible(False)
    title1.set_visible(False)
    #plt.show()
    plt.savefig('%s.png' % (output_prefix+'_hm')
                ,dpi=300
                ,transparent=True
    #            ,bbox_inches='tight'
               )
    plt.close('all')

def plot_map_label(correlation,res_list,r_dict, title, output_prefix):
    #M = np.ma.masked_invalid(correlation)
    M = np.array(correlation)
    res = np.array(res_list,dtype=int) 
    res_max = res.max()
    Z = np.empty((res_max+1,res_max+1))
    Z.fill(np.nan)  
    xx,yy = np.meshgrid(res,res)
    Z[xx,yy]=M
    ax = plt.subplots()[1]
    new_map = cm.Purples
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=0., vmax=1.35)

    fig = plt.gcf()
    fig.set_size_inches(6,5.4)
    fig.subplots_adjust(top=0.831,
                    bottom=0.152,
                    left=0.121,
                    right=0.977,
                    hspace=0.2,
                    wspace=0.2)
    ax.set_frame_on(False)
    ax.grid(False)

    plt.xticks(rotation=90)

    # Turn off all the ticks
    ax = plt.gca()
    for t in ax.xaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)

    for t in ax.yaxis.get_major_ticks():
        t.tick1line.set_visible = False
        t.tick2line.set_visible = False
        t.label.set_fontsize(16)
    
    mn = np.nan
    for di in r_dict:
        res_start = r_dict[di]['start']
        res_stop = r_dict[di]['stop']
        res_min = r_dict[di]['min'] 
        mn = min(res_min,mn)
        res_max = r_dict[di]['max']
        res_mid = (res_start + res_stop)/2
        if (r_dict[di]['txt'] == "mid"):
            txt_mid = res_mid
        else:
            txt_mid = r_dict[di]['txt']
        res_len = res_stop - res_start
        cl = r_dict[di]['color']
        ax.annotate(di,xy=(res_mid,res_max-90),xytext=(txt_mid,res_max+50),xycoords="data"
                    ,ha="center",va="center"
                    ,fontsize=10
                    ,c=cl
                    ,arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3",
                            facecolor=cl)
                    ,bbox=dict(fc='none',ec=cl)
                   )
        ax.annotate(di,xy=(res_max-90,res_mid),xytext=(res_max+40,txt_mid),xycoords="data"
                    ,ha="center",va="center"
                    ,c=cl
                    ,fontsize=10
                    ,rotation=-90
                    ,arrowprops=dict(arrowstyle="->",
                            connectionstyle="arc3")
                    ,bbox=dict(fc='none',ec=cl)
                   )
        recx = matplotlib.patches.Rectangle((res_start,res_min),res_len,res_max-mn
                                           ,lw=1
                                           ,fc='none',ec=cl)
        recy = matplotlib.patches.Rectangle((res_min,res_start),res_max-mn,res_len
                                           ,lw=1
                                           ,fc='none',ec=cl)
        ax.add_patch(recx)
        ax.add_patch(recy)

    ax.set_xlim(0,res_max+40)
    ax.set_ylim(0,res_max+40)
    title1 = plt.title(title, fontsize=30,y=1.08)
    xl1 = plt.xlabel('Residue Index', fontsize=24)
    xl2 = plt.ylabel("Residue Index", fontsize=24)

    cbar = plt.colorbar(heatmap, orientation="vertical")
    plt.tight_layout()

    # Turn off heatmap and plot labels
    heatmap.set_visible(False)
    #cbar.solids.set_visible(False)
    #plt.show()
    plt.savefig('%s.pdf' % (output_prefix+"_labels")
                ,dpi=300
                ,bbox_inches='tight'
                ,transparent=True)
    plt.close('all')

res_dict = {"6al2":{"PD":{"start":49,"stop":327,"min":49,"max":326,"color":"g","txt":"mid"}},
        "6al2-woloop":{"PD":{"start":49,"stop":327,"min":49,"max":326,"color":"g","txt":"mid"}},
         }
pdb = sys.argv[1]
pdb2 = sys.argv[2]
p1f = sys.argv[3]
p2f = sys.argv[4]
#title = sys.argv[5]
#out_pre = sys.argv[6]
pdbname = sys.argv[1]
pdbname2 = sys.argv[2]
fname = sys.argv[3]
fname2 = sys.argv[4]
name_split = pdbname.split('.')
name_split2 = pdbname2.split('.')

#pdb = mdt.load_pdb(str(pdbname))
corr = np.loadtxt(fname)
prt_name = name_split[0].split('/')
prt_name2 = name_split2[0].split('/') 
rep = int(name_split[1])
rep2 = int(name_split2[1])
chain_name = name_split[1]

P1 = np.loadtxt(p1f)
P2 = np.loadtxt(p2f)
#P3 = np.loadtxt(p3f)
topo = md.load(pdb)
topo2 = md.load(pdb2)
res_list = [int(str(topo.topology.residue(i))[3:]) for i in range(topo.n_residues)]
res_list2 = [int(str(topo2.topology.residue(i))[3:]) for i in range(topo2.n_residues)]
res = np.array(res_list,dtype=int)
res2 = np.array(res_list2,dtype=int) 
res_max = res.max()
res_min = res.min()
res_max2 = res2.max()
res_min2 = res2.min()

Z = np.empty((res_max+1,res_max+1))
Z2 = np.empty((res_max2+1,res_max2+1))
#Z = np.empty_like(P1)
Z.fill(np.nan)  
Z2.fill(np.nan)
xx,yy = np.meshgrid(res,res)
xx2,yy2 = np.meshgrid(res2,res2)
Z[xx,yy]=P1
Z2[xx2,yy2]=P2


#for x,y in combinations((P1,P2,P3),2):
sub = Z-Z2
P3 = np.abs(sub)
P3 = P3[49:,49:]
#title_dict = {"6vyb":,"6vxx","5x5b","5x58"}
plot_map(P3,res_list,res_dict[pdbname[:4]],f"{label_dict[prt_name[-1]][rep]}_{label_dict[prt_name2[-1]][rep2]}",f"{fname[:-4]}-{fname2[:-4]}")
#plot_map_label(P3,res_list,res_dict[title[:4]],title,out_pre)
#plot_map_nolabel(P3,res_list,res_dict[title[:4]],title,out_pre)
