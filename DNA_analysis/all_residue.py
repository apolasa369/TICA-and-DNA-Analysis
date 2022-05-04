import sys
import numpy as np
import mdtraj as mdt
import matplotlib
#matplotlib.rcParams['font.size'] = 9
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import matplotlib.pyplot as plt
from matplotlib import cm

pdbname = sys.argv[1]
fname = sys.argv[2]
name_split = pdbname.split('.')

pdb = mdt.load_pdb(str(pdbname))
corr = np.loadtxt(fname)
prt_name = name_split[0].split('/') 
#print(name_split)
print(prt_name[-1])
rep = int(name_split[1])
chain_name = name_split[1]
label_dict = {"6al2": {0:"YidC Set 1",1:"YidC Set 2"}
             ,"6al2TM": {0:"YidC ΔPD Set 1",1:"YidC ΔPD Set 2"}
             ,"6al2-woloop": {0:"YidC ΔC2 Set 1",1:"YidC ΔC2 Set 2"}
             ,"6al2-woloopTM": {0:"YidC ΔC2 ΔPD Set 1",1:"YidC ΔC2 ΔPD Set 2"}
            }
#p_dict = {"5x58":{"p1":"Protomer A","p2":"Protomer B","p3":"Protomer C"}
 #         ,"5x5b":{"p1":"Protomer A","p2":"Protomer B","p3":"Protomer C"}
  #        ,"6vxx":{"p1":"Protomer A","p2":"Protomer C","p3":"Protomer B"}
   #       ,"6vyb":{"p1":"Protomer A","p2":"Protomer C","p3":"Protomer B"}
    #     }



res_dict = {"6al2":{"PD":{"start":49,"stop":326,"min":49,"max":326,"color":"r","txt":"mid"}},
        "6al2-woloop":{"PD":{"start":49,"stop":326,"min":49,"max":326,"color":"r","txt":"mid"}},
    }

res_common = [int(str(p)[3:]) for p in pdb.topology.residues]

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
    new_map = cm.PuOr #[(cm.Spectral(i)) for i in range(0,250)]
    #colors = [('white')] + [(cm.coolwarm(i)) for i in range(0,250)]

    #new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=-1, vmax=1)

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
    title1 = plt.title(title, fontsize=30,y=1.08)
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
    new_map = cm.PuOr #[(cm.Spectral(i)) for i in range(0,250)]
    #colors = [('white')] + [(cm.coolwarm(i)) for i in range(0,250)]

    #new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=-1, vmax=1)

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
    
    ax.set_xlim(res_min,res_max)
    ax.set_ylim(res_min,res_max)
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
    new_map = cm.PuOr #[(cm.Spectral(i)) for i in range(0,250)]
    #colors = [('white')] + [(cm.coolwarm(i)) for i in range(0,250)]

    #new_map = matplotlib.colors.LinearSegmentedColormap.from_list('new_map', colors, N=300)
    heatmap = ax.pcolor(Z, cmap=new_map, vmin=-1, vmax=1)

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

    ax.set_xlim(res_min,res_max)
    ax.set_ylim(res_min,res_max)
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

plot_map(corr
         ,res_common
         ,res_dict[prt_name[-1]]
         ,f"{label_dict[prt_name[-1]][rep]}"
         ,f"{fname[:-4]}")
#plot_map_nolabel(corr
#         ,res_common
#         ,res_dict[prt_name[-1]]
#         ,f"{label_dict[prt_name[-1]][rep]} {p_dict[prt_name[-1]][chain_name]}"
#         ,f"{fname[:-4]}")
#plot_map_label(corr
#         ,res_common
#         ,res_dict[prt_name[-1]]
#         ,f"{label_dict[prt_name[-1]][rep]} {p_dict[prt_name[-1]][chain_name]}"
#         ,f"{fname[:-4]}")
