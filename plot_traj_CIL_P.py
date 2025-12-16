#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 08:46:23 2024

@author: atgrigas
"""

import glob
import numpy as np
import os
import sys
import math
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

pon = '0.01'
poff = '0.001'
target_phi = '0.4'
target_temp = '0'
gamma = '1000'
# %%
theta_cutoff = '2.0'

rr = '0.1'

seq = '0'
run = '0'

cutoff = 0

num_atoms = 256

frame_split = 10

file = 'traj_data/traj_seq_'+seq+'_'+run+'_pon_'+pon+'_poff_'+poff+'_phi'+target_phi+'_T_'+target_temp+'_gamma_'+gamma+'_theta_cutoff_'+theta_cutoff+'_rr_'+rr+'.txt'

sourcelines = np.loadtxt(file)
 
#%%       
 
Lx_list = []
Ly_list = []
coords_list = []
frame = 0
for i in range(0,len(sourcelines)):
    
    if i % (num_atoms + 1) == 0:
        Lx = sourcelines[i,0]
        Ly = sourcelines[i,1]
        
        coords = sourcelines[i+1:i+1+num_atoms,:2]
        avg = np.average(coords,axis=0)
        coords -= avg
        coords[:,0] = coords[:,0] - Lx*np.round(coords[:,0] / Lx)
        coords[:,1] = coords[:,1] - Ly*np.round(coords[:,1] / Ly)

        coords_list.append(coords)
        Lx_list.append(Lx)
        Ly_list.append(Ly)
    

    
#%% 

file = 'seqs/seq_'+seq+'.txt'

sigma_i_array =  np.loadtxt(file)


bonded_list_list = []
file = 'bonded_data/bonded_seq_'+seq+'_'+run+'_pon_'+pon+'_poff_'+poff+'_phi'+target_phi+'_T_'+target_temp+'_gamma_'+gamma+'_theta_cutoff_'+theta_cutoff+'_rr_'+rr+'.txt'
sourcelines = []
with open(file,'r') as f:
    for line in f:
        sourcelines.append(line)
count = 0
for line in sourcelines:
    bonded_list = []
    #if count < len(tension_list_list):
        #tension_list = tension_list_list[count]
    count += 1
    splits = line.split()
    bond = []
    
    bound_count = 0
    for i in range(0,len(splits)):
        
        if i % 2 == 0:
            
            
            bond = [int(splits[i]),int(splits[i+1])]
            bonded_list.append(bond)
                
            bound_count += 1

    bonded_list_list.append(bonded_list)




#%%

# Visualizing frames as singular mpl plots, can be viewed in an IDE like Spyder, or saved and loaded by imageJ

pbc_list = [[0,0],
            [0,1],
            [1,0],
            [1,1],
            [1,-1],
            [-1,1],
            [-1,0],
            [0,-1],
            [-1,-1]]

box = Lx_list[0]

rmsd_list = []
for idx in np.arange(0,len(coords_list),10):
    idx = int(idx)
    
    coords = coords_list[idx]
    Lx = Lx_list[idx]
    Ly = Ly_list[idx]
    bonded_list = bonded_list_list[idx]   
    #tension_list = tension_list_list[idx]    
    box = Lx
        
    fig,ax = plt.subplots(figsize=[7,7])     
    cmap = matplotlib.colormaps['winter']
    ax.set_xlim((-1.2*box/2,1.2*box/2))
    ax.set_ylim((-1.2*box/2,1.2*box/2)) 
    min_sigma = min(sigma_i_array)
    max_sigma = max(sigma_i_array)
    colors = []
    s_list = []
    s_list_alpha = []
    for j in range(0, num_atoms):
        r = sigma_i_array[j]
        r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
        marker_size = (1.*r_)**2
        s_list.append(marker_size)
        colors.append(cmap((r-min_sigma)/(max_sigma-min_sigma))) 
    
    for j in range(0,len(pbc_list)):
        pbc = pbc_list[j]
        
        x = coords[:,0] + Lx*pbc[0]
        y = coords[:,1] + Ly*pbc[1]
        
        if j == 0:
            plt.scatter(x, y, marker='o', facecolor=colors, edgecolor='k', linewidth=0.5, s=s_list)    
            #plt.quiver(x,y,arrows[:,0],arrows[:,1],color='r',zorder=300)
        else:
            plt.scatter(x, y, marker='o', facecolor='w', edgecolor='k', linewidth=0.5, s=s_list)    
    
        cmap = matplotlib.colormaps['bwr_r']
        norm = matplotlib.colors.Normalize(vmin=-0.6, vmax=0.6)
        
        cmap = matplotlib.colormaps['viridis_r']
        norm = matplotlib.colors.Normalize(vmin=-0.4, vmax=0.)
        xplot = []
        yplot = []
        for k in range(0,len(bonded_list)):
            bond = bonded_list[k]
            
            i = bond[0]
            j = bond[1]
            
            xi = x[i]
            xj = x[j]
            delta_x = xi-xj
            delta_x -= Lx * np.round(delta_x / Lx)
            
            yi = y[i]
            yj = y[j]
            delta_y = yi-yj
            delta_y -=  Ly * np.round(delta_y / Ly)
    
            xj = xi - delta_x
            yj = yi - delta_y        
    
        
            xplot.append(xi)
            xplot.append(xj)
            xplot.append(0)
            yplot.append(yi)
            yplot.append(yj)
            yplot.append(1)
            
            #tension = tension_list[k]
            
            #plt.plot([xi,xj],[yi,yj],c=cmap(norm(tension)),linestyle='-',zorder=200)
        
        mask = np.ma.array(yplot)
        for j in np.arange(2,len(xplot),3):
            mask[j] = np.ma.masked
    
        plt.plot(xplot,mask,color='grey',linestyle='-',zorder=200)
    
    
    plt.plot([-Lx/2,Lx/2],[-Ly/2,-Ly/2],color='k')
    plt.plot([-Lx/2,Lx/2],[Ly/2,Ly/2],color='k')
    
    plt.plot([-Lx/2,-Lx/2],[-Ly/2,Ly/2],color='k')
    plt.plot([Lx/2,Lx/2],[-Ly/2,Ly/2],color='k')
    
    
    ax.set_aspect('equal', adjustable='box')
    plt.axis('off')
   # plt.savefig('plot/plot_'+str(idx)+'.png',dpi=300)
    plt.show()
    plt.close()



#%%

# Rendering an mp4 directly


# =============================================================================
# from matplotlib.animation import FuncAnimation, FFMpegWriter
# 
# pbc_list = [[0,0],
#             [0,1],
#             [1,0],
#             [1,1],
#             [1,-1],
#             [-1,1],
#             [-1,0],
#             [0,-1],
#             [-1,-1]]
# 
# box = Lx_list[0]
# 
# print('Preparing PBC data...')
# pbc_coords_list = []
# pbc_xplot_list = []
# pbc_mask_list = []
# for idx in np.arange(0,int(len(coords_list)/5),10):
#     
#     pbc_coords = np.zeros((num_atoms*len(pbc_list),2))
#     pbc_arrows = np.zeros((num_atoms*len(pbc_list),2))
#     Lx = Lx_list[idx]
#     Ly = Ly_list[idx]
#     bonded_list = bonded_list_list[idx]        
#     
#     xplot_list = []
#     yplot_list = []
#     for pbc_idx in range(0,len(pbc_list)):
#         coords = np.copy(coords_list[idx])
#         pbc = pbc_list[pbc_idx]
#         
#         coords[:,0] = coords[:,0] + Lx*pbc[0]
#         coords[:,1] = coords[:,1] + Ly*pbc[1]
#         
#         pbc_coords[pbc_idx*num_atoms:(pbc_idx+1)*num_atoms,:] = coords
#         
#         x = np.copy(coords[:,0])
#         y = np.copy(coords[:,1])
#         xplot = []
#         yplot = []
#         for bond in bonded_list:
#             
#             i = bond[0]
#             j = bond[1]
#             
#             xi = x[i]
#             xj = x[j]
#             delta_x = xi-xj
#             delta_x -= Lx * np.round(delta_x / Lx)
#             
#             yi = y[i]
#             yj = y[j]
#             delta_y = yi-yj
#             delta_y -=  Ly * np.round(delta_y / Ly)
#     
#             xj = xi - delta_x
#             yj = yi - delta_y        
#     
#             xplot.append(xi)
#             xplot.append(xj)
#             xplot.append(0)
#             yplot.append(yi)
#             yplot.append(yj)
#             yplot.append(1)
#     
#         xplot_list += xplot
#         yplot_list += yplot
#     pbc_coords_list.append(pbc_coords)
#     mask = np.ma.array(yplot_list)
#     for j in np.arange(2,len(xplot_list),3):
#         mask[j] = np.ma.masked
#     pbc_xplot_list.append(xplot_list)
#     pbc_mask_list.append(mask)
# 
# #%%
# 
# fig,ax = plt.subplots(figsize=[7,7])     
# cmap = matplotlib.colormaps['winter']
# ax.set_xlim((-1.2*box/2,1.2*box/2))
# ax.set_ylim((-1.2*box/2,1.2*box/2))
# plt.tight_layout()
# 
# plt.close()
# 
# sigma_i_array_pbc = []
# for pbc_idx in range(0,len(pbc_list)):
#     sigma_i_array_pbc += list(sigma_i_array)
# 
# min_sigma = min(sigma_i_array_pbc)
# max_sigma = max(sigma_i_array_pbc)
# colors = []
# s_list = []
# s_list_alpha = []
# for j in range(0, len(sigma_i_array_pbc)):
#     r = sigma_i_array_pbc[j]
#     r_ = ax.transData.transform([r,0])[0] - ax.transData.transform([0,0])[0]
#     marker_size = (1.*r_)**2
#     s_list.append(marker_size)
#     if j < num_atoms:
#         colors.append(cmap((r-min_sigma)/(max_sigma-min_sigma)))
#     else:
#         colors.append((0,0,0,0))
# 
# #%%
# 
# print('Rendering...')
# 
# import matplotlib as mpl
# mpl.rcParams['animation.ffmpeg_path'] = '/opt/homebrew/bin/ffmpeg'
# 
# 
# fig,ax = plt.subplots(figsize=[7,7])     
# cmap = matplotlib.colormaps['winter']
# ax.set_xlim((-1.2*box/2,1.2*box/2))
# ax.set_ylim((-1.2*box/2,1.2*box/2))
# 
# ax.set_aspect('equal', adjustable='box')
# plt.axis('off')
# 
# #plt.tight_layout(pad=0.2)
# 
# coords = pbc_coords_list[0]
# 
# x = pbc_coords_list[0][:,0]
# y = pbc_coords_list[0][:,1]
# 
# scatter = ax.scatter(x, y, marker='o', facecolor=colors, edgecolor='k', linewidth=0.25, s=s_list) 
# line, = ax.plot(pbc_xplot_list[0],pbc_mask_list[0],color='grey',linestyle='-',zorder=200)
# 
# plt.plot([-Lx/2,Lx/2],[-Ly/2,-Ly/2],color='k',zorder=400)
# plt.plot([-Lx/2,Lx/2],[Ly/2,Ly/2],color='k',zorder=400)
# 
# plt.plot([-Lx/2,-Lx/2],[-Ly/2,Ly/2],color='k',zorder=400)
# plt.plot([Lx/2,Lx/2],[-Ly/2,Ly/2],color='k',zorder=400)
# 
# plt.tight_layout()
# 
# def update(frame):
#     scatter.set_offsets(list(zip(pbc_coords_list[frame][:,0], pbc_coords_list[frame][:,1])))  # Update scatter points
#     
#     line.set_data(pbc_xplot_list[frame], pbc_mask_list[frame])
#     
#     
#     if frame % (num_frames // 10) == 0 :
#         percent_done = (frame + 1) / num_frames * 100
#         print(f"Rendering... {percent_done:.1f}% completed")
#     
#     return scatter, line
# 
# num_frames = len(pbc_coords_list)
# 
# ani = FuncAnimation(fig, update, frames=num_frames, blit=True)
# 
# writer = FFMpegWriter(fps=15)
# ani.save('voro_at2_0_pon_'+pon+'_poff_'+poff+'_phi'+target_phi+'_T_'+target_temp+'_gamma_'+gamma+'_theta_cutoff_'+theta_cutoff+'_rr_'+rr+'.mp4', writer=writer, dpi=300)
# 
# plt.close(fig)
# =============================================================================
