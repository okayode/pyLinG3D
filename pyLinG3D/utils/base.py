import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# define a MATLAB-like patch function

def patch_b(ax, x, y, z, color='blue',alpha=0.25):
    pc = Poly3DCollection([list(zip(x,y,z))],linewidths=1,alpha=alpha)       
    pc.set_facecolor(color)                                                       
    ax.add_collection3d(pc)                         
    return pc

def patch_c(ax, x, y, z, color='cyan',alpha=0.25):
    pc = Poly3DCollection([list(zip(x,y,z))],linewidths=1,alpha=alpha)       
    pc.set_facecolor(color)                                                       
    ax.add_collection3d(pc)                         
    return pc

def patch_y(ax, x, y, z, color='yellow',alpha=0.25):
    pc = Poly3DCollection([list(zip(x,y,z))],linewidths=1,alpha=alpha)       
    pc.set_facecolor(color)                                                       
    ax.add_collection3d(pc)                         
    return pc

def patch_r(ax, x, y, z, color='red',alpha=0.25):
    pc = Poly3DCollection([list(zip(x,y,z))],linewidths=1,alpha=alpha)       
    pc.set_facecolor(color)                                                       
    ax.add_collection3d(pc)                         
    return pc

def DrawBackground(drug,tmax,timeStep,xmin,xmax,ymin,ymax,ax):
    drugmin=np.min(np.min(drug)); drugmax=np.max(np.max(drug)); drugstep=(drugmax-drugmin)/4
    kk=tmax/timeStep
    (Nx,Ny)=drug.shape; hgx=(xmax-xmin)/Nx; hgy=(ymax-ymin)/Ny;
    
    for ii in range(1,Nx+1):
        for jj in range(1,Ny+1):
            
            if drug[ii-1,jj-1]>=drugmin and drug[ii-1,jj-1]<(drugmin+drugstep):
                aa=[xmin+(ii-1)*hgx, xmin+ii*hgx, xmin+ii*hgx, xmin+(ii-1)*hgx, xmin+(ii-1)*hgx]
                bb=[kk, kk, kk, kk, kk]
                cc=[ymin+(jj-1)*hgy, ymin+(jj-1)*hgy, ymin+jj*hgy, ymin+jj*hgy, ymin+(jj-1)*hgy]
                patch_b(ax, aa, bb, cc)
                
            elif drug[ii-1,jj-1]>=(drugmin+drugstep) and drug[ii-1,jj-1]<(drugmin+2*drugstep):
                aa=[xmin+(ii-1)*hgx, xmin+ii*hgx, xmin+ii*hgx, xmin+(ii-1)*hgx, xmin+(ii-1)*hgx]
                bb=[kk, kk, kk, kk, kk]
                cc=[ymin+(jj-1)*hgy, ymin+(jj-1)*hgy, ymin+jj*hgy, ymin+jj*hgy, ymin+(jj-1)*hgy]
                patch_c(ax, aa, bb, cc)
                
            elif drug[ii-1,jj-1]>=(drugmin+2*drugstep) and drug[ii-1,jj-1]<(drugmin+3*drugstep):
                aa=[xmin+(ii-1)*hgx, xmin+ii*hgx, xmin+ii*hgx, xmin+(ii-1)*hgx, xmin+(ii-1)*hgx]
                bb=[kk, kk, kk, kk, kk]
                cc=[ymin+(jj-1)*hgy, ymin+(jj-1)*hgy, ymin+jj*hgy, ymin+jj*hgy, ymin+(jj-1)*hgy]
                patch_y(ax, aa, bb, cc)
                
            else:
                aa=[xmin+(ii-1)*hgx, xmin+ii*hgx, xmin+ii*hgx, xmin+(ii-1)*hgx, xmin+(ii-1)*hgx]
                bb=[kk, kk, kk, kk, kk]
                cc=[ymin+(jj-1)*hgy, ymin+(jj-1)*hgy, ymin+jj*hgy, ymin+jj*hgy, ymin+(jj-1)*hgy]
                patch_r(ax, aa, bb, cc)
                
def DefineColorPalette():
    col = np.array([[255,0,255],[255,0,0],[0,255,255],[0,0,255],[0,255,0],[0,0,0],[255,191,0],\
                    [255,255,0],[191,255,0],[128,128,0],[255,182,193],[0,191,255],[0,128,255],\
                    [250,235,215],[128,0,255],[154,205,50],[255,0,128],[102,0,0],[102,77,0],\
                    [0,102,102],[204,204,255],[255,204,255],[153,204,255],[255,153,153],[0,153,0],\
                    [0,153,153],[153,0,77],[255,228,225],[128,0,0],[102,102,153],[153,255,204],\
                    [218,112,214],[255,128,0],[192,192,192],[128,128,128],[75,0,130],[165,42,42],\
                    [216,191,216],[220,20,60],[245,222,179],[255,99,71],[255,127,80],[205,92,92],\
                    [240,128,128],[233,150,122],[250,128,114],[255,160,122],[255,69,0],\
                    [255,140,0],[255,165,0],[255,215,0],[184,134,11],[218,165,32],[0,100,0],\
                    [255,240,245],[188,143,143],[255,248,220],[50,205,50],[144,238,144],\
                    [152,251,152],[143,188,143],[0,250,154],[0,255,127],[46,139,87],[102,205,170],\
                    [60,179,113],[32,178,170],[47,79,79],[0,128,128],[0,139,139],[240,230,140],\
                    [245,245,220],[224,255,255],[0,206,209],[255,228,181],[255,20,147],\
                    [175,238,238],[127,255,212],[176,224,230],[95,158,160],[70,130,180],\
                    [100,149,237],[222,184,135],[30,144,255],[238,232,170],[189,183,107],\
                    [107,142,35],[124,252,0],[127,255,0],[173,255,47],\
                    [178,34,34],[221,160,221],[255,235,205]])/255
    return col

def find_a(condition):
    # res, = np.nonzero(np.ravel(condition))
    res = np.nonzero(np.ravel(condition))
    return res

def find_b(condition):
    res, = np.nonzero(np.ravel(condition))
    # res = np.nonzero(np.ravel(condition))
    return res