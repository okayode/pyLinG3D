import numpy as np
import matplotlib.pyplot as plt
import os
from os import path

from .utils import patch_b
from .utils import patch_c
from .utils import patch_y
from .utils import patch_r
from .utils import DrawBackground
from .utils import DefineColorPalette
from .utils import find_a

def linG3DClone(pathData,cloneNum,IsGradient,xmin,xmax,ymin,ymax,tmin,tmax,fileStep,toPrint):
    """
    This is a companion code for the paper "LinG3D: Visualizing the
    Spatio-Temporal Dynamics of Clonal Evolution" by A. Hu, A.M.E.
    Ojwang', K.D. Olumoyin, and K.A. Rejniak
    
    This code generates the 3D lineage tree of one clone of number
    specified in 'cloneNum' taking into account only the cells that
    survived to the end of simulation. It uses data from directory
    'pathData'.
    
    The following parameters need to be specified:
        pathData  -- directory with input data
        cloneNum  -- clone number to be drawn
        IsGradient -- 1 to draw drug in the background, 0 not to draw
        xmin,xmax,ymin,ymax -- dimensions of the spacial domain
        tmin, tmax          -- dimensions of the temporal domain
        fileStep            -- frequency of the sampled data
        toPrint    -- 1 to save the generated figure, 0 not to save
        
    It requires the following data in the pathData/data/ directory:
        cell_history.txt -- file with info about each cell
        cellID_##.txt    -- cell IDs in a file with index number ##
        cellXY_##.txt    -- cell coordinates in a file with index ##
        drug.txt         -- concentration of a drug for background
    
    for the examples discussed in the paper use:
        example 1: pathData='exampleB05';  cloneNum between 0 and 9
        example 2: pathData='exampleB005'; cloneNum between 0 and 147
        example 3: pathData='exampleExp';  numClones=10;
        
    January 10, 2024
    """
    
    dataDirectory = '/data/' # directory with cell and drug data
    timeStep=(tmax-tmin)/(2.5*(xmax-xmin))

    plt.grid()
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([tmin, tmax // timeStep])
    ax.set_zlim([ymin, ymax])
    
    col=DefineColorPalette()
    Ncol=col[:,0].shape[0]

    if toPrint==1:
        pathFigs=pathData+'/fig_clones';
        if path.exists(pathFigs):
            pass
        else:
            os.mkdir(pathFigs)
    
    # draw background with drug gradient
    if IsGradient == 1:
        drug = np.loadtxt(pathData + dataDirectory + 'drug.txt')
        DrawBackground(drug, tmax, timeStep, xmin, xmax, ymin, ymax, ax)
        
    # load cell history file
    hist = np.loadtxt(pathData + dataDirectory + 'cell_history.txt')
    # [cell ID, clone ID, mother ID, birth iter, div / death iter]
    
    print('clone=' + str(cloneNum))
    
    # load indices of all cells from cloneNum
    indLast, = find_a(hist[:,1]==cloneNum)
    tmp0 = find_a(hist[:,2]==tmax)
    
    np.concatenate((indLast,tmp0),axis=None)
    
    indLast = np.array(indLast)
    # print(indLast)
                   
    # define matrix of line segments (3D branches) to draw
    matrix_to_draw = np.zeros((1, 6))  # [x1,k1,y1,x2,k2,y2]
    Nmatrix = 0
                   
    for ii in range(len(indLast)):  # for every cell with index in indLast
        if ii % 100 == 0:
            print('... calculating')
                
        cellNum = int(hist[indLast[ii],0])  # cell ID
        # print(cellNum)
        mothNum = int(hist[indLast[ii],2])  # mother ID
        strtNum = int(hist[indLast[ii],3])  # cell birth
        endNum = int(hist[indLast[ii],4])
        endNum = max(tmax, min(endNum, tmax)) # cell div / death / tmax
        
        # find all appearances of the cellNum
        kkStart = fileStep * np.floor(strtNum / fileStep)  # initial file number
        kkEnd = fileStep * np.floor(endNum / fileStep)  # final file number
        kkStart = int(kkStart)
        kkEnd = int(kkEnd)
        # print(kkStart)
        
        for kk in range(kkEnd, kkStart, - fileStep):  # inspect all files
            # cell ID and cell XY from the first file
            fileMeID = np.loadtxt(pathData + dataDirectory + 'cellID_' + str(kk) + '.txt')  
            fileMeXY = np.loadtxt(pathData + dataDirectory + 'cellXY_' + str(kk) + '.txt')
                
            indMe, = find_a(fileMeID == cellNum) # find current indices of cellID
            # print(indMe)
                        
            fileMe2ID = np.loadtxt(pathData + dataDirectory + 'cellID_' + str(kk - fileStep) + '.txt')
            fileMe2XY = np.loadtxt(pathData + dataDirectory + 'cellXY_' + str(kk - fileStep) + '.txt')
                
            indMe2, = find_a(fileMe2ID == cellNum)  # find current indices of cellID
            # print(indMe2)
            
            if indMe.size == 0:
                pass
            elif indMe2.size == 0:
                #print(indMe)
                while kkStart < int(hist[mothNum-1,3]):  # find file with the grand-mother cell
                    mothNum = int(hist[mothNum-1,2])
                fileMe2ID = np.loadtxt(pathData + dataDirectory + 'cellID_' + str(kkStart) + '.txt')
                fileMe2XY = np.loadtxt(pathData + dataDirectory + 'cellXY_' + str(kkStart) + '.txt')
                # print(fileMe2XY)
                    
                indMe2, = find_a(fileMe2ID == mothNum)  # find current indices of mother cellID
                # print(indMe2)
                if indMe2.size == 0:
                    pass
                else:                     
                    matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                    tmp = np.array([0.0,0.0])
                    tmp2 = []
                    if fileMe2XY[indMe2][0].any() == 0:
                        tmp2.append(tmp)
                    else:
                        tmp2.append(fileMe2XY[indMe2][0])
                        
                    # print(tmp2)
                    # print(tmp2[0][1])
                    matrix_to_draw[Nmatrix,:] = [fileMeXY[indMe][0][0], kkStart + fileStep,\
                                                fileMeXY[indMe][0][1], tmp2[0][0],\
                                                kkStart, tmp2[0][1]]
                    Nmatrix += 1  # save branch to draw [x1,t1,y1,x2,t2,y2] 
                                
            #
            else:
                matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                tmp = np.array([0.0,0.0])
                tmp2 = []
                if fileMe2XY[indMe2][0].any() == 0:
                    tmp2.append(tmp)
                else:
                    tmp2.append(fileMe2XY[indMe2][0])
                # print(indMe2)
                # print(fileMe2XY)
                        
                matrix_to_draw[Nmatrix,:] = [tmp2[0][0], kk - fileStep,\
                                            tmp2[0][1], fileMeXY[indMe][0][0],\
                                            kk, fileMeXY[indMe][0][1]]
                Nmatrix += 1  # save branch to draw [x1,t1,y1,x2,t2,y2]
                   
    
    # drawing clones
    NumCol = cloneNum % Ncol  # clone color
    # plt.title('cells from clone=' + str(cloneNum))
    plt.ylabel('iterations/time x ' + str(timeStep),fontsize=8)
        
    for ii in range(Nmatrix):
        x = [matrix_to_draw[ii][0], matrix_to_draw[ii][3]]
        y = [matrix_to_draw[ii][1], matrix_to_draw[ii][4]]
        y = [yi / timeStep for yi in y]
        z = [matrix_to_draw[ii][2], matrix_to_draw[ii][5]]
        ax.plot(x, y, z, c=col[NumCol,:3],linewidth=1.1)
    tm = tmax/timeStep
    ax.view_init(40,-130)
    ax.set_box_aspect([1,2.5,1])
    # ax.set_xticks([-xmin,-xmin*0.5,0,0.5*xmax,xmax])
    ax.set_yticks([tm,tm-tm*.2,tm-tm*.4,tm-tm*.6,tm-tm*.8,0])
    # ax.set_zticks([-xmin,-xmin*0.5,0,0.5*xmax,xmax])
    ax.set_ylim(0,tm)
            
    if toPrint==1:
        plt.savefig(pathFigs + "/tree_clone_"+ str(cloneNum)+'.jpg', dpi=300)
    
    plt.show()
    
fig = plt.figure()
ax = plt.axes(projection='3d')






