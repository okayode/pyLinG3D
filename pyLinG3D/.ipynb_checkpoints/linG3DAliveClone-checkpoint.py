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
from .utils import find_b

class LinG3DAliveClone:
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
    def __init__(self, pathData, cloneNnum, IsGradient, xmin, xmax, ymin, ymax, tmin, tmax, fileStep, toPrint):
        self.pathData = pathData
        self.cloneNum = cloneNnum
        self.IsGradient = IsGradient
        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.tmin, self.tmax = tmin, tmax
        self.fileStep = fileStep
        self.toPrint = toPrint

        self.dataDirectory = '/data/'  # directory with cell and drug data
        self.timeStep = (tmax - tmin) / (2.5 * (xmax - xmin))

        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111, projection='3d')
        
        if self.toPrint == 1:
            self.pathFigs = os.path.join(self.pathData, 'fig_clonesAlive')
            if os.path.exists(self.pathFigs):
                pass
            else:
                os.mkdir(self.pathFigs)
    
    def linG3DAliveClone(self):
        self.ax.grid()
        self.ax.set_xlim([self.xmin, self.xmax])
        self.ax.set_ylim([self.tmin, self.tmax // self.timeStep])
        self.ax.set_zlim([self.ymin, self.ymax])

        col = DefineColorPalette()
        n_col = col[:, 0].shape[0]

        print('clone=' + str(self.cloneNum))

        # Draw background with drug gradient
        if self.IsGradient == 1:
            drug = np.loadtxt(self.pathData + self.dataDirectory + 'drug.txt')
            DrawBackground(drug, self.tmax, self.timeStep, self.xmin, self.xmax, self.ymin, self.ymax, self.ax)

        # Load cell history file
        hist = np.loadtxt(self.pathData + self.dataDirectory + 'cell_history.txt')

        # Load indices of all survived cells
        cell_id = np.loadtxt(self.pathData + self.dataDirectory + 'cellID_' + str(self.tmax) + '.txt')

        # Identify all survived cells from a given clone
        ind_last = []
        n_last = 0

        for ii in range(len(cell_id)):
            if hist[int(cell_id[ii]) - 1, 1] == self.cloneNum:
                ind_last.append(cell_id[ii])
                n_last += 1
        num_last = n_last
        for ii in range(num_last):
            mo_num = hist[int(ind_last[ii]) - 1, 2]
            while mo_num > 0 and (hist[int(mo_num) - 1, 1] == self.cloneNum):
                ind_last.append(mo_num)
                n_last += 1
                mo_num = hist[int(mo_num) - 1, 2]
        ind_last = np.unique(ind_last)

        # Define matrix of line segments (3D branches) to draw
        matrix_to_draw = np.zeros((1, 6))  # [x1, t1, y1, x2, t2, y2] X - time - Y axes
        n_matrix = 0

        for ii in range(len(ind_last)):
            if ii % 100 == 0:
                print('... calculating')

            cell_num = hist[int(ind_last[ii]) - 1, 0]
            mo_num = hist[int(ind_last[ii]) - 1, 2]
            strt_num = hist[int(ind_last[ii]) - 1, 3]
            end_num = hist[int(ind_last[ii]) - 1, 4]
            end_num = max(self.tmax, min(end_num, self.tmax))

            kk_start = self.fileStep * np.floor(strt_num / self.fileStep)
            kk_end = self.fileStep * np.floor(end_num / self.fileStep)
            kk_start = int(kk_start)
            kk_end = int(kk_end)

            for kk in range(kk_end, kk_start, -self.fileStep):
                file_me_id = np.loadtxt(self.pathData + self.dataDirectory + 'cellID_' + str(kk) + '.txt')
                file_me_xy = np.loadtxt(self.pathData + self.dataDirectory + 'cellXY_' + str(kk) + '.txt')

                ind_me = find_b(file_me_id == cell_num)

                file_me2_id = np.loadtxt(self.pathData + self.dataDirectory + 'cellID_' + str(kk - self.fileStep) + '.txt')
                file_me2_xy = np.loadtxt(self.pathData + self.dataDirectory + 'cellXY_' + str(kk - self.fileStep) + '.txt')

                ind_me2 = find_b(file_me2_id == cell_num)

                if ind_me.size == 0:
                    pass
                elif ind_me2.size == 0:
                    while kk_start < hist[int(mo_num) - 1, 3]:
                        mo_num = hist[int(mo_num) - 1, 2]
                    file_me2_id = np.loadtxt(self.pathData + self.dataDirectory + 'cellID_' + str(kk_start) + '.txt')
                    file_me2_xy = np.loadtxt(self.pathData + self.dataDirectory + 'cellXY_' + str(kk_start) + '.txt')

                    ind_me2 = find_b(file_me2_id == mo_num)
                    if ind_me2.size == 0:
                        pass
                    else:
                        matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                        matrix_to_draw[n_matrix, :] = [file_me_xy[ind_me][0][0], kk_start + self.fileStep,
                                                        file_me_xy[ind_me][0][1], file_me2_xy[ind_me2][0][0],
                                                        kk_start, file_me2_xy[ind_me2][0][1]]
                        n_matrix += 1
                else:
                    matrix_to_draw = np.row_stack((matrix_to_draw, np.zeros((1, 6))))
                    matrix_to_draw[n_matrix, :] = [file_me2_xy[ind_me2][0][0], kk - self.fileStep,
                                                    file_me2_xy[ind_me2][0][1], file_me_xy[ind_me][0][0],
                                                    kk, file_me_xy[ind_me][0][1]]
                    n_matrix += 1

        # Drawing clones
        num_col = self.cloneNum % n_col
        plt.ylabel('iterations/time x ' + str(self.timeStep), fontsize=8)

        for ii in range(n_matrix):
            x = [matrix_to_draw[ii, 0], matrix_to_draw[ii, 3]]
            y = [matrix_to_draw[ii, 1], matrix_to_draw[ii, 4]]
            y = [yi / self.timeStep for yi in y]
            z = [matrix_to_draw[ii, 2], matrix_to_draw[ii, 5]]
            self.ax.plot(x, y, z, c=col[num_col, :3], linewidth=1.3)
        tm = self.tmax / self.timeStep
        self.ax.view_init(40, -130)
        self.ax.set_box_aspect([1, 2.5, 1])
        self.ax.set_yticks([tm, tm - tm * 0.2, tm - tm * 0.4, tm - tm * 0.6, tm - tm * 0.8, 0])
        self.ax.set_ylim(0, tm)

        if self.toPrint == 1:
            plt.savefig(os.path.join(self.pathFigs, 'tree_alive_clone_' + str(self.cloneNum) + '.jpg'), dpi=300)

        plt.show()
           
