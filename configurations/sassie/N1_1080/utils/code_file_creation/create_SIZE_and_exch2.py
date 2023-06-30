
import os
import numpy as np
import argparse
import sys
sys.path.insert(1,os.path.join('..','..','utils'))
import sassie_functions as sf

def write_size_file(size_path,nPx,sNx,sNy):
    lines = 'C\n'
    lines += 'CBOP\n'
    lines += 'C\n'
    lines += 'C     /==========================================================\n'
    lines += 'C     | SIZE.h Declare size of underlying computational grid.    |\n'
    lines += 'C     |==========================================================|\n'
    lines += 'C     | The design here support a three-dimensional model grid   |\n'
    lines += 'C     | with indices I,J and K. The three-dimensional domain     |\n'
    lines += 'C     | is comprised of nPx*nSx blocks of size sNx along one axis|\n'
    lines += 'C     | nPy*nSy blocks of size sNy along another axis and one    |\n'
    lines += 'C     | block of size Nz along the final axis.                   |\n'
    lines += 'C     | Blocks have overlap regions of size OLx and OLy along the|\n'
    lines += 'C     | dimensions that are subdivided.                          |\n'
    lines += 'C     \==========================================================/\n'
    lines += 'C     \ev\n'
    lines += 'CEOP\n'
    lines += 'C     Voodoo numbers controlling data layout.\n'
    lines += 'C     sNx - No. X points in sub-grid.\n'
    lines += 'C     sNy - No. Y points in sub-grid.\n'
    lines += 'C     OLx - Overlap extent in X.\n'
    lines += 'C     OLy - Overlat extent in Y.\n'
    lines += 'C     nSx - No. sub-grids in X.\n'
    lines += 'C     nSy - No. sub-grids in Y.\n'
    lines += 'C     nPx - No. of processes to use in X.\n'
    lines += 'C     nPy - No. of processes to use in Y.\n'
    lines += 'C     Nx  - No. points in X for the total domain.\n'
    lines += 'C     Ny  - No. points in Y for the total domain.\n'
    lines += 'C     Nr  - No. points in Z for full process domain.\n'
    lines += '       INTEGER sNx\n'
    lines += '       INTEGER sNy\n'
    lines += '       INTEGER OLx\n'
    lines += '       INTEGER OLy\n'
    lines += '       INTEGER nSx\n'
    lines += '       INTEGER nSy\n'
    lines += '       INTEGER nPx\n'
    lines += '       INTEGER nPy\n'
    lines += '       INTEGER Nx\n'
    lines += '       INTEGER Ny\n'
    lines += '       INTEGER Nr\n'
    lines += '       PARAMETER (\n'
    lines += '     &           sNx =  '+str(sNx)+',\n'
    lines += '     &           sNy =  '+str(sNy)+',\n'
    lines += '     &           OLx =   8,\n'
    lines += '     &           OLy =   8,\n'
    lines += '     &           nSx =   1,\n'
    lines += '     &           nSy =   1,\n'
    lines += '     &           nPx =   '+str(nPx)+',\n'
    lines += '     &           nPy =   1,\n'
    lines += '     &           Nx  = sNx*nSx*nPx,\n'
    lines += '     &           Ny  = sNy*nSy*nPy,\n'
    lines += '     &           Nr  =  90 )\n'
    lines += ' \n'
    lines += 'C     MAX_OLX  - Set to the maximum overlap region size of any array\n'
    lines += 'C     MAX_OLY    that will be exchanged. Controls the sizing of exch\n'
    lines += 'C                routine buufers.\n'
    lines += '       INTEGER MAX_OLX\n'
    lines += '       INTEGER MAX_OLY\n'
    lines += '       PARAMETER ( MAX_OLX = OLx,\n'
    lines += '     &            MAX_OLY = OLy )\n'
    lines += ' \n'
    lines += '       integer     nobcs\n'
    lines += '       parameter ( nobcs = 4 )\n'

    f = open(size_path,'w')
    f.write(lines)
    f.close()

def write_exch2_file(exch2_path,llc,n,blank_list):
    lines = '# EXCH2 Package: Wrapper-2 User Choice\n'
    lines += '#--------------------\n'
    lines += '#  preDefTopol   :: pre-defined Topology selector:\n'
    lines += '#                :: = 0 : topology defined from processing "data.exch2";\n'
    lines += '#                :: = 1 : simple, single facet topology;\n'
    lines += '#                :: = 2 : customized topology (w2_set_myown_facets)\n'
    lines += '#                :: = 3 : 6-facet Cube (3 face-dims: nRed, nGreen, nBlue).\n'
    lines += '#  dimsFacets    :: facet pair of dimensions (n1x,n1y, n2x,n2y ...)\n'
    lines += '#  facetEdgeLink :: Face-Edge connectivity map:\n'
    lines += '#    facetEdgeLink(i,j)=XX.1 : face(j)-edge(i) (i=1,2,3,4 <==> N,S,E,W)\n'
    lines += '#    is connected to Northern edge of face "XX" ; similarly,\n'
    lines += '#    = XX.2 : to Southern.E, XX.3 = Eastern.E, XX.4 = Western.E of face "XX"\n'
    lines += '#  blankList     :: List of "blank" tiles\n'
    lines += '#  W2_mapIO      :: global map IO selector (-1 = old type ; 0 = 1 long line in X\n'
    lines += '#                :: 1 = compact, mostly in Y dir)\n'
    lines += '#  W2_printMsg   :: option for information messages printing\n'
    lines += '#                :: <0 : write to log file ; =0 : minimum print ;\n'
    lines += '#                :: =1 : no duplicated print ; =2 : all processes do print\n'
    lines += '#--------------------\n'
    lines += ' &W2_EXCH2_PARM01\n'
    lines += '# W2_printMsg= 1,\n'
    lines += '  W2_mapIO   = 1,\n'
    lines += '#\n'
    lines += '#-- 5 facets N0_1080 topology:\n'
    lines += '  preDefTopol = 0,\n'
    lines += '  dimsFacets = '+str(llc)+', '+str(n)+', '+str(llc)+', '+str(n)+', '+str(llc)+', '+str(llc)+', '+str(n)+', '+str(llc)+', '+str(n)+', '+str(llc)+', 0, 0,\n'
    lines += '#\n'
    lines += '# the north side of face 1 is connected to the west side of face 3\n'
    lines += '# the south side of face 1 is connected to nothing\n'
    lines += '# the east side of face 1 is connected to the west side of face 2\n'
    lines += '# the west side of face 1 is conneced to the north side of face 5\n'
    lines += '  facetEdgeLink(1,1)= 3.4, 0., 2.4, 5.1,\n'
    lines += '#\n'
    lines += '# the north side of face 2 is connected to the south side of face 3\n'
    lines += '# the south side of face 2 is connected to nothing\n'
    lines += '# the east side of face 2 is connected to the south side of face 4\n'
    lines += '# the west side of face 2 is conneced to the east side of face 1\n'
    lines += '  facetEdgeLink(1,2)= 3.2, 0. , 4.2 , 1.3,\n'
    lines += '#\n'
    lines += '# the north side of face 3 is connected to the west side of face 5\n'
    lines += '# the south side of face 3 is connected to the north side of face 2\n'
    lines += '# the east side of face 3 is connected to the west side of face 4\n'
    lines += '# the west side of face 3 is conneced to the north side of face 1\n'
    lines += '  facetEdgeLink(1,3)= 5.4, 2.1, 4.4, 1.1,\n'
    lines += '#\n'
    lines += '# the north side of face 4 is connected to the south side of face 5\n'
    lines += '# the south side of face 4 is connected to the east side of face 2\n'
    lines += '# the east side of face 4 is connected to nothing\n'
    lines += '# the west side of face 4 is conneced to the east side of face 3\n'
    lines += '  facetEdgeLink(1,4)= 5.2 , 2.3 , 0. , 3.3,\n'
    lines += '#\n'
    lines += '# the north side of face 5 is connected to the west side of face 1\n'
    lines += '# the south side of face 5 is connected to the north side of face 4\n'
    lines += '# the east side of face 5 is connected to nothing\n'
    lines += '# the west side of face 5 is conneced to the north side of face 3\n'
    lines += '  facetEdgeLink(1,5)= 1.4, 4.1 , 0. , 3.1,\n'
    lines += '#\n'
    lines += '  blankList = '+str(blank_list[0])+',\n'
    for tile in blank_list[1:]:
        lines += '   '+str(tile)+',\n'
    lines += ' &\n '
    lines += '\n '

    f = open(exch2_path, 'w')
    f.write(lines)
    f.close()
    

def create_blank_list_files(bathy_file,sNx,sNy,testing=False):

    llc = 1080
    n = 720

    n_rows = 4*n+llc
    n_cols = llc

    # step 1: get the bathy grid
    bathy_file = os.path.join('..', 'input', bathy_file)
    bathy_grid = np.fromfile(bathy_file, '>f4')
    bathy_grid_compact = np.reshape(bathy_grid, (n_rows, n_cols))
    bathy_grid_faces = sf.sassie_n1_compact_to_faces(bathy_grid_compact, llc, n)

    # step 2: run through all of the faces and count up the grid cells
    # if all the points are above sea level, then count em as blank
    tile_counter = 0
    blank_list = []

    if sNx == 180 and sNy == 180:
        manual_blank_list = [11, 17, 18, 29, 43, 74, 79, 105, 120]
    else:
        manual_blank_list = []

    if testing:
        print('Testing is true')

        n_test_cells = 2
        test_cell_counter = 0

        for face in range(1, 6):

            print('Face ' + str(face), np.shape(bathy_grid_faces[face]))

            nPx = int(np.shape(bathy_grid_faces[face])[1] / sNx)
            nPy = int(np.shape(bathy_grid_faces[face])[0] / sNy)

            for j in range(nPy):
                for i in range(nPx):
                    tile_counter += 1

                    bathy_subset = bathy_grid_faces[face][sNy * j:sNy * (j + 1), sNx * i:sNx * (i + 1)]

                    if test_cell_counter<n_test_cells and np.any(bathy_subset<0):
                        test_cell_counter+=1
                        print('   The blank tile is number '+str(tile_counter))
                    else:
                        blank_list.append(tile_counter)
    else:

        for face in range(1,6):

            print('Face '+str(face),np.shape(bathy_grid_faces[face]))

            nPx = int(np.shape(bathy_grid_faces[face])[1] / sNx)
            nPy = int(np.shape(bathy_grid_faces[face])[0] / sNy)

            for j in range(nPy):
                for i in range(nPx):
                    tile_counter += 1

                    bathy_subset = bathy_grid_faces[face][sNy * j:sNy * (j + 1), sNx * i:sNx * (i + 1)]

                    if np.all(bathy_subset>=0) or tile_counter in manual_blank_list:
                        blank_list.append(tile_counter)

    print('Counted '+str(tile_counter)+' total tiles')
    print('    of which '+str(len(blank_list))+' are blank')
    print('    and ' + str(tile_counter-len(blank_list)) + ' are not')

    size_path = os.path.join('..','code','SIZE.h')
    write_size_file(size_path, tile_counter - len(blank_list),sNx,sNy)

    exch2_path = os.path.join('..', 'input', 'data.exch2')
    write_exch2_file(exch2_path, llc, n, blank_list)




if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--bathy_file", action="store",
                        help="The bathymetry file to read in.", dest="bathy_file",
                        type=str, required=False, default='bathymetry.bin')

    parser.add_argument("-sNx", "--sNx", action="store",
                        help="The tile size in the x direction", dest="sNx",
                        type=int, required=True, default=30)

    parser.add_argument("-sNy", "--sNy", action="store",
                        help="The tile size in the y direction", dest="sNy",
                        type=int, required=True, default=30)

    parser.add_argument("-t", "--test", action="store",
                        help="Whether or not to configure for a local compile test (nPx = 1).", dest="testing",
                        type=int, required=False, default=0)

    args = parser.parse_args()
    bathy_file = args.bathy_file
    sNx = args.sNx
    sNy = args.sNy
    testing = args.testing
    if testing==0:
        testing = False
    else:
        testing = True

    create_blank_list_files(bathy_file, sNx, sNy, testing)