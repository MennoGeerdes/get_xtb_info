"""
Extracts key info from xTB thermodynamic (frequency) output files (and geometry optimization if xyz file is not specified).
 ~ Might be incorrect for coordination number, leading to an incorrect sum of atom nrs of coordinated species as well, 
 ~ check manually untill tested to that regard.
"""
#%% Imports

import math, argparse, os, fnmatch

#%% Functions

def getxtbinfo(outpath,X=None,H=None,xyzpath=None):
    info={}
    info['HOMO']='-'
    info['LUMO']='-'
    info['Escf']='-'
    info['Egibbs']='-'
    info['Xspecies'] ='-'
    info['Xnr'] ='-'
    info['Xcharge']='-'
    info['Xcoordination']='-'
    info['Xsumcoord']='-'
    info['Hcharge']='-'
    info['XH']='-'
   
    # read the xtb.out file
    with open(outpath, 'r', encoding=("utf-8")) as xtbout:
        text = [_.strip() for _ in xtbout.readlines()]
        
        # first find coordinates in output file
        start = 0
        for i, line in enumerate(text):
            if line == "final structure:":
                start = i+2
                break
        
        # find the final property printout
        for i, line in enumerate(text):
            if "Property Printout" in line:
                final = i
        
        # loop through the rest of the lines for thermodynamic results
        for i, line in enumerate(text):
            if i < final:
                continue
            if '(HOMO)' in line:
                HOMO = line.split()[-2]
            elif '(LUMO)' in line:
                LUMO = line.split()[-2]
            elif 'TOTAL ENERGY' in line:
                ESCF = line.split()[3]
            elif 'TOTAL FREE ENERGY' in line:
                EGIBBS = line.split()[4]
            elif 'total WBO' in line:
                wbo = i
            elif 'covCN' in line:
                covCN = i

        info['HOMO']=HOMO
        info['LUMO']=LUMO
        info['Escf']=ESCF
        info['Egibbs']=EGIBBS
        
        # if nr of X in input geometry is given, extract info regarding atom X
        if X:
            info['Xspecies'] = text[covCN+X].split()[2]
            info['Xnr'] = text[covCN+X].split()[1]
            info['Xcharge']=(text[covCN+X].split()[4])
            
            # the coordination number and coordinating species are extracted in the next steps,
            # this may yet have to be tested for correctness!
            wbolist = (text[wbo+X].split())[5::3]
            coordination=0
            for i in wbolist:
                if float(i) > 0.4:
                    coordination+=1
            info['Xcoordination']=coordination
            # atom nrs of species bonded to X, is correct if coordination number is correct.
            coordnrs=[]
            for i in range(coordination):
                coordnrs.append(text[wbo+X].split()[4+3*i])
            
            # sum of atomnrs of the coordinating species to X
            info['Xsumcoord']=str(sum([int(text[covCN+int(i)].split()[1]) for i in coordnrs]))         
            
        # if nr of H in input geometry is given, extract info regarding H
        if H:
            info['Hcharge']=(text[covCN+H].split()[4])

        # some selected distances are given directly in the output files, but not all of them. This is a more general method,
        # since the full geometry is given we can just use the Pythagorean theorem.
        # if X and H numbers are specified, but no path to xtbopt.xyz file is given: extract X-H distance from xtb.out file
        if X and H and xyzpath==None:
            x = [float(text[start+X].split()[0]),float(text[start+X].split()[1]),float(text[start+X].split()[2])]
            h = [float(text[start+H].split()[0]),float(text[start+H].split()[1]),float(text[start+H].split()[2])]
            
            # in the coordinates specified in xtb.out files, the distance between atoms is multiplied by the factor below, idk why.
            info['XH'] = (math.sqrt( (x[0]-h[0])**2 + (x[1]-h[1])**2 + (x[2]-h[2])**2 ))/1.88972594929721
        
        # The distance between atoms X and H can also be extracted from the xyz file that is written by the xTB program (xtbopt.xyz)        
        # if X, H and a path to xtbopt.xyz are given, extract X-H bond length from xtbopt.xyz file
        # I think this way might be better to extract distances, because of the weird factor in the xtb.out file.
        # differences are at most 1e-7
        if X and H and xyzpath:
            with open(xyzpath, 'r', encoding=("utf-8")) as xtbxyz:
                text = [_.strip() for _ in xtbxyz.readlines()]
                Xx=float(text[1+X].split()[1])
                Xy=float(text[1+X].split()[2])
                Xz=float(text[1+X].split()[3])
                Hx=float(text[1+H].split()[1])
                Hy=float(text[1+H].split()[2])
                Hz=float(text[1+H].split()[3])
                XH=math.sqrt((Xx-Hx)**2+(Xy-Hy)**2+(Xz-Hz)**2)
                info['XH']=str(XH)
    # return the dictionary with added info
    return(info)

# this callable returns the info from the getxtbinfo function as a csv line, which is printed out as output of this script
def linewriter(info):
    line = ''
    for key in info:
        line += ',' + str(info[key])
    return line

#%% Main code
'''
Extracts key info from xtb output files
Syntaxis:
    > python get_xtb_info.py xtb.out xtbopt.xyz=None X=None H=None
    - for help, enter: $ python get_xtb_info.py -h
This script prints a csv line with info. To see in what order info is printed add --header to print the header as well.
The start of the line is a comma, so you can loop over different paths to xtb.out,
and add the interation number before the line as index number.
'''

if __name__ == '__main__':

    # get arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('outpath', help='path to xTB output file')
    parser.add_argument('-xyz', help='path to xyz file from xTB geometry optimization')
    parser.add_argument('-X', help='number of atom X in xTB input', type=int)
    parser.add_argument('-H', help='number of atom H in xTB input', type=int)
    parser.add_argument('--header', help='if true: only print header. A valid path to xTB output is still required.', action='store_true')
    args = parser.parse_args()
    
    out = args.outpath
    x = args.X
    h = args.H
    xyz = args.xyz
        
    # retrieve output paths
    cwd = os.getcwd()
    if not os.path.isabs(out):
        out = os.path.join(cwd,out)
        folder, name = os.path.split(out)
        for f in os.listdir(folder):
            if not fnmatch.fnmatch(f,name):
                continue
            else:
                out = os.path.join(folder,f)
    
    if xyz:
        if not os.path.isabs(xyz):
            xyz = os.path.join(cwd,xyz)
            folder, name = os.path.split(xyz)
            for f in os.listdir(folder):
                if not fnmatch.fnmatch(f,name):
                    continue
                else:
                    xyz = os.path.join(folder,f)
    
    # get info and prepare line
    try:
        info = getxtbinfo(out,X=x,H=h,xyzpath=xyz)
    except (SystemExit, KeyboardInterrupt):
        raise
    except:
        info = None
    
    # print header if option --header is added, else print line with info, or error message if info == None
    if info:
        line = linewriter(info)
        header = 'index'
        for key in info:
            header += (',' + key)
        if args.header:
            print(header)
        else:
            print(line)
    else:
        print('Something went wrong. Check path to xTB output file?')