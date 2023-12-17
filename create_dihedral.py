import os
import numpy as np
import ase
from ase.visualize import view
import math
import pandas as pd
from ase import Atoms

def read_structures(filename):
    atoms = ase.io.read(filename)
    view(atoms)
    
def get_geometry(filename):
    atoms = []
    geometryx = []
    geometryy = []
    geometryz = []
    with open (filename, "r") as f:
        lines = f.readlines()
        for i in lines[2:]:
            data = i.split()
            atoms.append(data[0])
            geox = float(data[1]) 
            geoy = float(data[2]) 
            geoz = float(data[3])
            geometryx.append(geox)
            geometryy.append(geoy)
            geometryz.append(geoz)

    dic = {"Atom": atoms, "X": geometryx, "Y": geometryy, "Z": geometryz}    
    df = pd.DataFrame(dic)
    return df

def get_dihedral(filename, index1, index2, index3, index4):
### Create a dataframe to calculate the diihedral angle ###
    df = get_geometry(filename)
    
### Create 2 planes for the dihedral angle calculation ###
    atom1 = df["Atom"][int(index1)]
    coox = df["X"]
    atom2 = df["Atom"][int(index2)]
    cooy = df["Y"]
    atom3 = df["Atom"][int(index3)]
    cooz = df["Z"]
    atom4 = df["Atom"][int(index4)]
    x1 = coox[int(index1)]
    x2 = coox[int(index2)]
    x3 = coox[int(index3)]
    x4 = coox[int(index4)]
    y1 = cooy[int(index1)]
    y2 = cooy[int(index2)]
    y3 = cooy[int(index3)]
    y4 = cooy[int(index4)]
    z1 = cooz[int(index1)]
    z2 = cooz[int(index2)]
    z3 = cooz[int(index3)]
    z4 = cooz[int(index4)]

    #bond1 = (x2-x1)*x + (y2-y1)*y + (z2-z1)*z
    vector1 = np.array([(x2-x1), (y2-y1), (z2-z1)])
    #bond2 = (x3-x2)*x + (y3-y2)*y + (z3-z2)*z
    vector2 = np.array([(x3-x2), (y3-y2), (z3-z2)])
    #bond3 = (x4-x3)*x + (y4-y3)*y + (z4-z3)*z
    vector3 = np.array([(x4-x3), (y4-y3), (z4-z3)])

    pl1 = np.cross(vector1, vector2)
    pl2 = np.cross(vector2, vector3)
    ### Dihedral angle equation ###
    ### A1; A2; B1; B2; C1; C2 are coefficients
    ### A1x + B1y + C1z + D = 0 plane equation
    A1 = pl1[0]
    B1 = pl1[1]
    C1 = pl1[2]
    A2 = pl2[0]
    B2 = pl2[1]
    C2 = pl2[2]
    ### Numerator = A1*A2 + B1*B2 + C1*C2
    ### Denominator = (A1**2 + B1**2 + C1**2)**0.5 + (A2**2 + B2**2 + C2**2)**0.5
    numerator = A1*A2 + B1*B2 + C1*C2
    denominator = ((A1**2 + B1**2 + C1**2)**0.5) * ((A2**2 + B2**2 + C2**2)**0.5)

    theta = math.acos(numerator/denominator)
    theta_degree = math.degrees(theta)
    #print(A1, B1, C1)
    
    return theta_degree

def rotate_to_0(filename, rotate, non_rotate, D1, D2, D3, D4):
    ### First create a zero dihedral angle structure ###
    ### Rotate is the part you want to rotate ###
    ### Non_rotate is the part you wanna fix their geometry ###
    df = get_geometry(filename) # Please use initial file right here #
    new_df1 = df.iloc[rotate]
    new_df1.to_csv('temp_rotate.xyz', header=None, index=None, sep=' ', mode='a')
    new_df2 = df.iloc[non_rotate]
    new_df2.to_csv('temp_non_rotate.xyz', header=None, index=None, sep=' ', mode='a')
    
    with open ("temp_rotate.xyz", "r") as r:
        existing_content = r.read()

    with open ("rotate.xyz", "w") as r:  
        r.write(f"{len(rotate)}")
        r.write("\n\n")
        r.write(existing_content)

    with open ("temp_non_rotate.xyz", "r") as s:
        existing_content = s.read()

    with open ("non_rotate.xyz", "w") as s:  
        s.write(f"{len(non_rotate)}")
        s.write("\n\n")
        s.write(existing_content)
        
    ### Create 0 degree dihedral angle xyz file ###
    rotate_atoms = ase.io.read("rotate.xyz")
    non_rotate_atoms = ase.io.read("non_rotate.xyz")
    current_diangle = get_dihedral(filename, D1, D2, D3, D4) # D2 should be the rotated center (fixed) #
    
    D2x = df.iloc[D2]["X"]
    D2y = df.iloc[D2]["Y"]
    D2z = df.iloc[D2]["Z"]
    
    D3x = df.iloc[D3]["X"]
    D3y = df.iloc[D3]["Y"]
    D3z = df.iloc[D3]["Z"]
    
    rotate_atoms.rotate(current_diangle, (D3x-D2x, D3y-D2y, D3z-D2z), center=(D2x, D2y, D2z))
    rotate_atoms += non_rotate_atoms
    
    ase.io.write("0degree.xyz", rotate_atoms)
    # Delete all the intermediate file
    file_list = ["non_rotate.xyz","rotate.xyz","temp_rotate.xyz","temp_non_rotate.xyz"]
    for f in file_list:
        os.remove(f)
        
def check(D1, D2, D3, D4, N_non):
    ### Check the 0degree.xyz dihedral is correct or not ###
    ### N_non is the number of non-rotating atoms ###
    with open ("0degree.xyz", "r") as r:
        N_atoms = r.readline()
    rotate = [int(i) for i in range(int(N_atoms) - N_non)]
    non_rotate = [int(j) for j in range((int(N_atoms) - N_non), int(N_atoms))]
    #print(rotate)
    #print(non_rotate)
    not0 = True
    while not0:
        angle = get_dihedral("0degree.xyz", D1, D2, D3, D4)   
        if angle < 10e-6:
            not0 = False

        else:
            rotate_to_0("0degree.xyz", rotate, non_rotate, D1, D2, D3, D4)
            
    print("Successfully created 0 degree structure")
    
def rotate(i, D2, D3, N_non):
    ### Create structure with 5 degrees increasing ###
    ### i is the final angle +5 ex: 185 will create 5-180 file with gap 5 ###
    df = get_geometry("0degree.xyz")
    with open ("0degree.xyz", "r") as r:
        N_atoms = r.readline()
    rotate = [int(i) for i in range(int(N_atoms) - N_non)]
    non_rotate = [int(j) for j in range((int(N_atoms) - N_non), int(N_atoms))]
    for j in range(5, i, 5):
        new_df1 = df.iloc[rotate]
        new_df1.to_csv('temp_rotate.xyz', header=None, index=None, sep=' ', mode='a')
        new_df2 = df.iloc[non_rotate]
        new_df2.to_csv('temp_non_rotate.xyz', header=None, index=None, sep=' ', mode='a')
    
        with open ("temp_rotate.xyz", "r") as r:
            existing_content = r.read()

        with open ("rotate.xyz", "w") as r:  
            r.write(f"{len(rotate)}")
            r.write("\n\n")
            r.write(existing_content)

        with open ("temp_non_rotate.xyz", "r") as s:
            existing_content = s.read()

        with open ("non_rotate.xyz", "w") as s:  
            s.write(f"{len(non_rotate)}")
            s.write("\n\n")
            s.write(existing_content)
        
        rotate_atoms = ase.io.read("rotate.xyz")
        non_rotate_atoms = ase.io.read("non_rotate.xyz")
        
        D2x = df.iloc[D2]["X"]
        D2y = df.iloc[D2]["Y"]
        D2z = df.iloc[D2]["Z"]
    
        D3x = df.iloc[D3]["X"]
        D3y = df.iloc[D3]["Y"]
        D3z = df.iloc[D3]["Z"]
    
        rotate_atoms.rotate(-j, (D3x-D2x, D3y-D2y, D3z-D2z), center=(D2x, D2y, D2z))
        rotate_atoms += non_rotate_atoms

        ase.io.write(f"{j}degree.xyz", rotate_atoms)
        # Delete all the intermediate file
        file_list = ["non_rotate.xyz","rotate.xyz","temp_rotate.xyz","temp_non_rotate.xyz"]
        for f in file_list:
            os.remove(f)