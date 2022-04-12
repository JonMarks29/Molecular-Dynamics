# -*- coding: utf-8 -*-
"""
gro_handler.py 
last edit 24th April 2021

@author: Jonathan Machin
"""

import numpy as np
#import matplotlib.pyplot as plt


# atom class - DO NOT EXPLICITLY INSTANTIATE THIS CLASS IN SCRIPTS
class Atom():
    def __init__(self, resnumb, resname, atomname, atomnumb, x,y,z, vx=None,vy=None,vz=None):
        self.resnumb = resnumb
        self.resname = resname
        self.atomname = atomname
        self.atomnumb = atomnumb
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz

        self.coor = (self.x, self.y, self.z)
        
    def velocity(self):
        return (self.vx, self.vy, self.vz)

        
# residue class - DO NOT EXPLICITLY INSTANTIATE THIS CLASS IN SCRIPTS 
class Residue():
    def __init__(self, atoms, resnumb, resname):
        self.atoms = atoms
        self.resnumb = resnumb
        self.resname = resname
        
    def find_center(self, all_atoms):
        xc,yc,zc = 0,0,0
        for atom in self.atoms:
            x,y,z = all_atoms[atom].coor
            xc += x
            yc += y
            zc += z          
        x,y,z = xc/len(self.atoms), yc/len(self.atoms), zc/len(self.atoms)
        self.center = (x,y,z)
        return (x,y,z)
        
    def map_coor(self, new_center, replace, flip=False, dm=0):
        new_coor = {}
        xc,yc,zc = self.center
        xn,yn,zn = new_center
        for atom in self.atoms:
            x,y,z = replace.atoms[atom].coor
            atomname = replace.atoms[atom].atomname
            xd = round(xn + (xc - x), 3)
            yd = round(yn + (yc - y), 3)
            zd = round(zn + (zc - z), 3)
            
            #flip all atoms around z-axis
            if flip == True:
                zd = round(zd - 2*(zd - zn) +dm, 3)
            
            new_coor[atom] = (xd,yd,zd, atomname)
        
        return new_coor
        
        
        
    
class Gro():
    def __init__(self, file, exclude=[]):
        self.file = file
        self.atoms = {}
        self.residues = {}
        self.resnames = set()
        self.atomnames = set()

        
        self.read(exclude)
        
        
    def read(self, exclude):
        with open(self.file) as f:
            lines = f.readlines()
            
        xcent, ycent, zcent = 0,0,0
        self.file_title = lines[0][:-1]
        self.atom_count = int(lines[1][:-1].strip())
        
        l, box_found = 0, False
        while box_found == False:
            l -= 1
            lb = lines[l]
            if len(lb) > 5:
                box_found = True
                self.box = (float(lb[:-1].split()[0]), float(lb[:-1].split()[1]), float(lb[:-1].split()[2]))
        
        
        print('reading file with title: ', lines[0][:-1])
        print('number of atoms/beads in file: ', lines[1][:-1])
        print('boxsize (xyz):', self.box)
        
        residue = 0
        residue_atoms = []
        
        # initially tries to read gro files with velocity data, if fails reads it without velocity data
        for l in lines[2:-2]:
            if l[5:10].strip() in exclude:
                continue
            
            x,y,z = float(l[20:28].strip()), float(l[28:36].strip()), float(l[36:44].strip()) 
            try:
                self.atoms[int(l[15:20].strip())] = Atom(int(l[0:5].strip()), l[5:10].strip(), l[10:15].strip(), int(l[15:20].strip()), x,y,z, float(l[44:52].strip()), float(l[52:60].strip()), float(l[60:68].strip()))
            except ValueError:
                self.atoms[int(l[15:20].strip())] = Atom(int(l[0:5].strip()), l[5:10].strip(), l[10:15].strip(), int(l[15:20].strip()), x,y,z)
            self.resnames.add(l[5:10].strip())
            self.atomnames.add(l[10:15].strip()) 
            xcent, ycent, zcent = xcent+x, ycent+y, zcent+z  
            
            resnumber = int(l[0:5].strip())
            if residue != resnumber:
                residue = resnumber
                self.residues[residue] = Residue([int(l[15:20].strip())], residue, l[5:10].strip())
                
            elif residue == resnumber:
                self.residues[residue].atoms.append(int(l[15:20].strip()))
                
            
            # reports progress by considering atom number 
            # if int(l[15:20])%10000 == 0:
            #     print('at gro atom numb: ', int(l[15:20]))
                            
        
        #print('all atom bead names:', self.atomnames)
        self.center = (xcent/self.atom_count, ycent/self.atom_count, zcent/self.atom_count)
        for n, res in self.residues.items():
            res.find_center(self.atoms)
        print()
        
        
    def write(self, name, tidy=True):
        if tidy == True:
            self.tidy_numbers()
        
        
        title = self.file_title + ' MODIFIED\n'        
        with open(name, 'w') as f:
            f.write(title)
            f.write(str(len(self.atoms))+'\n')
            
            for n, atom in self.atoms.items():
                if atom.vx == None:
                    atom.vx, atom.vy, atom.vz = 0.0,0.0,0.0
                
                resnumb = "{:>5}".format(str(atom.resnumb))
                resname = "{:<5}".format(str(atom.resname))
                atomname = "{:>5}".format(str(atom.atomname))
                atomnumb = "{:>5}".format(str(atom.atomnumb))
                
                def pad_decimal(i, n=3):
                    i = str(i)
                    while len(i.split('.')[-1]) < n:
                        i += '0'    
                    while len(i.split('.')[-1]) > n:
                        i = i[:-1]
                    return i
 
                x,y,z = atom.coor
                x,y,z = pad_decimal(x), pad_decimal(y), pad_decimal(z)
                
                vx,vy,vz = atom.velocity()
                vx,vy,vz = pad_decimal(vx), pad_decimal(vy), pad_decimal(vz)
                       
                
                x = "{:>8}".format(x)[:8]
                y = "{:>8}".format(y)[:8]
                z = "{:>8}".format(z)[:8]
                vx = "{:>8}".format(vx)[:8]
                vy = "{:>8}".format(vy)[:8]
                vz = "{:>8}".format(vz)[:8]
                
                f.write(resnumb+resname+atomname+atomnumb+x+y+z+vx+vy+vz+'\n')
            
            boxx = "{:<10}".format(str(self.box[0]))
            boxy = "{:<10}".format(str(self.box[1]))
            boxz = "{:<10}".format(str(self.box[2]))
            
            f.write(boxx+boxy+boxz+'\n')
            f.write('\n')
            
        print('written current Gro object to', name)
                    
        
        
    # removes any check single-atom residues closer than cutoff to the given residue
    # IF MULTI-ATOM RESIDUES ARE SUPPLIED (in check) ONLY SOME ATOMS MAY BE RE REMOVED
    def solvent_clash(self, check=['W'], residue=None, cutoff=0.3):
        resatoms = self.residues[residue].atoms
        rxc, ryc, rzc = self.residues[residue].center
        
        remove_residues = set()
        for n, atom in self.atoms.items():
            if atom.resname not in check:
                continue
            ax,ay,az = atom.coor
            if abs(ax-rxc) < 3 and abs(ay-ryc) <  3:
                for resatom in resatoms:
                    x,y,z = self.atoms[resatom].coor
                    
                    if ((ax-x)**2 + (ay-y)**2 + (az-z)**2)**0.5 < cutoff:
                        remove_residues.add(atom.resnumb)
        print(len(remove_residues), 'residues of type', check, 'clash and will be removed')
        
        for ri in remove_residues:
            self.remove_residue(ri)
                        
    
    def residue_count(self):
        all_res = {'protein': 0}
        for resname in self.resnames:
            all_res[resname] = 0
        
        for n, res in self.residues.items():
            all_res[res.resname] += 1
          
        prot = ['TYR', 'ASN', 'VAL', 'LYS', 'ASP', 'MET', 'LEU', 'PRO', 'ARG', 'GLY', 
                'TRP', 'SER', 'THR', 'ILE', 'PHE', 'GLU', 'ALA', 'CYS', 'GLN', 'HIS',
                'protein']
        print()
        for res, count in all_res.items():
            if res in prot:
                all_res['protein'] += count
            else:
                print(res, count)
        print('protein', all_res['protein'],'\n')
        
        return all_res
            
    
    def tidy_numbers(self):
        new_atoms = {}
        new_residues = {}
        
        count = 1
        atom_sort = sorted(self.atoms.keys())
        for n in atom_sort:
            atom = self.atoms[n]
            atom.atomnumb = count
            resnumb = atom.resnumb
            new_atoms[count] = atom
            new_res_atoms = [count if x==n else x for x in self.residues[resnumb].atoms]
            self.residues[resnumb].atoms = new_res_atoms
            count += 1
        self.atoms = new_atoms
        
        count = 1
        residue_sort = sorted(self.residues.keys())
        for n in residue_sort:
            residue = self.residues[n]
            residue.resnumb = count
            for atom in residue.atoms:
                self.atoms[atom].resnumb = count
                
            new_residues[count] = residue
            count += 1
        self.residues = new_residues
        
        print('renumbered resnumbs and atomnumbs')
            
        
    def create_residue(self, new_coor, replace):
        resname = replace.residues[1].resname
        self.resnames.add(resname)        
        
        current_atom_count = len(self.atoms)
        resnumb = len(self.residues) + 2
        residue_atoms = []
        for atom, coor in new_coor.items():
            atomnumb = current_atom_count + atom
            x,y,z, atomname = coor
            self.atoms[atomnumb] = Atom(resnumb, resname, atomname, atomnumb, x,y,z)
            residue_atoms.append(atomnumb)
        self.residues[resnumb] = Residue(residue_atoms, resnumb, resname)
        center = self.residues[resnumb].find_center(self.atoms)
        
        print('created', resname, 'residue', resnumb, 'centered at', center)
        
        return resnumb
            
        
    def remove_residue(self, resnumb):
        atoms = self.residues[resnumb].atoms
        for a in atoms:
            del self.atoms[a]
        del self.residues[resnumb]
        
        
    def replace_residue(self, coor, original, replace, n=1, flip=False, dm=0, solvent_clash=True):
        res_remove = []
        x,y,z = coor
        for r in range(n):
            dist = np.inf
            closest = None
            for resnumb, res in self.residues.items():
                if res.resname != original:
                    continue
                
                xc,yc,zc = res.center
                dist_check = ((xc-x)**2 + (yc-y)**2 + (zc-z)**2)**0.5
                
                if dist_check < dist and resnumb not in res_remove:
                    closest = resnumb
                    dist = dist_check
            res_remove.append(closest)
        
        xc,yc,zc = 0,0,0
        for ri in res_remove:
            r = self.residues[ri]
            x,y,z = r.center
            xc += x
            yc += y
            zc += z
        
        xc, yc, zc = xc/len(res_remove), yc/len(res_remove), zc/len(res_remove)
            
        new_coor = replace.residues[1].map_coor((xc,yc,zc), replace, flip, dm)
        
        new_resnumb = self.create_residue(new_coor, replace)
        print('residues', res_remove, 'of type', original, 'will be removed')
        for ri in res_remove:
            self.remove_residue(ri)
            
        if solvent_clash == True:
            self.solvent_clash(check=['W', 'ION', 'NA+', 'CL-'], residue=new_resnumb, cutoff=0.4)
        
        self.tidy_numbers()
        print()
                    



