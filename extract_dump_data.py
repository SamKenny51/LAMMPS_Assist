import numpy as np
import random
import matplotlib.pyplot as plt
import os

def get_data_from_line(line, timeStep):
    splitLine = line[:-1].split(" ")
    return [timeStep]+[float(datum) for datum in splitLine]

class LAMMPSDumpDataStream:
    def __init__(self, fileName):
        self.fileName = fileName
        self.extract_header_info(verbose=False)
        self.currentData = []
        self.filePosition = dict()
        self.hookFile = fileName + "_hook.dat"
        if os.path.isfile(self.hookFile):
            self.retrieve_position_hooks(self.hookFile)
    
    def extract_header_info(self, verbose=True):
        with open(self.fileName, 'r') as f:
            counter = 0
            headerInfo = ''
            for line in f:
                if verbose:
                    if counter < 9:
                        headerInfo += line
                    else:
                        return headerInfo
                else:
                    if counter == 0:
                        self.firstLine = line
                    if counter == 1:
                        self.timeStep = int(line)
                    if counter == 3:
                        self.nAtoms = int(line)
                    if counter == 4:
                        self.boundaryConditions = line[:-1].split(" ")[3:]
                    if counter == 5:
                        self.xBounds = [float(num) for num in line.split(" ")]
                    if counter == 6:
                        self.yBounds = [float(num) for num in line.split(" ")]
                    if counter == 7:
                        self.zBounds = [float(num) for num in line.split(" ")]
                    if counter == 8:
                        headerKeys = ["time"]+line[:-1].split(" ")[2:]
                        self.headerMap = {headerKeys[i]:i for i in range(len(headerKeys))}
                    if counter > 0 and line == self.firstLine:
                        self.blockLength = counter - 1
                        break
                counter += 1

    def retrieve_position_hooks(self, hookFile):
        with open(hookFile, 'r') as f:
            for line in f:
                hook = [int(h) for h in line.split(",")]
                self.filePosition[hook[0]] = hook[1]

    def save_position_hook(self, timeStep, seekPosition):
        with open(self.hookFile, 'a') as f:
            f.write(f"{timeStep}, {seekPosition}\n")

    def extract_time_step(self, targetTimeStep, peekFlag=False, debug=False):
        self.currentData = None
        counter = 0
        dataArray = []
        timeStep = -1
        timeFlag = False
        searchMode = True
        
        print('Reading file...')
        with open(self.fileName, 'r') as f:
            if self.filePosition:
                print("Setting starting point...")
                if targetTimeStep in self.filePosition.keys():
                    print(f"Setting file position to cached timestep {targetTimeStep}: {self.filePosition[targetTimeStep]}")
                    timeStep = targetTimeStep - 1
                    f.seek(self.filePosition[targetTimeStep])
                    searchMode = False
                elif targetTimeStep == -1:
                    print("Reading entire file...")
                    timeStep = max(self.filePosition.keys())
                    f.seek(self.filePosition[timeStep])
                else:
                    for k in self.filePosition.keys():
                        if timeStep < targetTimeStep:
                            timeStep = k
                    print(f"Reading forward from cached timestep {timeStep}: {self.filePosition[timeStep]}")
                    f.seek(self.filePosition[timeStep])
            else:
                print('Searching file...')
                self.filePosition = dict()
                
            tempFilePos = f.tell()
            line = f.readline()
            while line:
                if debug: print(line)
                    
                if line == self.firstLine:
                    counter = 0
                    if debug: print(timeStep)
                    timeStep += 1
                    if searchMode:
                        # Save current position for future use
                        self.filePosition[timeStep] = tempFilePos
                        self.save_position_hook(timeStep, tempFilePos)
                        
                if targetTimeStep == -1:
                    tempFilePos = f.tell()
                    line = f.readline()
                    continue
                                        
                if timeStep == targetTimeStep and not timeFlag:
                    print(f"Timestep {timeStep} reached...")
                    timeFlag = True

                if timeFlag and counter == 1:
                    print(line)
                    
                if timeStep > targetTimeStep:
                    break
                
                if counter == 3:
                    self.nAtoms = int(line)
                                    
                if counter >= 9 and timeStep == targetTimeStep:
                    dataArray.append(get_data_from_line(line, timeStep))
                    
                counter += 1
                
                if searchMode:
                    tempFilePos = f.tell()
                    
                line = f.readline()
                
        self.currentData = np.array(dataArray)
        if peekFlag:
            return self.currentData[0,:]
        
    def bounds_by_type(self, typeIDs):
        bounds = {"x": [np.inf, -np.inf],
                  "y": [np.inf, -np.inf],
                  "z": [np.inf, -np.inf]}
        
        def expand_bounds(bounds, direction, dataline):
            # Update min
            if line[self.headerMap[direction]] < bounds[0]:
                bounds[0] = line[self.headerMap[direction]]
            # Update max
            if line[self.headerMap[direction]] > bounds[1]:
                bounds[1] = line[self.headerMap[direction]]
            return bounds
        
        if self.currentData is None:
            return bounds
        
        typeCheck = True
        for line in self.currentData:
            if any([int(line[self.headerMap["type"]]) == typeID for typeID in typeIDs]):
                typeCheck = False
                for direction in ["x","y","z"]:
                    bounds[direction] = expand_bounds(bounds[direction], direction, line)
        if typeCheck:
            raise ValueError("No atoms of given types. Exiting...")
        return bounds

    def bounding_vol_by_type(self, typeIDs):
        bounds = self.bounds_by_type(typeIDs)
        vol = (bounds['x'][1] - bounds['x'][0])*(bounds['y'][1] - bounds['y'][0])*(bounds['z'][1] - bounds['z'][0])
        if vol <= 0:
            raise ValueError("Non-positive volume computed. Exiting...")
        return vol

    def get_atom_properties_by_type(self, atomProperties, typeIDs):
        atomProps = {prop: [] for prop in atomProperties}
        for line in self.currentData:
            if any([int(line[self.headerMap["type"]]) == typeID for typeID in typeIDs]):
                for prop in atomProps.keys():
                    atomProps[prop].append(line[self.headerMap[prop]])
        for prop in atomProps.keys():
            atomProps[prop] = np.array(atomProps[prop])
        return atomProps

    def get_atom_properties_by_id(self, atomProperties, atomIDs):
        atomProps = {prop: [] for prop in atomProperties}
        for line in self.currentData:
            if any([int(line[self.headerMap["id"]]) == atomID for atomID in atomIDs]):
                for prop in atomProps.keys():
                    atomProps[prop].append(line[self.headerMap[prop]])
        for prop in atomProps.keys():
            atomProps[prop] = np.array(atomProps[prop])
        return atomProps
            
    def count_if_type(self, typeIDs):
        count = 0
        for line in self.currentData:
            if any([int(line[self.headerMap["type"]]) == typeID for typeID in typeIDs]):
                count += 1
        return count

    def sum_if(self, sumProp, prop, bounds):
        if sumProp not in headerMap.keys():
            raise ValueError("Conditional Sum of " + sumProp + " failed. No such property in header.")
        for p in props:
            if p not in headerMap.keys():
                raise ValueError("Conditional Sum of " + sumProp + " failed. " + p + " property not in header.")
        return np.sum(self.currentData[bounds[0]<=self.currentData[:,headerMap[prop]]<=bounds[1]][:,headerMap[sumProp]])

    

def if_between(array, minBound, maxBound):
    res = []
    for el in array:
        if minBound < el <= maxBound:
            res.append(True)
        else:
            res.append(False)
    return np.array(res)

def dir_bounds_by_type(dataStream, typeIDs, direction):
    bounds = [np.inf, -np.inf]
    for k in dataStream.filePosition.keys():
        print("Bounds")
        try:
            dataStream.extract_time_step(k)
        except TypeError as e:
            print("EOF reached...")
            break
        
        tempBounds = dataStream.bounds_by_type(typeIDs)[direction]
        
        if tempBounds[0] < bounds[0]:
                bounds[0] = tempBounds[0]
        if tempBounds[1] > bounds[1]:
                bounds[1] = tempBounds[1]
                
    return bounds

def box_volume(xBounds, yBounds, zBounds):
    return (xBounds[1]-xBounds[0])*(xBounds[1]-xBounds[0])*(zBounds[1]-zBounds[0])

def procedure1(dataStream, Nbins=50):
    zbounds = dir_bounds_by_type(dataStream, [2,3], 'z')
                
    zbins = np.linspace(zbounds[0], zbounds[1], num=Nbins)

    vol = box_volume(dataStream.xBounds, dataStream.yBounds, zbins[0:2])
    print("per bin volume: ",vol)

    VX, VY, VZ = [],[],[]
    FX, FY, FZ = [],[],[]
    NumArr = []
    
    for k in dataStream.filePosition.keys():
        print("Accumulator")
        sumVX, sumVY, sumVZ = np.zeros(zbins.shape[0]), np.zeros(zbins.shape[0]), np.zeros(zbins.shape[0])
        sumFX, sumFY, sumFZ = np.zeros(zbins.shape[0]), np.zeros(zbins.shape[0]), np.zeros(zbins.shape[0])
        Num =np.zeros(zbins.shape[0])
        try:
            dataStream.extract_time_step(k)
        except TypeError as e:
            break

        if dataStream.currentData is None:
            print("EOF reached...")
            break
        
        for line in dataStream.currentData:
            if 1.5<line[dataStream.headerMap['type']]<3.1:
                # Scan for the correct bin
                for i in range(zbins.shape[0]):
                    if i == 49:
                        break
                    if zbins[i]<=line[dataStream.headerMap['z']]<=zbins[i+1]:
                        break
                # accumulate in the bins
                sumVX[i] += line[dataStream.headerMap['vx']]
                sumVY[i] += line[dataStream.headerMap['vy']]
                sumVZ[i] += line[dataStream.headerMap['vz']]
                sumFX[i] += line[dataStream.headerMap['fx']]
                sumFY[i] += line[dataStream.headerMap['fy']]
                sumFZ[i] += line[dataStream.headerMap['fz']]
                Num[i] += 1
        print(sumVX)
        VX.append(sumVX)
        VY.append(sumVY)
        VZ.append(sumVZ)
        FX.append(sumFX)
        FY.append(sumFY)
        FZ.append(sumFZ)
        NumArr.append(Num)
    return zbins, VX, VY, VZ, FX, FY, FZ, NumArr, vol

def diffusion_coefficient(dataStream, nAtoms=100):
    nx = np.zeros(nAtoms)
    ny = np.zeros(nAtoms)
    nz = np.zeros(nAtoms)
    for k in dataStream.filePosition.keys():
        try:
            dataStream.extract_time_step(k)
        except TypeError as e:
            break

        if dataStream.currentData is None:
            print("EOF reached...")
            break

        if k == 0:
            atomIDs = dataStream.get_atom_properties_by_type(['id'], [2,3])
            atomIDs = list(map(int, atomIDs['id']))
            random.shuffle(atomIDs)
            atomIDs = atomIDs[:nAtoms]
            msds = np.zeros(len(dataStream.filePosition.keys()))
            r0 = dataStream.get_atom_properties_by_id(['x','y','z'], atomIDs)
            rprev = r0
            
            avgVel = np.zeros(len(dataStream.filePosition.keys()))
            v = dataStream.get_atom_properties_by_id(['vx', 'vy', 'vz'], atomIDs)
            avgVel[k] = np.mean(np.sqrt(v['vx']**2. + v['vy']**2. + v['vz']**2.))
        else:
            rnew = dataStream.get_atom_properties_by_id(['x','y','z'], atomIDs)
            for i in range(nAtoms):
                if rnew['x'][i] - rprev['x'][i] > dataStream.xBounds[1]/2:
                    nx[i] -= 1
                elif rnew['x'][i] - rprev['x'][i] < -dataStream.xBounds[1]/2:
                    nx[i] += 1

                if rnew['y'][i] - rprev['y'][i] > dataStream.yBounds[1]/2:
                    ny[i] -= 1
                elif rnew['y'][i] - rprev['y'][i] < -dataStream.yBounds[1]/2:
                    ny[i] += 1

                if rnew['z'][i] - rprev['z'][i] > dataStream.zBounds[1]/2:
                    nz[i] -= 1
                elif rnew['z'][i] - rprev['z'][i] < -dataStream.zBounds[1]/2:
                    nz[i] += 1

                rnew['x'][i] += nx[i] * dataStream.xBounds[1]
                rnew['y'][i] += ny[i] * dataStream.yBounds[1]
                rnew['z'][i] += nz[i] * dataStream.zBounds[1]
                
            msds[k] = np.mean((rnew['x'] - r0['x'])**2 + (rnew['y'] - r0['y'])**2 + (rnew['z'] - r0['z'])**2)
            rprev = rnew

            v = dataStream.get_atom_properties_by_id(['vx', 'vy', 'vz'], atomIDs)
            avgVel[k] = np.mean(np.sqrt(v['vx']**2. + v['vy']**2. + v['vz']**2.))
            
    return msds, avgVel
            
'''def save_bins_as_columns(arr):
    with open("outfile.dat", "w") as f:
	for row in arr:
            f.write((','.join(['%10.6f']*row.size)+'\n') % tuple(row))'''

def expression_selectionAND(dataStream, attrs, ops, vals):
    res = np.copy(dataStream.currentData)
    for i in range(len(attrs)):
        if ops[i] == '>':
            res = res[res[:, dataStream.headerMap[attrs[i]]] > vals[i], :]
        if ops[i] == '<':
            res = res[res[:, dataStream.headerMap[attrs[i]]] < vals[i], :]
        if ops[i] == '>=':
            res = res[res[:, dataStream.headerMap[attrs[i]]] >= vals[i], :]
        if ops[i] == '<=':
            res = res[res[:, dataStream.headerMap[attrs[i]]] <= vals[i], :]
        if ops[i] == '==':
            res = res[res[:, dataStream.headerMap[attrs[i]]] == vals[i], :]
        if ops[i] == '!=':
            res = res[res[:, dataStream.headerMap[attrs[i]]] != vals[i], :]
    return res

def atoms_from_molIDs(dataStream, molIDs):
    res = []
    for molID in molIDs:
        atom = expression_selectionAND(dataStream, ['mol'], ['=='], [molID])
        res.extend(list(atom))
    return np.array(res)

# specific to project =====================================================

potentialParams_Au_EtOH = {
5:[74534.5,-4,	0.0,	111996,	-2.2,	2.7],
4:[29633.0,	-2.6,	0.4,	17852.9,-1.2,	3.5],
3:[18594.0,	-3.7,	-1.2,	85932.8,	-3.6,	8.7],
2:[76040.5,	-3.3,	-0.5,	63929.7,	-1.4,	4.5]}

# eps_i, sigma_i, q_i // 2=ch3, 3=ch2, 4=O, 5=H
potentialParams_LJ_coul_long = {
    2:[0.008444935324472082, 3.75, 0],
    3:[0.0039639492339358755, 3.95, 0.265],
    4:[0.008014071277305138, 3.02, -0.7],
    5:[0, 0, 0.435]}

def LJCoulLong_Force(rvec, r, type1, type2, params):
    ep = np.sqrt(params[type1][0]*params[type2][0])
    sigma = 0.5 * (params[type1][1] + params[type2][1]) 
    C = 14.39964548 # eV * Angstrom / e^2 (Coulomb constant in metal units)
    if r < 14:
        return np.zeros((3,))
    return (C * params[type1] * params[type2] / r**2 + 48*sigma**12*ep*r**-13 - 24*sigma**6*ep*r**-7) * ( rvec/r )

def ModifiedMorseU(r, params):
    return 0.01036*params[0] * np.exp(params[1]*(r + params[2])) - params[3] * np.exp(params[4]*(r + params[5]))

def ModifiedMorseF(rvec, r, params):
    return 0.01036*(params[3]*params[4]*np.exp(params[4]*(r + params[5])) - params[0]*params[1]*np.exp(params[1]*(r + params[2]))) * (rvec/r)

def mod_pbc(rvec, xBounds, yBounds, zBounds):
    if rvec[0] > xBounds[1]/2:
        rvec[0] -= xBounds[1]
    elif rvec[0] < -xBounds[1]/2:
        rvec[0] += xBounds[1]
        
    if rvec[1] > yBounds[1]/2:
        rvec[1] -= yBounds[1]
    elif rvec[1] < -yBounds[1]/2:
        rvec[1] += yBounds[1]
        
    if rvec[2] > zBounds[1]/2:
        rvec[2] -= zBounds[1]
    elif rvec[2] < -zBounds[1]/2:
        rvec[2] += zBounds[1]
        
    return rvec

if __name__ == "__main__":
    dataStream = LAMMPSDumpDataStream('dump.AuEtOHThermalizedNoForce')
    surfDistCutoff = 6
    xbs = dataStream.xBounds
    ybs = dataStream.yBounds
    zbs = dataStream.zBounds
    #for k in dataStream.filePosition.keys():
    dataStream.extract_time_step(0)
    lowerSurfZ = np.mean(expression_selectionAND(dataStream, ['z','z','type'], ['>','<','=='], [22, 24, 1])[:,6])
    upperSurfZ = np.mean(expression_selectionAND(dataStream, ['z','z','type'], ['>','<','=='], [63, 65, 1])[:,6])
    #print(lowerSurfZ, upperSurfZ)
    
    lowerSurfaceMolIDs = set(expression_selectionAND(dataStream, ['z','z','type'], ['>','<','>'], [lowerSurfZ-surfDistCutoff, lowerSurfZ, 1])[:,3])
    upperSurfaceMolIDs = set(expression_selectionAND(dataStream, ['z','z','type'], ['>','<','>'], [upperSurfZ, upperSurfZ+surfDistCutoff, 1])[:,3])
    surfaceMolIDs = lowerSurfaceMolIDs.union(upperSurfaceMolIDs)
    surfaceMolIDs = list(map(int, surfaceMolIDs))
    
    print(f"Number of molecules in selection = {len(surfaceMolIDs)}")
    
    surfAtoms = atoms_from_molIDs(dataStream, surfaceMolIDs)
    print(surfAtoms.shape)
    #print(f"Number of atom IDs in selection = {len(atomIDs)}")
    #rvec = modPeriodicBox(rvec)
    # get atom ids of Au atoms
    AuAtoms = expression_selectionAND(dataStream, ['type'], ['=='], [1])

    forceSurfaceAtoms = np.zeros((surfAtoms.shape[0],3))
    for i in range(surfAtoms.shape[0]):
        for j in range(AuAtoms.shape[0]):
            rvec = surfAtoms[i,4:7] - AuAtoms[j,4:7]
            rvec = mod_pbc(rvec, xbs, ybs, zbs)
            r = np.linalg.norm(rvec)
            if r < 14:
                forceSurfaceAtoms[i,:] += ModifiedMorseF(rvec, r, potentialParams_Au_EtOH[surfAtoms[i,2]])
    
        
