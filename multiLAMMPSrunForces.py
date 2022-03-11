import os
import subprocess
import datetime

'''files = [file for file in os.listdir() if file[:2] == 'in']
print(files)
path = "D:\\MolecularDynamics\\handsOnFiles\\argonMultiTest\\"'''

logfile = 'D:\MolecularDynamics\AAAA NanoscaleSliding\AAAA EtOHonAu_molecularFlow\AAAA EtOHonAu_v2500\\force.log'
fLog = open(logfile, 'a')

path = "D:\\MolecularDynamics\\AAAA NanoscaleSliding\\"
subpaths = ["AAAA EtOHonAu_molecularFlow\\AAAA EtOHonAu_v2500\\",
            "AAAA EtOHonAu_v2500_fixedAdsorb\\",
            "AAAA NoImpulse\\"]

folders = ["v"+str(i) for i in range(1,7)]

files = ['AuSlab+etOH_SurfForce_RERUN.in',
        'AuSlab+etOH_MolecForce_RERUN.in']
dump = ['dump.AuEtOHSurf_Forces',
        'dump.AuEtOHMolec_Forces']

nfiles, limit = 0, 3

for i, sp in enumerate(subpaths):
    for j, folder in enumerate(folders):
        os.chdir(path+sp+folder)

        # Create directories to store data
        absPath = os.getcwd()
        p1 = os.path.join(absPath, "surfForce")
        p2 = os.path.join(absPath, "MolecForce")
        os.mkdir(p1)
        os.mkdir(p2)
        
        # first subpath folder v1 is a short test run.
        if i == 0 and folder == 'v1':
            continue
        
        for f, file in enumerate(files):
            logText = f"{file} started at {datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S')}"
            print(logText)
            fLog.write(logText+'\\n')
            
            subprocess.run(['lmp_serial', '-in', file])
            print(os.listdir())
            print(f"{file} finished.")
            # Move data to created folders
            if f == 0:
                os.rename(path+sp+folder+'\\'+dump[f], path+sp+folder+'\\surfForce\\'+dump[f])
                os.rename(path+sp+folder+'\\log.lammps', path+sp+folder+'\\surfForce\\log.lammps')
                os.chdir('.\\surfForce')
            elif f == 1:
                os.rename(path+sp+folder+'\\'+dump[f], path+sp+folder+'\\MolecForce\\'+dump[f])
                os.rename(path+sp+folder+'\\log.lammps', path+sp+folder+'\\MolecForce\\log.lammps')
                os.chdir('.\\MolecForce')
            print(os.getcwd(), os.listdir())
            nfiles += 1
            if nfiles > limit:
                break
        if nfiles > limit:
            break
    if nfiles > limit:
        break
            
    print(10*"#")

fLog.close()
