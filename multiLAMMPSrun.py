import os
import subprocess
import datetime

'''files = [file for file in os.listdir() if file[:2] == 'in']
print(files)
path = "D:\\MolecularDynamics\\handsOnFiles\\argonMultiTest\\"'''

path = "D:\\MolecularDynamics\\AAAA NanoscaleSliding\\"
#subpaths = ["AAAA EtOHonAu_molecularFlow\\AAAA EtOHonAu_v2500\\"]"AAAA EtOHonAu_v2500_fixedAdsorb\\",
subpaths = [
    "AAAA NoImpulse\\"]

folders = ["v"+str(i) for i in range(1,4)]

for i, sp in enumerate(subpaths):
    for j, folder in enumerate(folders):
        os.chdir(path+sp+folder)
        file = 'AuSlab+etOH_STRESS_RERUN.in'
        print(f"{file} started at {datetime.datetime.now().strftime('%m/%d/%Y %H:%M:%S')}")
        subprocess.run(['lmp_serial', '-in', file])
        print(os.listdir())
        print(f"{file} finished.")
        os.rename(path+sp+folder+'\\stress.out', path+sp+folder+'\\stress\\stress.out')
        os.rename(path+sp+folder+'\\log.lammps', path+sp+folder+'\\stress\\log.lammps')
        os.chdir('.\\stress')
        print(os.getcwd(), os.listdir())
            
    print(10*"#")
