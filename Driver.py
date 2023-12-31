"""
Written by Naveen Prakash
"""

mat2_thk = [0.5, 0.54, 0.58, 0.62, 0.66, 0.7, 0.74, 0.78, 0.82, 0.86, 0.9]
mat2_E   = [65e3, 70e3, 75e3, 80e3, 85e3, 90e3, 95e3, 100e3, 105e3]

from tempfile import mkstemp
from shutil import move, copymode
from os import fdopen, remove
from subprocess import check_output

def replace(file_path, pattern, subst):
    #Create temp file
    fh, abs_path = mkstemp()
    with fdopen(fh,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                if pattern in line:
                    new_file.write(subst)
                else:
                    new_file.write(line)
    #Copy the file permissions from the old file to the new file
    copymode(file_path, abs_path)
    #Remove original file
    remove(file_path)
    #Move new file
    move(abs_path, file_path)

file_path = "./Script_BiMaterial_Master.py"

count = 0
for i in mat2_thk:
    pattern = "mat2_thk="
    subst = "mat2_thk=" + str(i) + "\n"
    replace(file_path, pattern, subst)
    for j in mat2_E:
        pattern = "mat2_E="
        subst = "mat2_E=" + str(j) + "\n"
        replace(file_path, pattern, subst)

        new_cte = -1.3826E-19*pow(j,3) + 3.4898E-14*pow(j,2) - 2.9381E-09*j + 8.4091E-05 

        pattern = "mat2_cte="
        subst = "mat2_cte=" + str(new_cte) + "\n"
        replace(file_path, pattern, subst)
        
        count = count + 1
        print("Running case #" + str(count) + " ----> " + "Thk: " + str(i) + ", Mod: " +  str(j))
        
        check_output("abq2022le cae nogui=Script_BiMaterial_Master.py", shell=True)

print("Driver script done!")