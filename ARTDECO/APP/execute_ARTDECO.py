import os

file = "/mnt/csl/work/cristina.gil.diaz/ARTDECO-1.2.0/" \
"dates_ARTDECO.txt"

# Simulations GAC
directory_name = "/mnt/csl/work/cristina.gil.diaz/Dataset/" \
"ARTDECO/INPUTS/EXAMPLES/PAC/LW/LW_"

# Checking if the file exists
if not os.path.isfile(file):
   print("The file "+file+" does not exist.")
   exit(1)

# Reading line by line
file = open(file, "r")
Lines = file.readlines()
for line in Lines:
   # Removing unnecessary characters
   line = line.rstrip('\r\n')
   print('Case: '+line)

   # Full path of the Python command
   path = directory_name+line+'/config'

   # Checking the existence of the path
   if not os.path.exists(path):
      print("The input file "+path+" does not exist.")
      continue
   else:
      # Execution of the Python command
      os.system("python3 Artdeco.py execute "+path)
