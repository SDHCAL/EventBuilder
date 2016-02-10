import os
import subprocess
def post(run,daqname):
    subprocess.call(['./Test',"./1.xml",str(run),str(daqname)])
