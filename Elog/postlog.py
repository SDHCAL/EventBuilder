import os
import subprocess
def post(run,daqname):
    subprocess.call(['./Test',"./Geom.xml",str(run),str(daqname)])
