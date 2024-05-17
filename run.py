import subprocess
import time
import os
import sys

timeout = 1800
testdir = "test/" + sys.argv[1] + "/"
instdir = testdir + "instances/"
outpdir = testdir + "output/"
weberk = "build/weberknecht"

if __name__ == "__main__":
    succeded = []
    failed = []
    n = 0
    instances = os.listdir(instdir)
    instances.sort(key=lambda k: int(k.split(".")[0]))
    for filename in instances:
        print(f"running {filename}")
        n += 1
        with open(instdir + filename) as i, open(outpdir + filename[:-2] + "log", "w") as o:
            try:
                subprocess.call([weberk], stdin=i, stdout=o, timeout=timeout)
                succeded.append(filename)
            except subprocess.TimeoutExpired:
                failed.append(filename)
    
    output = f"number of instances: {n}\n"
    output += f"succeded {len(succeded) / n * 100 : 3.1f}%:\n{succeded}\n"
    output += f"failed   {len(failed) / n * 100 : 3.1f}%:\n{failed}\n"
    with open(outpdir + "run.log", "w") as o:
        o.write(output)
    print(output)
    