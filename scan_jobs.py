import getpass
import re
import subprocess
import os


all_jobs = subprocess.run(['scontrol', 'show', 'jobid', '-dd'], stdout=subprocess.PIPE).stdout.decode("utf-8")
all_jobs = all_jobs.split('\n\n')
outputs = []
jobids = []
for job in all_jobs:
    # only keep jobs that belong to current user
    userid = getpass.getuser()
    #userid = 'tczorro'
    if 'UserId={}'.format(userid) not in job:
        continue
    # filter out jobs by work directory
    if "WorkDir=/blue/rmirandaquintana/" not in job:
        continue
    # skip pending jobs
    #if "JobState=PENDING" not in job:
    #    continue
    # FIXME: if jobstate is pending, then an error will be raised (JobState=PENDING vs JobState=RUNNING) 
    #if not re.search(r'apg\w\.\w', job):
    #    continue
    # get output
    outputs.append(re.search(r'StdOut=(.+)\n', job).group(1))
    # get jobid
    jobids.append(re.search(r'JobId=(\d+) ', job).group(1))

# check output
def extract_sigma(line):
    re_sigma = re.search("^\s+\d+\s+\d+\s+[-\d\.\+e]+\s+[-\d\.\+e]+\s+([-\d\.\+e]+)\s+", line)
    if re_sigma:
        return re_sigma.group(1)

counter = 0
print(len(outputs))

#seen = {}
#dupes = []
#for x, y in zip(outputs, jobids):
#    if x not in seen:
#        seen[x] = 1
#    else:
#        if seen[x] >= 1:
#            dupes.append((x, y))
#        seen[x] += 1
#print(dupes)


#for output, jobid in sorted(zip(outputs, jobids), key=lambda x: x[0]):
for output, jobid in zip(outputs, jobids):
    pass
    
    print(output, jobid)
    subprocess.run(['scancel', jobid])
    #if 'ccsdt/' in output and 'h10' in output:
    #    print(output, jobid)
    #    #subprocess.run(['scancel', jobid])
    #with open(output, 'r') as f:
    #    content = f.read()
    #    if "Function evaluations" in content:
    #        print(output, jobid)
    #        # subprocess.run(['scancel', jobid])

    #subprocess.run(['scancel', jobid])

    #if 'oo-' in output and 'h10' not in output:
    #    print(output, jobid)
    #    subprocess.run(['scancel', jobid])
    #if 'ccsdt' in output:
    #    print(output, jobid)
    #    subprocess.run(['scancel', jobid])

    #i = os.path.join(os.path.split(output)[0], 'calculate.py')
    #continue
    #if 'apg' not in output:
    #    subprocess.run(['sed', '-r', 's/from wfns.solver.equation import cma/from wfns.solver.equation import cma\\nfrom wfns.upgrades import speedup_apg, speedup_sd, speedup_sign/', i, '-i'])
    #else:
    #    subprocess.run(['sed', '-r', 's/from wfns.solver.equation import cma/from wfns.solver.equation import cma\\nfrom wfns.upgrades import speedup_sd, speedup_sign/', i, '-i'])

    #with open(output, 'r') as f:
    #    lines = f.readlines()
    #if any('solutions rejected' in line for line in lines):
    #    print(output, jobid)
    #    #subprocess.run(['scancel', jobid])
        
    #sigma = None
    #for line in lines:
    #    temp_sigma = extract_sigma(line)
    #    if temp_sigma:
    #        sigma = temp_sigma
    ## cancel jobs with low enough sigma
    #if sigma is not None and float(sigma) < 1e-5:
    #    print(output)
    #    #print(jobid)
    #    #subprocess.run(['scancel', jobid])
    #    counter += 1
    #else:
    #    pass
    #    #print(output)
    #    #counter += 1
print(counter)
