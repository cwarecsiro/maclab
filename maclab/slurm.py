import subprocess, getpass
from subprocess import call
import datetime
import sys, os, re

class Jid:
    
    """  Class for slurm job ids as returned from an sbatch call """
    def __init__(self, jid):
        self.jid = jid.replace('\n', '')
        self.jid = int(str.split(jid.replace('\n', ''), ' ')[3])
        self.type = 'sbatch job ID'
    
    @property
    def scancel(self):
        """ simple wrapper around slurm cancel all jobs for user function """
        subprocess.run(["scancel", jobid])
    
    @property
    def log(self):
        jid = str(self.jid)
        scontrol = subprocess.run(['scontrol', 'show', 'job', jid],stdout=subprocess.PIPE).stdout.decode('utf-8')
        if len(scontrol):
            scontrol = scontrol.split('\n')
            stdout = []
            for i in scontrol:
                if(re.search('StdOut', i)):
                    stdout.append(i)
            tofilepath = re.sub(r'|'.join((r'StdOut=', r' ')), '', stdout[0])
            f = open(tofilepath, 'r')
            for line in f:
                print(line)
        else:
            try_default = os.getcwd() + '/slurm-' + jid + '.out'
            try:
                f = open(try_default, 'r')
                for line in f:
                    print(line)
            except:
                print('Job has finished and probably has a custom location for the log file which is on the TODO list...')
        
    @property
    def stats (self):
        JOBID = '--jobs=' + str(self.jid)
        return(subprocess.run(
                ['sacct', JOBID, '--format=User,JobID,account,Timelimit,elapsed,ReqMem,MaxRss,ExitCode'], 
                stdout=subprocess.PIPE).stdout.decode('utf-8'))


def sbatch(ex_script, jobname = None, sys_arg = None):
    """ 
    Summary
    -------
    Python interface to slurm sbatch cmd line function
    
    Args
    ----
    ex_script: executable file (e.g. .sh) to pass to sbacth
    jobname: string to pass as a job name
    sys_arg: dict of named arguments to pass to the program being run 
    
    Returns
    -------
    JOBID (class Jid)
    """
    call = []
    if jobname is not None:
        call.append("--job-name=" + jobname)
    if sys_arg is not None:
        args_list = [i for i in sys_arg.items()]
        arg_name = args_list[0][0]
        arg_val = args_list[0][1]
        call.append("--export=" + str(arg_name) + "=" + str(arg_val))
    call.append(ex_script)    
    if sys_arg is not None:
        call.append(str(arg_val))
    make_call = subprocess.run(["sbatch"] + call, stdout=subprocess.PIPE).stdout.decode('utf-8')
    return(Jid(make_call))


def squeue (user = None, jid = None):
    """ simple wrapper around slurm queue function """
    if user is None:
        user = getpass.getuser()
    else:
        user = str(user)
    if jid is not None:
        jid = str(jid)
        queue = subprocess.run(['squeue', '-u', user, '-j', jid], stdout=subprocess.PIPE).stdout.decode('utf-8')
        print(queue)
    else:
        queue = subprocess.run(['squeue', '-u', user], stdout=subprocess.PIPE).stdout.decode('utf-8')
        queue = queue.split('\n')
        for line in queue:
            print(line)
        
def scancel (JOBID = None):
    """ cancel all jobs for current user """
    if JOBID is None:
        cancel_ok = input('Warning! This will cancel ALL jobs for ' + 
                              getpass.getuser() + '. Type y/n to cancel or quit function')
        if cancel_ok == 'n':
            return
        else:
            subprocess.run(["scancel", "-u", "$USER"])
    else:
        subprocess.run(["scancel", JOBID])

def config(exfile, mem=48, time='00:59:00', ntasks=1, nodes=1, env='python/3.6.1', jobname=None, stdout = None, **args):
    if jobname is None:
        now = datetime.datetime.now()
        now = now.strftime("%Y-%m-%d")
        jobname = 'job_' + now
        temp_dst = '/home/' + getpass.getuser() + '/tmp_sh/' + jobname + '.sh'
    else:
        temp_dst = '/home/' + getpass.getuser() + '/tmp_sh/' + jobname + '.sh'
    try:
        f = open(temp_dst, 'w')
    except FileNotFoundError:
        os.mkdir('/home/' + getpass.getuser() + '/tmp_sh')
        f = open(temp_dst, 'w')
    
    f.write('#!/bin/bash' + '\n')
    f.write('#SBATCH --job-name=' + jobname + '\n')
    f.write('#SBATCH --nodes=' + str(nodes) + '\n')
    f.write('#SBATCH --ntasks-per-node=' + str(ntasks)+ '\n')
    f.write('#SBATCH --mem=' + str(mem) + 'GB' + '\n')
    f.write('#SBATCH --time=' + time + '\n')
    if stdout is not None:
        f.write('#SBATCH --output=%s' %stdout + '\n')
    f.write('\n')
    f.write('echo "Launching job"' + '\n')
    f.write('module load ' + env + '\n')
    if args:
        argv_args = []
        for key, val in args.items():
            argv_args.append(str(val))
        f.write('python ' + exfile + ' ' + ' '.join(argv_args) +'\n') 
    else:
        f.write('python ' + exfile + '\n') 
    f.write('\n')
    f.write('exit_code=$?' + '\n')
    f.write('if [ "$exit_code" -ne 0 ]; then ' + '\n')
    f.write('\t' + 'echo "----------------------------------------"' + '\n')
    f.write('\t' + 'echo "<< JOB FAIL! >> Exit code was $exit_code"' + '\n')
    f.write('else' + '\n')
    f.write('\t' + 'echo "Done!"' + '\n')
    f.write('fi')
    f.close()
    return(temp_dst)