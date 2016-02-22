# submitjob 2 -m 3 -c 4 ./script.pl
# runs for two hours, using 3 gigs of RAM, using 4 cores

import subprocess


def submit_new_job(script_name, list_of_args):

    command1 = ['submitjob',
                '48',
                '-m',
                '4',
                '-c',
                '4',
                script_name]
    p1 = subprocess.Popen(command1)
    exit_code = p1.wait()
    return exit_code
