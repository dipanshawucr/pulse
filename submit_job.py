# submitjob 2 -m 3 -c 4 ./script.pl
# runs for two hours, using 3 gigs of RAM, using 4 cores

import subprocess


def submit_new_job(script_name, list_of_args=[]):
    """
    Submits a new job to beagle with a max runtime for 48 hours,
    with 4 GB of memory, and 4 cores per job.
    Returns 0 upon task completion.

    :param script_name:
    :param list_of_args:
    :return:
    """

    command1 = ['submitjob',
                '48',
                '-m',
                '4',
                '-c',
                '4',
                script_name] + list_of_args
    p1 = subprocess.Popen(command1)
    exit_code = p1.wait()
    return exit_code
