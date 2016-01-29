import subprocess


def blast_events(uniprot_file_location, events_file_location, output_location):
    """

    :param uniprot_file_location:
    :param events_file_location:
    :param output_location:
    :return:
    """
    command1 = ['blastall',
                '-p',
                'blastx',
                '-d',
                uniprot_file_location,
                '-i',
                events_file_location,
                '-o',
                output_location,
                '-v1',
                '-b1',
                '-m8']
    p1 = subprocess.Popen(command1)
    exit_code = p1.wait()
    return exit_code
