import subprocess


def start_pfam_scan(pfam_scan_script_location, pfam_input_location, pfam_output_location, hmmer3_data_location):
    command1 = ['perl',
                pfam_scan_script_location,
                '-fasta',
                pfam_input_location,
                '-outfile',
                pfam_output_location,
                '-dir',
                hmmer3_data_location]
    p1 = subprocess.Popen(command1)
    exit_code = p1.wait()
    return exit_code
