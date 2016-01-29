import subprocess


def run_cufflinks(pulse_path, cell_line_name):
    output_file = open(pulse_path + '/output/for_preprocess/' + cell_line_name + '/transcripts.gtf', 'w')
    command1 = ['cufflinks', '-o', pulse_path + '/output/for_preprocess/' + cell_line_name + '/cufflinks_output', '-g',
                pulse_path + '/input/Homo_sapiens.GRCh37.70.gtf', pulse_path + '/output/for_preprocess/' +
                cell_line_name + '/no-rg/' + cell_line_name]
    p1 = subprocess.Popen(command1, stdout=output_file)
    exit_codes = p1.wait()
    return exit_codes
