import subprocess

def run_cufflinks(pulse_path, cell_line_name):
    output_file = open(pulse_path + '/output/for_preprocess/' + cell_line_name + '/transcript.gtf', 'w')
    command1 = ['cufflinks', '-o', pulse_path + '/output/for_preprocess/cufflinks_output/' + cell_line_name, '-g',
                pulse_path + '/input/Homo_sapiens.GRCh37.70.gtf', pulse_path + '/output/for_preprocess/' +
                cell_line_name + '/no-rg/' + cell_line_name]
    p1 = subprocess.Popen(command1, stdout=output_file)
    exit_codes = [p.wait() for p in p1]
    return exit_codes
