import os
import subprocess

PULSE_PATH = os.getcwd()


def run_samtools(cell_line_for_samtools):
    print(cell_line_for_samtools)
    paths_to_create = [r'./output/for_preprocess/' + cell_line_for_samtools,
                       r'./output/for_preprocess/' + cell_line_for_samtools + r'/no-rg',
                       r'./output/for_preprocess/' + cell_line_for_samtools + r'/cufflinks_output']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)

    command1 = ['samtools', 'view', '-h', PULSE_PATH + '/input/cell_lines/' + cell_line_for_samtools]
    command2 = ['grep',  '-v', '@RG']
    command3 = ['samtools', 'view', '-Sbh', '-']
    output_file = open(PULSE_PATH + '/output/for_preprocess/' +
                       cell_line_for_samtools + '/no-rg/' + cell_line_for_samtools, "w")

    print "Running: view"
    p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
    print "Running: grep"
    p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
    print "Running view2"
    p3 = subprocess.Popen(command3, stdin=p2.stdout, stdout=output_file)
    exit_codes = [p.wait() for p in p1, p2, p3]
    return exit_codes


def run_cufflinks(cell_line_name):
    print "Processing cell line: " + cell_line_name
    output_file = open(PULSE_PATH + '/output/for_preprocess/' + cell_line_name + '/transcript.gtf', 'w')
    command1 = ['cufflinks', '-o', PULSE_PATH + '/output/for_preprocess/cufflinks_output/' + cell_line_name, '-g',
                PULSE_PATH + '/input/Homo_sapiens.GRCh37.70.gtf', PULSE_PATH + '/output/for_preprocess/' +
                cell_line_name + '/no-rg/' + cell_line_name]
    p1 = subprocess.Popen(command1, stdout=output_file)
    p1.stdout.close()


if __name__ == "__main__":
    all_cell_lines = os.listdir("./input/cell_lines")
    for cell_line in all_cell_lines:
        print "Running samtools for: " + cell_line
        print run_samtools(cell_line)
        run_cufflinks(cell_line)
