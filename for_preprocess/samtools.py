import os
import subprocess

PULSE_PATH = os.getcwd()


def create_paths_for_cell_line(cell_line_for_samtools):
    """
    Creates all necessary paths for preprocessing:

    Base folder for cell line:
    - ../output/for_preprocess/<CELL_LINE_NAME>/

    samtools output:
    - ../output/for_preprocess/<CELL_LINE_NAME>/no-rg/

    cufflinks_output:
    - ../output/for_preprocess/<CELL_LINE_NAME>/cufflinks_output/

    :param cell_line_for_samtools:
    :return:
    """
    print(cell_line_for_samtools)
    paths_to_create = [r'../output/for_preprocess/' + cell_line_for_samtools,
                       r'../output/for_preprocess/' + cell_line_for_samtools + r'/no-rg',
                       r'../output/for_preprocess/' + cell_line_for_samtools + r'/cufflinks_output']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)


def run_samtools(pulse_path, cell_line_for_samtools):

    command1 = ['samtools', 'view', '-h', pulse_path + '/input/cell_lines/' + cell_line_for_samtools]
    command2 = ['grep',  '-v', '@RG']
    command3 = ['samtools', 'view', '-Sbh', '-']
    output_file = open(pulse_path + '/output/for_preprocess/' +
                       cell_line_for_samtools + '/no-rg/' + cell_line_for_samtools, "w")

    print "Running: view"
    p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
    print "Running: grep"
    p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
    print "Running view2"
    p3 = subprocess.Popen(command3, stdin=p2.stdout, stdout=output_file)
    exit_codes = [p.wait() for p in p1, p2, p3]
    return exit_codes
