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
    subprocess.Popen(command3, stdin=p2.stdout, stdout=output_file)


def run_cufflinks(no_rg_cell_line):
    print "Processing cell line: " + no_rg_cell_line
    output = "./output/transcripts" + no_rg_cell_line
    call("cufflinks", "-o", output, "-g",
         "./input/Homo_sapiens.GRCh37.70.gtf", "")


if __name__ == "__main__":
    all_cell_lines = os.listdir("./input/cell_lines")
    for cell_line in all_cell_lines:
        print "Running samtools for: " + cell_line
        run_samtools(cell_line)
