import os
from subprocess import call


def run_samtools(cell_line_for_samtools):
    print(cell_line_for_samtools)
    paths_to_create = [r'./output/for_preprocess/' + cell_line_for_samtools,
                       r'./output/for_preprocess/' + cell_line_for_samtools + r'/no-rg',
                       r'./output/for_preprocess/' + cell_line_for_samtools + r'/cufflinks_output']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)
    command = "samtools view -h ./input/cell_lines/" + cell_line_for_samtools + \
              " | grep -v @RG " + "| samtools view -Sbh - > ./output/for_preprocess/" + \
              cell_line_for_samtools + "/no-rg"
    print("Running: " + command)
    call(command)


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

