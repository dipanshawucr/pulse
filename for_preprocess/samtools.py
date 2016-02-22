import sys
import getopt
import subprocess


def run_samtools(pulse_path, cell_line_for_samtools):

    command1 = ['samtools', 'view', '-h', pulse_path + '/input/cell_lines/' + cell_line_for_samtools]
    command2 = ['grep',  '-v', '@RG']
    command3 = ['samtools', 'view', '-Sbh', '-']
    output_file = open(pulse_path + '/output/for_preprocess/' +
                       cell_line_for_samtools + '/no-rg/' + cell_line_for_samtools, "w")

    print("Running: view")
    p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
    print("Running: grep")
    p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
    print("Running view2")
    p3 = subprocess.Popen(command3, stdin=p2.stdout, stdout=output_file)
    exit_codes = [p.wait() for p in [p1, p2, p3]]
    return exit_codes


def main(argv):
    pulse_path = ''
    cell_line_for_samtools = ''

    try:
        opts, args = getopt.getopt(argv, "hp:c:", ["path=", "cell_line="])
    except getopt.GetoptError:
        print('samtools.py -p <pulse_path> -c <cell_line_for_samtools>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('samtools.py -p <pulse_path> -c <cell_line_for_samtools>')
            sys.exit()
        elif opt in ("-p", "--path"):
            pulse_path = arg
        elif opt in ("-c", "--cell-line"):
            cell_line_for_samtools = arg
    print(pulse_path, cell_line_for_samtools)
    run_samtools(pulse_path, cell_line_for_samtools)

if __name__ == "__main__":
    main(sys.argv[1:])
