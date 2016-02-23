import sys
import getopt
import subprocess


def run_cufflinks(pulse_path, cell_line_name):
    output_file = open(pulse_path + '/output/for_preprocess/' + cell_line_name + '/transcripts.gtf', 'w')
    command1 = ['cufflinks',
                '-o',
                pulse_path + '/output/for_preprocess/' + cell_line_name + '/cufflinks_output',
                '-g',
                pulse_path + '/input/Homo_sapiens.GRCh37.70.gtf',
                pulse_path + '/output/for_preprocess/' + cell_line_name + '/no-rg/' + cell_line_name]
    p1 = subprocess.Popen(command1, stdout=output_file)
    exit_codes = p1.wait()
    return exit_codes


def main(argv):
    pulse_path = ''
    cell_line_for_cufflinks = ''

    try:
        opts, args = getopt.getopt(argv, "hp:c:", ["path=", "cell_line="])
    except getopt.GetoptError:
        print('cufflinks.py -p <pulse_path> -c <cell_line_for_cufflinks>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('cufflinks.py -p <pulse_path> -c <cell_line_for_cufflinks>')
            sys.exit()
        elif opt in ("-p", "--path"):
            pulse_path = arg
        elif opt in ("-c", "--cell-line"):
            cell_line_for_cufflinks = arg
    run_cufflinks(pulse_path, cell_line_for_cufflinks)

if __name__ == "__main__":
    main(sys.argv[1:])
