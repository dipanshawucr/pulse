# used to call iupred and score all the seq in a fasta file. the ouput is suitable for the features script
import os


def delete_tmp_file(feature_extract_output_path):
    os.system("rm " + feature_extract_output_path + "/temp.seq")
    return


def write_tmp_file(content, feature_extract_output_path):
    write_to = open(feature_extract_output_path + '/temp.seq', 'w')
    print >> write_to, "temp seq"
    print >> write_to, content
    return


def exec_iupred(iupred_path, feature_extract_output_path):
    return os.popen(iupred_path + " " + feature_extract_output_path + "/temp.seq long").readlines()


def save_iupred_out(output, prot_id, output_file):
    for line in output:
        if line[0] == "#":
            pass
        else:
            pos = line[0:5].strip()
            res = line[6:8].strip()
            score = line[13:19].strip()
            binary_score = transform_score(score)
            print >> output_file, "%s\t%s\t%s" % (prot_id, str(pos), binary_score)
    return


def transform_score(score):
    if float(score) > .5:
        return "."
    else:
        return "*"


def generate_iupred_file(filename, feature_extract_output_path,
                          iupred_install_path, output_location):
    prot_seq = ''
    prot_id = ''
    output_file = open(output_location, 'w')
    with open(filename, 'r') as input_file:
        for line in input_file:
            if line[0] == ">":
                if len(prot_seq) > 1:
                    write_tmp_file(prot_seq, feature_extract_output_path)
                    output_iupred = exec_iupred(iupred_install_path, feature_extract_output_path)
                    save_iupred_out(output_iupred, prot_id, output_file)
                prot_id = line[1:].strip()
                prot_seq = ''
            else:
                prot_seq += line

        write_tmp_file(prot_seq, feature_extract_output_path)
        output_iupred = exec_iupred(iupred_install_path, feature_extract_output_path)
        save_iupred_out(output_iupred, prot_id, output_file)
