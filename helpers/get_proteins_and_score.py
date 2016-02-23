import os

try:
    PULSE_PATH = os.environ['PULSE_PATH']
except KeyError:
    print('Environment variable PULSE_PATH not setup')
    print('Using working directory')
    PULSE_PATH = os.getcwd()


def get_proteins(pulse_names_output_file):
    """
    Returns a list of all proteins that appear in order from
    original PULSE output.

    :param pulse_names_output_file:
    :return:
    """
    pulse_output = open(pulse_names_output_file, "r")
    protein_list = []
    for line in pulse_output:
        protein = line.split("//")[-1].strip()
        protein_list.append(protein)
    return protein_list


def get_scores(pulse_scores_output_file):
    """
    Returns a list of all proteins' scores that appear
    in the order from pulse output.

    :param pulse_scores_output_file:
    :return:
    """
    pulse_scores = open(pulse_scores_output_file, "r")
    score_list = []
    for line in pulse_scores:
        score = line.strip()
        score_list.append(score)
    return score_list

if __name__ == "__main__":
    proteins_list = get_proteins(PULSE_PATH+'/output/machine/names.txt')
    scores_list = get_scores(PULSE_PATH+'/output/machine/PULSE_Output.txt')

    # Combine two lists
    protein_to_score = zip(proteins_list, scores_list)

    for pair in protein_to_score:
        print(pair[0])

