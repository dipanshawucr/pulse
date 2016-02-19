import os


def create_paths_for_cell_line(cell_line):
    """
    Creates all necessary paths for for_preprocess:

    Base folder for cell line:
    - ../output/for_preprocess/<CELL_LINE_NAME>/

    samtools output:
    - ../output/for_preprocess/<CELL_LINE_NAME>/no-rg/

    cufflinks_output:
    - ../output/for_preprocess/<CELL_LINE_NAME>/cufflinks_output/

    :param cell_line:
    :return:
    """
    paths_to_create = [r'output/for_preprocess/' + cell_line,
                       r'output/for_preprocess/' + cell_line + r'/no-rg',
                       r'output/for_preprocess/' + cell_line + r'/cufflinks_output']
    for path in paths_to_create:
        if not os.path.exists(path):
            os.makedirs(path)
