import os


def create_paths_for_cell_line(pulse_path, cell_line):
    """
    Creates all necessary paths for features:

    Base folder for cell line:
    - ../output/features/<CELL_LINE_NAME>/

    :param cell_line:
    :param pulse_path:
    :return:
    """

    paths_to_create = [pulse_path + r'/output/features/' + cell_line]
    for path in paths_to_create:
        print path
        if not os.path.exists(path):
            os.makedirs(path)
