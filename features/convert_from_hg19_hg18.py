import subprocess


def use_remap_api(remap_api_path, mode, remap_from, remap_to, remap_input, remap_output):
    """

    :param remap_api_path:
    :param mode:
    :param remap_from:
    :param remap_to:
    :param remap_input:
    :param remap_output:
    :return:
    """
    command1 = ['perl',
                remap_api_path,
                '--mode',
                mode,
                '--from',
                remap_from,
                '--dest',
                remap_to,
                '--annotation',
                remap_input,
                '--report_out',
                remap_output]
    p1 = subprocess.Popen(command1)
    exit_code = p1.wait()
    return exit_code
