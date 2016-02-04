import unicodedata


def normalize_unicode_data(data):
    """
    Returns a Python-normalized String object from unicode.

    :param data:
    :return:
    """
    normalized_data = unicodedata.normalize('NFKD', data).encode('ascii', 'ignore')
    return normalized_data
