import os
import sys
import sqlite3


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


def score_differences(mapping, uniprot, start, end):
    count = 0
    try:
        if uniprot in mapping:
            if start <= end:
                for i in range(start, end + 1):
                    if i not in mapping[uniprot]:
                        pass
                    elif mapping[uniprot][i] == "*":
                        count += 1
                return (1.0 * count) / (end - start + 1)

    except KeyError, e:
        print 'WARNING!: Features not found '
        print str(e)
        print 'Please check ID tags or/and generate the missing features'
    return count


def score_differences_pfam(mapping, uniprot, start, end):
    try:
        count = 0
        enzymatic = 0
        if uniprot in mapping:
            if start <= end:
                for i in range(start, end + 1):
                    if i not in mapping[uniprot]:
                        pass
                    elif mapping[uniprot][i][0] == "*":
                        if mapping[uniprot][i][1]:
                            enzymatic = 1
                        count += 1
                return [(1.0 * count) / (end - start + 1), enzymatic]
    except KeyError, e:
        print 'WARNING!: Features not found for'
        print  str(e)
        print 'Please check ID tags or/and generate the missing features'
    return [count, 0]


# Used in features.uniprot_disorder
def link_db(db_path):
    try:
        db = sqlite3.connect(db_path)
    except sqlite3.Error, errmsg:
        print 'DB not available ' + str(errmsg)
        sys.exit()
    db_cursor = db.cursor()
    return db_cursor, db
