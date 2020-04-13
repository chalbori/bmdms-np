import os
import sys, inspect

CURRENT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, PARENT_DIR)
# print("sys.path: {}".format(sys.path))

from database import info_spectrum_mysql
from database.accounts import Database

def main(spec_id):
    """
    Basic example script to search query block against BMDMS-NP

    usage:

        (source code folder/sample):~$ python get_spectrum.py <spec_id>

    example:
        (source code folder/sample):~$ python get_spectrum.py 3457

    """
    uri_bmdms = "mysql+pymysql://{}:{}@{}/{}".format(
        Database.user_id, Database.user_password, Database.host, Database.database
    )
    repo_spectrum = info_spectrum_mysql.SpectrumRepository(db_uri=uri_bmdms)
    spectrum = repo_spectrum.get(spec_id=spec_id)
    print(spectrum)
    repo_spectrum.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(main.__doc__)
        exit()
    print(len(sys.argv))
    SPECTRUM_ID = sys.argv[1]

    main(spec_id=SPECTRUM_ID)
