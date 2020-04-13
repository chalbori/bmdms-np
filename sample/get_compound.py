import os
import sys, inspect

CURRENT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, PARENT_DIR)
# print("sys.path: {}".format(sys.path))

from database import info_chem_mysql
from database.accounts import Database

def main(compound_id):
    """
    Basic example script to get information of compound

    usage:

        (source code folder/sample):~$ python get_compound.py <compound_id>

    example:
        (source code folder/sample):~$ python get_compound.py "c20150100089"

    """
    uri_bmdms = "mysql+pymysql://{}:{}@{}/{}".format(
        Database.user_id, Database.user_password, Database.host, Database.database
    )
    repo_compound = info_chem_mysql.CompoundRepository(db_uri=uri_bmdms)
    compound = repo_compound.get(c_id=compound_id)
    print(compound)
    repo_compound.close()

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print(main.__doc__)
        exit()
    print(len(sys.argv))
    COMPOUND_ID = sys.argv[1]

    main(compound_id=COMPOUND_ID)
