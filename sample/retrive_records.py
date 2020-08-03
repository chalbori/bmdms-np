import os
import sys
import inspect
import subprocess

CURRENT_DIR = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
PARENT_DIR = os.path.dirname(CURRENT_DIR)
sys.path.insert(0, PARENT_DIR)
print("sys.path: {}".format(sys.path))
from database.accounts import Database
from database import info_chem_mysql
from database import info_spectrum_mysql
from tools import tool_file


def main(output_file, limit=288939):
    """
    Basic example script to retrieve all spectral records of BMDMS-NP

    usage:

        (source code folder/sample):~$ python retrieve_records.py <output_file> <num_of_records>

    example:
        (source code folder/sample):~$ python retrieve_records.py "bmdms-np.msp" "5"(max = 288939)

    parameters:
        <instrument_type> : "Orbitrap", "O", "QTOF", "Q-TOF", "Q"
        <ion_mode> : "[M+H]+", "H", "[M+Na]+", "Na"
        <precursor_mz> : float (see compounds.csv)
    """
    db_bmdms_uri = "mysql+pymysql://{}:{}@{}/{}".format(
        Database.user_id, Database.user_password, Database.host, Database.database
    )
    
    repo_spectrum = info_spectrum_mysql.SpectrumRepository(db_uri=db_bmdms_uri)
    repo_compound = info_chem_mysql.CompoundRepository(db_uri=db_bmdms_uri)
    
    stream_spectrum = repo_spectrum.get_stream(include_null=False, limit=limit)
    f = open(file=output_file, mode='w', encoding='utf-8')

    for _, spectrum in enumerate(stream_spectrum):
        chemical = repo_compound.get(c_id=spectrum.snapeaks_id)
        peaklist = spectrum.peaks_cut_noise("dynamic")
        if len(peaklist[0]) > 0:
            # print(tool_file.convert_spec_to_msp(chemical=chemical, spec=spectrum, peaklist=peaklist))
            file_text = tool_file.convert_spec_to_msp(chemical=chemical, spec=spectrum, peaklist=peaklist)

            f.write(file_text)
            f.write('\n')
    f.close()
    repo_compound.close()
    repo_spectrum.close()

if __name__ == "__main__":
    if len(sys.argv) != 2 or len(sys.argv) != 3:
        print(main.__doc__)
        exit()
    output_file = sys.argv[1]
    main(output_file=output_file, limit=sys.argv[2])
