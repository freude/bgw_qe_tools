from ase.io import cif
from ase.io import espresso


pseudopotentials = {'C': 'C.UPF', 'H': 'H.UPF'}


def cif2pwi(in_path, out_path):

    atoms = cif.read_cif(in_path, None)

    with open(out_path, 'w') as fd:
        espresso.write_espresso_in(fd, atoms,
                                   input_data=None,
                                   pseudopotentials=pseudopotentials,
                                   kspacing=(0, 0, 0),
                                   kpts=(7, 7, 5))




if __name__ == '__main__':

    in_path = '/home/mk/tetracene/tetracene.cif'
    out_path = 'tetracene.pwi'

    cif2pwi(in_path, out_path)
