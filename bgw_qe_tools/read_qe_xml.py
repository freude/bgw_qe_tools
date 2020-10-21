import numpy as np
import xml.etree.ElementTree as ET
from ase.io import Trajectory, cif
from ase.io import espresso
from ase.visualize import view
from ase import Atoms, Atom
from ase.calculators.singlepoint import (SinglePointDFTCalculator,
                                         SinglePointKPoint)
from ase.units import create_units
import matplotlib.pyplot as plt

# Quantum ESPRESSO uses CODATA 2006 internally
units = create_units('2006')
# atoms = cif.read_cif('/home/mk/tetracene/tetracene.cif', None)


def get_atoms_from_xml(xml_elem):
    atoms = []

    for item in xml_elem.find('atomic_structure').find('atomic_positions').findall('atom'):
        coords = np.array(item.text.split(), dtype=np.float) * units['Bohr']
        atom = Atom(symbol=item.items()[0][1],
                    position=coords)
        atoms.append(atom)

    a1 = np.array(xml_elem.find('atomic_structure').find('cell').find('a1').text.split(), dtype=np.float)
    a2 = np.array(xml_elem.find('atomic_structure').find('cell').find('a2').text.split(), dtype=np.float)
    a3 = np.array(xml_elem.find('atomic_structure').find('cell').find('a3').text.split(), dtype=np.float)

    cell = np.array([a1, a2, a3]) * units['Bohr']

    atoms = Atoms(atoms, cell=cell, pbc=True)

    return atoms


def _parse_xml(input_tag):
    atoms = get_atoms_from_xml(input_tag)

    parsed_data = {}

    for elem in input_tag:
        parsed_data.setdefault(elem.tag, {})
        for item in elem:
            parsed_data[elem.tag][item.tag] = item.text

    parsed_data['dft']['functional'] = input_tag.find('dft').find('functional').text

    parsed_data['atomic_species'] = {}
    for item in input_tag.find('atomic_species'):
        parsed_data['atomic_species'][item.items()[0][1]] = item.find('pseudo_file').text

    if input_tag.find('k_points_IBZ').find('monkhorst_pack') is not None:
        parsed_data['k_points'] = [int(item[1]) for item in input_tag.find('k_points_IBZ').find('monkhorst_pack').items()]

    if 'electric_field' in parsed_data:
        for item in input_tag.find('electric_field'):
            parsed_data['electric_field'][item.tag] = item.text

    return atoms, parsed_data


def parse_xml_input(xml_file):

    root = ET.parse(xml_file).getroot()
    input_tag = root.find('input')
    atoms, parsed_data = _parse_xml(input_tag)

    return atoms, parsed_data


def parse_xml_output(xml_file):

    root = ET.parse(xml_file).getroot()
    output_tag = root.find('output')
    atoms, parsed_data = _parse_xml(output_tag)

    return atoms, parsed_data


def xml2atoms(xml_elem):
    # symbols = [el.items()[0][1] for el in xml_elem.find('atomic_species').findall('species')]
    #
    # nat = int(xml_elem.find('atomic_structure').items()[0][1])
    # alat = float(xml_elem.find('atomic_structure').items()[1][1])

    atoms = []

    for item in xml_elem.find('atomic_structure').find('atomic_positions').findall('atom'):
        coords = np.array(item.text.split(), dtype=np.float) * units['Bohr']
        atom = Atom(symbol=item.items()[0][1],
                    position=coords)
        atoms.append(atom)

    a1 = np.array(xml_elem.find('atomic_structure').find('cell').find('a1').text.split(), dtype=np.float)
    a2 = np.array(xml_elem.find('atomic_structure').find('cell').find('a2').text.split(), dtype=np.float)
    a3 = np.array(xml_elem.find('atomic_structure').find('cell').find('a3').text.split(), dtype=np.float)

    cell = np.array([a1, a2, a3]) * units['Bohr']

    atoms = Atoms(atoms, cell=cell, pbc=True)

    # ------------------------------------------------------------------------------------------

    # Extract calculation results
    # Energy
    energy = 2 * float(xml_elem.find('total_energy').find('etot').text) * units['Ry']

    # Forces
    forces = []
    if xml_elem.find('forces') is not None:
        for item in xml_elem.find('forces').text.strip().split('\n'):
            forces.append(np.fromstring(item.strip(), sep=' ', dtype=np.float))

    forces = np.array(forces) * units['Ry'] / units['Bohr']

    # Stress
    stress = None

    # Magmoms
    magmoms = None

    # Fermi level / highest occupied level
    efermi = None
    if xml_elem.find('band_structure') is not None:
        efermi = 2 * float(xml_elem.find('band_structure').find('fermi_energy').text) * units['Ry']

    # # K-points
    ibzkpts = None
    # weights = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print them."
    #
    # for kpts_index in indexes[_PW_KPTS]:
    #     nkpts = int(pwo_lines[kpts_index].split()[4])
    #     kpts_index += 2
    #
    #     if pwo_lines[kpts_index].strip() == kpoints_warning:
    #         continue
    #
    #     # QE prints the k-points in units of 2*pi/alat
    #     # with alat defined as the length of the first
    #     # cell vector
    #     cell = structure.get_cell()
    #     alat = np.linalg.norm(cell[0])
    #     ibzkpts = []
    #     weights = []
    #     for i in range(nkpts):
    #         L = pwo_lines[kpts_index + i].split()
    #         weights.append(float(L[-1]))
    #         coord = np.array([L[-6], L[-5], L[-4].strip('),')],
    #                          dtype=float)
    #         coord *= 2 * np.pi / alat
    #         coord = kpoint_convert(cell, ckpts_kv=coord)
    #         ibzkpts.append(coord)
    #     ibzkpts = np.array(ibzkpts)
    #     weights = np.array(weights)
    #
    # # Bands
    # kpts = None
    # kpoints_warning = "Number of k-points >= 100: " + \
    #                   "set verbosity='high' to print the bands."
    #
    # for bands_index in indexes[_PW_BANDS] + indexes[_PW_BANDSTRUCTURE]:
    #     if image_index < bands_index < next_index:
    #         bands_index += 2
    #
    #         if pwo_lines[bands_index].strip() == kpoints_warning:
    #             continue
    #
    #         assert ibzkpts is not None
    #         spin, bands, eigenvalues = 0, [], [[], []]
    #
    #         while True:
    #             L = pwo_lines[bands_index].replace('-', ' -').split()
    #             if len(L) == 0:
    #                 if len(bands) > 0:
    #                     eigenvalues[spin].append(bands)
    #                     bands = []
    #             elif L == ['occupation', 'numbers']:
    #                 # Skip the lines with the occupation numbers
    #                 bands_index += len(eigenvalues[spin][0]) // 8 + 1
    #             elif L[0] == 'k' and L[1].startswith('='):
    #                 pass
    #             elif 'SPIN' in L:
    #                 if 'DOWN' in L:
    #                     spin += 1
    #             else:
    #                 try:
    #                     bands.extend(map(float, L))
    #                 except ValueError:
    #                     break
    #             bands_index += 1
    #
    #         if spin == 1:
    #             assert len(eigenvalues[0]) == len(eigenvalues[1])
    #         assert len(eigenvalues[0]) == len(ibzkpts), \
    #             (np.shape(eigenvalues), len(ibzkpts))
    #
    #         kpts = []
    #         for s in range(spin + 1):
    #             for w, k, e in zip(weights, ibzkpts, eigenvalues[s]):
    #                 kpt = SinglePointKPoint(w, s, k, eps_n=e)
    #                 kpts.append(kpt)

    # ------------------------------------------------------------------------------------------

    calc = SinglePointDFTCalculator(atoms, energy=energy,
                                    forces=forces, stress=stress,
                                    magmoms=magmoms, efermi=efermi,
                                    ibzkpts=ibzkpts)
    # calc.kpts = kpts
    atoms.calc = calc

    input_parameters = {}

    input_parameters['ecut'] = None
    if xml_elem.find('basis_set') is not None:
        input_parameters['ecut'] = 2 * float(xml_elem.find('basis_set').find('ecutwfc').text)

    input_parameters['input_dft'] = None
    if xml_elem.find('dft') is not None:
        input_parameters['input_dft'] = xml_elem.find('dft').find('functional').text.lower()

    input_parameters['k_points'] = None
    if xml_elem.find('band_structure') is not None:
        k_points = xml_elem.find('band_structure').find('starting_k_points').find('monkhorst_pack').items()
        k_points = [int(item[1]) for item in k_points]
        input_parameters['k_points'] = k_points

    input_parameters['fft_grid'] = None
    if xml_elem.find('basis_set') is not None:
        fft_points = xml_elem.find('basis_set').find('fft_grid').items()
        fft_points = [int(item[1]) for item in fft_points]
        input_parameters['fft_grid'] = fft_points

    return atoms, input_parameters


def traj_from_qe_xml(fileobj, index=-1, results_required=True):
    """Reads Quantum ESPRESSO output files.

    The atomistic configurations as well as results (energy, force, stress,
    magnetic moments) of the calculation are read for all configurations
    within the output file.

    Will probably raise errors for broken or incomplete files.

    Parameters
    ----------
    fileobj : file|str
        A file like object or filename
    index : slice
        The index of configurations to extract.
    results_required : bool
        If True, atomistic configurations that do not have any
        associated results will not be included. This prevents double
        printed configurations and incomplete calculations from being
        returned as the final configuration with no results data.

    Yields
    ------
    structure : Atoms
        The next structure from the index slice. The Atoms has a
        SinglePointCalculator attached with any results parsed from
        the file.


    """

    root = ET.parse(fileobj).getroot()
    output = root.find('output')
    steps = root.findall('step')

    atoms, input_parameters = xml2atoms(output)

    trajectory = None

    trajectory = Trajectory('t1.traj', 'a')
    atoms_list = []

    for step in steps:
        aaa, _ = xml2atoms(step)
        trajectory.write(aaa)
        atoms_list.append(aaa)

    trajectory.close()

    return atoms, input_parameters, atoms_list


def qe_xml_to_kgrid(file_name, kpoints=None, koffsets=None, q_shifts=None):
    atoms, input_parameters, _ = traj_from_qe_xml(file_name)
    kgrid = []

    # -------------------------------------------------------------------------
    # --------------------- collect data from xml file ------------------------
    # -------------------------------------------------------------------------

    elements = {}
    indices = []
    count = 0

    for item in atoms.get_chemical_symbols():
        if item in elements:
            indices.append(elements[item])
        else:
            count += 1
            elements[item] = count
            indices.append(elements[item])

    # cell in alat
    cell = atoms.get_cell() / atoms.get_cell().tolist()[0][0]

    # positions in alat
    positions = atoms.get_positions() / atoms.get_cell().tolist()[0][0]

    # fft grid parameters
    fft_grid = input_parameters['fft_grid']

    # k-grid parameters
    k_grid = input_parameters['k_points']

    # -------------------------------------------------------------------------
    # ------------------------ form the input files ---------------------------
    # -------------------------------------------------------------------------

    if k_points is not None:
        kgrid.append(' '.join(map(str, k_points)))
    else:
        kgrid.append(' '.join(map(str, k_grid[:3])))

    if k_offsets is not None:
        kgrid.append(' '.join(map(str, k_points)))
    else:
        kgrid.append(' '.join(map(str, 0.5 * np.array(k_grid[3:]))))

    if q_shifts is not None:
        kgrid.append(' '.join(map(str, q_shifts)))
    else:
        kgrid.append(' '.join(map(str, [0.0, 0.0, 0.0])))

    for item in cell:
        kgrid.append(' '.join(map(str, item)))

    kgrid.append(str(len(positions)))

    for j, item in enumerate(positions):
        kgrid.append(str(indices[j]) + ' ' + ' '.join(map(str, item)))

    kgrid.append(' '.join(map(str, fft_grid)))
    kgrid.append('.false.')
    kgrid.append('.false.')
    kgrid.append('.false.')

    return '\n'.join(kgrid)


def make_kgrids(file_name, kpoints=None, koffsets=None, q_shifts=None):

    kgrid = qe_xml_to_kgrid(file_name, k_points=k_points, k_offsets=k_offsets)
    qgrid = qe_xml_to_kgrid(file_name, k_points=k_points, k_offsets=k_offsets, q_shifts=q_shifts)

    with open("kgrid.in", "w") as text_file:
        text_file.write(kgrid)

    with open("kgrid_q.in", "w") as text_file:
        text_file.write(qgrid)


def make_pw2bgw(filename, outdir, prefix, real_or_complex=2):
    with open(filename, "r") as text_file:
        text = text_file.readlines()

    ans = []
    ans.append('&INPUT_PW2BGW')
    ans.append('prefix = \'' + prefix + '\'')
    if outdir is not None:
        ans.append('outdit = \'' + outdir + '\'')
    else:
        ans.append('outdit = \'tmp\'')

    ans.append('real_or_complex = ' + str(real_or_complex))
    ans.append('wfng_flag = ' + '.true.')
    ans.append('wfng_file = \'WFN\'')
    ans.append('wfng_kgrid = ' + '.true.')

    ans.append('wfng_nk1 = ' + text[0].split()[0])
    ans.append('wfng_nk2 = ' + text[0].split()[1])
    ans.append('wfng_nk3 = ' + text[0].split()[2])

    with open('k_grid.out', "r") as text_file:
        text = text_file.readlines()

    ans.append('wfng_dk1 = ' + text[1].split()[0])
    ans.append('wfng_dk2 = ' + text[1].split()[1])
    ans.append('wfng_dk3 = ' + text[1].split()[2])

    # -------------------------------------------------

    if abs(int(text[1].split()[0])) + abs(int(text[1].split()[0])) + abs(int(text[1].split()[0])) > 0.0:

        ans.append('rhog_flag = ' + '.true.')
        ans.append('rhog_file = \'' + 'RHO' + '\'')

        ans.append('vxcg_flag = ' + '.true.')
        ans.append('vxcg_file = \'' + 'VXC' + '\'')

        ans.append('vxc0_flag = ' + '.true.')
        ans.append('vxc0_file = \'' + 'vxc0.dat' + '\'')

        ans.append('vxc_flag = ' + '.true.')
        ans.append('vxc_file = \'' + 'vxc.dat' + '\'')

        ans.append('vxc_integral = \'g\'')
        ans.append('vxc_diag_nmin = ' + '180')
        ans.append('vxc_diag_nmax = ' + '220')

        ans.append('vxc_offdiag_nmin = ' + '0')
        ans.append('vxc_offdiag_nmax = ' + '0')

        ans.append('vxc_zero_rho_core = ' + '.true.')

    ans.append('/')
    aaa = '\n'.join(ans)

    return aaa


def dict_retyping(dic):

    for key, value in dic.items():
        if value == 'true':
            dic[key] = True
        elif value == 'false':
            dic[key] = False
        else:
            try:
                dic[key] = int(value)
            except (ValueError, TypeError):
                try:
                    dic[key] = float(value)
                except (ValueError, TypeError):
                    pass

    return dic


def _parsed_to_input_data(parsed_data):

    input_data = {'control': {}, 'system': {}, 'electrons': {}}

    # input_data = {'control': {'calculation': 'bands',
    #                           'restart_mode': 'from_scratch',
    #                           'tstress': True,
    #                           'tprnfor': True,
    #                           'outdir': 'tmp',
    #                           'prefix': 'slab_si_tc',
    #                           'forc_conv_thr': 1e-5,
    #                           'nstep': 100,
    #                           'tefield': True,
    #                           'dipfield': True
    #                           },
    #               'system': {'ecutwfc': 40.0,
    #                          'input_dft': 'vdw-df2-c09',
    #                          'edir': 3,
    #                          'emaxpos': 0.95,
    #                          'eopreg': 0.03,
    #                          'eamp': 0
    #                          },
    #               'electrons': {'mixing_beta': 0.1,
    #                             'conv_thr': 1e-8},
    #               'ions': {'ion_dynamics': 'bfgs'}}

    input_data['control'] = dict_retyping(parsed_data['control_variables'])

    input_data['control']['calculation'] = 'bands'
    input_data['control']['restart_mode'] = 'from_scratch'

    input_data['system']['ecutwfc'] = float(parsed_data['basis']['ecutwfc'])
    input_data['system']['input_dft'] = parsed_data['dft']['functional']

    if 'electric_field' in parsed_data:
        if parsed_data['electric_field']['dipole_correction']:
            input_data['control']['tefield'] = True
            input_data['control']['dipfield'] = True
            input_data['system']['edir'] = int(parsed_data['electric_field']['electric_field_direction'])
            input_data['system']['emaxpos'] = float(parsed_data['electric_field']['potential_max_position'])
            input_data['system']['eopreg'] = float(parsed_data['electric_field']['potential_decrease_width'])
            input_data['system']['eamp'] = float(parsed_data['electric_field']['electric_field_amplitude'])

    input_data['electrons'] = dict_retyping(parsed_data['electron_control'])

    return input_data


def make_input_bands(xml_file, kpoints=None, koffsets=None):

    atoms_in, data_in = parse_xml_input(xml_file)
    atoms_out, data_out = parse_xml_input(xml_file)
    output_file = 'espresso_40_5_4_1.pwi'

    input_data = _parsed_to_input_data(data_in)

    if kpoints is None:
        kpoints = tuple(data_in['k_points'][:3])

    if koffsets is None:
        koffsets = tuple(data_in['k_points'][3:])

    with open(output_file, 'w') as fd:
        espresso.write_espresso_in(fd, atoms_in,
                                   input_data=input_data,
                                   pseudopotentials=data_in['atomic_species'],
                                   kspacing=None,
                                   koffset=koffsets,
                                   kpts=kpoints)

    outdir = data_in['control_variables']['outdir']
    prefix = data_in['control_variables']['prefix']

    return outdir, prefix


if __name__ == '__main__':
    # file_name = '/home/mk/data-file-schema.xml'
    # file_name = '/home/mk/tetracene_opt2.xml'
    xml_file_name = '/home/mk//si_slab.xml'

    # make_kgrids(xml_file_name, k_points=[5, 4, 1], q_shifts=[0.001, 0.0, 0.0])
    outdir, prefix = make_input_bands(xml_file_name)
    make_pw2bgw(outdir, prefix)


    # kgrids_to_pw2bgw('kgrid.in', 'si_slab_70_5_4_1')
    # # file_name1 = '/home/mk/tetracene.xml'
    # # file_name2 = '/home/mk/tetracene1.xml'
    #
    # # atoms, atoms_list1 = read_qe_xml(file_name1)
    # # atoms, atoms_list2 = read_qe_xml(file_name2)
    # # atoms_list = atoms_list1+atoms_list2
    # atoms, ecut, atoms_list = traj_from_qe_xml(file_name)
    # atoms1, ecut, atoms_list1 = traj_from_qe_xml(file_name1)
    #
    # # print(len(atoms_list))
    # view(atoms_list1)
    #
    # etot = [item.get_total_energy() for item in atoms_list]
    # etot1 = [item.get_total_energy() for item in atoms_list1]
    #
    # etot = [np.linalg.norm(item.get_cell()[1]) for item in atoms_list]
    # etot1 = [np.linalg.norm(item.get_cell()[1]) for item in atoms_list1]
    #
    # plt.plot(etot)
    # plt.plot(etot1)
    # plt.show()
