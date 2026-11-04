from ase.spacegroup import crystal
from collections import Counter
import numpy as np
from ase.build import make_supercell
import math
from ase.neighborlist import neighbor_list
from itertools import combinations
a = 3.076
c = 25.3558

spacegroup = 164

DIRECTION = "0"
LENGTH_NM = 100.0
WIDTH_NM = 5.0
VACUUM_A = 20.0
BOND_CUTOFF = 2.45

basis = ["Ti", "Ti", "C", "F"]

positions = [
    (0.0, 0.0, 0.5000),      # Ti1 - wewnętrzny
    (1/3, 2/3, 0.5928),      # Ti2 - powierzchniowy
    (1/3, 2/3, 0.4492),      # C
    (0.0, 0.0, 0.3578)       # F
]

unit_cell = crystal(
    basis,
    positions,
    spacegroup=spacegroup,
    cellpar=[a, a, c, 90, 90, 120]
)
counts = Counter(unit_cell.get_chemical_symbols())

print(unit_cell)
print(counts)

assert counts == Counter({
    "Ti": 3,
    "C": 2,
    "F": 2
})
typing_cell = unit_cell.repeat(
    (2, 2, 1)
)
assert len(typing_cell) == 28
scaled = typing_cell.get_scaled_positions(
    wrap=True
)
def parity(value, base):

    delta = (value - base) % 1.0

    if (
        np.isclose(delta, 0.0, atol=1e-5)
        or np.isclose(delta, 1.0, atol=1e-5)
    ):
        return 0

    if np.isclose(delta, 0.5, atol=1e-5):
        return 1

    raise ValueError(
        f"Nie rozpoznano pozycji: "
        f"value={value}, base={base}"
    )
lower_C_map = {
    (0, 0): "cx3",
    (1, 0): "cx1",
    (0, 1): "cx4",
    (1, 1): "cx2"
}
upper_C_map = {
    (0, 0): "cx2",
    (1, 0): "cx4",
    (0, 1): "cx1",
    (1, 1): "cx3"
}
iff_types = []

for atom, frac in zip(
    typing_cell,
    scaled):
    x, y, z = frac
    if atom.symbol == "Ti":
        if np.isclose(
            z,
            0.5,
            atol=1e-5
        ):
            iff_types.append("ti1")
        else:
            iff_types.append("ti2")

    elif atom.symbol == "F":

        iff_types.append("fx")

    elif atom.symbol == "C":
        if z < 0.5:
            i = parity(x, 1/6)
            j = parity(y, 1/3)

            iff_types.append(
                lower_C_map[(i, j)]
            )
        else:
            i = parity(x, 1/3)
            j = parity(y, 1/6)

            iff_types.append(
                upper_C_map[(i, j)]
            )
    else:
        raise ValueError(
            f"Nieobsługiwany pierwiastek: "
            f"{atom.symbol}"
        )
typing_cell.new_array(
    "iff_type",
    np.array(
        iff_types,
        dtype="U3"
    )
)
type_counts = Counter(
    iff_types
)
print("IFF types:")
print(type_counts)

expected_types = Counter({
    "ti1": 4,
    "ti2": 8,

    "cx1": 2,
    "cx2": 2,
    "cx3": 2,
    "cx4": 2,

    "fx": 8
})
assert type_counts == expected_types
IFF_CHARGE = {
    "ti1": +0.350,
    "ti2": +0.675,

    "cx1": -0.350,
    "cx2": -0.350,
    "cx3": -0.350,
    "cx4": -0.350,

    "fx": -0.500
}
IFF_TYPE_ID = {
    "ti2": 1,
    "cx3": 2,
    "cx1": 3,
    "cx4": 4,
    "cx2": 5,
    "ti1": 6,
    "fx": 7
}
IFF_MASS = {
    1: (47.876000, "ti2"),
    2: (12.011150, "cx3"),
    3: (12.011150, "cx1"),
    4: (12.011150, "cx4"),
    5: (12.011150, "cx2"),
    6: (47.876000, "ti1"),
    7: (18.998400, "fx")
}
PAIR_COEFF = {
    1: (0.300000, 1.915444),
    2: (0.070000, 3.385417),
    3: (0.070000, 3.385417),
    4: (0.070000, 3.385417),
    5: (0.070000, 3.385417),
    6: (0.200000, 2.004510),
    7: (0.170000, 2.904326)
}
TI_C_K = 285.0
TI_C_R0 = 2.157
ANGLE_COEFF = {
    1:  (40.0,  90.0),
    2:  (40.0,  90.0),
    3:  (40.0,  90.0),
    4:  (40.0,  90.0),
    5:  (40.0,  90.0),
    6:  (40.0,  90.0),

    7:  (67.0,  93.5),
    8:  (0.0,   90.0),
    9:  (67.0,  91.6),

    10: (67.0,  93.5),
    11: (0.0,   90.0),
    12: (67.0,  91.6),

    13: (0.0,   90.0),
    14: (67.0,  93.5),
    15: (67.0,  91.6),

    16: (67.0,  93.5),
    17: (0.0,   90.0),
    18: (67.0,  91.6),

    19: (40.0,  90.0),
    20: (40.0,  90.0),
    21: (67.0, 180.0),

    22: (40.0,  90.0),
    23: (67.0, 180.0),
    24: (67.0, 180.0),

    25: (40.0,  90.0),
    26: (40.0,  90.0),
    27: (67.0, 180.0),
    28: (40.0,  90.0)
}

charges = np.array([
    IFF_CHARGE[t]
    for t in iff_types
])

total_charge = charges.sum()
print(
    f"Total charge = "
    f"{total_charge:.8f} e")
assert np.isclose(
    total_charge,
    0.0,
    atol=1e-8)
typing_cell.set_initial_charges(
    charges)
print("Etap typowania IFF zakończony poprawnie.")
P_MATRICES = {
    "0": [
        [1, 0, 0],
        [1, 2, 0],
        [0, 0, 1]
    ],

    "19.11": [
        [3, 1, 0],
        [1, 5, 0],
        [0, 0, 1]
    ],

    "30": [
        [2, 1, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]
}
P = P_MATRICES[DIRECTION]
h = P[0][0]
k = P[0][1]
angle = math.degrees(
    math.atan2(
        math.sqrt(3) * k,
        2 * h - k
    )
)
print(
    f"Transport angle = "
    f"{angle:.4f} deg"
)
seed = make_supercell(
    typing_cell,
    P
)
# Eksport ponizej zapisuje pudlo ortogonalne w osiach x, y, z.
# Obracamy zarowno atomy, jak i wektory komorki.
seed.rotate(-angle, "z", rotate_cell=True)
seed.wrap()
if not np.allclose(seed.cell.array, np.diag(seed.cell.lengths()), atol=1e-8):
    raise RuntimeError("Komorka nie jest ustawiona wzdluz osi kartezjanskich.")
seed.set_cell(np.diag(seed.cell.lengths()), scale_atoms=False)

cell_parameters = seed.cell.cellpar()
print(
    "Seed cell:",
    cell_parameters
)
gamma = cell_parameters[5]

assert np.isclose(
    gamma,
    90.0,
    atol=1e-6
)
seed_lx, seed_ly = seed.cell.lengths()[:2]
nx = max(
    1,
    int(round(
        LENGTH_NM * 10 / seed_lx
    ))
)
ny = max(
    1,
    int(round(
        WIDTH_NM * 10 / seed_ly
    ))
)
nanosheet = seed.repeat(
    (nx, ny, 1)
)
nanosheet.center(
    vacuum=VACUUM_A,
    axis=2
)
nanosheet.pbc = (True, True, False)
if "iff_type" not in nanosheet.arrays:

    raise RuntimeError(
        "ASE zgubiło typy IFF."
    )
types = nanosheet.get_array(
    "iff_type"
)
final_type_counts = Counter(
    types
)
n_formula = len(nanosheet) // 7
expected_final_types = Counter({
    "ti1": n_formula,
    "ti2": 2 * n_formula,

    "cx1": n_formula // 2,
    "cx2": n_formula // 2,
    "cx3": n_formula // 2,
    "cx4": n_formula // 2,

    "fx": 2 * n_formula
})
assert (
    final_type_counts
    == expected_final_types
)
final_charges = (
    nanosheet.get_initial_charges()
)
assert np.isclose(
    final_charges.sum(),
    0.0,
    atol=1e-8
)
symbols = nanosheet.get_chemical_symbols()
i_list, j_list, d_list = neighbor_list(
    "ijd",
    nanosheet,
    cutoff=BOND_CUTOFF
)
bond_dict = {}
for i, j, distance in zip(
    i_list,
    j_list,
    d_list
):
    pair = {
        symbols[i],
        symbols[j]
    }
    if pair == {"Ti", "C"}:
        atom1, atom2 = sorted(
            (int(i), int(j))
        )
        bond_dict[
            (atom1, atom2)
        ] = float(distance)
bond_pairs = sorted(
    bond_dict.keys()
)
bond_lengths = np.array(
    list(bond_dict.values())
)
expected_bonds = 12 * n_formula
print(
    f"Ti-C bonds = "
    f"{len(bond_pairs)}"
)
print(
    f"Expected bonds = "
    f"{expected_bonds}"
)
assert len(bond_pairs) == expected_bonds
print(
    f"Ti-C bond length: "
    f"{bond_lengths.min():.4f} - "
    f"{bond_lengths.max():.4f} A"
)
BOND_TYPE = {
    ("ti2", "cx3"): 1,
    ("ti2", "cx1"): 2,
    ("ti2", "cx4"): 3,
    ("ti2", "cx2"): 4,

    ("ti1", "cx3"): 5,
    ("ti1", "cx1"): 6,
    ("ti1", "cx4"): 7,
    ("ti1", "cx2"): 8
}
def get_bond_type(i, j):
    type1 = str(types[i])
    type2 = str(types[j])
    if type1.startswith("ti"):
        key = (type1, type2)
    elif type2.startswith("ti"):
        key = (type2, type1)
    else:
        raise ValueError(
            f"Nieznane wiązanie: "
            f"{type1}-{type2}"
        )
    return BOND_TYPE[key]
bonds = []
for bond_id, (i, j) in enumerate(
    bond_pairs,
    start=1
):

    bond_type = get_bond_type(
        i,
        j
    )

    bonds.append(
        (
            bond_id,
            bond_type,
            i + 1,
            j + 1
        )
    )
bond_type_counts = Counter(
    bond[1]
    for bond in bonds
)
print(
    f"Bond types = "
    f"{dict(bond_type_counts)}"
)
adjacency = [
    []
    for _ in range(len(nanosheet))
]
for i, j in bond_pairs:
    adjacency[i].append(j)
    adjacency[j].append(i)
expected_coordination = {
    "ti1": 6,
    "ti2": 3,

    "cx1": 6,
    "cx2": 6,
    "cx3": 6,
    "cx4": 6,

    "fx": 0
}
for i, neighbors in enumerate(adjacency):

    atom_type = str(types[i])

    actual = len(neighbors)

    expected = expected_coordination[
        atom_type
    ]

    if actual != expected:

        raise RuntimeError(
            f"Atom {i+1} ({atom_type}) "
            f"ma {actual} Bonds, "
            f"a powinien mieć {expected}"
        )
print(
    "Coordination test passed."
)
angle_count = 0
for neighbors in adjacency:

    n = len(neighbors)

    angle_count += (
        n * (n - 1) // 2
    )
expected_angles = (
    51 * n_formula
)
print(
    f"Possible angles = "
    f"{angle_count}"
)
print(
    f"Expected angles = "
    f"{expected_angles}"
)
assert (
    angle_count
    == expected_angles
)
TI2_ANGLE = {
    frozenset(("cx3", "cx1")): 1,
    frozenset(("cx3", "cx4")): 2,
    frozenset(("cx1", "cx4")): 3,
    frozenset(("cx1", "cx2")): 4,
    frozenset(("cx3", "cx2")): 5,
    frozenset(("cx4", "cx2")): 6
}
TI1_DIFFERENT_C = {
    frozenset(("cx3", "cx2")): 19,
    frozenset(("cx3", "cx4")): 20,
    frozenset(("cx4", "cx2")): 22,
    frozenset(("cx1", "cx2")): 25,
    frozenset(("cx1", "cx4")): 26,
    frozenset(("cx3", "cx1")): 28
}
TI1_SAME_C = {
    "cx3": 21,
    "cx2": 23,
    "cx4": 24,
    "cx1": 27
}
C_CENTER = {
    "cx3": {
        "ti2-ti2": 7,
        "mixed": 8,
        "ti1-ti1": 9
    },

    "cx1": {
        "ti2-ti2": 10,
        "mixed": 11,
        "ti1-ti1": 12
    },

    "cx4": {
        "mixed": 13,
        "ti2-ti2": 14,
        "ti1-ti1": 15
    },

    "cx2": {
        "ti2-ti2": 16,
        "mixed": 17,
        "ti1-ti1": 18
    }
}
def get_angle_type(i, center, k):

    left = str(types[i])
    middle = str(types[center])
    right = str(types[k])
    if middle == "ti2":

        key = frozenset(
            (left, right)
        )

        return TI2_ANGLE[key]
    if middle == "ti1":

        if left == right:

            return TI1_SAME_C[left]

        key = frozenset(
            (left, right)
        )

        return TI1_DIFFERENT_C[key]
    if middle.startswith("cx"):

        if left == "ti2" and right == "ti2":
            group = "ti2-ti2"

        elif left == "ti1" and right == "ti1":
            group = "ti1-ti1"

        else:
            group = "mixed"

        return C_CENTER[
            middle
        ][
            group
        ]
    raise ValueError(
        f"Nieznany kąt: "
        f"{left}-{middle}-{right}"
    )
angles = []
angle_id = 1
for center in range(len(nanosheet)):

    neighbors = sorted(
        adjacency[center]
    )
    for i, k in combinations(
        neighbors,
        2
    ):

        angle_type = get_angle_type(
            i,
            center,
            k
        )

        angles.append(
            (
                angle_id,
                angle_type,
                i + 1,
                center + 1,
                k + 1
            )
        )

        angle_id += 1
print(
    f"Generated angles = "
    f"{len(angles)}"
)
assert len(angles) == expected_angles
angle_type_counts = Counter(
    angle[1]
    for angle in angles
)
print(
    f"Angle types = "
    f"{dict(sorted(angle_type_counts.items()))}"
)

if not np.allclose(nanosheet.cell.array,
                   np.diag(nanosheet.cell.lengths()), atol=1e-8):
    raise RuntimeError("Writer wymaga pudla ortogonalnego zgodnego z osiami.")
Lx, Ly, Lz = nanosheet.cell.lengths()
nanosheet.wrap()
positions_final = nanosheet.get_positions()
direction_tag = DIRECTION.replace(".", "p")
output = (
    f"Ti3C2F2_IFF_"
    f"{direction_tag}deg_"
    f"{LENGTH_NM:.0f}x{WIDTH_NM:.0f}nm.data"
)
with open(
    output,
    "w",
    encoding="utf-8"
) as f:

    f.write(
        "Ti3C2F2 - IFF CVFF v1.6\n\n"
    )

    f.write(
        f"{len(nanosheet)} atoms\n"
    )

    f.write(
        f"{len(bonds)} bonds\n"
    )

    f.write(
        f"{len(angles)} angles\n\n"
    )

    f.write("7 atom types\n")
    f.write("8 bond types\n")
    f.write("28 angle types\n\n")

    f.write(
        f"0.0 {Lx:.10f} xlo xhi\n"
    )

    f.write(
        f"0.0 {Ly:.10f} ylo yhi\n"
    )

    f.write(
        f"0.0 {Lz:.10f} zlo zhi\n\n"
    )

    f.write("Masses\n\n")

    for type_id in range(1, 8):

        mass, name = IFF_MASS[
            type_id
        ]

        f.write(
            f"{type_id:4d} "
            f"{mass:12.6f} "
            f"# {name}\n"
        )

    f.write(
        "\nPair Coeffs # lj/cut/coul/long\n\n"
    )

    for type_id in range(1, 8):

        epsilon, sigma = PAIR_COEFF[
            type_id
        ]

        f.write(
            f"{type_id:4d} "
            f"{epsilon:12.6f} "
            f"{sigma:12.6f}\n"
        )

    f.write(
        "\nBond Coeffs # harmonic\n\n"
    )

    for bond_type in range(1, 9):

        f.write(
            f"{bond_type:4d} "
            f"{TI_C_K:12.6f} "
            f"{TI_C_R0:12.6f}\n"
        )

    f.write(
        "\nAngle Coeffs # harmonic\n\n"
    )

    for angle_type in range(1, 29):

        K, theta0 = ANGLE_COEFF[
            angle_type
        ]

        f.write(
            f"{angle_type:4d} "
            f"{K:12.6f} "
            f"{theta0:12.6f}\n"
        )

    f.write(
        "\nAtoms # full\n\n"
    )

    for i in range(len(nanosheet)):

        name = str(types[i])

        type_id = IFF_TYPE_ID[
            name
        ]

        charge = IFF_CHARGE[
            name
        ]

        x, y, z = positions_final[i]

        molecule_id = 1

        f.write(
            f"{i+1:8d} "
            f"{molecule_id:6d} "
            f"{type_id:4d} "
            f"{charge:10.6f} "
            f"{x:16.8f} "
            f"{y:16.8f} "
            f"{z:16.8f} "
            f"# {name}\n"
        )

    f.write(
        "\nBonds\n\n"
    )

    for (
        bond_id,
        bond_type,
        atom1,
        atom2
    ) in bonds:

        f.write(
            f"{bond_id:8d} "
            f"{bond_type:4d} "
            f"{atom1:8d} "
            f"{atom2:8d}\n"
        )

    f.write(
        "\nAngles\n\n"
    )

    for (
        angle_id,
        angle_type,
        atom1,
        center,
        atom3
    ) in angles:

        f.write(
            f"{angle_id:8d} "
            f"{angle_type:4d} "
            f"{atom1:8d} "
            f"{center:8d} "
            f"{atom3:8d}\n"
        )

print(
    f"Saved IFF file: {output}"
)

print("=" * 50)
print(
    f"Direction = {DIRECTION}")
print(
    f"Angle = {angle:.4f} deg")
print(
    f"Repetitions: "
    f"nx={nx}, ny={ny}")
print(
    f"Size = "
    f"{Lx/10:.4f} x "
    f"{Ly/10:.4f} nm")
print(
    f"Atoms = {len(nanosheet)}")
print(
    f"Formula units = {n_formula}")
print(
    f"IFF types = "
    f"{dict(final_type_counts)}")
print(
    f"Charge = "
    f"{final_charges.sum():.8f} e")
print("=" * 50)
print(f"Grubosc obwiedni atomowej = {np.ptp(positions_final[:, 2]):.6f} A")
print(f"Lz pudla z proznia = {Lz:.6f} A")
print("Struktura poczatkowa: wymaga minimalizacji i rownowazenia w LAMMPS.")
