//! Element data for H to Ne.

pub const BOHR_PER_ANGSTROM: f64 = 1.0 / 0.529_177_210_92;
pub const EV_PER_HARTREE: f64 = 27.211_386_245_988;
pub const DEBYE_PER_AU: f64 = 2.541_746_473;
/// Atomic mass unit in electron masses.
pub const AMU_IN_ME: f64 = 1_822.888_486;
/// Atomic unit of time in femtoseconds.
pub const FS_PER_AU_TIME: f64 = 0.024_188_843_265_857;

pub struct Element {
    pub symbol: &'static str,
    pub z: u32,
    /// Mass in amu.
    pub mass: f64,
    /// Covalent radius in Angstrom.
    pub covalent_radius: f64,
}

pub const ELEMENTS: [Element; 10] = [
    Element {
        symbol: "H",
        z: 1,
        mass: 1.008,
        covalent_radius: 0.31,
    },
    Element {
        symbol: "He",
        z: 2,
        mass: 4.0026,
        covalent_radius: 0.28,
    },
    Element {
        symbol: "Li",
        z: 3,
        mass: 6.94,
        covalent_radius: 1.28,
    },
    Element {
        symbol: "Be",
        z: 4,
        mass: 9.0122,
        covalent_radius: 0.96,
    },
    Element {
        symbol: "B",
        z: 5,
        mass: 10.81,
        covalent_radius: 0.84,
    },
    Element {
        symbol: "C",
        z: 6,
        mass: 12.011,
        covalent_radius: 0.76,
    },
    Element {
        symbol: "N",
        z: 7,
        mass: 14.007,
        covalent_radius: 0.71,
    },
    Element {
        symbol: "O",
        z: 8,
        mass: 15.999,
        covalent_radius: 0.66,
    },
    Element {
        symbol: "F",
        z: 9,
        mass: 18.998,
        covalent_radius: 0.57,
    },
    Element {
        symbol: "Ne",
        z: 10,
        mass: 20.18,
        covalent_radius: 0.58,
    },
];

pub fn element(symbol: &str) -> Option<&'static Element> {
    ELEMENTS
        .iter()
        .find(|e| e.symbol.eq_ignore_ascii_case(symbol))
}

pub fn element_by_z(z: u32) -> Option<&'static Element> {
    ELEMENTS.get((z as usize).wrapping_sub(1))
}
