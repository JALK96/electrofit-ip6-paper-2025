from __future__ import annotations

# GROMACS labels as they appear in hb.log for the default IP6 system.
GROMACS_OXYGEN_TO_GRO_PHOSPHATE = {
    "O": "P",
    "O6": "P",
    "O7": "P",
    "O8": "P",
    "O1": "P1",
    "O9": "P1",
    "O10": "P1",
    "O11": "P1",
    "O2": "P2",
    "O12": "P2",
    "O13": "P2",
    "O14": "P2",
    "O3": "P3",
    "O15": "P3",
    "O16": "P3",
    "O17": "P3",
    "O4": "P4",
    "O18": "P4",
    "O19": "P4",
    "O20": "P4",
    "O5": "P5",
    "O21": "P5",
    "O22": "P5",
    "O23": "P5",
}

GRO_TO_PAPER_PHOSPHATE = {
    "P": "P1",
    "P1": "P2",
    "P2": "P3",
    "P3": "P4",
    "P4": "P5",
    "P5": "P6",
}

# Alternate mapping for data already using paper-consistent labels.
LAUREN_GROMACS_OXYGEN_TO_GRO_PHOSPHATE = {
    "O1": "P1",
    "O7": "P1",
    "O8": "P1",
    "O9": "P1",
    "O2": "P2",
    "O10": "P2",
    "O11": "P2",
    "O12": "P2",
    "O3": "P3",
    "O13": "P3",
    "O14": "P3",
    "O15": "P3",
    "O4": "P4",
    "O16": "P4",
    "O17": "P4",
    "O18": "P4",
    "O5": "P5",
    "O19": "P5",
    "O20": "P5",
    "O21": "P5",
    "O6": "P6",
    "O22": "P6",
    "O23": "P6",
    "O24": "P6",
}

LAUREN_GRO_TO_PAPER_PHOSPHATE = {
    "P1": "P1",
    "P2": "P2",
    "P3": "P3",
    "P4": "P4",
    "P5": "P5",
    "P6": "P6",
}

BRIDGING_OXYGEN_TOKENS_DEFAULT = {"O", "O1", "O2", "O3", "O4", "O5"}
BRIDGING_OXYGEN_TOKENS_LAUREN = {"O1", "O2", "O3", "O4", "O5", "O6"}

# One-based oxygen labels used by parse_hbond_log_to_dataframe and comparison scripts.
PAPER_ONE_BASED_OXYGEN_TO_PHOSPHATE = {
    "O1": "P1",
    "O7": "P1",
    "O8": "P1",
    "O9": "P1",
    "O2": "P2",
    "O10": "P2",
    "O11": "P2",
    "O12": "P2",
    "O3": "P3",
    "O13": "P3",
    "O14": "P3",
    "O15": "P3",
    "O4": "P4",
    "O16": "P4",
    "O17": "P4",
    "O18": "P4",
    "O5": "P5",
    "O19": "P5",
    "O20": "P5",
    "O21": "P5",
    "O6": "P6",
    "O22": "P6",
    "O23": "P6",
    "O24": "P6",
}


def mapping_tables(*, lauren_labels: bool = False) -> tuple[dict[str, str], dict[str, str], set[str]]:
    if lauren_labels:
        return (
            LAUREN_GROMACS_OXYGEN_TO_GRO_PHOSPHATE,
            LAUREN_GRO_TO_PAPER_PHOSPHATE,
            BRIDGING_OXYGEN_TOKENS_LAUREN,
        )
    return (
        GROMACS_OXYGEN_TO_GRO_PHOSPHATE,
        GRO_TO_PAPER_PHOSPHATE,
        BRIDGING_OXYGEN_TOKENS_DEFAULT,
    )


def map_gromacs_oxygen_to_paper_phosphate(oxygen_token: str | None, *, lauren_labels: bool = False) -> str | None:
    if not oxygen_token:
        return None
    oxygen_to_gro, gro_to_paper, _ = mapping_tables(lauren_labels=lauren_labels)
    gro_pg = oxygen_to_gro.get(oxygen_token)
    if not gro_pg:
        return None
    return gro_to_paper.get(gro_pg, gro_pg)
