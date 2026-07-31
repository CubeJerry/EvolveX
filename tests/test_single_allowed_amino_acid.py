import statistics

import pandas as pd

from evolvex.model_generation import (
    _sample_variance_or_zero,
    get_acceptable_positions_mut_names_map,
    get_allowed_mutations_per_position_maps,
)
from evolvex.search_algorithms import get_random_mut_name


def test_sample_variance_or_zero_handles_singleton_without_changing_regular_variance():
    assert _sample_variance_or_zero([0.25]) == 0.0

    regular_values = [-0.5, 0.0, 0.5]
    assert _sample_variance_or_zero(regular_values) == statistics.variance(
        regular_values
    )


def test_single_d_allowed_amino_acid_survives_handoff_and_sampling(tmp_path):
    summary = pd.DataFrame(
        [
            {
                "position": "42",
                "binding_ddG": -0.4,
                "antibody_stability_ddG": 0.2,
                "complex_stability_ddG": 0.1,
                "original_residue": "A",
                "MakeAla": "N",
            }
        ],
        index=["AB42D"],
    )
    summary_path = tmp_path / "all_mutations_summary.csv"
    summary.to_csv(summary_path)

    allowed_mut_names, allowed_amino_acids, make_ala_positions = (
        get_allowed_mutations_per_position_maps(
            PDB_name="single_d_test",
            all_mutations_summary_file_path=summary_path,
        )
    )

    assert allowed_mut_names == {"42": ["AB42D"]}
    assert allowed_amino_acids == {"42": {"D"}}
    assert make_ala_positions == set()
    assert get_random_mut_name(["AB42"], allowed_amino_acids) == "AB42D"


def test_regular_multi_amino_acid_filtering_is_unchanged():
    summary = pd.DataFrame(
        [
            {
                "position": "53",
                "binding_ddG": 0.0,
                "antibody_stability_ddG": 0.0,
                "complex_stability_ddG": 0.0,
                "original_residue": "A",
            },
            {
                "position": "53",
                "binding_ddG": 0.1,
                "antibody_stability_ddG": 0.1,
                "complex_stability_ddG": 0.0,
                "original_residue": "A",
            },
            {
                "position": "53",
                "binding_ddG": 0.2,
                "antibody_stability_ddG": 0.2,
                "complex_stability_ddG": 0.0,
                "original_residue": "A",
            },
        ],
        index=["AB53D", "AB53F", "AB53N"],
    )

    acceptable = get_acceptable_positions_mut_names_map(
        summary,
        PDB_name="regular_position_test",
    )

    # The pre-existing low-variance large-hydrophobic filter still removes F,
    # while ordinary D and N choices remain available.
    assert acceptable["53"] == ["AB53D", "AB53N"]
