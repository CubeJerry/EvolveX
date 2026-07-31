import csv
from types import SimpleNamespace

from evolvex import search_algorithms


class ImmediateExecutor:
    def submit(self, function, *args):
        return function(*args)

    def cancel(self, futures):
        assert all(future is not None for future in futures)


def immediate_as_completed(futures, with_results=False):
    assert with_results is True
    for future in futures:
        yield future, future


def make_model(name):
    return SimpleNamespace(
        model_dir=SimpleNamespace(name=name),
        backbone_PDB_file_name="backbone",
    )


def count_iteration_rows(csv_path, iteration):
    if not csv_path.exists() or csv_path.stat().st_size == 0:
        return 0
    with csv_path.open(newline="") as handle:
        return sum(
            int(row["nth_iteration"]) == iteration
            for row in csv.DictReader(handle)
        )


def test_ga_search_streams_each_complete_iteration_before_starting_next(
    monkeypatch,
    tmp_path,
):
    csv_path = tmp_path / "generated_models_info.csv"

    def assert_previous_iteration_is_visible(nth_iteration):
        if nth_iteration > 1:
            assert count_iteration_rows(csv_path, nth_iteration - 1) == 2

    def fake_make_mc_step(
        model,
        nth_iteration,
        iteration_fraction,
        model_pdb_files_dir,
        globals_,
    ):
        assert_previous_iteration_is_visible(nth_iteration)
        return model, {
            "backbone_PDB_file_name": [model.backbone_PDB_file_name],
            "nth_model": [model.model_dir.name],
            "step": ["MC"],
            "nth_iteration": [nth_iteration],
        }

    def fake_make_recombination_step(
        model_1,
        model_2,
        nth_iteration,
        iteration_fraction,
        model_pdb_files_dir,
        globals_,
    ):
        assert_previous_iteration_is_visible(nth_iteration)
        return model_1, model_2, {
            "backbone_PDB_file_name": ["backbone", "backbone"],
            "nth_model": [model_1.model_dir.name, model_2.model_dir.name],
            "step": ["recombination", "recombination"],
            "nth_iteration": [nth_iteration, nth_iteration],
        }

    monkeypatch.setattr(search_algorithms, "as_completed", immediate_as_completed)
    monkeypatch.setattr(search_algorithms, "make_MC_step", fake_make_mc_step)
    monkeypatch.setattr(
        search_algorithms,
        "make_recombination_step",
        fake_make_recombination_step,
    )

    search_algorithms.GA_search(
        ImmediateExecutor(),
        [make_model("model_1"), make_model("model_2")],
        csv_path,
        tmp_path / "model_pdbs",
        SimpleNamespace(recombine_every_nth_iteration=3, max_iterations=6),
    )

    with csv_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert len(rows) == 12
    assert {
        iteration: count_iteration_rows(csv_path, iteration)
        for iteration in range(1, 7)
    } == {iteration: 2 for iteration in range(1, 7)}
