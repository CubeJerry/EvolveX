import sys
from pathlib import Path
from types import SimpleNamespace
import types


def install_import_stubs():
    if "pandas" not in sys.modules:
        pandas = types.ModuleType("pandas")
        pandas.read_csv = lambda *args, **kwargs: None
        sys.modules["pandas"] = pandas

    if "dask.distributed" not in sys.modules:
        dask = types.ModuleType("dask")
        distributed = types.ModuleType("dask.distributed")
        distributed.as_completed = lambda *args, **kwargs: None
        distributed.Client = object
        distributed.LocalCluster = object
        distributed.wait = lambda *args, **kwargs: None
        sys.modules["dask"] = dask
        sys.modules["dask.distributed"] = distributed

    if "dask_jobqueue" not in sys.modules:
        dask_jobqueue = types.ModuleType("dask_jobqueue")
        dask_jobqueue.SLURMCluster = object
        sys.modules["dask_jobqueue"] = dask_jobqueue

    if "Bio" not in sys.modules:
        bio = types.ModuleType("Bio")
        bio.SeqIO = types.SimpleNamespace(parse=lambda *args, **kwargs: [])
        bio.BiopythonParserWarning = Warning

        pdb = types.ModuleType("Bio.PDB")

        class PDBParser:
            def __init__(self, *args, **kwargs):
                pass

            def get_structure(self, *args, **kwargs):
                return types.SimpleNamespace(get_residues=lambda: [])

        pdb.PDBParser = PDBParser

        pdb_exceptions = types.ModuleType("Bio.PDB.PDBExceptions")
        pdb_exceptions.PDBConstructionWarning = Warning

        data = types.ModuleType("Bio.Data")
        iupac = types.ModuleType("Bio.Data.IUPACData")
        iupac.protein_letters_3to1 = {}

        sys.modules["Bio"] = bio
        sys.modules["Bio.PDB"] = pdb
        sys.modules["Bio.PDB.PDBExceptions"] = pdb_exceptions
        sys.modules["Bio.Data"] = data
        sys.modules["Bio.Data.IUPACData"] = iupac


install_import_stubs()


sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from evolvex import search_algorithms


class ImmediateExecutor:
    def submit(self, func, *args):
        return func(*args)

    def cancel(self, futures):
        assert all(future is not None for future in futures)


def immediate_as_completed(futures, with_results=False):
    assert with_results is True
    for future in futures:
        yield future, future


def make_model(name, backbone="backbone"):
    return SimpleNamespace(name=name, backbone_PDB_file_name=backbone)


def test_ga_search_writes_each_mc_iteration_before_starting_next(monkeypatch, tmp_path):
    events = []

    def fake_make_MC_step(model, nth_iteration, iteration_fraction, model_PDB_files_dir, GLOBALS):
        events.append(("submit_mc", nth_iteration, model.name))
        return model, {
            "step": ["MC"],
            "nth_iteration": [nth_iteration],
            "nth_model": [model.name],
        }

    def fake_make_recombination_step(model_1, model_2, nth_iteration, iteration_fraction, model_PDB_files_dir, GLOBALS):
        events.append(("submit_recombination", nth_iteration, model_1.name, model_2.name))
        return model_1, model_2, {
            "step": ["recombination", "recombination"],
            "nth_iteration": [nth_iteration, nth_iteration],
            "nth_model": [model_1.name, model_2.name],
        }

    def fake_write_generated_models_info(generated_models_info, generated_models_info_file_handle):
        for step, nth_iteration, nth_model in zip(
            generated_models_info["step"],
            generated_models_info["nth_iteration"],
            generated_models_info["nth_model"],
        ):
            events.append(("write", step, nth_iteration, nth_model))

    monkeypatch.setattr(search_algorithms, "as_completed", immediate_as_completed)
    monkeypatch.setattr(search_algorithms, "make_MC_step", fake_make_MC_step)
    monkeypatch.setattr(search_algorithms, "make_recombination_step", fake_make_recombination_step)
    monkeypatch.setattr(search_algorithms, "write_generated_models_info", fake_write_generated_models_info)

    search_algorithms.GA_search(
        ImmediateExecutor(),
        [make_model("model_1"), make_model("model_2")],
        tmp_path / "generated_models_info.csv",
        tmp_path / "model_pdbs",
        SimpleNamespace(recombine_every_nth_iteration=3, max_iterations=6),
    )

    assert [event[1] for event in events if event[0] == "submit_mc"] == [1, 1, 2, 2, 4, 4, 5, 5]
    assert [event[1] for event in events if event[0] == "submit_recombination"] == [3, 6]

    writes_for_iteration_1 = [i for i, event in enumerate(events) if event[:3] == ("write", "MC", 1)]
    submits_for_iteration_2 = [i for i, event in enumerate(events) if event[:2] == ("submit_mc", 2)]
    assert max(writes_for_iteration_1) < min(submits_for_iteration_2)

    writes_for_iteration_4 = [i for i, event in enumerate(events) if event[:3] == ("write", "MC", 4)]
    submits_for_iteration_5 = [i for i, event in enumerate(events) if event[:2] == ("submit_mc", 5)]
    assert max(writes_for_iteration_4) < min(submits_for_iteration_5)
