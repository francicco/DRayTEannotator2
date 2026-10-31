import hashlib
import io
import json
import os
import tarfile
from pathlib import Path

import pytest

from drayte.storage.archive import ArchiveStore


def make_tree(root):
    tree = root / "extension/extensionwork"
    (tree / "family").mkdir(parents=True)
    (tree / "family/repam-ranges.tsv").write_text("start\tend\n1\t4\n")
    (tree / "family/repam-repseq.fa").write_text(">f\nACGT\n")
    (tree / "empty").mkdir()
    os.utime(tree / "family/repam-ranges.tsv", ns=(1234567890123456789, 1234567890123456789))
    return tree


def test_roundtrip_preserves_bytes_mtimes_modes_and_internal_symlinks(tmp_path):
    tree = make_tree(tmp_path)
    (tree / "linked").symlink_to("family/repam-repseq.fa")
    os.chmod(tree / "family/repam-ranges.tsv", 0o640)
    store = ArchiveStore(tmp_path)
    original = store._inventory("extension/extensionwork")
    record = store.archive("extension/extensionwork")
    assert not tree.exists()
    assert store.records()["extension/extensionwork"] == record
    assert store.restore_all() == ["extension/extensionwork"]
    assert store._inventory("extension/extensionwork") == original
    assert store.records() == {}
    assert not (tmp_path / record["archive"]).exists()
    assert (tree / "linked").is_symlink()


def test_file_roundtrip_and_rearchive_after_restore(tmp_path):
    file = tmp_path / "Hmel.fa.align"
    file.write_bytes(b"alignment\n" * 200)
    store = ArchiveStore(tmp_path)
    first = store.archive(file.name)
    assert store.archive(file.name) == first
    store.restore_all()
    second = store.archive(file.name)
    assert second["archive"] != first["archive"]
    store.restore_all()
    assert file.read_bytes() == b"alignment\n" * 200


def test_corrupt_archive_refuses_restore_without_writing_source(tmp_path):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    record = store.archive("extension/extensionwork")
    (tmp_path / record["archive"]).write_bytes(b"corrupt")
    with pytest.raises(ValueError, match="checksum mismatch"):
        store.restore_all()
    assert not tree.exists()
    assert store.records()


def test_missing_archive_fails_explicitly(tmp_path):
    make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    record = store.archive("extension/extensionwork")
    (tmp_path / record["archive"]).unlink()
    with pytest.raises(FileNotFoundError, match="missing"):
        store.records()


def test_partial_source_removal_is_recoverable(tmp_path, monkeypatch):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    import drayte.storage.archive as module

    real_remove = module.shutil.rmtree

    def interrupted_remove(path, *args, **kwargs):
        if Path(path) == tree:
            (tree / "family/repam-repseq.fa").unlink()
            raise OSError("simulated interrupted removal")
        return real_remove(path, *args, **kwargs)

    interrupted_remove.avoids_symlink_attacks = True
    monkeypatch.setattr(module.shutil, "rmtree", interrupted_remove)
    with pytest.raises(OSError, match="interrupted"):
        store.archive("extension/extensionwork")
    assert store.records()
    monkeypatch.setattr(module.shutil, "rmtree", real_remove)
    store.restore_all()
    assert (tree / "family/repam-repseq.fa").read_text() == ">f\nACGT\n"


def test_changed_existing_file_is_never_overwritten(tmp_path):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    store.archive("extension/extensionwork")
    file = tree / "family/repam-ranges.tsv"
    file.parent.mkdir(parents=True)
    file.write_text("new data\n")
    with pytest.raises(ValueError, match="Conflicting"):
        store.restore_all()
    assert file.read_text() == "new data\n"
    assert store.records()


def test_archive_verification_failure_retains_originals(tmp_path, monkeypatch):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    original = store._inventory("extension/extensionwork")

    def fail(record):
        raise ValueError("simulated verification error")

    monkeypatch.setattr(store, "_verify_archive", fail)
    with pytest.raises(ValueError, match="verification"):
        store.archive("extension/extensionwork")
    assert store._inventory("extension/extensionwork") == original
    assert store.records() == {}
    assert list(store.directory.iterdir()) == []


@pytest.mark.parametrize("source", [".", "../outside", "/tmp/outside", ".storage", "extension/../outside"])
def test_unsafe_source_paths_are_rejected(tmp_path, source):
    store = ArchiveStore(tmp_path)
    with pytest.raises(ValueError):
        store.archive(source)


def test_external_or_symlinked_sources_are_rejected(tmp_path):
    tree = make_tree(tmp_path)
    (tree / "bad").symlink_to("/etc/passwd")
    store = ArchiveStore(tmp_path)
    with pytest.raises(ValueError, match="escapes"):
        store.archive("extension/extensionwork")
    assert tree.exists()
    (tree / "bad").unlink()
    (tmp_path / "linked-source").symlink_to(tree, target_is_directory=True)
    with pytest.raises(ValueError, match="symlink"):
        store.archive("linked-source")


def test_archive_catalog_cannot_use_external_payload_or_traversal(tmp_path):
    make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    record = store.archive("extension/extensionwork")
    catalog = json.loads(store.manifest.read_text())
    record = catalog["archives"]["extension/extensionwork"]
    record["members"][0]["name"] = "../outside"
    store.manifest.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="Unsafe"):
        store.restore_all()


def test_malicious_tar_member_rejected_even_with_updated_payload_hash(tmp_path):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    record = store.archive("extension/extensionwork")
    payload = tmp_path / record["archive"]
    with tarfile.open(payload, "w:gz") as archive:
        member = tarfile.TarInfo("../outside")
        member.size = 3
        archive.addfile(member, io.BytesIO(b"bad"))
    catalog = json.loads(store.manifest.read_text())
    catalog["archives"]["extension/extensionwork"]["archive_sha256"] = hashlib.sha256(payload.read_bytes()).hexdigest()
    store.manifest.write_text(json.dumps(catalog))
    with pytest.raises(ValueError, match="Unsafe"):
        store.restore_all()
    assert not tree.exists()
    assert not (tmp_path.parent / "outside").exists()


def test_restore_refuses_symlinked_parent(tmp_path):
    make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    store.archive("extension/extensionwork")
    (tmp_path / "extension").rmdir()
    (tmp_path / "elsewhere").mkdir()
    (tmp_path / "extension").symlink_to(tmp_path / "elsewhere", target_is_directory=True)
    with pytest.raises(ValueError, match="ancestor"):
        store.restore_all()
    assert not list((tmp_path / "elsewhere").iterdir())


def test_archive_retains_selected_inputs_and_restores_full_tree(tmp_path):
    tree = make_tree(tmp_path)
    scratch = tree / "family/intermediate.fa"
    scratch.write_text("ACGT" * 20000)
    store = ArchiveStore(tmp_path)
    before = store._inventory("extension/extensionwork")
    retained = "extension/extensionwork/family/repam-ranges.tsv"
    store.archive("extension/extensionwork", retain=[retained])
    assert (tmp_path / retained).is_file()
    assert not scratch.exists()
    assert not (tree / "family/repam-repseq.fa").exists()
    store.restore_all()
    assert store._inventory("extension/extensionwork") == before


@pytest.mark.parametrize("retained", ["elsewhere", "extension/extensionwork/family", "../outside"])
def test_invalid_retained_names_never_remove_originals(tmp_path, retained):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    before = store._inventory("extension/extensionwork")
    with pytest.raises(ValueError):
        store.archive("extension/extensionwork", retain=[retained])
    assert store._inventory("extension/extensionwork") == before
    assert not store.records()


def test_modified_retained_file_blocks_restore(tmp_path):
    tree = make_tree(tmp_path)
    store = ArchiveStore(tmp_path)
    retained = "extension/extensionwork/family/repam-ranges.tsv"
    store.archive("extension/extensionwork", retain=[retained])
    (tmp_path / retained).write_text("changed\n")
    with pytest.raises(ValueError, match="Conflicting"):
        store.restore_all()
    assert not (tree / "family/repam-repseq.fa").exists()
