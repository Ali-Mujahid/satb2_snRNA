"""Tiny invented transport fixtures only. No study analysis, files or ZIPs used."""
import importlib.util
import json
import tempfile
import unittest
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pack = load("pack", HERE / "package_archives.py")
restore = load("restore", HERE / "extract_reviewed_archive.py")


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.base = Path(self.tmp.name).resolve()
        self.native = self.base / "native"
        self.native.mkdir()
        self.input = self.native / "counts.bin"
        self.input.write_bytes(b"\x00\x01\xffEXACT\n" * 200)
        self.seal = self.native / "FINAL_VALIDATION.json"
        self.seal.write_text(json.dumps({"schema": "satb2-composition-final-validation/v1", "status": "COMPLETE_WITH_INCONCLUSIVE_OR_UNTESTABLE_SCOPES", "outputs": [{**pack.record(self.input), "path": "counts.bin"}]}))
        self.review = self.base / "actual_review.md"
        self.review.write_text("This is an invented unit-test fixture, not a study review.\n")
        self.cfg = {"schema": "satb2-followups-archive-delivery/v1", "approved_for_archives": True,
                    "allowed_source_roots": [str(self.base)], "allowed_output_roots": [str(self.base)], "destination_dir": str(self.base / "output"), "maximum_archive_bytes": 2_000_000,
                    "components": [{"id": "first", "archive_name": "component.zip", "role": "native_result", "readme": "INVENTED UNIT TEST: inconclusive, not biological data.\n", "gates": [{"kind": "composition", "root": str(self.native), "archive_prefix": "native", "seal": pack.record(self.seal), "review_receipts": [pack.record(self.review)]}]}]}
        self.config = self.base / "config.json"

    def tearDown(self):
        self.tmp.cleanup()

    def save(self):
        self.config.write_text(json.dumps(self.cfg))

    def test_exact_bytes_and_roundtrip(self):
        self.save()
        result = pack.build(self.config)
        archive = self.base / "output/component.zip"
        self.assertEqual(result["status"], "BUILT_AND_READBACK_VERIFIED_NOT_UPLOADED")
        self.assertEqual(self.input.read_bytes(), zipfile.ZipFile(archive).read("native/counts.bin"))
        out = restore.extract(archive, self.base / "restored", pack.sha(archive))
        self.assertEqual((out / "native/counts.bin").read_bytes(), self.input.read_bytes())

    def test_deterministic_archive_and_fresh_only(self):
        self.save(); pack.build(self.config)
        first = (self.base / "output/component.zip").read_bytes()
        with self.assertRaises(ValueError):
            pack.build(self.config)
        self.cfg["destination_dir"] = str(self.base / "other")
        self.save(); pack.build(self.config)
        self.assertEqual(first, (self.base / "other/component.zip").read_bytes())

    def test_unapproved_and_partial_fail_before_output(self):
        self.cfg["approved_for_archives"] = False; self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        self.cfg["approved_for_archives"] = True
        d = json.loads(self.seal.read_text()); d["status"] = "RUNNING"
        self.seal.write_text(json.dumps(d)); self.cfg["components"][0]["gates"][0]["seal"] = pack.record(self.seal); self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        self.assertFalse((self.base / "output").exists())

    def test_review_hash_and_native_hash_are_bound(self):
        self.save(); self.review.write_text("changed")
        with self.assertRaises(ValueError): pack.build(self.config)
        self.cfg["components"][0]["gates"][0]["review_receipts"] = [pack.record(self.review)]
        self.save(); self.input.write_bytes(b"changed")
        with self.assertRaises(ValueError): pack.build(self.config)
        self.assertFalse((self.base / "output").exists())

    def test_traversal_symlinks_and_outside_roots(self):
        for bad in ["../evil", "/evil", "a/../evil", "C:/evil", "a\\evil", "a//evil", "a/evil "]:
            with self.assertRaises(ValueError): pack.safe_member(bad)
        linked = self.base / "linked"; linked.symlink_to(self.input)
        self.cfg["components"][0]["files"] = [{"receipt": pack.record(linked), "member": "source/linked"}]; self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        with self.assertRaises(ValueError): pack.canonical(self.input, [self.base / "elsewhere"])

    def test_shared_dedup_and_alias_restore(self):
        component = self.cfg["components"][0]
        component["role"] = "shared_frozen_inputs_and_source"
        component["deduplicate_identical_files"] = True
        component["files"] = [{"receipt": pack.record(self.input), "member": "second/counts.bin"}]
        self.save(); result = pack.build(self.config)
        self.assertEqual(result["archives"][0]["requested_members"] - result["archives"][0]["stored_members"], 1)
        out = restore.extract(self.base / "output/component.zip", self.base / "restored")
        self.assertEqual((out / "second/counts.bin").read_bytes(), self.input.read_bytes())

    def test_native_models_cannot_dedup_paths(self):
        self.cfg["components"][0]["deduplicate_identical_files"] = True; self.save()
        with self.assertRaises(ValueError): pack.build(self.config)

    def test_display_receipt_required_and_component_budgets(self):
        gate = self.cfg["components"][0]["gates"][0]; gate["kind"] = "layout"
        self.seal.write_text(json.dumps({"schema": "satb2-cellchat-layout-only-final/v1", "status": "COMPLETE_PENDING_VISUAL_REVIEW", "outputs": [{**pack.record(self.input), "path": "counts.bin"}]}))
        gate["seal"] = pack.record(self.seal); self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        display = self.base / "display.json"
        display.write_text(json.dumps({"schema": "satb2-cellchat-layout-display-acceptance/v1", "approved_for_display": True, "layout_final_validation": {"sha256": gate["seal"]["sha256"]}}))
        gate["display_review_receipts"] = [pack.record(display)]
        self.cfg["maximum_archive_bytes"] = 100; self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        self.assertFalse((self.base / "output").exists())

    def test_mutable_environment_and_overlap_rejected(self):
        self.cfg["forbidden_source_roots"] = [str(self.native)]
        self.save()
        with self.assertRaises(ValueError): pack.build(self.config)
        self.cfg.pop("forbidden_source_roots")
        self.cfg["destination_dir"] = str(self.native / "archive")
        self.save()
        with self.assertRaises(ValueError): pack.build(self.config)

    def test_restore_rejects_zip_traversal_and_wrong_digest(self):
        self.save(); pack.build(self.config)
        archive = self.base / "output/component.zip"
        with self.assertRaises(ValueError): restore.extract(archive, self.base / "bad", "0" * 64)
        with zipfile.ZipFile(archive, "a") as z: z.writestr("../escape", "untrusted")
        with self.assertRaises(ValueError): restore.extract(archive, self.base / "bad")

    def test_selected_components_do_not_require_unselected_partial_seal(self):
        second = json.loads(json.dumps(self.cfg["components"][0]))
        second.update(id="second", archive_name="second.zip")
        second["gates"][0]["seal"]["path"] = str(self.base / "absent_pending_FINAL.json")
        self.cfg["components"].append(second); self.save()
        result = pack.build(self.config, ["first"])
        self.assertEqual(len(result["archives"]), 1)
        self.assertFalse((self.base / "output/second.zip").exists())

    def test_final_release_requires_all_components_across_fresh_invocations(self):
        second = json.loads(json.dumps(self.cfg["components"][0]))
        second.update(id="second", archive_name="second.zip")
        self.cfg["components"].append(second); self.save()
        pack.build(self.config, ["first"])
        first_receipt = self.base / "output/ARCHIVE_DELIVERY_RECEIPT.json"
        with self.assertRaises(ValueError): pack.verify_deliveries(self.config, [first_receipt])
        self.cfg["destination_dir"] = str(self.base / "second_output"); self.save()
        pack.build(self.config, ["second"])
        second_receipt = self.base / "second_output/ARCHIVE_DELIVERY_RECEIPT.json"
        verified = pack.verify_deliveries(self.config, [first_receipt, second_receipt])
        self.assertEqual(set(verified["components"]), {"first", "second"})
        with self.assertRaises(ValueError): pack.verify_deliveries(self.config, [first_receipt, first_receipt, second_receipt])

    def test_executable_wrapper_permission_restored_without_execution(self):
        wrapper = self.base / "run.sh"
        wrapper.write_text("#!/bin/sh\n# never executed by this test\n")
        wrapper.chmod(0o755)
        self.cfg["components"][0]["files"] = [{"member": "source/run.sh", "receipt": pack.record(wrapper)}]
        self.save(); pack.build(self.config)
        restored = restore.extract(self.base / "output/component.zip", self.base / "restored")
        self.assertEqual((restored / "source/run.sh").stat().st_mode & 0o777, 0o755)


if __name__ == "__main__":
    unittest.main()
