"""Legacy-run isolation tests."""

import json


def test_classify_legacy_and_v4_runs(tmp_path):
    from slingshot.v4.runs import classify_run, discover_v4_runs

    legacy = tmp_path / "results_old"
    legacy.mkdir()
    (legacy / "results.pkl").write_bytes(b"legacy")
    current = tmp_path / "v4_current"
    current.mkdir()
    (current / "manifest.json").write_text(
        json.dumps(
            {
                "schema_version": 4,
                "legacy_science_model": False,
            }
        ),
        encoding="utf-8",
    )
    assert classify_run(legacy)["legacy_science_model"] is True
    assert classify_run(current)["legacy_science_model"] is False
    assert discover_v4_runs(tmp_path) == [current]
