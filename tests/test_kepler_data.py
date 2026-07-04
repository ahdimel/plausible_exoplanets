"""Tests for kepler_data (DR25 ingestion): fiducial cut logic, KeplerTarget
parsing (14-duration CDPP order), KOI->host filtering, null-feh defaulting,
and score-cut variants. No network: runs against tests/fixtures/*, ~20 real
rows copied from the 2026-07-03 snapshots (two cells blanked on purpose —
kepid 10000785 feh and kepid 10000797 rrmscdpp06p0 — because real DR25 rows
never have those nulls in isolation)."""

import math
from pathlib import Path

import pytest

from exoverse.kepler_data import (
    CDPP_COLUMNS, CDPP_DURATIONS_HR, Cuts, FIDUCIAL, KOI_SNAPSHOT,
    STELLAR_SNAPSHOT, load_koi_systems, load_stellar_targets,
)

FIXTURES = Path(__file__).parent / "fixtures"

# kepids in the stellar fixture that pass the fiducial stellar cuts
PASSING_KEPIDS = {
    3128552, 1432789, 757450,             # KOI hosts (multi-3, multi-2, single)
    10811496, 3548044, 3112129, 8218274,  # hosts of KOI-cut-rejected KOIs
    10000785,                             # feh blank -> 0.0, otherwise passes
    10000800, 10000823, 10000827, 10000876,
}
# kepid -> the stellar cut that excludes it
EXCLUDED_KEPIDS = {
    10068797: "teff 3692 < 3900",
    10028433: "teff 7535 > 7300",
    10000981: "logg 3.490 < 4.0",
    9471974: "logg 3.907 < 4.0 (KOI K00119.01 host)",
    10001013: "dataspan 92.58 <= 365",
    10000672: "mass null",
    100000925: "dutycycle/dataspan/cdpp all null",
    10000797: "one cdpp column (rrmscdpp06p0) null",
}


@pytest.fixture
def data_dir(tmp_path):
    """Stage the fixtures under the canonical snapshot names, with a
    provenance comment line like the real snapshots carry."""
    for src, dst in [("dr25_stellar_fixture.csv", STELLAR_SNAPSHOT),
                     ("dr25_koi_fixture.csv", KOI_SNAPSHOT)]:
        text = (FIXTURES / src).read_text()
        (tmp_path / dst).write_text("# fixture rows from 2026-07-03 snapshot\n"
                                    + text)
    return tmp_path


# ------------------------------------------------------------ stellar cuts
def test_stellar_cuts_include_and_exclude(data_dir):
    kepids = {t.kepid for t in load_stellar_targets(data_dir)}
    assert kepids == PASSING_KEPIDS
    for kepid, reason in EXCLUDED_KEPIDS.items():
        assert kepid not in kepids, f"{kepid} should fail: {reason}"


def test_stellar_cut_edges(data_dir):
    # widening each cut readmits exactly the row it was excluding
    wide_teff = Cuts(teff_min=3600.0, teff_max=7600.0)
    kepids = {t.kepid for t in load_stellar_targets(data_dir, wide_teff)}
    assert {10068797, 10028433} <= kepids
    kepids = {t.kepid for t in load_stellar_targets(data_dir, Cuts(logg_min=3.4))}
    assert {10000981, 9471974} <= kepids
    kepids = {t.kepid for t in load_stellar_targets(data_dir, Cuts(dataspan_min=90.0))}
    assert 10001013 in kepids
    # null rows stay out no matter how loose the cuts are
    loose = Cuts(teff_min=0.0, teff_max=1e5, logg_min=-10.0, dataspan_min=0.0)
    kepids = {t.kepid for t in load_stellar_targets(data_dir, loose)}
    assert not {10000672, 100000925, 10000797} & kepids


# --------------------------------------------------------- target parsing
def test_kepler_target_parsing(data_dir):
    t = {x.kepid: x for x in load_stellar_targets(data_dir)}[3128552]
    assert (t.teff, t.logg, t.feh) == (5509.0, 4.308, 0.160)
    assert (t.radius, t.mass, t.kepmag) == (1.116, 0.9240, 14.523)
    assert (t.dutycycle, t.dataspan) == (0.6978, 1458.9310)
    assert len(t.cdpp_ppm) == len(CDPP_DURATIONS_HR) == len(CDPP_COLUMNS) == 14
    # documented duration order: index 0 = 1.5 h, 7 = 6.0 h, 13 = 15.0 h
    assert CDPP_DURATIONS_HR == (1.5, 2.0, 2.5, 3.0, 3.5, 4.5, 5.0, 6.0,
                                 7.5, 9.0, 10.5, 12.0, 12.5, 15.0)
    assert t.cdpp_ppm[0] == 229.670
    assert t.cdpp_ppm[7] == 141.746   # rrmscdpp06p0
    assert t.cdpp_ppm[13] == 115.737  # rrmscdpp15p0
    assert all(math.isfinite(c) for c in t.cdpp_ppm)
    assert CDPP_COLUMNS[0] == "rrmscdpp01p5"
    assert CDPP_COLUMNS[7] == "rrmscdpp06p0"
    assert CDPP_COLUMNS[13] == "rrmscdpp15p0"


def test_null_feh_defaults_to_solar(data_dir):
    t = {x.kepid: x for x in load_stellar_targets(data_dir)}[10000785]
    assert t.feh == 0.0
    assert (t.teff, t.radius, t.mass) == (5333.0, 0.650, 0.6350)


# ---------------------------------------------------------------- KOI cuts
def test_koi_systems_grouping_and_sorting(data_dir):
    systems = load_koi_systems(data_dir)
    assert set(systems) == {3128552, 1432789, 757450}
    assert [k.kepoi_name for k in systems[3128552]] == [
        "K02055.03", "K02055.02", "K02055.01"]  # sorted by period
    periods = [k.period for k in systems[3128552]]
    assert periods == sorted(periods)
    assert len(systems[1432789]) == 2
    single = systems[757450][0]
    assert single.kepoi_name == "K00889.01"
    assert (single.period, single.prad) == (8.884922995, 10.51)
    assert single.score == 0.999
    assert single.disposition == "CONFIRMED"


def test_koi_cut_exclusions(data_dir):
    names = {k.kepoi_name
             for kois in load_koi_systems(data_dir).values() for k in kois}
    assert "K00119.01" not in names   # host 9471974 fails stellar logg cut
    assert "K00755.01" not in names   # host absent from the stellar table
    assert "K00753.01" not in names   # FALSE POSITIVE (host passes stellar)
    assert "K02194.03" not in names   # score 0.496 <= 0.5
    assert "K04144.01" not in names   # period 0.4877 < 0.5 d
    assert "K01064.01" not in names   # prad 67.35 > 30 Re


def test_koi_score_variants(data_dir):
    # score > 0.0 readmits exactly the low-score candidate
    systems = load_koi_systems(data_dir, Cuts(koi_score_min=0.0))
    names = {k.kepoi_name for kois in systems.values() for k in kois}
    assert "K02194.03" in names
    assert 3548044 in systems
    assert "K00753.01" not in names   # FP score 0.0 still out (strict >)
    # score > 0.9 drops K00992.01 (0.848) but keeps its sibling
    systems = load_koi_systems(data_dir, Cuts(koi_score_min=0.9))
    assert [k.kepoi_name for k in systems[1432789]] == ["K00992.02"]
    assert len(systems[3128552]) == 3


def test_koi_period_and_prad_variants(data_dir):
    systems = load_koi_systems(data_dir, Cuts(koi_period_min=0.2))
    names = {k.kepoi_name for kois in systems.values() for k in kois}
    assert "K04144.01" in names
    systems = load_koi_systems(data_dir, Cuts(koi_prad_max=100.0))
    names = {k.kepoi_name for kois in systems.values() for k in kois}
    assert "K01064.01" in names


def test_missing_snapshot_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_stellar_targets(tmp_path)


def test_fiducial_constants():
    assert FIDUCIAL == Cuts()
    assert (FIDUCIAL.teff_min, FIDUCIAL.teff_max) == (3900.0, 7300.0)
    assert FIDUCIAL.logg_min == 4.0
    assert FIDUCIAL.dataspan_min == 365.0
    assert FIDUCIAL.koi_dispositions == ("CONFIRMED", "CANDIDATE")
    assert FIDUCIAL.koi_score_min == 0.5
    assert (FIDUCIAL.koi_period_min, FIDUCIAL.koi_period_max) == (0.5, 640.0)
    assert FIDUCIAL.koi_prad_max == 30.0
