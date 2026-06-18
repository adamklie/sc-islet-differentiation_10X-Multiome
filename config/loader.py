import csv
import yaml
from pathlib import Path

CONFIG_DIR = Path(__file__).parent


# ---------------------------------------------------------------------------
# YAML config loading
# ---------------------------------------------------------------------------

def load_yaml_config(name):
    """
    Load any YAML config file from this config directory by name.

    Parameters
    ----------
    name : str
        Config path relative to the config directory, with or without
        the ``.yaml`` extension.  Supports subdirectories
        (e.g., ``'colors/sc-islet-differentiation_10X-Multiome'``).

    Returns
    -------
    dict
        Parsed YAML contents.
    """
    if not name.endswith(".yaml"):
        name = name + ".yaml"
    yaml_path = CONFIG_DIR / name
    if not yaml_path.exists():
        raise FileNotFoundError(f"No config found: {yaml_path}")
    with open(yaml_path) as f:
        return yaml.safe_load(f)


def list_configs():
    """Return available config YAML paths relative to config dir (without extension)."""
    configs = []
    for f in sorted(CONFIG_DIR.rglob("*.yaml")):
        rel = f.relative_to(CONFIG_DIR)
        configs.append(str(rel.with_suffix("")))
    return configs


# ---------------------------------------------------------------------------
# TSV metadata loading (generic)
# ---------------------------------------------------------------------------

def _load_metadata_tsv(path):
    """Load a TSV metadata file as a list of row dicts."""
    rows = []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            if "display_order" in row:
                row["display_order"] = int(row["display_order"])
            rows.append(row)
    return rows


# ---------------------------------------------------------------------------
# Cell type metadata (single source of truth)
# ---------------------------------------------------------------------------

def load_cell_type_metadata(path=None):
    """
    Load the cell type metadata TSV as a list of row dicts.

    Returns
    -------
    list of dict
        Each dict has keys: ``cell_type``, ``grouping``, ``color``,
        ``display_order``.
    """
    path = Path(path) if path else CONFIG_DIR / "cell_type_metadata.tsv"
    return _load_metadata_tsv(path)


def load_cell_type_colors(path=None):
    """
    Return cell type colors and display order from the metadata TSV.

    Returns
    -------
    dict
        ``{"celltype_colors": {ct: hex, ...}, "celltype_order": [ct, ...]}``
    """
    rows = load_cell_type_metadata(path)
    rows.sort(key=lambda r: r["display_order"])
    return {
        "celltype_colors": {r["cell_type"]: r["color"] for r in rows},
        "celltype_order": [r["cell_type"] for r in rows],
    }


def load_grouping(path=None):
    """
    Return a dict mapping each cell type to its biological group.

    Returns
    -------
    dict
        ``{cell_type: group_name}``
    """
    rows = load_cell_type_metadata(path)
    return {r["cell_type"]: r["grouping"] for r in rows}


def get_group_members(grouping, group_name):
    """Return a list of cell types belonging to a group."""
    return [ct for ct, grp in grouping.items() if grp == group_name]


def get_group_names(grouping):
    """Return unique group names in insertion order."""
    seen = {}
    for grp in grouping.values():
        seen.setdefault(grp, None)
    return list(seen)


# ---------------------------------------------------------------------------
# Stage metadata
# ---------------------------------------------------------------------------

def load_stage_metadata(path=None):
    """Load the differentiation stage metadata TSV as a list of row dicts."""
    path = Path(path) if path else CONFIG_DIR / "stage_metadata.tsv"
    return _load_metadata_tsv(path)


def load_stage_colors(path=None):
    """
    Return differentiation stage colors and display order.

    Returns
    -------
    dict
        ``{"stage_colors": {stage: hex, ...}, "stage_order": [stage, ...]}``
    """
    rows = load_stage_metadata(path)
    rows.sort(key=lambda r: r["display_order"])
    return {
        "stage_colors": {r["stage"]: r["color"] for r in rows},
        "stage_order": [r["stage"] for r in rows],
    }


# ---------------------------------------------------------------------------
# Endocrine metadata
# ---------------------------------------------------------------------------

def load_endocrine_metadata(path=None):
    """Load the endocrine classification metadata TSV as a list of row dicts."""
    path = Path(path) if path else CONFIG_DIR / "endocrine_metadata.tsv"
    return _load_metadata_tsv(path)


def load_endocrine_colors(path=None):
    """
    Return endocrine classification colors and display order.

    Returns
    -------
    dict
        ``{"endocrine_colors": {label: hex, ...}, "endocrine_order": [label, ...]}``
    """
    rows = load_endocrine_metadata(path)
    rows.sort(key=lambda r: r["display_order"])
    return {
        "endocrine_colors": {r["endocrine"]: r["color"] for r in rows},
        "endocrine_order": [r["endocrine"] for r in rows],
    }


# ---------------------------------------------------------------------------
# Batch metadata
# ---------------------------------------------------------------------------

def load_batch_metadata(path=None):
    """Load the batch metadata TSV as a list of row dicts."""
    path = Path(path) if path else CONFIG_DIR / "batch_metadata.tsv"
    return _load_metadata_tsv(path)


def load_batch_colors(path=None):
    """
    Return batch colors and display order.

    Returns
    -------
    dict
        ``{"batch_colors": {batch: hex, ...}, "batch_order": [batch, ...]}``
    """
    rows = load_batch_metadata(path)
    rows.sort(key=lambda r: r["display_order"])
    return {
        "batch_colors": {r["batch"]: r["color"] for r in rows},
        "batch_order": [r["batch"] for r in rows],
    }
