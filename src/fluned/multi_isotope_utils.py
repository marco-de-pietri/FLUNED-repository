import pickle
import xml.etree.ElementTree as ET
from typing import Any

import h5py


def get_isotope_reactions_dicts(file_path):
    """
    Read depletion results and return mappings for isotopes and reactions.
    """

    with h5py.File(file_path, "r") as f:
        nuc_idx = pickle.loads(f["index_nuc"][...].tobytes())
        rx_idx = pickle.loads(f["index_rx"][...].tobytes())

    return nuc_idx, rx_idx


def map_targets_to_channels(
    xml_file: str, isotopes_index: dict[str, int], reaction_rates_index: dict[str, int]
) -> dict[str, list[dict[str, Any]]]:
    """
    Parse `xml_file` and map each target nuclide to its parent/reaction channels.
    """
    tree = ET.parse(xml_file)
    root = tree.getroot()
    result: dict[str, list[dict[str, Any]]] = {}

    for nuclide in root.findall("nuclide"):
        parent = nuclide.get("name")
        if parent not in isotopes_index:
            continue
        parent_idx = isotopes_index[parent]

        for rx in nuclide.findall("reaction"):
            rtype = rx.get("type")
            target = rx.get("target")

            if target is None:
                raise ValueError("target attribute missing in reaction element")

            if rtype not in reaction_rates_index:
                continue

            reaction_idx = reaction_rates_index[rtype]
            channel = {
                "parent_nuclide": parent,
                "reaction": rtype,
                "channel_index": [parent_idx, reaction_idx],
            }
            if target in result:
                result[target].append(channel)
            else:
                result[target] = [channel]

    return result


def filter_channels(
    xml_file: str, channels: dict[str, list[dict[str, Any]]]
) -> dict[str, list[dict[str, Any]]]:
    """
    Keep only targets that have a photon source in the chain XML.
    """

    tree = ET.parse(xml_file)
    root = tree.getroot()

    valid_sources = set()
    for nu in root.findall("nuclide"):
        name = nu.get("name")
        for src in nu.findall("source"):
            if src.get("particle") in ("photon",):
                valid_sources.add(name)
                break

    return {
        target: chans for target, chans in channels.items() if target in valid_sources
    }
