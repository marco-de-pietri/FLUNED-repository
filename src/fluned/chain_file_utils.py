import xml.etree.ElementTree as ET
from typing import Any


def map_targets_to_channels(
    xml_file: str, isotopes_index: dict[str, int], reaction_rates_index: dict[str, int]
) -> dict[str, list[dict[str, Any]]]:
    """
    Parses the XML at `xml_file` and returns a dict mapping each target atom
    to a list of channels. Each channel is a dict with:
      - 'parent_nuclide': name of the nuclide
      - 'reaction': reaction type
      - 'channel_index': [parent_isotope_index, reaction_rate_index]

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
    Given a mapping of target nuclides to reaction channels, remove any
    targets for which the corresponding <nuclide> element in the XML
    has no <source> child.
    """

    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Identify nuclides with a <source> of particle 'photon'
    valid_sources = set()
    for nu in root.findall("nuclide"):
        name = nu.get("name")
        for src in nu.findall("source"):
            if src.get("particle") in ("photon",):
                valid_sources.add(name)
                break

    # Filter out targets without a valid source
    return {
        target: chans for target, chans in channels.items() if target in valid_sources
    }
