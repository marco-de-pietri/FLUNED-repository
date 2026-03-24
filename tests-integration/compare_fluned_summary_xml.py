#!/usr/bin/env python3

import argparse
import math
import re
import sys
import xml.etree.ElementTree as ET


NUMERIC_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")


def normalize(text):
    return "" if text is None else text.strip()


def parse_number(text):
    if not text or not NUMERIC_RE.fullmatch(text):
        return None

    try:
        return float(text)
    except ValueError:
        return None


def compare_scalar(path, expected, actual, rel_tol, abs_tol, differences):
    expected_text = normalize(expected)
    actual_text = normalize(actual)

    expected_number = parse_number(expected_text)
    actual_number = parse_number(actual_text)

    if expected_number is not None and actual_number is not None:
        if not math.isclose(
            expected_number,
            actual_number,
            rel_tol=rel_tol,
            abs_tol=abs_tol,
        ):
            abs_diff = abs(expected_number - actual_number)
            scale = max(abs(expected_number), abs(actual_number), abs_tol)
            rel_diff = abs_diff / scale
            differences.append(
                f"{path}: numeric mismatch "
                f"expected={expected_text} actual={actual_text} "
                f"abs_diff={abs_diff:.6g} rel_diff={rel_diff:.6g}"
            )
        return

    if expected_text != actual_text:
        differences.append(
            f"{path}: value mismatch expected={expected_text!r} actual={actual_text!r}"
        )


def build_child_keys(children):
    keys = []
    anonymous_counts = {}

    for child in children:
        child_id = child.attrib.get("id")
        if child_id is not None:
            keys.append((child.tag, child_id))
            continue

        anonymous_counts[child.tag] = anonymous_counts.get(child.tag, 0) + 1
        keys.append((child.tag, anonymous_counts[child.tag]))

    return keys


def format_child_path(parent_path, key):
    tag, identifier = key
    if isinstance(identifier, str):
        return f"{parent_path}/{tag}[@id='{identifier}']"
    return f"{parent_path}/{tag}[{identifier}]"


def compare_elements(path, expected, actual, rel_tol, abs_tol, differences):
    if expected.tag != actual.tag:
        differences.append(
            f"{path}: tag mismatch expected={expected.tag!r} actual={actual.tag!r}"
        )
        return

    expected_attr_keys = set(expected.attrib)
    actual_attr_keys = set(actual.attrib)

    for attr_name in sorted(expected_attr_keys - actual_attr_keys):
        differences.append(f"{path}/@{attr_name}: missing attribute in actual XML")

    for attr_name in sorted(actual_attr_keys - expected_attr_keys):
        differences.append(f"{path}/@{attr_name}: unexpected attribute in actual XML")

    for attr_name in sorted(expected_attr_keys & actual_attr_keys):
        compare_scalar(
            f"{path}/@{attr_name}",
            expected.attrib[attr_name],
            actual.attrib[attr_name],
            rel_tol,
            abs_tol,
            differences,
        )

    expected_text = normalize(expected.text)
    actual_text = normalize(actual.text)
    if expected_text or actual_text:
        compare_scalar(path, expected_text, actual_text, rel_tol, abs_tol, differences)

    expected_children = list(expected)
    actual_children = list(actual)
    expected_keys = build_child_keys(expected_children)
    actual_keys = build_child_keys(actual_children)

    expected_map = {}
    for key, child in zip(expected_keys, expected_children):
        expected_map[key] = child

    actual_map = {}
    for key, child in zip(actual_keys, actual_children):
        actual_map[key] = child

    processed = set()

    for key in expected_keys:
        child_path = format_child_path(path, key)
        if key not in actual_map:
            differences.append(f"{child_path}: missing element in actual XML")
            continue

        compare_elements(
            child_path,
            expected_map[key],
            actual_map[key],
            rel_tol,
            abs_tol,
            differences,
        )
        processed.add(key)

    for key in actual_keys:
        if key in processed or key in expected_map:
            continue
        differences.append(
            f"{format_child_path(path, key)}: unexpected element in actual XML"
        )


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare two fluned_summary.xml files semantically."
    )
    parser.add_argument("expected_xml", help="Reference XML file")
    parser.add_argument("actual_xml", help="Produced XML file")
    parser.add_argument(
        "--rel-tol",
        type=float,
        default=1e-5,
        help="Relative tolerance for numeric values",
    )
    parser.add_argument(
        "--abs-tol",
        type=float,
        default=1e-12,
        help="Absolute tolerance for numeric values",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    try:
        expected_root = ET.parse(args.expected_xml).getroot()
        actual_root = ET.parse(args.actual_xml).getroot()
    except ET.ParseError as exc:
        print(f"XML parse error: {exc}", file=sys.stderr)
        return 2
    except OSError as exc:
        print(f"File error: {exc}", file=sys.stderr)
        return 2

    differences = []
    compare_elements(
        f"/{expected_root.tag}",
        expected_root,
        actual_root,
        args.rel_tol,
        args.abs_tol,
        differences,
    )

    if not differences:
        return 0

    print(
        "NOT PASSED - fluned_summary.xml differences below "
        f"(rel_tol={args.rel_tol}, abs_tol={args.abs_tol})"
    )
    for difference in differences:
        print(difference)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
