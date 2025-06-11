import linecache
import tracemalloc
import sys
import gzip


def check_int(check_str):
    """
    Check if a string can be converted to an integer
    """
    try:
        int(check_str)
        return True
    except ValueError:
        return False


def mod(vec1):
    sum = 0
    for el in vec1:
        sum += el**2

    return pow(sum, 0.5)


def check_float(str):
    """
    this function checks if a string is a float
    """
    try:
        float(str)
        return True
    except ValueError:
        return False


def open_utf8_or_gzip(file_path):
    """
    this function tries to open a file with utf-8 encoding, if it fails it
    tries to open it with gzip. It returns the file object
    """

    try:
        with open(file_path, "rb") as f:
            magic = f.read(2)
    except OSError:
        print("Couldn't open Volume V file")
        sys.exit(1)

    if magic == b"\x1f\x8b":
        # The file is gzip-compressed
        try:
            with gzip.open(file_path, "rt", encoding="utf-8") as inpFile:
                data = inpFile.read()
        except Exception as e:
            print(f"Error reading gzip file: {e}")
            sys.exit(1)

    else:
        # The file is a regular text file
        try:
            with open(file_path, "r", encoding="utf-8") as inpFile:
                data = inpFile.read()
        except Exception as e:
            print(f"Error reading text file: {e}")
            sys.exit(1)

    return data


def display_top(snapshot, key_type="lineno", limit=10):
    snapshot = snapshot.filter_traces(
        (
            tracemalloc.Filter(False, "<frozen importlib._bootstrap>"),
            tracemalloc.Filter(False, "<unknown>"),
        )
    )
    top_stats = snapshot.statistics(key_type)

    print("Top %s lines" % limit)
    for index, stat in enumerate(top_stats[:limit], 1):
        frame = stat.traceback[0]
        print(
            "#%s: %s:%s: %.1f KiB"
            % (index, frame.filename, frame.lineno, stat.size / 1024)
        )
        line = linecache.getline(frame.filename, frame.lineno).strip()
        if line:
            print("    %s" % line)

    other = top_stats[limit:]
    if other:
        size = sum(stat.size for stat in other)
        print("%s other: %.1f KiB" % (len(other), size / 1024))
    total = sum(stat.size for stat in top_stats)
    print("Total allocated size: %.1f KiB" % (total / 1024))

    return
