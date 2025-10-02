
import argparse
import subprocess

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='\n\n\033[1mGet S3 recursive bucket listing\033[0m\n\n',)
    parser.add_argument('s3link', type=str,
                        help='s3 bucket to query')
    parser.add_argument('--nameFilter', type=str, default=None,
                    help='Only include paths containing this substring')
    parser.add_argument('--excludeName', type=str, default=None,
                    help='Exclude paths containing this substring')
    parser.add_argument('--fileEndsWith', type=str, default=None,
                    help='Only include files ending with this string, e.g., ".h5"')
    parser.add_argument('--maxLevels', type=int, default=None,
                    help='Maximum recursion depth for directories; None = recurse fully')
    # Flat flag
    parser.add_argument('--flat', action='store_true',
                    help='Print files as a flat list of full paths instead of a tabbed hierarchy')
    args = parser.parse_args()
    printKeywords, searchKeywords = {}, {}
    for arg in {'nameFilter', 'fileEndsWith', 'excludeName'}:
        printKeywords[arg] = getattr(args, arg)
    #
    for arg in {'maxLevels'}:
        searchKeywords[arg] = getattr(args, arg)
    return args.s3link, searchKeywords, printKeywords, args.flat
    


def list_s3_tree(s3_uri, profile=None, maxLevels=None, level=0):
    """
    Recursively list S3 contents into a nested dictionary using `aws s3 ls`,
    starting from a given prefix in the bucket.

    Parameters
    ----------
    s3_uri : str
        S3 URI, e.g., 's3://my-bucket/xxx/yyy/zzz/'
    profile : str, optional
        AWS profile name.
    maxLevels : int, optional
        Maximum recursion depth. If None, recurse fully.
    level : int
        Current recursion depth (internal).

    Returns
    -------
    dict
        Structure like:
        {
          'files': [...],
          'directories': {
             'subdir/': { ... }
          }
        }
    """
    cmd_base = ["aws", "s3", "ls"]
    if profile:
        cmd_base += ["--profile", profile]

    tree = {"files": [], "directories": {}}

    try:
        result = subprocess.check_output(cmd_base + [s3_uri], text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error listing {s3_uri}: {e}")
        return tree

    for line in result.splitlines():
        if level == 2:
            print(line)
        line = line.strip()
        if line.startswith("PRE "):
            prefix = line.split()[-1]
            if maxLevels is None or level < maxLevels - 1:
                subtree = list_s3_tree(
                    s3_uri + prefix, profile=profile,
                    maxLevels=maxLevels, level=level+1
                )
                tree["directories"][prefix] = subtree
            else:
                tree["directories"][prefix] = {"files": [], "directories": {}}
        else:
            parts = line.split()
            if len(parts) >= 4:
                filename = parts[-1]
                tree["files"].append(s3_uri + filename)

    return tree


def print_s3_tree(tree, name="", level=0, 
                  nameFilter=None, fileEndsWith=None, excludeName=None):
    """
    Pretty-print the nested S3 dictionary with indentation.
    Supports filters:
      - nameFilter: include only names containing this substring (applies to dirs and files)
      - fileEndsWith: include only files ending with this string
      - excludeName: skip names containing this substring (applies to dirs and files)

    Parameters
    ----------
    tree : dict
        Nested dict returned by list_s3_tree.
    name : str
        Name of the current directory (optional).
    level : int
        Indentation level.
    """
    indent = "\t" * level

    # Check if this directory should be printed
    if name:
        if excludeName and excludeName in name:
            show_dir = False
        elif nameFilter is None or nameFilter in name:
            print(f"{indent}{name}/")
            show_dir = True
        else:
            show_dir = False
    else:
        show_dir = True  # top-level bucket

    # Print files
    for f in tree.get("files", []):
        fname = f.split("/")[-1]
        if excludeName and excludeName in fname:
            continue
        if nameFilter and nameFilter not in fname:
            continue
        if fileEndsWith and not fname.endswith(fileEndsWith):
            continue
        print(f"{indent}\t{fname}")

    # Recurse into subdirectories
    for dirname, subtree in tree.get("directories", {}).items():
        if excludeName and excludeName in dirname:
            continue
        if nameFilter and nameFilter not in dirname:
            continue
        print_s3_tree(subtree, dirname.rstrip("/"), level+1,
                      nameFilter=nameFilter, fileEndsWith=fileEndsWith, excludeName=excludeName)

def print_s3_files_flat(tree, nameFilter=None, fileEndsWith=None, excludeName=None):
    """
    Print full paths of all files in the nested S3 dictionary (flat list),
    with optional filters.

    Parameters
    ----------
    tree : dict
        Nested dict returned by list_s3_tree.
    nameFilter : str or None
        Include only names containing this substring (applies to dirs and files)
    fileEndsWith : str or None
        Include only files ending with this string
    excludeName : str or None
        Skip names containing this substring (applies to dirs and files)
    """
    # Print files at this level
    for f in tree.get("files", []):
        fname = f.split("/")[-1]
        if excludeName and excludeName in fname:
            continue
        if nameFilter and nameFilter not in fname:
            continue
        if fileEndsWith and not fname.endswith(fileEndsWith):
            continue
        print(f)

    # Recurse into subdirectories
    for dirname, subtree in tree.get("directories", {}).items():
        if excludeName and excludeName in dirname:
            continue
        if nameFilter and nameFilter not in dirname:
            continue
        print_s3_files_flat(subtree, nameFilter=nameFilter,
                            fileEndsWith=fileEndsWith, excludeName=excludeName)

def run():
    
    initialLink, searchKeywords, printKeywords, flat = parseCommandLine()
    s3files = list_s3_tree(initialLink, **searchKeywords)
    print(initialLink)
    if flat:
        print_s3_files_flat(s3files, **printKeywords)
    else:
        print_s3_tree(s3files, **printKeywords)

if __name__ == "__main__":
    run()