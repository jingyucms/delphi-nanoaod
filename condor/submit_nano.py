import os
import shutil
import glob
import sys
import subprocess
import multiprocessing
from pathlib import Path
import time
import re
import yaml

opendata_dir = Path("/eos/opendata/delphi")

def extract_version_info(version_str):
    """
    Extract year and short version from version string.
    Example: "v94c" -> year="1994", short_version="94c"
    """
    # Remove 'v' prefix
    short_version = version_str.lstrip('v')
    
    # Extract year prefix (first 2 digits)
    year_prefix = short_version[:2]
    
    # Map to full year
    year = "19" + year_prefix
    
    return year, short_version

def build_output_paths(cfg, date="251219"):
    """
    Build output paths for both TPCNtuple and DelphiNanoAOD.
    Returns: (tpcntuple_dir, delphinanoaod_dir, TYPE)
    """
    base = Path("/eos/experiment/eealliance/Samples/DELPHI")
    
    # Extract version info
    year, short_version = extract_version_info(cfg["version"])
    
    # Use copy_energy if exists, otherwise use energy
    energy = cfg.get("copy_energy") or cfg["energy"]
    
    if cfg["type"] == "data":
        # Data path (no stream subfolder)
        tpcntuple_dir = base / year / energy / "Data" / short_version / "TPCNtuple" / date
        delphinanoaod_dir = base / year / energy / "Data" / short_version / "DelphiNanoAOD" / date
        TYPE = "DATA"
    else:
        # MC path (with stream subfolder)
        stream = cfg["stream"]
        tpcntuple_dir = base / year / energy / "MC" / short_version / "TPCNtuple" / stream / date
        delphinanoaod_dir = base / year / energy / "MC" / short_version / "DelphiNanoAOD" / stream / date
        TYPE = "MC"
    
    return str(tpcntuple_dir), str(delphinanoaod_dir), TYPE

def build_patterns(nickname):
    cfg = config[nickname]
    
    if cfg["type"] == "data":
        patterns = [
            str(opendata_dir / "collision-data" / subdir / "**" / f"*{cfg['extension']}")
            for subdir in cfg["subdirs"]
        ]
    elif cfg["type"] == "sim" and "subdirs" in cfg:
        # Handle hybrid case: sim type but with data-like search path
        patterns = []
        
        if "file_patterns" in cfg:
            # Handle multiple subdirs with specific file patterns
            for subdir_config in cfg["file_patterns"]:
                subdir = subdir_config["subdir"]
                file_pattern = subdir_config["pattern"]
                
                # Check if pattern contains numeric range like [241-252]
                range_match = re.search(r'\[(\d+)-(\d+)\]', file_pattern)
                if range_match:
                    start, end = map(int, range_match.groups())
                    base_pattern = file_pattern.replace(range_match.group(0), '{}')
                    for num in range(start, end + 1):
                        pattern = str(opendata_dir / "collision-data" / subdir / base_pattern.format(num))
                        patterns.append(pattern)
                else:
                    pattern = str(opendata_dir / "collision-data" / subdir / file_pattern)
                    patterns.append(pattern)
        elif "file_pattern" in cfg:
            # Handle single file pattern for all subdirs
            for subdir in cfg["subdirs"]:
                pattern = str(opendata_dir / "collision-data" / subdir / cfg["file_pattern"])
                patterns.append(pattern)
        else:
            # Default to wildcard pattern
            for subdir in cfg["subdirs"]:
                pattern = str(opendata_dir / "collision-data" / subdir / f"*{cfg['extension']}")
                patterns.append(pattern)
    elif cfg["stream"] in ["pythia8", "pythia8_dire"]:
        # Pythia8 files are in the new EOS structure under SDST
        year, short_version = extract_version_info(cfg["version"])
        energy = cfg.get("copy_energy") or cfg["energy"]
        base = Path("/eos/experiment/eealliance/Samples/DELPHI")
        pattern = str(base / year / energy / "MC" / short_version / "SDST" / cfg["stream"] / "251219" / "simana*sdst")
        patterns = [pattern]
    else:
        # Original sim type logic
        pattern = str(
            opendata_dir / "simulated-data" / cfg["origin"] / "**" / cfg["version"] /
            "**" / f"{cfg['stream']}_*{cfg['energy']}*{cfg['extension']}"
        )
        patterns = [pattern]
    
    # Build output paths using new function
    tpcntuple_dir, delphinanoaod_dir, TYPE = build_output_paths(cfg)
    
    print("Search pattern:", patterns)
    return patterns, tpcntuple_dir, delphinanoaod_dir, TYPE

def find_matches(patterns):
    matches = []
    for pattern in patterns:
        matches.extend(glob.glob(pattern, recursive=True))
    return matches

def file_should_be_skipped(path, MIN_SIZE_BYTES, MAX_AGE_DAYS):
    if not os.path.exists(path):
        return False
    stat = os.stat(path)
    size_ok = stat.st_size > MIN_SIZE_BYTES
    age_ok = (time.time() - stat.st_mtime) < (MAX_AGE_DAYS * 86400)
    return size_ok and age_ok

def total_jobs(user):
    try:
        output = subprocess.check_output(
            ["condor_q", user, "-total"],
            stderr=subprocess.DEVNULL,
            text=True
        )
        for line in output.splitlines():
            if line.startswith("Total for query:"):
                # Example line: "Total for query: 200 jobs; 0 completed; 0 removed; 200 idle; 0 running; 0 held; 0 suspended"
                numbers = re.findall(r"\b\d+\b", line)
                # Only count idle + running jobs (as per previous context)
                if len(numbers) >= 5:
                    return int(numbers[3]) + int(numbers[4])
        return 0
    except Exception as e:
        print(f"Error checking condor queue: {e}")
        return 0

def load_config(yaml_path="condor/sample_list.yaml"):
    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)

if __name__ == "__main__":

    min_bytes = 100_000      # 100 KB
    max_days = 3

    MAX_QUEUE = 200
    USER = os.environ["USER"]

    config = load_config()
    
    #nickname = "short94_c2"
    #nickname = "sh_qqps_e91.25_c94_2l_c2"
    nickname = "sh_kk2f4146qqpy_e91.25_c94_2l_c2"
    #nickname = "sh_kk2f4146qqpydcy_e91.25_c94_2l_c2"
    #nickname = "sh_kk2f4146qqardcy_e91.25_r94_2l_c2"
    #nickname = "sh_apacic105_e91.25_w94_2l_c2"
    #nickname = "sh_zgpy_b94_2l_c2"
    #nickname = "sh_qqps_k94_2l_c2"
    #nickname = "sh_pythia8_94c"
    #nickname = "sh_pythia8_dire_94c"
    #nickname = "sh_kk2f4146tthl_e91.25_c94_2l_c2"

    #nickname = "short95_d2"
    #nickname = "sh_kk2f4146qqpy_e91.25_c95_1l_d2"
    #nickname = "sh_pythia8_95d"
    #nickname = "sh_pythia8_dire_95d"

    #nickname = "xsdst99_e192_e1"
    #nickname = "xsdst99_e196_e1"
    #nickname = "xsdst99_e200_e1"
    #nickname = "xsdst99_e202_e1"
    #nickname = "xs_kk2f4143qq_e191.6_r99_1l_e1"
    #nickname = "xs_kk2f4143qq_e195.5_l99_1l_e1"
    #nickname = "xs_kk2f4143qq_e199.5_c99_1l_e1"
    #nickname = "xs_kk2f4143qq_e201.6_l99_1l_e1"
    #nickname = "xs_wphact24cc_e191.6_m80.4_c99_1l_e1"
    #nickname = "xs_wphact24cc_e195.5_m80.4_c99_1l_e1"
    #nickname = "xs_wphact24cc_e199.5_m80.4_c99_1l_e1"
    #nickname = "xs_wphact24cc_e201.6_m80.4_c99_1l_e1"
    #nickname = "xs_wphact21nc4f_e191.6_m80.4_l99_1l_e1"
    #nickname = "xs_wphact21nc4f_e195.5_m80.4_l99_1l_e1"
    #nickname = "xs_wphact21nc4f_e199.5_m80.4_l99_1l_e1"
    #nickname = "xs_wphact21nc4f_e201.6_m80.4_l99_1l_e1"
    #nickname = "xs_wphact211ncgg_e191.6_m80.4_c99_1l_e1"
    #nickname = "xs_wphact211ncgg_e195.5_m80.4_c99_1l_e1"
    #nickname = "xs_wphact211ncgg_e199.5_m80.4_c99_1l_e1"
    #nickname = "xs_wphact211ncgg_e201.6_m80.4_c99_1l_e1"
    #nickname = "xs_kk2f4144tthl_e191.6_c99_1l_e1"
    #nickname = "xs_kk2f4144tthl_e195.5_c99_1l_e1"
    #nickname = "xs_kk2f4144tthl_e199.5_c99_1l_e1"
    #nickname = "xs_kk2f4144tthl_e201.6_c99_1l_e1"
    
    patterns, tpcntuple_dir, delphinanoaod_dir, run_type = build_patterns(nickname)
    matched_files = find_matches(patterns)

    print(f"Found {len(matched_files)} files for '{nickname}':")
    for f in matched_files:
        print(f)
    print(f"TPCNtuple directory: {tpcntuple_dir}")
    print(f"DelphiNanoAOD directory: {delphinanoaod_dir}")

    # Create both output directories
    os.makedirs(tpcntuple_dir, exist_ok=True)
    os.makedirs(delphinanoaod_dir, exist_ok=True)

    input("👉 Press Enter to start job submission...")
    
    for f in matched_files:
        fname = os.path.basename(f)
        output = f"nanoaod_{fname}.root"

        # Check if either output file exists and should be skipped
        tpcntuple_path = os.path.join(tpcntuple_dir, output)
        ttree_output = output.replace(".root", "_ttree.root")
        delphinanoaod_path = os.path.join(delphinanoaod_dir, ttree_output)
        
        if (file_should_be_skipped(tpcntuple_path, min_bytes, max_days) and 
            file_should_be_skipped(delphinanoaod_path, min_bytes, max_days)):
            print(f"⚠️ Skipping {output} (both outputs exist)")
            continue

        # To avoid ran out of disk space on lxplus
        while total_jobs(USER) >= MAX_QUEUE:
            print(f"⏳ Too many jobs in Condor queue for {USER}. Waiting 1 minutes...")
            time.sleep(60)
            
        print(f"✅ Submitting job for {output}")

        # Create a copy of current environment and add your variables
        env = os.environ.copy()
        env["IN"] = f
        env["OUT"] = output
        env["TPCNTUPLE_DIR"] = tpcntuple_dir
        env["DELPHINANOAOD_DIR"] = delphinanoaod_dir
        env["ISMC"] = run_type

        print(f, output)
        print(f"  -> TPCNtuple: {tpcntuple_dir}")
        print(f"  -> DelphiNanoAOD: {delphinanoaod_dir}")

        subprocess.run(["condor_submit", "condor/condor.sub"], env=env)
