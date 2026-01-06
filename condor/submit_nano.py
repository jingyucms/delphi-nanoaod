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
user_output_base = Path("/eos/user/z/zhangj/DELPHI")

def build_patterns(nickname):
    cfg = config[nickname]
    
    if cfg["type"] == "data":
        patterns = [
            str(opendata_dir / "collision-data" / subdir / "**" / f"*{cfg['extension']}")
            for subdir in cfg["subdirs"]
        ]
        copy_dir = user_output_base / "collision_data" / cfg["version"] / cfg["energy"]
        TYPE = "DATA"
    elif cfg["type"] == "sim" and "subdirs" in cfg:
        # Handle hybrid case: sim type but with data-like search path
        patterns = []
        
        if "file_patterns" in cfg:
            # Handle multiple subdirs with specific file patterns
            for subdir_config in cfg["file_patterns"]:
                subdir = subdir_config["subdir"]
                file_pattern = subdir_config["pattern"]
                
                # Check if pattern contains numeric range like [241-252]
                import re
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
        if cfg["copy_energy"]:
            copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["copy_energy"] / cfg["stream"]
        else:
            copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["energy"] / cfg["stream"]
        TYPE = "MC"
    elif cfg["stream"] == "pythia8":
        pattern = f"/eos/user/z/zhangj/DELPHI/simulation/{cfg['version']}/{cfg['energy']}/pythia8_sdst/simana*sdst"
        patterns = [pattern]
        copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["energy"] / cfg["stream"]
        TYPE = "MC"
    elif cfg["stream"] == "pythia8_dire":
        pattern = f"/eos/user/z/zhangj/DELPHI/simulation/{cfg['version']}/{cfg['energy']}/pythia8_dire_sdst/simana*sdst"
        patterns = [pattern]
        copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["energy"] / cfg["stream"]
        TYPE = "MC"
    else:
        # Original sim type logic
        pattern = str(
            opendata_dir / "simulated-data" / cfg["origin"] / "**" / cfg["version"] /
            "**" / f"{cfg['stream']}_*{cfg['energy']}*{cfg['extension']}"
        )
        patterns = [pattern]
        if cfg["copy_energy"]:
            copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["copy_energy"] / cfg["stream"]
        else:
            copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["energy"] / cfg["stream"]
        TYPE = "MC"
    
    print("Search pattern:", patterns)
    return patterns, str(copy_dir), TYPE

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
    
    patterns, copy_dir, run_type = build_patterns(nickname)
    matched_files = find_matches(patterns)

    copy_dir = copy_dir

    print(f"Found {len(matched_files)} files for '{nickname}':")
    for f in matched_files:
        print(f)
    print(f"Copy directory: {copy_dir}")

    os.makedirs(copy_dir, exist_ok=True)    # Pause for manual check
    input("👉 Press Enter to start job submission...")
    
    for f in matched_files:
        fname = os.path.basename(f)
        output = f"nanoaod_{fname}.root"

        output_path = os.path.join(copy_dir, output)

        if file_should_be_skipped(output_path, min_bytes, max_days):
            print(f"⚠️ Skipping {output}")
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
        env["COPYDIR"] = copy_dir
        env["ISMC"] = run_type

        print(f, output, copy_dir)

        subprocess.run(["condor_submit", "condor/condor.sub"], env=env)
