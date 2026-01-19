import os
import glob
import subprocess
import time
from pathlib import Path
import sys
import yaml

def load_config(yaml_path="../../condor/sample_list.yaml"):
    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)

opendata_dir = Path("/eos/opendata/delphi")
user_output_base = Path("/eos/user/z/zhangj/DELPHI")  # Keep old output location

# Load config early so it's available to all functions
config = None

def load_config_if_needed():
    """Lazy load config when first needed"""
    global config
    if config is None:
        config = load_config()

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

def build_input_pattern(nickname):
    """
    Build input file pattern from new EOS structure.
    Returns the pattern to find input files.
    """
    load_config_if_needed()
    cfg = config[nickname]
    base = Path("/eos/experiment/eealliance/Samples/DELPHI")
    
    # Extract version info
    year, short_version = extract_version_info(cfg["version"])
    
    # Use copy_energy if exists, otherwise use energy
    energy = cfg.get("copy_energy") or cfg["energy"]
    
    if cfg["type"] == "data":
        # Data: read from TPCNtuple (for *al.root files) or DelphiNanoAOD (for analysis)
        # Most analysis uses the _ttree.root files from DelphiNanoAOD
        input_dir = base / year / energy / "Data" / short_version / "TPCNtuple" / "251219"
    else:
        # MC: read from TPCNtuple directory
        stream = cfg["stream"]
        input_dir = base / year / energy / "MC" / short_version / "TPCNtuple" / stream / "251219"
    
    return str(input_dir)

def build_patterns(nickname):
    """
    Legacy function for compatibility - returns old output directory structure.
    Output directory stays the same (old EOS location).
    """
    load_config_if_needed()
    cfg = config[nickname]
    
    if cfg["type"] == "data":
        copy_dir = user_output_base / "collision_data" / cfg["version"] / cfg["energy"]
    else:
        energy = cfg.get("copy_energy") or cfg["energy"]
        copy_dir = user_output_base / "simulation" / cfg["version"] / energy / cfg["stream"]

    return str(copy_dir)

def is_mc_sample(nickname):
    """
    Determine if a sample is MC or data from config.
    Returns True for MC, False for data.
    """
    # Handle special legacy cases not in config
    if nickname == "ALEPHMC":
        return True
    elif nickname == "ALEPH":
        return False
    
    # For all other samples, check config
    load_config_if_needed()
    if nickname not in config:
        raise ValueError(f"Nickname '{nickname}' not found in config. Please add it to sample_list.yaml")
    
    return config[nickname]["type"] == "sim"

isGen = False

# Executable to run
#executable = "analysis_correlation.py"
#executable = "analysis_eec_lep2.py"
#executable = "analysis_eec.py"
#executable = "analysis_trk.py"
#executable = "create_response_matrices.py"
#executable = "create_response_matrices_thrust.py"
executable = "analysis_thrust.py"

nicknames = [
    "sh_kk2f4146qqpy_e91.25_c94_2l_c2",
    "sh_kk2f4146qqardcy_e91.25_r94_2l_c2",
    "sh_pythia8_94c",
    "sh_pythia8_dire_94c",
    "sh_kk2f4146qqpy_e91.25_c95_1l_d2",
    "sh_pythia8_95d",
    "sh_pythia8_dire_95d",
    "short94_c2",
    "short95_d2"
#    "ALEPHMC"
#    "ALEPH"
]

version = "v48"

# Loop through each nickname
for nickname in nicknames:
    print(f"\n{'='*60}")
    print(f"Processing nickname: {nickname}")
    print(f"{'='*60}")
    
    # Determine if this is MC
    is_mc = is_mc_sample(nickname)
    print(f"Sample type: {'MC' if is_mc else 'Data'}")
    
    # Handle special legacy cases (ALEPH only) that don't use new structure
    if nickname == "ALEPHMC":
        PATTERN = "/eos/user/z/zhangj/ALEPH/SamplesLEP1/ALEPHMC/LEP1MC1994_recons_aftercut-0*.root"
        JobFlavour = "workday"
    elif nickname == "ALEPH":
        PATTERN = "/eos/user/z/zhangj/ALEPH/SamplesLEP1/ALEPH/LEP1Data1994P*part*"
        JobFlavour = "longlunch"
    # All DELPHI samples use new structure with same pattern
    else:
        input_dir = build_input_pattern(nickname)
        PATTERN = f"{input_dir}/nanoaod_*.root"
        JobFlavour = "longlunch"

    print(f"Pattern: {PATTERN}")
    
    file_list = sorted(glob.glob(PATTERN))
    print(f"Found {len(file_list)} files for {nickname}")
    
    for file_path in file_list:
        file_name = os.path.basename(file_path)

        print(f"  Submitting: {file_name}")

        env = os.environ.copy()
        env["p"] = file_path
        env["f"] = file_name.replace("root", version+".root")
        env["JOB_FLAVOUR"] = JobFlavour
        env["EXE"] = executable
        env["ana"] = executable.replace("analysis_", "").replace(".py", "").replace("create", "")
        env["ana"] += '_'+nickname
        if isGen:
            env["ana"]+="gen"
        
        # Add --is-mc flag for MC samples
        if is_mc:
            env["EXTRA_ARGS"] = "--is-mc"
        else:
            env["EXTRA_ARGS"] = ""
        
        #subprocess.run(["condor_submit", "condor/condor.sub"], env=env, check=True)
        # Retry logic for transient HTCondor errors
        max_retries = 3
        for attempt in range(max_retries):
            result = subprocess.run(["condor_submit", "condor/condor.sub"], env=env)
            if result.returncode == 0:
                break
            if attempt < max_retries - 1:
                print(f"  ⚠️ HTCondor error, retry {attempt+1}/{max_retries}...")
                time.sleep(2)
            else:
                print(f"  ❌ Failed to submit {file_name} after {max_retries} attempts, skipping...")
    
    print(f"Completed {len(file_list)} submissions for {nickname}")

print(f"\n{'='*60}")
print("All submissions completed!")
print(f"{'='*60}")