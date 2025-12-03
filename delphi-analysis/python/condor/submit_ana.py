import os
import glob
import subprocess
from pathlib import Path
import sys
import yaml

def load_config(yaml_path="../../condor/sample_list.yaml"):
    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)

opendata_dir = Path("/eos/opendata/delphi")
user_output_base = Path("/eos/user/z/zhangj/DELPHI")

def build_patterns(nickname):
    cfg = config[nickname]
    
    if cfg["type"] == "data":
        copy_dir = user_output_base / "collision_data" / cfg["version"] / cfg["energy"]
    else:
        copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["copy_energy"] / cfg["stream"]

    return str(copy_dir)

isGen = False

#executable = "analysis_correlation.py"
#executable = "analysis_eec_lep2.py"
executable = "analysis_eec.py"
#executable = "analysis_trk.py"
#executable = "create_response_matrices.py"
#executable = "create_response_matrices_thrust.py"

nicknames = [
    "sh_kk2f4146qqpy_e91.25_c94_2l_c2",
#    "sh_kk2f4146qqardcy_e91.25_r94_2l_c2",
#    "Pythia8_94c",
#    "Pythia8_Dire_94c",
    "sh_kk2f4146qqpy_e91.25_c95_1l_d2",
#    "Pythia8_95d",
#    "Pythia8_Dire_95d"
    "short94_c2",
    "short95_d2"
#    "ALEPHMC"
#    "ALEPH"
]

#nicknames = [
    #"short94_c2"
#    "short95_d2"
#    "Pythia8_94c"
#]

#nicknames = [
    #"xs_kk2f4144tthl_e201.6_c99_1l_e1",
    #"xs_wphact211ncgg_e201.6_m80.4_c99_1l_e1",
    #"xs_wphact21nc4f_e201.6_m80.4_l99_1l_e1",
#    "xs_wphact24cc_e201.6_m80.4_c99_1l_e1",
    #"xs_kk2f4143qq_e201.6_l99_1l_e1",
#]

version = "v40"

config = load_config()

# Loop through each nickname
for nickname in nicknames:
    print(f"\n{'='*60}")
    print(f"Processing nickname: {nickname}")
    print(f"{'='*60}")
    
    if nickname == "ALEPHMC":
        PATTERN = "/eos/user/z/zhangj/ALEPH/SamplesLEP1/ALEPHMC/LEP1MC1994_recons_aftercut-0*.root"
        JobFlavour = "workday"
    elif nickname == "ALEPH":
        PATTERN = "/eos/user/z/zhangj/ALEPH/SamplesLEP1/ALEPH/LEP1Data1994P*part*"
        JobFlavour = "longlunch"
    elif nickname == "oldDELPHI":
        PATTERN  = "/eos/user/z/zhangj/DELPHI/simulation/1994_v2/qqps/qqps*.sdst.root"
        JobFlavour = "espresso"
    elif nickname == "Pythia8_94c":
        PATTERN = "/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/pythia8/nanoaod_simana_job_*.sdst.root"
        JobFlavour = "espresso"
    elif nickname == "Pythia8_Dire_94c":
        PATTERN = "/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/pythia8_dire/nanoaod_simana_job_*.sdst.root"
        JobFlavour = "espresso"
    elif nickname == "Pythia8_95d":
        PATTERN = "/eos/user/z/zhangj/DELPHI/simulation/v95d/91.25/pythia8/nanoaod_simana_job_*.sdst.root"
        JobFlavour = "espresso"
    elif nickname == "Pythia8_Dire_95d":
        PATTERN = "/eos/user/z/zhangj/DELPHI/simulation/v95d/91.25/pythia8_dire/nanoaod_simana_job_*.sdst.root"
        JobFlavour = "espresso"
    elif nickname == "sh_zgpy_b94_2l_c2" or nickname == "sh_qqps_k94_2l_c2":
        PATTERN = build_patterns(nickname)+"/nanoaod_*al.root"
        JobFlavour = "espresso"
    elif nickname == "short94_c2" or nickname == "short95_d2":
        PATTERN = build_patterns(nickname)+"/nanoaod_*al.root"
        JobFlavour = "longlunch"
    else:
        PATTERN = build_patterns(nickname)+"/nanoaod_*sdst.root"
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
        
        subprocess.run(["condor_submit", "condor/condor.sub"], env=env, check=True)
    
    print(f"Completed {len(file_list)} submissions for {nickname}")

print(f"\n{'='*60}")
print("All submissions completed!")
print(f"{'='*60}")
