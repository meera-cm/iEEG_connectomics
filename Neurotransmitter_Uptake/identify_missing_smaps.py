import pathlib
import os

alphaSeedsDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds'
alphaSmapDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds\HCP_MGH_32fold_groupconnectome (Horn 2017)'


alphaSeedsDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds'
alphaSmapDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\alpha_seeds\HCP_MGH_32fold_groupconnectome (Horn 2017)'

alphaSmapDir = pathlib.Path(alphaSmapDir)
parent_dir = pathlib.Path(alphaSmapDir).parent
for file in parent_dir.iterdir():
    if not file.is_file():
        continue
    basename = pathlib.Path(file).with_suffix("").name
    if basename == "ea_ui":
        continue
    new_file = alphaSmapDir / f"{basename}_struc_seed.nii"
    if not new_file.exists():
        print(basename, "missing")

betaSmapDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\beta_seeds\HCP_MGH_32fold_groupconnectome (Horn 2017)'

betaSmapDir = pathlib.Path(betaSmapDir)
parent_dir = pathlib.Path(betaSmapDir).parent
for file in parent_dir.iterdir():
    if not file.is_file():
        continue
    basename = pathlib.Path(file).with_suffix("").name
    if basename == "ea_ui":
        continue
    new_file = betaSmapDir / f"{basename}_struc_seed.nii"
    if not new_file.exists():
        print(basename, "missing")

thetaSmapDir = r'C:\Users\ICN_guest\Charité - Universitätsmedizin Berlin\Interventional Cognitive Neuromodulation - PROJECT ECoG Atlas Connectomics\data\original_data\seeds\theta_seeds\HCP_MGH_32fold_groupconnectome (Horn 2017)'

thetaSmapDir = pathlib.Path(thetaSmapDir)
parent_dir = pathlib.Path(thetaSmapDir).parent
for file in parent_dir.iterdir():
    if not file.is_file():
        continue
    basename = pathlib.Path(file).with_suffix("").name
    if basename == "ea_ui":
        continue
    new_file = thetaSmapDir / f"{basename}_struc_seed.nii"
    if not new_file.exists():
        print(basename, "missing")