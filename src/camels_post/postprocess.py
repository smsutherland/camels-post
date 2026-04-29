import shutil
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dirs",
        type=Path,
        nargs="*",
        help="Where to start looking for simulations to post-process",
    )
    parser.add_argument(
        "--target",
        type=Path,
        help="Directory to put sbatch scripts",
        default=Path.cwd(),
    )
    parser.add_argument(
        "--parallel",
        type=int,
        default=16,
        help="How much parallelization should each post-processing task use?",
    )
    parser.add_argument(
        "--ntasks",
        type=int,
        default=32,
        help="How many post-processing tasks should run at once?",
    )
    parser.add_argument(
        "--force",
        "-f",
        action="store_true",
        help="Make post-processing jobs start over",
    )
    args = parser.parse_args()
    dirs: list[Path] = args.dirs
    target: Path = args.target
    parallel: int = args.parallel
    ntasks: int = min(args.ntasks, len(dirs))
    force: bool = args.force
    force_flag: str = " -f" if force else ""

    post_script = Path(__file__).parent / "./data/post_processing.sh"

    target.mkdir(exist_ok=True)

    disbatch_script = target / "tasks"
    with open(disbatch_script, "w") as f:
        f.write("#DISBATCH PREFIX cd \n")
        f.write(
            f"#DISBATCH SUFFIX ; (bash post_processing.sh {force_flag}) &{'>' if force else '>>'} post_processing.log\n"
        )

        paths = set()
        for path in dirs:
            shutil.copy(post_script, path)
            path = str(path.resolve())
            if path not in paths:
                f.write(path + "\n")
                paths.add(path)

    sbatch_script = target / "job.sh"
    with open(sbatch_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --ntasks={ntasks}\n")
        f.write(f"#SBATCH --cpus-per-task={parallel}\n")
        f.write("#SBATCH --job-name=SWIMBA_postprocess\n")
        f.write("#SBATCH --partition=cmbas\n")
        f.write("module load disBatch\n")
        f.write("disBatch tasks\n")
