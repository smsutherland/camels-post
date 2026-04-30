import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "dirs",
        type=Path,
        nargs="+",
        help="Where to look for simulations to post-process",
    )
    parser.add_argument(
        "--target",
        type=Path,
        help="Directory to put sbatch scripts",
        default=Path.cwd(),
    )
    parser.add_argument(
        "--cpus-per-task",
        "-c",
        type=int,
        default=16,
        help="How many cpus should each post-processing task use?",
    )
    parser.add_argument(
        "--ntasks",
        "-n",
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
    parser.add_argument(
        "--preset",
        help="Which preset post-processing configuration to use?",
        choices=["swimba", "gadget-nbody"],
        default="swimba",
    )
    args = parser.parse_args()
    dirs: set[Path] = set(d.resolve() for d in args.dirs)
    target: Path = args.target.resolve()
    parallel: int = args.cpus_per_task
    ntasks: int = min(args.ntasks, len(dirs))
    force: bool = args.force
    force_flag: str = " -f" if force else ""
    preset: str = args.preset

    target.mkdir(exist_ok=True)

    with open(Path(__file__).parent / "post_process.sh") as f:
        post_process_file = f.read()

    with open(Path(__file__).parent / "presets" / (preset + ".sh")) as f:
        preset = f.read()

    post_process_file = post_process_file.replace("PRESET_OPTIONS", preset)

    with open(target / "post_process.sh", "w") as f:
        f.write(post_process_file)

    post_process_path = target / "post_process.sh"

    disbatch_script = target / "tasks"
    with open(disbatch_script, "w") as f:
        f.write("#DISBATCH PREFIX cd \n")
        f.write(
            f"#DISBATCH SUFFIX ; (bash {post_process_path} {force_flag}) &{'>' if force else '>>'} post_processing.log\n"
        )

        for path in dirs:
            f.write(str(path) + "\n")

    sbatch_script = target / "job.sh"
    with open(sbatch_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --ntasks={ntasks}\n")
        f.write(f"#SBATCH --cpus-per-task={parallel}\n")
        f.write("#SBATCH --job-name=SWIMBA_postprocess\n")
        f.write("#SBATCH --partition=cmbas\n")
        f.write("module load disBatch\n")
        f.write("disBatch tasks\n")

    disbatch_check = target / "check_tasks"
    with open(disbatch_check, "w") as f:
        f.write("#DISBATCH PREFIX cd \n")
        f.write(
            f"#DISBATCH SUFFIX ; (bash {post_process_path} check) &>> post_processing.log\n"
        )

        for path in dirs:
            f.write(str(path) + "\n")
