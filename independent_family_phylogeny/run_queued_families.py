#!/usr/bin/env python3
"""Wait for the original family batch, then extend the same run. Plan by default."""
import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import run_phylogeny

ROOT = Path(__file__).resolve().parent


def command(queue):
    options = queue["parameters"]
    argv = ["--resume", "--run-name", queue["run_name"], "--alignment", options["alignment"],
            "--trim", options["trim"], "--mafft-threads", str(options["mafft_threads"]),
            "--iqtree-threads", options["iqtree_threads"], "--seed", str(options["seed"])]
    for panel in options["excluded_panels"]:
        argv += ["--exclude-panel", panel]
    for family in queue["all_families"]:
        argv += ["--family", family]
    return argv


def readiness(queue):
    folder = ROOT / "results" / queue["run_name"]
    missing = [name for name in queue["wait_for_families"] if not (folder / name / "tree.complete.json").exists()]
    path = folder / "run_status.json"
    try:
        status = json.loads(path.read_text()) if path.exists() else None
    except json.JSONDecodeError:
        status = None  # Original runner can be finishing its status write.
    if missing:
        if status and status.get("status") == "failed":
            raise ValueError("The original batch reports failure. Recover its unfinished families before starting additions.")
        return False, "Waiting for original family trees: " + ", ".join(missing)
    if not status or not set(queue["wait_for_families"]) <= set(status.get("requested_families", [])):
        return False, "Original tree receipts are present; waiting for the original runner to publish its final status."
    return True, "Original family batch finished; completed stages will be verified and reused."


def verify(queue):
    if run_phylogeny.digest(ROOT / "library/library_manifest.json") != queue["library_manifest_sha256"]:
        raise ValueError("Prepared library changed after queue creation; refresh the queued extension before execution")
    for name in queue["wait_for_families"]:
        config = ROOT / "results" / queue["run_name"] / name / "run_config.json"
        if config.exists():
            values = json.loads(config.read_text())
            if any(values[key] != expected for key, expected in queue["parameters"].items()):
                raise ValueError(f"Queue parameters differ from the saved original run: {name}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execute", action="store_true", help="Start the extension only after the original batch finishes")
    parser.add_argument("--wait", action="store_true", help="Wait without launching analytical tools while the original batch is active")
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--queue", type=Path, default=ROOT / "pending_families.json")
    args = parser.parse_args(argv)
    try:
        if args.poll_seconds < 1:
            raise ValueError("--poll-seconds must be positive")
        queue = json.loads(args.queue.read_text())
        verify(queue)
        print("Queued additions: " + ", ".join(queue["added_families"]), flush=True)
        if not args.execute:
            _, message = readiness(queue)
            print(message, flush=True)
            return run_phylogeny.main(command(queue))
        if os.name == "nt":
            raise ValueError("Execute this continuation in the original WSL/Linux environment")
        import fcntl
        with (ROOT / "pending_families.lock").open("a") as lock:
            try:
                fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise ValueError("Another queued continuation is already running") from error
            last_message = None
            while True:
                ready, message = readiness(queue)
                if message != last_message:
                    print(message, flush=True)
                    last_message = message
                if ready:
                    break
                if not args.wait:
                    raise ValueError("Original batch is not finished; use --execute --wait to defer the extension safely")
                time.sleep(args.poll_seconds)
            verify(queue)
            result = run_phylogeny.main(["--execute"] + command(queue))
            if result == 0:
                run_phylogeny.dump(ROOT / "pending_families.complete.json", {
                    "completed_utc": datetime.now(timezone.utc).isoformat(), "run_name": queue["run_name"],
                    "added_families": queue["added_families"],
                    "queue_sha256": run_phylogeny.digest(args.queue)})
            return result
    except (OSError, ValueError, KeyError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
