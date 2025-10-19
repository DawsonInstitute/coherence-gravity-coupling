#!/usr/bin/env python3
import argparse
import csv
import hashlib
from pathlib import Path
from datetime import datetime, timezone


def sha256_of(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open('rb') as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', default='data_manifest.csv')
    p.add_argument('--roots', nargs='+', default=['results', 'papers/figures'])
    args = p.parse_args()

    rows = []
    for root in args.roots:
        root_path = Path(root)
        if not root_path.exists():
            continue
        for path in root_path.rglob('*'):
            if path.is_file():
                stat = path.stat()
                rows.append(
                    {
                        'path': str(path),
                        'sha256': sha256_of(path),
                        'size_bytes': stat.st_size,
                        'created_utc': datetime.fromtimestamp(stat.st_mtime, tz=timezone.utc).isoformat(),
                        'description': '',
                    }
                )

    with open(args.output, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['path', 'sha256', 'size_bytes', 'created_utc', 'description'])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} entries to {args.output}")


if __name__ == '__main__':
    main()
