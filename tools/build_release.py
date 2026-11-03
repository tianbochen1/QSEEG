"""Assemble and verify the complete computational bundle from pinned inputs."""
from pathlib import Path
from urllib.request import urlretrieve
import argparse
import hashlib
import json
import tempfile
import zipfile

ROOT = Path(__file__).resolve().parents[1]
BASE_URL = 'https://github.com/tianbochen1/QSEEG/releases/download/v1.0.2/QSEEG_code.zip'
BASE_SHA256 = '38bc82d95869c735edc1f12f845fae71f59386e9e0d59f64c0ab3617c44e656d'


def build(baseline, destination):
    assert hashlib.sha256(baseline.read_bytes()).hexdigest() == BASE_SHA256
    manifest = json.loads((ROOT / 'docs/code_manifest.json').read_text())
    destination.mkdir(parents=True, exist_ok=True)
    archive = destination / 'QSEEG_code.zip'
    with zipfile.ZipFile(baseline) as original, zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as output:
        for name, expected in sorted(manifest.items()):
            candidate = ROOT / ('docs/REPRODUCIBILITY.md' if name == 'README.md' else name)
            data = candidate.read_bytes() if candidate.is_file() else original.read('QSEEG_code/' + name)
            actual = hashlib.sha256(data).hexdigest()
            assert actual == expected, f'Content mismatch: {name}'
            info = zipfile.ZipInfo('QSEEG_code/' + name, date_time=(2026, 11, 3, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.create_system = 3
            info.external_attr = 0o100644 << 16
            output.writestr(info, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    with zipfile.ZipFile(archive) as output:
        assert output.testzip() is None
        assert len(output.namelist()) == len(manifest)
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    (destination / 'SHA256SUMS.txt').write_text(digest + '  QSEEG_code.zip\n')
    print(f'Verified {len(manifest)} files; archive SHA256 {digest}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--baseline', type=Path)
    parser.add_argument('--output', type=Path, default=ROOT / 'dist')
    args = parser.parse_args()
    if args.baseline is None:
        with tempfile.TemporaryDirectory() as directory:
            baseline = Path(directory) / 'baseline.zip'
            urlretrieve(BASE_URL, baseline)
            build(baseline, args.output)
    else:
        build(args.baseline, args.output)
