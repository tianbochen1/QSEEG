import json
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
import analyse_eeg
import eeg_features
from core import SPECTRAL_FLOOR, spectral_features


class SpectralEligibilityTests(unittest.TestCase):
    def test_one_empty_quantile_excludes_the_window_for_every_response(self):
        good = np.repeat(np.arange(20, dtype=float), 5)
        bad = np.r_[np.zeros(99), 100.]
        mixed = spectral_features(np.stack([good, bad]), [.1, .5, .9], [.06, .2], .02)
        reference = spectral_features(good, [.1, .5, .9], [.06, .2], .02)
        for values, expected in zip(mixed, reference):
            np.testing.assert_allclose(values[0], expected)
            self.assertTrue(np.isfinite(expected).all())
            self.assertTrue(np.isnan(values[1]).all())

    def test_screen_is_specific_to_window_length_and_channel(self):
        rng = np.random.default_rng(409)
        signals = rng.normal(size=(2, 2, 6000))
        signals[0, 0, :3000] = np.r_[np.zeros(2999), 100.]
        metadata = dict(subject=1, night=1, age=55, female=0, treatment=0,
                        study='SC', record='SC-test')
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            inputs = root / 'segments'
            outputs = root / 'features'
            inputs.mkdir()
            outputs.mkdir()
            path = inputs / 'SC-test.npz'
            np.savez(path, signals=signals, rows=[[2, 0, 1.], [2, 60, 2.]],
                     quality=np.ones((2, 2), dtype=bool), metadata=json.dumps(metadata))
            with patch.object(eeg_features, 'OUT', outputs), patch.object(analyse_eeg, 'FEATURES', outputs):
                eeg_features.one(path)
                for variant, channel, expected in [('main', 0, 2), ('short', 0, 1), ('short', 1, 2)]:
                    frame, quantile, ordinary, shape, candidates = analyse_eeg.load('SC', 2, channel, variant)
                    self.assertEqual(candidates, 2)
                    self.assertEqual(len(frame), expected)
                    for values in (quantile, ordinary, shape):
                        self.assertEqual(len(values), expected)
                        self.assertTrue(np.isfinite(values).all())
                with np.load(outputs / path.name) as feature:
                    legacy = {k: feature[k] for k in feature if not k.endswith('_valid')}
                legacy['short_q'][0, 0] = np.float32(np.log(SPECTRAL_FLOOR))
                np.savez(outputs / path.name, **legacy)
                with self.assertRaisesRegex(ValueError, 'rebuild features'):
                    analyse_eeg.load('SC', 2, 0, 'short')
                eeg_features.one(path)
                self.assertEqual(len(analyse_eeg.load('SC', 2, 0, 'short')[0]), 1)
                with np.load(outputs / path.name) as feature:
                    marked = {k: feature[k] for k in feature}
                marked['main_q'][0, 0, 0] = np.float32(np.log(SPECTRAL_FLOOR))
                np.savez(outputs / path.name, **marked)
                self.assertEqual(len(analyse_eeg.load('SC', 2, 0, 'main')[0]), 2)

    def test_incomplete_pairs_leave_the_fit_and_inclusion_counts_together(self):
        rng = np.random.default_rng(87)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            features = root / 'features'
            results = root / 'results'
            features.mkdir()
            results.mkdir()
            for subject in range(10):
                for night in [1, 2]:
                    treatment = int(night == 1 + subject % 2)
                    record = f'ST-{subject:02d}-{night}'
                    metadata = dict(subject=subject, night=night, age=35+subject, female=0,
                                    treatment=treatment, study='ST', record=record)
                    q = rng.normal(size=(1, 2, 9, 3))
                    ordinary = rng.normal(size=(1, 2, 3))
                    np.savez(features / (record+'.npz'), metadata=json.dumps(metadata),
                             rows=[[2, 0, rng.uniform()]], quality=[[True, True]],
                             short_valid=[[True, True]], short_q=q, short_o=ordinary, short_s=ordinary)
            with patch.object(analyse_eeg, 'FEATURES', features), patch.object(analyse_eeg, 'OUT', results):
                analyse_eeg.run('ST', 2, 0, 'short')
                meta = results / 'ST_2_0_short.json'
                first = json.loads(meta.read_text())
                self.assertEqual((first['N'], first['segments']), (10, 20))
                excluded = features / 'ST-00-2.npz'
                with np.load(excluded) as archive:
                    changed = {k: archive[k] for k in archive}
                changed['short_valid'][0, 0] = False
                np.savez(excluded, **changed)
                analyse_eeg.run('ST', 2, 0, 'short')
                after = json.loads(meta.read_text())
                self.assertNotEqual(first['input_sha256'], after['input_sha256'])
                self.assertEqual((after['N'], after['segments'], after['candidate_segments']), (9, 18, 20))
                included = pd.read_csv(results / 'ST_2_0_short_included.csv')
                self.assertNotIn(0, included.subject.to_numpy())
                self.assertEqual(included.groupby('subject').treatment.nunique().tolist(), [2] * 9)


if __name__ == '__main__':
    unittest.main()
