import hashlib
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch
from urllib.request import Request

import h5py

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import cloud_data


class CloudDataTests(unittest.TestCase):
    def test_redirect_does_not_forward_token_to_storage_host(self):
        request = Request('https://api.github.com/repos/a/b', headers={'Authorization': 'Bearer test'})
        redirected = cloud_data.SafeRedirect().redirect_request(
            request, None, 302, '', {}, 'https://release-assets.githubusercontent.com/file')
        self.assertIsNone(redirected.get_header('Authorization'))
        with self.assertRaises(ValueError):
            cloud_data.SafeRedirect().redirect_request(request, None, 302, '', {}, 'http://example.com/file')

    def test_download_validates_and_publishes_complete_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fixture = root / 'fixture.h5ad'
            with h5py.File(fixture, 'w') as handle:
                for key in ['X', 'obs', 'var']:
                    handle.create_group(key)
            payload = fixture.read_bytes()
            release = {'assets': [{'name': cloud_data.ASSET, 'id': 123, 'size': len(payload),
                                   'digest': 'sha256:' + hashlib.sha256(payload).hexdigest()}]}
            with patch.object(cloud_data, 'build_opener') as build:
                build.return_value.open.side_effect = [io.BytesIO(json.dumps(release).encode()), io.BytesIO(payload)]
                cloud_data.download(root)
                for call in build.return_value.open.call_args_list:
                    self.assertIsNone(call.args[0].get_header('Authorization'))
            self.assertEqual((root/'results_intermediate'/cloud_data.ASSET).read_bytes(), payload)
            self.assertFalse(list((root/'results_intermediate').glob('*.download')))

    def test_incomplete_download_is_not_exposed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            release = {'assets': [{'name': cloud_data.ASSET, 'id': 123, 'size': 100}]}
            with patch.object(cloud_data, 'build_opener') as build:
                build.return_value.open.side_effect = [io.BytesIO(json.dumps(release).encode()), io.BytesIO(b'partial')]
                with self.assertRaises(ValueError):
                    cloud_data.download(root, 'test')
            self.assertFalse((root/'results_intermediate'/cloud_data.ASSET).exists())
            self.assertEqual(list((root/'results_intermediate').iterdir()), [])
