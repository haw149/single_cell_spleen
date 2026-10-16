"""Download the GitHub release dataset into the app's local cache."""
import hashlib
import json
import os
from pathlib import Path
import threading
from urllib.error import HTTPError
from urllib.request import Request, HTTPRedirectHandler, build_opener
from urllib.parse import urlsplit

import h5py

REPOSITORY = "haw149/single_cell_spleen"
TAG = "data-v1"
ASSET = "clustered_filtered.h5ad"
_lock = threading.Lock()
_started = False
_status = ""


class SafeRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urlsplit(newurl).scheme != "https":
            raise ValueError("Dataset redirects must use HTTPS.")
        redirected = super().redirect_request(req, fp, code, msg, headers, newurl)
        if redirected is not None and urlsplit(req.full_url).netloc != urlsplit(newurl).netloc:
            redirected.remove_header("Authorization")
        return redirected


def download(root, token=None):
    target = Path(root) / "results_intermediate" / ASSET
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name("." + target.name + ".download")
    opener = build_opener(SafeRedirect())
    headers = {"Accept": "application/vnd.github+json",
               "X-GitHub-Api-Version": "2022-11-28", "User-Agent": "EDNA-Shiny"}
    if token:
        headers["Authorization"] = f"Bearer {token}"
    req = Request(f"https://api.github.com/repos/{REPOSITORY}/releases/tags/{TAG}", headers=headers)
    with opener.open(req, timeout=60) as response:
        release = json.load(response)
    asset = next((a for a in release.get("assets", []) if a["name"] == ASSET), None)
    if asset is None:
        raise ValueError(f"Release {TAG} has no {ASSET} attachment.")
    headers["Accept"] = "application/octet-stream"
    req = Request(f"https://api.github.com/repos/{REPOSITORY}/releases/assets/{int(asset['id'])}", headers=headers)
    digest = hashlib.sha256()
    try:
        with opener.open(req, timeout=120) as response, temporary.open("wb") as output:
            while chunk := response.read(1024 * 1024):
                output.write(chunk)
                digest.update(chunk)
        if temporary.stat().st_size != asset["size"]:
            raise ValueError("Dataset download was incomplete. Retry the download.")
        expected = asset.get("digest")
        if expected and expected.startswith("sha256:") and expected[7:] != digest.hexdigest():
            raise ValueError("Dataset checksum did not match the release attachment.")
        with h5py.File(temporary, "r") as handle:
            if not all(key in handle for key in ("X", "obs", "var")):
                raise ValueError("The attachment is not an AnnData file with X, obs and var.")
        temporary.replace(target)
    finally:
        temporary.unlink(missing_ok=True)


def ensure_dataset(root, retry=False):
    global _started, _status
    if (Path(root) / "results_intermediate" / ASSET).is_file():
        return
    with _lock:
        if _started and not (retry and not _status.startswith("Downloading")):
            return
        token = os.environ.get("EDNA_GITHUB_TOKEN")
        _started = True
        _status = "Downloading the dataset from GitHub. This may take several minutes; the page refreshes automatically."

    def worker():
        global _status
        try:
            download(root, token)
            _status = "Dataset downloaded. Click Refresh project if it has not appeared yet."
        except HTTPError as exc:
            _status = (f"GitHub dataset download failed (HTTP {exc.code}). Check that release data-v1 contains "
                       "clustered_filtered.h5ad. Private repositories also require EDNA_GITHUB_TOKEN with Contents: Read access. "
                       "Click Refresh project to retry.")
        except Exception:
            # Never expose request URLs, authorization headers, or token values in the UI/log.
            _status = "Dataset download failed validation or could not finish. Check the release attachment and click Refresh project to retry."
    threading.Thread(target=worker, daemon=True).start()


def status():
    return _status or "No AnnData file found. Add a saved .h5ad dataset to the deployment."
