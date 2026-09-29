"""
Base connector framework with 3-tier fallback architecture:
1. Live HTTP API with exponential backoff retries.
2. Local persistent disk cache.
3. Embedded curated JSON fixtures for offline / air-gapped reproducibility.
"""

from __future__ import annotations

import json
import os
import time
from pathlib import Path
from typing import Any, Dict, Optional

import requests


class ConnectorError(Exception):
    """Base exception for all multi-omic connector errors."""
    pass


class ConnectorNotFoundError(ConnectorError):
    """Raised when an accession or identifier is not found (HTTP 404)."""
    pass


class ConnectorNetworkError(ConnectorError):
    """Raised when network connectivity fails and no offline fallback is available."""
    pass


class BaseConnector:
    """Resilient connector base class with caching and fixture fallback."""

    def __init__(
        self,
        cache_dir: Optional[str] = None,
        timeout: float = 5.0,
        retries: int = 3,
        backoff_factor: float = 0.5,
    ) -> None:
        self.timeout = timeout
        self.retries = retries
        self.backoff_factor = backoff_factor

        if cache_dir is not None:
            self.cache_dir = Path(cache_dir)
        else:
            default_cache = Path.home() / ".cache" / "insilico_peptide"
            self.cache_dir = Path(os.getenv("INSILICO_CACHE_DIR", str(default_cache)))

        try:
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass

        self.fixtures_dir = Path(__file__).parent / "fixtures"

    def _get_cache_file(self, key: str) -> Path:
        """Derive cache file path from key."""
        safe_key = "".join(c if c.isalnum() or c in ("-", "_") else "_" for c in key)
        return self.cache_dir / f"{safe_key}.json"

    def load_from_cache(self, key: str) -> Optional[Dict[str, Any]]:
        """Load cached JSON data if file exists."""
        cache_file = self._get_cache_file(key)
        if cache_file.exists():
            try:
                with open(cache_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError):
                return None
        return None

    def save_to_cache(self, key: str, data: Dict[str, Any]) -> None:
        """Persist JSON data to disk cache."""
        try:
            cache_file = self._get_cache_file(key)
            with open(cache_file, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)
        except OSError:
            pass

    def load_fixture(self, filename: str) -> Optional[Dict[str, Any]]:
        """Load embedded curated fixture."""
        fixture_file = self.fixtures_dir / filename
        if fixture_file.exists():
            try:
                with open(fixture_file, "r", encoding="utf-8") as f:
                    return json.load(f)
            except (json.JSONDecodeError, OSError):
                return None
        return None

    def fetch_json(
        self,
        url: str,
        params: Optional[Dict[str, Any]] = None,
        cache_key: Optional[str] = None,
        fixture_name: Optional[str] = None,
        offline: bool = False,
    ) -> Dict[str, Any]:
        """
        Execute 3-tier resilient JSON retrieval.

        Args:
            url: API endpoint URL.
            params: Optional query parameters.
            cache_key: Cache key (defaults to sanitized URL if None).
            fixture_name: Embedded fixture filename to fallback to.
            offline: Force offline mode without making network requests.

        Returns:
            Decoded JSON dictionary.

        Raises:
            ConnectorNotFoundError: If endpoint returned 404 and no fixture available.
            ConnectorNetworkError: If network request failed and no cache/fixture exists.
        """
        key = cache_key or url

        # Tier 2: Check disk cache
        cached = self.load_from_cache(key)
        if cached is not None:
            return cached

        # If offline forced, attempt fixture fallback immediately
        if offline:
            if fixture_name:
                fixture_data = self.load_fixture(fixture_name)
                if fixture_data is not None:
                    return fixture_data
            raise ConnectorNotFoundError(
                f"Offline mode active: No local cache or fixture available for key '{key}'"
            )

        # Tier 1: Live HTTP API with exponential backoff
        last_error: Optional[Exception] = None
        for attempt in range(self.retries):
            try:
                headers = {"Accept": "application/json", "User-Agent": "insilico-peptide-toolkit/0.1.0"}
                resp = requests.get(url, params=params, headers=headers, timeout=self.timeout)

                if resp.status_code == 200:
                    data = resp.json()
                    self.save_to_cache(key, data)
                    return data
                elif resp.status_code == 404:
                    # If fixture exists, fallback
                    if fixture_name:
                        fixture_data = self.load_fixture(fixture_name)
                        if fixture_data is not None:
                            return fixture_data
                    raise ConnectorNotFoundError(f"Resource not found at {url} (HTTP 404)")
                else:
                    resp.raise_for_status()

            except requests.exceptions.RequestException as e:
                last_error = e
                # Exponential backoff
                if attempt < self.retries - 1:
                    time.sleep(self.backoff_factor * (2**attempt))

        # Network failed across all retries: Fallback to Tier 3 (Fixture)
        if fixture_name:
            fixture_data = self.load_fixture(fixture_name)
            if fixture_data is not None:
                return fixture_data

        raise ConnectorNetworkError(
            f"Failed to fetch data from {url} after {self.retries} attempts: {last_error}"
        )
