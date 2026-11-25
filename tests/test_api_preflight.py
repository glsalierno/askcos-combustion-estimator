"""Offline coverage for the ASKCOS reachability check."""

import sys
from unittest.mock import Mock

import pytest
import requests

import process_reactivity_oxygen as estimator


@pytest.mark.parametrize(
    "status,allow,expected",
    [
        (204, "", True),
        (405, "OPTIONS, POST", True),
        (405, "GET", False),
        (302, "", False),
        (403, "", False),
        (404, "", False),
        (503, "", False),
    ],
)
def test_api_status(monkeypatch, status, allow, expected):
    response = requests.Response()
    response.status_code = status
    response.headers["Allow"] = allow
    response._content_consumed = True
    methods = Mock(return_value=response)
    post = Mock(side_effect=AssertionError("preflight must not send a prediction"))
    monkeypatch.setattr(estimator.requests, "options", methods)
    monkeypatch.setattr(estimator.requests, "post", post)

    assert estimator.check_api("http://askcos.example/") is expected
    methods.assert_called_once_with(
        f"http://askcos.example{estimator.ASKCOS_ENDPOINT}",
        timeout=(3.05, 5),
        allow_redirects=False,
    )
    post.assert_not_called()


def test_check_api_failure_exits_without_reading_batch(monkeypatch, caplog):
    monkeypatch.setattr(sys, "argv", ["askcos-iterative", "--check-api"])
    monkeypatch.setattr(
        estimator.requests, "options", Mock(side_effect=requests.Timeout)
    )
    read = Mock(side_effect=AssertionError("preflight must not read the input"))
    monkeypatch.setattr(estimator.pd, "read_csv", read)

    with pytest.raises(SystemExit) as error:
        estimator.main()
    assert error.value.code == 1
    assert "Timeout" in caplog.text
    read.assert_not_called()
