class BarcodeKitError(Exception):
    """Base exception for user-facing toolkit failures."""


class ConfigError(BarcodeKitError):
    """Raised when configuration is missing or invalid."""


class GenBankError(BarcodeKitError):
    """Raised when GenBank retrieval or parsing fails."""


class GenBankSearchError(GenBankError):
    """Raised when GenBank search requests fail."""


class GenBankFetchError(GenBankError):
    """Raised when GenBank record fetch or response handling fails."""

    def __init__(
        self,
        message: str,
        *,
        error_type: str = "unexpected_response",
        retryable: bool = False,
    ):
        super().__init__(message)
        self.error_type = error_type
        self.retryable = retryable


class TaxonomyError(BarcodeKitError):
    """Raised when taxonomy standardization fails."""


class BuildError(BarcodeKitError):
    """Raised when a dataset cannot be built."""


class PhylogenyError(BarcodeKitError):
    """Raised when alignment, trimming, or tree reconstruction fails."""
