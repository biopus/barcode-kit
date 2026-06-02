from __future__ import annotations

import importlib
import pkgutil

import barcode_kit


def test_each_barcode_kit_module_declares_resolvable_exports() -> None:
    package_names = [barcode_kit.__name__]
    package_names.extend(
        f"{barcode_kit.__name__}.{module.name}"
        for module in pkgutil.iter_modules(barcode_kit.__path__)
    )

    for package_name in package_names:
        module = importlib.import_module(package_name)
        exports = getattr(module, "__all__", None)

        assert exports is not None, f"{package_name} must define __all__"
        assert isinstance(exports, list), f"{package_name}.__all__ must be a list"
        assert all(isinstance(name, str) for name in exports)
        assert len(exports) == len(set(exports)), f"{package_name}.__all__ has duplicates"
        for name in exports:
            assert not name.startswith("_"), f"{package_name} exports private name {name}"
            assert hasattr(module, name), f"{package_name}.__all__ includes missing {name}"
