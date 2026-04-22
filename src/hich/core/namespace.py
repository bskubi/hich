import tiledb as tdb
import numpy as np
import polars as pl
import pysam
from io import StringIO
import importlib
from typing import Protocol, Any, cast, Generic, TypeVar, Tuple

T = TypeVar("T")

def namespace_from_path(path, module_name: str, protocol: type[T]) -> T:
    """
        Import a namespace from a path and cast to a Protocol for IDE visibility.

        Args:
            path: Path to python script to import at runtime as module.
            module_name: Name of imported module.
            protocol: Protocol defining required attributes of imported module.

        Example:
            ```python
            # config.py
            def hello(): 
                print("Hello, world!")

            # main.py
            from typing import Protocol
            
            class HelloConfig(Protocol):
                def hello(self) -> None:
                    "Prints 'Hello, world!' to terminal."
                    ...
            
            # Hovering over 'hello' variable provides IDE hints from HelloConfig
            hello = namespace_from_path("config.py", "hello", HelloConfig)
            hello.hello()
            ```
    """
    # Create spec
    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load {path}")
    
    # Create module based on spec
    module = importlib.util.module_from_spec(spec)

    # Execute module
    spec.loader.exec_module(module)

    # Cast module to Protocol for IDE visibility
    module_T = cast(T, module)

    return module_T