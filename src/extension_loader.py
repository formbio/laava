#!/usr/bin/env python3
"""
LAAVA Extension Loader

CONFIGURATION-DRIVEN FUNCTION POINTER SYSTEM

This is NOT a traditional plugin system - it's a function pointer assignment system
that achieves zero runtime overhead by making decisions at import time.

How it works:
1. IMPORT TIME: Load config file and assign function pointers
2. RUNTIME: Direct function calls through pointers (zero overhead)

Example:
    # Import-time decision (happens once)
    cigar_func = get_extension('cigar_processor', default_function)
    
    # Runtime calls (zero overhead - direct function pointer invocation)
    result = cigar_func(data)  # Direct call, no dispatch overhead

Design principles:
- Zero runtime overhead (function pointer assignment at import time)
- Zero open source contamination (generic loading capability only)
- Configuration-driven (external YAML/JSON controls which functions get loaded)
- Backward compatible (graceful fallback to defaults when no config)
- Function pointer pattern (not plugin dispatch pattern)
"""

import os
import sys
import importlib
import logging
from typing import Dict, Any, Callable, Optional
from pathlib import Path

try:
    import yaml
    YAML_AVAILABLE = True
except ImportError:
    YAML_AVAILABLE = False

try:
    import json
    JSON_AVAILABLE = True
except ImportError:
    JSON_AVAILABLE = False


class ExtensionLoader:
    """
    Loads extensions based on configuration file.
    
    Supports both YAML and JSON configuration formats.
    """
    
    def __init__(self):
        self.extensions = {}
        self.config_loaded = False
        
    def load_config(self, config_path: Optional[str] = None) -> Dict[str, Any]:
        """
        Load extension configuration from file.
        
        Args:
            config_path: Path to configuration file. If None, uses environment variable.
            
        Returns:
            Dictionary of loaded extensions
        """
        if self.config_loaded:
            return self.extensions
            
        # Get config path from parameter or environment
        if config_path is None:
            config_path = os.environ.get('LAAVA_EXTENSIONS_CONFIG')
            
        if not config_path:
            logging.debug("No extension config specified")
            self.config_loaded = True
            return self.extensions
            
        config_file = Path(config_path)
        if not config_file.exists():
            logging.warning(f"Extension config file not found: {config_path}")
            self.config_loaded = True
            return self.extensions
            
        try:
            # Load configuration based on file extension
            if config_file.suffix.lower() in ['.yaml', '.yml']:
                if not YAML_AVAILABLE:
                    logging.error("YAML config specified but PyYAML not available")
                    self.config_loaded = True
                    return self.extensions
                    
                with open(config_file, 'r') as f:
                    config = yaml.safe_load(f)
                    
            elif config_file.suffix.lower() == '.json':
                if not JSON_AVAILABLE:
                    logging.error("JSON config specified but json module not available")
                    self.config_loaded = True
                    return self.extensions
                    
                with open(config_file, 'r') as f:
                    config = json.load(f)
                    
            else:
                logging.error(f"Unsupported config file format: {config_file.suffix}")
                self.config_loaded = True
                return self.extensions
                
            # Load extensions from config
            if 'extensions' in config:
                self._load_extensions_from_config(config['extensions'])
                logging.info(f"Loaded {len(self.extensions)} extensions from {config_path}")
            else:
                logging.warning(f"No 'extensions' section found in {config_path}")
                
        except Exception as e:
            logging.error(f"Error loading extension config {config_path}: {e}")
            
        self.config_loaded = True
        return self.extensions
    
    def _load_extensions_from_config(self, extensions_config: Dict[str, Any]):
        """Load extensions from configuration dictionary."""
        
        for extension_name, extension_spec in extensions_config.items():
            try:
                if isinstance(extension_spec, dict):
                    # Function replacement extension
                    if 'module' in extension_spec and 'function' in extension_spec:
                        func = self._load_function(extension_spec['module'], extension_spec['function'])
                        if func:
                            self.extensions[extension_name] = func
                            logging.debug(f"Loaded function extension: {extension_name}")
                    
                    # Class-based extension
                    elif 'module' in extension_spec and 'class' in extension_spec:
                        cls = self._load_class(extension_spec['module'], extension_spec['class'])
                        if cls:
                            # Instantiate with optional config
                            config = extension_spec.get('config', {})
                            instance = cls(config) if config else cls()
                            self.extensions[extension_name] = instance
                            logging.debug(f"Loaded class extension: {extension_name}")
                    
                    else:
                        logging.warning(f"Invalid extension spec for {extension_name}: missing module/function or module/class")
                        
                else:
                    logging.warning(f"Invalid extension spec for {extension_name}: expected dict, got {type(extension_spec)}")
                    
            except Exception as e:
                logging.error(f"Error loading extension {extension_name}: {e}")
    
    def _load_function(self, module_name: str, function_name: str) -> Optional[Callable]:
        """Load a function from a module."""
        try:
            module = importlib.import_module(module_name)
            if hasattr(module, function_name):
                func = getattr(module, function_name)
                if callable(func):
                    return func
                else:
                    logging.error(f"'{function_name}' in module '{module_name}' is not callable")
            else:
                logging.error(f"Function '{function_name}' not found in module '{module_name}'")
        except ImportError as e:
            logging.error(f"Could not import module '{module_name}': {e}")
        except Exception as e:
            logging.error(f"Error loading function '{function_name}' from '{module_name}': {e}")
        
        return None
    
    def _load_class(self, module_name: str, class_name: str) -> Optional[type]:
        """Load a class from a module."""
        try:
            module = importlib.import_module(module_name)
            if hasattr(module, class_name):
                cls = getattr(module, class_name)
                if isinstance(cls, type):
                    return cls
                else:
                    logging.error(f"'{class_name}' in module '{module_name}' is not a class")
            else:
                logging.error(f"Class '{class_name}' not found in module '{module_name}'")
        except ImportError as e:
            logging.error(f"Could not import module '{module_name}': {e}")
        except Exception as e:
            logging.error(f"Error loading class '{class_name}' from '{module_name}': {e}")
        
        return None
    
    def get_extension(self, name: str, default: Any = None) -> Any:
        """
        Get an extension by name.
        
        Args:
            name: Extension name
            default: Default value if extension not found
            
        Returns:
            Extension object or default value
        """
        return self.extensions.get(name, default)
    
    def has_extension(self, name: str) -> bool:
        """Check if an extension is loaded."""
        return name in self.extensions
    
    def list_extensions(self) -> Dict[str, str]:
        """List all loaded extensions with their types."""
        result = {}
        for name, ext in self.extensions.items():
            if callable(ext):
                result[name] = f"function: {ext.__name__}"
            else:
                result[name] = f"object: {type(ext).__name__}"
        return result


# Global extension loader instance
_extension_loader = None

def get_extension_loader() -> ExtensionLoader:
    """Get the global extension loader instance."""
    global _extension_loader
    if _extension_loader is None:
        _extension_loader = ExtensionLoader()
        _extension_loader.load_config()  # Load config on first access
    return _extension_loader

def get_extension(name: str, default: Any = None) -> Any:
    """
    Get an extension by name (FUNCTION POINTER ASSIGNMENT).
    
    This is the core of the zero-overhead function pointer system:
    - Returns optimized function if extension loaded
    - Returns default function if no extension found
    - Decision made once at import time for zero runtime overhead
    
    Args:
        name: Extension name
        default: Default function pointer to use as fallback
        
    Returns:
        Function pointer (either extension or default)
    """
    return get_extension_loader().get_extension(name, default)

def has_extension(name: str) -> bool:
    """Check if an extension is loaded (convenience function)."""
    return get_extension_loader().has_extension(name)
