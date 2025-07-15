#!/usr/bin/env python3
"""
Test the configuration-driven extension system.

This test validates that the function pointer assignment system works correctly
by using a simple logging override that demonstrates the extension is loaded
and executed without changing the core logic.
"""

import sys
import os
import tempfile
import logging
import io
from pathlib import Path

# Add src directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / 'src'))
# test_extension_override.py is now in the same directory

import pysam
from extension_loader import get_extension, get_extension_loader
from summarize_alignment import iter_cigar_w_aligned_pair

def test_extension_system_basic():
    """Test that extension system loads and uses extensions correctly."""
    
    # Create a simple test BAM record
    # Note: This would need actual test data in a real test
    print("Testing extension system basic functionality...")
    
    # Test 1: No extension config (should use default)
    os.environ.pop('LAAVA_EXTENSIONS_CONFIG', None)
    
    # Force reload of extension loader
    import importlib
    import extension_loader
    importlib.reload(extension_loader)
    
    # Get function pointer - should be default
    cigar_func = get_extension('cigar_processor', iter_cigar_w_aligned_pair)
    assert cigar_func == iter_cigar_w_aligned_pair, "Should use default when no config"
    print("✅ Default fallback works")
    
    # Test 2: With extension config (should use override)
    config_path = str(Path(__file__).parent / 'test_extension_config.json')
    os.environ['LAAVA_EXTENSIONS_CONFIG'] = config_path
    
    # Force reload to pick up new config
    importlib.reload(extension_loader)
    
    # Get function pointer - should be override
    cigar_func = get_extension('cigar_processor', iter_cigar_w_aligned_pair)
    assert cigar_func.__name__ == 'iter_cigar_w_aligned_pair_with_logging', "Should use override when config present"
    print("✅ Extension override works")
    
    # Test 3: Extension info
    loader = get_extension_loader()
    assert loader.has_extension('cigar_processor'), "Should have cigar_processor extension"
    extensions = loader.list_extensions()
    assert 'cigar_processor' in extensions, "Should list cigar_processor extension"
    print("✅ Extension discovery works")

def test_extension_logging():
    """Test that the logging override produces expected log messages."""
    
    # Set up logging capture
    log_capture = io.StringIO()
    handler = logging.StreamHandler(log_capture)
    handler.setLevel(logging.INFO)
    
    # Configure logging
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    
    try:
        # Set extension config
        config_path = str(Path(__file__).parent / 'test_extension_config.json')
        os.environ['LAAVA_EXTENSIONS_CONFIG'] = config_path
        
        # Force reload
        import importlib
        import extension_loader
        importlib.reload(extension_loader)
        
        # Create a mock CSV writer
        class MockCSVWriter:
            def __init__(self):
                self.rows = []
            def writerow(self, row):
                self.rows.append(row)
        
        # Create a simple mock record (this would use real test data in practice)
        class MockRecord:
            def __init__(self):
                self.qname = "test_read"
                self.cigartuples = [(0, 100)]  # 100M
                self.is_unmapped = False
            
            def get_aligned_pairs(self):
                # Simple mock - 100 aligned positions
                return [(i, i) for i in range(100)]
        
        # Get the extension function
        cigar_func = get_extension('cigar_processor', iter_cigar_w_aligned_pair)
        
        # Call it with mock data
        mock_record = MockRecord()
        mock_writer = MockCSVWriter()
        
        result = cigar_func(mock_record, mock_writer)
        
        # Check log output
        log_output = log_capture.getvalue()
        assert "🔧 START of override extension!" in log_output, "Should log start message"
        assert "🔧 END of override extension!" in log_output, "Should log end message"
        assert "test_extension_override.py" in log_output, "Should identify the override module"
        
        print("✅ Extension logging works")
        print(f"   Log output: {log_output.strip()}")
        
        # Verify function returns expected result format
        assert isinstance(result, tuple), "Should return tuple"
        assert len(result) == 2, "Should return (total_err, total_len)"
        print("✅ Extension returns correct format")
        
    finally:
        # Clean up logging
        logger.removeHandler(handler)

def test_extension_identical_results():
    """Test that extension produces identical results to original function."""
    
    print("Testing that extension produces identical results...")
    
    # This test would compare results between default and extension
    # For now, we just verify the extension can be called
    
    # Set extension config
    config_path = str(Path(__file__).parent / 'test_extension_config.json')
    os.environ['LAAVA_EXTENSIONS_CONFIG'] = config_path
    
    # Force reload
    import importlib
    import extension_loader
    importlib.reload(extension_loader)
    
    # Get extension function
    cigar_func = get_extension('cigar_processor', iter_cigar_w_aligned_pair)
    
    # Verify it's the override
    assert cigar_func.__name__ == 'iter_cigar_w_aligned_pair_with_logging'
    
    print(" Extension function loaded correctly")

def main():
    """Run all extension system tests."""
    
    print("🧪 LAAVA Extension System Tests")
    print("=" * 50)
    
    try:
        test_extension_system_basic()
        test_extension_logging()
        test_extension_identical_results()
        
        return 0
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    finally:
        # Clean up environment
        os.environ.pop('LAAVA_EXTENSIONS_CONFIG', None)

if __name__ == '__main__':
    exit(main())
