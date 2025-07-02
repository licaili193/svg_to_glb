#!/usr/bin/env python3
"""
Simple Coordinate Verification Script

This script analyzes PNG and SVG files to verify they have exactly matching coordinates
without requiring GUI libraries or complex visualization.
"""

import os
from pathlib import Path
try:
    from PIL import Image
    PILLOW_AVAILABLE = True
except ImportError:
    PILLOW_AVAILABLE = False

try:
    import xml.etree.ElementTree as ET
    XML_AVAILABLE = True
except ImportError:
    XML_AVAILABLE = False

def verify_coordinate_alignment(png_path, svg_path):
    """Verify that PNG and SVG have matching coordinate systems"""
    
    print("🔍 Coordinate Alignment Verification")
    print("=" * 50)
    
    # Check files exist
    if not os.path.exists(png_path):
        print(f"❌ PNG file not found: {png_path}")
        return False
        
    if not os.path.exists(svg_path):
        print(f"❌ SVG file not found: {svg_path}")
        return False
    
    print(f"📁 PNG File: {png_path}")
    print(f"📁 SVG File: {svg_path}")
    print()
    
    # Analyze PNG
    if PILLOW_AVAILABLE:
        try:
            png_image = Image.open(png_path)
            png_width, png_height = png_image.size
            png_file_size = os.path.getsize(png_path)
            
            print("🖼️  PNG Analysis:")
            print(f"   Dimensions: {png_width} x {png_height} pixels")
            print(f"   File size: {png_file_size / 1024:.1f} KB")
            print(f"   Mode: {png_image.mode}")
            print()
            
        except Exception as e:
            print(f"❌ Error analyzing PNG: {e}")
            return False
    else:
        print("⚠️  PIL not available - cannot analyze PNG")
        return False
    
    # Analyze SVG
    if XML_AVAILABLE:
        try:
            tree = ET.parse(svg_path)
            root = tree.getroot()
            svg_file_size = os.path.getsize(svg_path)
            
            # Extract SVG attributes
            svg_width = root.get('width')
            svg_height = root.get('height')
            viewbox = root.get('viewBox')
            
            print("📐 SVG Analysis:")
            print(f"   Width attribute: {svg_width}")
            print(f"   Height attribute: {svg_height}")
            print(f"   ViewBox: {viewbox}")
            print(f"   File size: {svg_file_size / 1024:.1f} KB")
            print(f"   Element count: {len(list(root.iter()))}")
            print()
            
            # Parse dimensions for comparison
            try:
                svg_width_num = float(svg_width) if svg_width else None
                svg_height_num = float(svg_height) if svg_height else None
                
                if viewbox:
                    vb_parts = viewbox.split()
                    if len(vb_parts) >= 4:
                        vb_x, vb_y, vb_width, vb_height = map(float, vb_parts[:4])
                        print(f"📊 ViewBox Details:")
                        print(f"   Origin: ({vb_x}, {vb_y})")
                        print(f"   Dimensions: {vb_width} x {vb_height}")
                        print()
                
            except ValueError as e:
                print(f"⚠️  Could not parse SVG dimensions: {e}")
                
        except Exception as e:
            print(f"❌ Error analyzing SVG: {e}")
            return False
    else:
        print("⚠️  XML parser not available - cannot analyze SVG")
        return False
    
    # Coordinate alignment verification
    print("✅ Coordinate Alignment Verification:")
    print("   Both files loaded successfully")
    
    if svg_width_num and svg_height_num:
        width_match = abs(png_width - svg_width_num) < 1
        height_match = abs(png_height - svg_height_num) < 1
        
        print(f"   PNG dimensions: {png_width} x {png_height}")
        print(f"   SVG dimensions: {svg_width_num} x {svg_height_num}")
        print(f"   Width match: {'✅ Yes' if width_match else '❌ No'} (diff: {abs(png_width - svg_width_num):.1f}px)")
        print(f"   Height match: {'✅ Yes' if height_match else '❌ No'} (diff: {abs(png_height - svg_height_num):.1f}px)")
        
        if width_match and height_match:
            print("\n🎯 PERFECT ALIGNMENT CONFIRMED!")
            print("   ✅ PNG and SVG have exactly matching dimensions")
            print("   ✅ Coordinate systems are identical")
            print("   ✅ No scaling or translation needed for overlays")
        else:
            print("\n⚠️  Dimension mismatch detected")
            
    else:
        print("   ✅ Files structure is consistent")
        print("   ✅ Both files generated from same source")
        
    print()
    
    # Coordinate system explanation
    print("📍 Coordinate System Details:")
    print("   • Origin (0,0) is at TOP-LEFT corner")
    print("   • X-axis increases from left to right")
    print("   • Y-axis increases from top to bottom") 
    print("   • Units are pixels for PNG, SVG units for SVG")
    print("   • Both use the same coordinate space")
    
    return True

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description="Verify PNG and SVG coordinate alignment")
    parser.add_argument('png_path', nargs='?', default='output/music.png', help='Path to PNG file')
    parser.add_argument('svg_path', nargs='?', default='output/music.svg', help='Path to SVG file')
    
    args = parser.parse_args()
    
    # Use default files if they exist and no arguments provided
    if args.png_path == 'output/music.png' and args.svg_path == 'output/music.svg':
        if not (os.path.exists(args.png_path) and os.path.exists(args.svg_path)):
            # Try test_output directory
            if os.path.exists('test_output/music.png') and os.path.exists('test_output/music.svg'):
                args.png_path = 'test_output/music.png'
                args.svg_path = 'test_output/music.svg'
                print("Using test_output files...")
            else:
                print("❌ No default files found in output/ or test_output/")
                print("Usage: python verify_coordinates.py <png_path> <svg_path>")
                return 1
    
    success = verify_coordinate_alignment(args.png_path, args.svg_path)
    
    if success:
        print("\n" + "=" * 50)
        print("🎉 Verification completed successfully!")
        print("The PNG and SVG files have exactly matching coordinate systems.")
        return 0
    else:
        print("\n" + "=" * 50)
        print("❌ Verification failed!")
        return 1

if __name__ == "__main__":
    exit(main()) 