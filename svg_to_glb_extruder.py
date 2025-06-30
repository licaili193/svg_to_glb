#!/usr/bin/env python3
"""
SVG to GLB Extruder

A comprehensive tool for extruding SVG files (especially music scores) into 3D GLB meshes.
Handles complex challenges like:
- Pure line/path strokes appearing in 3D mesh
- Self-intersecting regions (like treble clef)
- Shapes with holes (sharp symbols, half notes)

Usage:
    python svg_to_glb_extruder.py input.svg [--output output.glb] [--extrude-height 5.0] [--stroke-width 0.5]
"""

import argparse
import sys
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

try:
    import trimesh
    import trimesh.path
    from trimesh.path import Path2D
    import trimesh.creation
except ImportError:
    print("Error: trimesh not found. Install with: pip install trimesh")
    sys.exit(1)

try:
    import cairosvg
    from PIL import Image
    import io
except ImportError:
    print("Warning: cairosvg or PIL not found. PNG rendering will be disabled.")
    print("Install with: pip install cairosvg pillow")
    cairosvg = None

try:
    import shapely
    from shapely.geometry import Polygon, MultiPolygon
    from shapely.ops import unary_union
    from shapely import buffer, simplify
except ImportError:
    print("Error: shapely not found. Install with: pip install shapely")
    sys.exit(1)

try:
    import svgpathtools
    from svgpathtools import svg2paths, Path as SVGPath
except ImportError:
    print("Warning: svgpathtools not found. SVG parsing will be limited.")
    print("Install with: pip install svgpathtools")
    svgpathtools = None

try:
    import xml.etree.ElementTree as ET
except ImportError:
    print("Warning: xml.etree.ElementTree not available")
    ET = None


class SVGToGLBExtruder:
    """
    Main class for converting SVG files to extruded GLB meshes.
    """
    
    def __init__(self, stroke_width=0.5, extrude_height=5.0, resolution=0.1):
        """
        Initialize the extruder.
        
        Args:
            stroke_width (float): Width to give to stroke paths when converting to filled shapes
            extrude_height (float): Height to extrude the 2D shapes
            resolution (float): Resolution for path discretization
        """
        self.stroke_width = stroke_width
        self.extrude_height = extrude_height
        self.resolution = resolution
        
    def render_svg_to_png(self, svg_file, png_file, target_width=1600):
        """
        Render SVG to PNG for verification purposes with proper dimensions.
        
        Args:
            svg_file (str): Path to input SVG file
            png_file (str): Path to output PNG file
            target_width (int): Target width in pixels for high resolution
        """
        if not cairosvg:
            print("Warning: Cannot render PNG - cairosvg not available")
            return False
            
        try:
            # First, get SVG dimensions
            svg_width, svg_height = self._get_svg_dimensions(svg_file)
            
            if svg_width and svg_height:
                # Calculate aspect ratio and target dimensions
                aspect_ratio = svg_height / svg_width
                width = target_width
                height = int(target_width * aspect_ratio)
            else:
                # Fallback to default dimensions
                width, height = target_width, int(target_width * 0.75)
            
            # Convert SVG to PNG with proper dimensions
            cairosvg.svg2png(
                url=svg_file,
                write_to=png_file,
                output_width=width,
                output_height=height
            )
            print(f"✓ Rendered PNG: {png_file} ({width}x{height})")
            return True
        except Exception as e:
            print(f"Warning: Failed to render PNG: {e}")
            return False
    
    def _get_svg_dimensions(self, svg_file):
        """
        Extract SVG dimensions from the file.
        
        Args:
            svg_file (str): Path to SVG file
            
        Returns:
            tuple: (width, height) or (None, None) if not found
        """
        try:
            if not ET:
                return None, None
                
            tree = ET.parse(svg_file)
            root = tree.getroot()
            
            # Try to get width and height attributes
            width = root.get('width')
            height = root.get('height')
            
            if width and height:
                # Parse numeric values (remove units like 'px', 'pt', etc.)
                import re
                width_num = re.findall(r'[\d.]+', str(width))
                height_num = re.findall(r'[\d.]+', str(height))
                
                if width_num and height_num:
                    return float(width_num[0]), float(height_num[0])
            
            # Fallback: try viewBox
            viewbox = root.get('viewBox')
            if viewbox:
                values = viewbox.split()
                if len(values) >= 4:
                    return float(values[2]), float(values[3])
                    
        except Exception as e:
            print(f"Could not extract SVG dimensions: {e}")
            
        return None, None
    
    def load_svg_paths(self, svg_file):
        """
        Load SVG and extract 2D paths using svgpathtools.
        
        Args:
            svg_file (str): Path to SVG file
            
        Returns:
            list: List of shapely polygons extracted from SVG
        """
        if not svgpathtools:
            print("Error: svgpathtools not available for SVG parsing")
            return None
            
        try:
            # Parse SVG file
            paths, attributes = svg2paths(svg_file)
            
            if not paths:
                print("Warning: No paths found in SVG file")
                return None
            
            print(f"✓ Found {len(paths)} paths in SVG")
            
            polygons = []
            skipped_paths = []
            path_stats = {'filled': 0, 'stroked': 0, 'empty': 0, 'complex': 0}
            
            # Process each path
            for i, path in enumerate(paths):
                attr = attributes[i] if i < len(attributes) else {}
                
                # Track path characteristics for debugging
                fill = attr.get('fill', '')
                stroke = attr.get('stroke', '')
                path_length = path.length()
                
                if fill and fill.lower() not in ['none', 'transparent', '']:
                    path_stats['filled'] += 1
                elif stroke and stroke.lower() not in ['none', 'transparent', '']:
                    path_stats['stroked'] += 1
                elif path_length > 0:
                    path_stats['complex'] += 1
                else:
                    path_stats['empty'] += 1
                
                try:
                    # Convert SVG path to polygon
                    polygon = self._svg_path_to_polygon(path, attr)
                    if polygon is not None:
                        polygons.append(polygon)
                    else:
                        skipped_paths.append({
                            'index': i,
                            'length': path_length,
                            'fill': fill,
                            'stroke': stroke,
                            'id': attr.get('id', f'path_{i}')
                        })
                except Exception as e:
                    print(f"Warning: Could not process path {i}: {e}")
                    continue
            
            print(f"✓ Successfully processed {len(polygons)} paths into polygons")
            print(f"  Path types: {path_stats['filled']} filled, {path_stats['stroked']} stroked, {path_stats['complex']} complex, {path_stats['empty']} empty")
            
            if len(skipped_paths) > 0:
                print(f"⚠ Skipped {len(skipped_paths)} paths during conversion")
                # Show details for first few skipped paths
                for skip in skipped_paths[:5]:
                    print(f"  - {skip['id']}: length={skip['length']:.1f}, fill='{skip['fill']}', stroke='{skip['stroke']}'")
                if len(skipped_paths) > 5:
                    print(f"  ... and {len(skipped_paths) - 5} more")
            
            return polygons
            
        except Exception as e:
            print(f"Error loading SVG: {e}")
            return None
    
    def _svg_path_to_polygon(self, svg_path, attributes):
        """
        Convert an SVG path to a shapely polygon.
        
        Args:
            svg_path: svgpathtools Path object
            attributes (dict): SVG path attributes
            
        Returns:
            shapely.geometry.Polygon or None
        """
        try:
            # Sample the path to get discrete points with higher resolution for complex shapes
            path_length = svg_path.length()
            if path_length == 0:
                return None
                
            # Use higher sampling for complex paths
            num_samples = max(200, int(path_length / (self.resolution * 0.5)))
            
            points = []
            for i in range(num_samples + 1):
                t = i / num_samples
                try:
                    point = svg_path.point(t)
                    points.append([point.real, point.imag])
                except:
                    continue
            
            if len(points) < 3:
                return None
            
            # Remove duplicate consecutive points
            filtered_points = [points[0]]
            for point in points[1:]:
                if np.linalg.norm(np.array(point) - np.array(filtered_points[-1])) > 1e-3:
                    filtered_points.append(point)
            points = filtered_points
            
            if len(points) < 3:
                return None
            
            # Determine rendering style - be more permissive
            fill = attributes.get('fill', '')
            stroke = attributes.get('stroke', '')
            stroke_width_attr = attributes.get('stroke-width', str(self.stroke_width))
            
            # Parse stroke width safely
            try:
                stroke_width = float(stroke_width_attr)
            except (ValueError, TypeError):
                stroke_width = self.stroke_width
            
            # Check if path is closed (forms a loop)
            is_closed = np.linalg.norm(np.array(points[0]) - np.array(points[-1])) < 1e-2
            
            # Strategy 1: Try as filled shape if fill is specified or path is closed
            if (fill and fill.lower() not in ['none', 'transparent', '']) or is_closed:
                try:
                    # Ensure path is closed
                    if not is_closed and len(points) >= 3:
                        points.append(points[0])
                    
                    polygon = Polygon(points)
                    if polygon.is_valid and polygon.area > 1e-6:
                        return polygon
                    
                    # Try to fix invalid polygon
                    if not polygon.is_valid:
                        fixed_polygon = polygon.buffer(0)
                        if fixed_polygon.is_valid and not fixed_polygon.is_empty:
                            return fixed_polygon
                except Exception as e:
                    pass
            
            # Strategy 2: Try as stroke if stroke is specified or not filled
            if (stroke and stroke.lower() not in ['none', 'transparent', '']) or not fill or fill.lower() in ['none', 'transparent', '']:
                try:
                    from shapely.geometry import LineString
                    line = LineString(points)
                    if line.is_valid and line.length > 0:
                        # Use appropriate stroke width
                        width = max(stroke_width/2, 0.01)  # Minimum width to ensure visibility
                        buffered = buffer(line, distance=width)
                        if buffered.is_valid and not buffered.is_empty:
                            return buffered
                except Exception as e:
                    pass
            
            # Strategy 3: Fallback - treat as default stroke for any remaining valid paths
            try:
                if len(points) >= 3:
                    if is_closed:
                        # Try as polygon first
                        polygon = Polygon(points)
                        if polygon.is_valid and polygon.area > 1e-6:
                            return polygon
                    
                    # Treat as stroke line
                    from shapely.geometry import LineString
                    line = LineString(points)
                    if line.is_valid and line.length > 0:
                        buffered = buffer(line, distance=self.stroke_width/2)
                        if buffered.is_valid and not buffered.is_empty:
                            return buffered
            except Exception as e:
                pass
                
            return None
            
        except Exception as e:
            # Don't print error for every failed path as it clutters output
            return None
    
    def convert_strokes_to_fills(self, polygons):
        """
        Since we're already handling stroke conversion in _svg_path_to_polygon,
        this method just validates and cleans up the polygons.
        
        Args:
            polygons (list): List of shapely polygons
            
        Returns:
            list: Cleaned list of polygons
        """
        if not polygons:
            print("Warning: No polygons provided")
            return None
        
        valid_polygons = []
        for poly in polygons:
            if poly is not None and hasattr(poly, 'is_valid') and poly.is_valid and not poly.is_empty:
                valid_polygons.append(poly)
        
        print(f"✓ Validated {len(valid_polygons)} polygons")
        return valid_polygons
    
    def handle_complex_geometry(self, polygons):
        """
        Handle complex geometry issues like self-intersections and holes.
        
        Args:
            polygons (list): List of shapely polygons
            
        Returns:
            list: Processed polygons ready for extrusion
        """
        try:
            processed_polygons = []
            invalid_count = 0
            fixed_count = 0
            
            for i, poly in enumerate(polygons):
                if isinstance(poly, (Polygon, MultiPolygon)):
                    original_poly = poly
                    
                    # Handle self-intersections by making geometry valid
                    if not poly.is_valid:
                        invalid_count += 1
                        try:
                            poly = poly.buffer(0)  # This often fixes self-intersections
                            if poly.is_valid and not poly.is_empty:
                                fixed_count += 1
                        except:
                            # If buffer fails, try a different approach
                            try:
                                from shapely import make_valid
                                poly = make_valid(original_poly)
                                if poly.is_valid and not poly.is_empty:
                                    fixed_count += 1
                            except:
                                # Skip this polygon if we can't fix it
                                continue
                    
                    # Only simplify if the polygon is relatively large to avoid losing detail
                    if poly.is_valid and hasattr(poly, 'area') and poly.area > 10:
                        # Use gentler simplification to preserve details
                        simplified = simplify(poly, tolerance=0.005)
                        if simplified.is_valid and not simplified.is_empty:
                            poly = simplified
                    
                    if poly.is_valid and not poly.is_empty:
                        if isinstance(poly, MultiPolygon):
                            # Split multipolygons into individual polygons
                            for sub_poly in poly.geoms:
                                if sub_poly.is_valid and not sub_poly.is_empty and hasattr(sub_poly, 'area') and sub_poly.area > 1e-6:
                                    processed_polygons.append(sub_poly)
                        else:
                            if hasattr(poly, 'area') and poly.area > 1e-6:
                                processed_polygons.append(poly)
            
            if invalid_count > 0:
                print(f"✓ Fixed {fixed_count}/{invalid_count} invalid polygons")
            
            print(f"✓ Processed complex geometry: {len(processed_polygons)} valid polygons")
            return processed_polygons
            
        except Exception as e:
            print(f"Error handling complex geometry: {e}")
            return polygons
    
    def polygons_to_path2d(self, polygons):
        """
        Convert shapely polygons back to trimesh Path2D for extrusion.
        
        Args:
            polygons (list): List of shapely polygons
            
        Returns:
            trimesh.path.Path2D: Path2D object ready for extrusion
        """
        try:
            # Combine all polygons into a single path
            vertices = []
            entities = []
            
            for poly in polygons:
                if isinstance(poly, Polygon):
                    # Add exterior coordinates
                    coords = list(poly.exterior.coords)
                    if len(coords) >= 3:
                        start_idx = len(vertices)
                        vertices.extend(coords[:-1])  # Exclude duplicate last point
                        
                        # Create line entities connecting the vertices
                        for i in range(len(coords) - 1):
                            next_i = (i + 1) % (len(coords) - 1)
                            line = trimesh.path.entities.Line([start_idx + i, start_idx + next_i])
                            entities.append(line)
                    
                    # Add holes (interiors)
                    for interior in poly.interiors:
                        coords = list(interior.coords)
                        if len(coords) >= 3:
                            start_idx = len(vertices)
                            vertices.extend(coords[:-1])  # Exclude duplicate last point
                            
                            # Create line entities for the hole (reverse direction)
                            for i in range(len(coords) - 2, -1, -1):
                                next_i = (i - 1) % (len(coords) - 1)
                                line = trimesh.path.entities.Line([start_idx + i, start_idx + next_i])
                                entities.append(line)
            
            if not vertices:
                print("Error: No vertices found in polygons")
                return None
            
            # Convert to numpy array and ensure 2D
            vertices_array = np.array(vertices)
            if vertices_array.shape[1] > 2:
                vertices_array = vertices_array[:, :2]
            
            path_2d = Path2D(entities=entities, vertices=vertices_array)
            print(f"✓ Created Path2D with {len(vertices)} vertices and {len(entities)} entities")
            
            return path_2d
            
        except Exception as e:
            print(f"Error converting polygons to Path2D: {e}")
            return None
    
    def extrude_to_3d(self, path_2d):
        """
        Extrude 2D path to 3D mesh.
        
        Args:
            path_2d (trimesh.path.Path2D): 2D path to extrude
            
        Returns:
            trimesh.Trimesh: 3D mesh
        """
        try:
            # Extrude the path
            mesh = path_2d.extrude(height=self.extrude_height)
            
            if mesh is None:
                print("Error: Extrusion failed")
                return None
            
            # Handle the case where extrusion returns a list of meshes
            if isinstance(mesh, list):
                if len(mesh) == 0:
                    print("Error: Extrusion returned empty list")
                    return None
                elif len(mesh) == 1:
                    mesh = mesh[0]
                else:
                    # Combine multiple meshes
                    print(f"Combining {len(mesh)} meshes...")
                    mesh = trimesh.util.concatenate(mesh)
            
            # Ensure the mesh is valid
            if hasattr(mesh, 'is_watertight'):
                print(f"✓ Extruded to 3D mesh: {len(mesh.vertices)} vertices, {len(mesh.faces)} faces")
                print(f"  Watertight: {mesh.is_watertight}")
                print(f"  Volume: {mesh.volume:.2f}")
            else:
                print(f"✓ Extruded to 3D mesh")
            
            return mesh
            
        except Exception as e:
            print(f"Error extruding to 3D: {e}")
            return None
    
    def save_glb(self, mesh, output_file):
        """
        Save mesh as GLB file.
        
        Args:
            mesh (trimesh.Trimesh): 3D mesh to save
            output_file (str): Output GLB file path
        """
        try:
            # Export as GLB
            mesh.export(output_file)
            print(f"✓ Saved GLB: {output_file}")
            
            # Print some statistics
            file_size = Path(output_file).stat().st_size
            print(f"  File size: {file_size / 1024:.1f} KB")
            
        except Exception as e:
            print(f"Error saving GLB: {e}")
    
    def process_svg(self, svg_file, output_glb, render_png=True):
        """
        Complete pipeline to process SVG to GLB.
        
        Args:
            svg_file (str): Input SVG file path
            output_glb (str): Output GLB file path
            render_png (bool): Whether to render a PNG for verification
        """
        print(f"Processing SVG: {svg_file}")
        print(f"Output GLB: {output_glb}")
        print(f"Extrude height: {self.extrude_height}")
        print(f"Stroke width: {self.stroke_width}")
        print("-" * 50)
        
        # Step 1: Render PNG for verification
        if render_png:
            png_file = str(Path(output_glb).with_suffix('.png'))
            self.render_svg_to_png(svg_file, png_file)
        
        # Step 2: Load SVG paths
        polygons = self.load_svg_paths(svg_file)
        if not polygons:
            print("Failed to load SVG paths")
            return False
        
        # Step 3: Convert strokes to fills and handle complex geometry
        polygons = self.convert_strokes_to_fills(polygons)
        if not polygons:
            print("Failed to process SVG shapes")
            return False
        
        # Step 4: Handle complex geometry (self-intersections, holes)
        polygons = self.handle_complex_geometry(polygons)
        if not polygons:
            print("Failed to handle complex geometry")
            return False
        
        # Step 5: Convert back to Path2D for extrusion
        final_path_2d = self.polygons_to_path2d(polygons)
        if final_path_2d is None:
            print("Failed to create final Path2D")
            return False
        
        # Step 6: Extrude to 3D
        mesh = self.extrude_to_3d(final_path_2d)
        if mesh is None:
            print("Failed to extrude to 3D")
            return False
        
        # Step 7: Save as GLB
        self.save_glb(mesh, output_glb)
        
        print("-" * 50)
        print("✓ Successfully completed SVG to GLB conversion!")
        return True


def main():
    parser = argparse.ArgumentParser(
        description="Convert SVG files to extruded GLB meshes",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python svg_to_glb_extruder.py music_score.svg
  python svg_to_glb_extruder.py input.svg --output custom.glb --extrude-height 10 --stroke-width 1.0
  python svg_to_glb_extruder.py sheet.svg --no-png
        """
    )
    
    parser.add_argument('input_svg', help='Input SVG file path')
    parser.add_argument('--output', '-o', help='Output GLB file path (default: input name with .glb extension)')
    parser.add_argument('--extrude-height', type=float, default=5.0, help='Height to extrude (default: 5.0)')
    parser.add_argument('--stroke-width', type=float, default=0.5, help='Width for stroke paths (default: 0.5)')
    parser.add_argument('--resolution', type=float, default=0.1, help='Path discretization resolution (default: 0.1)')
    parser.add_argument('--no-png', action='store_true', help='Skip PNG rendering')
    
    args = parser.parse_args()
    
    # Validate input file
    input_path = Path(args.input_svg)
    if not input_path.exists():
        print(f"Error: Input file '{args.input_svg}' not found")
        sys.exit(1)
    
    # Determine output path
    if args.output:
        output_path = args.output
    else:
        output_path = str(input_path.with_suffix('.glb'))
    
    # Create extruder and process
    extruder = SVGToGLBExtruder(
        stroke_width=args.stroke_width,
        extrude_height=args.extrude_height,
        resolution=args.resolution
    )
    
    success = extruder.process_svg(
        svg_file=str(input_path),
        output_glb=output_path,
        render_png=not args.no_png
    )
    
    if not success:
        print("Conversion failed!")
        sys.exit(1)


if __name__ == "__main__":
    main() 