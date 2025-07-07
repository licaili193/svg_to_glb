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
        Convert an SVG path to a shapely polygon, handling discontinuities properly.
        
        Args:
            svg_path: svgpathtools Path object
            attributes (dict): SVG path attributes
            
        Returns:
            shapely.geometry.Polygon or MultiPolygon or None
        """
        try:
            # Check path continuity and split if needed
            continuous_subpaths = self._split_path_at_discontinuities(svg_path)
            
            if not continuous_subpaths:
                return None
            
            # Process each continuous subpath separately
            valid_polygons = []
            for subpath_segments in continuous_subpaths:
                polygon = self._process_continuous_subpath(subpath_segments, attributes)
                if polygon is not None:
                    valid_polygons.append(polygon)
            
            if not valid_polygons:
                return None
            elif len(valid_polygons) == 1:
                return valid_polygons[0]
            else:
                # Return MultiPolygon for multiple discontinuous parts
                from shapely.geometry import MultiPolygon
                return MultiPolygon(valid_polygons)
                
        except Exception as e:
            # Don't print error for every failed path as it clutters output
            return None
    
    def _split_path_at_discontinuities(self, svg_path):
        """
        Split an SVG path into continuous subpaths at discontinuity points.
        
        Args:
            svg_path: svgpathtools Path object
            
        Returns:
            list: List of continuous subpath segments
        """
        if len(svg_path) == 0:
            return []
        
        continuous_subpaths = []
        current_subpath = []
        discontinuity_count = 0
        
        for i, segment in enumerate(svg_path):
            if i == 0:
                # First segment always starts a new subpath
                current_subpath = [segment]
            else:
                # Check if this segment is continuous with the previous
                prev_segment = svg_path[i-1]
                gap_distance = abs(segment.start - prev_segment.end)
                
                if gap_distance > 1e-6:  # Significant gap = discontinuity
                    discontinuity_count += 1
                    
                    # Save the current subpath and start a new one
                    if current_subpath:
                        continuous_subpaths.append(current_subpath)
                    current_subpath = [segment]
                else:
                    # Continuous - add to current subpath
                    current_subpath.append(segment)
        
        # Don't forget the last subpath
        if current_subpath:
            continuous_subpaths.append(current_subpath)
        
        if discontinuity_count > 0:
            print(f"  ⚠️  Split path with {discontinuity_count} discontinuities into {len(continuous_subpaths)} continuous subpaths")
        
        return continuous_subpaths
    
    def _process_continuous_subpath(self, segments, attributes):
        """
        Process a continuous subpath into a polygon.
        
        Args:
            segments: List of continuous path segments
            attributes (dict): SVG path attributes
            
        Returns:
            shapely.geometry.Polygon or None
        """
        if not segments:
            return None
            
        try:
            # Calculate total length of this subpath
            subpath_length = sum(segment.length() for segment in segments)
            
            if subpath_length == 0:
                return None
            
            # Use adaptive sampling based on subpath complexity and length
            base_samples = max(20, int(subpath_length * 2))  # 2 samples per unit length minimum
            max_samples = 500  # Cap per subpath to prevent excessive computation
            num_samples = min(base_samples, max_samples)
            
            points = []
            last_point = None
            
            # Sample each segment proportionally
            for segment in segments:
                segment_length = segment.length()
                if segment_length == 0:
                    continue
                
                # Number of samples for this segment (proportional to length)
                segment_samples = max(3, int((segment_length / subpath_length) * num_samples))
                
                for i in range(segment_samples):
                    t = i / (segment_samples - 1) if segment_samples > 1 else 0
                    try:
                        complex_point = segment.point(t)
                        point = [float(complex_point.real), float(complex_point.imag)]
                        
                        # Skip invalid points
                        if any(not np.isfinite(coord) for coord in point):
                            continue
                            
                        # Skip duplicate points (more precise threshold)
                        if last_point is not None:
                            distance = np.linalg.norm(np.array(point) - np.array(last_point))
                            if distance < 1e-6:  # Very small threshold for duplicates
                                continue
                        
                        points.append(point)
                        last_point = point
                        
                    except Exception as e:
                        # Skip problematic points but continue sampling
                        continue
            
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
            
            # Check if subpath is closed (forms a loop)
            is_closed = np.linalg.norm(np.array(points[0]) - np.array(points[-1])) < 1e-3
            
            # For very small paths, use a relative threshold instead
            if not is_closed:
                path_bbox_size = max(
                    max(p[0] for p in points) - min(p[0] for p in points),
                    max(p[1] for p in points) - min(p[1] for p in points)
                )
                if path_bbox_size > 0:
                    relative_threshold = path_bbox_size * 0.01  # 1% of bounding box
                    is_closed = np.linalg.norm(np.array(points[0]) - np.array(points[-1])) < relative_threshold
            
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
            return None
    
    def convert_strokes_to_fills(self, polygons):
        """
        Since we're already handling stroke conversion in _svg_path_to_polygon,
        this method just validates and cleans up the polygons.
        Now also handles MultiPolygon objects from discontinuous paths.
        
        Args:
            polygons (list): List of shapely polygons (may include MultiPolygons)
            
        Returns:
            list: Cleaned list of individual polygons
        """
        if not polygons:
            print("Warning: No polygons provided")
            return None
        
        valid_polygons = []
        multipolygon_count = 0
        
        for poly in polygons:
            if poly is not None:
                # Handle MultiPolygon objects from discontinuous paths
                if isinstance(poly, MultiPolygon):
                    multipolygon_count += 1
                    for sub_poly in poly.geoms:
                        if sub_poly is not None and hasattr(sub_poly, 'is_valid') and sub_poly.is_valid and not sub_poly.is_empty:
                            valid_polygons.append(sub_poly)
                # Handle regular Polygon objects
                elif hasattr(poly, 'is_valid') and poly.is_valid and not poly.is_empty:
                    valid_polygons.append(poly)
        
        if multipolygon_count > 0:
            print(f"✓ Expanded {multipolygon_count} discontinuous path(s) into {len(valid_polygons)} individual polygons")
        else:
            print(f"✓ Validated {len(valid_polygons)} polygons")
        return valid_polygons
    
    def handle_complex_geometry(self, polygons):
        """
        Handle complex geometry issues like self-intersections and holes.
        Now also handles discontinuous paths that were split into MultiPolygons.
        
        Args:
            polygons (list): List of shapely polygons (including MultiPolygons from discontinuous paths)
            
        Returns:
            list: Processed polygons ready for extrusion
        """
        try:
            # First, try to resolve overlapping simple strokes (like staff lines)
            # This helps with musical notation where staff lines from adjacent measures overlap
            resolved_polygons = self._resolve_overlapping_simple_strokes(polygons)
            
            processed_polygons = []
            invalid_count = 0
            fixed_count = 0
            discontinuous_count = 0
            
            for i, poly in enumerate(resolved_polygons):
                if isinstance(poly, (Polygon, MultiPolygon)):
                    original_poly = poly
                    
                    # First, handle MultiPolygons from discontinuous paths
                    if isinstance(poly, MultiPolygon):
                        discontinuous_count += 1
                        # Process each sub-polygon separately
                        for sub_poly in poly.geoms:
                            processed_sub, was_invalid, was_fixed = self._process_single_polygon(sub_poly)
                            if was_invalid:
                                invalid_count += 1
                            if was_fixed:
                                fixed_count += 1
                            if processed_sub is not None:
                                if isinstance(processed_sub, MultiPolygon):
                                    # Recursively handle nested MultiPolygons
                                    for nested_poly in processed_sub.geoms:
                                        if nested_poly.is_valid and not nested_poly.is_empty and hasattr(nested_poly, 'area') and nested_poly.area > 1e-6:
                                            processed_polygons.append(nested_poly)
                                else:
                                    if processed_sub.is_valid and not processed_sub.is_empty and hasattr(processed_sub, 'area') and processed_sub.area > 1e-6:
                                        processed_polygons.append(processed_sub)
                    else:
                        # Handle regular Polygons
                        processed_poly, was_invalid, was_fixed = self._process_single_polygon(poly)
                        if was_invalid:
                            invalid_count += 1
                        if was_fixed:
                            fixed_count += 1
                        if processed_poly is not None:
                            if isinstance(processed_poly, MultiPolygon):
                                # Split MultiPolygons into individual polygons
                                for sub_poly in processed_poly.geoms:
                                    if sub_poly.is_valid and not sub_poly.is_empty and hasattr(sub_poly, 'area') and sub_poly.area > 1e-6:
                                        processed_polygons.append(sub_poly)
                            else:
                                if processed_poly.is_valid and not processed_poly.is_empty and hasattr(processed_poly, 'area') and processed_poly.area > 1e-6:
                                    processed_polygons.append(processed_poly)
            
            if invalid_count > 0:
                print(f"✓ Fixed {fixed_count}/{invalid_count} invalid polygons")
            
            if discontinuous_count > 0:
                print(f"✓ Processed {discontinuous_count} discontinuous paths (split at gaps)")
            
            print(f"✓ Processed complex geometry: {len(processed_polygons)} valid polygons")
            return processed_polygons
            
        except Exception as e:
            print(f"Error handling complex geometry: {e}")
            return polygons
    
    def _resolve_overlapping_simple_strokes(self, polygons):
        """
        Resolve overlapping simple rectangular strokes (like staff lines) while preserving 
        complex shapes with holes and self-intersections.
        
        Args:
            polygons (list): List of shapely polygons
            
        Returns:
            list: Processed polygons with simple overlaps resolved
        """
        try:
            if not polygons:
                return polygons
            
            # Separate simple strokes from complex shapes
            simple_strokes = []
            complex_shapes = []
            
            for poly in polygons:
                if isinstance(poly, MultiPolygon):
                    # MultiPolygons are generally complex
                    complex_shapes.append(poly)
                elif isinstance(poly, Polygon):
                    if self._is_simple_rectangular_stroke(poly):
                        simple_strokes.append(poly)
                    else:
                        complex_shapes.append(poly)
                else:
                    complex_shapes.append(poly)
            
            # If we have simple strokes, try to union overlapping ones
            if simple_strokes:
                print(f"  Detected {len(simple_strokes)} simple strokes, {len(complex_shapes)} complex shapes")
                
                # Group simple strokes by approximate orientation and position
                stroke_groups = self._group_similar_strokes(simple_strokes)
                
                # Union overlapping strokes within each group
                resolved_strokes = []
                for group in stroke_groups:
                    if len(group) > 1:
                        # Try to union overlapping strokes in this group
                        try:
                            unioned = unary_union(group)
                            if isinstance(unioned, MultiPolygon):
                                resolved_strokes.extend(unioned.geoms)
                            else:
                                resolved_strokes.append(unioned)
                        except Exception as e:
                            # If union fails, keep individual strokes
                            resolved_strokes.extend(group)
                    else:
                        resolved_strokes.extend(group)
                
                if len(resolved_strokes) < len(simple_strokes):
                    print(f"  ✓ Resolved {len(simple_strokes)} simple strokes into {len(resolved_strokes)} unified strokes")
                
                # Combine resolved strokes with complex shapes
                return resolved_strokes + complex_shapes
            else:
                # No simple strokes detected, return original
                return polygons
                
        except Exception as e:
            print(f"  Warning: Could not resolve overlapping strokes: {e}")
            return polygons
    
    def _is_simple_rectangular_stroke(self, poly):
        """
        Check if a polygon is a simple rectangular stroke (like a staff line).
        This includes both simple rectangles and buffered line strings.
        
        Args:
            poly: Shapely Polygon
            
        Returns:
            bool: True if it's a simple rectangular stroke
        """
        try:
            if not isinstance(poly, Polygon) or not poly.is_valid:
                return False
            
            # Check if it has holes (complex shapes usually have holes)
            if len(poly.interiors) > 0:
                return False
            
            # Get basic characteristics
            bbox = poly.bounds
            width = bbox[2] - bbox[0]
            height = bbox[3] - bbox[1]
            
            if width <= 0 or height <= 0:
                return False
                
            aspect_ratio = max(width, height) / min(width, height)
            
            # Staff lines typically have high aspect ratio (long and thin)
            # This is the key indicator for line-like shapes
            if aspect_ratio > 10:  # Lower threshold to catch more line-like shapes
                return True
            
            # For polygons with exactly 4 vertices, check if it's a proper rectangle
            coords = list(poly.exterior.coords)[:-1]  # Remove duplicate closing point
            if len(coords) == 4:
                # Check if it's roughly rectangular by examining the angles
                # For a rectangle, consecutive edges should be roughly perpendicular
                edges = []
                for i in range(4):
                    p1 = np.array(coords[i])
                    p2 = np.array(coords[(i + 1) % 4])
                    edge = p2 - p1
                    edges.append(edge)
                
                # Check if edges are roughly perpendicular
                tolerance = 0.1  # Allow some tolerance for angles
                for i in range(4):
                    edge1 = edges[i]
                    edge2 = edges[(i + 1) % 4]
                    
                    # Normalize edges
                    edge1_norm = edge1 / np.linalg.norm(edge1) if np.linalg.norm(edge1) > 0 else edge1
                    edge2_norm = edge2 / np.linalg.norm(edge2) if np.linalg.norm(edge2) > 0 else edge2
                    
                    # Check if dot product is close to 0 (perpendicular)
                    dot_product = np.dot(edge1_norm, edge2_norm)
                    if abs(dot_product) > tolerance:
                        return False
                
                # If we reach here, it's a proper rectangle
                return True
            
            return False
            
        except Exception as e:
            return False
    
    def _group_similar_strokes(self, strokes):
        """
        Group similar strokes that might be overlapping (like staff lines).
        
        Args:
            strokes (list): List of simple stroke polygons
            
        Returns:
            list: List of groups, where each group contains similar strokes
        """
        try:
            if not strokes:
                return []
            
            groups = []
            remaining_strokes = strokes.copy()
            
            while remaining_strokes:
                # Start a new group with the first remaining stroke
                current_stroke = remaining_strokes.pop(0)
                current_group = [current_stroke]
                current_bounds = current_stroke.bounds
                
                # Find other strokes that might be in the same "line" (e.g., same staff line)
                i = 0
                while i < len(remaining_strokes):
                    other_stroke = remaining_strokes[i]
                    other_bounds = other_stroke.bounds
                    
                    # Check if strokes are roughly aligned (same y-coordinate range for horizontal lines)
                    y_overlap = min(current_bounds[3], other_bounds[3]) - max(current_bounds[1], other_bounds[1])
                    y_size = max(current_bounds[3] - current_bounds[1], other_bounds[3] - other_bounds[1])
                    
                    # Check if strokes potentially overlap or are adjacent
                    x_gap = max(current_bounds[0], other_bounds[0]) - min(current_bounds[2], other_bounds[2])
                    
                    # If they have significant y-overlap and are close in x, they might be the same staff line
                    if y_size > 0 and y_overlap / y_size > 0.5 and x_gap < 5:  # Allow small gaps
                        current_group.append(other_stroke)
                        remaining_strokes.pop(i)
                        
                        # Update bounds to include this stroke
                        current_bounds = (
                            min(current_bounds[0], other_bounds[0]),
                            min(current_bounds[1], other_bounds[1]),
                            max(current_bounds[2], other_bounds[2]),
                            max(current_bounds[3], other_bounds[3])
                        )
                    else:
                        i += 1
                
                groups.append(current_group)
            
            return groups
            
        except Exception as e:
            print(f"  Warning: Could not group similar strokes: {e}")
            return [[stroke] for stroke in strokes]  # Return each stroke as its own group
    
    def _process_single_polygon(self, poly):
        """
        Process a single polygon to fix issues and simplify if needed.
        
        Args:
            poly: Single Polygon object
            
        Returns:
            tuple: (processed_polygon_or_None, was_invalid, was_fixed)
        """
        try:
            original_poly = poly
            was_invalid = False
            was_fixed = False
            
            # Handle self-intersections by making geometry valid
            if not poly.is_valid:
                was_invalid = True
                try:
                    poly = poly.buffer(0)  # This often fixes self-intersections
                    if poly.is_valid and not poly.is_empty:
                        was_fixed = True
                except:
                    # If buffer fails, try a different approach
                    try:
                        from shapely import make_valid
                        poly = make_valid(original_poly)
                        if poly.is_valid and not poly.is_empty:
                            was_fixed = True
                    except:
                        # Skip this polygon if we can't fix it
                        return None, was_invalid, was_fixed
            
            # Only simplify if the polygon is relatively large to avoid losing detail
            if poly.is_valid and hasattr(poly, 'area') and poly.area > 10:
                # Use gentler simplification to preserve details
                simplified = simplify(poly, tolerance=0.005)
                if simplified.is_valid and not simplified.is_empty:
                    poly = simplified
            
            result = poly if poly.is_valid and not poly.is_empty else None
            return result, was_invalid, was_fixed
            
        except Exception as e:
            return None, False, False
    
    def polygons_to_path2d(self, polygons):
        """
        Convert shapely polygons back to trimesh Path2D for extrusion.
        
        Args:
            polygons (list): List of shapely polygons
            
        Returns:
            trimesh.path.Path2D: Path2D object ready for extrusion
        """
        try:
            # Use trimesh's built-in polygon handling for better reliability
            from trimesh.path.entities import Line
            from shapely.geometry import MultiPolygon
            
            # Combine all polygons for unified processing
            all_vertices = []
            all_entities = []
            
            for poly in polygons:
                if isinstance(poly, Polygon):
                    # Process exterior boundary
                    ext_coords = list(poly.exterior.coords)
                    if len(ext_coords) >= 4:  # At least 3 unique points + 1 closing point
                        # Remove the duplicate closing point
                        ext_coords = ext_coords[:-1]
                        
                        if len(ext_coords) >= 3:
                            start_idx = len(all_vertices)
                            all_vertices.extend([[pt[0], pt[1]] for pt in ext_coords])
                            
                            # Create proper closed loop - connect each point to the next
                            for i in range(len(ext_coords)):
                                next_i = (i + 1) % len(ext_coords)
                                line = Line([start_idx + i, start_idx + next_i])
                                all_entities.append(line)
                    
                    # Process holes (interior boundaries) 
                    for interior in poly.interiors:
                        hole_coords = list(interior.coords)
                        if len(hole_coords) >= 4:  # At least 3 unique points + 1 closing point
                            # Remove the duplicate closing point
                            hole_coords = hole_coords[:-1]
                            
                            if len(hole_coords) >= 3:
                                start_idx = len(all_vertices)
                                all_vertices.extend([[pt[0], pt[1]] for pt in hole_coords])
                                
                                # For holes, maintain the same orientation as given in shapely
                                # (shapely handles the correct orientation internally)
                                for i in range(len(hole_coords)):
                                    next_i = (i + 1) % len(hole_coords)
                                    line = Line([start_idx + i, start_idx + next_i])
                                    all_entities.append(line)
            
            if not all_vertices:
                print("Error: No vertices found in polygons")
                return None
            
            # Convert to numpy array
            vertices_array = np.array(all_vertices, dtype=np.float64)
            
            # Create Path2D with proper entity list
            path_2d = Path2D(entities=all_entities, vertices=vertices_array)
            
            print(f"✓ Created Path2D with {len(all_vertices)} vertices and {len(all_entities)} entities")
            
            # Validate the path and provide debugging info
            try:
                if hasattr(path_2d, 'polygons_closed'):
                    closed_polys = path_2d.polygons_closed
                    if closed_polys is not None and len(closed_polys) > 0:
                        total_area = sum([p.area for p in closed_polys])
                        print(f"  Found {len(closed_polys)} closed polygon(s), total area: {total_area:.2f}")
                    else:
                        print("  Warning: No closed polygons found - check path connectivity")
                        
                if hasattr(path_2d, 'polygons_full'):
                    all_polys = path_2d.polygons_full
                    if all_polys is not None:
                        print(f"  Total polygons (including open): {len(all_polys)}")
                    
            except Exception as e:
                print(f"  Warning: Could not validate path geometry: {e}")
            
            return path_2d
            
        except Exception as e:
            print(f"Error converting polygons to Path2D: {e}")
            print("Attempting alternative approach...")
            
            # Fallback: Use simpler approach without holes
            try:
                vertices = []
                entities = []
                
                for poly in polygons:
                    if isinstance(poly, Polygon):
                        # Only process exterior (skip holes for this fallback)
                        coords = list(poly.exterior.coords)[:-1]  # Remove duplicate
                        if len(coords) >= 3:
                            start_idx = len(vertices)
                            vertices.extend([[pt[0], pt[1]] for pt in coords])
                            
                            # Simple sequential connection
                            for i in range(len(coords)):
                                next_i = (i + 1) % len(coords)
                                entities.append(Line([start_idx + i, start_idx + next_i]))
                
                if vertices:
                    vertices_array = np.array(vertices, dtype=np.float64)
                    path_2d = Path2D(entities=entities, vertices=vertices_array)
                    print(f"✓ Created fallback Path2D with {len(vertices)} vertices")
                    return path_2d
                    
            except Exception as e2:
                print(f"Fallback also failed: {e2}")
                
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

            # Flip y-axis so that SVG y-down becomes y-up in 3D
            mesh.apply_scale([1, -1, 1])
            
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