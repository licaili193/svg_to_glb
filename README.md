# Music Score Processing Suite

A comprehensive suite of tools for processing and visualizing music scores, including SVG to 3D conversion and OSMD-based music rendering with coordinate verification.

## Features

### 🎵 Music Rendering (OSMD)
- **MusicXML to SVG/PNG**: Render MusicXML files using OpenSheetMusicDisplay
- **Infinite Width**: No line wrapping - continuous horizontal music layout
- **Coordinate Alignment**: SVG and PNG with exactly matching coordinate systems
- **Customizable Output**: Configurable height while maintaining aspect ratio

### 🎯 Coordinate Verification
- **Perfect Alignment**: Verify PNG and SVG coordinate system matching
- **Visual Comparison**: Side-by-side visualization and overlay analysis
- **Detailed Metrics**: Precise measurements and difference calculations

### 🎪 SVG to 3D Conversion
- **SVG to GLB conversion**: Complete pipeline from 2D vector graphics to 3D meshes
- **PNG verification**: Renders PNG from SVG for visual verification
- **GLB export**: Industry-standard 3D format compatible with web browsers, 3D viewers, and game engines

### Advanced Geometry Handling
- **Stroke to fill conversion**: Converts line strokes (staff lines, stems, etc.) to filled shapes for proper 3D extrusion
- **Self-intersection handling**: Automatically fixes self-intersecting polygons (like treble clefs)
- **Hole support**: Properly handles shapes with holes (sharp symbols, half notes, etc.)
- **Complex path processing**: Handles intricate music notation symbols and paths

### Music Notation Specific
- **Staff lines**: Converts thin stroke lines to extruded 3D elements
- **Note symbols**: Handles complex note shapes including stems, flags, and beams
- **Clefs and symbols**: Processes complex musical symbols like treble clefs, time signatures
- **Text elements**: Converts musical text and numbers to 3D

## Installation

### Prerequisites
- **Python 3.8+** for SVG processing and coordinate verification
- **Node.js 16+** for OSMD music rendering
- **npm** package manager

### Install Dependencies

```bash
# Install Python dependencies
pip install -r requirements.txt

# Install Node.js dependencies for OSMD rendering
npm install
```

### Development Setup

```bash
# Clone repository
git clone <repository-url>
cd music-score-processing-suite

# Install all dependencies
pip install -r requirements.txt
npm install

# Test OSMD rendering
node render_music.js --height 600

# Test coordinate verification
python verify_coordinates.py
```

### Alternative Installation (Minimal)

```bash
# Minimal installation (PNG rendering disabled)
pip install trimesh[easy] shapely numpy
```

## Usage

### 🎵 Music Rendering (OSMD)

```bash
# Render MusicXML to SVG and PNG with default settings (600px height)
node render_music.js

# Custom height and output directory
node render_music.js --height 800 --output custom_output

# Specify custom MusicXML file
node render_music.js --input path/to/music.musicxml --height 400
```

**Options:**
- `--input`, `-i`: Input MusicXML file (default: `resources/MozaVeilSample.musicxml`)
- `--height`, `-h`: Desired PNG height in pixels (default: 600)
- `--output`, `-o`: Output directory (default: `output`)

### 🎯 Coordinate Verification

```bash
# Quick verification with default files
python verify_coordinates.py

# Verify custom files
python verify_coordinates.py path/to/music.png path/to/music.svg

# Advanced visualization (requires matplotlib)
python visualize_alignment.py

# Save comparison plot
python visualize_alignment.py --save alignment_check.png
```

### 🎪 SVG to 3D Conversion

```bash
# Convert SVG to GLB with default settings
python svg_to_glb_extruder.py music_score.svg

# Specify custom output file
python svg_to_glb_extruder.py input.svg --output my_model.glb
```

### Advanced Usage

```bash
# Custom extrusion height and stroke width
python svg_to_glb_extruder.py sheet.svg --extrude-height 10.0 --stroke-width 1.5

# Skip PNG rendering for faster processing
python svg_to_glb_extruder.py input.svg --no-png

# Fine-tune path resolution
python svg_to_glb_extruder.py complex.svg --resolution 0.05
```

### Command Line Options

- `input_svg`: Input SVG file path (required)
- `--output`, `-o`: Output GLB file path (default: input name with .glb extension)
- `--extrude-height`: Height to extrude the 2D shapes (default: 5.0)
- `--stroke-width`: Width for stroke paths when converting to fills (default: 0.5)
- `--resolution`: Path discretization resolution (default: 0.1, smaller = more detail)
- `--no-png`: Skip PNG rendering for verification

## Programming Interface

You can also use the tool programmatically:

```python
from svg_to_glb_extruder import SVGToGLBExtruder

# Create extruder with custom settings
extruder = SVGToGLBExtruder(
    stroke_width=1.0,
    extrude_height=8.0,
    resolution=0.05
)

# Process SVG file
success = extruder.process_svg(
    svg_file="music_score.svg",
    output_glb="output.glb",
    render_png=True
)
```

## 🎯 Coordinate Verification Results

The coordinate verification system demonstrates **perfect alignment** between generated PNG and SVG files:

### ✅ Verified Alignment
- **Width match**: 12890 × 12890 pixels (0.0px difference)
- **Height match**: 960 × 960 pixels (0.0px difference)  
- **Origin consistency**: Both use (0,0) at top-left corner
- **Scale matching**: No transformation needed for overlays

### 📊 Technical Details
- **Coordinate System**: Origin at top-left, X→right, Y→bottom
- **SVG Structure**: `<svg width="12890" height="960" viewBox="0 0 5136.735... 382.565...">`
- **PNG Properties**: RGBA format, 12890×960 pixels, infinite width layout
- **Element Count**: 3248 SVG elements processed

### 🔬 Verification Methods
```bash
# Text-based verification (no GUI required)
python verify_coordinates.py
# Output: 🎯 PERFECT ALIGNMENT CONFIRMED!

# Visual verification with matplotlib
python visualize_alignment.py
# Shows side-by-side comparison and overlay analysis
```

## How It Works

### Processing Pipeline

1. **SVG Loading**: Parses SVG file and extracts vector paths using trimesh
2. **PNG Rendering**: Creates PNG for visual verification using cairosvg
3. **Stroke Processing**: Converts stroke paths to filled shapes using shapely buffering
4. **Geometry Cleanup**: Fixes self-intersections and handles complex polygons
5. **Path Reconstruction**: Converts processed shapes back to trimesh Path2D format
6. **3D Extrusion**: Extrudes 2D paths to 3D meshes with specified height
7. **GLB Export**: Saves final 3D mesh in GLB format

### Handling Complex Shapes

#### Self-Intersecting Polygons
- Uses shapely's `buffer(0)` operation to automatically fix self-intersections
- Applies polygon simplification to reduce complexity while preserving shape

#### Shapes with Holes
- Preserves polygon interiors (holes) during processing
- Correctly handles nested shapes and complex topologies

#### Stroke Conversion
- Converts line strokes to filled polygons using buffering
- Ensures thin elements like staff lines appear in the final 3D mesh
- Configurable stroke width for different line thicknesses

## Music Notation Challenges Addressed

### 1. Staff Lines and Stems
Staff lines in music notation are typically drawn as strokes. The tool converts these to filled rectangles that can be properly extruded to 3D.

### 2. Complex Symbols
Musical symbols like treble clefs often have self-intersecting paths. The geometry cleanup process fixes these automatically.

### 3. Note Heads with Holes
Symbols like half notes and whole notes have holes in the center. The tool preserves these holes in the 3D extrusion.

### 4. Overlapping Elements
When multiple musical elements overlap, the tool handles the complex boolean operations required to create a valid 3D mesh.

## Output Formats

### GLB Files
- Binary format, optimized for size and loading speed
- Compatible with web browsers, Three.js, Babylon.js
- Supported by 3D viewers like Blender, 3D Viewer apps
- Can be imported into game engines (Unity, Unreal Engine)

### PNG Files (Verification)
- Generated automatically for visual verification
- Helps ensure SVG was processed correctly
- Can be disabled with `--no-png` flag

## Troubleshooting

### OSMD Rendering Issues

**"Module not found" errors**
- Run `npm install` to install Node.js dependencies
- Ensure you're using Node.js 16+ (check with `node --version`)

**"MusicXML parsing failed"**
- Verify the MusicXML file is valid and properly formatted
- Check file encoding (UTF-8 or UTF-16 supported)

**"Canvas/Graphics errors"**
- On macOS: May need to install Cairo (`brew install cairo`)
- On Linux: Install system graphics libraries
- Try using different Node.js versions if issues persist

### Coordinate Verification Issues

**"PIL not available"**
- Install Pillow: `pip install pillow`

**"matplotlib not available"** 
- Install matplotlib: `pip install matplotlib`
- For text-only verification, use `verify_coordinates.py` instead

**"Cairo graphics library not available"**
- On macOS: `brew install cairo`
- On Ubuntu: `sudo apt-get install libcairo2-dev`
- The basic verification still works without Cairo

### SVG to 3D Conversion Issues

**"Failed to load SVG paths"**
- Ensure SVG file is valid and contains vector paths
- Some SVG files may need preprocessing to simplify complex elements

**"No valid polygons found"**
- SVG may contain only text or raster images
- Try increasing `--stroke-width` for thin line elements

**"Extrusion failed"**
- Input polygons may be too complex
- Try simplifying the SVG or adjusting `--resolution`

### Performance Tips

- Use `--no-png` to skip PNG rendering for faster processing
- Increase `--resolution` value (e.g., 0.2) for faster processing of complex files
- Simplify SVG files beforehand if they contain excessive detail

### Quality Tips

- Decrease `--resolution` value (e.g., 0.05) for higher quality output
- Adjust `--stroke-width` based on the thickness of lines in your SVG
- Use appropriate `--extrude-height` for your intended 3D printing or display scale

## Dependencies

### Python Dependencies
- **trimesh**: 3D mesh processing and GLB export
- **shapely**: 2D geometry operations and polygon handling  
- **numpy**: Numerical computing
- **matplotlib**: Visualization and plotting (coordinate verification)
- **pillow**: Image processing for PNG analysis
- **cairosvg**: SVG to PNG rendering for verification
- **svgpathtools**: Advanced SVG path processing
- **scipy**: Scientific computing for mesh operations

### Node.js Dependencies  
- **opensheetmusicdisplay**: Music notation rendering from MusicXML
- **jsdom**: DOM environment for server-side rendering
- **canvas**: Canvas API implementation for PNG generation
- **sharp**: High-performance image processing

## License

This project is provided as-is for educational and research purposes. Please ensure you have appropriate rights to process any SVG files you use with this tool.

## Contributing

Contributions are welcome! Areas for improvement:
- Support for additional SVG features
- Performance optimizations for large files
- Additional output formats
- Better handling of text elements
- GUI interface

## Examples

The tool has been tested with various types of SVG files:
- Music scores from notation software
- Technical drawings with thin lines
- Logo designs with complex shapes
- Architectural drawings
- Scientific diagrams

For music notation specifically, it handles:
- Staff notation
- Chord symbols
- Dynamic markings
- Articulation marks
- Time signatures and key signatures
- Complex orchestral scores 