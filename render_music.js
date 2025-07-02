#!/usr/bin/env node

import fs from 'fs';
import path from 'path';
import { JSDOM } from 'jsdom';
import { createCanvas, loadImage } from 'canvas';
import sharp from 'sharp';
import pkg from 'opensheetmusicdisplay';
const { OpenSheetMusicDisplay } = pkg;

/**
 * Renders MusicXML to SVG and PNG with infinite width and exact coordinate matching
 */
class MusicRenderer {
    constructor(options = {}) {
        this.desiredHeight = options.height || 400;
        this.outputDir = options.outputDir || 'output';
        this.musicXmlPath = options.musicXmlPath || 'resources/MozaVeilSample.musicxml';
        
        // Ensure output directory exists
        if (!fs.existsSync(this.outputDir)) {
            fs.mkdirSync(this.outputDir, { recursive: true });
        }
    }

    /**
     * Set up JSDOM environment for server-side rendering
     */
    setupDOM() {
        const dom = new JSDOM(`
            <!DOCTYPE html>
            <html>
            <head><title>OSMD Renderer</title></head>
            <body>
                <div id="osmdContainer"></div>
            </body>
            </html>
        `, {
            pretendToBeVisual: true,
            resources: "usable"
        });

        global.window = dom.window;
        global.document = dom.window.document;
        global.HTMLElement = dom.window.HTMLElement;
        global.HTMLDivElement = dom.window.HTMLDivElement;
        global.Node = dom.window.Node;
        global.navigator = dom.window.navigator;
        global.DOMParser = dom.window.DOMParser;
        global.XMLHttpRequest = dom.window.XMLHttpRequest;

        return dom;
    }

    /**
     * Configure OSMD for infinite width rendering
     */
    configureOSMD(osmd) {
        const options = {
            // Disable line wrapping for infinite width
            autoResize: false,
            
            // Disable title and subtitle rendering
            drawTitle: false,
            drawSubtitle: false,
            drawComposer: false,
            drawLyricist: false,
            drawCredits: false,
            
            // Ensure we get a single long line
            drawingParameters: "compacttight",
            
            // Set backend to SVG
            backend: "svg",
            
            // Disable various elements we don't need
            drawMeasureNumbers: false,
            drawMeasureNumbersOnlyAtSystemStart: false,
            
            // Ensure single staff line rendering
            renderSingleHorizontalStaffline: true,
            
            // Disable cursor
            disableCursor: true,
            
            // Set large page format to avoid wrapping
            pageFormat: "Endless"
        };

        osmd.setOptions(options);
        
        // Set custom page format for infinite width (much larger values)
        osmd.setCustomPageFormat(99999, 200); // Very wide, reasonable height in OSMD units
    }

    /**
     * Load and render the MusicXML file
     */
    async loadAndRender() {
        try {
            console.log('Setting up DOM environment...');
            const dom = this.setupDOM();
            
            console.log('Initializing OSMD...');
            const container = document.getElementById('osmdContainer');
            const osmd = new OpenSheetMusicDisplay(container);
            
            console.log('Configuring OSMD for infinite width...');
            this.configureOSMD(osmd);
            
            console.log('Loading MusicXML file...');
            // Read the file as buffer first to detect encoding
            const fileBuffer = fs.readFileSync(this.musicXmlPath);
            let musicXmlContent;
            
            // Check if it's UTF-16 by looking for BOM or encoding declaration
            const bufferStr = fileBuffer.toString('utf8');
            if (bufferStr.includes('encoding="UTF-16"') || fileBuffer[0] === 0xFF || fileBuffer[0] === 0xFE) {
                console.log('Detected UTF-16 encoding, converting...');
                musicXmlContent = fileBuffer.toString('utf16le');
            } else {
                musicXmlContent = bufferStr;
            }
            
            console.log('Rendering music...');
            await osmd.load(musicXmlContent);
            osmd.render();
            
            console.log('Extracting SVG...');
            const svgElement = container.querySelector('svg');
            if (!svgElement) {
                throw new Error('No SVG element found after rendering');
            }
            
            // Get the actual rendered dimensions from SVG attributes
            // getBBox() is not available in JSDOM, so we use viewBox or width/height
            let actualWidth, actualHeight;
            
            const viewBox = svgElement.getAttribute('viewBox');
            if (viewBox) {
                const [x, y, w, h] = viewBox.split(' ').map(Number);
                actualWidth = w;
                actualHeight = h;
            } else {
                actualWidth = parseFloat(svgElement.getAttribute('width')) || 800;
                actualHeight = parseFloat(svgElement.getAttribute('height')) || 600;
            }
            
            // If dimensions are still not reasonable, calculate from content
            if (!actualWidth || !actualHeight || actualWidth < 100 || actualHeight < 50) {
                // Find all elements with x,y positions to estimate bounds
                const allElements = svgElement.querySelectorAll('*[x], *[y], *[transform]');
                let maxX = 0, maxY = 0;
                
                allElements.forEach(el => {
                    const x = parseFloat(el.getAttribute('x')) || 0;
                    const y = parseFloat(el.getAttribute('y')) || 0;
                    maxX = Math.max(maxX, x + 100); // Add some padding
                    maxY = Math.max(maxY, y + 50);
                });
                
                actualWidth = Math.max(actualWidth || 0, maxX, 1000);
                actualHeight = Math.max(actualHeight || 0, maxY, 200);
            }
            
            console.log(`Rendered dimensions: ${actualWidth} x ${actualHeight}`);
            
            // Calculate scale factor to match desired height
            const scaleFactor = this.desiredHeight / actualHeight;
            const scaledWidth = Math.round(actualWidth * scaleFactor);
            
            console.log(`Scaled dimensions: ${scaledWidth} x ${this.desiredHeight}`);
            console.log(`Scale factor: ${scaleFactor}`);
            
            // Update SVG to match PNG dimensions exactly
            svgElement.setAttribute('width', scaledWidth);
            svgElement.setAttribute('height', this.desiredHeight);
            svgElement.setAttribute('viewBox', `0 0 ${actualWidth} ${actualHeight}`);
            
            // Save SVG
            const svgContent = svgElement.outerHTML;
            const svgPath = path.join(this.outputDir, 'music.svg');
            fs.writeFileSync(svgPath, svgContent);
            console.log(`SVG saved to: ${svgPath}`);
            
            // Convert to PNG with exact matching coordinates
            console.log('Converting to PNG...');
            const pngPath = path.join(this.outputDir, 'music.png');
            await this.convertSVGToPNG(svgContent, pngPath, scaledWidth, this.desiredHeight);
            console.log(`PNG saved to: ${pngPath}`);
            
            return {
                svgPath,
                pngPath,
                dimensions: {
                    width: scaledWidth,
                    height: this.desiredHeight,
                    originalWidth: actualWidth,
                    originalHeight: actualHeight,
                    scaleFactor
                }
            };
            
        } catch (error) {
            console.error('Error rendering music:', error);
            throw error;
        }
    }

    /**
     * Convert SVG to PNG using Sharp with exact coordinate matching
     */
    async convertSVGToPNG(svgContent, outputPath, width, height) {
        try {
            const pngBuffer = await sharp(Buffer.from(svgContent))
                .resize(width, height, {
                    fit: 'fill',
                    background: { r: 255, g: 255, b: 255, alpha: 1 }
                })
                .png()
                .toBuffer();
            
            fs.writeFileSync(outputPath, pngBuffer);
        } catch (error) {
            console.error('Error converting SVG to PNG:', error);
            throw error;
        }
    }
}

/**
 * Main execution function
 */
async function main() {
    const args = process.argv.slice(2);
    
    // Parse command line arguments
    let height = 400;
    let inputFile = 'resources/MozaVeilSample.musicxml';
    let outputDir = 'output';
    
    for (let i = 0; i < args.length; i++) {
        if (args[i] === '--height' || args[i] === '-h') {
            height = parseInt(args[i + 1]);
            i++;
        } else if (args[i] === '--input' || args[i] === '-i') {
            inputFile = args[i + 1];
            i++;
        } else if (args[i] === '--output' || args[i] === '-o') {
            outputDir = args[i + 1];
            i++;
        } else if (args[i] === '--help') {
            console.log(`
Usage: node render_music.js [options]

Options:
  --height, -h <number>   Desired height of the PNG output (default: 400)
  --input, -i <path>      Path to MusicXML file (default: resources/MozaVeilSample.musicxml)
  --output, -o <path>     Output directory (default: output)
  --help                  Show this help message

Examples:
  node render_music.js --height 600
  node render_music.js --height 800 --output ./rendered
  node render_music.js -h 500 -i ./my-music.musicxml -o ./results
            `);
            return;
        }
    }
    
    console.log('OSMD Music Renderer');
    console.log('===================');
    console.log(`Input file: ${inputFile}`);
    console.log(`Output directory: ${outputDir}`);
    console.log(`Desired height: ${height}px`);
    console.log('');
    
    // Check if input file exists
    if (!fs.existsSync(inputFile)) {
        console.error(`Error: Input file '${inputFile}' not found`);
        process.exit(1);
    }
    
    try {
        const renderer = new MusicRenderer({
            height,
            musicXmlPath: inputFile,
            outputDir
        });
        
        const result = await renderer.loadAndRender();
        
        console.log('\n✅ Rendering completed successfully!');
        console.log(`📁 Files saved to: ${outputDir}/`);
        console.log(`📄 SVG: ${result.svgPath}`);
        console.log(`🖼️  PNG: ${result.pngPath}`);
        console.log(`📏 Final dimensions: ${result.dimensions.width} x ${result.dimensions.height}px`);
        console.log(`🔍 Scale factor: ${result.dimensions.scaleFactor.toFixed(4)}`);
        
    } catch (error) {
        console.error('❌ Rendering failed:', error.message);
        process.exit(1);
    }
}

// Run the main function
main().catch(console.error); 