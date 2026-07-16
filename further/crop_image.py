#!/usr/bin/env python
import os
from PIL import Image

def crop_center_1000(input_path, output_path):
    if not os.path.exists(input_path):
        raise FileNotFoundError(f"Input image not found: {input_path}")
        
    print(f"Loading image from {input_path}...")
    img = Image.open(input_path)
    w, h = img.size
    print(f"Original image dimensions: {w}x{h}")
    
    # Target crop size
    target_w, target_h = 1000, 1000
    
    if w < target_w or h < target_h:
        raise ValueError(
            f"Image dimensions ({w}x{h}) are smaller than target crop size ({target_w}x{target_h})"
        )
        
    # Calculate crop box coordinates
    left = (w - target_w) // 2
    top = (h - target_h) // 2
    right = left + target_w
    bottom = top + target_h
    
    print(f"Cropping box (left, top, right, bottom): ({left}, {top}, {right}, {bottom})")
    cropped_img = img.crop((left, top, right, bottom))
    
    # Double check crop size
    cw, ch = cropped_img.size
    print(f"Cropped image dimensions: {cw}x{ch}")
    assert cw == target_w and ch == target_h, f"Cropped dimensions {cw}x{ch} do not match target {target_w}x{target_h}!"
    
    cropped_img.save(output_path)
    print(f"Saved cropped image to {output_path}")

if __name__ == "__main__":
    input_img = "/Users/Igniz/Desktop/ICRAR/further/NGC4254_observed_VRI.png"
    output_img = "/Users/Igniz/Desktop/ICRAR/further/NGC4254_observed_VRI_1000_1000.png"
    crop_center_1000(input_img, output_img)
